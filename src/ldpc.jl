using SparseArrays

# Short teaching LDPC: Gallager's regular construction, Fan's array code, and
# flooding sum-product / min-sum. Not a 3GPP NR / LTE encoder. The size cap
# exists because encoding is dense Gaussian elimination over GF(2).
const _LDPC_MAX_DIM = 4096

@doc raw"""
    LDPCCode

Binary linear code stored by a sparse parity-check matrix, for short teaching simulations.

The code is the null space of `H` over GF(2). `n` is the block length, `m` is the number of rows of `H` (including any redundant checks), `rank` is ``\mathrm{rank}_2(H)``, and `k = n - rank` is the dimension. Information bits occupy `info_pos` (increasing column indices) and parity bits occupy `parity_pos` (the pivot columns, also increasing). Row `r` of `P` (`rank` × `k`, entries in ``\{0,1\}``) is the parity equation

```math
c_{\texttt{parity\_pos}[r]} = \sum_{t=1}^{k} P_{r t}\, u_t \pmod{2}.
```

`column_weight` is the common variable-node degree, or `0` when the columns are not regular. `row_weight` is the same quantity for checks. `edge_var[e]` and `edge_check[e]` are the 1-based endpoints of Tanner-graph edge `e`; `var_edges[j]` and `check_edges[i]` list the edges incident on variable `j` and check `i`. Edges are created by scanning the columns of `H` from left to right.

`H` may have more rows than `rank`. Redundant checks stay in the Tanner graph and in [`ldpc_syndrome`](@ref). They do not add parity bits.

Constructors: [`gallager_ldpc`](@ref), [`array_ldpc`](@ref), [`hamming74_ldpc`](@ref), [`ldpc_from_H`](@ref). Encoding is [`ldpc_encode`](@ref). This is a reference implementation for short codes, not a 3GPP NR base graph and not the Richardson–Urbanke linear-time encoder.
"""
struct LDPCCode
    H::SparseMatrixCSC{Int8,Int}
    n::Int
    m::Int
    k::Int
    rank::Int
    info_pos::Vector{Int}
    parity_pos::Vector{Int}
    P::Matrix{UInt8}
    edge_var::Vector{Int}
    edge_check::Vector{Int}
    var_edges::Vector{Vector{Int}}
    check_edges::Vector{Vector{Int}}
    column_weight::Int
    row_weight::Int
end

@doc raw"""
    LDPCResult

Hard decision and diagnostics from [`ldpc_decode_details`](@ref).

# Fields
- `bits`: length-`k` information bits, read from `codeword` at `info_pos`. Positive posterior LLR maps to bit 0.
- `codeword`: length-`n` hard decisions on every variable.
- `posterior`: length-`n` a-posteriori LLRs, ``L_j = \log P(c_j=0\mid \mathbf{y})/P(c_j=1\mid \mathbf{y})`` in exact sum-product on a tree. A positive value favours bit 0. These are not extrinsic LLRs.
- `iterations`: number of flooding rounds actually run. Zero when the channel hard decision is already a codeword and early stopping is on, or when `max_iter` is 0.
- `converged`: true when `codeword` has syndrome 0. Then `codeword == ldpc_encode(code, bits)`. If decoding stops on the iteration limit, the hard vector may not be a codeword, and `bits` are only the systematic coordinates of that vector.
"""
struct LDPCResult
    bits::Vector{Int}
    codeword::Vector{Int}
    posterior::Vector{Float64}
    iterations::Int
    converged::Bool
end

function Base.show(io::IO, code::LDPCCode)
    print(io, "LDPCCode(n=", code.n, ", k=", code.k, ", m=", code.m, ", rank=", code.rank)
    if code.column_weight > 0 && code.row_weight > 0
        print(io, ", regular=(", code.column_weight, ",", code.row_weight, ")")
    end
    print(io, ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", code::LDPCCode)
    return show(io, code)
end

function Base.show(io::IO, result::LDPCResult)
    print(io, "LDPCResult(k=", length(result.bits),
          ", iterations=", result.iterations,
          ", converged=", result.converged, ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", result::LDPCResult)
    return show(io, result)
end

function _check_shape(m::Int, n::Int)
    (m >= 1 && n >= 1) || throw(ArgumentError(
        "parity-check matrix must be non-empty, got $(m)×$(n)",
    ))
    (m <= _LDPC_MAX_DIM && n <= _LDPC_MAX_DIM) || throw(ArgumentError(
        "parity-check matrix is $(m)×$(n); this encoder accepts at most " *
        "$(_LDPC_MAX_DIM) rows and columns (dense Gaussian elimination)",
    ))
    return nothing
end

function _is_prime(p::Int)
    p >= 2 || return false
    (p == 2 || p == 3) && return true
    (p % 2 == 0 || p % 3 == 0) && return false
    d = 5
    while d * d <= p
        (p % d == 0 || p % (d + 2) == 0) && return false
        d += 6
    end
    return true
end

# SplitMix64 (Steele, Lea, and Flood, OOPSLA 2014). Arithmetic is mod 2^64.
function _splitmix64!(state::Ref{UInt64})
    state[] += 0x9E3779B97F4A7C15
    z = state[]
    z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
    z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
    return z ⊻ (z >> 31)
end

function _fisher_yates!(state::Ref{UInt64}, n::Int)
    perm = collect(0:(n - 1))
    for i in n:-1:2
        j = Int(_splitmix64!(state) % UInt64(i)) + 1
        perm[i], perm[j] = perm[j], perm[i]
    end
    return perm
end

function _band_supports(perm::Vector{Int}, μ::Int, dc::Int)
    supports = Vector{Vector{Int}}(undef, μ)
    for r in 0:(μ - 1)
        cols = Vector{Int}(undef, dc)
        for s in 0:(dc - 1)
            cols[s + 1] = perm[r * dc + s + 1] + 1
        end
        sort!(cols)
        supports[r + 1] = cols
    end
    return supports
end

function _gallager_permutation(n::Int, μ::Int, dc::Int, seed::UInt64, used::Set{Vector{Int}})
    state = Ref(seed)
    for _ in 1:10_000
        perm = _fisher_yates!(state, n)
        supports = _band_supports(perm, μ, dc)
        any(s -> s in used, supports) && continue
        length(Set(supports)) == length(supports) || continue
        return perm, supports
    end
    throw(ArgumentError(
        "could not build a Gallager band with distinct row supports for n=$n, " *
        "row weight $dc; choose a larger block length or a smaller column weight",
    ))
end

function _append_band!(rows, cols, perm, t, μ, dc)
    for r in 0:(μ - 1)
        row = (t - 1) * μ + r + 1
        base = r * dc
        for s in 0:(dc - 1)
            push!(rows, row)
            push!(cols, perm[base + s + 1] + 1)
        end
    end
    return nothing
end

function _sparse_ones(rows::Vector{Int}, cols::Vector{Int}, m::Int, n::Int)
    H = sparse(rows, cols, ones(Int8, length(rows)), m, n)
    if nnz(H) != length(rows) || any(!isone, nonzeros(H))
        throw(ArgumentError("duplicate parity-check edges"))
    end
    return H
end

function _constant_degree(degrees::Vector{Int})
    d = degrees[1]
    return all(==(d), degrees) ? d : 0
end

function _gf2_systematic(H::SparseMatrixCSC{Int8,Int})
    m, n = size(H)
    A = zeros(UInt8, m, n)
    rows = rowvals(H)
    vals = nonzeros(H)
    for j in 1:n
        for idx in nzrange(H, j)
            A[rows[idx], j] = vals[idx]
        end
    end
    pivot_col = Int[]
    h = 1
    for col in 1:n
        h > m && break
        piv = 0
        for row in h:m
            if A[row, col] != 0x00
                piv = row
                break
            end
        end
        piv == 0 && continue
        if piv != h
            for j in 1:n
                A[h, j], A[piv, j] = A[piv, j], A[h, j]
            end
        end
        for row in 1:m
            row == h && continue
            A[row, col] == 0x00 && continue
            for j in 1:n
                A[row, j] = xor(A[row, j], A[h, j])
            end
        end
        push!(pivot_col, col)
        h += 1
    end
    rank = length(pivot_col)
    for row in (rank + 1):m
        for j in 1:n
            A[row, j] == 0x00 || error("LDPC Gaussian elimination left a nonzero row")
        end
    end
    for r in 1:rank
        A[r, pivot_col[r]] == 0x01 || error("LDPC pivot is not 1")
        for r2 in 1:rank
            r2 == r && continue
            A[r2, pivot_col[r]] == 0x00 || error("LDPC pivot column was not cleared")
        end
    end
    is_pivot = falses(n)
    is_pivot[pivot_col] .= true
    info_pos = [j for j in 1:n if !is_pivot[j]]
    k = length(info_pos)
    k == n - rank || error("LDPC dimension does not match the pivot count")
    if k < 1
        throw(ArgumentError(
            "parity-check matrix has full column rank; the only word is zero, " *
            "which this constructor rejects",
        ))
    end
    P = Matrix{UInt8}(undef, rank, k)
    for r in 1:rank
        for t in 1:k
            P[r, t] = A[r, info_pos[t]]
        end
    end
    return info_pos, pivot_col, P, rank
end

function _ldpc_code(H::SparseMatrixCSC{Int8,Int})
    m, n = size(H)
    _check_shape(m, n)
    nnz(H) > 0 || throw(ArgumentError("parity-check matrix has no ones"))
    all(isone, nonzeros(H)) || throw(ArgumentError("parity-check entries must be 0 or 1"))
    edge_var = Int[]
    edge_check = Int[]
    var_edges = [Int[] for _ in 1:n]
    check_edges = [Int[] for _ in 1:m]
    rows = rowvals(H)
    for j in 1:n
        for idx in nzrange(H, j)
            i = rows[idx]
            e = length(edge_var) + 1
            push!(edge_var, j)
            push!(edge_check, i)
            push!(var_edges[j], e)
            push!(check_edges[i], e)
        end
    end
    info_pos, parity_pos, P, rank = _gf2_systematic(H)
    coldeg = [length(edges) for edges in var_edges]
    rowdeg = [length(edges) for edges in check_edges]
    return LDPCCode(
        H, n, m, n - rank, rank, info_pos, parity_pos, P,
        edge_var, edge_check, var_edges, check_edges,
        _constant_degree(coldeg), _constant_degree(rowdeg),
    )
end

@doc raw"""
    gallager_ldpc(n::Integer, column_weight::Integer, row_weight::Integer) -> LDPCCode

Deterministic ``(d_v, d_c)``-regular LDPC code of length `n` in Gallager's block form.

# Arguments
- `n`: block length. `row_weight` must divide `n`. At most 4096, and the resulting number of checks must be at most 4096.
- `column_weight`: variable-node degree ``d_v \ge 2``, the number of bands.
- `row_weight`: check-node degree ``d_c \ge 2``.

# Construction
Let ``\mu = n/d_c`` and stack ``d_v`` bands of ``\mu`` rows (Gallager, IRE Trans. Inf. Theory, 8(1):21–28, 1962, Section II). Band 1 is the band matrix: row ``r`` has ones in columns ``(r-1)d_c+1`` through ``r d_c``.

Band ``t = 2,\ldots,d_v`` is that band after a column permutation. The permutation is Fisher–Yates driven by SplitMix64 (Steele, Lea, and Flood, OOPSLA 2014) started at `UInt64(t)`, not by Julia's global RNG. Each draw consumes ``n-1`` mixer outputs; index ``j`` in the swap of length ``i`` is `(z mod i) + 1` (1-based). A draw is rejected, and the stream continues, if any of its row supports is already present. The first accepted draw is kept. At most 10 000 draws are tried. The same `(n, d_v, d_c)` always yields the same matrix.

Gallager's 1962 ensemble does not fix the permutations. The rule above is the definition used by this function, so the matrix is reproducible. It does not claim a designed girth (this is not PEG). Rows are not duplicated.

# Returns
An [`LDPCCode`](@ref) with `column_weight = d_v`, `row_weight = d_c`, and `m = n d_v / d_c` checks. `k = n - \mathrm{rank}_2(H)` may exceed ``n - m`` when the bands are linearly dependent over GF(2). The constructor rejects a matrix of dimension 0.

# References
Gallager, IRE Trans. Inf. Theory, 8(1):21–28, 1962; Gallager, *Low-Density Parity-Check Codes*, MIT Press, 1963. Systematic encoding is Gaussian elimination, Lin and Costello, *Error Control Coding*, 2nd ed., §3.2.
"""
function gallager_ldpc(n::Integer, column_weight::Integer, row_weight::Integer)
    n = Int(n)
    dv = Int(column_weight)
    dc = Int(row_weight)
    n >= 2 || throw(ArgumentError("block length must be ≥ 2, got $n"))
    dv >= 2 || throw(ArgumentError("column_weight must be ≥ 2, got $dv"))
    dc >= 2 || throw(ArgumentError("row_weight must be ≥ 2, got $dc"))
    dc <= n || throw(ArgumentError("row_weight $dc exceeds block length $n"))
    n % dc == 0 || throw(ArgumentError(
        "row_weight=$dc must divide block length n=$n",
    ))
    μ = n ÷ dc
    m = dv * μ
    if n > _LDPC_MAX_DIM || m > _LDPC_MAX_DIM
        throw(ArgumentError(
            "Gallager matrix is $(m)×$(n); this encoder accepts at most " *
            "$(_LDPC_MAX_DIM) rows and columns",
        ))
    end
    rows = Int[]
    cols = Int[]
    sizehint!(rows, m * dc)
    sizehint!(cols, m * dc)
    perm0 = collect(0:(n - 1))
    used = Set{Vector{Int}}()
    for support in _band_supports(perm0, μ, dc)
        push!(used, support)
    end
    _append_band!(rows, cols, perm0, 1, μ, dc)
    for t in 2:dv
        perm, supports = _gallager_permutation(n, μ, dc, UInt64(t), used)
        for support in supports
            push!(used, support)
        end
        _append_band!(rows, cols, perm, t, μ, dc)
    end
    code = _ldpc_code(_sparse_ones(rows, cols, m, n))
    (code.column_weight == dv && code.row_weight == dc) ||
        error("internal: Gallager matrix is not ($dv, $dc)-regular")
    return code
end

@doc raw"""
    array_ldpc(column_weight::Integer, row_weight::Integer, prime::Integer) -> LDPCCode

Fan's array LDPC code: a ``(d_v, d_c)``-regular quasi-cyclic parity-check matrix built from circulant permutation blocks.

# Arguments
- `column_weight`: ``d_v``, with ``2 \le d_v \le p``.
- `row_weight`: ``d_c``, with ``2 \le d_c \le p``.
- `prime`: a prime ``p`` with ``2 \le p`` and ``d_v p \le 4096``, ``d_c p \le 4096``.

# Construction
`H` is a ``d_v \times d_c`` array of ``p \times p`` blocks. Block ``(r, c)``, indices starting at 0, is the circulant permutation matrix ``P^{rc \bmod p}``. The generator ``P`` has a one in column ``(i+1) \bmod p`` of row ``i`` (0-based), so ``P^0 = I``. Block length ``n = d_c p`` and `H` has ``m = d_v p`` rows. Every column has weight ``d_v`` and every row has weight ``d_c``.

This matrix has no 4-cycles (Fan, 2000). The constructor still rejects dimension 0. Rank over GF(2) may be strictly less than ``m``, so the rate ``k/n`` may exceed ``1 - d_v/d_c``.

# References
Fan, J. L., "Array codes as LDPC codes," Proc. 2nd Int. Symp. Turbo Codes and Related Topics, Brest, pp. 543–546, 2000. The same circulant-permutation form is the quasi-cyclic construction in Fossorier, IEEE Trans. Inf. Theory, 50(8):1788–1793, 2004. This is not a 3GPP NR base graph.
"""
function array_ldpc(column_weight::Integer, row_weight::Integer, prime::Integer)
    dv = Int(column_weight)
    dc = Int(row_weight)
    p = Int(prime)
    p >= 2 || throw(ArgumentError("prime must be ≥ 2, got $p"))
    p <= _LDPC_MAX_DIM || throw(ArgumentError(
        "prime must be ≤ $(_LDPC_MAX_DIM), got $p",
    ))
    _is_prime(p) || throw(ArgumentError("array length factor must be prime, got $p"))
    dv >= 2 || throw(ArgumentError("column_weight must be ≥ 2, got $dv"))
    dc >= 2 || throw(ArgumentError("row_weight must be ≥ 2, got $dc"))
    dv <= p || throw(ArgumentError("column_weight $dv exceeds prime $p"))
    dc <= p || throw(ArgumentError("row_weight $dc exceeds prime $p"))
    m = dv * p
    n = dc * p
    if m > _LDPC_MAX_DIM || n > _LDPC_MAX_DIM
        throw(ArgumentError(
            "array matrix is $(m)×$(n); this encoder accepts at most " *
            "$(_LDPC_MAX_DIM) rows and columns",
        ))
    end
    rows = Int[]
    cols = Int[]
    sizehint!(rows, m * dc)
    sizehint!(cols, m * dc)
    for r in 0:(dv - 1)
        for c in 0:(dc - 1)
            shift = mod(r * c, p)
            for i in 0:(p - 1)
                push!(rows, r * p + i + 1)
                push!(cols, c * p + mod(i + shift, p) + 1)
            end
        end
    end
    code = _ldpc_code(_sparse_ones(rows, cols, m, n))
    (code.column_weight == dv && code.row_weight == dc) ||
        error("internal: array matrix is not ($dv, $dc)-regular")
    return code
end

@doc raw"""
    hamming74_ldpc() -> LDPCCode

The binary ``[7, 4]`` Hamming code as a parity-check fixture.

Columns of `H`, from left to right, are the integers ``1,\ldots,7`` written in three bits with the least significant bit in row 1:

```math
H = \begin{bmatrix}
1&0&1&0&1&0&1\\
0&1&1&0&0&1&1\\
0&0&0&1&1&1&1
\end{bmatrix}.
```

The minimum distance is 3. Column weights are not constant (`column_weight == 0`); every row has weight 4. The Tanner graph contains 4-cycles, so flooding sum-product on this matrix is not maximum-likelihood. The matrix is here so tests can enumerate all 16 codewords. It is not an LDPC ensemble. Hard-decision syndrome decoding of the same codewords is [`hamming74_encode`](@ref) / [`hamming74_decode`](@ref), not this belief-propagation fixture.

# References
Lin and Costello, *Error Control Coding*, 2nd ed., §3.2: the parity-check matrix of the length-``2^3-1`` Hamming code has every nonzero vector of ``\mathrm{GF}(2)^3`` as a column.
"""
function hamming74_ldpc()
    H = Int8[
        1 0 1 0 1 0 1
        0 1 1 0 0 1 1
        0 0 0 1 1 1 1
    ]
    return ldpc_from_H(H)
end

@doc raw"""
    ldpc_from_H(H::AbstractMatrix) -> LDPCCode

Build an [`LDPCCode`](@ref) from an explicit parity-check matrix over GF(2).

# Arguments
`H` has integer or Boolean entries in ``\{0,1\}``. Size at most 4096 × 4096, at least one one, and dimension ``k = n - \mathrm{rank}_2(H) \ge 1``. A floating-point matrix is rejected. Duplicate ones in a sparse matrix are rejected.

Zero rows are kept as idle checks. A column of zeros is an unprotected information bit. Degree-1 checks are allowed: the decoder pins that variable to bit 0. Regular Gallager and array codes do not contain degree-1 checks.

Gaussian elimination scans columns from the left and takes the first available pivot, so `info_pos` and `parity_pos` are deterministic. This is the systematic form of Lin and Costello, 2nd ed., §3.2, with parity bits left in the pivot columns rather than packed at the end of the block.
"""
function ldpc_from_H(H::SparseMatrixCSC{<:Integer})
    m, n = size(H)
    _check_shape(m, n)
    I, J, V = findnz(H)
    for t in eachindex(V)
        V[t] == 1 || throw(ArgumentError(
            "parity-check entries must be 0 or 1, got $(V[t]) at ($(I[t]), $(J[t]))",
        ))
    end
    isempty(I) && throw(ArgumentError("parity-check matrix has no ones"))
    return _ldpc_code(_sparse_ones(I, J, m, n))
end

function ldpc_from_H(H::AbstractMatrix{<:Integer})
    m, n = size(H)
    _check_shape(m, n)
    I = Int[]
    J = Int[]
    for j in 1:n
        for i in 1:m
            v = H[i, j]
            if v == 0
                continue
            elseif v == 1
                push!(I, i)
                push!(J, j)
            else
                throw(ArgumentError(
                    "parity-check entries must be 0 or 1, got $v at ($i, $j)",
                ))
            end
        end
    end
    isempty(I) && throw(ArgumentError("parity-check matrix has no ones"))
    return _ldpc_code(_sparse_ones(I, J, m, n))
end

function ldpc_from_H(H::AbstractMatrix{Bool})
    return ldpc_from_H(Int8.(H))
end

function ldpc_from_H(::AbstractMatrix{<:AbstractFloat})
    throw(ArgumentError("parity-check matrix must have integer entries 0 or 1"))
end

@doc raw"""
    ldpc_column_degrees(code::LDPCCode) -> Vector{Int}

Variable-node degrees, column by column of `code.H`. For a regular construction every entry equals `code.column_weight`. A value `code.column_weight == 0` means the columns are not all the same degree; the vector returned here is then the degree distribution as a list of `n` integers.
"""
function ldpc_column_degrees(code::LDPCCode)
    return [length(edges) for edges in code.var_edges]
end

@doc raw"""
    ldpc_row_degrees(code::LDPCCode) -> Vector{Int}

Check-node degrees, one per row of `code.H`, including redundant and degree-0 rows. See [`ldpc_column_degrees`](@ref).
"""
function ldpc_row_degrees(code::LDPCCode)
    return [length(edges) for edges in code.check_edges]
end

@doc raw"""
    ldpc_encode(code::LDPCCode, bits::AbstractVector{<:Integer}) -> Vector{Int}

Systematic encoding. `bits` has length `code.k` and entries in ``\{0,1\}``. Bit `t` is written to column `code.info_pos[t]` and the parity bits are `P * bits` over GF(2).

The map is deterministic. Every returned vector has syndrome 0 under [`ldpc_syndrome`](@ref), including when `H` has redundant rows. Encoding does not allocate a generator matrix; [`ldpc_generator_matrix`](@ref) is the same linear map, with `c = u G` in the row-vector convention.

# References
Lin and Costello, 2nd ed., §3.2. The Richardson–Urbanke ``O(n)`` encoder (IEEE Trans. Inf. Theory, 47(2):638–656, 2001) is not implemented; Gaussian elimination is done once, when the code is constructed.
"""
function ldpc_encode(code::LDPCCode, bits::AbstractVector{<:Integer})
    k = code.k
    length(bits) == k || throw(ArgumentError(
        "expected $(k) information bits, got $(length(bits))",
    ))
    for b in bits
        (b == 0 || b == 1) || throw(ArgumentError("bits must be 0 or 1"))
    end
    c = zeros(Int, code.n)
    @inbounds for i in 1:k
        c[code.info_pos[i]] = Int(bits[i])
    end
    P = code.P
    @inbounds for r in 1:code.rank
        s = 0
        for t in 1:k
            s ⊻= Int(P[r, t]) & Int(bits[t])
        end
        c[code.parity_pos[r]] = s
    end
    return c
end

@doc raw"""
    ldpc_generator_matrix(code::LDPCCode) -> Matrix{Int}

Binary generator matrix, size `k` × `n`. Row-vector convention over GF(2): the codeword for information row `u` is `u * G`. Columns `info_pos` form ``I_k``. Columns `parity_pos` form ``P^{\mathsf T}``. Rows are a basis of the code.
"""
function ldpc_generator_matrix(code::LDPCCode)
    G = zeros(Int, code.k, code.n)
    @inbounds for i in 1:code.k
        G[i, code.info_pos[i]] = 1
        for r in 1:code.rank
            G[i, code.parity_pos[r]] = Int(code.P[r, i])
        end
    end
    return G
end

@doc raw"""
    ldpc_syndrome(code::LDPCCode, bits::AbstractVector{<:Integer}) -> Vector{Int}

`H` times `bits` over GF(2), length `m`. `bits` is a length-`n` vector with entries in ``\{0,1\}``. The result is the zero vector if and only if `bits` is a codeword. Redundant rows of `H` are included.
"""
function ldpc_syndrome(code::LDPCCode, bits::AbstractVector{<:Integer})
    length(bits) == code.n || throw(ArgumentError(
        "length must be $(code.n), got $(length(bits))",
    ))
    for b in bits
        (b == 0 || b == 1) || throw(ArgumentError("bits must be 0 or 1"))
    end
    s = zeros(Int, code.m)
    @inbounds for e in eachindex(code.edge_var)
        if bits[code.edge_var[e]] == 1
            s[code.edge_check[e]] ⊻= 1
        end
    end
    return s
end

function _syndrome_is_zero(code::LDPCCode, hard::Vector{Int})
    @inbounds for i in 1:code.m
        s = 0
        for e in code.check_edges[i]
            s ⊻= hard[code.edge_var[e]]
        end
        iszero(s) || return false
    end
    return true
end

# Binary box-plus: 2 atanh(tanh(x/2) tanh(y/2)), in the Jacobian form
# sign(x) sign(y) min(|x|,|y|) + log1p(exp(-(|x|+|y|))) - log1p(exp(-||x|-|y||)).
# Richardson and Urbanke, Modern Coding Theory, Cambridge, 2008, §2.5 / the
# check-node update; Ryan and Lin, Channel Codes, Cambridge, 2009, Ch. 5.
function _llr_boxplus(x::Float64, y::Float64)
    x_inf = isinf(x)
    y_inf = isinf(y)
    if x_inf || y_inf
        if x_inf && y_inf
            return (x < 0) ⊻ (y < 0) ? -Inf : Inf
        elseif x_inf
            return x < 0 ? -y : y
        else
            return y < 0 ? -x : x
        end
    end
    ax = abs(x)
    ay = abs(y)
    mag = min(ax, ay) + log1p(exp(-(ax + ay))) - log1p(exp(-abs(ax - ay)))
    mag < 0 && (mag = 0.0)
    return (x < 0) ⊻ (y < 0) ? -mag : mag
end

function _check_update!(incoming, outgoing, d::Int, algorithm::Symbol,
                        alpha::Float64, beta::Float64, prefix, suffix)
    if d == 1
        # Empty box-plus. A degree-1 check forces its variable to 0.
        outgoing[1] = Inf
        return nothing
    end
    if algorithm === :sum_product
        prefix[1] = incoming[1]
        for i in 2:d
            prefix[i] = _llr_boxplus(prefix[i - 1], incoming[i])
        end
        suffix[d] = incoming[d]
        for i in (d - 1):-1:1
            suffix[i] = _llr_boxplus(incoming[i], suffix[i + 1])
        end
        outgoing[1] = suffix[2]
        outgoing[d] = prefix[d - 1]
        for i in 2:(d - 1)
            outgoing[i] = _llr_boxplus(prefix[i - 1], suffix[i + 1])
        end
        return nothing
    end
    total_neg = false
    amin = Inf
    second = Inf
    imin = 1
    @inbounds for i in 1:d
        x = incoming[i]
        x < 0 && (total_neg = !total_neg)
        a = abs(x)
        if a < amin
            second = amin
            amin = a
            imin = i
        elseif a < second
            second = a
        end
    end
    @inbounds for i in 1:d
        mag = i == imin ? second : amin
        if beta != 0
            mag = mag > beta ? mag - beta : 0.0
        end
        mag *= alpha
        neg = total_neg ⊻ (incoming[i] < 0)
        outgoing[i] = neg ? -mag : mag
    end
    return nothing
end

function _pack_result(code::LDPCCode, hard::Vector{Int}, posterior::Vector{Float64},
                      iterations::Int, converged::Bool)
    return LDPCResult(hard[code.info_pos], hard, posterior, iterations, converged)
end

@doc raw"""
    ldpc_decode_details(code::LDPCCode, llr::AbstractVector{<:AbstractFloat};
                        algorithm::Symbol=:normalized_min_sum, max_iter::Integer=20,
                        alpha::Real=0.75, beta::Real=0.15, early_stop::Bool=true) -> LDPCResult

Flooding iterative decoding of an [`LDPCCode`](@ref) from channel LLRs.

# Arguments
- `llr`: length `code.n`, finite, one LLR per code bit,

```math
L_j = \log\frac{p(y_j\mid c_j=0)}{p(y_j\mid c_j=1)}.
```

  A positive value favours bit 0, the same convention as [`demodulate_soft`](@ref CommToolkit.Modulation.demodulate_soft). Integer vectors are rejected: hard bits are not LLRs. For BPSK from [`modulate`](@ref CommToolkit.Modulation.modulate) (``0\mapsto +1``, ``1\mapsto -1``) in AWGN,

```math
L = \frac{2y}{\sigma^2}, \qquad \sigma^2 = \frac{N_0}{2},
```

  with the noise variance `awgn` uses. Coded BPSK of symbol energy 1 has information-bit energy ``E_b = n/k``, so

```math
\sigma^2 = \frac{n/k}{2\,E_b/N_0}
```

  and the call is `awgn(symbols, snr_db; energy_per_bit = n/k)`. There is no separate a-priori input; fold a prior into `llr` before calling.
- `algorithm`:
  - `:sum_product` — log-domain sum-product. The check-to-variable message is the box-plus of the other incoming messages, ``2\operatorname{atanh}\prod\tanh(L/2)``.
  - `:min_sum` — the min-sum (BP-based) approximation of that box-plus.
  - `:normalized_min_sum` — min-sum scaled by `alpha` (default). `alpha` must lie in ``(0, 1]``.
  - `:offset_min_sum` — min-sum with an offset `beta` subtracted inside the magnitude, then clipped at 0. `beta` must be ``\ge 0``.
- `max_iter`: maximum number of flooding rounds. `0` returns the channel hard decision and does not iterate.
- `early_stop`: when true, return as soon as the hard-decision syndrome is zero, including before the first round.

`alpha` is used only by `:normalized_min_sum`. `beta` is used only by `:offset_min_sum`. The other algorithm ignores both.

# Schedule
One round updates every check from the current variable-to-check messages, then every variable. Variable-to-check messages are initialized to the channel LLR. This is the flooding schedule on the Tanner graph (Kschischang, Frey, and Loeliger, IEEE Trans. Inf. Theory, 47(2):498–519, 2001; the parity-check form goes back to Gallager, 1962). It is not a layered or shuffled schedule.

On a tree, `:sum_product` equals bitwise MAP after enough rounds (the graph diameter). With cycles it is approximate. Min-sum replaces the box-plus magnitude by a minimum (Fossorier, Mihaljević, and Imai, IEEE Trans. Commun., 47(5):673–680, 1999). Normalized min-sum multiplies that minimum by ``\alpha`` (Chen and Fossorier, IEEE Trans. Commun., 50(3):406–414, 2002). Offset min-sum subtracts ``\beta`` (Chen and Fossorier, IEEE Commun. Lett., 6(5):208–210, 2002; Chen, Dholakia, Eleftheriou, Fossorier, and Hu, IEEE Trans. Commun., 53(8):1288–1299, 2005).

The defaults `alpha = 0.75` and `beta = 0.15` are fixed constants in the range those papers select by density evolution for regular ensembles. They are not re-optimized for each matrix.

If `early_stop` is true and the channel hard decision is already a codeword, `iterations` is 0 and `posterior` equals the channel LLR. Callers who want a-posteriori LLRs should pass `early_stop=false`. Hard decisions use `posterior .< 0` (a tie at 0 maps to bit 0).

# References
Also Ryan and Lin, *Channel Codes: Classical and Modern*, Cambridge University Press, 2009, Chapter 5, for the syndrome stopping rule and the log-domain check update. Not implemented: NR base graphs, rate matching, CRC attachment, HARQ, puncturing, shortening, and Turbo iteration.
"""
function ldpc_decode_details(code::LDPCCode, llr::AbstractVector{<:AbstractFloat};
                             algorithm::Symbol=:normalized_min_sum,
                             max_iter::Integer=20,
                             alpha::Real=0.75,
                             beta::Real=0.15,
                             early_stop::Bool=true)
    length(llr) == code.n || throw(ArgumentError(
        "LLR length must be $(code.n), got $(length(llr))",
    ))
    all(isfinite, llr) || throw(ArgumentError("channel LLRs must be finite"))
    max_iter = Int(max_iter)
    max_iter >= 0 || throw(ArgumentError("max_iter must be ≥ 0, got $max_iter"))
    if algorithm === :sum_product || algorithm === :min_sum
        alpha = 1.0
        beta = 0.0
    elseif algorithm === :normalized_min_sum
        alpha_ok = isfinite(alpha) && 0 < alpha <= 1
        alpha_ok || throw(ArgumentError(
            "alpha must lie in (0, 1] for :normalized_min_sum, got $alpha",
        ))
        alpha = Float64(alpha)
        beta = 0.0
    elseif algorithm === :offset_min_sum
        beta_ok = isfinite(beta) && beta >= 0
        beta_ok || throw(ArgumentError(
            "beta must be ≥ 0 for :offset_min_sum, got $beta",
        ))
        alpha = 1.0
        beta = Float64(beta)
    else
        throw(ArgumentError(
            "algorithm must be :sum_product, :min_sum, :normalized_min_sum, " *
            "or :offset_min_sum, got $(repr(algorithm))",
        ))
    end

    channel = Float64.(llr)
    posterior = copy(channel)
    hard = Vector{Int}(undef, code.n)
    @inbounds for j in 1:code.n
        hard[j] = channel[j] < 0 ? 1 : 0
    end
    converged = _syndrome_is_zero(code, hard)
    if max_iter == 0 || (early_stop && converged)
        return _pack_result(code, hard, posterior, 0, converged)
    end

    nedges = length(code.edge_var)
    v2c = Vector{Float64}(undef, nedges)
    c2v = zeros(Float64, nedges)
    @inbounds for e in 1:nedges
        v2c[e] = channel[code.edge_var[e]]
    end
    maxdeg = 1
    for edges in code.check_edges
        maxdeg = max(maxdeg, length(edges))
    end
    incoming = Vector{Float64}(undef, maxdeg)
    outgoing = Vector{Float64}(undef, maxdeg)
    prefix = Vector{Float64}(undef, maxdeg)
    suffix = Vector{Float64}(undef, maxdeg)

    iterations = 0
    @inbounds for it in 1:max_iter
        iterations = it
        for i in 1:code.m
            edges = code.check_edges[i]
            d = length(edges)
            d == 0 && continue
            for t in 1:d
                incoming[t] = v2c[edges[t]]
            end
            _check_update!(incoming, outgoing, d, algorithm, alpha, beta, prefix, suffix)
            for t in 1:d
                c2v[edges[t]] = outgoing[t]
            end
        end
        for j in 1:code.n
            edges = code.var_edges[j]
            s = channel[j]
            finite = isfinite(s)
            for e in edges
                msg = c2v[e]
                s += msg
                finite &= isfinite(msg)
            end
            posterior[j] = s
            hard[j] = s < 0 ? 1 : 0
            if finite
                for e in edges
                    v2c[e] = s - c2v[e]
                end
            else
                for e in edges
                    acc = channel[j]
                    for e2 in edges
                        e2 == e && continue
                        acc += c2v[e2]
                    end
                    v2c[e] = acc
                end
            end
        end
        converged = _syndrome_is_zero(code, hard)
        if early_stop && converged
            break
        end
    end
    return _pack_result(code, hard, posterior, iterations, converged)
end

function ldpc_decode_details(::LDPCCode, llr::AbstractVector{<:Integer}; kwargs...)
    _ = kwargs
    throw(ArgumentError(
        "LDPC decoding expects channel LLRs (positive towards bit 0), not hard bits; " *
        "got $(typeof(llr))",
    ))
end

function ldpc_decode_details(::LDPCCode, llr::AbstractVector{<:Real}; kwargs...)
    _ = kwargs
    throw(ArgumentError(
        "channel LLRs must be floating-point values, got $(eltype(llr))",
    ))
end

@doc raw"""
    ldpc_decode(code::LDPCCode, llr::AbstractVector{<:AbstractFloat}; kwargs...) -> Vector{Int}

Information bits from [`ldpc_decode_details`](@ref). Keyword arguments are the same. The default algorithm is `:normalized_min_sum`.

Positive channel LLR means bit 0. The returned bits are the systematic coordinates of the hard decision, length `code.k`.
"""
function ldpc_decode(code::LDPCCode, llr::AbstractVector; kwargs...)
    return ldpc_decode_details(code, llr; kwargs...).bits
end
