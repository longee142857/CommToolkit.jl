# Short teaching Arikan polar codes: Bhattacharyya information-set
# construction, non-systematic and systematic encoding, and
# successive-cancellation (SC) decoding. Not a 3GPP NR polar code: no CRC
# attachment, no list (SCL) decoding, no rate matching, no NR reliability
# sequence. The block-length cap keeps frames short for teaching.

const _POLAR_MAX_N = 256

@doc raw"""
    PolarCode

Short teaching Arikan polar code of block length `n` (a power of two, at most 256) and dimension `k`.

`info_pos` lists the `k` information positions (1-based, increasing) of the pre-transform vector `u`; the complement is frozen to 0. `frozen[i]` is true when position `i` of `u` is frozen. `design_snr_db` is the binary-input AWGN ``E_s/N_0`` in decibels at symbol energy 1 that the information-set construction was designed for, and `z` holds the Bhattacharyya parameter of each of the `n` synthesized sub-channels, in the natural order the successive-cancellation decoder walks.

The codeword for `u` is `x = u G_N` over GF(2) with `G_N = F^{\otimes n}`, `F = [[1, 0], [1, 1]]`, in natural order, without the bit-reversal permutation. Arikan's original `G_N = B_N F^{\otimes n}` differs only by a permutation of the codeword coordinates, which a memoryless channel does not see.

Constructor: [`polar_code`](@ref). Encoding is [`polar_encode`](@ref); decoding is [`polar_sc_decode`](@ref). This is a teaching Arikan polar code with SC decoding, not a 3GPP NR polar code: no CRC attachment, no list decoding, no rate matching, and the information set is not the NR reliability sequence.
"""
struct PolarCode
    n::Int
    k::Int
    info_pos::Vector{Int}
    frozen::BitVector
    design_snr_db::Float64
    z::Vector{Float64}
end

function Base.show(io::IO, code::PolarCode)
    print(io, "PolarCode(n=", code.n, ", k=", code.k,
          ", design_snr_db=", code.design_snr_db, ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", code::PolarCode)
    return show(io, code)
end

# Bhattacharyya parameters of the n synthesized sub-channels, in the natural
# order the SC decoder walks. Each stage replaces every z by its minus and
# plus children, Z^- = 2Z - Z^2 and Z^+ = Z^2 (Arikan 2009, Proposition 5).
function _polar_bhattacharyya(n::Int, z0::Float64)
    z = [z0]
    while length(z) < n
        m = length(z)
        nxt = Vector{Float64}(undef, 2 * m)
        @inbounds for i in 1:m
            zi = z[i]
            nxt[2 * i - 1] = 2 * zi - zi * zi
            nxt[2 * i] = zi * zi
        end
        z = nxt
    end
    return z
end

@doc raw"""
    polar_code(n::Integer, k::Integer; design_snr_db::Real=0.0) -> PolarCode

Construct a short teaching Arikan polar code of block length `n` and dimension `k`.

# Arguments
- `n`: block length, a power of two with ``2 \le n \le 256``. The cap keeps frames short for teaching. There is no rate matching, so block lengths that are not powers of two do not exist here.
- `k`: number of information bits per codeword, ``1 \le k \le n``.
- `design_snr_db`: binary-input AWGN ``E_s/N_0`` in decibels at symbol energy 1, used only to pick the information set. Must be finite. It is a design point, not an operating point: the encoder and decoder accept any channel.

# Construction
The Bhattacharyya parameter of the binary-input AWGN channel at ``E_s/N_0`` (``\sigma^2 = N_0/2`` per dimension, the same convention as [`awgn`](@ref CommToolkit.Channel.awgn)) is

```math
Z_0 = \int_{-\infty}^{\infty} \sqrt{p(y \mid 0)\,p(y \mid 1)}\;dy = e^{-E_s/N_0}.
```

Arikan's recursion propagates it through the ``\log_2 N`` combining stages,

```math
Z^- = 2Z - Z^2, \qquad Z^+ = Z^2,
```

producing one parameter per synthesized sub-channel, in the natural order the successive-cancellation decoder walks. The `k` indices with the smallest `z` become `info_pos`; exact ties go to the lower index. `Z^+ = Z^2` is exact for every symmetric binary-input memoryless channel; `Z^- = 2Z - Z^2` holds with equality on the binary erasure channel and is otherwise an upper bound on the exact Bhattacharyya parameter (Arikan, IEEE Trans. Inf. Theory, 55(7):3051–3073, 2009, Proposition 5). Running the recursion on the AWGN ``Z_0`` is the standard teaching design rule: exact for a BEC, a heuristic elsewhere. It is not the 3GPP NR polar reliability sequence of TS 38.212 §5.3.1.2, and no CRC, interleaver, or parity check is added on top.

The construction is deterministic and uses no random numbers: the same `(n, k, design_snr_db)` always yields the same `info_pos` and `z`.

# Returns
A [`PolarCode`](@ref). See also [`polar_encode`](@ref) and [`polar_sc_decode`](@ref).
"""
function polar_code(n::Integer, k::Integer; design_snr_db::Real=0.0)
    n = Int(n)
    k = Int(k)
    n >= 2 || throw(ArgumentError("block length must be ≥ 2, got $n"))
    ispow2(n) || throw(ArgumentError(
        "block length must be a power of two, got $n",
    ))
    n <= _POLAR_MAX_N || throw(ArgumentError(
        "block length is $n; this teaching construction accepts at most " *
        "$_POLAR_MAX_N",
    ))
    1 <= k <= n || throw(ArgumentError(
        "dimension must satisfy 1 ≤ k ≤ n=$n, got $k",
    ))
    snr = Float64(design_snr_db)
    isfinite(snr) || throw(ArgumentError(
        "design_snr_db must be finite, got $design_snr_db",
    ))
    z = _polar_bhattacharyya(n, exp(-10^(snr / 10)))
    order = sortperm(z; alg=MergeSort) # stable: exact ties keep the lower index
    info_pos = sort(order[1:k])
    frozen = trues(n)
    frozen[info_pos] .= false
    return PolarCode(n, k, info_pos, frozen, snr, z)
end

# In-place polar transform x = u F^{⊗n} over GF(2), F = [1 0; 1 1], natural
# order (no bit-reversal permutation). An involution: two applications return
# the input, since F^2 = I over GF(2).
function _polar_transform!(x::AbstractVector{<:Integer})
    n = length(x)
    len = 1
    while len < n
        step = 2 * len
        @inbounds for i in 1:step:n
            for j in 0:(len - 1)
                x[i + j] ⊻= x[i + j + len]
            end
        end
        len = step
    end
    return x
end

@doc raw"""
    polar_encode(code::PolarCode, bits::AbstractVector{<:Integer}; systematic::Bool=false) -> Vector{Int}

Map `code.k` information bits to a length-`code.n` codeword `x = u G_N` over GF(2), with `G_N = F^{\otimes n}`, `F = [[1, 0], [1, 1]]`, in natural order.

`bits` has length `code.k` and entries in ``\{0,1\}``.

- `systematic = false` (default): `u` is zero except `u[code.info_pos] = bits`. This is Arikan's original `x = u G_N`.
- `systematic = true`: the frozen positions of `u` are still zero, but the information positions are solved so that the codeword itself carries the message, `x[code.info_pos] == bits`. The solve is the triangular system `u_A T_{A,A} = \mathrm{bits}` over GF(2), where `A` is `code.info_pos` and `T = F^{\otimes n}`. `T` is lower triangular with unit diagonal (`T[r, c] = 1` exactly when the 1-bits of `c - 1` are a subset of those of `r - 1`), so every principal submatrix `T_{A,A}` is invertible and back-substitution runs from the last information position down. This is the systematic-encoding idea of Arikan (IEEE Commun. Lett., 15(8):860–862, 2011) written for the natural-order transform; the decoder reads the message from the systematic positions of the re-encoded estimate.

The transform is the ``\mathcal{O}(N \log N)`` butterfly: for stage `s = 0, \ldots, \log_2 N - 1`, every pair `(x_i, x_{i+2^s})` with bit `s` of `i` clear becomes `(x_i \oplus x_{i+2^s},\; x_{i+2^s})`. It is an involution, `G_N G_N = I`. Both forms produce words of the same linear code; only the message-to-codeword map differs. There is no CRC attachment and no rate matching.
"""
function polar_encode(code::PolarCode, bits::AbstractVector{<:Integer}; systematic::Bool=false)
    k = code.k
    length(bits) == k || throw(ArgumentError(
        "expected $(k) information bits, got $(length(bits))",
    ))
    for b in bits
        (b == 0 || b == 1) || throw(ArgumentError("bits must be 0 or 1"))
    end
    A = code.info_pos
    u = zeros(Int, code.n)
    if systematic
        # u_A T_{A,A} = bits over GF(2), T = F^{⊗n}. Column c of the product
        # sees only rows r ≥ c (T is lower triangular), so solving from the
        # largest information position down is back-substitution. Positions
        # outside A stay zero.
        @inbounds for t in k:-1:1
            c = A[t]
            s = Int(bits[t])
            c0 = c - 1
            for r in (t + 1):k
                ri = A[r]
                u[ri] == 0 && continue
                (c0 & ~(ri - 1)) == 0 && (s ⊻= 1)
            end
            u[c] = s
        end
    else
        @inbounds for t in 1:k
            u[A[t]] = Int(bits[t])
        end
    end
    return _polar_transform!(u)
end

# Recursive SC on the tree of the natural-order transform. L[ℓ] and B[ℓ] are
# per-level workspaces; the node at level ℓ covers m positions starting at
# offset (0-based). A leaf decides its bit (frozen to 0), an internal node
# computes the box-plus LLRs for its left child, saves the left decisions,
# applies the g update (1 - 2û)α_i + α_{i+m}, recurses right, and assembles
# its own partial-sum segment (left ⊕ right, right).
function _polar_sc!(uhat::Vector{Int}, L, B, level::Int, m::Int, offset::Int,
                    frozen::BitVector)
    if m == 1
        bit = frozen[offset + 1] ? 0x00 : (L[level][1] < 0 ? 0x01 : 0x00)
        uhat[offset + 1] = Int(bit)
        B[level][1] = bit
        return nothing
    end
    h = m >> 1
    α = L[level]
    β = L[level + 1]
    @inbounds for i in 1:h
        β[i] = _llr_boxplus(α[i], α[i + h])
    end
    _polar_sc!(uhat, L, B, level + 1, h, offset, frozen)
    bl = B[level]
    bc = B[level + 1]
    @inbounds for i in 1:h
        bl[i] = bc[i]
    end
    @inbounds for i in 1:h
        β[i] = (bl[i] == 0x01 ? -α[i] : α[i]) + α[i + h]
    end
    _polar_sc!(uhat, L, B, level + 1, h, offset + h, frozen)
    @inbounds for i in 1:h
        bl[i] ⊻= bc[i]
        bl[i + h] = bc[i]
    end
    return nothing
end

@doc raw"""
    polar_sc_decode(code::PolarCode, llr::AbstractVector{<:AbstractFloat}; systematic::Bool=false) -> Vector{Int}

Successive-cancellation (SC) decoding of a [`PolarCode`](@ref) from channel LLRs.

# Arguments
- `llr`: length `code.n`, finite, one LLR per code bit,

```math
L_j = \log\frac{p(y_j \mid x_j = 0)}{p(y_j \mid x_j = 1)}.
```

  A positive value favours bit 0, the same convention as [`demodulate_soft`](@ref CommToolkit.Modulation.demodulate_soft). Integer vectors are rejected: hard bits are not LLRs. For BPSK from [`modulate`](@ref CommToolkit.Modulation.modulate) (``0 \mapsto +1``, ``1 \mapsto -1``) in AWGN, `L = 2y/\sigma^2` with `\sigma^2 = N_0/2`. Coded BPSK of symbol energy 1 has information-bit energy ``E_b = n/k``, so the channel call is `awgn(symbols, snr_db; energy_per_bit = n/k)` and `\sigma^2 = (n/k)/(2\,E_b/N_0)`, exactly as in the LDPC interface.
- `systematic`: must match the encoder. With `false` (default) the returned bits are `\hat{u}[\texttt{info\_pos}]`. With `true` the estimate is re-encoded, `\hat{x} = \hat{u} G_N`, and `\hat{x}[\texttt{info\_pos}]` is returned. `\hat{u}` is zero on the frozen positions, so `\hat{x}` is always a codeword and the two readings of a correct frame are identical.

# Algorithm
Arikan's SC decoder on the tree of the natural-order transform. A node covering `2m` positions with LLRs `\alpha` computes, for `i = 1, \ldots, m`,

```math
\beta_i = \alpha_i \boxplus \alpha_{i+m},
```

for its left child, where `a \boxplus b = 2\,\mathrm{atanh}(\tanh(a/2)\tanh(b/2))` is the exact log-domain check update, the same box-plus the LDPC decoder uses. After the left subtree is decided,

```math
\beta'_i = (1 - 2\hat{u}^{\mathrm{left}}_i)\,\alpha_i + \alpha_{i+m}
```

feeds the right child. Frozen leaves are decided as 0; a tie at LLR 0 maps to bit 0. One frame costs `\mathcal{O}(N \log N)` LLR operations. SC is not maximum-likelihood: decided bits are never revisited, and an early error propagates down the tree. There is no list (SCL) decoding, no CRC-aided path selection, and no belief-propagation variant here.

# References
Arikan, IEEE Trans. Inf. Theory, 55(7):3051–3073, 2009. This is a teaching SC decoder, not a 3GPP NR polar decoder.
"""
function polar_sc_decode(code::PolarCode, llr::AbstractVector{<:AbstractFloat};
                         systematic::Bool=false)
    n = code.n
    length(llr) == n || throw(ArgumentError(
        "LLR length must be $(n), got $(length(llr))",
    ))
    all(isfinite, llr) || throw(ArgumentError("channel LLRs must be finite"))
    levels = trailing_zeros(n) + 1
    L = [Vector{Float64}(undef, n >> (ℓ - 1)) for ℓ in 1:levels]
    B = [Vector{UInt8}(undef, n >> (ℓ - 1)) for ℓ in 1:levels]
    L[1] .= llr
    uhat = Vector{Int}(undef, n)
    _polar_sc!(uhat, L, B, 1, n, 0, code.frozen)
    systematic && _polar_transform!(uhat)
    return uhat[code.info_pos]
end

function polar_sc_decode(::PolarCode, llr::AbstractVector{<:Integer}; kwargs...)
    _ = kwargs
    throw(ArgumentError(
        "polar SC decoding expects channel LLRs (positive towards bit 0), not hard bits; " *
        "got $(typeof(llr))",
    ))
end

function polar_sc_decode(::PolarCode, llr::AbstractVector{<:Real}; kwargs...)
    _ = kwargs
    throw(ArgumentError(
        "channel LLRs must be floating-point values, got $(eltype(llr))",
    ))
end
