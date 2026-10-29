# GF(2^m) for m = 2, 3, 4, the [7, 4] Hamming syndrome decoder, the binary
# primitive BCH code of length 15, dimension 7, designed distance 5, and
# narrow-sense primitive Reed–Solomon codes over those fields.
# Not a general-n BCH constructor and not a long (CCSDS) Reed–Solomon library.

const _GF_PRIMS = Dict{Int,Int}(
    2 => 0b111,      # x^2 + x + 1
    3 => 0b1011,     # x^3 + x + 1
    4 => 0b10011,    # x^4 + x + 1
)

@doc raw"""
    GF2m(m::Integer)

Galois field ``\mathrm{GF}(2^m)`` for ``m\in\{2,3,4\}``.

Elements are integers ``0,\ldots,2^m-1``. The integer is the bitmask of a polynomial over GF(2): bit 0 is the constant term, bit ``m-1`` is the coefficient of ``x^{m-1}``. Addition is [`gf_add`](@ref). Multiplication uses the discrete log of a fixed primitive element ``\alpha=x`` (integer 2) and one of the primitive polynomials

```math
x^2+x+1,\qquad x^3+x+1,\qquad x^4+x+1.
```

``m=4`` is the field under [`bch157_encode`](@ref), and these fields are the symbol fields of [`rs_code`](@ref). There is no ``\mathrm{GF}(2^8)`` table and no Conway polynomial search. Build the field with `GF2m(m)`; the fields of the struct are the log and antilog tables and are not a second constructor.

# References
Lidl and Niederreiter, *Finite Fields*, 2nd ed., Cambridge, 1997, §2.5, for the polynomial representation. The three polynomials are the standard primitive ones of those degrees (Lin and Costello, *Error Control Coding*, 2nd ed., Appendix A).
"""
struct GF2m
    m::Int
    prim::Int
    order::Int
    exp::Vector{Int}
    log::Vector{Int}
end

function GF2m(m::Integer)
    mm = Int(m)
    prim = get(_GF_PRIMS, mm, 0)
    prim == 0 && throw(ArgumentError(
        "GF(2^m) is implemented for m = 2, 3, and 4, got $mm",
    ))
    order = 1 << mm
    n = order - 1
    exp = Vector{Int}(undef, 2n)
    logt = fill(-1, order)
    x = 1
    hi = order
    for i in 0:(n - 1)
        logt[x] < 0 || error("internal: polynomial for m=$mm is not primitive")
        exp[i + 1] = x
        logt[x] = i
        x <<= 1
        if (x & hi) != 0
            x ⊻= prim
        end
        x &= n
    end
    x == 1 || throw(ArgumentError("primitive polynomial for m=$mm did not cycle"))
    for i in 0:(n - 1)
        exp[n + i + 1] = exp[i + 1]
    end
    return GF2m(mm, prim, order, exp, logt)
end

function Base.show(io::IO, F::GF2m)
    print(io, "GF(2^", F.m, ")")
end

function _gf_elem(F::GF2m, a::Integer)
    aa = Int(a)
    (0 <= aa < F.order) || throw(ArgumentError(
        "GF(2^$(F.m)) element out of range: $aa",
    ))
    return aa
end

@doc raw"""
    gf_add(F::GF2m, a::Integer, b::Integer) -> Int

Sum of two elements of `F`. The sum is bitwise XOR. `F` is used only to reject integers outside ``0,\ldots,2^m-1``.
"""
function gf_add(F::GF2m, a::Integer, b::Integer)
    return _gf_elem(F, a) ⊻ _gf_elem(F, b)
end

@doc raw"""
    gf_mul(F::GF2m, a::Integer, b::Integer) -> Int

Product of two elements of `F`. A zero factor returns 0. The nonzero product is ``\alpha^{\log a+\log b}``.
"""
function gf_mul(F::GF2m, a::Integer, b::Integer)
    aa = _gf_elem(F, a)
    bb = _gf_elem(F, b)
    (aa == 0 || bb == 0) && return 0
    return F.exp[F.log[aa] + F.log[bb] + 1]
end

@doc raw"""
    gf_div(F::GF2m, a::Integer, b::Integer) -> Int

Quotient ``a \cdot b^{-1}`` in `F`. Dividing by zero throws `ArgumentError`. The numerator 0 returns 0.
"""
function gf_div(F::GF2m, a::Integer, b::Integer)
    aa = _gf_elem(F, a)
    bb = _gf_elem(F, b)
    bb == 0 && throw(ArgumentError("division by zero in GF(2^$(F.m))"))
    aa == 0 && return 0
    n = F.order - 1
    e = mod(F.log[aa] - F.log[bb], n)
    return F.exp[e + 1]
end

@doc raw"""
    gf_inv(F::GF2m, a::Integer) -> Int

Multiplicative inverse of a nonzero element of `F`. The inverse of 0 throws `ArgumentError`.
"""
function gf_inv(F::GF2m, a::Integer)
    aa = _gf_elem(F, a)
    aa == 0 && throw(ArgumentError("inverse of zero in GF(2^$(F.m))"))
    n = F.order - 1
    return F.exp[n - F.log[aa] + 1]
end

@doc raw"""
    gf_pow(F::GF2m, a::Integer, e::Integer) -> Int

Power ``a^e`` in `F`. The exponent may be negative; it is reduced modulo ``2^m-1`` for ``a \neq 0``. ``0^0 = 1`` and ``0^e = 0`` for ``e > 0``. A negative power of 0 throws `ArgumentError`.
"""
function gf_pow(F::GF2m, a::Integer, e::Integer)
    aa = _gf_elem(F, a)
    ee = Int(e)
    if aa == 0
        ee > 0 && return 0
        ee == 0 && return 1
        throw(ArgumentError("negative power of zero in GF(2^$(F.m))"))
    end
    n = F.order - 1
    return F.exp[mod(F.log[aa] * ee, n) + 1]
end

const _GF16 = GF2m(4)

# g(x) = x^8 + x^7 + x^6 + x^4 + 1, MSB first (coefficient of x^8 first).
const _BCH157_G = [1, 1, 1, 0, 1, 0, 0, 0, 1]
const _BCH157_N = 15
const _BCH157_K = 7

function _gf2_mod(dividend::Vector{Int}, g::Vector{Int})
    r = copy(dividend)
    deg = length(g) - 1
    n = length(r)
    @inbounds for i in 1:(n - deg)
        if r[i] != 0
            for j in 1:length(g)
                r[i + j - 1] ⊻= g[j]
            end
        end
    end
    return r[(n - deg + 1):n]
end

function _horner(F::GF2m, coefs::AbstractVector{<:Integer}, point::Int)
    acc = 0
    for c in coefs
        acc = xor(gf_mul(F, acc, point), Int(c) & 1)
    end
    return acc
end

@doc raw"""
    hamming74_encode(data::AbstractVector{<:Integer}) -> Vector{Int}

Encode 4 bits into the binary ``[7,4]`` Hamming code, minimum distance 3.

The parity-check matrix is the one in [`hamming74_ldpc`](@ref): column ``j`` (1-based) is the integer ``j`` with its least significant bit in row 1. Data bits `(d_1,d_2,d_3,d_4)` occupy codeword positions `(3,5,6,7)`. The parity bits are

```math
c_1 = d_1+d_2+d_4,\qquad
c_2 = d_1+d_3+d_4,\qquad
c_4 = d_2+d_3+d_4
```

over GF(2), and the codeword is ``(c_1,c_2,d_1,c_4,d_2,d_3,d_4)``. `hamming74_encode([1, 0, 1, 1])` is `[0, 1, 1, 0, 0, 1, 1]`.

`data` has length 4 and entries in ``\{0,1\}``. This is a hard algebraic encoder. It does not call `ldpc_encode`, and the systematic positions are not the pivot order of that Gaussian elimination. The set of 16 codewords is the same.
"""
function hamming74_encode(data::AbstractVector{<:Integer})
    length(data) == 4 || throw(ArgumentError(
        "Hamming(7, 4) takes 4 information bits, got $(length(data))",
    ))
    d = _binary_bits(data)
    d1, d2, d3, d4 = d[1], d[2], d[3], d[4]
    c1 = d1 ⊻ d2 ⊻ d4
    c2 = d1 ⊻ d3 ⊻ d4
    c4 = d2 ⊻ d3 ⊻ d4
    return [c1, c2, d1, c4, d2, d3, d4]
end

@doc raw"""
    hamming74_decode(received::AbstractVector{<:Integer}) -> Vector{Int}

Syndrome decode one word of the binary ``[7,4]`` Hamming code. Corrects one error and returns the 4 data bits in the order [`hamming74_encode`](@ref) accepts.

The syndrome is ``H\mathbf{r}`` with the matrix of [`hamming74_ldpc`](@ref). Packed as an integer with row 1 in the least significant bit, a nonzero syndrome equals the 1-based index of the single column that matches it, and that bit is flipped. The zero syndrome leaves the word unchanged.

Two errors do not decode to the transmitted data: their syndrome equals some third column, that bit is flipped, and the result is a different codeword. This decoder always returns a codeword. It is not an erasure decoder and it has no failure flag. `received` has length 7 and entries in ``\{0,1\}``.
"""
function hamming74_decode(received::AbstractVector{<:Integer})
    length(received) == 7 || throw(ArgumentError(
        "Hamming(7, 4) word length must be 7, got $(length(received))",
    ))
    r = _binary_bits(received)
    s1 = r[1] ⊻ r[3] ⊻ r[5] ⊻ r[7]
    s2 = r[2] ⊻ r[3] ⊻ r[6] ⊻ r[7]
    s3 = r[4] ⊻ r[5] ⊻ r[6] ⊻ r[7]
    syn = s1 | (s2 << 1) | (s3 << 2)
    if syn != 0
        r[syn] ⊻= 1
    end
    return [r[3], r[5], r[6], r[7]]
end

@doc raw"""
    bch157_encode(msg::AbstractVector{<:Integer}) -> Vector{Int}

Encode 7 bits into the binary primitive BCH code of length 15, dimension 7, and designed distance 5.

The code is the narrow-sense BCH code over [`GF2m`](@ref)`(4)` with roots ``\alpha,\alpha^2,\alpha^3,\alpha^4``, ``\alpha`` a root of ``x^4+x+1``. The generator is

```math
g(x)=(x^4+x+1)(x^4+x^3+x^2+x+1)=x^8+x^7+x^6+x^4+1.
```

The codeword is systematic. Entry 1 of the returned vector is the coefficient of ``x^{14}``, and entry 15 is the coefficient of ``x^0``. The first 7 entries are `msg`; the last 8 are the remainder of ``m(x)\,x^{8}`` divided by ``g(x)``. `msg` has length 7 and entries in ``\{0,1\}``.

The coefficients of ``g(x)`` itself, padded on the high-degree side, are the codeword `[0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1]`. The all-zero message encodes to the all-zero word.

Designed distance 5 corrects 2 errors. The true minimum distance of this code is 5 (Lin and Costello, 2nd ed., §6.1, the length-15 double-error-correcting BCH code). This is not a Reed–Solomon code and not a BCH code of any other length.
"""
function bch157_encode(msg::AbstractVector{<:Integer})
    length(msg) == _BCH157_K || throw(ArgumentError(
        "BCH(15, 7) takes $(_BCH157_K) information bits, got $(length(msg))",
    ))
    bits = _binary_bits(msg)
    shifted = vcat(bits, zeros(Int, _BCH157_N - _BCH157_K))
    parity = _gf2_mod(shifted, _BCH157_G)
    return vcat(bits, parity)
end

function _bch157_locate(rx::Vector{Int})
    F = _GF16
    α = 2
    S1 = _horner(F, rx, gf_pow(F, α, 1))
    S3 = _horner(F, rx, gf_pow(F, α, 3))
    if S1 == 0 && S3 == 0
        return Int[], 0
    end
    if S1 == 0
        return Int[], -1
    end
    Λ1 = S1
    Λ2 = gf_div(F, xor(S3, gf_pow(F, S1, 3)), S1)
    locs = Int[]
    for i in 0:(_BCH157_N - 1)
        x = i == 0 ? 1 : gf_pow(F, α, _BCH157_N - i)
        val = xor(1, gf_mul(F, Λ1, x), gf_mul(F, Λ2, gf_mul(F, x, x)))
        if val == 0
            push!(locs, _BCH157_N - i)
        end
    end
    deg = Λ2 == 0 ? 1 : 2
    if length(locs) != deg
        return Int[], -1
    end
    return locs, length(locs)
end

@doc raw"""
    bch157_decode(received::AbstractVector{<:Integer}) -> Vector{Int}

Decode one word of the binary BCH code from [`bch157_encode`](@ref). Corrects every pattern of at most 2 errors and returns the 7 systematic bits.

Syndromes are ``S_j=r(\alpha^j)`` for ``j=1,3`` (the even syndromes of a binary code are powers of ``S_1``). The locator for at most two errors is the Peterson solution

```math
\Lambda(x)=1+S_1 x+\frac{S_3+S_1^3}{S_1}x^2
```

when ``S_1\neq 0``, and ``\Lambda(x)=1`` when ``S_1=S_3=0`` (Lin and Costello, 2nd ed., §6.2–§6.3; Peterson, IRE Trans. Inf. Theory, 1960). Roots are found by a Chien search. A root ``\alpha^{-i}`` marks an error in the coefficient of ``x^i``.

If ``S_1=0`` but ``S_3\neq 0``, or the number of roots is not the degree of ``\Lambda``, the pattern is not a correctable pattern of weight at most 2. The function then returns the 7 received systematic bits with no bits flipped. There is no separate failure flag. A pattern of 3 or more errors may also be miscorrected to a different codeword; designed distance 5 does not promise to correct it.

`received` has length 15 and entries in ``\{0,1\}``, in the same coefficient order as `bch157_encode`.
"""
function bch157_decode(received::AbstractVector{<:Integer})
    length(received) == _BCH157_N || throw(ArgumentError(
        "BCH(15, 7) word length must be $(_BCH157_N), got $(length(received))",
    ))
    rx = _binary_bits(received)
    locs, nerr = _bch157_locate(rx)
    if nerr > 0
        for i in locs
            rx[i] ⊻= 1
        end
    end
    return rx[1:_BCH157_K]
end

# ---------- Narrow-sense primitive Reed–Solomon codes ----------

# Evaluate a polynomial over F whose coefficients are given highest degree
# first (the codeword order of rs_encode and bch157_encode). Unlike _horner,
# the coefficients are arbitrary field elements.
function _gf_eval(F::GF2m, coefs::AbstractVector{Int}, point::Int)
    acc = 0
    for c in coefs
        acc = xor(gf_mul(F, acc, point), c)
    end
    return acc
end

# Same evaluation for a polynomial stored lowest degree first (Berlekamp–Massey
# output and derived polynomials).
function _gf_eval_low(F::GF2m, coefs::Vector{Int}, point::Int)
    acc = 0
    for i in length(coefs):-1:1
        acc = xor(gf_mul(F, acc, point), coefs[i])
    end
    return acc
end

# Formal derivative over GF(2^m): p x^p differentiates to c x^(p-1) for odd p
# and to 0 for even p.
function _gf_deriv_low(coefs::Vector{Int})
    length(coefs) == 1 && return Int[0]
    out = zeros(Int, length(coefs) - 1)
    for p in 1:(length(coefs) - 1)
        if isodd(p)
            out[p] = coefs[p + 1]
        end
    end
    return out
end

# Remainder of a high-degree-first polynomial divided by the monic g over F.
function _gf_poly_mod(F::GF2m, dividend::Vector{Int}, g::Vector{Int})
    r = copy(dividend)
    deg = length(g) - 1
    n = length(r)
    @inbounds for i in 1:(n - deg)
        c = r[i]
        if c != 0
            for j in 2:(deg + 1)
                r[i + j - 1] ⊻= gf_mul(F, c, g[j])
            end
        end
    end
    return r[(n - deg + 1):n]
end

# g(x) = prod_{j=1}^{nroots} (x + α^j), α = x, stored highest degree first.
function _rs_generator(F::GF2m, nroots::Int)
    g = Int[1]
    for j in 1:nroots
        root = gf_pow(F, 2, j)
        h = Vector{Int}(undef, length(g) + 1)
        h[1] = g[1]
        for i in 2:length(g)
            h[i] = xor(g[i], gf_mul(F, g[i - 1], root))
        end
        h[length(g) + 1] = gf_mul(F, g[end], root)
        g = h
    end
    return g
end

# S_j = r(α^j) for j = 1..d, r stored highest degree first.
function _rs_syndromes(F::GF2m, word::Vector{Int}, d::Int)
    S = Vector{Int}(undef, d)
    x = 2
    for j in 1:d
        S[j] = _gf_eval(F, word, x)
        x = gf_mul(F, x, 2)
    end
    return S
end

# Massey's shift-register synthesis on the syndrome sequence S_1..S_d. Returns
# the error-locator polynomial Λ (lowest degree first) and its degree L.
function _berlekamp_massey(F::GF2m, S::Vector{Int})
    Λ = Int[1]
    B = Int[1]
    L = 0
    m = 1
    b = 1
    for r in 0:(length(S) - 1)
        d = S[r + 1]
        for i in 1:L
            i + 1 > length(Λ) && break
            d ⊻= gf_mul(F, Λ[i + 1], S[r + 1 - i])
        end
        if d == 0
            m += 1
            continue
        end
        coef = gf_div(F, d, b)
        T = copy(Λ)
        for _ in 1:(m + length(B) - length(Λ))
            push!(Λ, 0)
        end
        for j in 1:length(B)
            B[j] == 0 && continue
            Λ[j + m] ⊻= gf_mul(F, coef, B[j])
        end
        if 2L <= r
            L = r + 1 - L
            B = T
            b = d
            m = 1
        else
            m += 1
        end
    end
    return Λ, L
end

# Ω(x) ≡ S(x) Λ(x) (mod x^d) with S(x) = S_1 + S_2 x + ..., lowest degree first.
function _rs_key_equation(F::GF2m, S::Vector{Int}, Λ::Vector{Int})
    Ω = zeros(Int, length(S))
    for i in 1:length(S)
        S[i] == 0 && continue
        for j in 1:length(Λ)
            Λ[j] == 0 && continue
            deg = i + j - 2
            deg >= length(S) && continue
            Ω[deg + 1] ⊻= gf_mul(F, S[i], Λ[j])
        end
    end
    return Ω
end

@doc raw"""
    RSCode

Narrow-sense primitive Reed–Solomon code over [`GF2m`](@ref). Built by [`rs_code`](@ref), shortened by [`rs_shorten`](@ref), and used by [`rs_encode`](@ref) and [`rs_decode`](@ref).

The fields are the parameter list of the code, not a second constructor. `F` is the symbol field; `n` and `k` are the effective length and dimension after shortening; `t` is the guaranteed correction radius ``\lfloor (n-k)/2 \rfloor``; `g` is the generator polynomial with the coefficient of the highest power first; `s` is the number of shortened symbols.
"""
struct RSCode
    F::GF2m
    n::Int
    k::Int
    t::Int
    g::Vector{Int}
    s::Int
end

function Base.show(io::IO, C::RSCode)
    print(io, "RS(", C.n, ", ", C.k, ") over GF(2^", C.F.m, ")")
    C.s == 0 || print(io, ", shortened by ", C.s)
    return nothing
end

@doc raw"""
    rs_code(m::Integer, k::Integer) -> RSCode

Narrow-sense primitive Reed–Solomon ``(n, k)`` code over ``\mathrm{GF}(2^m)`` with ``n = 2^m - 1``.

The field is [`GF2m`](@ref)`(m)` with ``m \in \{2, 3, 4\}``, so ``n \in \{3, 7, 15\}``. The generator has the ``n - k`` consecutive roots ``\alpha, \alpha^2, \ldots, \alpha^{n-k}``, ``\alpha = x``:

```math
g(x) = \prod_{j=1}^{n-k} \bigl(x + \alpha^j\bigr).
```

By the BCH bound the minimum distance is at least ``n - k + 1``, and a Reed–Solomon code meets the Singleton bound with equality, so ``d = n - k + 1`` and the guaranteed correction radius is

```math
t = \left\lfloor \frac{n - k}{2} \right\rfloor .
```

`k` satisfies ``1 \le k \le n - 1``. ``k = n - 1`` gives ``t = 0``: [`rs_decode`](@ref) then checks the syndrome but corrects nothing. These are short teaching codes. They are not the CCSDS ``(255, 223)`` or ``(255, 239)`` codes, there is no ``\mathrm{GF}(2^8)`` table, and [`rs_decode`](@ref) is a hard-decision bounded-distance decoder, not a soft or list decoder.

# References
Reed and Solomon, J. SIAM 8(2):300–304, 1960; Lin and Costello, *Error Control Coding*, 2nd ed., ch. 7.
"""
function rs_code(m::Integer, k::Integer)
    F = GF2m(m)
    n = F.order - 1
    kk = Int(k)
    (1 <= kk <= n - 1) || throw(ArgumentError(
        "RS($n, k) over GF(2^$(F.m)) needs 1 <= k <= $(n - 1), got $kk",
    ))
    t = (n - kk) >> 1
    g = _rs_generator(F, n - kk)
    return RSCode(F, n, kk, t, g, 0)
end

@doc raw"""
    rs_shorten(C::RSCode, s::Integer) -> RSCode

Shorten `C` by `s` information symbols: the ``(n, k)`` code becomes the ``(n - s, k - s)`` code with the same generator and the same correction radius ``t``.

[`rs_encode`](@ref) prepends `s` zero symbols to the message, encodes the parent, and drops the `s` leading zeros of the systematic word. [`rs_decode`](@ref) re-inserts the zeros before decoding and refuses any error locator that falls in the shortened positions: those symbols were never transmitted, so such a locator marks an uncorrectable pattern. ``0 \le s \le k - 1`` so the effective dimension stays positive. Shortening does not reduce the minimum distance (Lin and Costello, 2nd ed., ch. 4).
"""
function rs_shorten(C::RSCode, s::Integer)
    ss = Int(s)
    (0 <= ss <= C.k - 1) || throw(ArgumentError(
        "shortening RS($(C.n), $(C.k)) needs 0 <= s <= $(C.k - 1), got $ss",
    ))
    return RSCode(C.F, C.n - ss, C.k - ss, C.t, C.g, C.s + ss)
end

@doc raw"""
    rs_encode(C::RSCode, msg::AbstractVector{<:Integer}) -> Vector{Int}

Systematic encoding of `k` information symbols into the Reed–Solomon codeword of length `n` of [`RSCode`](@ref).

Symbols are integers in ``0, \ldots, 2^m - 1``, the element representation of [`GF2m`](@ref). Entry 1 of the returned vector is the coefficient of ``x^{n-1}`` and entry `n` the coefficient of ``x^0``, the same order as [`bch157_encode`](@ref). The first `k` entries are `msg`; the last ``n - k`` entries are the remainder of ``m(x)\,x^{n-k}`` divided by the generator:

```math
c(x) = m(x)\,x^{n-k} + \bigl(m(x)\,x^{n-k} \bmod g(x)\bigr),
```

so ``c(\alpha^j) = 0`` for ``j = 1, \ldots, n - k``. A code from [`rs_shorten`](@ref) takes the effective ``k`` symbols and returns the effective ``n`` symbols.
"""
function rs_encode(C::RSCode, msg::AbstractVector{<:Integer})
    length(msg) == C.k || throw(ArgumentError(
        "RS($(C.n), $(C.k)) takes $(C.k) information symbols, got $(length(msg))",
    ))
    info = [_gf_elem(C.F, a) for a in msg]
    d = length(C.g) - 1
    shifted = vcat(zeros(Int, C.s), info, zeros(Int, d))
    parity = _gf_poly_mod(C.F, shifted, C.g)
    return vcat(info, parity)
end

@doc raw"""
    rs_decode(C::RSCode, rx::AbstractVector{<:Integer}) -> Vector{Int}

Bounded-distance decoding of one received word of [`RSCode`](@ref). Corrects every pattern of at most `t` symbol errors and returns the `k` systematic symbols.

Syndromes are ``S_j = r(\alpha^j)`` for ``j = 1, \ldots, n - k``. The error-locator polynomial ``\Lambda(x)`` comes from the Berlekamp–Massey algorithm, its roots from a Chien search over ``\alpha^1, \ldots, \alpha^n``, and the error values from the Forney formula

```math
e_J = \frac{\Omega(X_J^{-1})}{\Lambda'(X_J^{-1})},
\qquad
\Omega(x) \equiv S(x)\,\Lambda(x) \pmod{x^{\,n-k}},
```

with ``S(x) = S_1 + S_2 x + \cdots``, ``X_J`` the locator of the error at entry ``J`` (the coefficient of ``x^{\,n-J}``), and ``\Lambda'`` the formal derivative over ``\mathrm{GF}(2^m)``. Entry ``J`` is corrected by XOR with ``e_J``.

The failure contract is the one of [`bch157_decode`](@ref): there is no failure flag. If the locator degree exceeds `t`, the number of roots is not the locator degree, a root falls in the shortened positions of [`rs_shorten`](@ref), or the corrected word still has a nonzero syndrome, the pattern is not a correctable pattern of weight at most `t` and the function returns the `k` received systematic symbols unchanged. A pattern of more than `t` errors may instead be miscorrected to a different codeword; ``d = n - k + 1`` does not promise to correct it.

This is the hard-decision bounded-distance decoder. It is not an erasure decoder, not the Guruswami–Sudan list decoder, and not a soft algebraic decoder such as Kötter–Vardy. `rx` has length `n` and entries in ``0, \ldots, 2^m - 1``, in the same coefficient order as [`rs_encode`](@ref).

# References
Massey, IEEE Trans. Inf. Theory 15(1):122–127, 1969; Chien, IEEE Trans. Inf. Theory 10(4):357–363, 1964; Forney, IEEE Trans. Inf. Theory 11(4):549–557, 1965; Lin and Costello, 2nd ed., ch. 7.
"""
function rs_decode(C::RSCode, rx::AbstractVector{<:Integer})
    length(rx) == C.n || throw(ArgumentError(
        "RS($(C.n), $(C.k)) word length must be $(C.n), got $(length(rx))",
    ))
    r = [_gf_elem(C.F, a) for a in rx]
    F = C.F
    s = C.s
    N = C.n + s
    word = s == 0 ? copy(r) : vcat(zeros(Int, s), r)
    d = length(C.g) - 1
    S = _rs_syndromes(F, word, d)
    out = r[1:C.k]
    all(iszero, S) && return out
    Λ, L = _berlekamp_massey(F, S)
    L > C.t && return out
    locs = Vector{Int}(undef, L)
    nroots = 0
    x = 2
    for J in 1:N
        if _gf_eval_low(F, Λ, x) == 0
            nroots += 1
            nroots <= L && (locs[nroots] = J)
        end
        x = gf_mul(F, x, 2)
    end
    nroots == L || return out
    # A root in the shortened positions would correct a symbol that was never
    # transmitted; the pattern is not a correctable pattern of this code.
    for i in 1:L
        locs[i] > s || return out
    end
    Ω = _rs_key_equation(F, S, Λ)
    Λd = _gf_deriv_low(Λ)
    for i in 1:L
        J = locs[i]
        xinv = gf_pow(F, 2, mod(J, N))
        den = _gf_eval_low(F, Λd, xinv)
        den == 0 && return out
        word[J] ⊻= gf_div(F, _gf_eval_low(F, Ω, xinv), den)
    end
    all(iszero, _rs_syndromes(F, word, d)) || return out
    return word[(s + 1):(s + C.k)]
end
