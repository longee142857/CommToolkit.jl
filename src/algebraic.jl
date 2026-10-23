# GF(2^m) for m = 2, 3, 4, the [7, 4] Hamming syndrome decoder, and the
# binary primitive BCH code of length 15, dimension 7, designed distance 5.
# Not a Reed–Solomon library and not a general-n BCH constructor.

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

``m=4`` is the field under [`bch157_encode`](@ref). There is no ``\mathrm{GF}(2^8)`` table and no Conway polynomial search. Build the field with `GF2m(m)`; the fields of the struct are the log and antilog tables and are not a second constructor.

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
