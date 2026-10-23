"""
    Scrambling

Additive (synchronous) LFSR scrambler. A block is XORed with a linear-feedback sequence that restarts from a seed on every call. The same call with the same polynomial and seed is the descrambler.

This is not a self-synchronizing (multiplicative) scrambler, not a Gold code, and not a 3GPP NR scrambler. It is not wired into the MCS / HARQ loop. Named degrees are classical connection polynomials; the bits are the raw shift-register output, not an inverted ITU-T O.150 test pattern.
"""
module Scrambling

export prbs_polynomial, lfsr_sequence, scramble, descramble

@doc raw"""
    prbs_polynomial(degree::Integer) -> Int

Connection polynomial of a classical binary m-sequence, as an integer bitmask. Bit ``k`` is the coefficient of ``x^k``, so the value is odd and its highest set bit is ``x^{\mathrm{degree}}``.

| degree | polynomial |
| --- | --- |
| 7 | ``x^7 + x^6 + 1`` |
| 9 | ``x^9 + x^5 + 1`` |
| 11 | ``x^{11} + x^9 + 1`` |
| 15 | ``x^{15} + x^{14} + 1`` |
| 20 | ``x^{20} + x^3 + 1`` |
| 23 | ``x^{23} + x^{18} + 1`` |
| 29 | ``x^{29} + x^{27} + 1`` |
| 31 | ``x^{31} + x^{28} + 1`` |

Degree 7 is the usual non-standard trinomial (not the IEEE 802.11 scrambler ``x^7 + x^4 + 1``). Degrees 9 and above are connection polynomials from ITU-T Recommendation O.150. Degree 20 is ``x^{20} + x^3 + 1``, not the other O.150 polynomial ``x^{20} + x^{17} + 1``. O.150 inverts the length-``2^{15}-1``, ``2^{23}-1``, ``2^{29}-1``, and ``2^{31}-1`` patterns, and the ``x^{20} + x^{17} + 1`` variant stuffs extra ones. Those modifications are not applied. Any other polynomial is passed straight to [`lfsr_sequence`](@ref) as a bitmask; this function only names the eight above.
"""
function prbs_polynomial(degree::Integer)
    d = Int(degree)
    d == 7 && return (1 << 7) | (1 << 6) | 1
    d == 9 && return (1 << 9) | (1 << 5) | 1
    d == 11 && return (1 << 11) | (1 << 9) | 1
    d == 15 && return (1 << 15) | (1 << 14) | 1
    d == 20 && return (1 << 20) | (1 << 3) | 1
    d == 23 && return (1 << 23) | (1 << 18) | 1
    d == 29 && return (1 << 29) | (1 << 27) | 1
    d == 31 && return (1 << 31) | (1 << 28) | 1
    throw(ArgumentError(
        "no named PRBS polynomial for degree $d; supported degrees are 7, 9, 11, 15, 20, 23, 29, and 31",
    ))
end

function _lfsr_params(poly::Integer)
    poly > 0 || throw(ArgumentError("polynomial bitmask must be positive, got $poly"))
    n = ndigits(poly, base=2) - 1
    n >= 1 || throw(ArgumentError("LFSR degree must be ≥ 1, got $n"))
    n <= 62 || throw(ArgumentError("LFSR degree must be ≤ 62, got $n"))
    isodd(poly) || throw(ArgumentError(
        "connection polynomial must include the x^0 term (bit 0 set), got $poly",
    ))
    tapmask = UInt64(poly) & ((UInt64(1) << n) - UInt64(1))
    return n, tapmask
end

function _seed_state(seed::Union{Integer,Nothing}, degree::Int)
    maxs = (UInt64(1) << degree) - UInt64(1)
    if seed === nothing
        return maxs
    end
    seed >= 1 || throw(ArgumentError("seed must be in 1:$(Int(maxs)), got $seed"))
    # maxs fits in Int for degree ≤ 62, but seed may be a wider Integer.
    if seed > maxs
        throw(ArgumentError("seed must be in 1:$(Int(maxs)), got $seed"))
    end
    return UInt64(seed)
end

@doc raw"""
    lfsr_sequence(poly::Integer, nbits::Integer; seed=nothing) -> Vector{Int}

First `nbits` bits of a Fibonacci LFSR, each entry `0` or `1`.

# Polynomial
`poly` is a bitmask: bit ``k`` is the coefficient of ``x^k``. The degree ``n`` is the highest set bit, and bit 0 must be set (the ``x^0`` term). ``1 \le n \le 62``. [`prbs_polynomial`](@ref) supplies the classical trinomials. Primitivity is not checked.

# Register
The state is an integer whose bits ``0, \ldots, n-1`` are the shift register. Bit 0 is the bit emitted on this step. If the polynomial is ``x^n + a_{n-1} x^{n-1} + \cdots + a_0`` and the output is ``s_0, s_1, \ldots``, the state at step ``t`` holds ``s_t`` in bit 0 through ``s_{t+n-1}`` in bit ``n-1``, and the feedback is

```math
s_{t+n} = a_{n-1} s_{t+n-1} + \cdots + a_0 s_t \pmod{2}.
```

The new bit is written into bit ``n-1`` after the register shifts toward bit 0. `seed = nothing` loads every stage with one, so the first ``n`` outputs are `1` and the seed value is ``2^n - 1``. Any other seed is that integer state and must lie in ``1:(2^n - 1)`` (the all-zero register is rejected). Each call starts again at `seed`; nothing is retained between calls.

# Notes
`nbits = 0` returns an empty `Vector{Int}` after the polynomial is checked. The sequence is exactly periodic with period ``2^n - 1`` only when `poly` is primitive. Named polynomials from [`prbs_polynomial`](@ref) are primitive. The bits are not inverted and are not stuffed: this is not an ITU-T O.150 test-pattern generator.
"""
function lfsr_sequence(poly::Integer, nbits::Integer; seed::Union{Integer,Nothing}=nothing)
    nbits >= 0 || throw(ArgumentError("nbits must be ≥ 0, got $nbits"))
    degree, tapmask = _lfsr_params(poly)
    state = _seed_state(seed, degree)
    nb = Int(nbits)
    out = Vector{Int}(undef, nb)
    @inbounds for i in 1:nb
        out[i] = Int(state & 0x1)
        feedback = count_ones(state & tapmask) & 1
        state = (state >> 1) | (UInt64(feedback) << (degree - 1))
    end
    return out
end

@doc raw"""
    scramble(bits::AbstractVector{<:Integer}, poly::Integer; seed=nothing) -> Vector{Int}

Additive scrambler. Each bit of `bits` is XORed with the corresponding bit of [`lfsr_sequence`](@ref)`(poly, length(bits); seed)`.

`bits` must contain only `0` and `1` (`Bool` is accepted). The result is a `Vector{Int}` of the same length and does not share storage with `bits`. `seed = nothing` is the all-ones register. The same `poly` and `seed` make [`descramble`](@ref) the inverse, including when `poly` is not primitive, because the sequence is restarted rather than free-running.
"""
function scramble(bits::AbstractVector{<:Integer}, poly::Integer; seed::Union{Integer,Nothing}=nothing)
    seq = lfsr_sequence(poly, length(bits); seed=seed)
    out = Vector{Int}(undef, length(bits))
    i = 0
    for b in bits
        i += 1
        b == 0 || b == 1 || throw(ArgumentError("bits must be 0 or 1, got $b"))
        @inbounds out[i] = Int(b) ⊻ seq[i]
    end
    return out
end

@doc raw"""
    descramble(bits::AbstractVector{<:Integer}, poly::Integer; seed=nothing) -> Vector{Int}

Inverse of [`scramble`](@ref) for the same `poly` and `seed`.

Additive scrambling is an involution when the LFSR restarts: XORing the same sequence a second time restores the block. `descramble` is that second call. A different seed, or a different polynomial, does not invert the first call.
"""
function descramble(bits::AbstractVector{<:Integer}, poly::Integer; seed::Union{Integer,Nothing}=nothing)
    return scramble(bits, poly; seed=seed)
end

end # module Scrambling
