module Modulation

export ModulationScheme, BPSK, QPSK, QAM, modulate, demodulate, demodulate_soft,
    constellation, bits_per_symbol

@doc raw"""
    abstract type ModulationScheme end

Marker for a digital constellation. Concrete schemes are [`BPSK`](@ref), [`QPSK`](@ref), and square [`QAM`](@ref).

Symbols are normalised to unit average energy ``E_s = 1``. With ``k`` bits per symbol, ``E_b = 1/k``.
"""
abstract type ModulationScheme end

@doc raw"""
    BPSK()

Binary phase-shift keying. The bit map is ``0 \mapsto +1``, ``1 \mapsto -1``, so ``E_s = E_b = 1``.

In AWGN the bit error probability is ``Q(\sqrt{2 E_b/N_0})``
(Proakis and Salehi, *Digital Communications*, 5th ed., McGraw-Hill, 2008, §5.2).
"""
struct BPSK <: ModulationScheme end

# Gray QPSK, unit average energy. Index 2*b0 + b1 + 1:
# 00 → (+1+j)/√2, 01 → (+1−j)/√2, 10 → (−1+j)/√2, 11 → (−1−j)/√2.
const _QPSK = ComplexF64[1 + 1im, 1 - 1im, -1 + 1im, -1 - 1im] ./ sqrt(2)

@doc raw"""
    QPSK()

Gray-labelled quadrature phase-shift keying with unit average symbol energy.

Bit pairs ``(b_0, b_1)``, ``b_0`` first, map to

| bits | symbol |
|------|--------|
| 00 | ``(+1+j)/\sqrt{2}`` |
| 01 | ``(+1-j)/\sqrt{2}`` |
| 10 | ``(-1+j)/\sqrt{2}`` |
| 11 | ``(-1-j)/\sqrt{2}`` |

Adjacent quadrants differ by one bit. Each component is an independent BPSK decision of amplitude ``1/\sqrt{2}``, so the bit error probability in AWGN equals the BPSK formula ``Q(\sqrt{2 E_b/N_0})`` with ``E_b = E_s/2`` (Proakis and Salehi, 5th ed., §5.2).
"""
struct QPSK <: ModulationScheme end

@doc raw"""
    QAM(m::Integer)

Square Gray-labelled ``M``-QAM, ``M \in \{4, 16, 64, 256, \ldots\}``.

Each axis is a Gray-labelled ``\sqrt{M}``-PAM alphabet ``\{-(\sqrt{M}-1), \ldots, \sqrt{M}-1\}`` scaled by ``\sqrt{3 / (2(M-1))}``, which gives ``E_s = 1``. The two most-significant bits of each ``\log_2 M``-bit block select the in-phase Gray label; the remaining bits select the quadrature label (MSB first).

The exact average bit error rate of this labelling in AWGN is the square-QAM formula of Cho and Yoon, IEEE Trans. Commun., 50(7):1074–1080, 2002. A common high-SNR approximation is

```math
P_b \approx \frac{4}{\log_2 M}\left(1 - \frac{1}{\sqrt{M}}\right)
Q\left(\sqrt{\frac{3 \log_2 M}{M-1}\,\frac{E_b}{N_0}}\right)
```

(Proakis and Salehi, 5th ed., §5.2). ``M = 4`` is Gray 4-QAM, distinct in type from [`QPSK`](@ref) but with the same geometry.
"""
struct QAM <: ModulationScheme
    m::Int
    function QAM(m::Integer)
        M = Int(m)
        if M < 4 || !ispow2(M) || isqrt(M)^2 != M
            throw(ArgumentError("M must be a square power of two ≥ 4 (4, 16, 64, …), got $m"))
        end
        return new(M)
    end
end

@doc raw"""
    bits_per_symbol(scheme::ModulationScheme) -> Int

Number of bits carried by one symbol: 1 for BPSK, 2 for QPSK, and ``\log_2 M`` for `QAM(M)`.
"""
bits_per_symbol(::BPSK) = 1
bits_per_symbol(::QPSK) = 2
bits_per_symbol(scheme::QAM) = Int(log2(scheme.m))

function _qam_constellation(M::Int)
    n = isqrt(M)
    scale = sqrt(2 * (M - 1) / 3)
    pam = [-(n - 1):2:(n - 1);] ./ scale
    gray = [xor(i, i >> 1) for i in 0:(n - 1)]
    constell = Vector{ComplexF64}(undef, M)
    for i in 1:n, j in 1:n
        constell[gray[i] * n + gray[j] + 1] = pam[i] + pam[j] * im
    end
    return constell
end

@doc raw"""
    constellation(scheme::ModulationScheme) -> Vector{ComplexF64}

Ideal constellation points in the modulator's index order, with average energy 1.

For `QAM`, index ``i`` (1-based) is the symbol for the ``\log_2 M`` bits of ``i - 1`` written MSB first.
"""
constellation(::BPSK) = ComplexF64[1, -1]
constellation(::QPSK) = copy(_QPSK)
constellation(scheme::QAM) = _qam_constellation(scheme.m)

function _check_bits(bits)
    all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("bits must be 0 or 1"))
    return nothing
end

@doc raw"""
    modulate(bits::AbstractVector{<:Integer}, scheme::ModulationScheme)

Map bits to unit-energy symbols.

# Arguments
- `bits`: entries in ``\{0, 1\}``. Length must be a multiple of [`bits_per_symbol`](@ref)`(scheme)`.
- `scheme`: [`BPSK`](@ref), [`QPSK`](@ref), or [`QAM`](@ref).

# Returns
- `BPSK`: real vector, ``0 \mapsto +1``, ``1 \mapsto -1``.
- `QPSK` and `QAM`: `Vector{ComplexF64}` using the Gray maps documented on those types.

# References
Proakis and Salehi, *Digital Communications*, 5th ed., §5.2 (constellation normalisation and Gray labelling).
"""
function modulate end

function modulate(bits::AbstractVector{<:Integer}, ::BPSK)
    _check_bits(bits)
    return 1 .- 2 .* bits
end

function modulate(bits::AbstractVector{<:Integer}, ::QPSK)
    _check_bits(bits)
    length(bits) % 2 == 0 || throw(ArgumentError("bit length must be even for QPSK"))
    n_sym = length(bits) ÷ 2
    syms = Vector{ComplexF64}(undef, n_sym)
    @inbounds for i in 1:n_sym
        b0 = bits[2i - 1]
        b1 = bits[2i]
        syms[i] = _QPSK[2 * b0 + b1 + 1]
    end
    return syms
end

function modulate(bits::AbstractVector{<:Integer}, scheme::QAM)
    _check_bits(bits)
    M = scheme.m
    k = bits_per_symbol(scheme)
    length(bits) % k == 0 ||
        throw(ArgumentError("bit length must be a multiple of log2(M) = $k"))
    constell = _qam_constellation(M)
    n_sym = length(bits) ÷ k
    syms = Vector{ComplexF64}(undef, n_sym)
    @inbounds for i in 1:n_sym
        idx = 0
        base = (i - 1) * k
        for j in 1:k
            idx = (idx << 1) | bits[base + j]
        end
        syms[i] = constell[idx + 1]
    end
    return syms
end

@doc raw"""
    demodulate(symbols, scheme::ModulationScheme) -> Vector{Int}

Hard decisions, inverse to [`modulate`](@ref) in the noiseless case.

# Arguments
- `symbols`: real or complex samples for BPSK (only the real part is tested). Complex samples for QPSK and QAM.
- `scheme`: the modulator that produced the symbols.

# Returns
Bits in ``\{0, 1\}``, MSB first for QPSK and QAM. QPSK decides each component by its sign. QAM uses minimum-distance decoding, which is maximum likelihood in isotropic AWGN (Proakis and Salehi, 5th ed., §4.3 and §5.2).
"""
function demodulate end

function demodulate(syms::AbstractVector{<:Real}, ::BPSK)
    return Int.(syms .< 0)
end

function demodulate(syms::AbstractVector{<:Complex}, ::BPSK)
    return Int.(real.(syms) .< 0)
end

function demodulate(syms::AbstractVector{<:Complex}, ::QPSK)
    n_sym = length(syms)
    bits = Vector{Int}(undef, 2 * n_sym)
    @inbounds for i in 1:n_sym
        bits[2i - 1] = Int(real(syms[i]) < 0)
        bits[2i] = Int(imag(syms[i]) < 0)
    end
    return bits
end

function demodulate(syms::AbstractVector{<:Complex}, scheme::QAM)
    M = scheme.m
    k = bits_per_symbol(scheme)
    constell = _qam_constellation(M)
    bits = Vector{Int}(undef, length(syms) * k)
    @inbounds for i in eachindex(syms)
        _, min_idx = findmin(abs2.(syms[i] .- constell))
        label = min_idx - 1
        base = (i - 1) * k
        for j in 1:k
            bits[base + j] = (label >> (k - j)) & 1
        end
    end
    return bits
end

@doc raw"""
    demodulate_soft(symbols::AbstractVector{<:Real}, ::BPSK, noise_var::Real) -> Vector{Float64}

BPSK log-likelihood ratio for the bit map ``0 \mapsto +1``, ``1 \mapsto -1``,

```math
L(y) = \log\frac{p(y \mid b=0)}{p(y \mid b=1)} = \frac{2y}{\sigma^2},
```

where `noise_var` is the real noise variance ``\sigma^2 = N_0/2`` (Proakis and Salehi, 5th ed., §4.2). A positive LLR favours bit 0. `noise_var` must be positive.
"""
function demodulate_soft(syms::AbstractVector{<:Real}, ::BPSK, noise_var::Real)
    noise_var > 0 || throw(ArgumentError("noise_var must be positive, got $noise_var"))
    return (2 .* syms) ./ float(noise_var)
end

end # module Modulation
