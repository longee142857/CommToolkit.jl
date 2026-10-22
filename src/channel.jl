module Channel

import ..Modulation
using ..Metrics: snr_db_to_linear

export awgn, rayleigh_channel, block_rayleigh, block_rayleigh_channel,
    apply_fir, fir_channel

@doc raw"""
    awgn(symbols, snr_db::Real; energy_per_bit::Real=1.0)
    awgn(symbols, snr_db::Real, bits_per_symbol::Integer)

Add real or circularly symmetric complex AWGN.

# Arguments
- `symbols`: a real vector is given one real noise component; a complex vector is given independent real and imaginary components.
- `snr_db`: ``E_b/N_0`` in decibels. `Inf` adds no noise.
- `energy_per_bit`: ``E_b`` in the same units as the symbol energy. The keyword form uses this value directly.
- `bits_per_symbol`: positional form. Assumes unit average symbol energy ``E_s = 1``, so ``E_b = 1/\texttt{bits\_per\_symbol}``.

# Returns
A noisy copy of `symbols`. The noise on each real dimension has variance

```math
\sigma^2 = \frac{N_0}{2} = \frac{E_b}{2\,E_b/N_0}.
```

For a complex symbol the total noise power is ``N_0``. This is the baseband AWGN model used with the BPSK / QPSK / ``M``-QAM error formulas in Proakis and Salehi, *Digital Communications*, 5th ed., §4.2 and §5.2.

# Notes
The positional method does not estimate ``E_s`` from `symbols`. Callers that use [`Modulation.modulate`](@ref) already have ``E_s = 1``.
"""
function awgn(syms::AbstractVector{T}, snr_db::Real; energy_per_bit::Real=1.0) where {T}
    energy_per_bit > 0 || throw(ArgumentError("energy_per_bit must be positive"))
    ρ = snr_db_to_linear(snr_db)
    noise_var = energy_per_bit / (2 * ρ)
    noise_std = sqrt(noise_var)
    n = length(syms)
    if T <: Complex
        noise = noise_std .* (randn(n) .+ im .* randn(n))
    else
        noise = noise_std .* randn(n)
    end
    return syms .+ noise
end

function awgn(syms::AbstractVector, snr_db::Real, bits_per_symbol::Integer)
    bits_per_symbol >= 1 || throw(ArgumentError("bits_per_symbol must be ≥ 1"))
    return awgn(syms, snr_db; energy_per_bit=1.0 / bits_per_symbol)
end

@doc raw"""
    rayleigh_channel(symbols, snr_db::Real, bits_per_symbol::Integer) -> (rx, h)
    rayleigh_channel(symbols, snr_db::Real, bits_per_symbol::Integer, perfect_csi::Bool)

Flat Rayleigh fading, independent from symbol to symbol,

```math
y = h x + n, \qquad h \sim \mathcal{CN}(0, 1),
```

so ``\mathbb{E}[|h|^2] = 1``. Real `symbols` are promoted to complex. The noise follows the same ``E_b/N_0`` convention as [`awgn`](@ref), with ``E_s = 1`` and ``E_b = 1/\texttt{bits\_per\_symbol}``. ``\mathrm{snr\_db} = +\infty`` disables the noise.

# Returns
- The two-argument (plus SNR) method returns `(rx, h)`: the faded observation and the channel coefficients.
- `perfect_csi = true` returns the zero-forced samples `rx ./ h` and drops `h`.
- `perfect_csi = false` returns `rx` only, with no phase compensation. Coherent BPSK demodulation of that output is not the perfect-CSI receiver.

# References
With perfect CSI, coherent BPSK has the closed-form bit error rate

```math
P_b = \frac12\left(1 - \sqrt{\frac{\bar\gamma_b}{1 + \bar\gamma_b}}\right),
\qquad \bar\gamma_b = E_b/N_0,
```

because the instantaneous SNR is ``|h|^2 \bar\gamma_b``. See Goldsmith, *Wireless Communications*, Cambridge University Press, 2005, §6.3, and Proakis and Salehi, *Digital Communications*, 5th ed., §13.3.
"""
function rayleigh_channel(syms::AbstractVector, snr_db::Real, bits_per_symbol::Integer)
    n = length(syms)
    syms_c = ComplexF64.(syms)
    h = (randn(n) .+ im .* randn(n)) ./ sqrt(2)
    faded = syms_c .* h
    return awgn(faded, snr_db, bits_per_symbol), h
end

function rayleigh_channel(syms::AbstractVector, snr_db::Real, bits_per_symbol::Integer, perfect_csi::Bool)
    rx, h = rayleigh_channel(syms, snr_db, bits_per_symbol)
    return perfect_csi ? rx ./ h : rx
end

@doc raw"""
    block_rayleigh(n::Integer, block_length::Integer) -> Vector{ComplexF64}

Block-constant Rayleigh coefficients of length `n`. Samples inside one block are equal, and distinct blocks are independent:

```math
h \sim \mathcal{CN}(0, 1), \qquad \mathbb{E}[|h|^2] = 1.
```

Block ``b = 0, 1, \ldots`` occupies indices ``b \cdot \texttt{block\_length} + 1`` through the end of that block. The last block is shorter when `n` is not a multiple of `block_length`. `block_length == 1` is symbol-by-symbol fading, the same marginal law as [`rayleigh_channel`](@ref).

Blocks are independent. There is no Doppler filter and no Jakes spectrum. The draws use `randn`, on the same global stream as `rayleigh_channel`.
"""
function block_rayleigh(n::Integer, block_length::Integer)
    n = Int(n)
    block_length = Int(block_length)
    n >= 0 || throw(ArgumentError("n must be ≥ 0, got $n"))
    block_length >= 1 || throw(ArgumentError("block_length must be ≥ 1, got $block_length"))
    h = Vector{ComplexF64}(undef, n)
    i = 1
    while i <= n
        g = (randn() + im * randn()) / sqrt(2)
        stop = min(i + block_length - 1, n)
        @inbounds for j in i:stop
            h[j] = g
        end
        i = stop + 1
    end
    return h
end

@doc raw"""
    block_rayleigh_channel(symbols, snr_db::Real, bits_per_symbol::Integer, block_length::Integer) -> (rx, h)

Block Rayleigh fading, ``y = h x + n``, with `h` from [`block_rayleigh`](@ref).

# Arguments
- `symbols`: real or complex. Real inputs are promoted to complex. ``E_s = 1`` is assumed and is not estimated.
- `snr_db`: ``E_b/N_0`` in decibels, the same positional convention as [`rayleigh_channel`](@ref). ``E_b = 1/\texttt{bits\_per\_symbol}``. `+Inf` adds no noise.
- `bits_per_symbol`: passed to [`awgn`](@ref). Must be ≥ 1.
- `block_length`: number of consecutive symbols that share one coefficient. This is not a perfect-CSI flag.

# Returns
The faded observation `rx` and the per-symbol coefficients `h` (constant on each block).

# Coherent use
Hard BPSK decisions with perfect CSI are `demodulate(real.(rx ./ h), BPSK())`. The division scales the noise on symbol ``i`` to variance ``\sigma^2/|h_i|^2``, with ``\sigma^2 = N_0/2`` from [`awgn`](@ref), so one scalar variance is the wrong input to [`Modulation.demodulate_soft`](@ref). The coherent BPSK LLR, in the sign convention of that function, is

```math
L_i = \frac{2\,\Re(h_i^* y_i)}{\sigma^2}.
```

Symbol-wise Rayleigh fading (`block_length == 1`) gives an interleaver nothing to spread. A block longer than one symbol is the case in which a block interleaver changes the burst pattern seen by a short convolutional code.
"""
function block_rayleigh_channel(
    syms::AbstractVector,
    snr_db::Real,
    bits_per_symbol::Integer,
    block_length::Integer,
)
    bits_per_symbol >= 1 || throw(ArgumentError("bits_per_symbol must be ≥ 1"))
    syms_c = ComplexF64.(syms)
    h = block_rayleigh(length(syms_c), block_length)
    return awgn(syms_c .* h, snr_db, bits_per_symbol), h
end

@doc raw"""
    apply_fir(signal, taps) -> Vector

Full linear convolution of `signal` with the FIR taps `taps`,

```math
y[k] = \sum_{j=1}^{L} h_j\, x[k-j+1],
```

with ``x[i] = 0`` outside `1:length(signal)`. The result has length ``n + L - 1``.

This is the same model as the tall Toeplitz matrix in [`CommToolkit.Equalization.zf_equalize`](@ref): `zf_equalize(apply_fir(s, h), h)` recovers `s` when `h` is not identically zero, up to linear-algebra roundoff. Nothing is truncated and the convolution is not circular. The first and last ``L - 1`` samples are the filter warmup and the tail.

`signal` and `taps` must be non-empty. The element type is `promote_type(eltype(signal), eltype(taps), Float64)`. Taps are used as given: they are not normalised, not drawn at random, and not a fading process.
"""
function apply_fir(signal::AbstractVector{<:Number}, taps::AbstractVector{<:Number})
    n = length(signal)
    m = length(taps)
    n >= 1 || throw(ArgumentError("signal must be non-empty"))
    m >= 1 || throw(ArgumentError("taps must be non-empty"))
    T = promote_type(eltype(signal), eltype(taps), Float64)
    out = zeros(T, n + m - 1)
    @inbounds for i in 1:n
        xi = convert(T, signal[i])
        for j in 1:m
            out[i + j - 1] += xi * convert(T, taps[j])
        end
    end
    return out
end

function _entry_noise_var(::Type{T}, snr_db::Real, energy_per_bit::Real) where {T}
    per_dim = energy_per_bit / (2 * snr_db_to_linear(snr_db))
    return T <: Complex ? 2 * per_dim : per_dim
end

@doc raw"""
    fir_channel(signal, taps, snr_db::Real; energy_per_bit::Real=1.0) -> (rx, noise_var)
    fir_channel(signal, taps, snr_db::Real, bits_per_symbol::Integer) -> (rx, noise_var)

[`apply_fir`](@ref), then AWGN on the convolved samples.

# Arguments
- `signal`, `taps`: passed to [`apply_fir`](@ref). Taps are deterministic. This is not [`CommToolkit.OFDM.TDLChannel`](@ref), not a 3GPP delay profile, and not a Jakes spectrum.
- `snr_db`: ``E_b/N_0`` in decibels, the same convention as [`awgn`](@ref). It is **not** the signal-to-noise ratio of the samples after the FIR, and the tap energy is not folded into ``E_b``. `+Inf` adds no noise and does not draw from the global random stream.
- `energy_per_bit`: ``E_b`` in the keyword method. Must be positive. It is not estimated from `signal`.
- `bits_per_symbol`: positional method. Assumes ``E_s = 1``, so ``E_b = 1/\texttt{bits\_per\_symbol}``, and must be ≥ 1. Callers that used [`Modulation.modulate`](@ref) already have ``E_s = 1``.

# Returns
- `rx`: the convolved sequence, plus noise when `snr_db` is finite. Same length as [`apply_fir`](@ref).
- `noise_var`: variance of one entry of that additive noise, in the sense [`CommToolkit.Equalization.mmse_equalize`](@ref) uses. A real convolution has variance ``N_0/2``; a complex convolution (complex `signal` or complex `taps`) has variance ``N_0 = E_b / (E_b/N_0)``. `+Inf` returns `noise_var = 0`.

Pass `noise_var` to `mmse_equalize` only when the equalizer sees this same `rx`. A later matched filter colors the noise; the scalar above is then not the variance at the symbol sampler.
"""
function fir_channel(
    signal::AbstractVector{<:Number},
    taps::AbstractVector{<:Number},
    snr_db::Real;
    energy_per_bit::Real=1.0,
)
    energy_per_bit > 0 || throw(ArgumentError("energy_per_bit must be positive"))
    y = apply_fir(signal, taps)
    noise_var = _entry_noise_var(eltype(y), snr_db, energy_per_bit)
    if isinf(snr_db) && snr_db > 0
        return y, noise_var
    end
    return awgn(y, snr_db; energy_per_bit=energy_per_bit), noise_var
end

function fir_channel(
    signal::AbstractVector{<:Number},
    taps::AbstractVector{<:Number},
    snr_db::Real,
    bits_per_symbol::Integer,
)
    bits_per_symbol >= 1 || throw(ArgumentError(
        "bits_per_symbol must be ≥ 1, got $bits_per_symbol",
    ))
    return fir_channel(signal, taps, snr_db; energy_per_bit=1.0 / bits_per_symbol)
end

end # module Channel
