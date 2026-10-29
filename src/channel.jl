module Channel

import ..Modulation
using ..Metrics: snr_db_to_linear

import FFTW: ifft

export awgn, rayleigh_channel, block_rayleigh, block_rayleigh_channel,
    jakes_fading, jakes_channel,
    rice_fading, rice_channel,
    nakagami_envelope, nakagami_fading, nakagami_channel,
    mrc, mrc_channel,
    apply_fir, fir_channel,
    CDLProfile, cdl_profile, CDLTaps, cdl_taps, apply_cdl, cdl_channel

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

Blocks are independent. There is no Doppler filter. The Clarke spectrum is [`jakes_fading`](@ref), not this function. The draws use `randn`, on the same global stream as `rayleigh_channel`.
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

function _rice_k(K::Real)
    Kf = float(K)
    (Kf >= 0 && !isnan(Kf)) || throw(ArgumentError(
        "K must be the linear Rice factor ≥ 0, not decibels, got $K",
    ))
    return Kf
end

@doc raw"""
    rice_fading(n::Integer, K::Real) -> Vector{ComplexF64}

`n` i.i.d. flat Rice coefficients with unit average power and line-of-sight phase 0,

```math
h = \sqrt{\frac{K}{K+1}} + \sqrt{\frac{1}{K+1}}\, g, \qquad g \sim \mathcal{CN}(0, 1).
```

`K` is the **linear** power ratio ``P_{\mathrm{LOS}} / P_{\mathrm{scatter}}``. It is not in decibels; pass ``10^{K_{\mathrm{dB}}/10}``. ``\mathbb{E}[|h|^2] = 1`` for every finite `K`.

- `K = 0` is ``\mathcal{CN}(0, 1)``, the same law and the same `randn` pattern as the coefficient draw inside [`rayleigh_channel`](@ref): all real parts, then all imaginary parts.
- `K = Inf` is the constant `1`. No scatter is drawn, and the global random stream is not advanced.
- `K < 0` or `NaN` is rejected. `n < 0` is rejected.

There is no Doppler filter and no Jakes spectrum. Draws use `randn`.
"""
function rice_fading(n::Integer, K::Real)
    n = Int(n)
    n >= 0 || throw(ArgumentError("n must be ≥ 0, got $n"))
    Kf = _rice_k(K)
    if isinf(Kf)
        return fill(ComplexF64(1), n)
    end
    los = sqrt(Kf / (Kf + 1))
    sc = sqrt(1 / (2 * (Kf + 1)))
    return los .+ sc .* (randn(n) .+ im .* randn(n))
end

@doc raw"""
    rice_channel(symbols, snr_db::Real, bits_per_symbol::Integer; K::Real) -> (rx, h)
    rice_channel(symbols, snr_db::Real, bits_per_symbol::Integer, perfect_csi::Bool; K::Real)

Flat Rice fading, independent from symbol to symbol, then AWGN,

```math
y = h x + n,
```

with `h` from [`rice_fading`](@ref). Real `symbols` are promoted to complex. ``E_s = 1`` is assumed and is not estimated. The noise follows the same ``E_b/N_0`` convention as [`awgn`](@ref), with ``E_b = 1/\texttt{bits\_per\_symbol}``. Because ``\mathbb{E}[|h|^2] = 1``, `snr_db` is the average ``E_b/N_0``. ``\mathrm{snr\_db} = +\infty`` disables the noise but still draws `h` (except when `K = Inf`, where `h` is identically 1 and no scatter is drawn).

`K` is a required keyword. It is the linear Rice factor documented on [`rice_fading`](@ref), not a decibel value.

# Returns
- Without `perfect_csi`: `(rx, h)`, the faded observation and the coefficients.
- `perfect_csi = true`: the zero-forced samples `rx ./ h`, and `h` is dropped. A coefficient that lands on 0 yields `NaN` in that entry.
- `perfect_csi = false`: `rx` only, with no phase compensation.

# Coherent BPSK
With perfect CSI the conditional bit error rate is ``Q(\sqrt{2 |h|^2 \bar\gamma_b})``. Averaging over this `h` is the Rice integral used by the tests (Craig's representation of `Q` against the Rician moment-generating function; Proakis and Salehi, 5th ed., §13.3). `K = 0` reduces to the Rayleigh formula on [`rayleigh_channel`](@ref). `K = Inf` is ordinary AWGN.

Zero-forcing scales the noise on symbol ``i`` to variance ``\sigma^2 / |h_i|^2``. One scalar variance is the wrong input to [`Modulation.demodulate_soft`](@ref). The coherent BPSK LLR is ``2\Re(h_i^* y_i) / \sigma^2``.
"""
function rice_channel(syms::AbstractVector, snr_db::Real, bits_per_symbol::Integer; K::Real)
    bits_per_symbol >= 1 || throw(ArgumentError("bits_per_symbol must be ≥ 1"))
    syms_c = ComplexF64.(syms)
    h = rice_fading(length(syms_c), K)
    return awgn(syms_c .* h, snr_db, bits_per_symbol), h
end

function rice_channel(
    syms::AbstractVector,
    snr_db::Real,
    bits_per_symbol::Integer,
    perfect_csi::Bool;
    K::Real,
)
    rx, h = rice_channel(syms, snr_db, bits_per_symbol; K=K)
    return perfect_csi ? rx ./ h : rx
end

@doc raw"""
    mrc(y, h) -> Vector{ComplexF64}

Maximum-ratio combining with ideal channel knowledge.

`y` and `h` are `L × N` matrices, one column per symbol and one row per branch:

```math
y_{\ell,n} = h_{\ell,n}\, x_n + \text{noise}.
```

The combined sample is

```math
z_n = \frac{\sum_{\ell=1}^{L} h_{\ell,n}^{*}\, y_{\ell,n}}{\sum_{\ell=1}^{L} |h_{\ell,n}|^2}.
```

When the branch noises are i.i.d. circular Gaussian of per-dimension variance ``\sigma^2`` and the model above holds, ``z_n = x_n + n'_n`` with

```math
\mathrm{Var}(\Re n'_n) = \mathrm{Var}(\Im n'_n) = \sigma^2 \Big/ \sum_{\ell} |h_{\ell,n}|^2.
```

A hard demodulator written for the original constellation can take `z` directly. [`Modulation.demodulate_soft`](@ref) cannot: it expects one noise variance for the whole vector. The coherent BPSK LLR, in that function's sign convention, is

```math
L_n = \frac{2}{\sigma^2}\,\Re\sum_{\ell} h_{\ell,n}^{*}\, y_{\ell,n},
```

which uses the numerator before the normalisation above.

`size(y)` must equal `size(h)`, and `L ≥ 1`. A column whose coefficients are all 0 yields `NaN`. Real matrices are promoted to complex. This is receive diversity with known gains. It is not Alamouti decoding, not a MIMO detector, and it does not estimate `h`.
"""
function mrc(y::AbstractMatrix{<:Number}, h::AbstractMatrix{<:Number})
    size(y) == size(h) || throw(ArgumentError(
        "y and h must have the same size, got $(size(y)) and $(size(h))",
    ))
    L, N = size(y)
    L >= 1 || throw(ArgumentError("MRC needs at least one branch, got $L"))
    z = Vector{ComplexF64}(undef, N)
    @inbounds for n in 1:N
        acc = 0.0 + 0.0im
        den = 0.0
        for ℓ in 1:L
            hl = ComplexF64(h[ℓ, n])
            acc += conj(hl) * ComplexF64(y[ℓ, n])
            den += abs2(hl)
        end
        z[n] = acc / den
    end
    return z
end

@doc raw"""
    mrc_channel(symbols, snr_db::Real, bits_per_symbol::Integer, L::Integer; K::Real=0) -> (z, h)

`L`-branch receive diversity, flat fading, ideal-CSI maximum-ratio combining.

Each branch observes the same unit-energy symbol through an independent coefficient from [`rice_fading`](@ref), then independent AWGN at the stated average ``E_b/N_0``. `K` is the linear Rice factor, default `0` (Rayleigh). It is not in decibels. ``\mathbb{E}[|h_{\ell}|^2] = 1`` on every branch, so `snr_db` is the average per-branch ``E_b/N_0``, not the SNR after combining. The symbol energy is not split across branches.

# Returns
- `z`: length-`N` vector from [`mrc`](@ref), so a noiseless call returns the symbols.
- `h`: the `L × N` coefficients that were applied and then used as CSI. There is no estimation error.

`L ≥ 1` and `bits_per_symbol ≥ 1`. ``\mathrm{snr\_db} = +\infty`` adds no noise. For coherent BPSK and i.i.d. Rayleigh branches (`K = 0`), the tests compare against

```math
P_b = \left(\frac{1-\mu}{2}\right)^L
\sum_{k=0}^{L-1}\binom{L-1+k}{k}\left(\frac{1+\mu}{2}\right)^k,
\qquad
\mu = \sqrt{\frac{\bar\gamma_b}{1+\bar\gamma_b}},
```

with ``\bar\gamma_b`` the per-branch average (Proakis and Salehi, 5th ed., §13.4). `L = 1` is the single-branch Rayleigh or Rice receiver, up to the normalisation in [`mrc`](@ref), which for one branch is `y ./ h`.

Branches are independent. There is no spatial correlation, no Jakes spectrum, and no transmit diversity.
"""
function mrc_channel(
    syms::AbstractVector{<:Number},
    snr_db::Real,
    bits_per_symbol::Integer,
    L::Integer;
    K::Real=0.0,
)
    bits_per_symbol >= 1 || throw(ArgumentError("bits_per_symbol must be ≥ 1"))
    L = Int(L)
    L >= 1 || throw(ArgumentError("L must be ≥ 1, got $L"))
    x = ComplexF64.(syms)
    N = length(x)
    h = Matrix{ComplexF64}(undef, L, N)
    y = Matrix{ComplexF64}(undef, L, N)
    for ℓ in 1:L
        hl = rice_fading(N, K)
        @inbounds for n in 1:N
            h[ℓ, n] = hl[n]
        end
        y[ℓ, :] = awgn(x .* hl, snr_db, bits_per_symbol)
    end
    return mrc(y, h), h
end

function _jakes_params(fd::Real, Ts::Real)
    f = float(fd)
    t = float(Ts)
    (isfinite(f) && f >= 0) || throw(ArgumentError(
        "fd must be a finite maximum Doppler frequency ≥ 0, in hertz, got $fd",
    ))
    (isfinite(t) && t > 0) || throw(ArgumentError(
        "Ts must be a finite sample interval > 0, in seconds, got $Ts",
    ))
    f * t < 0.5 || throw(ArgumentError(
        "fd*Ts must be < 1/2 so the Clarke spectrum does not alias, got fd*Ts=$(f * t)",
    ))
    return f, t
end

"""Power of the unit-area Clarke spectrum on `[f_lo, f_hi] ∩ [-fd, fd]`."""
function _clarke_interval_power(f_lo::Float64, f_hi::Float64, fd::Float64)
    a = max(f_lo, -fd)
    b = min(f_hi, fd)
    b <= a && return 0.0
    return (asin(b / fd) - asin(a / fd)) / π
end

"""
Clarke bin powers in FFT order, summing to 1.

Bin `k` (0-based) covers one DFT cell of width `1/(N Ts)`. For even `N` the
Nyquist cell is split across `±1/(2 Ts)`. `fd` is strictly positive;
`jakes_fading` handles `fd == 0` before calling this.
"""
function _clarke_bin_powers(N::Int, fd::Float64, Ts::Float64)
    P = zeros(N)
    Δf = 1 / (N * Ts)
    fs = 1 / Ts
    for k in 0:(N - 1)
        if iseven(N) && k == N ÷ 2
            P[k + 1] = _clarke_interval_power(fs / 2 - Δf / 2, fs / 2, fd) +
                       _clarke_interval_power(-fs / 2, -fs / 2 + Δf / 2, fd)
        elseif k <= (N - 1) ÷ 2
            P[k + 1] = _clarke_interval_power((k - 0.5) * Δf, (k + 0.5) * Δf, fd)
        else
            kk = k - N
            P[k + 1] = _clarke_interval_power((kk - 0.5) * Δf, (kk + 0.5) * Δf, fd)
        end
    end
    s = sum(P)
    s > 0 || error("Clarke bin powers summed to $s")
    s == 1 && return P
    return P .* inv(s)
end

@doc raw"""
    jakes_fading(n::Integer, fd::Real, Ts::Real) -> Vector{ComplexF64}

`n` samples of a unit-power, flat, time-correlated Rayleigh process at spacing `Ts` seconds, with maximum Doppler `fd` hertz.

The autocorrelation the construction targets is Clarke's,

```math
\mathbb{E}[h_t h_{t-m}^{*}] = J_0(2\pi f_d m T_s),
\qquad \mathbb{E}[|h|^2] = 1,
```

which is the inverse Fourier transform of the Jakes spectrum

```math
S(f) = \frac{1}{\pi f_d \sqrt{1 - (f/f_d)^2}}, \qquad |f| < f_d,
```

and ``S(f) = 0`` otherwise (Clarke, Bell System Technical Journal, 1968; the spectrum is the one named in Jakes, *Microwave Mobile Communications*, 1974). ``\mathbb{E}[h_t h_{t-m}] = 0``: the process is proper, so the in-phase and quadrature components are uncorrelated and each has variance ``1/2``.

# Implementation
This is **not** the classical Jakes sum of sinusoids. It is one period of a circularly symmetric complex Gaussian process synthesised on an `n`-point DFT:

1. FFT bin ``k`` is given power ``P_k`` equal to the integral of ``S`` over that bin. The antiderivative is ``\arcsin(f/f_d)/\pi``, so the band-edge singularity is integrable and no bin is sampled at ``f = \pm f_d``. For even `n` the Nyquist bin is the pair of half-cells at ``\pm 1/(2 T_s)``.
2. ``P_k`` is divided by ``\sum P`` so the discrete spectrum has unit power. The sum differs from 1 only by rounding in `asin`.
3. Bins with ``P_k > 0`` get independent ``\mathcal{CN}(0, P_k)`` coefficients. Zero-power bins stay zero and do not draw. The time series is `n * ifft` of that vector (`FFTW.ifft` divides by `n`).

The ensemble autocorrelation is exactly the inverse DFT of ``\{P_k\}``, including the periodic replicas of an `n`-point DFT. It approaches ``J_0(2\pi f_d m T_s)`` when the bins are fine, that is when `n * fd * Ts` is large. A short record, or a Doppler so small that only the DC bin is occupied, does not resolve ``J_0``. The sequence is exactly periodic with period `n`. It is not an IIR filter, not a streaming generator, and not a tapped-delay line.

`fd = 0` skips the FFT. One ``\mathcal{CN}(0, 1)`` draw is repeated on every sample, which is ``J_0(0) = 1``.

# Arguments
- `n`: number of samples. `n ≥ 0`. `n = 0` returns an empty vector and does not draw.
- `fd`: maximum Doppler in hertz. Finite, `≥ 0`.
- `Ts`: sample interval in seconds. Finite and `> 0`.

`fd * Ts` must lie in `[0, 1/2)`. A symbol rate of 10 kHz and a 50 Hz Doppler is `jakes_fading(n, 50, 1e-4)`.

# Random draws
For `fd = 0`: `randn()` for the real part, then `randn()` for the imaginary part. For `fd > 0`: the same pair, in increasing FFT-bin index, and only for bins with positive Clarke power. [`jakes_channel`](@ref) then passes the faded symbols through [`awgn`](@ref), which draws its own noise.

Each sample is exactly ``\mathcal{CN}(0, 1)`` in ensemble, for every admissible `fd`. Ideal-CSI coherent BPSK therefore has the same **average** bit error rate as symbol-wise [`rayleigh_channel`](@ref). The correlation changes the burst structure, not that average. This is not frequency-selective fading, not Rice, and not a 3GPP delay profile.
"""
function jakes_fading(n::Integer, fd::Real, Ts::Real)
    n = Int(n)
    n >= 0 || throw(ArgumentError("n must be ≥ 0, got $n"))
    f, t = _jakes_params(fd, Ts)
    n == 0 && return ComplexF64[]
    if f == 0
        c = (randn() + im * randn()) / sqrt(2)
        return fill(c, n)
    end
    P = _clarke_bin_powers(n, f, t)
    spec = zeros(ComplexF64, n)
    @inbounds for k in 1:n
        Pk = P[k]
        Pk > 0 || continue
        σ = sqrt(Pk / 2)
        spec[k] = σ * randn() + im * σ * randn()
    end
    return n .* ifft(spec)
end

@doc raw"""
    jakes_channel(symbols, snr_db::Real, bits_per_symbol::Integer, fd::Real, Ts::Real) -> (rx, h)
    jakes_channel(symbols, snr_db::Real, bits_per_symbol::Integer, fd::Real, Ts::Real, perfect_csi::Bool)

Flat time-correlated Rayleigh fading, then AWGN,

```math
y = h x + n,
```

with `h` from [`jakes_fading`](@ref). Real `symbols` are promoted to complex. ``E_s = 1`` is assumed and is not estimated. The noise follows the same ``E_b/N_0`` convention as [`awgn`](@ref), with ``E_b = 1/\texttt{bits\_per\_symbol}``. Because ``\mathbb{E}[|h|^2] = 1``, `snr_db` is the average ``E_b/N_0``. ``\mathrm{snr\_db} = +\infty`` adds a zero noise vector and still draws that vector through [`awgn`](@ref).

`fd` and `Ts` are the maximum Doppler in hertz and the sample interval in seconds, as on [`jakes_fading`](@ref). `fd * Ts` must be `< 1/2`.

# Returns
- Without `perfect_csi`: `(rx, h)`, the faded observation and the coefficients that were applied.
- `perfect_csi = true`: the zero-forced samples `rx ./ h`, and `h` is dropped. A coefficient that lands on 0 yields `NaN` in that entry.
- `perfect_csi = false`: `rx` only, with no phase compensation.

# Coherent use
Ideal-CSI coherent BPSK has conditional bit error rate ``Q(\sqrt{2 |h|^2 \bar\gamma_b})``. The expectation over the marginal of `h` is the Rayleigh formula on [`rayleigh_channel`](@ref), at every `fd`. The time correlation does not enter that average. Errors arrive in bursts of length on the order of ``1/(f_d T_s)`` samples, so a bit-i.i.d. confidence interval computed on one drop is miscalibrated.

Zero-forcing scales the noise on symbol ``i`` to variance ``\sigma^2 / |h_i|^2``. One scalar variance is the wrong input to [`Modulation.demodulate_soft`](@ref). The coherent BPSK LLR is ``2\Re(h_i^* y_i) / \sigma^2``.

This is one flat coefficient per symbol. It is not a frequency-selective Jakes tap model, not [`block_rayleigh`](@ref), and not [`CommToolkit.OFDM.TDLChannel`](@ref).
"""
function jakes_channel(
    syms::AbstractVector,
    snr_db::Real,
    bits_per_symbol::Integer,
    fd::Real,
    Ts::Real,
)
    bits_per_symbol >= 1 || throw(ArgumentError("bits_per_symbol must be ≥ 1"))
    syms_c = ComplexF64.(syms)
    h = jakes_fading(length(syms_c), fd, Ts)
    return awgn(syms_c .* h, snr_db, bits_per_symbol), h
end

function jakes_channel(
    syms::AbstractVector,
    snr_db::Real,
    bits_per_symbol::Integer,
    fd::Real,
    Ts::Real,
    perfect_csi::Bool,
)
    rx, h = jakes_channel(syms, snr_db, bits_per_symbol, fd, Ts)
    return perfect_csi ? rx ./ h : rx
end

function _nakagami_m(m::Real)
    mf = float(m)
    if isinf(mf) && mf > 0
        return mf
    end
    (isfinite(mf) && mf >= 0.5) || throw(ArgumentError(
        "m must be the Nakagami parameter ≥ 1/2 (or +Inf), got $m",
    ))
    return mf
end

"""Degrees of freedom `2m` when that value is an integer in `Float64`; otherwise `nothing`."""
function _nakagami_chi2_df(m::Float64)
    tw = 2 * m
    k = round(Int, tw)
    (k >= 1 && tw == k) || return nothing
    return k
end

"""Gamma(shape ≥ 1, scale = 1) by Marsaglia and Tsang, ACM TOMS 26(3), 2000."""
function _gamma_marsaglia_tsang(a::Float64)
    d = a - 1 / 3
    c = inv(sqrt(9 * d))
    for _ in 1:10_000
        x = randn()
        v = 1 + c * x
        while v <= 0
            x = randn()
            v = 1 + c * x
        end
        v = v * v * v
        u = rand()
        if u < 1 - 0.0331 * (x^2)^2 || log(u) < 0.5 * x^2 + d * (1 - v + log(v))
            return d * v
        end
    end
    error("Marsaglia–Tsang gamma sampler exceeded 10000 rejections")
end

"""Gamma(shape = `a` > 0, scale = 1). Shapes below 1 use the Marsaglia–Tsang reduction."""
function _gamma_scale1(a::Float64)
    a >= 1 && return _gamma_marsaglia_tsang(a)
    return _gamma_marsaglia_tsang(a + 1) * rand()^(1 / a)
end

function _nakagami_chi2_envelope(n::Int, k::Int)
    z = randn(k * n)
    R = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        s = 0.0
        base = (i - 1) * k
        for j in 1:k
            s += z[base + j]^2
        end
        R[i] = sqrt(s / k)
    end
    return R
end

@doc raw"""
    nakagami_envelope(n::Integer, m::Real) -> Vector{Float64}

`n` i.i.d. Nakagami-m amplitudes with ``\mathbb{E}[R^2] = 1`` and ``m \ge 1/2``,

```math
p(r) = \frac{2 m^m}{\Gamma(m)} r^{2m-1} \exp(-m r^2), \qquad r \ge 0.
```

``R^2`` is Gamma with shape `m` and scale `1/m`. There is no Doppler filter: samples are independent. This is the envelope only. The complex coefficient is [`nakagami_fading`](@ref).

# Generator
- `m = 1`: ``R = \sqrt{(Z_1^2 + Z_2^2)/2}`` with `randn(n)` then `randn(n)`. The modulus of [`nakagami_fading`](@ref) at `m = 1` consumes the same stream.
- `m = Inf`: the constant `1`. No random draw. This is a deterministic amplitude, not a deterministic complex coefficient.
- `2m` an integer value in `Float64` (so `m = 1/2, 1, 3/2, 2, …`, and the `m = 1` case above is special-cased out of this branch): ``R^2`` is the mean of `2m` squared standard normals, drawn as one `randn(2m * n)` vector, sample by sample. `m = 1/2` is ``R = |Z|``.
- Every other finite `m`: Marsaglia and Tsang's rejection sampler for the gamma power (ACM Trans. Math. Softw. 26(3):363–372, 2000), then the positive square root. Shapes below 1 use their ``\Gamma(m) = \Gamma(m+1)\,U^{1/m}`` reduction. The number of `randn` / `rand` calls depends on the accept/reject path.

`n ≥ 0`. `n = 0` returns an empty vector and does not draw. `m < 1/2`, `NaN`, and a negative infinity are rejected.
"""
function nakagami_envelope(n::Integer, m::Real)
    n = Int(n)
    n >= 0 || throw(ArgumentError("n must be ≥ 0, got $n"))
    mf = _nakagami_m(m)
    n == 0 && return Float64[]
    isinf(mf) && return ones(n)
    if mf == 1
        return hypot.(randn(n), randn(n)) ./ sqrt(2)
    end
    k = _nakagami_chi2_df(mf)
    if k !== nothing
        return _nakagami_chi2_envelope(n, k)
    end
    scale = inv(mf)
    R = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        R[i] = sqrt(_gamma_scale1(mf) * scale)
    end
    return R
end

@doc raw"""
    nakagami_fading(n::Integer, m::Real) -> Vector{ComplexF64}

`n` i.i.d. complex Nakagami-m coefficients with ``\mathbb{E}[|h|^2] = 1``.

For `m = 1` the law is ``\mathcal{CN}(0, 1)``, and the draws match the coefficient inside [`rayleigh_channel`](@ref) and [`rice_fading`](@ref) at `K = 0`: `randn(n)` for the real part, then `randn(n)` for the imaginary part. No separate phase draw.

For every other `m`, the amplitude is [`nakagami_envelope`](@ref) and the phase is an independent ``\mathrm{Unif}[0, 2\pi)`` from `2π * rand(n)`, taken after the amplitude stream. The phase is **not** tied to a Gaussian in-phase/quadrature pair when `m ≠ 1`. `m = Inf` is therefore a uniform phase on the unit circle, not the constant `1` returned by [`rice_fading`](@ref) at `K = Inf`.

`m` is the Nakagami parameter, ``\ge 1/2`` or `+Inf`, not in decibels. `n ≥ 0`. `n = 0` returns an empty vector and does not draw. Samples are independent across time. This is not [`jakes_fading`](@ref).
"""
function nakagami_fading(n::Integer, m::Real)
    n = Int(n)
    n >= 0 || throw(ArgumentError("n must be ≥ 0, got $n"))
    mf = _nakagami_m(m)
    n == 0 && return ComplexF64[]
    if mf == 1
        return (randn(n) .+ im .* randn(n)) ./ sqrt(2)
    end
    return nakagami_envelope(n, mf) .* cis.(2π .* rand(n))
end

@doc raw"""
    nakagami_channel(symbols, snr_db::Real, bits_per_symbol::Integer; m::Real) -> (rx, h)
    nakagami_channel(symbols, snr_db::Real, bits_per_symbol::Integer, perfect_csi::Bool; m::Real)

i.i.d. flat Nakagami-m fading, then AWGN,

```math
y = h x + n,
```

with `h` from [`nakagami_fading`](@ref). Real `symbols` are promoted to complex. ``E_s = 1`` is assumed and is not estimated. The noise follows the same ``E_b/N_0`` convention as [`awgn`](@ref), with ``E_b = 1/\texttt{bits\_per\_symbol}``. Because ``\mathbb{E}[|h|^2] = 1``, `snr_db` is the average ``E_b/N_0``. ``\mathrm{snr\_db} = +\infty`` adds a zero noise vector and still draws that vector through [`awgn`](@ref), except that the fading draw follows [`nakagami_fading`](@ref) (`m = Inf` draws only the phase).

`m` is a required keyword. It is the Nakagami parameter on [`nakagami_fading`](@ref), not a decibel value.

# Returns
- Without `perfect_csi`: `(rx, h)`, the faded observation and the coefficients.
- `perfect_csi = true`: the zero-forced samples `rx ./ h`, and `h` is dropped. A coefficient that lands on 0 yields `NaN` in that entry.
- `perfect_csi = false`: `rx` only, with no phase compensation.

# Coherent BPSK
With perfect CSI the conditional bit error rate is ``Q(\sqrt{2 |h|^2 \bar\gamma_b})``. Averaging over this envelope is the Gamma-MGF integral used by the tests,

```math
P_b = \frac{1}{\pi}\int_0^{\pi/2}
\left(\frac{m\sin^2\theta}{m\sin^2\theta + \bar\gamma_b}\right)^{m} d\theta.
```

`m = 1` reduces to the Rayleigh formula on [`rayleigh_channel`](@ref). For positive integer `m` the same integral equals `m`-branch i.i.d. Rayleigh MRC with per-branch average ``\bar\gamma_b / m`` (the Gamma(`m`) power is the sum of `m` unit-mean exponential powers, scaled back to average 1). `m = Inf` has ``|h| = 1``, so perfect CSI removes only a per-symbol phase and the bit error rate is ordinary AWGN. That phase is still random: `m = Inf` is not [`rice_channel`](@ref) at `K = Inf`.

Zero-forcing scales the noise on symbol ``i`` to variance ``\sigma^2 / |h_i|^2``. One scalar variance is the wrong input to [`Modulation.demodulate_soft`](@ref). The coherent BPSK LLR is ``2\Re(h_i^* y_i) / \sigma^2``.

Coefficients are independent from symbol to symbol. There is no Jakes spectrum and no spatial correlation.
"""
function nakagami_channel(syms::AbstractVector, snr_db::Real, bits_per_symbol::Integer; m::Real)
    bits_per_symbol >= 1 || throw(ArgumentError("bits_per_symbol must be ≥ 1"))
    syms_c = ComplexF64.(syms)
    h = nakagami_fading(length(syms_c), m)
    return awgn(syms_c .* h, snr_db, bits_per_symbol), h
end

function nakagami_channel(
    syms::AbstractVector,
    snr_db::Real,
    bits_per_symbol::Integer,
    perfect_csi::Bool;
    m::Real,
)
    rx, h = nakagami_channel(syms, snr_db, bits_per_symbol; m=m)
    return perfect_csi ? rx ./ h : rx
end

include("cdl.jl")

end # module Channel
