"""
    Sync

Teaching-grade symbol timing and carrier synchronisation. Not a production
synchroniser, and not a 3GPP NR, LTE, or O-RAN receiver. There is no PSS,
SSS, DMRS, lock detector, or AGC.

Single-carrier timing is a Gardner detector on an even number of samples per
symbol, plus a linear interpolator. Carrier recovery is a second-order Costas
loop for [`BPSK`](@ref) and [`QPSK`](@ref) at one sample per symbol. OFDM
frame detection is a Schmidl–Cox repeated-half preamble. It sits beside
[`OFDM.ofdm_sync`](@ref) and does not change that function: `ofdm_sync` is
still the van de Beek cyclic-prefix correlator, unambiguous only for a
carrier offset in ``(-1/2, 1/2)`` subcarrier spacings, with no integer-CFO
estimate and no shift of the buffer by the timing index.
"""
module Sync

using FFTW
import ..Modulation: BPSK, QPSK, PSK
import ..Pulse: matched_filter, downsample_symbols
import ..OFDM
import ..Visualization: pulse_shape

export fractional_delay,
    gardner_ted, gardner_loop, gardner_open_loop, single_carrier_timing,
    costas_detector, costas_loop,
    schmidl_preamble, schmidl_cox, integer_cfo_estimate, schmidl_sync

# Same Park–Miller constants as the OFDM pilots, different seed. Preamble
# signs do not touch the global RNG.
function _lcg_unit(seed::Int)
    r = (seed * 1103515245 + 12345) & 0x7fffffff
    return Float64(r) / 2147483648.0
end

function _as_float_vector(samples::AbstractVector{<:Number})
    T = eltype(samples) <: Complex ? ComplexF64 : Float64
    return T.(samples)
end

function _check_finite_sample(y::Number, what::AbstractString)
    if !(isfinite(real(y)) && isfinite(imag(y)))
        throw(ArgumentError("$what must be finite"))
    end
    return nothing
end

@doc raw"""
    fractional_delay(samples, delay::Real) -> Vector

Delay `samples` by `delay` samples with a two-point linear interpolator.

Sample `i` of the result (1-based) is the input at time `(i - 1) + delay` on the clock where input index 1 is time 0. An integer `delay` is a slice and keeps the element type. A non-integer delay promotes the result to `Float64` or `ComplexF64` and drops every output whose right-hand neighbour would fall past the end of `samples`, so the result is shorter than the input.

`delay` must be finite and ≥ 0. This is not a bandlimited fractional-delay filter and not a Farrow structure. It exists so [`gardner_open_loop`](@ref) and [`gardner_loop`](@ref) can move a sampling phase that is not an integer number of samples.
"""
function fractional_delay(samples::AbstractVector{<:Number}, delay::Real)
    δ = float(delay)
    isfinite(δ) || throw(ArgumentError("delay must be finite, got $delay"))
    δ >= 0 || throw(ArgumentError("delay must be ≥ 0, got $delay"))
    n = length(samples)
    n >= 1 || throw(ArgumentError("samples must be non-empty"))
    nout = floor(Int, n - δ)
    nout >= 1 || throw(ArgumentError(
        "delay $delay leaves no samples in a length-$n vector",
    ))
    if δ == floor(δ)
        shift = Int(δ)
        return Vector{eltype(samples)}(samples[(1 + shift):(nout + shift)])
    end
    out = Vector{promote_type(eltype(samples), Float64)}(undef, nout)
    @inbounds for i in 1:nout
        p = (i - 1) + δ          # 0-based input time
        i0 = floor(Int, p)       # 0-based left index
        μ = p - i0
        out[i] = (1 - μ) * samples[i0 + 1] + μ * samples[i0 + 2]
    end
    return out
end

"""
    _interp(samples, index::Int, mu::Float64)

Value at 1-based position `index + mu`, with `mu` in ``[0, 1)``. `mu == 0` does not read `index + 1`.
"""
function _interp(samples::AbstractVector, index::Int, mu::Float64)
    mu == 0 && return samples[index]
    return (1 - mu) * samples[index] + mu * samples[index + 1]
end

@doc raw"""
    gardner_ted(samples; sps, phase=1) -> Vector{Float64}

Gardner timing-error detector on a uniform grid (Gardner, IEEE Trans. Commun., 34(5):423–429, 1986).

# Arguments
- `samples`: matched-filter output, or any pulse-shaped waveform. One sample stream, real or complex.
- `sps`: samples per symbol. Must be even and ≥ 2, so the midpoint between two symbol strobes lies on a sample. This is the same `sps` as [`pulse_shape`](@ref) / [`downsample_symbols`](@ref). It is not a symbol period in seconds.
- `phase`: 1-based index of the first symbol strobe. Symbol strobe ``k`` (``k ≥ 1``) is index `phase + (k - 1) * sps`. `phase` must lie in `1:sps`. The first error uses strobe 2, because strobe 1 has no previous symbol.

# Error
For each complete triple the midpoint is `sps ÷ 2` samples before the late strobe (halfway between the early and late strobes):

```math
e = \Re\!\left(y_{\mathrm{mid}}^{*}\,(y_{\mathrm{late}} - y_{\mathrm{early}})\right).
```

That is ``y_I(\mathrm{mid})\,(y_I(\mathrm{late}) - y_I(\mathrm{early}))`` plus the same product on the quadrature component.

# Sign
On a symmetric Nyquist pulse (a root-raised-cosine matched to itself), the average of `e` is an odd function of the timing error about the autocorrelation peak. A strobe that is **late** — the kept sample sits after the peak, toward the next symbol — produces a **positive** average. A strobe that is early produces a negative average. The slope through that zero is positive, so a loop that moves the strobe later when its fractional interval `mu` grows must apply `mu ← mu - alpha * e`.

There is a second zero about half a symbol away, with the opposite slope. It is not a stable lock point of that update. [`gardner_open_loop`](@ref) keeps the positive-slope zero.

# What this is not
The detector does not resolve an integer number of symbols of delay; only the sampling phase inside one symbol period moves the average. A run of identical symbols makes `late - early ≈ 0`, so the average carries no timing information. A residual carrier that is large compared with the symbol rate biases the average. No interpolator runs inside this function: a fractional offset has to be removed with [`fractional_delay`](@ref) first, or absorbed by [`gardner_loop`](@ref). Odd `sps` is rejected rather than interpolated. This is not a Mueller–Müller detector and not a production timing-recovery loop.
"""
function gardner_ted(samples::AbstractVector{<:Number}; sps::Integer, phase::Integer=1)
    sps = Int(sps)
    phase = Int(phase)
    (sps >= 2 && iseven(sps)) || throw(ArgumentError(
        "sps must be even and ≥ 2 so the midpoint lies on a sample, got $sps",
    ))
    1 <= phase <= sps || throw(ArgumentError(
        "phase must lie in 1:sps = 1:$sps, got $phase",
    ))
    n = length(samples)
    first = phase + sps
    first <= n || throw(ArgumentError(
        "need at least $first samples for one Gardner error at phase=$phase, sps=$sps; got $n",
    ))
    half = sps ÷ 2
    nerr = (n - first) ÷ sps + 1
    errs = Vector{Float64}(undef, nerr)
    @inbounds for k in 1:nerr
        late_i = first + (k - 1) * sps
        mid = samples[late_i - half]
        early = samples[late_i - sps]
        late = samples[late_i]
        _check_finite_sample(mid, "samples")
        _check_finite_sample(early, "samples")
        _check_finite_sample(late, "samples")
        errs[k] = real(conj(mid) * (late - early))
    end
    return errs
end

@doc raw"""
    gardner_loop(samples; sps, alpha=0.15, mu0=0.0, basepoint=sps+1, n_symbols=nothing)

First-order Gardner loop with a linear interpolator.

# Clock
`basepoint` is the 1-based index of the sample at or just before the first **late** strobe. The strobe sits at position `basepoint + mu` on that same clock (`mu` in ``[0, 1)``). The early strobe is `sps` samples earlier and the midpoint is `sps ÷ 2` samples earlier, both with the same `mu`. The default `basepoint = sps + 1` is simply the second symbol of a `phase = 1` grid; it is not an estimate of the autocorrelation peak. Pass the peak index when the bulk delay of the pulse is known.

The error is the same formula as [`gardner_ted`](@ref), evaluated on the interpolated triple. The update is

```math
\mu \leftarrow \mu - \alpha\, e,
```

then ``\mu`` is wrapped into ``[0, 1)`` and `basepoint` advances by `sps` plus the number of samples that wrap carried. Because a late strobe makes `e > 0`, the update pulls `mu` backward. `alpha` is samples per error-unit. The error scales like the square of the waveform amplitude, so a gain that locks unit-energy root-raised-cosine BPSK does not transfer unchanged to a different pulse or to QPSK (the QPSK S-curve of this detector is flatter). The default `0.15` is a teaching value for unit-energy BPSK on `rrc_filter`, not a loop-bandwidth design.

# Returns
A named tuple `(symbols, error, mu, basepoint)`. Entry `k` is one symbol: `symbols[k]` is the interpolated late strobe, `error[k]` is the Gardner error, and the sampling instant is `basepoint[k] + mu[k]`. `mu[k]` is the interval **used** for that symbol, before the update.

# Limits
There is no symbol-clock frequency state. A constant sampling-phase error is what the loop is for. A clock that is fast or slow by a large fraction of a sample per symbol walks `mu` through its wrap and is not acquired cleanly. If one update fails to advance the strobe (`alpha * e` on the order of a symbol or more), the function throws rather than run backwards. Linear interpolation biases the lock point by a fraction of a sample; it is not a bandlimited interpolator. No lock detector. Not production timing recovery.

`n_symbols` stops after that many symbols. `nothing` runs until the next strobe would leave the buffer.
"""
function gardner_loop(samples::AbstractVector{<:Number};
                      sps::Integer,
                      alpha::Real=0.15,
                      mu0::Real=0.0,
                      basepoint::Integer=sps + 1,
                      n_symbols::Union{Nothing,Integer}=nothing)
    sps = Int(sps)
    (sps >= 2 && iseven(sps)) || throw(ArgumentError(
        "sps must be even and ≥ 2, got $sps",
    ))
    α = float(alpha)
    (isfinite(α) && α > 0) || throw(ArgumentError("alpha must be finite and > 0, got $alpha"))
    mu = float(mu0)
    (isfinite(mu) && 0 <= mu < 1) || throw(ArgumentError(
        "mu0 must lie in [0, 1), got $mu0",
    ))
    if n_symbols !== nothing
        n_symbols = Int(n_symbols)
        n_symbols >= 1 || throw(ArgumentError("n_symbols must be ≥ 1, got $n_symbols"))
    end
    base = Int(basepoint)
    base >= sps + 1 || throw(ArgumentError(
        "basepoint must be ≥ sps + 1 so the previous symbol exists, got $base",
    ))
    x = _as_float_vector(samples)
    n = length(x)
    half = sps ÷ 2
    T = eltype(x)
    symbols = T[]
    errors = Float64[]
    mus = Float64[]
    bases = Int[]
    while n_symbols === nothing || length(symbols) < n_symbols
        last_needed = max(base + (mu == 0 ? 0 : 1),
                          base - half + (mu == 0 ? 0 : 1))
        if base - sps < 1 || last_needed > n
            break
        end
        late = _interp(x, base, mu)
        early = _interp(x, base - sps, mu)
        mid = _interp(x, base - half, mu)
        _check_finite_sample(late, "samples")
        e = real(conj(mid) * (late - early))
        push!(symbols, late)
        push!(errors, e)
        push!(mus, mu)
        push!(bases, base)

        mu_next = mu - α * e
        shift = floor(Int, mu_next)
        mu_next -= shift
        base_next = base + sps + shift
        base_next > base || throw(ArgumentError(
            "timing update did not advance the strobe; reduce alpha",
        ))
        mu = mu_next
        base = base_next
    end
    isempty(symbols) && throw(ArgumentError(
        "buffer does not contain a full Gardner strobe at basepoint=$base, sps=$sps",
    ))
    return (symbols=symbols, error=errors, mu=mus, basepoint=bases)
end

function _positive_slope_zero(grid::Vector{Float64}, means::Vector{Float64}, period::Float64)
    n = length(grid)
    best = nothing
    best_slope = 0.0
    for i in 1:n
        j = i == n ? 1 : i + 1
        a = means[i]
        b = means[j]
        span = i == n ? (grid[1] + period - grid[end]) : (grid[j] - grid[i])
        rise = b - a
        # Strict rise: an all-zero curve (no symbol transitions) is not a lock.
        if a <= 0 <= b && rise > 0 && span > 0
            frac = -a / rise
            z = mod(grid[i] + frac * span, period)
            slope = rise / span
            if slope > best_slope
                best_slope = slope
                best = z
            end
        end
    end
    return best, best_slope
end

@doc raw"""
    gardner_open_loop(samples; sps, n_frac=8, phase=1) -> NamedTuple

Grid search for the positive-slope zero of [`gardner_ted`](@ref).

`n_frac` is the number of trial delays per sample, so the grid on one symbol period has `sps * n_frac` points: `0, 1/n_frac, …, sps - 1/n_frac`. At delay `d` the waveform is passed through [`fractional_delay`](@ref) and the detector is run at `phase` (default 1, strobes at indices `1, 1+sps, …` of the **delayed** waveform). The returned `delay` lies in ``[0, \texttt{sps})`` and is the linearly interpolated zero at which the average error changes from negative to positive.

After `fractional_delay(samples, delay)`, [`downsample_symbols`](@ref) with `offset = phase - 1` keeps that sampling phase. `delay` is defined only modulo `sps`. A value just below `sps` and a value of 0 are the same phase; the non-negative representative is what [`fractional_delay`](@ref) can apply. The integer symbol delay is not estimated. See [`single_carrier_timing`](@ref), which fills it in from `length(taps)` when the pulse is known.

`mean_error` and `grid` are the curve that was searched. The average includes the start and the end of the burst, where [`pulse_shape`](@ref) has truncated the pulse, so the burst needs to be long next to the filter span.

Throws if the curve has no positive-slope zero, or if the steepest such slope is below `0.01` times the average sample power ``\sum |r|^2 / N``. A burst with no symbol transitions sits near that floor: `late - early` vanishes in the steady state, and the truncation at the ends of [`pulse_shape`](@ref) leaves a slope on the order of ``3 \times 10^{-4}`` per sample for constant unit-energy BPSK. A mixed unit-energy root-raised-cosine QPSK burst is near ``0.1`` per sample, and alternating BPSK is steeper. The factor `0.01` is a teaching guard against locking onto truncation noise. It is not a detector threshold derived from a false-alarm rate.

The search covers one symbol period and one stable zero. It does not follow a drifting clock. Linear interpolation of both the waveform and the zero crossing leaves a fraction of a sample of bias. Not production timing recovery.
"""
function gardner_open_loop(samples::AbstractVector{<:Number};
                           sps::Integer, n_frac::Integer=8, phase::Integer=1)
    sps = Int(sps)
    n_frac = Int(n_frac)
    phase = Int(phase)
    (sps >= 2 && iseven(sps)) || throw(ArgumentError(
        "sps must be even and ≥ 2, got $sps",
    ))
    n_frac >= 1 || throw(ArgumentError("n_frac must be ≥ 1, got $n_frac"))
    ngrid = sps * n_frac
    grid = Vector{Float64}(undef, ngrid)
    means = Vector{Float64}(undef, ngrid)
    step = 1 / n_frac
    for i in 1:ngrid
        d = (i - 1) * step
        grid[i] = d
        y = d == 0 ? samples : fractional_delay(samples, d)
        e = gardner_ted(y; sps=sps, phase=phase)
        means[i] = sum(e) / length(e)
    end
    delay, slope = _positive_slope_zero(grid, means, float(sps))
    power = sum(abs2, samples) / length(samples)
    if delay === nothing || slope < 0.01 * power
        throw(ArgumentError(
            "Gardner open loop found no usable positive-slope zero; the burst may have too few transitions",
        ))
    end
    return (delay=delay, mean_error=means, grid=grid)
end

@doc raw"""
    single_carrier_timing(symbols, taps, sps; sample_delay=0, n_frac=8) -> NamedTuple

Glue [`pulse_shape`](@ref), [`matched_filter`](@ref), [`gardner_open_loop`](@ref), [`fractional_delay`](@ref), and [`downsample_symbols`](@ref).

The chain is: shape `symbols` with `taps` at `sps` samples per symbol, matched-filter with `conj.(reverse(taps))`, optionally delay the matched output by `sample_delay` samples (a known extra offset, so a caller can watch the estimate move), search the Gardner positive-slope zero, fractional-delay by that estimate, and keep every `sps`-th sample starting at index 1.

# Returns
`(symbols, delay, symbol_index)`.

- `symbols`: one sample per symbol at the corrected phase. Not sliced to the input length.
- `delay`: the open-loop delay in samples, in ``[0, \texttt{sps})``.
- `symbol_index`: 1-based index in the returned vector where **input symbol 1** lands, using the known bulk delay of the pulse rather than anything Gardner measured:

```math
\texttt{symbol\_index}
  = \operatorname{round}\!\left(
      \frac{\texttt{length(taps)} - \texttt{sample\_delay} - \texttt{delay} - 1}{\texttt{sps}}
    \right) + 1.
```

Input symbol ``k`` then sits at `symbol_index + k - 1`, when that index is inside the vector. [`pulse_shape`](@ref) discards the filter tail, so symbols within about half a filter span of the **end** of the burst are not the full matched-filter response and will not match. Gardner still does not estimate this integer delay; the formula uses `length(taps)`.

`sps` must be even and ≥ 2. `taps` are real, as [`pulse_shape`](@ref) requires. `sample_delay` must be finite and ≥ 0.

Not a receiver. No carrier loop, no equaliser, no noise model, no tracking of a symbol-clock frequency.
"""
function single_carrier_timing(symbols::AbstractVector{<:Number},
                               taps::AbstractVector{<:Real},
                               sps::Integer;
                               sample_delay::Real=0.0,
                               n_frac::Integer=8)
    sps = Int(sps)
    (sps >= 2 && iseven(sps)) || throw(ArgumentError(
        "sps must be even and ≥ 2, got $sps",
    ))
    δ = float(sample_delay)
    (isfinite(δ) && δ >= 0) || throw(ArgumentError(
        "sample_delay must be finite and ≥ 0, got $sample_delay",
    ))
    isempty(taps) && throw(ArgumentError("taps must be non-empty"))
    isempty(symbols) && throw(ArgumentError("symbols must be non-empty"))
    wave = matched_filter(pulse_shape(symbols, taps, sps), taps)
    if δ != 0
        wave = fractional_delay(wave, δ)
    end
    est = gardner_open_loop(wave; sps=sps, n_frac=n_frac, phase=1)
    aligned = est.delay == 0 ? wave : fractional_delay(wave, est.delay)
    sampled = downsample_symbols(aligned, sps; offset=0)
    peak = float(length(taps)) - δ - est.delay
    symbol_index = round(Int, (peak - 1) / sps) + 1
    return (symbols=sampled, delay=est.delay, symbol_index=symbol_index)
end

@doc raw"""
    costas_detector(sample, scheme) -> Float64

One phase-detector output of a Costas loop.

`scheme` is [`BPSK`](@ref) or [`QPSK`](@ref). No other scheme is accepted. In particular [`PSK`](@ref)`(4)` sits on the axes, a ``\pi/4`` rotation away from [`QPSK`](@ref), and is rejected; the QPSK detector below assumes the ``\pi/4`` constellation.

Let ``y = I + jQ`` be one symbol-rate sample, already downconverted by the current phase estimate.

- BPSK, bits on the real axis (``0 \mapsto +1``, ``1 \mapsto -1``):

```math
e = \operatorname{sign}(I)\, Q,
```

with ``\operatorname{sign}(0)`` taken positive. For a symbol ``\pm e^{j\theta}`` and ``|\theta| < \pi/2``, ``e = \sin\theta``.

- QPSK, the four points ``(\pm 1 \pm j)/\sqrt{2}``:

```math
e = \operatorname{sign}(I)\, Q - \operatorname{sign}(Q)\, I.
```

For a constellation point times ``e^{j\theta}`` and ``|\theta| < \pi/4``, ``e = \sqrt{2}\,\sin\theta``. The same expression is ``\sqrt{2}`` times the imaginary part of ``d^{*} y`` where ``d`` is the unit-energy hard decision ``(\operatorname{sign} I + j\operatorname{sign} Q)/\sqrt{2}``.

The sign is the sign of a small residual carrier phase. A loop that multiplies the next sample by ``e^{-j\hat\theta}`` must move ``\hat\theta`` in the direction of `e`.

BPSK is blind to a rotation of ``\pi`` (every bit flips, and `e` is still 0 at the rotated point). QPSK is blind to a rotation of ``\pi/2``. Those are stable false locks once ``|\theta|`` leaves the sector above. The detector is evaluated at one sample per symbol. It is not a carrier loop by itself, and it is not defined for OQPSK, ``\pi/4``-DQPSK, or QAM.
"""
function costas_detector(sample::Number, ::BPSK)
    _check_finite_sample(sample, "sample")
    i = real(sample)
    return (i >= 0 ? 1.0 : -1.0) * float(imag(sample))
end

function costas_detector(sample::Number, ::QPSK)
    _check_finite_sample(sample, "sample")
    i = float(real(sample))
    q = float(imag(sample))
    si = i >= 0 ? 1.0 : -1.0
    sq = q >= 0 ? 1.0 : -1.0
    return si * q - sq * i
end

function costas_detector(sample::Number, scheme)
    _check_finite_sample(sample, "sample")
    hint = scheme isa PSK ? " PSK(M) lies on a different grid from QPSK()." : ""
    throw(ArgumentError(
        "Costas detector is defined for BPSK and QPSK only, got $(typeof(scheme)).$hint",
    ))
end

@doc raw"""
    costas_loop(samples, scheme; alpha=0.06, beta=0.0015, phase0=0.0, freq0=0.0) -> NamedTuple

Second-order Costas loop at one sample per symbol.

`scheme` is [`BPSK`](@ref) or [`QPSK`](@ref), passed to [`costas_detector`](@ref). `samples` are symbol-rate complex baseband samples **after** timing recovery. A fractional timing error is not estimated here.

The state is a phase in radians and a frequency in radians per symbol. On sample ``y``:

```math
\begin{aligned}
y_k &= y\, e^{-j\hat\theta}, \\
e_k &= \text{Costas detector}(y_k), \\
\hat\omega &\leftarrow \hat\omega + \beta\, e_k, \\
\hat\theta &\leftarrow \hat\theta + \hat\omega + \alpha\, e_k.
\end{aligned}
```

``\hat\theta`` is then wrapped to ``(-\pi, \pi]``. ``\beta = 0`` is a first-order loop. ``\beta < 0`` is rejected.

# Returns
`(symbols, error, phase, frequency)`, one entry per input sample.

- `symbols[k]` is ``y_k``, derotated by the phase **before** the update.
- `phase[k]` is that same phase (radians removed from sample ``k``).
- `frequency[k]` is ``\hat\omega`` before the update on sample ``k``.
- `error[k]` is the detector output.

A positive residual phase produces a positive error, so `phase` moves up when the input was rotated by a positive angle ``e^{j\theta}``, ``\theta > 0``. A positive frequency offset (radians per symbol, the same sign as ``e^{j\omega (k-1)}``) drives `frequency` upward.

# Capture
The BPSK detector's linear region in this sign convention is ``|\theta| < \pi/2``. The QPSK region is ``|\theta| < \pi/4``. Outside those sectors the loop can lock to the opposite bit polarity (BPSK, period ``\pi``) or to a neighbouring quadrant (QPSK, period ``\pi/2``) without `error` growing. The default gains ``\alpha = 0.06``, ``\beta = 0.0015`` are teaching values that settle a frequency offset on the order of ``10^{-2}`` radians per symbol on unit-energy symbols. They are not a design for a stated loop bandwidth, and a offset of a large fraction of a radian per symbol is outside their acquisition range. There is no lock detector, no AGC, and no frequency sweep. Noise bandwidth is not optimised. Not a production carrier recovery, and not an OFDM equaliser.
"""
function costas_loop(samples::AbstractVector{<:Number}, scheme;
                     alpha::Real=0.06,
                     beta::Real=0.0015,
                     phase0::Real=0.0,
                     freq0::Real=0.0)
    α = float(alpha)
    β = float(beta)
    (isfinite(α) && α > 0) || throw(ArgumentError("alpha must be finite and > 0, got $alpha"))
    (isfinite(β) && β >= 0) || throw(ArgumentError("beta must be finite and ≥ 0, got $beta"))
    ϕ0 = float(phase0)
    ω0 = float(freq0)
    (isfinite(ϕ0) && isfinite(ω0)) || throw(ArgumentError("phase0 and freq0 must be finite"))
    isempty(samples) && throw(ArgumentError("samples must be non-empty"))
    phase = mod(ϕ0 + π, 2π) - π
    freq = ω0
    n = length(samples)
    symbols = Vector{ComplexF64}(undef, n)
    errors = Vector{Float64}(undef, n)
    phases = Vector{Float64}(undef, n)
    freqs = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        y = ComplexF64(samples[k])
        yrot = y * exp(-im * phase)
        e = costas_detector(yrot, scheme)
        symbols[k] = yrot
        errors[k] = e
        phases[k] = phase
        freqs[k] = freq
        freq += β * e
        phase += freq + α * e
        phase = mod(phase + π, 2π) - π
    end
    return (symbols=symbols, error=errors, phase=phases, frequency=freqs)
end

@doc raw"""
    schmidl_preamble(n_fft::Integer; n_cp::Integer=0) -> NamedTuple

One Schmidl–Cox training symbol (Schmidl and Cox, IEEE Trans. Commun., 45(12):1613–1621, 1997).

`n_fft` must be even and ≥ 4. In the 0-based FFT bin index ``k``, odd bins are zero and even bins (including DC) carry a deterministic ``\pm 1`` from a Park–Miller LCG seeded by `900 + j` (`j` is the 1-based Julia index, so the occupied bins are the odd Julia indices). An IFFT of that grid is two identical halves of length `n_fft ÷ 2`. Julia's `ifft` divides by `n_fft`, the same normalisation as [`OFDM.generate_ofdm_symbol`](@ref).

`n_cp` copies the last `n_cp` IFFT samples to the front. `0 ≤ n_cp ≤ n_fft`.

Returns `(time, freq, n_cp)`. `time` includes the prefix. `freq` is the length-`n_fft` grid and does not include a prefix. The signs do not use the global RNG.

[`OFDM.ofdm_tx`](@ref) does **not** prepend this symbol. [`OFDM.ofdm_sync`](@ref) does not look for it. The two estimators are separate: the cyclic-prefix correlator remains the van de Beek estimator documented on [`OFDM.cfo_estimate`](@ref).
"""
function schmidl_preamble(n_fft::Integer; n_cp::Integer=0)
    N = Int(n_fft)
    (N >= 4 && iseven(N)) || throw(ArgumentError(
        "n_fft must be even and ≥ 4, got $n_fft",
    ))
    n_cp = Int(n_cp)
    0 <= n_cp <= N || throw(ArgumentError(
        "n_cp must lie in 0:n_fft = 0:$N, got $n_cp",
    ))
    freq = zeros(ComplexF64, N)
    @inbounds for j in 1:2:N
        freq[j] = _lcg_unit(900 + j) > 0.5 ? 1.0 + 0.0im : -1.0 + 0.0im
    end
    time = ifft(freq)
    cp = n_cp == 0 ? ComplexF64[] : time[(end - n_cp + 1):end]
    return (time=vcat(cp, time), freq=freq, n_cp=n_cp)
end

@doc raw"""
    schmidl_cox(rx; half::Integer, search_range::Union{Nothing,Integer}=nothing) -> NamedTuple

Repeated-half frame metric and fractional carrier-frequency offset.

`half` is the length of one repeated half (``L = N/2`` for a preamble from [`schmidl_preamble`](@ref) with `n_fft = 2 * half`). At each 1-based start ``d``,

```math
\begin{aligned}
P(d) &= \sum_{m=0}^{L-1} r_{d+m}\, r_{d+m+L}^{*}, \\
R(d) &= \sum_{m=0}^{L-1} |r_{d+m+L}|^2, \\
M(d) &= \frac{|P(d)|^2}{R(d)^2 + 10^{-10}}.
\end{aligned}
```

`d` is the earliest maximiser of ``M``. Ties within a relative `1e-8` keep the earliest index, the same rule as [`OFDM.sto_estimate`](@ref). The default searches every start that can hold `2 * half` samples. `search_range` caps the number of starts.

The fractional CFO, in subcarrier spacings of an FFT of length `2 * half`, is

```math
\hat\varepsilon_f = -\angle P(d) / \pi.
```

The sign matches [`OFDM.cfo_correct`](@ref): a tone ``\exp(j 2\pi \varepsilon (i-1) / N)`` with ``N = 2L`` produces ``\hat\varepsilon_f = \varepsilon`` when ``\varepsilon \in (-1, 1)``. The lag is ``L`` samples, so the phase is ``\pi\varepsilon`` and the unambiguous interval is ``(-1, 1)``, twice the van de Beek interval. An even integer number of subcarrier spacings (``\varepsilon \leftarrow \varepsilon + 2z``) does not change ``P(d)`` and is invisible here. [`integer_cfo_estimate`](@ref) searches for that ``2z`` on a known `freq` grid. This function does not call [`OFDM.cfo_estimate`](@ref) and does not modify [`OFDM.ofdm_sync`](@ref).

A cyclic prefix makes ``M(d)`` a plateau rather than a single peak: every window that still sees the repeated halves scores like the true start. On a noiseless plateau the tie rule returns the earliest index, which is the first sample of the preamble. In noise the maximiser can sit anywhere on the plateau. Schmidl and Cox average along the plateau; this function does not. Treating `d + n_cp` as the FFT window is therefore only reliable when `d` really is the first sample. That is a teaching detector, not production frame sync.

`half` must be ≥ 2. The input must be finite and at least `2 * half` samples long.
"""
function schmidl_cox(rx::AbstractVector{<:Complex};
                     half::Integer,
                     search_range::Union{Nothing,Integer}=nothing)
    L = Int(half)
    L >= 2 || throw(ArgumentError("half must be ≥ 2, got $half"))
    x = ComplexF64.(rx)
    max_d = length(x) - 2 * L + 1
    max_d >= 1 || throw(ArgumentError(
        "need at least $(2 * L) samples for half=$L; got $(length(x))",
    ))
    nsearch = max_d
    if search_range !== nothing
        nsearch = Int(search_range)
        nsearch >= 1 || throw(ArgumentError("search_range must be ≥ 1, got $search_range"))
        nsearch = min(nsearch, max_d)
    end
    metric = Vector{Float64}(undef, nsearch)
    correlation = Vector{ComplexF64}(undef, nsearch)
    @inbounds for d in 1:nsearch
        s = 0.0 + 0.0im
        R = 0.0
        for m in 0:(L - 1)
            a = x[d + m]
            b = x[d + m + L]
            _check_finite_sample(a, "rx")
            _check_finite_sample(b, "rx")
            s += a * conj(b)
            R += abs2(b)
        end
        correlation[d] = s
        metric[d] = abs2(s) / (R^2 + 1e-10)
    end
    peak, idx = findmax(metric)
    while idx > 1 && metric[idx - 1] >= peak * (1 - 1e-8)
        idx -= 1
    end
    return (
        d=idx,
        epsilon_f=-angle(correlation[idx]) / π,
        metric=metric,
        correlation=correlation[idx],
    )
end

@doc raw"""
    integer_cfo_estimate(rx, freq; sto::Integer=1, z_max::Integer=4) -> Int

Even-integer carrier offset, in subcarrier spacings, after the fractional part has been removed.

`freq` is the length-`N` frequency grid of the repeated-half preamble (`N` even), the `freq` field of [`schmidl_preamble`](@ref). `sto` is the 1-based index of the first sample of an ISI-free FFT window of length `N` (the useful symbol, not a sample inside the cyclic prefix). The window is transformed with Julia's unnormalised `fft`.

The search is the coherent correlation against even shifts only,

```math
\hat m = \arg\max_{z \in \{-z_{\max},\ldots,z_{\max}\}}
\left|\sum_{k=1}^{N} Y[(k + 2z - 1) \bmod N + 1]\, X[k]^{*}\right|^2,
```

and the return value is ``\hat m = 2z``. A positive value means energy moved to higher bin indices, which is the same convention as [`OFDM.cfo_correct`](@ref): the time-domain tone ``\exp(j 2\pi \varepsilon (i-1) / N)`` shifts bin ``k`` toward bin ``k + \varepsilon``.

# Limits
``\hat\varepsilon_f`` from [`schmidl_cox`](@ref) must already have been removed with [`OFDM.cfo_correct`](@ref). A leftover fractional offset smears the bins and the maximiser is then not the integer CFO. The correlation is coherent, so it assumes a flat channel over the occupied bins. A linear phase from a timing error inside the window (an FFT that starts in the cyclic prefix, or past the symbol) reduces the peak and can move it. Only even integers are identifiable once the fractional estimator has absorbed the odd part into ``(-1, 1)``; the search does not try odd ``m``. Offsets with ``|2z| > 2 z_{\max}`` are outside the capture range and the maximiser is whatever even integer inside the window scores best, which need not be the truth. ``z_{\max} = 4`` covers ``\{-8,-6,\ldots,8\}``. The function does not change [`OFDM.cfo_estimate`](@ref), whose capture range stays ``(-1/2, 1/2)``.
"""
function integer_cfo_estimate(rx::AbstractVector{<:Complex},
                              freq::AbstractVector{<:Complex};
                              sto::Integer=1,
                              z_max::Integer=4)
    N = length(freq)
    (N >= 4 && iseven(N)) || throw(ArgumentError(
        "freq length must be even and ≥ 4, got $N",
    ))
    z_max = Int(z_max)
    z_max >= 0 || throw(ArgumentError("z_max must be ≥ 0, got $z_max"))
    sto = Int(sto)
    sto >= 1 || throw(ArgumentError("sto must be ≥ 1, got $sto"))
    last = sto + N - 1
    last <= length(rx) || throw(ArgumentError(
        "FFT window sto=$sto length $N runs past the end of a length-$(length(rx)) vector",
    ))
    Y = fft(ComplexF64.(rx[sto:last]))
    X = ComplexF64.(freq)
    best_m = 0
    best_ρ = -1.0
    for z in -z_max:z_max
        m = 2z
        ρ = 0.0 + 0.0im
        @inbounds for k in 1:N
            ρ += Y[mod1(k + m, N)] * conj(X[k])
        end
        val = abs2(ρ)
        if val > best_ρ
            best_ρ = val
            best_m = m
        end
    end
    return best_m
end

@doc raw"""
    schmidl_sync(rx, freq; n_cp::Integer=0, z_max::Integer=4, search_range::Union{Nothing,Integer}=nothing) -> NamedTuple

Fractional Schmidl–Cox correction followed by the even-integer search.

`freq` is the preamble grid from [`schmidl_preamble`](@ref). `half` is `length(freq) ÷ 2`. The steps are:

1. [`schmidl_cox`](@ref) on `rx`.
2. [`OFDM.cfo_correct`](@ref) by `epsilon_f`, with FFT length `length(freq)`. The phase reference is sample 1 of `rx`, the same convention as `cfo_correct`.
3. [`integer_cfo_estimate`](@ref) at `sto = d + n_cp`.
4. A second `cfo_correct` by that even integer, unless the integer is 0.

The total offset reported as `epsilon` is `epsilon_f + epsilon_i`. It equals a tone ``\exp(j 2\pi \varepsilon (i-1) / N)`` when ``\varepsilon \in (-1, 1) + 2\{-z_{\max},\ldots,z_{\max}\}``, the fractional part is not on the branch cut ``\pm 1``, the channel is flat, and `d` is the first sample of the preamble so that `d + n_cp` is the first sample of the useful symbol.

`n_cp` is **not** read from the waveform. Pass the same value that built the preamble. On a noiseless plateau [`schmidl_cox`](@ref) returns the earliest maximum, and `d + n_cp` is then the useful symbol. If noise moves `d` along the plateau, `sto` moves with it and the integer estimate is no longer the estimate documented above. This function does not average the plateau.

The returned `signal` is the whole input after both corrections, including any leading delay and the cyclic prefix. It is not shifted by `d`. That matches [`OFDM.ofdm_sync`](@ref), which also returns an unshifted buffer — and which is still the van de Beek estimator, not this one. `ofdm_sync` is not called.

Returns `(signal, d, sto, epsilon_f, epsilon_i, epsilon, metric)`.
"""
function schmidl_sync(rx::AbstractVector{<:Complex},
                      freq::AbstractVector{<:Complex};
                      n_cp::Integer=0,
                      z_max::Integer=4,
                      search_range::Union{Nothing,Integer}=nothing)
    N = length(freq)
    (N >= 4 && iseven(N)) || throw(ArgumentError(
        "freq length must be even and ≥ 4, got $N",
    ))
    n_cp = Int(n_cp)
    n_cp >= 0 || throw(ArgumentError("n_cp must be ≥ 0, got $n_cp"))
    sc = schmidl_cox(rx; half=N ÷ 2, search_range=search_range)
    corrected = OFDM.cfo_correct(rx, sc.epsilon_f, N)
    sto = sc.d + n_cp
    ε_i = integer_cfo_estimate(corrected, freq; sto=sto, z_max=z_max)
    signal = ε_i == 0 ? corrected : OFDM.cfo_correct(corrected, float(ε_i), N)
    return (
        signal=signal,
        d=sc.d,
        sto=sto,
        epsilon_f=sc.epsilon_f,
        epsilon_i=ε_i,
        epsilon=sc.epsilon_f + ε_i,
        metric=sc.metric,
    )
end

end # module Sync
