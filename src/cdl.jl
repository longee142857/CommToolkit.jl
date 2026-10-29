# Teaching clustered delay line (CDL): the delay–power columns of
# 3GPP TR 38.901 V16.1.0 Tables 7.7.1-1/2/3, delay scaling by Eq. (7.7-1),
# sample quantization with tap merging, one Rayleigh gain per tap (reusing
# jakes_fading for the optional Doppler), and a time-varying convolution.
# This file is included from channel.jl, inside module Channel.

# TR 38.901 V16.1.0, Table 7.7.1-1 (CDL-A, NLOS): normalized delays and powers.
const _CDL_A_DELAYS = [
    0.0000, 0.3819, 0.4025, 0.5868, 0.4610, 0.5375, 0.6708, 0.5750, 0.7618,
    1.5375, 1.8978, 2.2242, 2.1718, 2.4942, 2.5119, 3.0582, 4.0810, 4.4579,
    4.5695, 4.7966, 5.0066, 5.3043, 9.6586,
]
const _CDL_A_POWERS_DB = [
    -13.4, 0.0, -2.2, -4.0, -6.0, -8.2, -9.9, -10.5, -7.5, -15.9, -6.6,
    -16.7, -12.4, -15.2, -10.8, -11.3, -12.7, -16.2, -18.3, -18.9, -16.6,
    -19.9, -29.7,
]

# TR 38.901 V16.1.0, Table 7.7.1-2 (CDL-B, NLOS).
const _CDL_B_DELAYS = [
    0.0000, 0.1072, 0.2155, 0.2095, 0.2870, 0.2986, 0.3752, 0.5055, 0.3681,
    0.3697, 0.5700, 0.5283, 1.1021, 1.2756, 1.5474, 1.7842, 2.0169, 2.8294,
    3.0219, 3.6187, 4.1067, 4.2790, 4.7834,
]
const _CDL_B_POWERS_DB = [
    0.0, -2.2, -4.0, -3.2, -9.8, -1.2, -3.4, -5.2, -7.6, -3.0, -8.9, -9.0,
    -4.8, -5.7, -7.5, -1.9, -7.6, -12.2, -9.8, -11.4, -14.9, -9.2, -11.3,
]

# TR 38.901 V16.1.0, Table 7.7.1-3 (CDL-C, NLOS).
const _CDL_C_DELAYS = [
    0.0000, 0.2099, 0.2219, 0.2329, 0.2176, 0.6366, 0.6448, 0.6560, 0.6584,
    0.7935, 0.8213, 0.9336, 1.2285, 1.3083, 2.1704, 2.7105, 4.2589, 4.6003,
    5.4902, 5.6077, 6.3065, 6.6374, 7.0427, 8.6523,
]
const _CDL_C_POWERS_DB = [
    -4.4, -1.2, -3.5, -5.2, -2.5, 0.0, -2.2, -3.9, -7.4, -7.1, -10.7, -11.1,
    -5.1, -6.8, -8.7, -13.2, -13.9, -13.9, -15.8, -17.1, -16.0, -15.7,
    -21.6, -22.8,
]

const _CDL_TABLES = (
    A=(delays=_CDL_A_DELAYS, powers_db=_CDL_A_POWERS_DB),
    B=(delays=_CDL_B_DELAYS, powers_db=_CDL_B_POWERS_DB),
    C=(delays=_CDL_C_DELAYS, powers_db=_CDL_C_POWERS_DB),
)

@doc raw"""
    CDLProfile(delays, powers_db; name=:custom)

A clustered-delay-line delay–power profile. `delays[i]` is the normalized delay of cluster `i` (unitless, in units of the RMS delay spread, as tabulated in 3GPP TR 38.901), and `powers_db[i]` is its average power in decibels. Cluster order is kept as given; the tabulated profiles are not monotonic in delay.

# Fields
- `name::Symbol`: profile label (`:A`, `:B`, `:C` from [`cdl_profile`](@ref), `:custom` otherwise).
- `delays::Vector{Float64}`: normalized delays. Finite and ``\ge 0``. [`cdl_taps`](@ref) multiplies them by its `delay_spread` to get seconds, the scaling of TR 38.901, eq. (7.7-1).
- `powers_db::Vector{Float64}`: cluster powers in decibels, finite. They are normalized to sum 1 in linear scale inside [`cdl_taps`](@ref), so a constant offset in dB changes nothing.

The two vectors must have the same non-zero length. This is the delay–power part of a CDL only: no angles of departure or arrival, no cluster angle spreads, no cross-polarization ratios, no antenna geometry. The LOS profiles CDL-D and CDL-E are not part of the teaching subset.

# References
3GPP TR 38.901, V16.1.0, §7.7.1, Tables 7.7.1-1 to 7.7.1-3 (delay and power columns only), and §7.7.3, eq. (7.7-1).
"""
struct CDLProfile
    name::Symbol
    delays::Vector{Float64}
    powers_db::Vector{Float64}
    function CDLProfile(name::Symbol, delays::Vector{Float64}, powers_db::Vector{Float64})
        length(delays) == length(powers_db) || throw(ArgumentError(
            "delays and powers_db must have the same length",
        ))
        isempty(delays) && throw(ArgumentError("a CDL profile needs at least one cluster"))
        all(d -> isfinite(d) && d >= 0, delays) || throw(ArgumentError(
            "delays must be finite and ≥ 0",
        ))
        all(isfinite, powers_db) || throw(ArgumentError("powers_db must be finite"))
        return new(name, delays, powers_db)
    end
end

function CDLProfile(delays::AbstractVector{<:Real}, powers_db::AbstractVector{<:Real};
                    name::Symbol=:custom)
    return CDLProfile(name, Float64.(delays), Float64.(powers_db))
end

@doc raw"""
    cdl_profile(name::Symbol) -> CDLProfile

The tabulated NLOS clustered-delay-line profiles of 3GPP TR 38.901: `cdl_profile(:A)`, `cdl_profile(:B)`, and `cdl_profile(:C)` return Tables 7.7.1-1, 7.7.1-2, and 7.7.1-3 of V16.1.0 (23, 23, and 24 clusters), restricted to the normalized-delay and power columns.

The tabulated delays are normalized so that the power-weighted RMS delay spread

```math
\mathrm{DS} = \sqrt{\frac{\sum_i p_i \tau_i^2}{\sum_i p_i}
  - \left(\frac{\sum_i p_i \tau_i}{\sum_i p_i}\right)^2}, \qquad p_i = 10^{P_i/10},
```

is 1 for each table; [`cdl_taps`](@ref) recovers absolute delays in seconds by multiplying by the wanted RMS delay spread (TR 38.901, §7.7.3, eq. (7.7-1)).

The LOS profiles CDL-D and CDL-E are rejected: this teaching subset is NLOS-only and has no specular K-factor path. Angles, cluster spreads, and XPR are not part of a [`CDLProfile`](@ref).
"""
function cdl_profile(name::Symbol)
    name in (:D, :E) && throw(ArgumentError(
        "CDL-$name is an LOS profile; the teaching subset is the NLOS CDL-A/B/C only",
    ))
    haskey(_CDL_TABLES, name) || throw(ArgumentError(
        "unknown CDL profile $name; expected :A, :B, or :C",
    ))
    table = _CDL_TABLES[name]
    return CDLProfile(name, copy(table.delays), copy(table.powers_db))
end

"""Power-weighted RMS delay spread of `delays` with linear `powers`. Internal."""
function _rms_delay_spread(delays::Vector{Float64}, powers::Vector{Float64})
    s = sum(powers)
    μ = sum(d .* p for (d, p) in zip(delays, powers)) / s
    m2 = sum(d^2 * p for (d, p) in zip(delays, powers)) / s
    return sqrt(max(m2 - μ^2, 0.0))
end

@doc raw"""
    CDLTaps(delays, gains, powers)

Sample-quantized taps of a clustered delay line, as produced by [`cdl_taps`](@ref).

# Fields
- `delays::Vector{Int}`: tap delays in samples, ascending, unique, ``\ge 0``.
- `gains::Matrix{ComplexF64}`: `gains[i, t]` is the complex coefficient of tap `i` at input-sample index `t`. One column is a static drop; [`cdl_taps`](@ref) returns `n` columns, one per input sample, and a single column when `fd = 0`.
- `powers::Vector{Float64}`: average tap powers, ``\mathbb{E}[|\texttt{gains}[i, t]|^2]``, normalized so ``\sum_i \texttt{powers}[i] = 1``.

A static drop plugs into the OFDM chain: `TDLChannel(taps.gains[:, 1], taps.delays)` is a valid [`CommToolkit.OFDM.TDLChannel`](@ref) for `apply_channel` and `apply_multipath`.
"""
struct CDLTaps
    delays::Vector{Int}
    gains::Matrix{ComplexF64}
    powers::Vector{Float64}
    function CDLTaps(delays::Vector{Int}, gains::Matrix{ComplexF64}, powers::Vector{Float64})
        length(delays) == length(powers) == size(gains, 1) || throw(ArgumentError(
            "delays, powers, and rows of gains must have the same length",
        ))
        isempty(delays) && throw(ArgumentError("CDLTaps needs at least one tap"))
        all(d -> d >= 0, delays) || throw(ArgumentError("delays must be ≥ 0"))
        issorted(delays) && allunique(delays) || throw(ArgumentError(
            "delays must be ascending and unique",
        ))
        all(p -> isfinite(p) && p >= 0, powers) || throw(ArgumentError(
            "powers must be finite and ≥ 0",
        ))
        return new(delays, gains, powers)
    end
end

@doc raw"""
    cdl_taps(profile::CDLProfile, n::Integer, Ts::Real;
             delay_spread=100e-9, fd=0.0, n_taps=nothing) -> CDLTaps

Generate `n` samples of a teaching clustered delay line at sample interval `Ts` seconds.

# Procedure
1. Absolute delays: ``\tau_i = \texttt{delays}[i] \times \texttt{delay\_spread}`` seconds (TR 38.901, §7.7.3, eq. (7.7-1); the tabulated normalized delays have unit RMS delay spread, so `delay_spread` is the RMS delay spread of the result). `delay_spread` is finite and positive, in seconds; the default `100e-9` is the "nominal delay spread" of TR 38.901, Table 7.7.3-1. With `delay_spread = 1.0` the profile delays are read as seconds directly.
2. Sample quantization: cluster ``i`` lands on sample delay ``d_i = \texttt{round}(\tau_i / T_s)`` (ties to even). Clusters landing on the same sample **merge**: their powers add, because independent Rayleigh gains summing on one tap have the summed power. When `Ts` is large compared with the delay spread, everything merges onto sample 0 and the channel is flat.
3. Optional thinning: `n_taps` keeps only that many strongest merged taps (largest power, ties broken by the smaller delay) and renormalizes. This is a teaching convenience for short equalizers and demos; it is not part of TR 38.901. The default `nothing` keeps every merged tap.
4. Gains: tap ``i`` gets ``g_i(t) = \sqrt{p_i}\, h_i(t)`` with ``p_i`` its normalized power (``\sum_i p_i = 1``) and `h_i` an independent [`jakes_fading`](@ref)`(n, fd, Ts)` draw, in ascending delay order. `fd = 0` is a static drop: one ``\mathcal{CN}(0, 1)`` draw per tap and **one gain column** (the same draws as the first column of the general case), which [`apply_cdl`](@ref) broadcasts over all input samples. `fd > 0` returns `n` columns, a flat Clarke spectrum per tap with `fd * Ts < 1/2`. There is one gain per tap; the full 38.901 coefficient generation (20 rays per cluster from the angle tables, Steps 10–11 of §7.5) is **not** implemented.

`n ≥ 0`. `n = 0` returns taps with zero gain columns and does not draw. `Ts` is finite and positive.

# Random draws
The global `randn` stream, one [`jakes_fading`](@ref) call per tap in ascending delay order (that function's own draw order applies inside each call; `fd = 0` is one `randn` pair per tap). The same `Random.seed!` reproduces the same taps. `n = 0` does not draw.

# What this is not
Not the geometric TR 38.901 CDL: no departure/arrival angles, no cluster angle spreads, no XPR, no antenna patterns, no spatial consistency, and no MIMO extension. Not an NR compliance model. The per-tap Doppler is Clarke's spectrum, not the ray-sum Doppler of §7.5.
"""
function cdl_taps(profile::CDLProfile, n::Integer, Ts::Real;
                  delay_spread::Real=100e-9, fd::Real=0.0, n_taps=nothing)
    n = Int(n)
    n >= 0 || throw(ArgumentError("n must be ≥ 0, got $n"))
    t = float(Ts)
    (isfinite(t) && t > 0) || throw(ArgumentError(
        "Ts must be a finite sample interval > 0, in seconds, got $Ts",
    ))
    ds = float(delay_spread)
    (isfinite(ds) && ds > 0) || throw(ArgumentError(
        "delay_spread must be a finite RMS delay spread > 0, in seconds, got $delay_spread",
    ))
    f, _ = _jakes_params(fd, t)
    if n_taps !== nothing
        n_taps = Int(n_taps)
        n_taps >= 1 || throw(ArgumentError("n_taps must be ≥ 1, got $n_taps"))
    end

    # Normalized power per cluster.
    powers = exp10.(profile.powers_db ./ 10)
    powers ./= sum(powers)

    # Quantize cluster delays to samples and merge collisions (powers add).
    merged = Dict{Int,Float64}()
    @inbounds for i in eachindex(profile.delays)
        d = round(Int, profile.delays[i] * ds / t)
        merged[d] = get(merged, d, 0.0) + powers[i]
    end
    order = sort!(collect(keys(merged)))
    delays = order
    tap_powers = [merged[d] for d in order]

    if n_taps !== nothing && length(delays) > n_taps
        keep = sortperm(tap_powers; rev=true, alg=Base.DEFAULT_STABLE)[1:n_taps]
        # Stable descending power keeps the earlier (smaller) delay on ties.
        keep = sort(keep; by=i -> delays[i])
        delays = delays[keep]
        tap_powers = tap_powers[keep]
    end
    tap_powers ./= sum(tap_powers)

    ncols = f == 0 ? min(n, 1) : n
    gains = Matrix{ComplexF64}(undef, length(delays), ncols)
    for i in eachindex(delays)
        gains[i, :] .= sqrt(tap_powers[i]) .* jakes_fading(ncols, f, t)
    end
    return CDLTaps(delays, gains, tap_powers)
end

@doc raw"""
    apply_cdl(signal, taps::CDLTaps) -> Vector{ComplexF64}

Time-varying tapped-delay-line convolution of `signal` with `taps`,

```math
y[k] = \sum_{i:\, 1 \le k - d_i \le n} g_i[k - d_i]\; x[k - d_i],
\qquad k = 1, \ldots, n + d_{\max},
```

where ``d_i`` is `taps.delays[i]` and ``n`` is `length(signal)`. The gain is indexed by the **input** sample, so the tail beyond sample ``n`` uses only columns that exist. The result has length ``n + \texttt{maximum(taps.delays)}``, the full linear convolution, matching [`CommToolkit.OFDM.apply_multipath`](@ref) and [`apply_fir`](@ref).

`taps.gains` must have one column (a static drop, broadcast over all input samples) or exactly `length(signal)` columns. `signal` must be non-empty; real inputs are promoted to complex. Taps are used as given: this function does not draw from the random stream and does not normalize anything.
"""
function apply_cdl(signal::AbstractVector{<:Number}, taps::CDLTaps)
    n = length(signal)
    n >= 1 || throw(ArgumentError("signal must be non-empty"))
    ncols = size(taps.gains, 2)
    (ncols == 1 || ncols == n) || throw(ArgumentError(
        "taps.gains must have 1 column (static) or length(signal) = $n columns, got $ncols",
    ))
    x = ComplexF64.(signal)
    max_delay = maximum(taps.delays)
    out = zeros(ComplexF64, n + max_delay)
    static = ncols == 1
    @inbounds for i in eachindex(taps.delays)
        d = taps.delays[i]
        for j in 1:n
            g = static ? taps.gains[i, 1] : taps.gains[i, j]
            out[j + d] += g * x[j]
        end
    end
    return out
end

@doc raw"""
    cdl_channel(symbols, snr_db::Real, bits_per_symbol::Integer;
                Ts, profile=cdl_profile(:A), delay_spread=100e-9, fd=0.0, n_taps=nothing)
        -> (rx, taps)

Teaching clustered-delay-line fading, then AWGN,

```math
y = \mathrm{CDL}(x) + n,
```

with `taps` from [`cdl_taps`](@ref)`(profile, length(symbols), Ts; ...)` and the convolution from [`apply_cdl`](@ref). Real `symbols` are promoted to complex. ``E_s = 1`` is assumed and is not estimated. The noise follows the same ``E_b/N_0`` convention as [`awgn`](@ref), with ``E_b = 1/\texttt{bits\_per\_symbol}``, added on **every** sample of the convolution tail as well. Because the tap powers are normalized, ``\sum_i \mathbb{E}[|g_i|^2] = 1``, `snr_db` is the average ``E_b/N_0`` at the channel output, in the same sense as [`rayleigh_channel`](@ref). ``\mathrm{snr\_db} = +\infty`` disables the noise.

`Ts` is a **required keyword**: the sample interval in seconds (the symbol period for symbol-rate input). Together with `delay_spread` it decides how many clusters resolve onto separate taps; a slow symbol rate merges the whole profile onto tap 0 and the channel degenerates to flat Rayleigh. `profile`, `delay_spread`, `fd`, and `n_taps` are forwarded to [`cdl_taps`](@ref).

# Returns
- `rx`: the faded observation of length `length(symbols) + maximum(taps.delays)`, plus noise.
- `taps`: the [`CDLTaps`](@ref) that were applied. `taps.gains[:, 1]` and `taps.delays` of a static drop (`fd = 0`, or any single-column slice) build a [`CommToolkit.OFDM.TDLChannel`](@ref) for the OFDM receivers.

# Closed form
Only the fully merged case has one: when every cluster quantizes to sample 0, the composite coefficient is ``\mathcal{CN}(0, 1)`` and ideal-CSI coherent BPSK has the flat-Rayleigh bit error rate on [`rayleigh_channel`](@ref). The dispersive case has no golden BER curve in the tests; OFDM checks are perfect-recovery smoke tests above the cyclic prefix.

This is not the geometric TR 38.901 CDL (no angles, no ray sums, no MIMO), not an NR compliance channel, and the per-tap Doppler is the Clarke spectrum of [`jakes_fading`](@ref), not the 38.901 ray-sum Doppler.
"""
function cdl_channel(syms::AbstractVector, snr_db::Real, bits_per_symbol::Integer;
                     Ts::Real, profile::CDLProfile=cdl_profile(:A),
                     delay_spread::Real=100e-9, fd::Real=0.0, n_taps=nothing)
    bits_per_symbol >= 1 || throw(ArgumentError("bits_per_symbol must be ≥ 1"))
    syms_c = ComplexF64.(syms)
    taps = cdl_taps(profile, length(syms_c), Ts;
                    delay_spread=delay_spread, fd=fd, n_taps=n_taps)
    isempty(syms_c) && return ComplexF64[], taps
    return awgn(apply_cdl(syms_c, taps), snr_db, bits_per_symbol), taps
end
