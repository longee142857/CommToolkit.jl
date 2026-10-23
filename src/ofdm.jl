module OFDM

using FFTW
using ..Modulation: modulate, demodulate, QPSK, QAM

export OFDMParams, TDLChannel, modulation,
    generate_ofdm_symbol, ofdm_tx,
    apply_multipath, apply_channel,
    extract_ofdm_symbols, ofdm_rx_known, ofdm_rx_ls, ofdm_rx_mmse,
    interpolate_channel,
    sto_estimate, cfo_estimate, cfo_correct, ofdm_sync

@doc raw"""
    OFDMParams(; N_fft=64, N_cp=16, mod_order=4, n_guard=6, pilot_spacing=6)

Parameters of one baseband CP-OFDM transmitter.

# Arguments
- `N_fft`: FFT length.
- `N_cp`: cyclic-prefix length in samples. A static multipath of delay at most `N_cp` samples is absorbed by the prefix (cyclic convolution on the FFT window).
- `mod_order`: `4` selects QPSK; any larger square power of two selects `QAM(mod_order)`.
- `n_guard`: number of guard subcarriers at each FFT edge. The DC bin (`N_fft ÷ 2 + 1`, 1-based) is also nulled.
- `pilot_spacing`: pilots are placed on every `pilot_spacing`-th *used* subcarrier (used = active, excluding DC), not on every `pilot_spacing`-th FFT bin.

# Fields
`N_data`, `N_pilot`, `bits_per_symbol`, `pilot_indices`, `data_indices`, `null_indices` (DC only; the guards are simply unused). Indices are 1-based FFT bins in Julia's `fft` order (bin 1 is DC of a conventional zero-based FFT only after an `fftshift`; here bin `k` has digital frequency `(k - 1) / N_fft`).

This is a textbook subcarrier map. It is not a 3GPP NR resource grid and it does not implement DMRS, PTRS, or any O-RAN fronthaul mapping.
"""
struct OFDMParams
    N_fft::Int
    N_data::Int
    N_pilot::Int
    N_cp::Int
    mod_order::Int
    bits_per_symbol::Int
    pilot_spacing::Int
    pilot_indices::Vector{Int}
    data_indices::Vector{Int}
    null_indices::Vector{Int}
end

function OFDMParams(; N_fft::Int=64, N_cp::Int=16, mod_order::Int=4,
                    n_guard::Int=6, pilot_spacing::Int=6)
    N_fft >= 4 || throw(ArgumentError("N_fft must be ≥ 4"))
    N_cp >= 0 || throw(ArgumentError("N_cp must be ≥ 0"))
    n_guard >= 0 || throw(ArgumentError("n_guard must be ≥ 0"))
    pilot_spacing >= 1 || throw(ArgumentError("pilot_spacing must be ≥ 1"))
    if mod_order == 4
        bps = 2
    else
        bps = Int(log2(mod_order))
        # QAM constructor validates the square power-of-two constraint.
        QAM(mod_order)
    end
    2 * n_guard < N_fft - 1 || throw(ArgumentError("guard bands leave no active subcarriers"))

    active = collect((n_guard + 1):(N_fft - n_guard))
    dc_idx = N_fft ÷ 2 + 1
    used = filter(!=(dc_idx), active)
    pilot_positions = collect(pilot_spacing:pilot_spacing:length(used))
    pilot_indices = sort(used[pilot_positions])
    data_indices = sort(setdiff(used, pilot_indices))
    isempty(data_indices) && throw(ArgumentError("pilot map leaves no data subcarriers"))

    return OFDMParams(N_fft, length(data_indices), length(pilot_indices),
                      N_cp, mod_order, bps, pilot_spacing,
                      pilot_indices, data_indices, [dc_idx])
end

@doc raw"""
    modulation(params::OFDMParams) -> QPSK or QAM

`QPSK()` when `params.mod_order == 4`, otherwise `QAM(params.mod_order)`.
"""
modulation(params::OFDMParams) = params.mod_order == 4 ? QPSK() : QAM(params.mod_order)

# 31-bit LCG (Park–Miller / glibc-style constants). Used so pilots and TDL
# phases do not depend on the global RNG or on the Julia version.
function _lcg_rand(seed::Int)
    r = (seed * 1103515245 + 12345) & 0x7fffffff
    return Float64(r) / 2147483648.0
end

function _generate_pilots(n::Int)
    return [_lcg_rand(42 + i) > 0.5 ? 1.0 + 0.0im : -1.0 + 0.0im for i in 1:n]
end

@doc raw"""
    generate_ofdm_symbol(params::OFDMParams, data_syms, pilot_syms) -> Vector{ComplexF64}

Build one CP-OFDM symbol. `data_syms` and `pilot_syms` fill `params.data_indices` and `params.pilot_indices`. Every other bin, including DC and the guards, is zero.

The IFFT is Julia's `FFTW.ifft` (forward `fft` is unnormalised, `ifft` divides by `N_fft`), so `fft(ifft(X)) ≈ X`. The cyclic prefix is the last `N_cp` time-domain samples, prepended. The returned length is `N_fft + N_cp`.
"""
function generate_ofdm_symbol(params::OFDMParams,
                              data_syms::AbstractVector{<:Complex},
                              pilot_syms::AbstractVector{<:Complex})
    length(data_syms) == params.N_data || throw(ArgumentError(
        "expected $(params.N_data) data symbols, got $(length(data_syms))",
    ))
    length(pilot_syms) == params.N_pilot || throw(ArgumentError(
        "expected $(params.N_pilot) pilots, got $(length(pilot_syms))",
    ))
    N = params.N_fft
    freq = zeros(ComplexF64, N)
    freq[params.data_indices] .= data_syms
    freq[params.pilot_indices] .= pilot_syms
    time_domain = ifft(freq)
    cp = time_domain[(end - params.N_cp + 1):end]
    return vcat(cp, time_domain)
end

@doc raw"""
    ofdm_tx(params::OFDMParams, bits::AbstractVector{<:Integer}) -> (signal, pilots)

Modulate `bits` onto successive OFDM symbols.

# Arguments
- `bits`: entries in ``\{0, 1\}``. The vector is zero-padded up to a whole number of OFDM symbols. Those padding bits are not part of the caller's payload; receivers return them and the caller should drop them.

# Returns
- `signal`: concatenated CP-OFDM symbols, `ComplexF64`.
- `pilots`: the deterministic length-`N_pilot` BPSK pilot sequence, reused on every symbol. The sequence is an LCG, not a standard gold or Zadoff–Chu sequence.
"""
function ofdm_tx(params::OFDMParams, bits::AbstractVector{<:Integer})
    all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("bits must be 0 or 1"))
    bits_per_ofdm = params.N_data * params.bits_per_symbol
    n_sym = max(ceil(Int, length(bits) / bits_per_ofdm), 1)
    padded = Vector{Int}(undef, n_sym * bits_per_ofdm)
    fill!(padded, 0)
    copyto!(padded, 1, bits, 1, length(bits))
    pilots = _generate_pilots(params.N_pilot)
    scheme = modulation(params)
    tx = ComplexF64[]
    for s in 1:n_sym
        idx = (s - 1) * bits_per_ofdm
        data_bits = padded[(idx + 1):(idx + bits_per_ofdm)]
        append!(tx, generate_ofdm_symbol(params, modulate(data_bits, scheme), pilots))
    end
    return tx, pilots
end

@doc raw"""
    TDLChannel(taps, delays)
    TDLChannel(; delays_samples=[0, 1, 3], gain_dB=[0.0, -2.0, -5.0])

A static tapped-delay-line channel. `taps[i]` is the complex coefficient of a delay of `delays[i]` samples.

The keyword constructor sets the tap magnitude to ``10^{\mathrm{gain\_dB}/20}`` and a deterministic phase ``e^{j 2\pi u_i}`` from the same LCG used for pilots (`u_i = _lcg(123 + i)`). Tap powers are **not** renormalised to one. This is a static TDL, not a 3GPP TDL-A/B/C/D/E profile and not a fading process.
"""
struct TDLChannel
    taps::Vector{ComplexF64}
    delays::Vector{Int}
    function TDLChannel(taps::AbstractVector{<:Complex}, delays::AbstractVector{<:Integer})
        length(taps) == length(delays) ||
            throw(ArgumentError("taps and delays must have the same length"))
        isempty(taps) && throw(ArgumentError("TDL channel needs at least one tap"))
        any(d -> d < 0, delays) && throw(ArgumentError("delays must be ≥ 0"))
        return new(ComplexF64.(taps), Int.(delays))
    end
end

function TDLChannel(; delays_samples::Vector{Int}=[0, 1, 3],
                    gain_dB::Vector{<:Real}=[0.0, -2.0, -5.0])
    length(delays_samples) == length(gain_dB) ||
        throw(ArgumentError("delays_samples and gain_dB must have the same length"))
    taps = ComplexF64[
        10^(g / 20) * cispi(2 * _lcg_rand(123 + i))
        for (i, g) in enumerate(gain_dB)
    ]
    return TDLChannel(taps, delays_samples)
end

@doc raw"""
    apply_multipath(signal, channel::TDLChannel) -> Vector{ComplexF64}

Linear convolution of `signal` with the TDL. The result is longer than `signal` by `maximum(channel.delays)`.
"""
function apply_multipath(signal::AbstractVector{<:Complex}, channel::TDLChannel)
    L = length(signal)
    max_delay = maximum(channel.delays)
    output = zeros(ComplexF64, L + max_delay)
    @inbounds for (tap, delay) in zip(channel.taps, channel.delays)
        for i in 1:L
            output[i + delay] += tap * signal[i]
        end
    end
    return output
end

function _channel_freq_response(channel::TDLChannel, n::Int)
    H = zeros(ComplexF64, n)
    @inbounds for k in 1:n
        f = (k - 1) / n
        for (tap, delay) in zip(channel.taps, channel.delays)
            H[k] += tap * exp(-im * 2π * f * delay)
        end
    end
    return H
end

function _awgn_ofdm(signal::AbstractVector{ComplexF64}, snr_db::Real)
    P = sum(abs2, signal) / length(signal)
    snr_lin = exp10(float(snr_db) / 10)
    noise_power = P / snr_lin
    noise_std = sqrt(noise_power / 2)
    n = length(signal)
    return signal .+ noise_std .* randn(n) .+ im .* noise_std .* randn(n)
end

@doc raw"""
    apply_channel(signal, channel::TDLChannel, snr_db::Real, params::OFDMParams) -> (rx, H_freq)

Apply the static TDL and then AWGN.

# Arguments
- `snr_db`: signal-to-noise ratio in decibels defined on the **time-domain samples after fading**, ``\mathrm{SNR} = P_{\mathrm{sample}} / \mathbb{E}[|n|^2]``. It is not ``E_b/N_0``. The cyclic prefix, the null subcarriers, and the tap powers all sit inside ``P_{\mathrm{sample}}``. `Inf` disables the noise.
- The multipath tail past the original signal length is discarded. When `maximum(delays) ≤ N_cp` this does not disturb the cyclic convolution inside each kept OFDM symbol, except that energy leaking past the last sample of the burst is dropped.

# Returns
- `rx`: noisy waveform, same length as the input measured in whole OFDM symbols (`N_fft + N_cp` each).
- `H_freq`: length-`N_fft` frequency response ``H[k] = \sum_\ell a_\ell e^{-j 2\pi (k-1) \ell / N}``, which is exactly `fft` of the zero-padded impulse response (Julia's unnormalised FFT). Ideal equalisation divides the FFT bin by this value.
"""
function apply_channel(signal::AbstractVector{<:Complex}, channel::TDLChannel,
                       snr_db::Real, params::OFDMParams)
    faded = apply_multipath(ComplexF64.(signal), channel)
    sym_len = params.N_fft + params.N_cp
    sym_len > 0 || throw(ArgumentError("OFDM symbol length must be positive"))
    n_sym = length(signal) ÷ sym_len
    expected = n_sym * sym_len
    rx = if length(faded) > expected
        faded[1:expected]
    elseif length(faded) < expected
        vcat(faded, zeros(ComplexF64, expected - length(faded)))
    else
        faded
    end
    rx = _awgn_ofdm(rx, snr_db)
    return rx, _channel_freq_response(channel, params.N_fft)
end

@doc raw"""
    extract_ofdm_symbols(rx_signal, params::OFDMParams) -> Vector{Vector{ComplexF64}}

Drop the cyclic prefix of each symbol and return `fft` of the `N_fft` samples. Incomplete trailing samples are ignored. The FFT window starts `N_cp` samples into each symbol, which is correct when the burst is aligned to sample 1 and the delay spread fits in the prefix.
"""
function extract_ofdm_symbols(rx_signal::AbstractVector{<:Complex}, params::OFDMParams)
    N = params.N_fft
    N_cp = params.N_cp
    sym_len = N + N_cp
    n_sym = length(rx_signal) ÷ sym_len
    rx = ComplexF64.(rx_signal)
    return map(1:n_sym) do i
        idx = (i - 1) * sym_len + N_cp + 1
        fft(rx[idx:(idx + N - 1)])
    end
end

@doc raw"""
    ofdm_rx_known(rx_signal, params, H_freq, pilot_syms) -> Vector{Int}

Equalise data subcarriers by the known frequency response `H_freq ./` and hard-demodulate. `pilot_syms` is accepted for signature symmetry with the estimators and is not used. Zero-padding bits added by [`ofdm_tx`](@ref) are included in the result.
"""
function ofdm_rx_known(rx_signal::AbstractVector{<:Complex},
                       params::OFDMParams,
                       H_freq::AbstractVector{<:Complex},
                       pilot_syms::AbstractVector{<:Complex})
    length(H_freq) == params.N_fft ||
        throw(ArgumentError("H_freq must have length N_fft = $(params.N_fft)"))
    freq_syms = extract_ofdm_symbols(rx_signal, params)
    scheme = modulation(params)
    rx_bits = Int[]
    for fd in freq_syms
        data_eq = fd[params.data_indices] ./ H_freq[params.data_indices]
        append!(rx_bits, demodulate(data_eq, scheme))
    end
    return rx_bits
end

@doc raw"""
    interpolate_channel(H_pilot, params::OFDMParams) -> Vector{ComplexF64}

Copy pilot estimates onto `params.pilot_indices` and linearly interpolate, in the FFT-bin index, onto `params.data_indices`. Bins outside the pilot span are held at the nearest pilot (zero-order). Guard and DC bins stay zero.

Linear interpolation in the bin index is exact only for a frequency response that is linear between pilots. A multipath response is a sum of complex exponentials, so this estimator has a model error and an error floor that does not vanish as the SNR grows. It is not the Wiener interpolator of Edfors et al., IEEE Trans. Commun., 46(7):931–939, 1998.
"""
function interpolate_channel(H_pilot::AbstractVector{<:Complex}, params::OFDMParams)
    length(H_pilot) == params.N_pilot || throw(ArgumentError(
        "expected $(params.N_pilot) pilot estimates, got $(length(H_pilot))",
    ))
    N = params.N_fft
    H_est = zeros(ComplexF64, N)
    pilot_idx = params.pilot_indices
    @inbounds for (i, idx) in enumerate(pilot_idx)
        H_est[idx] = H_pilot[i]
    end
    @inbounds for d_idx in params.data_indices
        left = findlast(<(d_idx), pilot_idx)
        right = findfirst(>(d_idx), pilot_idx)
        if left === nothing
            H_est[d_idx] = H_pilot[right]
        elseif right === nothing
            H_est[d_idx] = H_pilot[left]
        else
            α = (d_idx - pilot_idx[left]) / (pilot_idx[right] - pilot_idx[left])
            H_est[d_idx] = H_pilot[left] * (1 - α) + H_pilot[right] * α
        end
    end
    return H_est
end

@doc raw"""
    ofdm_rx_ls(rx_signal, params, pilot_syms) -> Vector{Int}

Least-squares pilot estimates ``\hat H_p = Y_p / X_p``, then [`interpolate_channel`](@ref), zero-forcing on the data bins, and hard demodulation. `pilot_syms` must be the sequence returned by [`ofdm_tx`](@ref).
"""
function ofdm_rx_ls(rx_signal::AbstractVector{<:Complex},
                    params::OFDMParams,
                    pilot_syms::AbstractVector{<:Complex})
    freq_syms = extract_ofdm_symbols(rx_signal, params)
    scheme = modulation(params)
    rx_bits = Int[]
    for fd in freq_syms
        H_pilot = fd[params.pilot_indices] ./ pilot_syms
        H_est = interpolate_channel(H_pilot, params)
        data_eq = fd[params.data_indices] ./ H_est[params.data_indices]
        append!(rx_bits, demodulate(data_eq, scheme))
    end
    return rx_bits
end

@doc raw"""
    ofdm_rx_mmse(rx_signal, params, pilot_syms, snr_db::Real) -> Vector{Int}

Shrink each least-squares pilot estimate and then do the same linear interpolation as [`ofdm_rx_ls`](@ref),

```math
\hat H_p^{\mathrm{shrunk}}
  = \frac{|\hat H_p^{\mathrm{LS}}|^2}{|\hat H_p^{\mathrm{LS}}|^2 + 10^{-\mathrm{snr\_db}/10}}
    \,\hat H_p^{\mathrm{LS}}.
```

`snr_db` is used only as this shrinkage parameter. The formula is a per-pilot heuristic: it treats ``|\hat H|^2`` as if it were the prior power and it does **not** divide the noise variance by the FFT length or by the pilot energy. It is not the linear MMSE estimator that uses the channel frequency-correlation matrix (Edfors et al., 1998). At high `snr_db` the shrinkage tends to 1 and the receiver matches LS.
"""
function ofdm_rx_mmse(rx_signal::AbstractVector{<:Complex},
                      params::OFDMParams,
                      pilot_syms::AbstractVector{<:Complex},
                      snr_db::Real)
    freq_syms = extract_ofdm_symbols(rx_signal, params)
    scheme = modulation(params)
    noise_var = exp10(-float(snr_db) / 10)
    rx_bits = Int[]
    for fd in freq_syms
        H_ls = fd[params.pilot_indices] ./ pilot_syms
        power = abs2.(H_ls)
        H_shrunk = power ./ (power .+ noise_var) .* H_ls
        H_est = interpolate_channel(H_shrunk, params)
        data_eq = fd[params.data_indices] ./ H_est[params.data_indices]
        append!(rx_bits, demodulate(data_eq, scheme))
    end
    return rx_bits
end

@doc raw"""
    sto_estimate(rx_signal, params::OFDMParams; search_range=N_fft+N_cp) -> Int

Coarse symbol timing from the cyclic-prefix correlation (van de Beek, Sandell, and Börjesson, IEEE Trans. Signal Process., 45(7):1800–1805, 1997, the magnitude form of their correlation metric),

```math
\Lambda(t) = \frac{\bigl|\sum_{i=0}^{N_{cp}-1} r_{t+i} r_{t+i+N}^{*}\bigr|^2}
{\bigl(\sum_{i=0}^{N_{cp}-1} |r_{t+i}|^2 + |r_{t+i+N}|^2\bigr)^2}.
```

# Arguments
- `search_range`: number of candidate starts. The default is one OFDM symbol, so the window contains a single burst start when the delay is shorter than one symbol.

# Returns
The 1-based index of the maximiser. A noiseless burst has an identical peak on every symbol; when several maximisers agree up to roundoff (`1e-8` relative), the **earliest** one is returned. In noise the global maximiser inside the window is kept, which need not be the first symbol if `search_range` covers more than one symbol.
"""
function sto_estimate(rx_signal::AbstractVector{<:Complex}, params::OFDMParams;
                      search_range::Int=params.N_fft + params.N_cp)
    N = params.N_fft
    N_cp = params.N_cp
    sym = N + N_cp
    search_range >= 1 || throw(ArgumentError("search_range must be ≥ 1"))
    L = min(length(rx_signal), search_range + sym)
    nmetric = L - sym
    nmetric >= 1 || throw(ArgumentError("signal is shorter than one OFDM symbol plus the search window"))
    rx = ComplexF64.(rx_signal)
    metric = Vector{Float64}(undef, nmetric)
    @inbounds for t in 1:nmetric
        s = 0.0 + 0.0im
        e = 0.0
        for i in 0:(N_cp - 1)
            a = rx[t + i]
            b = rx[t + i + N]
            s += a * conj(b)
            e += abs2(a) + abs2(b)
        end
        metric[t] = abs2(s) / (e^2 + 1e-10)
    end
    peak, idx = findmax(metric)
    while idx > sym && metric[idx - sym] >= peak * (1 - 1e-8)
        idx -= sym
    end
    return idx
end

@doc raw"""
    cfo_estimate(rx_signal, params::OFDMParams; sto=sto_estimate(rx_signal, params)) -> Float64

Carrier-frequency offset, in units of the subcarrier spacing, from the cyclic-prefix correlation at a known start `sto` (van de Beek et al., 1997),

```math
\hat\varepsilon = -\frac{1}{2\pi}\arg\sum_{i=0}^{N_{cp}-1} r_{t+i}\, r_{t+i+N}^{*}.
```

The estimate is unambiguous for ``\varepsilon \in (-1/2, 1/2)``. `sto` is 1-based.

This is the van de Beek cyclic-prefix estimator. It does not estimate an integer number of subcarrier spacings: an offset of ``0.7`` is reported near ``-0.3``. A repeated-half preamble estimator with fractional range ``(-1, 1)`` and a separate even-integer search lives in `Sync` (`schmidl_cox`, `integer_cfo_estimate`). Those functions do not call this one, and this function does not look for that preamble.
"""
function cfo_estimate(rx_signal::AbstractVector{<:Complex}, params::OFDMParams;
                      sto::Int=sto_estimate(rx_signal, params))
    N = params.N_fft
    N_cp = params.N_cp
    sto >= 1 || throw(ArgumentError("sto must be ≥ 1, got $sto"))
    last = sto + N + N_cp - 1
    last <= length(rx_signal) || throw(ArgumentError(
        "signal length $(length(rx_signal)) is too short for CFO estimation at sto=$sto",
    ))
    rx = ComplexF64.(rx_signal)
    corr = zero(ComplexF64)
    @inbounds for i in 0:(N_cp - 1)
        corr += rx[sto + i] * conj(rx[sto + N + i])
    end
    return -angle(corr) / (2π)
end

@doc raw"""
    cfo_correct(rx_signal, cfo_norm::Real, n_fft::Integer) -> Vector{ComplexF64}

Multiply sample ``i`` (1-based) by ``e^{-j 2\pi \varepsilon (i-1) / N}``, where ``\varepsilon`` is `cfo_norm` in subcarrier spacings and ``N`` is `n_fft`. This inverts a residual frequency ``\varepsilon / N`` cycles per sample. Symbol timing is not modified.
"""
function cfo_correct(rx_signal::AbstractVector{<:Complex}, cfo_norm::Real, n_fft::Integer)
    n_fft >= 1 || throw(ArgumentError("n_fft must be ≥ 1"))
    ε = float(cfo_norm)
    N = float(n_fft)
    return ComplexF64[rx_signal[i] * exp(-im * 2π * ε * (i - 1) / N) for i in eachindex(rx_signal)]
end

@doc raw"""
    ofdm_sync(rx_signal, params::OFDMParams; search_range=N_fft+N_cp) -> (corrected, sto, cfo)

Run [`sto_estimate`](@ref) and [`cfo_estimate`](@ref), then [`cfo_correct`](@ref) the whole vector.

# Returns
The CFO-corrected waveform (still including any leading delay), the 1-based timing index, and the normalised CFO. The returned waveform is **not** shifted by `sto`; a caller that needs the FFT window to start at sample 1 has to discard `sto - 1` samples.

The CFO half of this function is [`cfo_estimate`](@ref): capture ``(-1/2, 1/2)``, no integer offset. `Sync.schmidl_sync` is a different estimator that needs a repeated-half preamble and also returns an unshifted buffer. Calling `ofdm_sync` does not run that estimator, and `schmidl_sync` does not call `ofdm_sync`.
"""
function ofdm_sync(rx_signal::AbstractVector{<:Complex}, params::OFDMParams;
                   search_range::Int=params.N_fft + params.N_cp)
    sto = sto_estimate(rx_signal, params; search_range=search_range)
    cfo = cfo_estimate(rx_signal, params; sto=sto)
    return cfo_correct(rx_signal, cfo, params.N_fft), sto, cfo
end

end # module OFDM
