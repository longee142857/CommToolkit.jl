# OFDM channel models: TDL multipath + AWGN

using FFTW

struct TDLChannel
    taps::Vector{ComplexF64}
    delays::Vector{Int}
end

function TDLChannel(; delays_samples::Vector{Int}=[0, 1, 3],
                      gain_dB::Vector{<:Real}=[0.0, -2.0, -5.0])
    gains_linear = [10^(g/20) for g in gain_dB]
    taps = [g * cispi(2 * _lcg_rand(123 + i)) for (i, g) in enumerate(gains_linear)]
    TDLChannel(taps, delays_samples)
end

function apply_multipath(signal::Vector{ComplexF64}, channel::TDLChannel)
    L = length(signal)
    max_delay = maximum(channel.delays)
    output = zeros(ComplexF64, L + max_delay)

    for (tap, delay) in zip(channel.taps, channel.delays)
        for i in 1:L
            output[i + delay] += tap * signal[i]
        end
    end
    output
end

function _channel_freq_response(channel::TDLChannel, n::Int)
    H = zeros(ComplexF64, n)
    for k in 1:n
        f = (k - 1) / n
        for (tap, delay) in zip(channel.taps, channel.delays)
            H[k] += tap * exp(-1im * 2π * f * delay)
        end
    end
    H
end

# Apply channel: multipath + AWGN.
# Returns (rx_signal, channel_freq_response)
function apply_channel(signal::Vector{ComplexF64}, channel::TDLChannel,
                       snr_db::Real, params::OFDMParams)
    N = params.N_fft
    N_cp = params.N_cp

    faded = apply_multipath(signal, channel)
    sym_len = N + N_cp
    total_syms = length(signal) ÷ sym_len
    expected = total_syms * sym_len

    rx = if length(faded) > expected
        faded[1:expected]
    elseif length(faded) < expected
        vcat(faded, zeros(ComplexF64, expected - length(faded)))
    else
        faded
    end

    # AWGN using signal power based SNR (OFDM convention)
    rx = _awgn_ofdm(rx, snr_db)

    H = _channel_freq_response(channel, N)
    rx, H
end

function _awgn_ofdm(signal::Vector{ComplexF64}, snr_db::Real)
    P = sum(abs2, signal) / length(signal)
    snr_lin = 10^(snr_db / 10)
    noise_power = P / snr_lin
    noise_std = sqrt(noise_power / 2)
    signal + noise_std * randn(length(signal)) + 1im * noise_std * randn(length(signal))
end
