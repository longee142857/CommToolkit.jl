# OFDM receiver: known channel, LS, MMSE

using FFTW

function extract_ofdm_symbols(rx_signal::Vector{ComplexF64}, params::OFDMParams)
    N = params.N_fft
    N_cp = params.N_cp
    sym_len = N + N_cp
    n_syms = length(rx_signal) ÷ sym_len

    map(1:n_syms) do i
        idx = (i-1) * sym_len + N_cp + 1
        fft(rx_signal[idx:idx+N-1])
    end
end

# Known channel receiver (benchmark upper bound)
function ofdm_rx_known(rx_signal::Vector{ComplexF64},
                        params::OFDMParams,
                        H_freq::Vector{ComplexF64},
                        pilot_syms::Vector{ComplexF64})
    freq_syms = extract_ofdm_symbols(rx_signal, params)
    scheme = modulation(params)

    rx_bits = Int[]
    for fd in freq_syms
        data_rx = fd[params.data_indices]
        H_data = H_freq[params.data_indices]
        data_eq = data_rx ./ H_data
        append!(rx_bits, demodulate(data_eq, scheme))
    end
    rx_bits
end

# LS channel estimation receiver
function ofdm_rx_ls(rx_signal::Vector{ComplexF64},
                     params::OFDMParams,
                     pilot_syms::Vector{ComplexF64})
    freq_syms = extract_ofdm_symbols(rx_signal, params)
    scheme = modulation(params)

    rx_bits = Int[]
    for fd in freq_syms
        H_pilot = fd[params.pilot_indices] ./ pilot_syms
        H_est = interpolate_channel(H_pilot, params)

        data_rx = fd[params.data_indices]
        data_eq = data_rx ./ H_est[params.data_indices]
        append!(rx_bits, demodulate(data_eq, scheme))
    end
    rx_bits
end

# MMSE channel estimation receiver
function ofdm_rx_mmse(rx_signal::Vector{ComplexF64},
                       params::OFDMParams,
                       pilot_syms::Vector{ComplexF64},
                       snr_db::Real)
    freq_syms = extract_ofdm_symbols(rx_signal, params)
    scheme = modulation(params)
    snr_lin = 10^(snr_db / 10)

    rx_bits = Int[]
    for fd in freq_syms
        H_pilot_ls = fd[params.pilot_indices] ./ pilot_syms

        # Simplified MMSE: Wiener-like scaling
        noise_var = 1.0 / snr_lin
        H_mmse = abs2.(H_pilot_ls) ./ (abs2.(H_pilot_ls) .+ noise_var) .* H_pilot_ls
        H_est = interpolate_channel(H_mmse, params)

        data_rx = fd[params.data_indices]
        data_eq = data_rx ./ H_est[params.data_indices]
        append!(rx_bits, demodulate(data_eq, scheme))
    end
    rx_bits
end

# Linear interpolation of channel estimates
function interpolate_channel(H_pilot::Vector{ComplexF64}, params::OFDMParams)
    N = params.N_fft
    H_est = zeros(ComplexF64, N)
    pilot_idx = params.pilot_indices
    data_idx = params.data_indices

    for (i, idx) in enumerate(pilot_idx)
        H_est[idx] = H_pilot[i]
    end

    for d_idx in data_idx
        left = findlast(x -> x < d_idx, pilot_idx)
        right = findfirst(x -> x > d_idx, pilot_idx)

        if left === nothing
            H_est[d_idx] = H_pilot[right]
        elseif right === nothing
            H_est[d_idx] = H_pilot[left]
        else
            α = (d_idx - pilot_idx[left]) / (pilot_idx[right] - pilot_idx[left])
            H_est[d_idx] = H_pilot[left] * (1 - α) + H_pilot[right] * α
        end
    end
    H_est
end
