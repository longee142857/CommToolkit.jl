# OFDM synchronization: STO and CFO estimation/compensation

# CP-based STO estimation
# Returns the estimated start index of the first OFDM symbol
function sto_estimate(rx_signal::Vector{ComplexF64}, params::OFDMParams;
                      search_range::Int=params.N_fft + params.N_cp)
    N = params.N_fft
    N_cp = params.N_cp
    L = min(length(rx_signal), search_range + N + N_cp)

    metric = zeros(Float64, L - N - N_cp)
    for t in 1:length(metric)
        s = 0.0
        e = 0.0
        for i in 0:N_cp-1
            s += rx_signal[t+i] * conj(rx_signal[t+i+N])
            e += abs2(rx_signal[t+i]) + abs2(rx_signal[t+i+N])
        end
        metric[t] = abs(s)^2 / (e^2 + 1e-10)
    end

    _, est_start = findmax(metric)
    est_start
end

# CFO estimation from CP correlation
# Returns normalized CFO (relative to subcarrier spacing)
function cfo_estimate(rx_signal::Vector{ComplexF64}, params::OFDMParams;
                      sto::Int=sto_estimate(rx_signal, params))
    N = params.N_fft
    N_cp = params.N_cp
    corr = sum(rx_signal[sto:sto+N_cp-1] .* conj(rx_signal[sto+N:sto+N+N_cp-1]))
    -angle(corr) / (2π)
end

# CFO correction
function cfo_correct(rx_signal::Vector{ComplexF64}, cfo_norm::Real, n_fft::Int)
    [rx_signal[i] * exp(-1im * 2π * cfo_norm * (i-1) / n_fft) for i in 1:length(rx_signal)]
end

# Combined sync: STO + CFO estimation + CFO correction
function ofdm_sync(rx_signal::Vector{ComplexF64}, params::OFDMParams;
                   search_range::Int=params.N_fft + params.N_cp)
    coarse_sto = sto_estimate(rx_signal, params, search_range=search_range)
    cfo = cfo_estimate(rx_signal, params, sto=coarse_sto)
    corrected = cfo_correct(rx_signal, cfo, params.N_fft)
    corrected, coarse_sto, cfo
end
