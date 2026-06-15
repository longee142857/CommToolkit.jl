# OFDM transmitter

using FFTW

function _generate_pilots(n::Int)
    [_lcg_rand(42 + i) > 0.5 ? 1.0+0.0im : -1.0+0.0im for i in 1:n]
end

function generate_ofdm_symbol(params::OFDMParams,
                              data_syms::Vector{ComplexF64},
                              pilot_syms::Vector{ComplexF64})
    N = params.N_fft
    freq_domain = zeros(ComplexF64, N)
    freq_domain[params.data_indices] .= data_syms
    freq_domain[params.pilot_indices] .= pilot_syms

    time_domain = ifft(freq_domain)

    cp = time_domain[end-params.N_cp+1:end]
    vcat(cp, time_domain)
end

function ofdm_tx(params::OFDMParams, bits::Vector{Int})
    bits_per_ofdm_sym = params.N_data * params.bits_per_symbol
    n_ofdm_syms = ceil(Int, length(bits) / bits_per_ofdm_sym)
    total_bits_needed = n_ofdm_syms * bits_per_ofdm_sym
    padded_bits = vcat(bits, zeros(Int, total_bits_needed - length(bits)))

    pilot_syms = _generate_pilots(params.N_pilot)
    scheme = modulation(params)

    tx_signal = ComplexF64[]
    for sym_idx in 1:n_ofdm_syms
        idx = (sym_idx - 1) * bits_per_ofdm_sym
        data_bits = padded_bits[idx+1:idx+bits_per_ofdm_sym]
        data_syms = modulate(data_bits, scheme)
        ofdm_sym = generate_ofdm_symbol(params, data_syms, pilot_syms)
        append!(tx_signal, ofdm_sym)
    end

    tx_signal, pilot_syms
end
