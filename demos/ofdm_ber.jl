import CommToolkit: OFDMParams, TDLChannel, ofdm_tx, apply_channel
import CommToolkit: ofdm_rx_known, ofdm_rx_ls, ofdm_rx_mmse, ber
import Plots

params = OFDMParams(N_fft=64, N_cp=16, mod_order=4)
channel = TDLChannel(delays_samples=[0, 1, 3], gain_dB=[0, -2, -5])
snr_range = 0:2:18
n_syms = 200

ber_known = Float64[]
ber_ls = Float64[]
ber_mmse = Float64[]

for snr in snr_range
    bits = rand(0:1, params.N_data * params.bits_per_symbol * n_syms)
    tx, pilots = ofdm_tx(params, bits)
    rx, H = apply_channel(tx, channel, snr, params)

    bits_k = ofdm_rx_known(rx, params, H, pilots)
    push!(ber_known, ber(bits, bits_k[1:length(bits)]))

    bits_ls = ofdm_rx_ls(rx, params, pilots)
    push!(ber_ls, ber(bits, bits_ls[1:length(bits)]))

    bits_m = ofdm_rx_mmse(rx, params, pilots, snr)
    push!(ber_mmse, ber(bits, bits_m[1:length(bits)]))

    println("SNR=$(snr)dB: known=$(round(ber_known[end], sigdigits=3)), " *
            "LS=$(round(ber_ls[end], sigdigits=3)), " *
            "MMSE=$(round(ber_mmse[end], sigdigits=3))")
end

p = Plots.plot(snr_range, ber_known,
    yaxis=:log, marker=:circle, lw=2,
    label="Known channel", legend=:bottomleft)
Plots.plot!(p, snr_range, ber_ls,
    yaxis=:log, marker=:square, lw=2,
    label="LS estimation")
Plots.plot!(p, snr_range, ber_mmse,
    yaxis=:log, marker=:diamond, lw=2,
    label="MMSE estimation")
Plots.xlabel!("SNR (dB)")
Plots.ylabel!("Bit Error Rate")
Plots.title!("OFDM BER: QPSK over 3-tap TDL Channel")

Plots.savefig(joinpath(@__DIR__, "ofdm_ber.png"))
Plots.display(p)
println("OFDM BER simulation complete.")
