import CommToolkit: modulate, demodulate, ber, awgn, rayleigh_channel
import CommToolkit: BPSK, QPSK, QAM, snr_db_to_linear
import Plots

n_bits = 200000
snr_range = 0:2:14

ber_bpsk_awgn = Float64[]
ber_qpsk_awgn = Float64[]
ber_16qam_awgn = Float64[]
ber_bpsk_rayleigh = Float64[]

for snr in snr_range
    tx = rand(0:1, n_bits)

    # BPSK AWGN
    rx = awgn(modulate(tx, BPSK()), snr, 1)
    push!(ber_bpsk_awgn, ber(tx, demodulate(rx, BPSK())))

    # QPSK AWGN
    rx = awgn(modulate(tx, QPSK()), snr, 2)
    push!(ber_qpsk_awgn, ber(tx, demodulate(rx, QPSK())))

    # 16QAM AWGN
    n_qam = n_bits ÷ 4 * 4
    tx_qam = tx[1:n_qam]
    rx = awgn(modulate(tx_qam, QAM(16)), snr, 4)
    push!(ber_16qam_awgn, ber(tx_qam, demodulate(rx, QAM(16))))

    # BPSK Rayleigh (perfect CSIR)
    rx, _ = rayleigh_channel(modulate(tx, BPSK()), snr, 1)
    rx_eq = rayleigh_channel(modulate(tx, BPSK()), snr, 1, true)
    push!(ber_bpsk_rayleigh, ber(tx, demodulate(rx_eq, BPSK())))
end

p = Plots.plot(snr_range, ber_bpsk_awgn,
    yaxis=:log, marker=:circle, lw=2,
    label="BPSK AWGN", legend=:bottomleft)
Plots.plot!(p, snr_range, ber_qpsk_awgn,
    yaxis=:log, marker=:square, lw=2,
    label="QPSK AWGN")
Plots.plot!(p, snr_range, ber_16qam_awgn,
    yaxis=:log, marker=:diamond, lw=2,
    label="16QAM AWGN")
Plots.plot!(p, snr_range, ber_bpsk_rayleigh,
    yaxis=:log, marker=:x, lw=2, ls=:dash,
    label="BPSK Rayleigh", color=:red)
Plots.xlabel!("Eb/N0 (dB)")
Plots.ylabel!("Bit Error Rate")
Plots.title!("BER Comparison: Modulation Schemes over AWGN and Rayleigh")

Plots.savefig(joinpath(@__DIR__, "ber_comparison.png"))
Plots.display(p)
println("Comparison BER simulation complete.")
