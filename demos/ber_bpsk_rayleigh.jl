import CommToolkit: modulate, demodulate, ber, rayleigh_channel, BPSK
import Plots

n_bits = 2 * 10^5
snr_db_range = 0:5:40

ber_sim = zeros(length(snr_db_range))

for (i, snr) in enumerate(snr_db_range)
    tx_bits = rand(0:1, n_bits)
    syms = modulate(tx_bits, BPSK())
    rx = rayleigh_channel(syms, snr, 1, true)  # perfect CSI
    rx_bits = demodulate(rx, BPSK())
    ber_sim[i] = ber(tx_bits, rx_bits)
end

# Theoretical BER for BPSK over flat Rayleigh with perfect CSIR:
# Pb = 0.5 * (1 - sqrt(SNR/(1+SNR)))
snr_lin = 10 .^(snr_db_range ./ 10)
ber_theory = 0.5 .* (1 .- sqrt.(snr_lin ./ (1 .+ snr_lin)))

p = Plots.plot(snr_db_range, ber_sim,
    yaxis=:log, marker=:circle,
    label="Simulated (perfect CSIR)", legend=:bottomleft)
Plots.plot!(p, snr_db_range, ber_theory,
    yaxis=:log, line=:dash,
    label="Theory")
Plots.xlabel!("Eb/N0 (dB)")
Plots.ylabel!("Bit Error Rate")
Plots.title!("BPSK over Flat Rayleigh Fading (perfect CSIR)")

Plots.savefig(joinpath(@__DIR__, "ber_bpsk_rayleigh.png"))
Plots.display(p)
println("Rayleigh BER simulation complete.")
