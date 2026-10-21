import CommToolkit: modulate, demodulate, ber, awgn, QPSK
import Plots
using SpecialFunctions

n_bits = 10^6
snr_db_range = 0:2:14

ber_sim = zeros(length(snr_db_range))

for (i, snr) in enumerate(snr_db_range)
    tx_bits = rand(0:1, n_bits)
    syms = modulate(tx_bits, QPSK())
    rx_syms = awgn(syms, snr, 2)  # 2 bits per symbol
    rx_bits = demodulate(rx_syms, QPSK())
    ber_sim[i] = ber(tx_bits, rx_bits)
end

# Theoretical BER for QPSK ≈ BPSK (same Eb/N0)
ber_theory = 0.5 * erfc.(sqrt.(10 .^(snr_db_range ./ 10)))

p = Plots.plot(snr_db_range, ber_sim,
    yaxis=:log, marker=:circle,
    label="Simulated", legend=:bottomleft)
Plots.plot!(p, snr_db_range, ber_theory,
    yaxis=:log, line=:dash,
    label="Theory")
Plots.xlabel!("Eb/N0 (dB)")
Plots.ylabel!("Bit Error Rate")
Plots.title!("QPSK over AWGN")

Plots.savefig(joinpath(@__DIR__, "ber_qpsk.png"))
Plots.display(p)
println("QPSK BER simulation complete.")
