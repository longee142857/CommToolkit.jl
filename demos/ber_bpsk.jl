# BPSK BER over AWGN
# Compares simulated BER with theoretical Eb/N0 curve

import CommToolkit: modulate, demodulate, ber, awgn, BPSK
import Plots
using SpecialFunctions

# Simulation parameters
n_bits = 10^6
snr_db_range = 0:2:12

ber_sim = zeros(length(snr_db_range))

for (i, snr) in enumerate(snr_db_range)
    tx_bits = rand(0:1, n_bits)
    syms = modulate(tx_bits, BPSK())
    rx_syms = awgn(syms, snr, 1)
    rx_bits = demodulate(rx_syms, BPSK())
    ber_sim[i] = ber(tx_bits, rx_bits)
end

# Theoretical BER for BPSK: 0.5 * erfc(sqrt(Eb/N0))
ber_theory = 0.5 * erfc.(sqrt.(10 .^(snr_db_range ./ 10)))

# Plot
p = Plots.plot(snr_db_range, ber_sim,
    yaxis=:log, marker=:circle,
    label="Simulated", legend=:bottomleft)
Plots.plot!(p, snr_db_range, ber_theory,
    yaxis=:log, line=:dash,
    label="Theory")
Plots.xlabel!("Eb/N0 (dB)")
Plots.ylabel!("Bit Error Rate")
Plots.title!("BPSK over AWGN")

Plots.savefig(joinpath(@__DIR__, "ber_bpsk.png"))
Plots.display(p)
println("BER simulation complete. Plot saved to demos/ber_bpsk.png")
