import CommToolkit: modulate, demodulate, ber, awgn, QAM
import Plots
using SpecialFunctions

n_bits = 2 * 10^6
snr_db_range = 0:2:16

ber_sim = zeros(length(snr_db_range))

for (i, snr) in enumerate(snr_db_range)
    tx_bits = rand(0:1, n_bits)
    syms = modulate(tx_bits, QAM(16))
    rx_syms = awgn(syms, snr, 4)  # 4 bits per symbol
    rx_bits = demodulate(rx_syms, QAM(16))
    ber_sim[i] = ber(tx_bits, rx_bits)
end

# Approximate theoretical BER for 16QAM:
# Pb ≈ (3/4) * erfc(sqrt(4*Eb/(5*N0)))  ... simplified
# More precisely using 3/(2k)*Q(sqrt(3k*SNR/(M-1))) for square QAM
k = 4  # bits per symbol
ber_theory = (3/(2*k)) * erfc.(sqrt.(3 * k * 10 .^(snr_db_range ./ 10) / (2 * (16 - 1))))

p = Plots.plot(snr_db_range, ber_sim,
    yaxis=:log, marker=:circle,
    label="Simulated", legend=:bottomleft)
Plots.plot!(p, snr_db_range, ber_theory,
    yaxis=:log, line=:dash,
    label="Theory (approx)")
Plots.xlabel!("Eb/N0 (dB)")
Plots.ylabel!("Bit Error Rate")
Plots.title!("16QAM over AWGN")

Plots.savefig(joinpath(@__DIR__, "ber_16qam.png"))
Plots.display(p)
println("16QAM BER simulation complete.")
