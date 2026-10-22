# Hard-decision BER against elementary closed forms.
# DBPSK: (1/2) exp(−Eb/N0). Rayleigh and MRC: Proakis and Salehi, 5th ed., §13.3–§13.4.
# The Gray M-PSK sector integral and the DQPSK Marcum-Q expression are in test/reference.jl.
# No Plots.

using CommToolkit
using Printf

function mc_ber(nbits, simulate)
    bits = rand(0:1, nbits)
    return count(!=(0), bits .- simulate(bits)) / nbits
end

function ber_mrc(L, γ)
    μ = sqrt(γ / (1 + γ))
    α = (1 - μ) / 2
    β = (1 + μ) / 2
    s = sum(binomial(L - 1 + k, k) * β^k for k in 0:(L - 1); init=0.0)
    return α^L * s
end

nbits = 100_000
snr_db = 6.0
γ = 10^(snr_db / 10)

cases = [
    "DBPSK" => (mc_ber(nbits, bits -> demodulate(awgn(modulate(bits, DBPSK()), snr_db, 1), DBPSK())),
                0.5 * exp(-γ)),
    "Rice K=0 BPSK" => (mc_ber(nbits, bits -> demodulate(
                    rice_channel(modulate(bits, BPSK()), snr_db, 1, true; K=0), BPSK())),
                0.5 * (1 - sqrt(γ / (1 + γ)))),
    "MRC L=2 Rayleigh BPSK" => (mc_ber(nbits, bits -> begin
                    z, _ = mrc_channel(modulate(bits, BPSK()), snr_db, 1, 2; K=0)
                    demodulate(z, BPSK())
                end),
                ber_mrc(2, γ)),
]

println("Eb/N0 = ", snr_db, " dB, bits = ", nbits, " (stream not seeded)")
for (name, (sim, formula)) in cases
    @printf("%-24s  sim %.4e   formula %.4e\n", name, sim, formula)
end

bits = rand(0:1, 63)
@assert demodulate(modulate(bits, PSK(8)), PSK(8)) == bits
dbits = rand(0:1, 64)
@assert demodulate(modulate(dbits, DBPSK()) .* cis(0.2), DBPSK()) == dbits
@assert demodulate(modulate(dbits, DQPSK()) .* cis(-0.4), DQPSK()) == dbits
