# Teaching link: punctured convolutional MCS, Chase of one punctured word,
# and one BICM receiver. Not an NR/LTE scheduler, not incremental redundancy.
#
# Es = 1, σ² = N0/2, seed 0x0c0de1162.
# Throughput is delivered payload bits per transmitted symbol.
# The MCS thresholds below were read off a Julia 1.10.12 sweep of this script's
# fixed-MCS columns (80 blocks, integer Es/N0 from -2 dB to 16 dB). The 5 dB
# bin is a near-tie between QPSK 2/3 and QPSK 3/4; do not treat it as stable.
# Run from the repo root:
#   julia --project=. demos/link_puncture_bicm.jl

using Random
using Printf
using CommToolkit

const SEED = 0x0c0de1162
const N_BLOCKS = 80

table = teaching_punctured_mcs_table()
policy = ThresholdAMC(teaching_punctured_amc_thresholds())
snrs = [0.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0]

println("CommToolkit punctured MCS link (not NR, not IR-HARQ).")
println("Seed ", string(SEED, base=16), ", ", N_BLOCKS,
        " blocks, receiver = Viterbi unless noted.")
println("Payload 54 bits, mother length 144, CRC-16 and a 2-bit tail are overhead.")
println("Spectral efficiency (payload bits / transmitted symbol):")
for row in table
    @printf("  %-12s  rate %-4s  %6.3f\n", row.name, string(row.rate),
            spectral_efficiency(row))
end
println()
println("Fixed columns follow teaching_punctured_mcs_table. adaptive uses")
println("teaching_punctured_amc_thresholds and ideal Es/N0. Chase×3 repeats the")
println("QPSK-3/4 word and sums bit LLRs; deleted mother bits stay at LLR 0.")
println()

header = "Es/N0  " * join([@sprintf("%8s", string(i)) for i in 1:length(table)]) *
         "   adaptive    QPSK-3/4 Chase×3"
println(header)
println("dB      " * " fixed throughput (bit/use)")

for snr in snrs
    cells = String[]
    for i in eachindex(table)
        Random.seed!(SEED)
        report = simulate_link(snr, N_BLOCKS; table=table, mcs=i, max_tx=1)
        push!(cells, @sprintf("%8.3f", report.throughput))
    end
    Random.seed!(SEED)
    adapt = simulate_link(snr, N_BLOCKS; table=table, policy=policy, max_tx=1)
    Random.seed!(SEED)
    chase = simulate_link(snr, N_BLOCKS; table=table, mcs=4, max_tx=3)
    @printf("%5.1f  %s   %6.3f/%4.2f   %6.3f/%4.2f\n",
            snr, join(cells), adapt.throughput, adapt.bler,
            chase.throughput, chase.bler)
end

println()
println("BICM on one 8-PSK row, set-partitioning labels, 12-row interleaver,")
println("max-log, max_tx = 1. Iterations 4 against one demapper pass, same seed.")
println("This is not the Gray Viterbi receiver, and it is not a capacity claim.")

psk = [mcs_convolutional("8PSK-1/2", PSK(8), 48)]
labels = bicm_psk8_id_labels()
for iterations in (1, 4)
    Random.seed!(SEED)
    report = simulate_link(6.0, N_BLOCKS; table=psk, mcs=1, receiver=:bicm,
                          bicm_iterations=iterations, bicm_n_rows=12,
                          bicm_labels=labels, approx=:maxlog)
    @printf("  iterations %d   BER %.4f   BLER %.3f   T %.3f bit/use   bit errors %d/%d\n",
            iterations, report.ber, report.bler, report.throughput,
            report.n_bit_errors, report.n_payload_bits)
end
