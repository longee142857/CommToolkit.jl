# Teaching link: fixed MCS, threshold AMC, and Chase combining on AWGN.
#
# This is not an NR/LTE scheduler. Es = 1, σ² = N0/2, seed 0x0c0de1162.
# Throughput is delivered payload bits per transmitted symbol.
# Run from the repo root:
#   julia --project=. demos/link_mcs_harq.jl

using Random
using Printf
using CommToolkit

const SEED = 0x0c0de1162
const N_BLOCKS = 80

table = teaching_mcs_table()
policy = ThresholdAMC(teaching_amc_thresholds())
snrs = [0.0, 2.0, 5.0, 8.0, 12.0, 16.0]

println("CommToolkit teaching link (not NR). Seed ", string(SEED, base=16),
        ", ", N_BLOCKS, " blocks, max_iter LDPC = 20.")
println("MCS spectral efficiency (payload bits / symbol):")
for row in table
    @printf("  %-18s  %6.3f\n", row.name, spectral_efficiency(row))
end
println()
header = "Es/N0  " * join([@sprintf("%8s", string(i)) for i in 1:length(table)]) *
         "   adaptive      16QAM Chase×3"
println(header)
println("dB      " * " fixed throughput (bit/use)" * " "^18 * "T / BLER")

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
println("Columns 1–5 are fixed MCS rows in teaching_mcs_table order, one transmission.")
println("adaptive uses teaching_amc_thresholds and ideal Es/N0.")
println("16QAM Chase×3 repeats the same coded block and sums bit LLRs.")
