# Teaching multi-user KPI→action mini xApp: three UEs on independent single
# links, sharing one transport-block budget per round and one share-weighted
# mean power budget. A rule table maps cumulative link KPIs onto discrete
# actions; a fixed allocation is the baseline it is compared against.
#
# This is a teaching mini xApp, not O-RAN: no E2 interface, no near-RT RIC, no
# fronthaul or functional split, no NR scheduler, no NR CQI, and no O-RAN
# compliance claim. Multi-user here means independent single links plus the
# two shared resource constraints — no inter-UE interference model, no cell
# geometry. The PHY is the existing single-link teaching loop, untouched.
#
# Es = 1, σ² = N0/2, seed 0x0c0de1162.
# Run from the repo root:
#   julia --project=. demos/xapp_multiuser.jl

using Random
using Printf
using CommToolkit

const SEED = 0x0c0de1162
const N_ROUNDS = 20
const BLOCKS_PER_ROUND = 30

ues = [XAppUE("UE-0dB", 0.0), XAppUE("UE-3dB", 3.0), XAppUE("UE-6dB", 6.0)]

println("CommToolkit teaching mini xApp (KPI→action, not O-RAN).")
println("Seed ", string(SEED, base=16), ", ", N_ROUNDS, " rounds × ",
        BLOCKS_PER_ROUND, " blocks, AWGN, max_tx = 1, Chase.")
println("UEs: ", join(("$(u.name) at $(u.esn0_db) dB" for u in ues), ", "), ".")
println()
println("KPIs observed per UE (cumulative, no sliding window): BLER and")
println("delivered payload bits. Action space per UE per round: an MCS row of")
println("teaching_punctured_mcs_table, a share of the round's block budget, and")
println("a power offset on the discrete ladder [-3, 0, +3] dB. Constraint:")
println("share-weighted mean linear power ≤ 1 (the Es = 1 baseline).")
println()
println("Rule table (RuleXAppPolicy): share ∝ 1/(delivered bits + 1), offset")
println("steps one rung up above BLER 0.1 and one rung down below 0.01, MCS")
println("from the teaching thresholds at base Es/N0 + offset. Baseline")
println("(FixedXAppPolicy): QPSK-1/2 for everyone, uniform shares, 0 dB.")
println()

function print_report(label, report)
    println(label)
    println("  UE        blocks  delivered    BLER   T (bit/use)")
    for kpi in report.kpis
        @printf("  %-8s %6d  %9d  %7.3f   %9.3f\n",
                kpi.name, kpi.n_blocks, kpi.n_delivered, kpi.bler,
                kpi.throughput)
    end
    @printf("  system: throughput %.3f bit/use, BLER %.3f, Jain %.3f\n",
            report.system_throughput, report.system_bler, report.jain)
    println("  last actions: ", join((string(a) for a in report.last_actions), "  "))
    println()
    return nothing
end

Random.seed!(SEED)
baseline = simulate_xapp(ues, FixedXAppPolicy(mcs=2), N_ROUNDS;
                         blocks_per_round=BLOCKS_PER_ROUND)
print_report("Baseline (fixed QPSK-1/2, uniform shares):", baseline)

Random.seed!(SEED)
ruled = simulate_xapp(ues, RuleXAppPolicy(length(ues)), N_ROUNDS;
                      blocks_per_round=BLOCKS_PER_ROUND)
print_report("KPI→action rule table:", ruled)

println("Same seed for both arms: the rule table separates the delivered-block")
println("Clopper–Pearson intervals at α = 1e-4 in test/xapp.jl; the numbers")
println("above are point estimates on this one stream, not a capacity claim.")
println()
println("Honest corner: the weak UE's cumulative BLER keeps early failures")
println("forever (no forgetting factor), so its offset rung can stay up after")
println("the link has recovered. That is the documented teaching simplification,")
println("not a sliding-window outer loop. Full geometry CDL, NR LDPC, and a real")
println("O-RAN E2SM remain outside the toolkit.")
