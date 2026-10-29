# Teaching multi-user KPI→action mini xApp. Not O-RAN: no E2 interface, no
# near-RT RIC, no NR scheduler, no NR CQI. Policies are checked
# deterministically on synthetic KPIs and on noiseless links; the closed-loop
# comparison against the fixed baseline uses seed 0x0c0de1162 with
# Clopper–Pearson at BER_ALPHA on delivered-block counts. Counts are not
# locked: a different randn stream must still separate the intervals below.

@testset "XAppUE / XAppAction / UEKPI" begin
    ue = XAppUE("ue1", 3.0)
    @test ue.name == "ue1"
    @test ue.esn0_db == 3.0
    @test XAppUE("inf", Inf).esn0_db == Inf
    @test occursin("ue1", sprint(show, ue))
    @test_throws ArgumentError XAppUE("", 0.0)
    @test_throws ArgumentError XAppUE("bad", NaN)
    @test_throws ArgumentError XAppUE("bad", -Inf)

    action = XAppAction(2, 1 / 3, -3.0)
    @test action.mcs == 2
    @test action.share ≈ 1 / 3
    @test action.power_offset_db == -3.0
    @test occursin("mcs=2", sprint(show, action))
    @test_throws ArgumentError XAppAction(0, 0.5, 0.0)
    @test_throws ArgumentError XAppAction(1, -0.1, 0.0)
    @test_throws ArgumentError XAppAction(1, 1.1, 0.0)
    @test_throws ArgumentError XAppAction(1, NaN, 0.0)
    @test_throws ArgumentError XAppAction(1, 0.5, Inf)

    kpi = UEKPI("ue1", 20, 18, 17, 17 * 54, 20 * 72)
    @test kpi.bler ≈ 0.1
    @test kpi.throughput ≈ 17 * 54 / (20 * 72)
    @test occursin("BLER", sprint(show, kpi))
    empty_kpi = UEKPI("ue1", 0, 0, 0, 0, 0)
    @test empty_kpi.bler == 0.0
    @test empty_kpi.throughput == 0.0
    @test_throws ArgumentError UEKPI("", 1, 1, 1, 54, 72)
    @test_throws ArgumentError UEKPI("ue1", -1, 0, 0, 0, 0)
    @test_throws ArgumentError UEKPI("ue1", 10, 11, 11, 594, 720)
    @test_throws ArgumentError UEKPI("ue1", 10, 9, 10, 540, 720)
end

@testset "Jain fairness index" begin
    @test jain_index([1.0, 1.0, 1.0]) ≈ 1.0
    @test jain_index([1.0, 0.0, 0.0]) ≈ 1 / 3
    @test jain_index([3.0, 1.0]) ≈ 0.8
    @test jain_index([2.0, 1.0, 1.0]) ≈ 16 / (3 * 6)
    @test jain_index([5]) ≈ 1.0
    @test isnan(jain_index([0.0, 0.0]))
    @test_throws ArgumentError jain_index(Float64[])
    @test_throws ArgumentError jain_index([1.0, -1.0])
    @test_throws ArgumentError jain_index([1.0, Inf])
end

@testset "XApp policies" begin
    fixed = FixedXAppPolicy()
    @test fixed.mcs == 1
    @test fixed.shares === nothing
    @test fixed.power_offset_db == 0.0
    @test occursin("FixedXAppPolicy", sprint(show, fixed))
    @test_throws ArgumentError FixedXAppPolicy(mcs=0)
    @test_throws ArgumentError FixedXAppPolicy(shares=Float64[])
    @test_throws ArgumentError FixedXAppPolicy(shares=[0.5, -0.5])
    @test_throws ArgumentError FixedXAppPolicy(shares=[NaN, 1.0])
    @test_throws ArgumentError FixedXAppPolicy(shares=[0.0, 0.0])
    @test_throws ArgumentError FixedXAppPolicy(power_offset_db=NaN)

    rule = RuleXAppPolicy(3)
    @test rule.n_ues == 3
    @test rule.target_bler == 0.1
    @test rule.offsets_db == [-3.0, 0.0, 3.0]
    @test rule.thresholds_db == teaching_punctured_amc_thresholds()
    @test rule.offset_index == [2, 2, 2] # the rung closest to 0 dB
    @test occursin("RuleXAppPolicy", sprint(show, rule))
    @test_throws ArgumentError RuleXAppPolicy(0)
    @test_throws ArgumentError RuleXAppPolicy(2; target_bler=0.0)
    @test_throws ArgumentError RuleXAppPolicy(2; target_bler=1.0)
    @test_throws ArgumentError RuleXAppPolicy(2; offsets_db=Float64[])
    @test_throws ArgumentError RuleXAppPolicy(2; offsets_db=[0.0, Inf])
    @test_throws ArgumentError RuleXAppPolicy(2; offsets_db=[0.0, 0.0])
    @test_throws ArgumentError RuleXAppPolicy(2; offsets_db=[3.0, 0.0])
    @test_throws ArgumentError RuleXAppPolicy(2; thresholds_db=[NaN])
    @test_throws ArgumentError RuleXAppPolicy(2; thresholds_db=[1.0, 0.0])
end

@testset "xapp_actions rule table" begin
    ues = [XAppUE("ue1", 0.0), XAppUE("ue2", 3.0), XAppUE("ue3", 6.0)]
    zero_kpis = [UEKPI(ue.name, 0, 0, 0, 0, 0) for ue in ues]

    # Fixed baseline: uniform shares, constant row and offset; KPIs ignored.
    fixed = FixedXAppPolicy(mcs=2)
    actions = xapp_actions(fixed, ues, zero_kpis)
    @test [a.mcs for a in actions] == [2, 2, 2]
    @test [a.share for a in actions] ≈ fill(1 / 3, 3)
    @test [a.power_offset_db for a in actions] == [0.0, 0.0, 0.0]
    weighted = FixedXAppPolicy(mcs=1, shares=[2.0, 1.0, 1.0], power_offset_db=-3.0)
    wactions = xapp_actions(weighted, ues, zero_kpis)
    @test [a.share for a in wactions] ≈ [0.5, 0.25, 0.25]
    @test all(a -> a.power_offset_db == -3.0, wactions)
    @test_throws ArgumentError xapp_actions(fixed, ues, zero_kpis[1:2])
    @test_throws ArgumentError xapp_actions(fixed, ues, zero_kpis; power_budget=0.0)
    @test_throws ArgumentError xapp_actions(
        FixedXAppPolicy(mcs=1, shares=[1.0, 1.0]), ues, zero_kpis)
    # The fixed policy does no power control: +3 dB everywhere is mean power 2.
    @test_throws ArgumentError xapp_actions(
        FixedXAppPolicy(mcs=1, power_offset_db=3.0), ues, zero_kpis;
        power_budget=1.0)

    # Rule table on zero KPIs: equal weights, initial rung, genie thresholds.
    rule = RuleXAppPolicy(3)
    ractions = xapp_actions(rule, ues, zero_kpis)
    @test [a.share for a in ractions] ≈ fill(1 / 3, 3)
    @test [a.power_offset_db for a in ractions] == [0.0, 0.0, 0.0]
    @test [a.mcs for a in ractions] == [1, 2, 4] # 0 dB, 3 dB, 6 dB on [-Inf,2,4,6,8]

    # Less delivered so far, larger share next round.
    kpis = [
        UEKPI("ue1", 10, 10, 10, 100, 1440),
        UEKPI("ue2", 10, 10, 10, 200, 1440),
        UEKPI("ue3", 10, 10, 10, 400, 1440),
    ]
    wactions = xapp_actions(rule, ues, kpis)
    shares = [a.share for a in wactions]
    @test sum(shares) ≈ 1.0
    @test shares[1] > shares[2] > shares[3]
    @test shares ≈ [1 / 101, 1 / 201, 1 / 401] ./ sum([1 / 101, 1 / 201, 1 / 401])
    @test_throws ArgumentError xapp_actions(rule, ues[1:2], zero_kpis[1:2])
    @test_throws ArgumentError xapp_actions(rule, ues, zero_kpis[1:2])
    @test_throws ArgumentError xapp_actions(RuleXAppPolicy(2), ues, zero_kpis)

    # Budget projection: +3 dB on one UE and -3 dB on another at equal shares
    # costs 0.5*2 + 0.5*0.5 = 1.25 > 1, so the highest rung steps down.
    two = [XAppUE("ue1", 3.0), XAppUE("ue2", 3.0)]
    twokpi = [UEKPI(ue.name, 10, 10, 10, 540, 720) for ue in two]
    rule2 = RuleXAppPolicy(2)
    rule2.offset_index .= [3, 1]
    projected = xapp_actions(rule2, two, twokpi)
    @test [a.power_offset_db for a in projected] == [0.0, -3.0]
    cost = sum(a.share * snr_db_to_linear(a.power_offset_db) for a in projected)
    @test cost <= 1.0

    # A budget infeasible even at the minimum rungs is an ArgumentError.
    @test_throws ArgumentError xapp_actions(
        RuleXAppPolicy(2; offsets_db=[1.0, 2.0]), two, twokpi; power_budget=0.5)

    # The MCS row follows the projected offset, not the raw one: +3 dB would
    # select the higher row at 5 dB, the budget pulls it back to 0 dB and 2 dB.
    one = [XAppUE("ue1", 2.0)]
    onekpi = [UEKPI("ue1", 10, 10, 10, 540, 720)]
    rule3 = RuleXAppPolicy(1; offsets_db=[0.0, 3.0], thresholds_db=[-Inf, 4.0])
    rule3.offset_index .= [2]
    pulled = xapp_actions(rule3, one, onekpi; power_budget=1.0)
    @test only(pulled).power_offset_db == 0.0
    @test only(pulled).mcs == 1
    loose = xapp_actions(rule3, one, onekpi; power_budget=2.0)
    @test only(loose).power_offset_db == 3.0
    @test only(loose).mcs == 2
end

@testset "observe_kpi! offset walk" begin
    rule = RuleXAppPolicy(2)
    @test rule.offset_index == [2, 2]
    high = UEKPI("ue1", 10, 5, 5, 270, 720)         # BLER 0.5 > 0.1
    low = UEKPI("ue2", 100, 100, 100, 5400, 7200)   # BLER 0 < 0.01
    mid = UEKPI("ue1", 100, 95, 95, 5130, 7200)     # BLER 0.05, in between
    @test observe_kpi!(rule, 1, high) == 3
    @test observe_kpi!(rule, 1, high) == 3          # clamped at the top rung
    @test observe_kpi!(rule, 2, low) == 1
    @test observe_kpi!(rule, 2, low) == 1           # clamped at the bottom rung
    @test observe_kpi!(rule, 1, mid) == 3           # between the rails: hold
    @test observe_kpi!(rule, 2, UEKPI("ue2", 0, 0, 0, 0, 0)) == 1 # no blocks: hold
    @test_throws ArgumentError observe_kpi!(rule, 3, mid)
    @test observe_kpi!(FixedXAppPolicy(), 1, mid) === nothing
end

function _assert_xapp_report(report; n_ues::Int)
    @test length(report.ues) == n_ues
    @test length(report.kpis) == n_ues
    @test length(report.last_actions) == n_ues
    @test report.total_blocks == sum(k -> k.n_blocks, report.kpis)
    @test report.total_ack == sum(k -> k.n_ack, report.kpis)
    @test report.total_delivered == sum(k -> k.n_delivered, report.kpis)
    @test report.total_success_bits == sum(k -> k.n_success_bits, report.kpis)
    @test report.total_channel_uses == sum(k -> k.n_channel_uses, report.kpis)
    @test report.system_bler ≈ 1 - report.total_ack / report.total_blocks
    @test report.system_throughput ≈ report.total_success_bits / report.total_channel_uses
    @test sum(a.share for a in report.last_actions) ≈ 1.0
    return nothing
end

@testset "Noiseless xApp loop" begin
    table = teaching_punctured_mcs_table()
    ues = [XAppUE("ue1", Inf), XAppUE("ue2", Inf), XAppUE("ue3", Inf)]

    # Fixed baseline: 9 blocks a round split 3/3/3, every block delivered at
    # the row's spectral efficiency. Jain is exactly 1 on equal deliveries.
    Random.seed!(RNG_SEED)
    fixed = simulate_xapp(ues, FixedXAppPolicy(mcs=2), 2; blocks_per_round=9)
    _assert_xapp_report(fixed; n_ues=3)
    @test [k.n_blocks for k in fixed.kpis] == [6, 6, 6]
    @test fixed.total_delivered == 18
    @test fixed.system_bler == 0.0
    @test fixed.system_throughput ≈ spectral_efficiency(table[2])
    @test all(k -> k.throughput ≈ spectral_efficiency(table[2]), fixed.kpis)
    @test fixed.jain == 1.0
    @test fixed.total_channel_uses == 18 * table[2].n_symbols
    @test [a.mcs for a in fixed.last_actions] == [2, 2, 2]
    @test occursin("XAppReport", sprint(show, fixed))
    @test occursin("Jain", sprint(show, MIME("text/plain"), fixed))

    # A share of zero is never scheduled; its KPI stays empty and drags Jain
    # to 2/3 on two equal deliveries. 8 blocks a round split 4/4/0 with no
    # largest-remainder tie.
    Random.seed!(RNG_SEED)
    skewed = simulate_xapp(ues, FixedXAppPolicy(mcs=2, shares=[0.5, 0.5, 0.0]), 2;
                           blocks_per_round=8)
    _assert_xapp_report(skewed; n_ues=3)
    @test [k.n_blocks for k in skewed.kpis] == [8, 8, 0]
    @test skewed.kpis[3].bler == 0.0
    @test skewed.kpis[3].throughput == 0.0
    @test skewed.jain ≈ 2 / 3

    # Largest-remainder ties go to the lower UE index: one block a round with
    # equal shares lands on UE1 every round.
    Random.seed!(RNG_SEED)
    single = simulate_xapp(ues, FixedXAppPolicy(mcs=2), 4; blocks_per_round=1)
    _assert_xapp_report(single; n_ues=3)
    @test [k.n_blocks for k in single.kpis] == [4, 0, 0]

    # Rule table on a noiseless link: the genie thresholds pick the fastest
    # row, the BLER rail walks the offset down to the bottom rung, and the
    # delivered bits stay equal, so Jain stays 1.
    Random.seed!(RNG_SEED)
    ruled = simulate_xapp(ues, RuleXAppPolicy(3), 3; blocks_per_round=9)
    _assert_xapp_report(ruled; n_ues=3)
    @test ruled.system_bler == 0.0
    @test ruled.system_throughput ≈ spectral_efficiency(table[5])
    @test ruled.jain == 1.0
    @test [a.mcs for a in ruled.last_actions] == [5, 5, 5]
    @test [a.power_offset_db for a in ruled.last_actions] == [-3.0, -3.0, -3.0]

    # harq=:ir forwards to the link loop: noiseless first transmissions pass
    # the CRC on pattern A, so channel uses match one transmission per block.
    Random.seed!(RNG_SEED)
    ir = simulate_xapp(ues, FixedXAppPolicy(mcs=4), 2; blocks_per_round=6,
                       max_tx=2, harq=:ir)
    _assert_xapp_report(ir; n_ues=3)
    @test ir.system_bler == 0.0
    @test ir.total_channel_uses == 12 * table[4].n_symbols
end

@testset "xApp reproducibility and fading smoke" begin
    ues = [XAppUE("ue1", 0.0), XAppUE("ue2", 3.0), XAppUE("ue3", 6.0)]
    fields = fieldnames(XAppReport)
    for chan in (:awgn, :rayleigh)
        Random.seed!(RNG_SEED)
        first = simulate_xapp(ues, RuleXAppPolicy(3), 3; channel=chan)
        Random.seed!(RNG_SEED)
        second = simulate_xapp(ues, RuleXAppPolicy(3), 3; channel=chan)
        for f in fields
            @test getfield(first, f) == getfield(second, f)
        end
        _assert_xapp_report(first; n_ues=3)
        @test 0 <= first.system_bler <= 1
    end

    # Ideal-CSI round fading of a noiseless base link stays noiseless:
    # Inf + a finite fade is still Inf.
    Random.seed!(RNG_SEED)
    calm = simulate_xapp([XAppUE("ue1", Inf)], FixedXAppPolicy(mcs=2), 3;
                         channel=:rayleigh)
    @test calm.system_bler == 0.0
    @test calm.jain == 1.0
end

@testset "Rule policy versus fixed baseline" begin
    # Three UEs at 0, 3, and 6 dB on the punctured table. The baseline forces
    # QPSK-1/2 for everyone: the weak UE fails more than half its blocks, the
    # strong UE wastes capacity. The rule table gives the weak UE BPSK-1/2 and
    # the strong UE QPSK-3/4, so it delivers a clearly higher fraction of the
    # 600 offered blocks, with higher system throughput and a Jain index much
    # closer to 1. Seed 0x0c0de1162, Clopper–Pearson α = 1e-4 on the delivered
    # counts; throughput and fairness are point-estimate comparisons on the
    # same draws. Not an NR scheduler comparison, not O-RAN.
    ues = [XAppUE("ue1", 0.0), XAppUE("ue2", 3.0), XAppUE("ue3", 6.0)]
    n_rounds = 20
    Random.seed!(RNG_SEED)
    base = simulate_xapp(ues, FixedXAppPolicy(mcs=2), n_rounds)
    Random.seed!(RNG_SEED)
    rule = simulate_xapp(ues, RuleXAppPolicy(3), n_rounds)
    _assert_xapp_report(base; n_ues=3)
    _assert_xapp_report(rule; n_ues=3)
    @test base.total_blocks == rule.total_blocks == n_rounds * 30
    nack_b = base.total_blocks - base.total_ack
    nack_r = rule.total_blocks - rule.total_ack
    @test nack_b >= 30
    @test nack_r >= 30
    _, hi_b = clopper_pearson(base.total_delivered, base.total_blocks, BER_ALPHA)
    lo_r, _ = clopper_pearson(rule.total_delivered, rule.total_blocks, BER_ALPHA)
    @test lo_r > hi_b
    @test rule.system_throughput > base.system_throughput
    @test rule.jain > base.jain
    @test base.last_actions == [XAppAction(2, 1 / 3, 0.0) for _ in 1:3]
end

@testset "simulate_xapp argument checks" begin
    ues = [XAppUE("ue1", 3.0), XAppUE("ue2", 6.0)]
    rule = RuleXAppPolicy(2)
    @test_throws ArgumentError simulate_xapp(XAppUE[], FixedXAppPolicy(), 1)
    @test_throws ArgumentError simulate_xapp(ues, rule, 0)
    @test_throws ArgumentError simulate_xapp(ues, rule, 1; blocks_per_round=0)
    @test_throws ArgumentError simulate_xapp(ues, rule, 1; table=MCS[])
    @test_throws ArgumentError simulate_xapp(ues, rule, 1; max_tx=0)
    @test_throws ArgumentError simulate_xapp(ues, rule, 1; harq=:typeii)
    @test_throws ArgumentError simulate_xapp(ues, rule, 1; approx=:maxlog2)
    @test_throws ArgumentError simulate_xapp(ues, rule, 1; channel=:geometry)
    @test_throws ArgumentError simulate_xapp(ues, rule, 1; power_budget=0.0)
    @test_throws ArgumentError simulate_xapp(ues, rule, 1; power_budget=Inf)
    @test_throws ArgumentError simulate_xapp(ues, RuleXAppPolicy(3), 1)
    @test_throws ArgumentError simulate_xapp(
        ues, RuleXAppPolicy(2; thresholds_db=[-Inf, 0.0]), 1)
    @test_throws ArgumentError simulate_xapp(ues, FixedXAppPolicy(mcs=99), 1)
    @test_throws ArgumentError simulate_xapp(
        ues, FixedXAppPolicy(shares=[1.0, 1.0, 1.0]), 1)
    # A fixed offset that breaks the budget is rejected on the first round.
    @test_throws ArgumentError simulate_xapp(
        ues, FixedXAppPolicy(power_offset_db=3.0), 1; power_budget=1.0)
end
