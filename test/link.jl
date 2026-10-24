# Teaching CRC, MCS table, Chase combining, and threshold AMC.
# Noise is the package convention σ² = N0/2 with Es = 1, so the link SNR is Es/N0.
# Intervals are Clopper–Pearson at BER_ALPHA (reference.jl). Counts are not locked:
# a different randn stream must still separate the intervals below.

@testset "CRC-16-CCITT" begin
    @test CRC16_CCITT_POLY == 0x1021
    @test crc16_ccitt(Int[]) == 0x0000
    empty_word = crc_append(Int[])
    @test empty_word == zeros(Int, 16)
    @test crc_check(empty_word)

    # CRC-16/XMODEM check value: ASCII "123456789", each byte MSB first.
    ascii_bits = Int[]
    for byte in codeunits("123456789")
        for shift in 7:-1:0
            push!(ascii_bits, (Int(byte) >> shift) & 1)
        end
    end
    @test crc16_ccitt(ascii_bits) == 0x31c3
    @test crc_check(crc_append(ascii_bits))
    @test crc16_ccitt(crc_append(ascii_bits)) == 0x0000

    # g(x) = x^16 + x^12 + x^5 + 1. A shift of g is an undetectable error.
    generator = zeros(Int, 17)
    generator[1] = 1
    generator[5] = 1
    generator[12] = 1
    generator[17] = 1
    word = crc_append(zeros(Int, 24))
    forged = copy(word)
    for i in eachindex(generator)
        forged[8 + i] ⊻= generator[i]
    end
    @test crc_check(forged)
    @test crc16_ccitt(forged) == 0x0000

    # Every odd-weight pattern is detectable because g(1) = 0. One-bit flips
    # cover the link payloads; three bits is another odd weight.
    Random.seed!(RNG_SEED)
    protected = crc_append(rand(0:1, 48))
    for i in eachindex(protected)
        bad = copy(protected)
        bad[i] ⊻= 1
        @test !crc_check(bad)
    end
    odd = copy(protected)
    odd[1] ⊻= 1
    odd[2] ⊻= 1
    odd[4] ⊻= 1
    @test !crc_check(odd)

    @test_throws ArgumentError crc16_ccitt([0, 2])
    @test_throws ArgumentError crc_append([0, 1, 3])
    @test_throws ArgumentError crc_check([0, 1, 0])
end

@testset "Teaching MCS table" begin
    table = teaching_mcs_table()
    @test [row.name for row in table] == [
        "BPSK-conv-1/2", "QPSK-conv-1/2", "QPSK-array-LDPC",
        "16QAM-conv-1/2", "64QAM-conv-1/2",
    ]
    ses = spectral_efficiency.(table)
    @test issorted(ses, lt=<)
    @test ses ≈ [48 / 132, 48 / 66, 41 / 44, 48 / 33, 48 / 22]
    for row in table
        @test row.crc_bits == 16
        @test row.n_encoder_in == row.payload_bits + 16
        @test row.n_coded % bits_per_symbol(row.scheme) == 0
        @test row.n_symbols == row.n_coded ÷ bits_per_symbol(row.scheme)
        @test spectral_efficiency(row) ≈ row.payload_bits / row.n_symbols
        @test occursin(row.name, sprint(show, row))
    end
    ldpc_row = table[3]
    @test ldpc_row.family === :ldpc
    @test ldpc_row.ldpc.n == 88
    @test ldpc_row.ldpc.k == 57
    @test ldpc_row.ldpc.m == 33
    @test ldpc_row.ldpc.rank == 31
    @test ldpc_row.payload_bits == 41
    @test table[1].family === :convolutional
    @test table[1].ldpc === nothing
    # A second call returns another vector and the same code.
    @test teaching_mcs_table()[3].ldpc === ldpc_row.ldpc

    @test_throws ArgumentError mcs_convolutional("empty", BPSK(), 0)
    @test_throws ArgumentError mcs_convolutional("not-a-multiple", QAM(16), 1)
    @test_throws ArgumentError mcs_convolutional("not-a-multiple", QAM(64), 1)
    @test_throws ArgumentError mcs_ldpc("short", BPSK(), hamming74_ldpc())

    edges = teaching_amc_thresholds()
    @test length(edges) == length(table)
    @test edges[1] == -Inf
    @test edges[2:end] == [2.0, 3.0, 7.0, 11.0]
    policy = ThresholdAMC(edges)
    @test select_mcs(policy, -5.0) == 1
    @test select_mcs(policy, 0.0) == 1
    @test select_mcs(policy, 2.0) == 2
    @test select_mcs(policy, 2.9) == 2
    @test select_mcs(policy, 3.0) == 3
    @test select_mcs(policy, 6.9) == 3
    @test select_mcs(policy, 7.0) == 4
    @test select_mcs(policy, 10.9) == 4
    @test select_mcs(policy, 11.0) == 5
    @test select_mcs(policy, Inf) == 5
    # Equal thresholds: the higher index wins, so a dominated row is skipped.
    tied = ThresholdAMC([-Inf, 5.0, 5.0])
    @test select_mcs(tied, 5.0) == 3
    @test select_mcs(tied, 4.9) == 1
    @test select_mcs(ThresholdAMC([0.0, 5.0]), -10.0) == 1
    @test occursin("ThresholdAMC", sprint(show, policy))

    @test_throws ArgumentError ThresholdAMC(Float64[])
    @test_throws ArgumentError ThresholdAMC([1.0, 0.0])
    @test_throws ArgumentError ThresholdAMC([NaN])
    @test_throws ArgumentError select_mcs(policy, NaN)
    @test_throws ArgumentError select_mcs(policy, -Inf)
end

@testset "Outer-loop offset" begin
    policy = OuterLoopAMC([-Inf, 0.0]; target_bler=0.1, offset_db=0.0, step_db=0.5,
                           min_offset_db=-10.0, max_offset_db=10.0)
    @test policy.step_up_db ≈ 0.5
    @test policy.step_down_db ≈ 0.5 * 0.9 / 0.1
    @test observe_acknack!(policy, true) ≈ -0.5
    @test observe_acknack!(policy, false) ≈ -0.5 + 4.5
    @test select_mcs(policy, 4.0) == select_mcs(
        ThresholdAMC(policy.thresholds_db), 4.0 - policy.offset_db,
    )
    clamped = OuterLoopAMC([-Inf]; target_bler=0.5, offset_db=0.0, step_db=1.0,
                            min_offset_db=-2.0, max_offset_db=2.0)
    observe_acknack!(clamped, false)
    observe_acknack!(clamped, false)
    observe_acknack!(clamped, false)
    @test clamped.offset_db == 2.0
    observe_acknack!(clamped, true)
    observe_acknack!(clamped, true)
    observe_acknack!(clamped, true)
    observe_acknack!(clamped, true)
    @test clamped.offset_db == -2.0
    @test occursin("OuterLoopAMC", sprint(show, clamped))
    @test_throws ArgumentError OuterLoopAMC([-Inf]; target_bler=0.0)
    @test_throws ArgumentError OuterLoopAMC([-Inf]; step_db=0.0)
    @test_throws ArgumentError OuterLoopAMC([-Inf]; min_offset_db=2.0, max_offset_db=1.0)
end

function _assert_link_report(report; max_tx::Int)
    @test 0 <= report.n_delivered <= report.n_ack <= report.n_blocks
    @test report.n_undetected == report.n_ack - report.n_delivered
    @test 0 <= report.n_bit_errors <= report.n_payload_bits
    @test report.n_payload_bits > 0
    @test report.ber ≈ report.n_bit_errors / report.n_payload_bits
    @test report.bler ≈ (report.n_blocks - report.n_ack) / report.n_blocks
    @test report.n_transmissions >= report.n_blocks
    @test report.n_transmissions <= report.n_blocks * max_tx
    @test report.n_channel_uses >= report.n_transmissions
    @test report.throughput ≈ report.n_success_bits / report.n_channel_uses
    @test report.mean_transmissions ≈ report.n_transmissions / report.n_blocks
    @test sum(report.mcs_histogram) == report.n_blocks
    @test length(report.mcs_histogram) == length(report.mcs_names)
    @test report.spectral_efficiency >= 0
    @test report.mean_mcs_efficiency >= 0
    return nothing
end

@testset "Noiseless and high-SNR link" begin
    table = teaching_mcs_table()
    for (i, row) in enumerate(table)
        Random.seed!(RNG_SEED)
        report = simulate_link(Inf, 2; table=table, mcs=i, max_tx=3)
        _assert_link_report(report; max_tx=3)
        @test report.bler == 0.0
        @test report.ber == 0.0
        @test report.n_undetected == 0
        @test report.n_delivered == 2
        @test report.mean_transmissions == 1.0
        @test report.mcs_histogram[i] == 2
        @test report.throughput ≈ spectral_efficiency(row)
        @test report.spectral_efficiency ≈ report.throughput
        @test report.mean_mcs_efficiency ≈ report.throughput
        @test report.n_channel_uses == 2 * row.n_symbols
    end
    # 16 dB is above the 64-QAM row. The lower Clopper–Pearson end of the
    # delivered-block rate, times that row's SE, still clears the 16-QAM SE.
    Random.seed!(RNG_SEED)
    high = simulate_link(16.0, 40; mcs=5, max_tx=1)
    _assert_link_report(high; max_tx=1)
    lo, _ = clopper_pearson(high.n_delivered, high.n_blocks, BER_ALPHA)
    @test lo * spectral_efficiency(table[5]) > spectral_efficiency(table[4])
    @test high.mcs_histogram == [0, 0, 0, 0, 40]
    @test occursin("BLER", sprint(show, MIME("text/plain"), high))
    @test occursin("LinkReport", sprint(show, high))

    Random.seed!(RNG_SEED)
    soft = simulate_link(Inf, 1; mcs=3, ldpc_algorithm=:sum_product, approx=:maxlog)
    @test soft.ber == 0.0
    Random.seed!(RNG_SEED)
    qam = simulate_link(Inf, 1; mcs=4, approx=:maxlog, ldpc_max_iter=0)
    @test qam.ber == 0.0
end

@testset "Chase combining versus a single transmission" begin
    # 16-QAM convolutional row at 5 dB Es/N0. One transmission is mostly CRC
    # failures; three Chase transmissions are mostly ACKs. n = 100 keeps the
    # test inside CI. α = 10^-4.
    n = 100
    Random.seed!(RNG_SEED)
    single = simulate_link(5.0, n; mcs=4, max_tx=1)
    Random.seed!(RNG_SEED)
    chase = simulate_link(5.0, n; mcs=4, max_tx=3)
    _assert_link_report(single; max_tx=1)
    _assert_link_report(chase; max_tx=3)
    @test single.mcs_histogram[4] == n
    @test chase.mcs_histogram[4] == n
    @test single.n_channel_uses == n * teaching_mcs_table()[4].n_symbols
    @test chase.n_channel_uses == chase.n_transmissions * teaching_mcs_table()[4].n_symbols
    @test chase.mean_transmissions > 1
    @test chase.n_channel_uses > single.n_channel_uses
    nack_s = single.n_blocks - single.n_ack
    nack_c = chase.n_blocks - chase.n_ack
    @test nack_s >= 30
    lo_s, _ = clopper_pearson(nack_s, n, BER_ALPHA)
    _, hi_c = clopper_pearson(nack_c, n, BER_ALPHA)
    @test hi_c < lo_s
    @test chase.throughput > single.throughput
    @test chase.ber < single.ber
end

@testset "Adaptive MCS versus fixed rows" begin
    table = teaching_mcs_table()
    policy = ThresholdAMC(teaching_amc_thresholds())
    n = 80
    # 5 dB selects the QPSK array LDPC on the published thresholds.
    # BPSK is the conservative row. 64-QAM is the aggressive row.
    Random.seed!(RNG_SEED)
    adapt = simulate_link(5.0, n; table=table, policy=policy)
    Random.seed!(RNG_SEED)
    conservative = simulate_link(5.0, n; table=table, mcs=1)
    Random.seed!(RNG_SEED)
    aggressive = simulate_link(5.0, n; table=table, mcs=5)
    _assert_link_report(adapt; max_tx=1)
    _assert_link_report(conservative; max_tx=1)
    _assert_link_report(aggressive; max_tx=1)
    @test adapt.mcs_histogram == [0, 0, n, 0, 0]
    @test conservative.mcs_histogram[1] == n
    @test aggressive.mcs_histogram[5] == n
    @test adapt.n_undetected == 0
    @test conservative.n_undetected == 0
    @test aggressive.n_undetected == 0
    se_a = spectral_efficiency(table[3])
    se_c = spectral_efficiency(table[1])
    se_g = spectral_efficiency(table[5])
    lo_a, _ = clopper_pearson(adapt.n_delivered, n, BER_ALPHA)
    _, hi_c = clopper_pearson(conservative.n_delivered, n, BER_ALPHA)
    _, hi_g = clopper_pearson(aggressive.n_delivered, n, BER_ALPHA)
    # max_tx = 1 and one MCS: throughput = (delivered blocks) / n * SE.
    @test adapt.throughput ≈ adapt.n_delivered / n * se_a
    @test lo_a * se_a > hi_c * se_c
    @test lo_a * se_a > hi_g * se_g

    # 16 dB selects 64-QAM. Its throughput lower bound clears the BPSK row.
    Random.seed!(RNG_SEED)
    high = simulate_link(16.0, n; table=table, policy=policy)
    Random.seed!(RNG_SEED)
    stuck = simulate_link(16.0, n; table=table, mcs=1)
    @test high.mcs_histogram[5] == n
    @test stuck.mcs_histogram[1] == n
    lo_h, _ = clopper_pearson(high.n_delivered, n, BER_ALPHA)
    _, hi_s = clopper_pearson(stuck.n_delivered, n, BER_ALPHA)
    @test lo_h * spectral_efficiency(table[5]) > hi_s * spectral_efficiency(table[1])

    # 0 dB stays on BPSK. Selection does not depend on the noise draw.
    Random.seed!(RNG_SEED)
    low = simulate_link(0.0, 12; table=table, policy=policy)
    @test low.mcs_histogram[1] == 12
end

@testset "Outer loop and moment SNR" begin
    rows = [
        mcs_convolutional("BPSK-conv-1/2", BPSK(), 48),
        mcs_convolutional("QPSK-conv-1/2", QPSK(), 48),
    ]
    # Offset 20 dB puts the first block below the QPSK threshold (24 − 20 = 4).
    # Each ACK removes 8 dB, so the next three blocks sit on QPSK.
    policy = OuterLoopAMC([-Inf, 5.0]; target_bler=0.1, offset_db=20.0, step_db=8.0,
                           min_offset_db=-30.0, max_offset_db=40.0)
    Random.seed!(RNG_SEED)
    climbed = simulate_link(24.0, 4; table=rows, policy=policy, max_tx=1)
    _assert_link_report(climbed; max_tx=1)
    @test climbed.n_ack == 4
    @test climbed.n_delivered == 4
    @test climbed.mcs_histogram == [1, 3]
    @test policy.offset_db ≈ 20.0 - 4 * 8.0
    uses = rows[1].n_symbols + 3 * rows[2].n_symbols
    @test climbed.n_channel_uses == uses
    @test climbed.throughput ≈ (4 * 48) / uses
    @test climbed.spectral_efficiency ≈
          (spectral_efficiency(rows[1]) + 3 * spectral_efficiency(rows[2])) / 4
    @test climbed.mean_mcs_efficiency ≈ climbed.spectral_efficiency
    @test climbed.throughput < climbed.spectral_efficiency

    frozen = OuterLoopAMC([-Inf, 5.0]; offset_db=1.25, step_db=0.5)
    before = frozen.offset_db
    Random.seed!(RNG_SEED)
    fixed = simulate_link(Inf, 2; table=rows, mcs=1, policy=frozen, max_tx=2)
    @test fixed.ber == 0.0
    @test fixed.mean_transmissions == 1.0
    @test frozen.offset_db == before

    # Both rows are BPSK, so every noiseless symbol has energy exactly 1 and the
    # moment estimate is +∞. A higher-order constellation in one short block does
    # not, because the empirical energy is not the ensemble energy.
    bpsk_rows = [
        mcs_convolutional("BPSK-low", BPSK(), 48),
        mcs_convolutional("BPSK-high", BPSK(), 96),
    ]
    Random.seed!(RNG_SEED)
    tracked = simulate_link(Inf, 4; table=bpsk_rows, policy=ThresholdAMC([-Inf, 0.0]),
                            channel_quality=:moment, initial_esn0_db=-20.0, snr_window=2)
    _assert_link_report(tracked; max_tx=1)
    @test tracked.mcs_histogram == [1, 3]
    @test tracked.bler == 0.0
    tracked_uses = bpsk_rows[1].n_symbols + 3 * bpsk_rows[2].n_symbols
    tracked_bits = bpsk_rows[1].payload_bits + 3 * bpsk_rows[2].payload_bits
    @test tracked.n_channel_uses == tracked_uses
    @test tracked.throughput ≈ tracked_bits / tracked_uses

    Random.seed!(RNG_SEED)
    finite = simulate_link(8.0, 4; channel_quality=:moment, snr_window=3)
    _assert_link_report(finite; max_tx=1)
    @test isfinite(finite.throughput)
end

@testset "Es/N0 moment estimator" begin
    # Real BPSK: E[y^2] = 1 + N0/2. Power 1.25 ⇒ N0 = 0.5 ⇒ Es/N0 = 3.0103 dB.
    @test estimate_esn0_db(fill(sqrt(1.25), 8), BPSK()) ≈ 10 * log10(2) rtol = 1e-12
    # Complex: E[|y|^2] = 1 + N0. The same N0.
    @test estimate_esn0_db(fill(ComplexF64(sqrt(1.5)), 8), QPSK()) ≈ 10 * log10(2) rtol = 1e-12
    @test estimate_esn0_db(ones(4), BPSK()) == Inf
    @test estimate_esn0_db(constellation(QPSK()), QPSK()) == Inf
    @test estimate_esn0_db([0.0], BPSK()) == Inf
    @test_throws ArgumentError estimate_esn0_db(Float64[], BPSK())
end

@testset "Link argument checks" begin
    @test_throws ArgumentError simulate_link(NaN, 1)
    @test_throws ArgumentError simulate_link(0.0, 0)
    @test_throws ArgumentError simulate_link(0.0, 1; max_tx=0)
    @test_throws ArgumentError simulate_link(0.0, 1; approx=:approx)
    @test_throws ArgumentError simulate_link(0.0, 1; channel_quality=:pilot)
    @test_throws ArgumentError simulate_link(0.0, 1; mcs=0)
    @test_throws ArgumentError simulate_link(0.0, 1; mcs=9)
    @test_throws ArgumentError simulate_link(0.0, 1; table=teaching_mcs_table()[1:2])
    @test_throws ArgumentError simulate_link(0.0, 1; table=MCS[])
    @test_throws ArgumentError simulate_link(0.0, 1; ldpc_algorithm=:turbo)
    @test_throws ArgumentError simulate_link(0.0, 1; ldpc_max_iter=-1)
    @test_throws ArgumentError simulate_link(0.0, 1; initial_esn0_db=0.0)
    @test_throws ArgumentError simulate_link(0.0, 1; channel_quality=:moment, snr_window=0)
    @test_throws ArgumentError simulate_link(0.0, 1; channel_quality=:moment,
                                             initial_esn0_db=NaN)
    @test_throws ArgumentError simulate_link(0.0, 1; receiver=:nr)
    @test_throws ArgumentError simulate_link(0.0, 1; bicm_iterations=4)
    @test_throws ArgumentError simulate_link(0.0, 1; bicm_n_rows=2)
    @test_throws ArgumentError simulate_link(0.0, 1; bicm_labels=bicm_psk8_id_labels())
    @test_throws ArgumentError simulate_link(0.0, 1; bicm_llr_clip=10)
    @test_throws ArgumentError simulate_link(Inf, 1; receiver=:bicm)
    @test_throws ArgumentError simulate_link(Inf, 1; receiver=:bicm, mcs=3)
    @test_throws ArgumentError simulate_link(Inf, 1; table=teaching_punctured_mcs_table(),
                                             receiver=:bicm, bicm_n_rows=5)
    @test_throws ArgumentError simulate_link(Inf, 1; receiver=:bicm, mcs=1, bicm_iterations=0)
end

@testset "Punctured convolutional MCS" begin
    table = teaching_punctured_mcs_table()
    @test [row.name for row in table] == [
        "BPSK-1/2", "QPSK-1/2", "QPSK-2/3", "QPSK-3/4", "16QAM-1/2",
    ]
    @test [row.rate for row in table] == [1//2, 1//2, 2//3, 3//4, 1//2]
    @test all(row -> row.family === :convolutional && row.ldpc === nothing, table)
    @test all(row -> row.payload_bits == 54 && row.n_mother == 144, table)
    @test [row.n_coded for row in table] == [144, 144, 108, 96, 144]
    @test [row.n_symbols for row in table] == [144, 72, 54, 48, 36]
    ses = spectral_efficiency.(table)
    @test issorted(ses, lt=<)
    @test ses ≈ [54 / 144, 54 / 72, 1.0, 54 / 48, 54 / 36]
    @test teaching_punctured_mcs_table()[3].rate == 2//3
    @test occursin("rate=2//3", sprint(show, table[3]))
    @test occursin("QPSK-2/3", sprint(show, table[3]))

    # The original five-row table is still the unpunctured list.
    old = teaching_mcs_table()
    @test old[1].rate == 1//2
    @test old[1].n_mother == old[1].n_coded == 132
    @test old[1].n_symbols == 132
    @test old[3].rate === nothing
    @test old[3].n_mother == 88
    @test spectral_efficiency.(old) ≈ [48 / 132, 48 / 66, 41 / 44, 48 / 33, 48 / 22]

    @test_throws ArgumentError mcs_convolutional("float", QPSK(), 54; rate=0.5)
    @test_throws ArgumentError mcs_convolutional("other", QPSK(), 54; rate=4//5)
    @test_throws ArgumentError mcs_convolutional("partial", QPSK(), 1; rate=2//3)

    for (i, row) in enumerate(table)
        Random.seed!(RNG_SEED)
        report = simulate_link(Inf, 2; table=table, mcs=i, max_tx=3)
        _assert_link_report(report; max_tx=3)
        @test report.bler == 0.0
        @test report.ber == 0.0
        @test report.n_delivered == 2
        @test report.mean_transmissions == 1.0
        @test report.throughput ≈ spectral_efficiency(row)
        @test report.n_channel_uses == 2 * row.n_symbols
        @test report.mcs_histogram[i] == 2
    end

    # Same CRC word, Gray labels, one demapper pass. Log-MAP on a noiseless
    # punctured block still passes the CRC. n_rows divides the transmitted length.
    for (idx, n_rows, iters) in ((3, 6, 2), (4, 8, 2), (2, 12, 1))
        Random.seed!(RNG_SEED)
        report = simulate_link(Inf, 1; table=table, mcs=idx, max_tx=2,
                               receiver=:bicm, bicm_iterations=iters, bicm_n_rows=n_rows)
        _assert_link_report(report; max_tx=2)
        @test report.ber == 0.0
        @test report.n_delivered == 1
        @test report.mean_transmissions == 1.0
        @test report.throughput ≈ spectral_efficiency(table[idx])
    end

    # 8-PSK set partitioning is not a Gray Viterbi row. The BICM receiver owns the labels.
    psk = mcs_convolutional("8PSK-1/2", PSK(8), 48; rate=1//2)
    @test psk.n_mother == 132
    @test psk.n_coded == 132
    @test psk.n_symbols == 44
    Random.seed!(RNG_SEED)
    psk_report = simulate_link(Inf, 1; table=[psk], mcs=1, receiver=:bicm,
                               bicm_iterations=3, bicm_n_rows=12,
                               bicm_labels=bicm_psk8_id_labels(), approx=:maxlog)
    @test psk_report.ber == 0.0
    @test psk_report.throughput ≈ spectral_efficiency(psk)

    # Julia normalises 2//4 to 1//2, so this is the mother code, not a new mask.
    normalised = mcs_convolutional("norm", BPSK(), 48; rate=2//4)
    @test normalised.rate == 1//2
    @test normalised.n_coded == normalised.n_mother
end

@testset "Punctured adaptive MCS versus fixed rows" begin
    # Design sweep: Julia 1.10.12, seed RNG_SEED, 80 blocks, no HARQ.
    # 3 dB is inside the QPSK-1/2 bin. 10 dB is inside the 16-QAM bin.
    # 5 dB is a near-tie and is not used. α = 1e-4, same interval as the
    # unpunctured adaptive test. Not an NR CQI comparison.
    table = teaching_punctured_mcs_table()
    edges = teaching_punctured_amc_thresholds()
    @test length(edges) == length(table)
    @test edges[1] == -Inf
    @test edges[2:end] == [2.0, 4.0, 6.0, 8.0]
    policy = ThresholdAMC(edges)
    @test select_mcs(policy, 0.0) == 1
    @test select_mcs(policy, 1.9) == 1
    @test select_mcs(policy, 2.0) == 2
    @test select_mcs(policy, 3.0) == 2
    @test select_mcs(policy, 3.9) == 2
    @test select_mcs(policy, 4.0) == 3
    @test select_mcs(policy, 5.9) == 3
    @test select_mcs(policy, 6.0) == 4
    @test select_mcs(policy, 7.9) == 4
    @test select_mcs(policy, 8.0) == 5
    @test select_mcs(policy, 10.0) == 5

    n = 80
    Random.seed!(RNG_SEED)
    adapt = simulate_link(3.0, n; table=table, policy=policy)
    Random.seed!(RNG_SEED)
    conservative = simulate_link(3.0, n; table=table, mcs=1)
    Random.seed!(RNG_SEED)
    aggressive = simulate_link(3.0, n; table=table, mcs=5)
    _assert_link_report(adapt; max_tx=1)
    _assert_link_report(conservative; max_tx=1)
    _assert_link_report(aggressive; max_tx=1)
    @test adapt.mcs_histogram == [0, n, 0, 0, 0]
    @test conservative.mcs_histogram[1] == n
    @test aggressive.mcs_histogram[5] == n
    se_a = spectral_efficiency(table[2])
    se_c = spectral_efficiency(table[1])
    se_g = spectral_efficiency(table[5])
    @test adapt.throughput ≈ adapt.n_delivered / n * se_a
    lo_a, _ = clopper_pearson(adapt.n_delivered, n, BER_ALPHA)
    _, hi_c = clopper_pearson(conservative.n_delivered, n, BER_ALPHA)
    _, hi_g = clopper_pearson(aggressive.n_delivered, n, BER_ALPHA)
    @test lo_a * se_a > hi_c * se_c
    @test lo_a * se_a > hi_g * se_g

    Random.seed!(RNG_SEED)
    high = simulate_link(10.0, n; table=table, policy=policy)
    Random.seed!(RNG_SEED)
    stuck = simulate_link(10.0, n; table=table, mcs=1)
    @test high.mcs_histogram[5] == n
    @test stuck.mcs_histogram[1] == n
    lo_h, _ = clopper_pearson(high.n_delivered, n, BER_ALPHA)
    _, hi_s = clopper_pearson(stuck.n_delivered, n, BER_ALPHA)
    @test lo_h * spectral_efficiency(table[5]) > hi_s * spectral_efficiency(table[1])

    Random.seed!(RNG_SEED)
    low = simulate_link(0.0, 12; table=table, policy=policy)
    @test low.mcs_histogram[1] == 12
end

@testset "Chase combining of one punctured word" begin
    # QPSK rate 3/4 at 2 dB Es/N0. One transmission is mostly CRC failures.
    # Three Chase transmissions of that same punctured word are mostly ACKs.
    # The repeat does not send the bits the mask deleted. n = 80, α = 1e-4.
    table = teaching_punctured_mcs_table()
    n = 80
    Random.seed!(RNG_SEED)
    single = simulate_link(2.0, n; table=table, mcs=4, max_tx=1)
    Random.seed!(RNG_SEED)
    chase = simulate_link(2.0, n; table=table, mcs=4, max_tx=3)
    _assert_link_report(single; max_tx=1)
    _assert_link_report(chase; max_tx=3)
    @test single.mcs_histogram[4] == n
    @test chase.mcs_histogram[4] == n
    @test single.n_channel_uses == n * table[4].n_symbols
    @test chase.n_channel_uses == chase.n_transmissions * table[4].n_symbols
    @test chase.mean_transmissions > 1
    nack_s = single.n_blocks - single.n_ack
    nack_c = chase.n_blocks - chase.n_ack
    @test nack_s >= 40
    lo_s, _ = clopper_pearson(nack_s, n, BER_ALPHA)
    _, hi_c = clopper_pearson(nack_c, n, BER_ALPHA)
    @test hi_c < lo_s
    @test chase.throughput > single.throughput
    @test chase.ber < single.ber
end

@testset "BICM receiver versus one demapper pass" begin
    # Terminated (7, 5), payload 48, CRC-16, 8-PSK set partitioning, 12-row
    # interleaver, max-log. Es/N0 = 6 dB, 80 blocks, seed RNG_SEED.
    # iterations = 4 against iterations = 1 on the same AWGN draws.
    # Not Gray, not NR, not a capacity claim. Clopper–Pearson α = 1e-4 on
    # payload bit errors and on CRC failures, the same interval as test/bicm.jl.
    # max_tx = 1, so this is not a HARQ gain.
    row = mcs_convolutional("8PSK-1/2", PSK(8), 48)
    table = [row]
    labels = bicm_psk8_id_labels()
    n = 80
    function run_bicm(iterations)
        Random.seed!(RNG_SEED)
        return simulate_link(6.0, n; table=table, mcs=1, receiver=:bicm,
                             bicm_iterations=iterations, bicm_n_rows=12,
                             bicm_labels=labels, approx=:maxlog)
    end
    one = run_bicm(1)
    four = run_bicm(4)
    _assert_link_report(one; max_tx=1)
    _assert_link_report(four; max_tx=1)
    @test one.n_payload_bits == four.n_payload_bits == n * row.payload_bits
    @test one.n_channel_uses == four.n_channel_uses == n * row.n_symbols
    @test one.n_bit_errors >= 200
    @test four.n_bit_errors >= 40
    @test four.n_bit_errors < one.n_bit_errors
    lo1, _ = clopper_pearson(one.n_bit_errors, one.n_payload_bits, BER_ALPHA)
    _, hi4 = clopper_pearson(four.n_bit_errors, four.n_payload_bits, BER_ALPHA)
    @test hi4 < lo1
    nack1 = one.n_blocks - one.n_ack
    nack4 = four.n_blocks - four.n_ack
    @test nack1 >= 60
    lo_b, _ = clopper_pearson(nack1, n, BER_ALPHA)
    _, hi_b = clopper_pearson(nack4, n, BER_ALPHA)
    @test hi_b < lo_b
    @test four.throughput > one.throughput
    @test four.ber < one.ber
end
