# Gardner timing, Costas carrier recovery, Schmidl–Cox preamble sync.
# Direction and known-offset checks. Not a BER curve, and not a claim that
# the loops are production synchronisers.

function _circ(a, b, period)
    d = mod(a - b, period)
    return min(d, period - d)
end

function _mixed_bits(n)
    return [((i * 17 + 3) % 7) > 3 ? 1 : 0 for i in 1:n]
end

@testset "Fractional delay" begin
    ramp = [0.0, 1.0, 2.0, 3.0, 4.0]
    @test fractional_delay(ramp, 0) == ramp
    @test eltype(fractional_delay(ramp, 0)) == Float64
    @test fractional_delay(1:5, 2) == [3, 4, 5]
    @test eltype(fractional_delay(1:5, 2)) == Int
    @test fractional_delay(ramp, 0.5) ≈ [0.5, 1.5, 2.5, 3.5]
    z = ComplexF64[1, 2, 3, 4]
    @test fractional_delay(z, 1) == ComplexF64[2, 3, 4]
    @test fractional_delay(z, 0.25)[1] ≈ 1.25

    @test_throws ArgumentError fractional_delay(ramp, -0.1)
    @test_throws ArgumentError fractional_delay(ramp, Inf)
    @test_throws ArgumentError fractional_delay(Float64[], 0)
    @test_throws ArgumentError fractional_delay(ramp, 5)
end

@testset "Gardner detector" begin
    # Hand triple, sps = 2, phase = 1. Late index 3: early=1, mid=0.5, late=-1.
    y = [1.0, 0.5, -1.0, -0.5, 1.0]
    e = gardner_ted(y; sps=2, phase=1)
    @test e[1] ≈ 0.5 * ((-1.0) - 1.0)
    @test e[2] ≈ (-0.5) * (1.0 - (-1.0))
    yc = ComplexF64[1 + 0im, 0.25 + 0.5im, -0.5 + 0.25im, 0, 1]
    ec = gardner_ted(yc; sps=2, phase=1)
    mid = yc[2]
    @test ec[1] ≈ real(conj(mid) * (yc[3] - yc[1]))

    @test_throws ArgumentError gardner_ted(y; sps=3, phase=1)
    @test_throws ArgumentError gardner_ted(y; sps=1, phase=1)
    @test_throws ArgumentError gardner_ted(y; sps=2, phase=0)
    @test_throws ArgumentError gardner_ted(y; sps=2, phase=3)
    @test_throws ArgumentError gardner_ted([1.0, 2.0]; sps=2, phase=1)
    @test_throws ArgumentError gardner_ted([1.0, NaN, 3.0, 4.0, 5.0]; sps=2, phase=1)

    # S-curve on a root-raised-cosine matched filter. Late strobes (positive
    # delta) average positive; early strobes average negative.
    sps = 4
    taps = rrc_filter(0.35, 10, sps)
    bits = _mixed_bits(160)
    syms = modulate(bits, BPSK())
    wave = matched_filter(pulse_shape(syms, taps, sps), taps)
    p = length(taps)
    ideal = mod(p - 1, sps)
    function mean_ted(delta)
        d = mod(ideal + delta, sps)
        z = d == 0 ? wave : fractional_delay(wave, d)
        err = gardner_ted(z; sps=sps, phase=1)
        return sum(err) / length(err)
    end
    @test mean_ted(0.0) > mean_ted(-0.5)
    @test mean_ted(0.5) > 0
    @test mean_ted(-0.5) < 0
    @test abs(mean_ted(0.0)) < abs(mean_ted(0.5))

    qwave = matched_filter(pulse_shape(modulate(_mixed_bits(160), QPSK()), taps, sps), taps)
    qideal = mod(length(taps) - 1, sps)
    function qmean(delta)
        d = mod(qideal + delta, sps)
        z = d == 0 ? qwave : fractional_delay(qwave, d)
        err = gardner_ted(z; sps=sps, phase=1)
        return sum(err) / length(err)
    end
    @test qmean(0.5) > 0
    @test qmean(-0.5) < 0
end

@testset "Gardner loop and open loop" begin
    sps = 4
    taps = rrc_filter(0.35, 10, sps)
    n_sym = 180
    alt = [isodd(i) ? 0 : 1 for i in 1:n_sym]
    syms = modulate(alt, BPSK())
    wave = matched_filter(pulse_shape(syms, taps, sps), taps)
    peak0 = length(taps)          # symbol 1
    peak = peak0 + sps            # symbol 2, so a previous symbol exists
    late = gardner_loop(wave; sps=sps, alpha=0.2, mu0=0.35, basepoint=peak, n_symbols=100)
    late_err = [late.basepoint[k] + late.mu[k] - (peak + (k - 1) * sps) for k in eachindex(late.mu)]
    @test late_err[1] ≈ 0.35 atol = 1e-12
    @test sum(late.error[1:20]) / 20 > 0
    @test late_err[end] < late_err[1]
    @test abs(late_err[end]) < 0.08
    @test abs(late_err[end]) < abs(late_err[1])

    early = gardner_loop(wave; sps=sps, alpha=0.2, mu0=0.65, basepoint=peak - 1, n_symbols=100)
    early_err = [early.basepoint[k] + early.mu[k] - (peak + (k - 1) * sps) for k in eachindex(early.mu)]
    @test early_err[1] ≈ -0.35 atol = 1e-12
    @test sum(early.error[1:20]) / 20 < 0
    @test early_err[end] > early_err[1]
    @test abs(early_err[end]) < 0.08

    # Open loop: an extra 0.4 sample delay. The estimate is modulo sps.
    # Two linear interpolators do not reproduce the original symbols; the
    # claim is that the delay matches the true phase, and that the samples
    # match a single interpolation at that true phase.
    extra = 0.4
    shifted = fractional_delay(wave, extra)
    est = gardner_open_loop(shifted; sps=sps, n_frac=8, phase=1)
    truth = mod((peak0 - extra) - 1, sps)
    @test _circ(est.delay, truth, sps) < 0.02
    @test 0 <= est.delay < sps
    @test length(est.grid) == sps * 8
    oracle = downsample_symbols(fractional_delay(shifted, truth), sps; offset=0)
    sampled = downsample_symbols(fractional_delay(shifted, est.delay), sps; offset=0)
    nkeep = min(length(oracle), length(sampled))
    @test maximum(abs, sampled[10:(nkeep - 10)] .- oracle[10:(nkeep - 10)]) < 1e-3
    peak_after = peak0 - extra - est.delay
    sym_index = round(Int, (peak_after - 1) / sps) + 1

    glued = single_carrier_timing(syms, taps, sps; sample_delay=extra, n_frac=8)
    @test _circ(glued.delay, truth, sps) < 0.02
    @test glued.symbol_index == sym_index

    # No injected fractional delay: the interior matches the transmitted symbols.
    # The pulse tail that pulse_shape drops is outside this window.
    plain = single_carrier_timing(syms, taps, sps; n_frac=8)
    @test _circ(plain.delay, 0.0, sps) < 0.02
    window = 30:(n_sym - 30)
    @test maximum(abs, plain.symbols[plain.symbol_index .+ (window .- 1)] .- syms[window]) < 0.03
    @test abs(plain.symbols[plain.symbol_index] - syms[1]) < 0.03

    qsyms = modulate(_mixed_bits(140), QPSK())
    qest = gardner_open_loop(matched_filter(pulse_shape(qsyms, taps, sps), taps); sps=sps, n_frac=8)
    @test _circ(qest.delay, 0.0, sps) < 0.05

    # No transitions: the detector stays at zero and the open loop refuses.
    flat = modulate(zeros(Int, 80), BPSK())
    flat_wave = matched_filter(pulse_shape(flat, taps, sps), taps)
    @test_throws ArgumentError gardner_open_loop(flat_wave; sps=sps, n_frac=4)

    # A gain that steps by more than a symbol must not walk backwards forever.
    sps2 = 4
    wild = zeros(32)
    wild[1] = -10.0
    wild[3] = 10.0
    wild[5] = 10.0
    @test_throws ArgumentError gardner_loop(wild; sps=sps2, alpha=1.0, mu0=0.0, basepoint=5)

    # Default basepoint is the second symbol of a phase-1 grid, and mu = 0
    # takes the sample with no interpolation.
    stopped = gardner_loop(wave; sps=sps, alpha=0.2, mu0=0.0)
    @test stopped.basepoint[1] == sps + 1
    @test stopped.mu[1] == 0
    @test length(stopped.symbols) >= 2
    @test gardner_ted(wave; sps=sps, phase=3) isa Vector{Float64}

    @test_throws ArgumentError gardner_loop(wave; sps=3)
    @test_throws ArgumentError gardner_loop(wave; sps=sps, alpha=0)
    @test_throws ArgumentError gardner_loop(wave; sps=sps, alpha=Inf)
    @test_throws ArgumentError gardner_loop(wave; sps=sps, mu0=1)
    @test_throws ArgumentError gardner_loop(wave; sps=sps, mu0=-0.1)
    @test_throws ArgumentError gardner_loop(wave; sps=sps, basepoint=sps)
    @test_throws ArgumentError gardner_loop(wave; sps=sps, n_symbols=0)
    @test_throws ArgumentError gardner_loop([1.0, 2.0, 3.0]; sps=2, basepoint=5)
    @test_throws ArgumentError gardner_open_loop(wave; sps=sps, n_frac=0)
    @test_throws ArgumentError single_carrier_timing(syms, taps, 3)
    @test_throws ArgumentError single_carrier_timing(syms, taps, sps; sample_delay=-1)
    @test_throws ArgumentError single_carrier_timing(syms, Float64[], sps)
    @test_throws ArgumentError single_carrier_timing(Float64[], taps, sps)
end

@testset "Costas loop" begin
    bbits = _mixed_bits(400)
    bpsk = ComplexF64.(modulate(bbits, BPSK()))
    θ = 0.35
    y = bpsk .* exp(im * θ)
    @test costas_detector(exp(im * θ), BPSK()) > 0
    @test costas_detector(exp(im * -θ), BPSK()) < 0
    @test costas_detector(-1.0 + 0im, BPSK()) ≈ 0 atol = 1e-12
    # A rotation of π is a bit flip and a stable false lock of the BPSK detector.
    @test demodulate(bpsk .* -1, BPSK()) == (1 .- bbits)
    @test abs(costas_detector(-bpsk[1], BPSK())) < 1e-12

    locked = costas_loop(y, BPSK(); alpha=0.08, beta=0.0)
    @test sum(locked.error[1:20]) / 20 > 0
    @test locked.phase[end] > locked.phase[1]
    @test abs(locked.phase[end] - θ) < 0.02
    @test abs(locked.phase[end] - θ) < abs(locked.phase[1] - θ)
    tail = locked.symbols[(end - 40):end]
    @test demodulate(tail, BPSK()) == bbits[(end - 40):end]

    yneg = bpsk .* exp(im * (-0.4))
    neg = costas_loop(yneg, BPSK(); alpha=0.08, beta=0.0)
    @test sum(neg.error[1:20]) / 20 < 0
    @test neg.phase[end] < neg.phase[1]
    @test abs(neg.phase[end] - (-0.4)) < 0.02

    ω = 0.012
    yf = bpsk .* exp.(im * ω .* (0:(length(bpsk) - 1)))
    tracked = costas_loop(yf, BPSK())
    f_early = sum(tracked.frequency[1:30]) / 30
    f_late = sum(tracked.frequency[(end - 80):end]) / 80
    @test f_late > f_early
    @test f_late > 0
    @test abs(f_late - ω) < 0.004
    ang = angle.(tracked.symbols[(end - 80):end] .* conj.(bpsk[(end - 80):end]))
    @test sqrt(sum(abs2, ang) / length(ang)) < 0.05

    qbits = _mixed_bits(800)
    qpsk = modulate(qbits, QPSK())
    θq = 0.22
    yq = qpsk .* exp(im * θq)
    s0 = (1 + 1im) / sqrt(2)
    @test costas_detector(s0 * exp(im * θq), QPSK()) > 0
    @test costas_detector(s0 * exp(im * -θq), QPSK()) < 0
    @test costas_detector(-s0 * exp(im * θq), QPSK()) > 0
    qlock = costas_loop(yq, QPSK(); alpha=0.08, beta=0.0)
    @test sum(qlock.error[1:20]) / 20 > 0
    @test qlock.phase[end] > qlock.phase[1]
    @test abs(qlock.phase[end] - θq) < 0.02
    qtail_syms = qlock.symbols[(end - 30):end]
    qtail_bits = qbits[(end - 61):end]
    @test demodulate(qtail_syms, QPSK()) == qtail_bits

    @test_throws ArgumentError costas_detector(0.2 + 0.1im, PSK(4))
    @test_throws ArgumentError costas_detector(0.2 + 0.1im, QAM(16))
    @test_throws ArgumentError costas_detector(NaN, BPSK())
    @test_throws ArgumentError costas_loop(ComplexF64[], BPSK())
    @test_throws ArgumentError costas_loop(bpsk, BPSK(); alpha=0)
    @test_throws ArgumentError costas_loop(bpsk, BPSK(); beta=-0.1)
    @test_throws ArgumentError costas_loop(bpsk, BPSK(); phase0=Inf)
    @test_throws ArgumentError costas_loop(bpsk, PSK(8))
end

@testset "Schmidl-Cox and integer CFO" begin
    N = 64
    pre = schmidl_preamble(N; n_cp=16)
    @test length(pre.freq) == N
    @test pre.n_cp == 16
    @test length(pre.time) == N + 16
    @test all(iszero, pre.freq[2:2:end])
    @test all(x -> x == 1 || x == -1, pre.freq[1:2:end])
    body = pre.time[(pre.n_cp + 1):end]
    L = N ÷ 2
    @test body[1:L] ≈ body[(L + 1):end] atol = 1e-12
    @test pre.time[1:16] ≈ pre.time[(end - 15):end] atol = 1e-12

    # Noiseless plateau: the earliest maximum is the first preamble sample.
    delay = 11
    pad = vcat(zeros(ComplexF64, delay), pre.time, zeros(ComplexF64, 24))
    sc = schmidl_cox(pad; half=L)
    @test sc.d == delay + 1
    @test sc.epsilon_f ≈ 0 atol = 1e-9
    plateau = sc.metric[sc.d:(sc.d + pre.n_cp)]
    @test minimum(plateau) > 0.99
    @test maximum(abs, plateau .- plateau[1]) < 1e-8

    function tone(x, ε)
        return [x[i] * exp(im * 2π * ε * (i - 1) / N) for i in eachindex(x)]
    end

    for ε in (0.2, -0.35, 0.7)
        y = tone(pre.time, ε)
        hat = schmidl_cox(y; half=L)
        @test hat.d == 1
        @test hat.epsilon_f ≈ ε atol = 1e-8
    end
    pure = schmidl_sync(tone(pre.time, 0.2), pre.freq; n_cp=pre.n_cp)
    @test pure.epsilon_i == 0
    @test pure.epsilon_f ≈ 0.2 atol = 1e-8
    @test pure.d == 1
    @test length(pure.signal) == length(pre.time)

    # Even-integer ambiguity of the half-symbol correlator, then the search.
    # 2.25 = 0.25 + 2, -3.4 = 0.6 + (-4), 1.3 = -0.7 + 2.
    for (ε, εf, εi) in ((2.25, 0.25, 2), (-3.4, 0.6, -4), (1.3, -0.7, 2), (-1.7, 0.3, -2))
        y = vcat(zeros(ComplexF64, 5), tone(pre.time, ε))
        sync = schmidl_sync(y, pre.freq; n_cp=pre.n_cp, z_max=4)
        @test sync.d == 6
        @test sync.sto == sync.d + pre.n_cp
        @test sync.epsilon_f ≈ εf atol = 1e-8
        @test sync.epsilon_i == εi
        @test sync.epsilon ≈ ε atol = 1e-8
        # Correcting the reported total puts the two halves back in phase.
        undone = cfo_correct(y, sync.epsilon, N)
        b = undone[(sync.sto):(sync.sto + N - 1)]
        @test b[1:L] ≈ b[(L + 1):end] atol = 1e-8
    end

    # z_max = 0 cannot see the even integer. The fractional piece is still there.
    buried = schmidl_sync(tone(pre.time, 2.25), pre.freq; n_cp=pre.n_cp, z_max=0)
    @test buried.epsilon_f ≈ 0.25 atol = 1e-8
    @test buried.epsilon_i == 0

    # No prefix: the metric peak is the first useful sample, and a delay shifts it.
    bare = schmidl_preamble(N; n_cp=0)
    ybare = vcat(zeros(ComplexF64, 5), tone(bare.time, 2.3))
    sbare = schmidl_sync(ybare, bare.freq; n_cp=0, z_max=4)
    @test sbare.d == 6
    @test sbare.epsilon ≈ 2.3 atol = 1e-8

    # van de Beek is unchanged: |ε| = 0.7 is outside (-1/2, 1/2) and folds.
    # Schmidl–Cox on its own preamble still reports 0.7.
    params = OFDMParams(N_fft=32, N_cp=8, mod_order=4, n_guard=2, pilot_spacing=4)
    obits = _mixed_bits(params.N_data * params.bits_per_symbol * 4)
    tx, pilots = ofdm_tx(params, obits)
    Nn = params.N_fft
    rxc = [tx[i] * exp(im * 2π * 0.7 * (i - 1) / Nn) for i in eachindex(tx)]
    @test cfo_estimate(rxc, params) ≈ -0.3 atol = 1e-8
    corrected, sto, cfo_hat = ofdm_sync(rxc, params)
    @test sto == 1
    @test cfo_hat ≈ -0.3 atol = 1e-8
    @test length(corrected) == length(rxc)
    @test schmidl_cox(tone(pre.time, 0.7); half=L).epsilon_f ≈ 0.7 atol = 1e-8

    # search_range stops the scan before a late preamble, so the peak is not that start.
    late_only = vcat(zeros(ComplexF64, 30), bare.time)
    full = schmidl_cox(late_only; half=L)
    @test full.d == 31
    capped = schmidl_cox(late_only; half=L, search_range=10)
    @test capped.d <= 10
    @test length(capped.metric) == 10
    @test schmidl_cox(pre.time; half=L, search_range=10^6).d == 1

    @test_throws ArgumentError schmidl_preamble(63)
    @test_throws ArgumentError schmidl_preamble(2)
    @test_throws ArgumentError schmidl_preamble(N; n_cp=-1)
    @test_throws ArgumentError schmidl_preamble(N; n_cp=N + 1)
    @test_throws ArgumentError schmidl_cox(pre.time; half=1)
    @test_throws ArgumentError schmidl_cox(pre.time[1:10]; half=L)
    @test_throws ArgumentError schmidl_cox(pre.time; half=L, search_range=0)
    @test_throws ArgumentError schmidl_cox(vcat(pre.time[1:10], [NaN * im]); half=4)
    @test_throws ArgumentError integer_cfo_estimate(pre.time, pre.freq[1:3]; sto=1)
    @test_throws ArgumentError integer_cfo_estimate(pre.time, pre.freq; sto=0)
    @test_throws ArgumentError integer_cfo_estimate(pre.time, pre.freq; z_max=-1)
    @test_throws ArgumentError integer_cfo_estimate(pre.time, pre.freq; sto=length(pre.time))
    @test_throws ArgumentError schmidl_sync(pre.time, pre.freq[1:5]; n_cp=0)
    @test_throws ArgumentError schmidl_sync(pre.time, pre.freq; n_cp=-1)
end
