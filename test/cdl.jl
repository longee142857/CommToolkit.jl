# Teaching CDL (TR 38.901 CDL-A/B/C delay–power profiles). Table values are
# frozen against TR 38.901 V16.1.0 Tables 7.7.1-1/2/3 (normalized delay and
# power columns). The only closed form is the fully merged flat case, which
# degenerates to rayleigh_channel; the dispersive checks are power/shape
# statistics and a perfect-recovery OFDM smoke test. No golden BER curve.

@testset "CDL profiles" begin
    pa = cdl_profile(:A)
    pb = cdl_profile(:B)
    pc = cdl_profile(:C)
    @test pa.name == :A && pb.name == :B && pc.name == :C
    @test length(pa.delays) == 23
    @test length(pb.delays) == 23
    @test length(pc.delays) == 24

    # Frozen corners of TR 38.901 V16.1.0 Tables 7.7.1-1/2/3.
    @test pa.delays[1] == 0.0 && pa.powers_db[1] == -13.4
    @test pa.delays[2] == 0.3819 && pa.powers_db[2] == 0.0
    @test pa.delays[23] == 9.6586 && pa.powers_db[23] == -29.7
    @test pa.delays[10] == 1.5375 && pa.powers_db[10] == -15.9
    @test pb.delays[1] == 0.0 && pb.powers_db[1] == 0.0
    @test pb.delays[23] == 4.7834 && pb.powers_db[23] == -11.3
    @test pb.delays[16] == 1.7842 && pb.powers_db[16] == -1.9
    @test pc.delays[1] == 0.0 && pc.powers_db[1] == -4.4
    @test pc.delays[6] == 0.6366 && pc.powers_db[6] == 0.0
    @test pc.delays[24] == 8.6523 && pc.powers_db[24] == -22.8

    @test all(>=(0), pa.delays) && all(>=(0), pb.delays) && all(>=(0), pc.delays)

    # The tabulated normalized delays have unit power-weighted RMS delay
    # spread (TR 38.901 §7.7.3), so delay_spread scales them linearly.
    for p in (pa, pb, pc)
        w = exp10.(p.powers_db ./ 10)
        w ./= sum(w)
        @test CommToolkit.Channel._rms_delay_spread(p.delays, w) ≈ 1.0 atol = 2e-3
    end

    # cdl_profile returns copies: mutating one does not corrupt the table.
    pa.delays[1] = 99.0
    @test cdl_profile(:A).delays[1] == 0.0

    @test_throws ArgumentError cdl_profile(:D)
    @test_throws ArgumentError cdl_profile(:E)
    @test_throws ArgumentError cdl_profile(:Z)

    # Custom profiles keep caller order and validate.
    pcust = CDLProfile([0.0, 0.5, 0.25], [0.0, -3.0, -10.0])
    @test pcust.name == :custom
    @test pcust.delays == [0.0, 0.5, 0.25]
    named = CDLProfile([0.0], [0.0]; name=:indoor)
    @test named.name == :indoor
    @test_throws ArgumentError CDLProfile([0.0, 1.0], [0.0])
    @test_throws ArgumentError CDLProfile(Float64[], Float64[])
    @test_throws ArgumentError CDLProfile([0.0, -1.0], [0.0, 0.0])
    @test_throws ArgumentError CDLProfile([0.0, NaN], [0.0, 0.0])
    @test_throws ArgumentError CDLProfile([0.0], [Inf])
end

@testset "CDL tap quantization and merging" begin
    pa = cdl_profile(:A)

    # CDL-A, fs = 100 MHz (Ts = 10 ns), nominal DS = 100 ns: frozen merged taps.
    taps = cdl_taps(pa, 4, 1e-8; delay_spread=100e-9)
    @test taps.delays == [0, 4, 5, 6, 7, 8, 15, 19, 22, 25, 31, 41, 45, 46, 48, 50, 53, 97]
    @test size(taps.gains) == (18, 1)   # fd = 0: one static column
    @test sum(taps.powers) ≈ 1.0
    @test issorted(taps.delays) && allunique(taps.delays)

    # Merged power equals the sum of the colliding clusters. Tap 4 holds
    # clusters 2 and 3 (0.3819 and 0.4025 both quantize to 4 at fs = 100 MHz).
    w = exp10.(pa.powers_db ./ 10)
    w ./= sum(w)
    @test taps.powers[2] ≈ w[2] + w[3] rtol = 1e-12
    @test taps.powers[1] ≈ w[1] rtol = 1e-12
    @test taps.powers[end] ≈ w[23] rtol = 1e-12

    # delay_spread = 1.0 reads the normalized delays as seconds: Ts = 0.5 s
    # quantizes CDL-B to round.(2τ).
    tb = cdl_taps(cdl_profile(:B), 4, 0.5; delay_spread=1.0)
    @test tb.delays == [0, 1, 2, 3, 4, 6, 7, 8, 9, 10]
    @test tb.delays == sort(unique(round.(Int, 2 .* cdl_profile(:B).delays)))

    # Slow sampling merges the whole profile onto tap 0: a flat channel.
    flat = cdl_taps(pa, 8, 1.0; delay_spread=100e-9)
    @test flat.delays == [0]
    @test flat.powers == [1.0]
    @test size(flat.gains) == (1, 1)

    # Thinning keeps the strongest merged taps and renormalizes.
    thin = cdl_taps(pa, 4, 1e-8; delay_spread=100e-9, n_taps=3)
    @test thin.delays == [4, 5, 6]
    @test sum(thin.powers) ≈ 1.0
    @test thin.powers ≈ [0.642995, 0.161513, 0.195492] rtol = 1e-4
    @test thin.powers ≈ taps.powers[2:4] ./ sum(taps.powers[2:4]) rtol = 1e-12
    # Ties in power break toward the smaller delay.
    tie = CDLProfile([0.3, 0.1, 0.2], [0.0, 0.0, 0.0])
    tied = cdl_taps(tie, 4, 0.05; delay_spread=1.0, n_taps=2)
    @test tied.delays == [2, 4]

    # n = 0 returns empty gains and does not draw from the global stream.
    Random.seed!(1)
    draw = randn()
    Random.seed!(1)
    empty_taps = cdl_taps(pa, 0, 1e-8)
    @test size(empty_taps.gains) == (18, 0)
    @test randn() == draw

    # fd = 0: one static column, same per-tap draws as jakes_fading(1, 0, Ts).
    Random.seed!(2)
    st = cdl_taps(pa, 64, 1e-8; delay_spread=100e-9)
    Random.seed!(2)
    href = [only(jakes_fading(1, 0.0, 1e-8)) for _ in 1:18]
    @test st.gains[:, 1] ≈ sqrt.(st.powers) .* href
    Random.seed!(2)
    @test cdl_taps(pa, 64, 1e-8; delay_spread=100e-9).gains == st.gains

    # fd > 0: n columns, reproducible, ergodic power per tap.
    Random.seed!(RNG_SEED)
    tv = cdl_taps(pa, 8192, 1e-8; delay_spread=100e-9, fd=2e6)   # fd*Ts = 0.02
    @test size(tv.gains) == (18, 8192)
    Random.seed!(RNG_SEED)
    @test cdl_taps(pa, 8192, 1e-8; delay_spread=100e-9, fd=2e6).gains == tv.gains
    means = [sum(abs2, tv.gains[i, :]) / 8192 for i in 1:18]
    # Single-record ergodic estimates; the flat marginal statistics of the
    # same generator are pinned down in the jakes_fading tests.
    @test sum(means) ≈ 1.0 rtol = 0.05
    for i in 1:18
        tv.powers[i] >= 0.1 || continue
        @test means[i] ≈ tv.powers[i] rtol = 0.25
    end

    @test_throws ArgumentError cdl_taps(pa, -1, 1e-8)
    @test_throws ArgumentError cdl_taps(pa, 4, 0.0)
    @test_throws ArgumentError cdl_taps(pa, 4, -1e-8)
    @test_throws ArgumentError cdl_taps(pa, 4, Inf)
    @test_throws ArgumentError cdl_taps(pa, 4, 1e-8; delay_spread=0.0)
    @test_throws ArgumentError cdl_taps(pa, 4, 1e-8; delay_spread=-1e-9)
    @test_throws ArgumentError cdl_taps(pa, 4, 1e-8; delay_spread=Inf)
    @test_throws ArgumentError cdl_taps(pa, 4, 1e-8; fd=-1.0)
    @test_throws ArgumentError cdl_taps(pa, 4, 1e-8; fd=6e7)   # fd*Ts = 0.6 ≥ 1/2
    @test_throws ArgumentError cdl_taps(pa, 4, 1e-8; n_taps=0)
end

@testset "apply_cdl" begin
    # Hand-computed time-varying convolution: gains indexed by input sample.
    taps = CDLTaps([0, 1],
                   ComplexF64[1 2 3; 10 20 30],
                   [0.5, 0.5])
    x = ComplexF64[1, -1, 2]
    y = apply_cdl(x, taps)
    @test y ≈ ComplexF64[1 * 1, 2 * -1 + 10 * 1, 3 * 2 + 20 * -1, 30 * 2]

    # A static single column is the FIR / TDL convolution with the same taps.
    Random.seed!(3)
    st = cdl_taps(cdl_profile(:A), 16, 1e-8; delay_spread=100e-9)
    xs = randn(50) .+ im .* randn(50)
    ycdl = apply_cdl(xs, st)
    fir_taps = zeros(ComplexF64, maximum(st.delays) + 1)
    for (d, g) in zip(st.delays, st.gains[:, 1])
        fir_taps[d + 1] = g
    end
    @test ycdl ≈ apply_fir(xs, fir_taps)
    tdl = TDLChannel(st.gains[:, 1], st.delays)
    @test ycdl ≈ CommToolkit.OFDM.apply_multipath(xs, tdl)

    # Real input promotes to complex.
    @test apply_cdl(ones(4), CDLTaps([0], reshape([2.0 + 0im], 1, 1), [1.0])) == fill(2.0 + 0im, 4)

    @test_throws ArgumentError apply_cdl(ComplexF64[], st)
    wrong = CDLTaps([0], reshape([1.0 + 0im, 2.0 + 0im], 1, 2), [1.0])
    @test_throws ArgumentError apply_cdl(ones(3), wrong)

    # CDLTaps field validation.
    @test_throws ArgumentError CDLTaps([1, 0], zeros(ComplexF64, 2, 1), [0.5, 0.5])
    @test_throws ArgumentError CDLTaps([0, 0], zeros(ComplexF64, 2, 1), [0.5, 0.5])
    @test_throws ArgumentError CDLTaps([-1], zeros(ComplexF64, 1, 1), [1.0])
    @test_throws ArgumentError CDLTaps([0], zeros(ComplexF64, 1, 1), [1.0, 0.0])
    @test_throws ArgumentError CDLTaps([0], zeros(ComplexF64, 2, 1), [1.0])
    @test_throws ArgumentError CDLTaps(Int[], zeros(ComplexF64, 0, 1), Float64[])
    @test_throws ArgumentError CDLTaps([0], zeros(ComplexF64, 1, 1), [-1.0])
end

@testset "cdl_channel" begin
    Random.seed!(RNG_SEED)
    syms = modulate(rand(0:1, 2 * 400), QPSK())
    rx, taps = cdl_channel(syms, Inf, 2; Ts=1e-8, profile=cdl_profile(:A))
    @test length(rx) == length(syms) + maximum(taps.delays)
    @test rx ≈ apply_cdl(syms, taps)

    # The same seed reproduces the same bits, taps, and observation.
    Random.seed!(RNG_SEED)
    syms2 = modulate(rand(0:1, 2 * 400), QPSK())
    rx2, taps2 = cdl_channel(syms2, Inf, 2; Ts=1e-8, profile=cdl_profile(:A))
    @test taps2.gains == taps.gains
    @test rx2 == rx

    # Empty input: no output samples, and cdl_taps drew zero gain columns.
    rx_e, taps_e = cdl_channel(ComplexF64[], 0.0, 1; Ts=1e-8)
    @test isempty(rx_e) && size(taps_e.gains, 2) == 0
    @test_throws ArgumentError cdl_channel(syms, 0.0, 0; Ts=1e-8)
    @test_throws UndefKeywordError cdl_channel(syms, 0.0, 2)

    # Interior output power of one static drop tracks the realized tap power
    # sum |g_i|² (unit-energy symbols, cross terms average out). The drop
    # average of that sum is 1 because sum p_i = 1 and E|h_i|² = 1; the
    # ergodic check of that statement is in the quantization testset.
    Random.seed!(RNG_SEED)
    n = 200_000
    xs = modulate(rand(0:1, n), BPSK())
    tD = cdl_taps(cdl_profile(:A), n, 1e-8; delay_spread=100e-9)
    yD = apply_cdl(xs, tD)
    dmax = maximum(tD.delays)
    interior = yD[(dmax + 1):n]
    @test sum(abs2, interior) / length(interior) ≈ sum(abs2, tD.gains[:, 1]) rtol = 0.02
end

@testset "CDL static drop on the OFDM chain" begin
    Random.seed!(RNG_SEED)
    params = OFDMParams(N_fft=32, N_cp=8, mod_order=4, n_guard=2, pilot_spacing=4)
    bits = rand(0:1, params.N_data * params.bits_per_symbol * 6)
    tx, pilots = ofdm_tx(params, bits)

    # fs = 25 MHz (Ts = 40 ns), DS = 30 ns: taps resolve within the prefix.
    drop = cdl_taps(cdl_profile(:A), 1, 40e-9; delay_spread=30e-9)
    @test drop.delays == [0, 1, 2, 3, 4, 7]
    @test maximum(drop.delays) <= params.N_cp
    tdl = TDLChannel(drop.gains[:, 1], drop.delays)
    rx, H = apply_channel(tx, tdl, Inf, params)
    @test ber(bits, ofdm_rx_known(rx, params, H, pilots)[1:length(bits)]) == 0.0

    # Same waveform through apply_cdl matches the TDL path before truncation.
    n_sym = length(tx) ÷ (params.N_fft + params.N_cp)
    keep = n_sym * (params.N_fft + params.N_cp)
    @test apply_cdl(tx, drop)[1:keep] ≈ CommToolkit.OFDM.apply_multipath(tx, tdl)[1:keep]
end

@testset "CDL flat degenerate is flat Rayleigh" begin
    # delay_spread far below one sample: every cluster merges onto tap 0, so
    # the composite coefficient is CN(0, 1) and the closed form of
    # rayleigh_channel applies. This is the only BER formula claimed here.
    γb = 10^(8 / 10)
    target = 300
    k = 0
    n = 0
    Random.seed!(RNG_SEED)
    while k < target && n < 200_000
        bit = rand(0:1)
        s = modulate([bit], BPSK())
        rx, taps = cdl_channel(s, 8.0, 1; Ts=1e-8, delay_spread=1e-12)
        @assert length(taps.delays) == 1
        zf = rx[1] / taps.gains[1, 1]
        k += demodulate([zf], BPSK())[1] != bit
        n += 1
    end
    @test n < 200_000
    assert_theory_inside(k, n, ber_bpsk_rayleigh(γb); min_errors=300)

    # Correlated flat record: the conditional-BER average tracks Rayleigh,
    # mirroring the jakes_fading marginal test (no i.i.d. interval here).
    Random.seed!(RNG_SEED)
    tv = cdl_taps(cdl_profile(:C), 80_000, 1.0; delay_spread=1e-12, fd=0.05)
    @test length(tv.delays) == 1
    g = tv.gains[1, :]
    pb = sum(qfunc(sqrt(2 * abs2(gk) * γb)) for gk in g) / length(g)
    @test pb ≈ ber_bpsk_rayleigh(γb) atol = 0.006
    @test sum(abs2, g) / length(g) ≈ 1 rtol = 0.04
end
