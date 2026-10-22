# Single-carrier FIR, matched filter, and symbol-rate sampling.
#
# apply_fir is the tall Toeplitz model used by zf_equalize / mmse_equalize.
# The pulse chain is pulse_shape → sample-rate FIR → matched_filter →
# downsample_symbols. A one-sample pulse makes that cascade exactly
# symbol-spaced, so ZF/MMSE recover the symbols with no noise.
# Root-raised-cosine with no multipath is only approximately Nyquist:
# pulse_shape drops its tail, and the taps themselves are truncated.

@testset "FIR convolution" begin
    h = [1.0, 0.5, -0.25]
    syms = [1.0, -1.0, 0.5, 0.25]
    rx = apply_fir(syms, h)
    @test length(rx) == length(syms) + length(h) - 1
    @test eltype(rx) == Float64
    manual = zeros(length(rx))
    for i in eachindex(rx)
        for j in eachindex(h)
            src = i - j + 1
            if 1 <= src <= length(syms)
                manual[i] += syms[src] * h[j]
            end
        end
    end
    @test rx ≈ manual atol = 1e-12
    @test zf_equalize(rx, h) ≈ syms atol = 1e-10
    @test mmse_equalize(rx, h, 0.0) ≈ syms atol = 1e-10

    stored = copy(syms)
    apply_fir(syms, h)
    @test syms == stored

    hc = [1.0 + 0.2im, 0.4 - 0.1im, -0.3 + 0.05im]
    qsyms = modulate([0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0], QPSK())
    rxc = apply_fir(qsyms, hc)
    @test eltype(rxc) <: Complex
    @test zf_equalize(rxc, hc) ≈ qsyms atol = 1e-8
    @test mmse_equalize(rxc, hc, 0.0) ≈ qsyms atol = 1e-8

    # Real symbols through complex taps stay on the complex equaliser model.
    mixed = apply_fir([1.0, -1.0, 1.0], [0.5 + 0.0im, 0.2im])
    @test eltype(mixed) <: Complex
    @test zf_equalize(mixed, [0.5 + 0.0im, 0.2im]) ≈ [1.0, -1.0, 1.0] atol = 1e-10

    @test apply_fir([1.0, -1.0], [1.0]) ≈ [1.0, -1.0]
    @test_throws ArgumentError apply_fir(Float64[], [1.0])
    @test_throws ArgumentError apply_fir([1.0], Float64[])
end

@testset "FIR noise" begin
    h = [1.0, -0.4]
    syms = ones(8)
    clean, nv_inf = fir_channel(syms, h, Inf)
    @test clean ≈ apply_fir(syms, h)
    @test nv_inf == 0.0

    Random.seed!(1)
    before = randn()
    Random.seed!(1)
    fir_channel(syms, h, Inf; energy_per_bit=1.0)
    @test randn() == before

    snr_db = 3.0
    eb = 2.0
    Random.seed!(11)
    rx, nv = fir_channel(syms, h, snr_db; energy_per_bit=eb)
    Random.seed!(11)
    @test rx ≈ awgn(apply_fir(syms, h), snr_db; energy_per_bit=eb)
    ρ = snr_db_to_linear(snr_db)
    @test nv ≈ eb / (2 * ρ)

    hc = [1.0 + 0im, 0.3]
    Random.seed!(12)
    rxc, nvc = fir_channel(ComplexF64.(syms), hc, snr_db; energy_per_bit=eb)
    @test nvc ≈ eb / ρ
    @test eltype(rxc) <: Complex
    Random.seed!(12)
    @test rxc ≈ awgn(apply_fir(ComplexF64.(syms), hc), snr_db; energy_per_bit=eb)

    # Positional Eb/N0 assumes Es = 1. Real BPSK at high SNR is error-free
    # on this channel, and the returned variance is the MMSE input.
    Random.seed!(0x0c0de1162)
    bits = rand(0:1, 64)
    bpsk = modulate(bits, BPSK())
    y, n_bpsk = fir_channel(bpsk, h, 40.0, 1)
    @test n_bpsk ≈ 1.0 / (2 * snr_db_to_linear(40.0))
    @test ber(bits, demodulate(mmse_equalize(y, h, n_bpsk), BPSK())) == 0.0
    @test ber(bits, demodulate(zf_equalize(y, h), BPSK())) == 0.0

    @test_throws ArgumentError fir_channel(syms, h, 0.0; energy_per_bit=0.0)
    @test_throws ArgumentError fir_channel(syms, h, 0.0; energy_per_bit=-1.0)
    @test_throws ArgumentError fir_channel(syms, h, 0.0, 0)
end

@testset "Matched filter and downsampling" begin
    g = [0.2, -0.4 + 0.1im, 0.7]
    x = [1.0 + 0.0im, 0.0, -0.5]
    @test matched_filter(x, g) ≈ apply_fir(x, conj.(reverse(g))) atol = 1e-12

    rrc = rrc_filter(0.35, 6, 4)
    @test rrc ≈ reverse(rrc) atol = 1e-12
    # An impulse comes back as the pulse. The pulse peak is the centre tap.
    mf_delta = matched_filter([1.0], rrc)
    @test mf_delta ≈ rrc atol = 1e-12
    @test argmax(abs.(mf_delta)) == (length(rrc) + 1) ÷ 2
    # Zero lag of the autocorrelation is index length(taps), and equals the energy.
    ac = matched_filter(rrc, rrc)
    @test argmax(abs.(ac)) == length(rrc)
    @test ac[length(rrc)] ≈ sum(abs2, rrc) atol = 1e-12
    @test ac[length(rrc)] ≈ 1 atol = 1e-12

    @test downsample_symbols([1, 0, -1, 0, 1, 0], 2; offset=0) == [1, -1, 1]
    @test eltype(downsample_symbols([1, 0, -1, 0], 2)) == Int
    @test downsample_symbols(1:7, 3; offset=0) == [1, 4, 7]
    # offset is not taken modulo sps: 5 skips one symbol and keeps phase 2.
    @test downsample_symbols(collect(1:10), 3; offset=5) == [6, 9]
    @test downsample_symbols([10, 20, 30], 1; offset=1) == [20, 30]
    @test downsample_symbols(ComplexF64[1, 2, 3, 4], 2; offset=1) == ComplexF64[2, 4]

    @test_throws ArgumentError downsample_symbols([1.0, 2.0], 0)
    @test_throws ArgumentError downsample_symbols([1.0, 2.0], 2; offset=-1)
    @test_throws ArgumentError downsample_symbols([1.0, 2.0], 2; offset=2)
    @test_throws ArgumentError matched_filter(Float64[], [1.0])
    @test_throws ArgumentError matched_filter([1.0], Float64[])
end

@testset "Pulse, FIR, matched filter, equalizer" begin
    # Causal one-sample pulse at tap 2. pulse_shape's discarded tail is zero,
    # and a multipath delay of one symbol stays on the sampling grid.
    sps = 2
    g = [0.0, 1.0]
    h_sample = [1.0 + 0.0im, 0.0 + 0.0im, 0.25 - 0.5im]
    h_sym = [1.0 + 0.0im, 0.25 - 0.5im]
    bits = [0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1]
    syms = modulate(bits, QPSK())

    shaped = pulse_shape(syms, g, sps)
    faded = apply_fir(shaped, h_sample)
    filtered = matched_filter(faded, g)
    # Transmit taps [0, 1] place each symbol on the next sample.
    # The matched taps are [1, 0], which add no further delay.
    sampled = downsample_symbols(filtered, sps; offset=1)

    @test sampled ≈ apply_fir(syms, h_sym) atol = 1e-12
    zf = zf_equalize(sampled, h_sym)
    mmse = mmse_equalize(sampled, h_sym, 0.0)
    @test zf ≈ syms atol = 1e-8
    @test mmse ≈ syms atol = 1e-8
    @test demodulate(zf, QPSK()) == bits
    @test demodulate(mmse, QPSK()) == bits

    # No multipath and no noise: RRC sampled on the autocorrelation peak.
    β, span, sps_rc = 0.25, 16, 8
    taps = rrc_filter(β, span, sps_rc)
    n_sym = 64
    isolated = zeros(n_sym)
    isolated[32] = 1.0
    burst = Float64[isodd(i) ? 1.0 : -1.0 for i in 1:n_sym]
    offset = length(taps) - 1

    function sample_rrc(symbols)
        wave = apply_fir(pulse_shape(symbols, taps, sps_rc), [1.0])
        mf = matched_filter(wave, taps)
        return downsample_symbols(mf, sps_rc; offset=offset)
    end

    isolated_rx = sample_rrc(isolated)
    @test isolated_rx[32] ≈ 1 atol = 1e-8
    for lag in 1:4
        @test abs(isolated_rx[32 + lag]) < 2e-3
        @test abs(isolated_rx[32 - lag]) < 2e-3
    end

    # Full autocorrelation support of these symbols lies inside the burst,
    # so the transmit truncation does not reach them.
    dense = sample_rrc(burst)
    window = (span + 1):(n_sym - span)
    @test maximum(abs, dense[window] .- burst[window]) < 2e-3
    # pulse_shape has already dropped the tail, so the last symbol is not recovered.
    @test abs(dense[n_sym] - burst[n_sym]) > 0.5
end
