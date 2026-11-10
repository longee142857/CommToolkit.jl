# Textbook noise quantities, not a transcription of the implementation.
#
# Es = 1 for modulate. k bits per symbol, Eb = Es/k.
# Linear ρ = Eb/N0. N0 = Eb/ρ. One real dimension has variance σ² = N0/2.
# A complex symbol has two independent dimensions, so E[|n|²] = N0.
# A real symbol has one dimension, so E[n²] = σ² = N0/2.
#
# Antipodal likelihood, symbols ±A, noise variance σ² on that axis:
#   L = log N(y; +A, σ²) / N(y; -A, σ²) = 2 A y / σ².
# BPSK uses A = 1 on the real axis. Gray QPSK uses A = 1/√2 on each axis
# (00 → (+1+j)/√2, and bit 0 is the sign of the real part).
#
# The terminated (7, 5) code is not measured here. conv_encode maps n
# information bits to 2(n+2) coded bits, effective rate n/(2(n+2)).

@testset "Noise contract: power and soft LLR" begin
    # Six-standard-error band. Sample second moment of N(0, σ²) has
    # standard error σ² √(2/n). Complex power |n|² is exponential with
    # mean N0, standard error N0/√n.
    real_ok(s2, σ2, n) = abs(s2 - σ2) <= 6 * σ2 * sqrt(2 / n)
    complex_power_ok(p, N0, n) = abs(p - N0) <= 6 * N0 / sqrt(n)

    ρ0 = exp10(0.0)
    # BPSK, Eb = 1, Eb/N0 = 0 dB → N0 = 1, σ² = 1/2.
    N0_bpsk = 1 / ρ0
    σ2_bpsk = N0_bpsk / 2
    nreal = 16_384
    Random.seed!(11)
    bits = rand(0:1, nreal)
    x = modulate(bits, BPSK())
    y = awgn(x, 0.0, 1)
    e = y .- x
    @test real_ok(sum(abs2, e) / nreal, σ2_bpsk, nreal)
    @test demodulate_soft(y, BPSK(), σ2_bpsk) ≈ (2 .* y ./ σ2_bpsk) atol = 1e-8

    # Same Eb/N0 on a complex BPSK vector: two dimensions, total power N0.
    Random.seed!(12)
    xc = ComplexF64.(modulate(rand(0:1, nreal), BPSK()))
    yc = awgn(xc, 0.0, 1)
    ec = yc .- xc
    @test complex_power_ok(sum(abs2, ec) / nreal, N0_bpsk, nreal)
    @test real_ok(sum(abs2, real.(ec)) / nreal, σ2_bpsk, nreal)
    @test real_ok(sum(abs2, imag.(ec)) / nreal, σ2_bpsk, nreal)
    yhand = 0.2 + 0.9im
    @test demodulate_soft([yhand], BPSK(), σ2_bpsk) ≈ [2 * real(yhand) / σ2_bpsk] atol = 1e-8

    # QPSK, Eb = 1/2, Eb/N0 = 0 dB → N0 = 1/2, σ² = 1/4, A = 1/√2.
    N0_qpsk = (1 / 2) / ρ0
    σ2_qpsk = N0_qpsk / 2
    nq = 8_192
    Random.seed!(13)
    xq = modulate(rand(0:1, 2 * nq), QPSK())
    yq = awgn(xq, 0.0, 2)
    eq = yq .- xq
    @test complex_power_ok(sum(abs2, eq) / nq, N0_qpsk, nq)
    @test real_ok(sum(abs2, real.(eq)) / nq, σ2_qpsk, nq)
    @test real_ok(sum(abs2, imag.(eq)) / nq, σ2_qpsk, nq)
    yq1 = ComplexF64(0.3 + 0.4im)
    A = inv(sqrt(2))
    @test demodulate_soft([yq1], QPSK(), σ2_qpsk) ≈
          [2 * A * real(yq1) / σ2_qpsk, 2 * A * imag(yq1) / σ2_qpsk] atol = 1e-8

    # MIMO: pass the textbook N0, not N0/2. Identity channel, unit-energy BPSK.
    N0_mimo = 1 / exp10(5 / 10)
    H = ComplexF64[1 0; 0 1]
    xs = ComplexF64[1.0, -1.0]
    Random.seed!(14)
    ntrial = 4_000
    psum = rsum = isum = 0.0
    for _ in 1:ntrial
        err = apply_mimo(xs, H, N0_mimo) .- xs
        psum += sum(abs2, err)
        rsum += sum(abs2, real.(err))
        isum += sum(abs2, imag.(err))
    end
    nent = ntrial * length(xs)
    @test complex_power_ok(psum / nent, N0_mimo, nent)
    @test real_ok(rsum / nent, N0_mimo / 2, nent)
    @test real_ok(isum / nent, N0_mimo / 2, nent)
    Random.seed!(15)
    y1 = apply_mimo(ComplexF64[1.0], ones(ComplexF64, 1, 1), N0_mimo)
    @test demodulate_soft(y1, BPSK(), N0_mimo / 2) ≈ [2 * real(y1[1]) / (N0_mimo / 2)] atol = 1e-8

    # FIR. Real entry variance is σ² = N0/2. Complex entry variance is N0.
    nt = 16_384
    Random.seed!(16)
    xr = modulate(rand(0:1, nt), BPSK())
    rxr, nvr = fir_channel(xr, [1.0], 0.0, 1)
    @test nvr == σ2_bpsk
    @test real_ok(sum(abs2, rxr .- xr) / nt, σ2_bpsk, nt)
    Random.seed!(17)
    xcf = ComplexF64.(modulate(rand(0:1, nt), BPSK()))
    rxc, nvc = fir_channel(xcf, ComplexF64[1.0], 0.0, 1)
    @test nvc == N0_bpsk
    @test complex_power_ok(sum(abs2, rxc .- xcf) / nt, N0_bpsk, nt)
    @test real_ok(sum(abs2, real.(rxc .- xcf)) / nt, σ2_bpsk, nt)
    @test real_ok(sum(abs2, imag.(rxc .- xcf)) / nt, σ2_bpsk, nt)

    # OFDM sample SNR: E[|n|²] = P_sample / 10^(snr_db/10), each dimension half.
    # That total is not the Eb/N0 reading N0 = Eb/ρ with Eb taken from Es = 1.
    params = OFDMParams(N_fft=32, N_cp=8, mod_order=4, n_guard=2, pilot_spacing=4)
    obits = zeros(Int, params.N_data * params.bits_per_symbol * 20)
    tx, _ = ofdm_tx(params, obits)
    flat = TDLChannel(delays_samples=[0], gain_dB=[0.0])
    clean, _ = apply_channel(tx, flat, Inf, params)
    Random.seed!(18)
    noisy, _ = apply_channel(tx, flat, 10.0, params)
    onoise = noisy .- clean
    non = length(onoise)
    P = sum(abs2, clean) / non
    snr_lin = exp10(1.0)
    Nsamp = P / snr_lin
    @test complex_power_ok(sum(abs2, onoise) / non, Nsamp, non)
    @test real_ok(sum(abs2, real.(onoise)) / non, Nsamp / 2, non)
    @test real_ok(sum(abs2, imag.(onoise)) / non, Nsamp / 2, non)
    eb_style = (1 / 2) / snr_lin
    @test abs(Nsamp - eb_style) / eb_style > 0.5
end

@testset "Infinite SNR bypasses additive noise only" begin
    function next_randn(f; seed=0x0c0de1162)
        Random.seed!(seed)
        f()
        return randn()
    end
    baseline = next_randn(() -> nothing)

    y = awgn(Float64[1.0, -1.0, 0.25], Inf, 1)
    @test y == Float64[1.0, -1.0, 0.25]
    @test y isa Vector{Float64}
    yc = awgn(ComplexF64[1.0, -1.0im], Inf; energy_per_bit=1.0)
    @test yc == ComplexF64[1.0, -1.0im]
    @test yc isa Vector{ComplexF64}
    ye = awgn(Float64[], Inf, 1)
    @test ye isa Vector{Float64} && isempty(ye)
    yce = awgn(ComplexF64[], Inf)
    @test yce isa Vector{ComplexF64} && isempty(yce)
    @test next_randn(() -> awgn(Float64[1.0, -1.0], Inf, 1)) == baseline
    @test next_randn(() -> awgn(ComplexF64[1.0, -1.0], Inf, 1)) == baseline
    @test next_randn(() -> awgn(Float64[], Inf, 1)) == baseline
    @test next_randn(() -> awgn(Float64[1.0, -1.0], 0.0, 1)) != baseline

    @test_throws ArgumentError awgn([1.0], Inf; energy_per_bit=0)
    @test_throws ArgumentError awgn([1.0], Inf; energy_per_bit=-2)
    @test_throws ArgumentError awgn([1.0, -1.0], Inf, 0)
    @test next_randn(() -> (try awgn([1.0], Inf; energy_per_bit=0) catch end)) == baseline
    @test next_randn(() -> (try awgn([1.0], Inf, 0) catch end)) == baseline

    x = ones(6)
    Random.seed!(3)
    rx, h = rayleigh_channel(x, Inf, 1)
    @test rx ≈ ComplexF64.(x) .* h
    @test eltype(rx) <: Complex
    Random.seed!(3)
    (randn(6) .+ im .* randn(6)) ./ sqrt(2)
    tail_h = randn()
    @test next_randn(() -> rayleigh_channel(x, Inf, 1); seed=3) == tail_h
    @test next_randn(() -> rayleigh_channel(x, 0.0, 1); seed=3) != tail_h

    Random.seed!(4)
    jrx, jh = jakes_channel(ones(8), Inf, 1, 10.0, 1e-3)
    @test jrx ≈ complex.(ones(8)) .* jh
    Random.seed!(4)
    jakes_fading(8, 10.0, 1e-3)
    tail_j = randn()
    @test next_randn(() -> jakes_channel(ones(8), Inf, 1, 10.0, 1e-3); seed=4) == tail_j

    syms = ones(ComplexF64, 16)
    Random.seed!(5)
    crx, ct = cdl_channel(syms, Inf, 1; Ts=1e-8, profile=cdl_profile(:A))
    @test crx ≈ apply_cdl(syms, ct)
    Random.seed!(5)
    cdl_taps(cdl_profile(:A), length(syms), 1e-8)
    tail_c = randn()
    @test next_randn(() -> cdl_channel(syms, Inf, 1; Ts=1e-8, profile=cdl_profile(:A)); seed=5) == tail_c
    @test next_randn(() -> (try
        cdl_channel(ones(4), Inf, 0; Ts=1e-8)
    catch end)) == baseline

    params = OFDMParams(N_fft=32, N_cp=8, mod_order=4)
    obits = zeros(Int, params.N_data * params.bits_per_symbol * 4)
    tx, _ = ofdm_tx(params, obits)
    flat = TDLChannel(delays_samples=[0], gain_dB=[0.0])
    clean, _ = apply_channel(tx, flat, Inf, params)
    @test clean isa Vector{ComplexF64}
    @test length(clean) == length(tx)
    @test next_randn(() -> apply_channel(tx, flat, Inf, params)) == baseline
    @test next_randn(() -> apply_channel(tx, flat, 0.0, params)) != baseline

    @test_throws ArgumentError fir_channel([1.0, -1.0], [1.0], Inf; energy_per_bit=-1)
    @test next_randn(() -> (try fir_channel([1.0], [1.0], Inf; energy_per_bit=0) catch end)) == baseline
end

@testset "Hard Viterbi rejects non-bits; float [0, 1] stays an LLR" begin
    info = zeros(Int, 16)
    coded = conv_encode(info)
    @test viterbi_decode(coded) == info
    # 1.0 favours bit 0. It is not a hard 1. The all-zero word is the coded image of these zeros.
    @test viterbi_decode(fill(1.0, length(coded))) == info
    @test viterbi_decode(fill(0.25, length(coded))) == info
    soft01 = Float64.(coded)
    @test length(viterbi_decode(soft01)) == length(info)
    erased = zeros(length(coded))
    @test length(viterbi_decode(erased)) == length(info)

    long_bad = fill(0, length(coded))
    long_bad[end] = 2
    err = try
        viterbi_decode(long_bad)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("0 or 1", err.msg)

    short_bad = try
        viterbi_decode(Int[2, 2])
        nothing
    catch e
        e
    end
    @test short_bad isa ArgumentError
    @test occursin("termination", short_bad.msg)
    @test !occursin("0 or 1", short_bad.msg)

    empty_i = try
        viterbi_decode(Int[])
        nothing
    catch e
        e
    end
    @test empty_i isa ArgumentError
    @test occursin("termination", empty_i.msg)

    odd = try
        viterbi_decode(Int[0, 1, 0])
        nothing
    catch e
        e
    end
    @test odd isa ArgumentError
    @test occursin("even", odd.msg)

    # Non-finite soft values fail before the length check. Empty floats are finite.
    nan_err = try
        viterbi_decode(Float64[NaN])
        nothing
    catch e
        e
    end
    @test nan_err isa ArgumentError
    @test occursin("finite", nan_err.msg)
    @test !occursin("even", nan_err.msg)

    inf_err = try
        viterbi_decode(Float64[0.0, Inf])
        nothing
    catch e
        e
    end
    @test inf_err isa ArgumentError
    @test occursin("finite", inf_err.msg)

    empty_f = try
        viterbi_decode(Float64[])
        nothing
    catch e
        e
    end
    @test empty_f isa ArgumentError
    @test occursin("termination", empty_f.msg)

    # Exactly the tail: 4 coded bits, 0 information bits. Legal, including LLRs in [0, 1].
    @test viterbi_decode(Int[0, 0, 0, 0]) == Int[]
    @test viterbi_decode(Float64[0.0, 1.0, 0.0, 1.0]) == Int[]
end
