using Test
using Random
using LinearAlgebra
using CommToolkit

include("reference.jl")
include("contract.jl")

# Headless GR, in case the Plots extension is exercised.
ENV["GKSwstype"] = get(ENV, "GKSwstype", "nul")

@testset "Package hygiene" begin
    text = read(joinpath(pkgdir(CommToolkit), "Project.toml"), String)
    @test occursin("b0120378-1646-4a57-aae9-1dec7113c692", text)
    @test pkgversion(CommToolkit) == v"1.2.2"
    @test !occursin("a1b2c3d4-e5f6-7890-abcd-ef1234567890", text)
    deps = match(r"(?s)\[deps\](.*?)(?:\n\[|\z)", text).captures[1]
    @test !occursin("Plots", deps)
    @test occursin("[weakdeps]", text)
    @test occursin("CommToolkitPlotsExt", text)
    license = read(joinpath(pkgdir(CommToolkit), "LICENSE"), String)
    @test occursin("MIT License", license)
    @test occursin("刘恒宇", license)
    @test occursin("longee142857", license)
    plots_uuid = Base.UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80")
    @test all(id -> id.uuid != plots_uuid, keys(Base.loaded_modules))
end

@testset "Reference formulas" begin
    # M = 4 Gray QAM is exactly Gray QPSK / BPSK.
    for db in (0.0, 5.0, 10.0)
        γ = 10^(db / 10)
        @test ber_square_mqam(4, γ) ≈ ber_bpsk_awgn(γ) rtol = 1e-12
    end
    # Zero observed errors: closed upper endpoint.
    α = 0.05
    n = 1000
    _, hi = clopper_pearson(0, n, α)
    @test hi ≈ 1 - (α / 2)^(1 / n) rtol = 1e-12
    lo, _ = clopper_pearson(0, n, α)
    @test lo == 0.0
    @test clopper_pearson(n, n, α)[2] == 1.0
end

@testset "Modulation maps" begin
    @test modulate([0, 0, 1, 1], BPSK()) == [1, 1, -1, -1]
    @test demodulate([1.0, -1.0, 0.5, -0.5], BPSK()) == [0, 1, 0, 1]
    @test demodulate(ComplexF64[1, -0.2], BPSK()) == [0, 1]
    @test demodulate_soft([1.0, -1.0], BPSK(), 0.5) ≈ [4.0, -4.0]
    @test bits_per_symbol(BPSK()) == 1
    @test bits_per_symbol(QPSK()) == 2
    @test bits_per_symbol(QAM(64)) == 6
    @test snr_db_to_linear(10) ≈ 10.0
    @test snr_db_to_linear(0) ≈ 1.0

    @test_throws ArgumentError QAM(8)
    @test_throws ArgumentError QAM(2)
    @test_throws ArgumentError modulate([0, 1, 2], BPSK())
    @test_throws ArgumentError modulate([0, 1, 0], QPSK())

    for scheme in (BPSK(), QPSK(), QAM(16), QAM(64))
        Random.seed!(RNG_SEED)
        bps = bits_per_symbol(scheme)
        tx = rand(0:1, 2000 * bps)
        @test demodulate(modulate(tx, scheme), scheme) == tx
    end

    # Unit average energy on a complete constellation.
    for scheme in (QPSK(), QAM(16), QAM(64), QAM(256))
        pts = constellation(scheme)
        @test sum(abs2, pts) / length(pts) ≈ 1.0 rtol = 1e-12
    end
    @test sum(abs2, constellation(BPSK())) / 2 ≈ 1.0

    # Gray labelling: every minimum-distance pair differs by one bit.
    for M in (4, 16, 64)
        pts = constellation(QAM(M))
        k = Int(log2(M))
        dmin2 = 6 / (M - 1)
        n_neighbors = 0
        for i in 1:M, j in (i + 1):M
            abs2(pts[i] - pts[j]) ≈ dmin2 || continue
            n_neighbors += 1
            @test count_ones(xor(i - 1, j - 1)) == 1
        end
        @test n_neighbors > 0
        @test k == bits_per_symbol(QAM(M))
    end
end

@testset "Metrics" begin
    @test ber([0, 0, 0], [0, 0, 0]) == 0.0
    @test ber([0, 0, 0], [1, 1, 1]) == 1.0
    @test ber([0, 1, 0], [1, 1, 0]) == 1 / 3
    @test_throws ArgumentError ber([0, 1], [0])
end

function _mc_awgn(scheme, snr_db; target=250, max_bits=2_000_000, batch=100_000)
    bps = bits_per_symbol(scheme)
    k = 0
    n = 0
    while k < target && n < max_bits
        nb = min(batch, max_bits - n)
        nb -= mod(nb, bps)
        nb > 0 || break
        bits = rand(0:1, nb)
        rx = awgn(modulate(bits, scheme), snr_db, bps)
        k += count(!=(0), bits .- demodulate(rx, scheme))
        n += nb
    end
    return k, n
end

function _mc_rayleigh(snr_db; target=400, max_bits=1_000_000, batch=100_000)
    k = 0
    n = 0
    while k < target && n < max_bits
        nb = min(batch, max_bits - n)
        bits = rand(0:1, nb)
        rx = rayleigh_channel(modulate(bits, BPSK()), snr_db, 1, true)
        k += count(!=(0), bits .- demodulate(rx, BPSK()))
        n += nb
    end
    return k, n
end

@testset "AWGN BER vs closed form" begin
    Random.seed!(RNG_SEED)
    for snr_db in (0.0, 4.0, 6.0, 8.0)
        γ = 10^(snr_db / 10)
        k, n = _mc_awgn(BPSK(), snr_db)
        assert_theory_inside(k, n, ber_bpsk_awgn(γ))
    end
    Random.seed!(RNG_SEED)
    for snr_db in (0.0, 5.0, 8.0)
        γ = 10^(snr_db / 10)
        k, n = _mc_awgn(QPSK(), snr_db)
        assert_theory_inside(k, n, ber_qpsk_awgn(γ))
    end
    Random.seed!(RNG_SEED)
    for snr_db in (2.0, 6.0, 10.0, 12.0)
        γ = 10^(snr_db / 10)
        k, n = _mc_awgn(QAM(16), snr_db)
        assert_theory_inside(k, n, ber_square_mqam(16, γ))
    end
    Random.seed!(RNG_SEED)
    for snr_db in (6.0, 10.0, 14.0)
        γ = 10^(snr_db / 10)
        k, n = _mc_awgn(QAM(64), snr_db)
        assert_theory_inside(k, n, ber_square_mqam(64, γ))
    end
end

@testset "Rayleigh BPSK BER vs closed form" begin
    Random.seed!(RNG_SEED)
    for snr_db in (0.0, 8.0, 15.0, 20.0)
        γ = 10^(snr_db / 10)
        k, n = _mc_rayleigh(snr_db)
        assert_theory_inside(k, n, ber_bpsk_rayleigh(γ); min_errors=300)
    end
end

@testset "Convolutional code" begin
    Random.seed!(RNG_SEED)
    tx = rand(0:1, 200)
    coded = conv_encode(tx)
    @test length(coded) == 2 * (length(tx) + 2)
    @test viterbi_decode(coded) == tx

    coded1 = copy(coded)
    coded1[3] ⊻= 1
    @test viterbi_decode(coded1) == tx
    coded2 = copy(coded)
    coded2[10] ⊻= 1
    coded2[11] ⊻= 1
    @test viterbi_decode(coded2) == tx

    ad, cd = distance_spectrum(16)
    # Lin and Costello, 2nd ed., transfer function of the (2,1,2) code;
    # same coefficients as the classical (7, 5) K = 3 spectrum.
    @test ad[5:12] == [1, 2, 4, 8, 16, 32, 64, 128]
    @test cd[5:12] == [1, 4, 12, 32, 80, 192, 448, 1024]
    @test all(iszero, cd[1:4])

    # Hard-decision Viterbi on AWGN. Info-bit Eb/N0 uses Es = 1 on coded BPSK
    # and R = 1/2, i.e. awgn(..., energy_per_bit = 2). The BSC crossover is
    # Q(√(2 R Eb/N0)) = Q(√(Eb/N0)).
    for snr_db in (5.0, 6.0)
        Random.seed!(RNG_SEED)
        γb = 10^(snr_db / 10)
        p_bsc = qfunc(sqrt(γb))
        bound = conv_union_bound(p_bsc; dmax=24)
        uncoded = ber_bpsk_awgn(γb)
        target = 250
        k = 0
        n = 0
        frame = 2000
        while k < target && n < 800_000
            bits = rand(0:1, frame)
            rx = awgn(modulate(conv_encode(bits), BPSK()), snr_db; energy_per_bit=2.0)
            dec = viterbi_decode(demodulate(rx, BPSK()))
            k += count(!=(0), dec .- bits)
            n += frame
        end
        lo, hi = clopper_pearson(k, n, BER_ALPHA)
        @test k >= target
        # Significant coding gain against uncoded coherent BPSK at the same Eb/N0.
        @test hi < uncoded
        # Not worse than the hard-decision union bound (Proakis §8.2).
        @test hi <= bound
        @test lo < bound
    end
end

@testset "Soft demodulation" begin
    # Closed form, Proakis and Salehi, 5th ed., §4.2: L = 2 Re(y) / σ², σ² = N0/2.
    for approx in (:exact, :maxlog)
        @test demodulate_soft([1.0, -0.5, 0.0], BPSK(), 0.5; approx=approx) ≈ [4.0, -2.0, 0.0]
        @test demodulate_soft(ComplexF64[1.0 - 3im, -0.25 + 2im], BPSK(), 0.25; approx=approx) ≈
              [8.0, -2.0]
    end
    @test demodulate_soft(Float64[], BPSK(), 1.0) == Float64[]

    # Sign follows the bit map 0 ↦ +1, 1 ↦ −1, for every finite positive noise variance.
    yb = [-1.5, -0.2, 0.2, 1.5]
    Lb = demodulate_soft(yb, BPSK(), 0.8)
    @test all(Lb[yb .< 0] .< 0)
    @test all(Lb[yb .> 0] .> 0)

    # |L| scales as 1/σ². Same convention as awgn: σ² = Eb / (2 Eb/N0) with Eb = 1.
    y0 = [0.8]
    vars = (1.0, 0.5, 0.2)
    mags = [only(demodulate_soft(y0, BPSK(), σ2)) for σ2 in vars]
    @test mags ≈ [1.6, 3.2, 8.0]
    @test issorted(abs.(mags))

    # Gray QPSK: exact = max-log = (√2 / σ²) · (Re y, Im y).
    σ2 = 0.25
    yq = [0.3 - 0.4im, -0.2 + 0.1im]
    closed = Float64[]
    for y in yq
        append!(closed, (sqrt(2) / σ2) .* (real(y), imag(y)))
    end
    Lex = demodulate_soft(yq, QPSK(), σ2; approx=:exact)
    Lmx = demodulate_soft(yq, QPSK(), σ2; approx=:maxlog)
    @test Lex ≈ closed rtol = 1e-10
    @test Lmx ≈ closed rtol = 1e-10

    # Frozen 16-QAM point. Independent enumeration of the Gray lattice at σ² = 1/2
    # (N0 = 1): see the max-log / log-sum-exp formulas in the demodulate_soft docstring.
    y16 = ComplexF64[0.2 + 0.5im]
    @test demodulate_soft(y16, QAM(16), 0.5; approx=:exact) ≈
          [-0.41028054307236383, -0.7376431735811475, -1.0312550627687205, -0.4539223388071142] atol = 1e-9
    @test demodulate_soft(y16, QAM(16), 0.5; approx=:maxlog) ≈
          [-0.25298221281347033, -0.5470177871865296, -0.632455532033676, -0.16754446796632413] atol = 1e-9

    # Square Gray QAM factors across I and Q. The first k/2 LLRs ignore Im(y).
    for (M, σq) in ((16, 0.3), (64, 0.2))
        k = bits_per_symbol(QAM(M))
        half = k ÷ 2
        for approx in (:exact, :maxlog)
            L1 = demodulate_soft([0.2 + 0.5im], QAM(M), σq; approx=approx)
            L2 = demodulate_soft([0.2 - 1.1im], QAM(M), σq; approx=approx)
            L3 = demodulate_soft([-0.7 + 0.5im], QAM(M), σq; approx=approx)
            @test L1[1:half] ≈ L2[1:half] rtol = 1e-10
            @test L1[(half + 1):k] ≈ L3[(half + 1):k] rtol = 1e-10
            @test !(L1[(half + 1):k] ≈ L2[(half + 1):k])
        end
    end

    # At a constellation point the nearest neighbour dominates: hard slice matches demodulate.
    for scheme in (BPSK(), QPSK(), QAM(4), QAM(16), QAM(64), QAM(256))
        pts = constellation(scheme)
        for approx in (:exact, :maxlog)
            L = demodulate_soft(pts, scheme, 1e-6; approx=approx)
            @test all(isfinite, L)
            @test Int.(L .< 0) == demodulate(pts, scheme)
        end
    end

    # Max-log hard decisions are minimum-distance decisions whenever the ML point is unique.
    Random.seed!(RNG_SEED)
    yml = randn(200) .+ im .* randn(200)
    for scheme in (QPSK(), QAM(16), QAM(64))
        pts = constellation(scheme)
        L = demodulate_soft(yml, scheme, 0.4; approx=:maxlog)
        hard = demodulate(yml, scheme)
        @test Int.(L .< 0) == hard
    end

    # 16-QAM exact |LLR| grows as σ² falls, for a symbol sitting on a constellation point.
    p0 = constellation(QAM(16))[1:1]
    prev = zeros(4)
    for σ2 in (0.2, 0.1, 0.05, 0.02)
        mag = abs.(demodulate_soft(p0, QAM(16), σ2; approx=:exact))
        @test all(mag .> prev)
        prev = mag
    end
    # Max-log magnitude is exactly proportional to 1/N0.
    Lm1 = demodulate_soft(y16, QAM(16), 0.4; approx=:maxlog)
    Lm2 = demodulate_soft(y16, QAM(16), 0.1; approx=:maxlog)
    @test Lm2 ≈ 4 .* Lm1 rtol = 1e-10

    # Conditional law of the BPSK LLR. Bit 0: y ∼ N(1, σ²), L = 2y/σ² ∼ N(2/σ², 4/σ²).
    Random.seed!(RNG_SEED)
    σb = 0.5
    nb = 20_000
    ybit = 1 .+ sqrt(σb) .* randn(nb)
    Lbmc = demodulate_soft(ybit, BPSK(), σb)
    μb = 2 / σb
    σLb = 2 / sqrt(σb)
    kmed = count(<(μb), Lbmc)
    ktail = count(<(μb - σLb), Lbmc)
    lo_m, hi_m = clopper_pearson(kmed, nb, BER_ALPHA)
    lo_t, hi_t = clopper_pearson(ktail, nb, BER_ALPHA)
    @test lo_m <= 0.5 <= hi_m
    @test lo_t <= qfunc(1.0) <= hi_t

    # Gray QPSK in-phase bit, transmitted 0: L ∼ N(1/σ², 2/σ²).
    Random.seed!(RNG_SEED)
    nq = 20_000
    a = 1 / sqrt(2)
    yI = a .+ sqrt(σb) .* randn(nq) .+ im .* (a .+ sqrt(σb) .* randn(nq))
    LI = demodulate_soft(yI, QPSK(), σb)[1:2:end]
    μI = 1 / σb
    σLI = sqrt(2 / σb)
    lo_q, hi_q = clopper_pearson(count(<(μI), LI), nq, BER_ALPHA)
    lo_qs, hi_qs = clopper_pearson(count(<(μI - σLI), LI), nq, BER_ALPHA)
    @test lo_q <= 0.5 <= hi_q
    @test lo_qs <= qfunc(1.0) <= hi_qs

    # awgn's σ² = Eb / (2 Eb/N0). Positional BPSK uses Eb = 1. The LLR of a
    # noise-free +1 is 2/σ²; symbols that actually passed through awgn must
    # follow the same conditional law, which locks the two calibrations together.
    snr_db = 0.0
    σ_awgn = 1 / (2 * 10^(snr_db / 10))
    @test only(demodulate_soft([1.0], BPSK(), σ_awgn)) ≈ 2 / σ_awgn
    Random.seed!(RNG_SEED)
    ncal = 20_000
    rxcal = awgn(modulate(zeros(Int, ncal), BPSK()), snr_db, 1)
    Lcal = demodulate_soft(rxcal, BPSK(), σ_awgn)
    μcal = 2 / σ_awgn
    σLcal = 2 / sqrt(σ_awgn)
    lo_c, hi_c = clopper_pearson(count(<(μcal), Lcal), ncal, BER_ALPHA)
    lo_cs, hi_cs = clopper_pearson(count(<(μcal - σLcal), Lcal), ncal, BER_ALPHA)
    @test lo_c <= 0.5 <= hi_c
    @test lo_cs <= qfunc(1.0) <= hi_cs

    @test_throws ArgumentError demodulate_soft([1.0], BPSK(), 0.0)
    @test_throws ArgumentError demodulate_soft([1.0], BPSK(), -0.1)
    @test_throws ArgumentError demodulate_soft([1.0], BPSK(), Inf)
    @test_throws ArgumentError demodulate_soft([1.0 + 0im], QPSK(), 0.2; approx=:approx)
    @test_throws ArgumentError demodulate_soft([0.1, -0.2], QPSK(), 0.2)
    @test_throws ArgumentError demodulate_soft([0.1], QAM(16), 0.2)
end

@testset "Soft-input decoding" begin
    Random.seed!(RNG_SEED)
    # L = 1 − 2ĉ reproduces Hamming Viterbi, including the tie rule.
    for n in (1, 7, 40)
        bits = rand(0:1, n)
        coded = conv_encode(bits)
        for flips in (Int[], [1], [2, 5], [1, 4, 8, 9])
            rx = copy(coded)
            for i in flips
                i <= length(rx) || continue
                rx[i] ⊻= 1
            end
            llr = Float64.(1 .- 2 .* rx)
            @test viterbi_decode(llr) == viterbi_decode(rx)
            @test viterbi_decode(1000 .* llr) == viterbi_decode(rx)
        end
    end

    bits = rand(0:1, 32)
    coded = conv_encode(bits)
    strong = 30 .* Float64.(1 .- 2 .* coded)
    @test viterbi_decode(strong) == bits
    Le = bcjr_decode(strong)
    @test length(Le) == length(bits)
    @test Int.(Le .< 0) == bits
    @test all(abs.(Le) .> 20)

    # Tie on every branch: input 0 is written first, so the survivor is the zero word.
    @test viterbi_decode(zeros(2 * (15 + 2))) == zeros(Int, 15)

    # Free-distance event of the (7, 5) code: input 1,0,0 produces coded 11 10 11,
    # and the two tail bits stay 0 once the encoder has returned to state 0.
    event = conv_encode([1, 0, 0])
    @test event == [1, 1, 1, 0, 1, 1, 0, 0, 0, 0]
    @test viterbi_decode(Float64.(1 .- 2 .* event)) == [1, 0, 0]

    # Uninformative channel, uniform prior: every information bit has APP LLR 0.
    @test bcjr_decode(zeros(2 * (6 + 2))) ≈ zeros(6) atol = 1e-9

    # Brute-force MAP on the terminated code. n = 8 information bits, 256 codewords.
    function brute_force_app(llr::Vector{Float64}, apriori::Vector{Float64})
        n = length(apriori)
        acc0 = fill(-Inf, n)
        acc1 = fill(-Inf, n)
        for mask in 0:(2^n - 1)
            info = [(mask >> (k - 1)) & 1 for k in 1:n]
            coded = conv_encode(info)
            logp = 0.0
            @inbounds for i in eachindex(coded)
                logp += 0.5 * (1 - 2 * coded[i]) * llr[i]
            end
            @inbounds for k in 1:n
                logp += 0.5 * (1 - 2 * info[k]) * apriori[k]
            end
            @inbounds for k in 1:n
                if info[k] == 0
                    acc0[k] = acc0[k] == -Inf ? logp : max(acc0[k], logp) + log1p(exp(-abs(acc0[k] - logp)))
                else
                    acc1[k] = acc1[k] == -Inf ? logp : max(acc1[k], logp) + log1p(exp(-abs(acc1[k] - logp)))
                end
            end
        end
        return acc0 .- acc1
    end

    for trial in 1:4
        info_n = 8
        La = trial == 1 ? zeros(info_n) : randn(info_n)
        ch = 3 .* randn(2 * (info_n + 2))
        app = brute_force_app(ch, La)
        @test bcjr_decode(ch; apriori=La) ≈ app .- La atol = 1e-6 rtol = 1e-6
    end
    ch0 = randn(2 * (4 + 2))
    @test bcjr_decode(ch0) ≈ bcjr_decode(ch0; apriori=zeros(4)) atol = 1e-12

    @test viterbi_decode(conv_encode(Int[])) == Int[]
    @test bcjr_decode(zeros(4)) == Float64[]

    @test_throws ArgumentError viterbi_decode(Float64[0.1, 0.2, 0.3])
    @test_throws ArgumentError viterbi_decode(Float64[0.0, NaN])
    @test_throws ArgumentError viterbi_decode(Float64[0.0, Inf])
    @test_throws ArgumentError bcjr_decode(Float64[0.1])
    @test_throws ArgumentError bcjr_decode(zeros(8); apriori=[0.0])
    @test_throws ArgumentError bcjr_decode(zeros(8); apriori=[0.0, NaN])

    # Soft union bound: free-distance term, monotone in Eb/N0, finite at the SNRs used below.
    b4 = conv_soft_union_bound(10^(4 / 10))
    b6 = conv_soft_union_bound(10^(6 / 10))
    @test b4 > qfunc(sqrt(5 * 10^(4 / 10)))
    @test b6 > qfunc(sqrt(5 * 10^(6 / 10)))
    @test b6 < b4 < 1
    @test b6 < ber_bpsk_awgn(10^(6 / 10))

    # Same Eb/N0, same frames. Soft-decision PEP is Q(√(d γb)); hard decisions quantise first.
    # Frame length 500 keeps the terminated-code BCJR cheap; the tail is 2 bits.
    function coded_counts(snr_db::Float64; max_info::Int, frame::Int)
        σ2 = 2.0 / (2 * 10^(snr_db / 10))
        k_hard = 0
        k_soft = 0
        k_map = 0
        n = 0
        while n < max_info
            bits = rand(0:1, frame)
            rx = awgn(modulate(conv_encode(bits), BPSK()), snr_db; energy_per_bit=2.0)
            hard_bits = demodulate(rx, BPSK())
            llr = demodulate_soft(rx, BPSK(), σ2)
            k_hard += count(!=(0), viterbi_decode(hard_bits) .- bits)
            k_soft += count(!=(0), viterbi_decode(llr) .- bits)
            k_map += count(!=(0), Int.(bcjr_decode(llr) .< 0) .- bits)
            n += frame
        end
        return k_hard, k_soft, k_map, n
    end

    # 3 dB and 4 dB. Clopper–Pearson intervals (α = 10⁻⁴) of soft Viterbi and
    # BCJR sit strictly under the hard-decision interval and under uncoded BPSK.
    # The soft point estimate is compared with the truncated union bound: below
    # it, and within a factor of ten. The interval itself is not required to lie
    # under the bound — at 4 dB the bound is only about 1.5× the estimate, so a
    # different `randn` stream can push the upper endpoint over it.
    # Hard-decision coding gain at 5 dB and 6 dB is covered by the testset above.
    for snr_db in (3.0, 4.0)
        Random.seed!(RNG_SEED)
        γb = 10^(snr_db / 10)
        kh, ks, km, nbits = coded_counts(snr_db; max_info=200_000, frame=500)
        lo_h, hi_h = clopper_pearson(kh, nbits, BER_ALPHA)
        _, hi_s = clopper_pearson(ks, nbits, BER_ALPHA)
        _, hi_m = clopper_pearson(km, nbits, BER_ALPHA)
        bound = conv_soft_union_bound(γb)
        hard_bound = conv_union_bound(qfunc(sqrt(γb)))
        uncoded = ber_bpsk_awgn(γb)
        p_soft = ks / nbits
        p_map = km / nbits
        @test kh >= 100
        @test ks >= 20
        @test km >= 20
        @test hi_s < lo_h
        @test hi_m < lo_h
        @test hi_s < uncoded
        @test hi_m < uncoded
        @test p_soft < bound
        @test p_map < bound
        @test p_soft > bound / 10
        @test p_map > bound / 10
        @test hi_h <= hard_bound
    end
end

@testset "Equalizers" begin
    h = [1.0, 0.5]
    syms = [1.0, -1.0, 1.0, -1.0]
    L = length(h)
    rx = zeros(length(syms) + L - 1)
    for i in eachindex(rx)
        for j in 1:L
            src = i - j + 1
            if 1 <= src <= length(syms)
                rx[i] += syms[src] * h[j]
            end
        end
    end
    @test zf_equalize(rx, h) ≈ syms atol = 1e-10
    @test mmse_equalize(rx, h, 0.0) ≈ syms atol = 1e-10

    hc = [1.0 + 0.2im, 0.4 - 0.1im]
    qsyms = modulate([0, 1, 1, 0, 1, 1, 0, 0, 1, 0], QPSK())
    rxc = zeros(ComplexF64, length(qsyms) + length(hc) - 1)
    for i in eachindex(rxc)
        for j in eachindex(hc)
            src = i - j + 1
            if 1 <= src <= length(qsyms)
                rxc[i] += qsyms[src] * hc[j]
            end
        end
    end
    @test zf_equalize(rxc, hc) ≈ qsyms atol = 1e-8
    @test mmse_equalize(rxc, hc, 0.0) ≈ qsyms atol = 1e-8
end

@testset "Pulse shaping" begin
    beta, span, sps = 0.5, 6, 4
    rc = raised_cosine(beta, span, sps)
    @test length(rc) == span * sps + 1
    @test sum(rc) ≈ 1 atol = 1e-10
    rrc = rrc_filter(beta, span, sps)
    @test length(rrc) == span * sps + 1
    @test sum(abs2, rrc) ≈ 1 atol = 1e-10
    @test sum(raised_cosine(0.0, 8, 4)) ≈ 1 atol = 1e-10

    # Two truncated RRC filters in cascade are nearly Nyquist at symbol spacing.
    β, span2, sps2 = 0.25, 16, 8
    g = rrc_filter(β, span2, sps2)
    ntap = length(g)
    casc = zeros(2ntap - 1)
    for i in 1:ntap, j in 1:ntap
        casc[i + j - 1] += g[i] * g[j]
    end
    peak = argmax(casc)
    @test casc[peak] ≈ 1 atol = 1e-8
    for lag in 1:4
        @test abs(casc[peak + lag * sps2]) < 2e-3
    end

    syms = [1.0, -1.0, 1.0]
    @test length(pulse_shape(syms, rc, sps)) == length(syms) * sps
end

@testset "OFDM" begin
    Random.seed!(RNG_SEED)
    params = OFDMParams(N_fft=32, N_cp=8, mod_order=4, n_guard=2, pilot_spacing=4)
    @test CommToolkit.modulation(params) isa QPSK
    bits = rand(0:1, params.N_data * params.bits_per_symbol * 6)
    tx, pilots = ofdm_tx(params, bits)

    multipath = TDLChannel(delays_samples=[0, 2, 5], gain_dB=[0.0, -3.0, -6.0])
    @test maximum(multipath.delays) <= params.N_cp
    rx, H = apply_channel(tx, multipath, Inf, params)
    got = ofdm_rx_known(rx, params, H, pilots)
    @test ber(bits, got[1:length(bits)]) == 0.0

    p16 = OFDMParams(N_fft=32, N_cp=8, mod_order=16, n_guard=2, pilot_spacing=4)
    b16 = rand(0:1, p16.N_data * p16.bits_per_symbol * 4)
    tx16, pil16 = ofdm_tx(p16, b16)
    rx16, H16 = apply_channel(tx16, multipath, Inf, p16)
    @test ber(b16, ofdm_rx_known(rx16, p16, H16, pil16)[1:length(b16)]) == 0.0

    flat = TDLChannel(delays_samples=[0], gain_dB=[0.0])
    rxf, Hf = apply_channel(tx, flat, Inf, params)
    @test ber(bits, ofdm_rx_known(rxf, params, Hf, pilots)[1:length(bits)]) == 0.0
    @test ber(bits, ofdm_rx_ls(rxf, params, pilots)[1:length(bits)]) == 0.0
    @test ber(bits, ofdm_rx_mmse(rxf, params, pilots, 40.0)[1:length(bits)]) == 0.0

    # Timing: 1-based index of the first sample of the burst.
    for delay in (0, 1, 5, 7)
        delayed = vcat(zeros(ComplexF64, delay), tx)
        @test sto_estimate(delayed, params) == delay + 1
        @test sto_estimate(delayed, params; search_range=3 * (params.N_fft + params.N_cp)) == delay + 1
    end

    N = params.N_fft
    for cfo in (0.0, 0.17, -0.23)
        rxc = [tx[i] * exp(im * 2π * cfo * (i - 1) / N) for i in eachindex(tx)]
        @test cfo_estimate(rxc, params) ≈ cfo atol = 1e-8
        corrected, sto, cfo_hat = ofdm_sync(rxc, params)
        @test sto == 1
        @test cfo_hat ≈ cfo atol = 1e-8
        Hflat = ones(ComplexF64, params.N_fft)
        @test ber(bits, ofdm_rx_known(corrected, params, Hflat, pilots)[1:length(bits)]) == 0.0
    end

    # Known-CSI BER falls as the time-domain SNR rises. Intervals must separate.
    function ofdm_errors(snr_db, n_sym)
        b = rand(0:1, params.N_data * params.bits_per_symbol * n_sym)
        sig, pil = ofdm_tx(params, b)
        y, Hf = apply_channel(sig, multipath, snr_db, params)
        dec = ofdm_rx_known(y, params, Hf, pil)
        k = count(!=(0), b .- dec[1:length(b)])
        return k, length(b)
    end
    Random.seed!(RNG_SEED)
    k_lo, n_lo = ofdm_errors(4.0, 40)
    k_hi, n_hi = ofdm_errors(18.0, 40)
    lo_lo, lo_hi = clopper_pearson(k_lo, n_lo, BER_ALPHA)
    _, hi_hi = clopper_pearson(k_hi, n_hi, BER_ALPHA)
    @test k_lo >= 20
    @test hi_hi < lo_lo
end

@testset "MIMO" begin
    Random.seed!(RNG_SEED)
    Hrand = mimo_channel(2, 2)
    @test size(Hrand) == (2, 2)

    H = ComplexF64[1.2+0.1im 0.3-0.2im; 0.15+0.05im -1.0+0.4im]
    @test cond(H) < 10
    bpsk = ComplexF64[1, -1]
    qpsk = constellation(QPSK())
    x = ComplexF64[1, -1]
    y = apply_mimo(x, H, 0.0)
    @test y ≈ H * x
    @test ml_detect(y, H, bpsk) ≈ x
    @test zf_detect(y, H) ≈ x atol = 1e-10
    @test mmse_detect(y, H, 0.0) ≈ x atol = 1e-10
    @test sic_detect(y, H, bpsk) ≈ x
    xq = qpsk[[1, 4]]
    yq = apply_mimo(xq, H, 0.0)
    @test ml_detect(yq, H, qpsk) ≈ xq
    @test sic_detect(yq, H, qpsk) ≈ xq

    Ht = mimo_channel(4, 2) # n_tx = 4, n_rx = 2
    P = zf_precoder(Ht)
    @test size(P) == (4, 2)
    @test Ht * P ≈ I(2) atol = 1e-8
    Hs = ComplexF64.(H)
    @test Hs * zf_precoder(Hs) ≈ I(2) atol = 1e-8

    V, S = svd_precoder(H)
    F = svd(H)
    @test V ≈ F.V
    @test S ≈ F.S
    @test H * V ≈ F.U * Diagonal(S) atol = 1e-8
    # The pre-0.2 return value was Vᴴ, which is not a precoder.
    @test size(V) == (size(H, 2), length(S))

    code = alamouti_encode(0.3 - 0.4im, 0.8 + 0.1im)
    @test code[1, 2] ≈ -conj(code[2, 1])
    @test code[2, 2] ≈ conj(code[1, 1])
    h1, h2 = 0.5 + 0.3im, -0.2 + 0.8im
    r1 = h1 * code[1, 1] + h2 * code[2, 1]
    r2 = h1 * code[1, 2] + h2 * code[2, 2]
    x1, x2 = alamouti_decode_2x1(r1, r2, h1, h2)
    @test x1 ≈ code[1, 1] atol = 1e-10
    @test x2 ≈ code[2, 1] atol = 1e-10

    Ha = [0.5+0.3im -0.2+0.8im; 0.1-0.4im 0.7+0.2im]
    r1v = Ha * code[:, 1]
    r2v = Ha * code[:, 2]
    y1, y2 = alamouti_decode_2x2(r1v, r2v, Ha)
    @test y1 ≈ code[1, 1] atol = 1e-10
    @test y2 ≈ code[2, 1] atol = 1e-10

    # ZF hard decisions: high-SNR interval sits strictly below the low-SNR interval.
    function mimo_zf_errors(noise_var, n_vec)
        k = 0
        for _ in 1:n_vec
            bits = rand(0:1, 2)
            s = ComplexF64.(modulate(bits, BPSK()))
            z = zf_detect(apply_mimo(s, H, noise_var), H)
            k += count(!=(0), bits .- demodulate(real.(z), BPSK()))
        end
        return k, 2 * n_vec
    end
    Random.seed!(RNG_SEED)
    k0, n0 = mimo_zf_errors(1.0, 4000)       # complex N0 = 1, Es = 1
    k1, n1 = mimo_zf_errors(1e-3, 4000)      # 30 dB above that
    lo0, hi0 = clopper_pearson(k0, n0, BER_ALPHA)
    _, hi1 = clopper_pearson(k1, n1, BER_ALPHA)
    @test k0 >= 20
    @test hi1 < lo0
    @test hi0 > hi1
end

include("interleave.jl")
include("ldpc.jl")
include("polar.jl")
include("link.jl")
include("xapp.jl")
include("sova.jl")
include("turbo.jl")
include("algebraic.jl")
include("pulse_fir.jl")
include("lms_dfe.jl")
include("psk_rice_mrc.jl")
include("evm.jl")
include("scrambling.jl")
include("waterfill.jl")
include("kronecker_mimo.jl")
include("sync.jl")
include("jakes_nakagami.jl")
include("cdl.jl")
include("bicm.jl")
include("exit.jl")

@testset "Plots extension" begin
    @test_throws ErrorException plot_constellation(ComplexF64[1.0, -1.0], BPSK())
    @test_throws ErrorException plot_exit_chart(EXITCurve([0.0, 1.0], [0.0, 1.0]))
    @eval using Plots
    Random.seed!(RNG_SEED)
    rx = awgn(modulate(rand(0:1, 400), QPSK()), 12.0, 2)
    plt = plot_constellation(rx, QPSK(); title="QPSK")
    @test plt isa Plots.Plot
    g = rrc_filter(0.35, 6, 4)
    shaped = real.(pulse_shape(modulate(rand(0:1, 80), BPSK()), g, 4))
    eye = plot_eye_diagram(shaped, 4; n_segments=10)
    @test eye isa Plots.Plot
    pref = plot_pulse_response(raised_cosine(0.35, 6, 4))
    @test pref isa Plots.Plot
    demap_curve = EXITCurve(collect(0.0:0.5:1.0), [0.5, 0.6, 0.7])
    decoder_curve = EXITCurve(collect(0.0:0.5:1.0), [0.4, 0.7, 0.9])
    @test plot_exit_chart(demap_curve) isa Plots.Plot
    @test plot_exit_chart(demap_curve, decoder_curve) isa Plots.Plot
end
