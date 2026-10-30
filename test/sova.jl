@testset "SOVA" begin
    # Noiseless terminated word: ML = MAP, and the LLR signs agree.
    bits = [1, 0, 1, 1, 0, 0, 1, 0, 1, 1]
    coded = conv_encode(bits)
    strong = 40.0 .* (1.0 .- 2.0 .* coded)
    Ls = sova_decode(strong)
    Lb = bcjr_decode(strong)
    @test length(Ls) == length(bits)
    @test all(isfinite, Ls)
    @test all(Ls .* Lb .> 0)
    @test Int.(Ls .< 0) == bits
    @test Int.(Lb .< 0) == bits
    @test viterbi_decode(strong) == bits
    @test all(abs.(Ls) .> 5)

    # Same boundary, same hard decision, on a correctable free-distance event.
    event = conv_encode([1, 0, 0])
    @test event == [1, 1, 1, 0, 1, 1, 0, 0, 0, 0]
    Le = 25.0 .* (1.0 .- 2.0 .* event)
    @test Int.(sova_decode(Le) .< 0) == [1, 0, 0]
    @test all(sova_decode(Le) .* bcjr_decode(Le) .> 0)
    @test viterbi_decode(Le) == [1, 0, 0]

    # One and two coded-bit flips stay inside dfree. Signs still match BCJR.
    for flips in ([3], [4, 7])
        rx = copy(coded)
        for i in flips
            rx[i] ⊻= 1
        end
        llr = 20.0 .* (1.0 .- 2.0 .* rx)
        @test Int.(sova_decode(llr) .< 0) == bits
        @test Int.(bcjr_decode(llr) .< 0) == bits
        @test viterbi_decode(llr) == bits
        @test all(sova_decode(llr) .* bcjr_decode(llr) .> 0)
    end

    # Uninformative channel: every information bit has LLR 0, as in BCJR.
    z = zeros(2 * (10 + 2))
    @test sova_decode(z) ≈ zeros(10) atol = 1e-8
    @test bcjr_decode(z) ≈ zeros(10) atol = 1e-8

    # Branch metric is linear in the channel LLRs, so a positive scale doubles the SOVA LLRs.
    # Negating the channel is not an odd map: the tail inputs are fixed at 0.
    llr = Float64[
        0.2, -0.4, 1.1, 0.3, -0.7, 0.5, 0.8, -0.2,
        0.4, 1.2, -0.6, 0.9, 0.1, -1.0, 0.3, 0.6,
    ]
    @test sova_decode(2.0 .* llr) ≈ 2.0 .* sova_decode(llr) atol = 1e-8 rtol = 1e-8

    # High-SNR perturbation of a real codeword: signs still agree.
    word = Float64.(conv_encode([1, 0, 1, 0, 1, 1, 0, 0, 1, 0]))
    tilted = 8.0 .* (1.0 .- 2.0 .* word) .+ 0.15 .* sin.(eachindex(word))
    @test all(sova_decode(tilted) .* bcjr_decode(tilted) .> 0)
    @test Int.(sova_decode(tilted) .< 0) == [1, 0, 1, 0, 1, 1, 0, 0, 1, 0]
    @test viterbi_decode(tilted) == [1, 0, 1, 0, 1, 1, 0, 0, 1, 0]

    # A priori of zero is the default. A mild prior does not flip a decisive channel,
    # and the extrinsic signs still agree with BCJR.
    ninfo = length(bits)
    @test sova_decode(strong; apriori=zeros(ninfo)) ≈ Ls atol = 1e-9
    La = fill(0.4, ninfo)
    @test all(sova_decode(strong; apriori=La) .* bcjr_decode(strong; apriori=La) .> 0)
    @test Int.((sova_decode(strong; apriori=La) .+ La) .< 0) == bits

    # Erasure from the rate-2/3 mask: a strong kept-bit channel still agrees.
    mask = puncture_mask(2//3)
    kept = puncture(coded, mask)
    mother = depuncture(30.0 .* (1.0 .- 2.0 .* Float64.(kept)), mask, length(coded); fill=0.0)
    @test Int.(sova_decode(mother) .< 0) == bits
    @test all(sova_decode(mother) .* bcjr_decode(mother) .> 0)

    @test sova_decode(zeros(4)) == Float64[]
    @test_throws ArgumentError sova_decode(Float64[0.1, 0.2, 0.3])
    @test_throws ArgumentError sova_decode(Float64[0.0, NaN])
    @test_throws ArgumentError sova_decode(Float64[0.0, Inf])
    @test_throws ArgumentError sova_decode(zeros(8); apriori=[0.0])
    @test_throws ArgumentError sova_decode(zeros(8); apriori=[0.0, NaN])
end

@testset "SOVA coded extrinsic" begin
    # Default keyword is the 1.2.1 vector. return_coded=true keeps that vector
    # as the first element and adds one extrinsic LLR per coded bit.
    bits = [1, 0, 1, 1, 0, 0, 1, 0, 1, 1]
    coded = conv_encode(bits)
    strong = 40.0 .* (1.0 .- 2.0 .* coded)
    info, Lc = sova_decode(strong; return_coded=true)
    @test sova_decode(strong; return_coded=false) == sova_decode(strong)
    @test info == sova_decode(strong)
    @test length(Lc) == length(coded)
    @test Int.(info .< 0) == bits
    @test viterbi_decode(strong) == bits
    # Equal-amplitude noiseless word: every competing gap that fires is the
    # free-distance event, ρ = 5·40. Coded extrinsic removes the channel LLR,
    # so |L_e(c)| = 160 and its sign matches the coded bit. Not a BCJR marginal.
    @test all(abs.(info) .≈ 200.0)
    @test all(abs.(Lc) .≈ 160.0)
    @test Int.(Lc .< 0) == coded
    @test Int.((Lc .+ strong) .< 0) == coded

    La = fill(0.4, length(bits))
    info_a, Lc_a = sova_decode(strong; apriori=La, return_coded=true)
    @test info_a ≈ sova_decode(strong; apriori=La) atol = 1e-9
    @test Int.((info_a .+ La) .< 0) == bits
    @test Int.((Lc_a .+ strong) .< 0) == coded

    # Fixed channel, deterministic traceback, including the coded vector.
    llr = Float64[
        0.2, -0.4, 1.1, 0.3, -0.7, 0.5, 0.8, -0.2,
        0.4, 1.2, -0.6, 0.9, 0.1, -1.0, 0.3, 0.6,
    ]
    @test sova_decode(llr) ≈ [-1.4, -1.4, -2.3, 1.6, -2.3, 1.6] atol = 1e-9
    _, coded_llr = sova_decode(llr; return_coded=true)
    @test coded_llr ≈ [
        -1.6, -1.0, 1.4, -1.7, -1.6, 0.9, 0.6, -1.2,
        1.2, 1.2, -1.7, 1.4, -1.7, -1.3, 1.3, 1.0,
    ] atol = 1e-9
    @test sova_decode(llr; return_coded=true)[1] == sova_decode(llr)
    # Branch metrics are linear in the channel LLRs, and coded extrinsic
    # subtracts the channel LLR, so the whole pair scales.
    info2, coded2 = sova_decode(2.0 .* llr; return_coded=true)
    @test info2 ≈ 2.0 .* sova_decode(llr) atol = 1e-8
    @test coded2 ≈ 2.0 .* coded_llr atol = 1e-8

    # Tie: an all-zero channel gives every path the same metric. The survivor
    # prefers the smaller input bit, then the smaller previous state, so the
    # traceback is the all-zero word. Disagreements have Δ = 0, and a zero
    # gap cancels the prior, so every extrinsic LLR is 0. `L .< 0` reads that
    # zero as bit 0, which matches this survivor. The choice is deterministic.
    z = zeros(2 * (10 + 2))
    zi, zc = sova_decode(z; return_coded=true)
    @test zi ≈ zeros(10) atol = 1e-8
    @test zc ≈ zeros(length(z)) atol = 1e-8
    @test zi == sova_decode(z)

    # Truncation: a one-bit terminated word. Five coded bits see a weight-5
    # competitor (ρ = 5A). The remaining coded bit never splits along the
    # survivor, so its reliability is the finite cap Σ|L|, not ±∞.
    one = [1]
    one_coded = conv_encode(one)
    A = 12.0
    Lone = A .* (1.0 .- 2.0 .* one_coded)
    info_one, Lc_one = sova_decode(Lone; return_coded=true)
    bound = sum(abs, Lone)
    ρc = abs.(Lc_one .+ Lone)
    @test length(one_coded) == 6
    @test info_one ≈ [-5A]
    @test abs(info_one[1]) < bound
    @test count(r -> isapprox(r, 5A; atol=1e-8), ρc) == 5
    @test count(r -> isapprox(r, bound; atol=1e-8), ρc) == 1
    @test ρc[4] ≈ bound atol = 1e-8
    @test all(isfinite, Lc_one)
    @test Int.((Lc_one .+ Lone) .< 0) == one_coded

    # Tail only: no information bit, four coded LLRs, still a tie at LLR 0.
    empty_info, tail_llr = sova_decode(zeros(4); return_coded=true)
    @test empty_info == Float64[]
    @test tail_llr ≈ zeros(4) atol = 1e-8
    @test_throws ArgumentError sova_decode(Float64[0.1, 0.2, 0.3]; return_coded=true)
end

@testset "SOVA versus BCJR" begin
    # Es = 1 on coded BPSK, Eb/N0 through energy_per_bit = 2 (nominal R = 1/2),
    # the same bookkeeping as the convolutional section of test/runtests.jl,
    # so σ² = 1/γb. Seed 0x0c0de1162. Teaching SOVA is not MAP: the same point
    # may favour bcjr_decode. Counts below were checked on Julia 1.10.12 and
    # 1.13.3; assertions are the stable inequalities, not the raw counts.
    function _sova_bcjr_counts(snr_db, frames, ninfo)
        γ = 10^(snr_db / 10)
        σ2 = 1 / γ
        es = eb = vit_dis = hard_dis = 0
        coded_app = coded_ext = 0
        mis = mib = 0.0
        Random.seed!(0x0c0de1162)
        for _ in 1:frames
            bits = rand(0:1, ninfo)
            coded = conv_encode(bits)
            y = awgn(modulate(coded, BPSK()), snr_db; energy_per_bit=2.0)
            L = demodulate_soft(y, BPSK(), σ2)
            Ls = sova_decode(L)
            Lb = bcjr_decode(L)
            info, Lc = sova_decode(L; return_coded=true)
            @test info == Ls
            hs = Int.(Ls .< 0)
            hb = Int.(Lb .< 0)
            es += count(!=(0), hs .- bits)
            eb += count(!=(0), hb .- bits)
            vit_dis += count(!=(0), hs .- viterbi_decode(L))
            hard_dis += count(!=(0), hs .- hb)
            mis += exit_mi_soft(Ls, bits)
            mib += exit_mi_soft(Lb, bits)
            coded_app += count(!=(0), Int.((Lc .+ L) .< 0) .- coded)
            coded_ext += count(!=(0), Int.(Lc .< 0) .- coded)
        end
        return (; es, eb, vit_dis, hard_dis, mis=mis / frames, mib=mib / frames,
                coded_app, coded_ext, nbits=frames * ninfo)
    end

    # 0 dB, 48 frames × 64 information bits. Point estimate: SOVA makes more
    # information-bit errors and its soft mutual information is lower
    # (over-confident gaps). Clopper–Pearson α = 1e-4 intervals still overlap,
    # so this is not a claim of a significant separation.
    low = _sova_bcjr_counts(0.0, 48, 64)
    lo_s, hi_s = clopper_pearson(low.es, low.nbits, BER_ALPHA)
    lo_b, hi_b = clopper_pearson(low.eb, low.nbits, BER_ALPHA)
    @test low.es > low.eb
    @test low.mib > low.mis + 0.03
    @test lo_s < hi_b && lo_b < hi_s
    # Slicing the coded extrinsic alone is not the survivor label: where the
    # gap is smaller than the channel LLR, L_e(c) can point the other way.
    @test low.coded_ext > low.coded_app

    # 2 dB, 40 × 64. Both upper bounds sit under uncoded BPSK. The point
    # estimates are close and may rank SOVA slightly better; the intervals
    # overlap. Termination makes SOVA and hard Viterbi disagree on a bit.
    mid = _sova_bcjr_counts(2.0, 40, 64)
    lo_sm, hi_sm = clopper_pearson(mid.es, mid.nbits, BER_ALPHA)
    lo_bm, hi_bm = clopper_pearson(mid.eb, mid.nbits, BER_ALPHA)
    uncoded = ber_bpsk_awgn(10^(2.0 / 10))
    @test hi_sm < uncoded
    @test hi_bm < uncoded
    @test lo_sm < hi_bm && lo_bm < hi_sm
    @test mid.coded_ext > mid.coded_app
    @test mid.es < 40 && mid.eb < 40

    # 6 dB, 40 × 64. High SNR: information-bit hard decisions match Viterbi
    # and BCJR, and match the transmitted bits. Coded APP matches the codeword.
    high = _sova_bcjr_counts(6.0, 40, 64)
    @test high.es == 0 && high.eb == 0
    @test high.vit_dis == 0 && high.hard_dis == 0
    @test high.coded_app == 0 && high.coded_ext == 0
    @test high.mis > 0.99 && high.mib > 0.99
end
