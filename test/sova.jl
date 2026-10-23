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
