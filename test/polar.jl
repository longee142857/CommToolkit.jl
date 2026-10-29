# Teaching Arikan polar code checks.
#
# Construction: Arikan, IEEE Trans. Inf. Theory 55(7):3051–3073, 2009 (channel
#   polarization; Proposition 5 for the Bhattacharyya recursion Z^- = 2Z − Z²,
#   Z^+ = Z², exact on the BEC and a teaching heuristic on AWGN).
# Systematic encoding: Arikan, IEEE Commun. Lett. 15(8):860–862, 2011 (message
#   read from fixed codeword positions), written for the natural-order
#   transform, whose principal submatrices are triangular with unit diagonal.
# Decoding: successive cancellation on the transform tree; the box-plus is the
#   exact log-domain check update shared with the LDPC decoder.
# Noise uses the same σ² = N0/2 as awgn / demodulate_soft. Not a 3GPP NR
# polar code: no CRC, no list decoding, no rate matching.

@testset "Polar construction" begin
    # Frozen construction vectors (deterministic, no randomness involved).
    c84 = polar_code(8, 4)
    @test c84.n == 8 && c84.k == 4 && c84.design_snr_db == 0.0
    @test c84.info_pos == [4, 6, 7, 8]
    @test c84.frozen == BitVector([1, 1, 1, 0, 1, 0, 0, 0])
    @test polar_code(16, 8; design_snr_db=2.0).info_pos ==
          [8, 10, 11, 12, 13, 14, 15, 16]
    @test polar_code(16, 5; design_snr_db=-1.0).info_pos == [8, 12, 14, 15, 16]
    @test polar_code(32, 16; design_snr_db=2.0).info_pos ==
          [12, 14, 15, 16, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]

    # Bhattacharyya recursion frozen at design 0 dB, N = 4: Z0 = e⁻¹,
    # z = [M(M(Z0)), P(M(Z0)), M(P(Z0)), P(P(Z0))], M(z) = 2z − z², P(z) = z².
    c4 = polar_code(4, 2)
    z0 = exp(-1.0)
    zm = 2 * z0 - z0^2
    zp = z0^2
    @test c4.z ≈ [2 * zm - zm^2, zm^2, 2 * zp - zp^2, zp^2] rtol = 1e-12
    @test c4.z ≈ [0.8403386998488147, 0.36050849836372917,
                  0.25235492758449124, 0.018315638888734182] rtol = 1e-12
    @test c4.info_pos == [3, 4]

    # The information set is the k smallest z; ties keep the lower index.
    for (n, k, dsnr) in ((8, 4, 0.0), (16, 8, 2.0), (64, 32, 2.0), (256, 128, 2.0))
        c = polar_code(n, k; design_snr_db=dsnr)
        zmax_in = maximum(c.z[c.info_pos])
        zmin_out = minimum(c.z[c.frozen])
        @test zmax_in <= zmin_out
        @test length(c.info_pos) == k
        @test issorted(c.info_pos)
        @test all(!iszero, c.z)
        # Reproducible: same arguments, same construction.
        c2 = polar_code(n, k; design_snr_db=dsnr)
        @test c2.info_pos == c.info_pos
        @test c2.z == c.z
        @test c2.frozen == c.frozen
    end

    # Rate-1 degenerate code is legal: nothing is frozen.
    call = polar_code(8, 8)
    @test call.info_pos == collect(1:8)
    @test !any(call.frozen)

    # show() smoke.
    @test sprint(show, c84) == "PolarCode(n=8, k=4, design_snr_db=0.0)"
    @test sprint(show, MIME"text/plain"(), c84) == sprint(show, c84)

    @test_throws ArgumentError polar_code(3, 2)          # not a power of two
    @test_throws ArgumentError polar_code(12, 4)         # not a power of two
    @test_throws ArgumentError polar_code(0, 1)
    @test_throws ArgumentError polar_code(512, 4)        # teaching cap is 256
    @test_throws ArgumentError polar_code(8, 0)
    @test_throws ArgumentError polar_code(8, 9)
    @test_throws ArgumentError polar_code(8, 4; design_snr_db=NaN)
    @test_throws ArgumentError polar_code(8, 4; design_snr_db=Inf)
end

@testset "Polar transform and encoding" begin
    # N = 2 by hand: x = u F, F = [1 0; 1 1], row-vector convention.
    c22 = polar_code(2, 2)
    @test polar_encode(c22, [0, 0]) == [0, 0]
    @test polar_encode(c22, [1, 0]) == [1, 0]
    @test polar_encode(c22, [0, 1]) == [1, 1]
    @test polar_encode(c22, [1, 1]) == [0, 1]

    # N = 4 rate 1: basis vectors map to the rows of F ⊗ F.
    c44 = polar_code(4, 4)
    @test polar_encode(c44, [1, 0, 0, 0]) == [1, 0, 0, 0]
    @test polar_encode(c44, [0, 1, 0, 0]) == [1, 1, 0, 0]
    @test polar_encode(c44, [0, 0, 1, 0]) == [1, 0, 1, 0]
    @test polar_encode(c44, [0, 0, 0, 1]) == [1, 1, 1, 1]

    # The butterfly is an involution, G_N G_N = I over GF(2).
    Random.seed!(RNG_SEED)
    for n in (2, 4, 8, 16, 32, 64, 128, 256)
        u = rand(0:1, n)
        v = CommToolkit.Coding._polar_transform!(copy(u))
        @test CommToolkit.Coding._polar_transform!(v) == u
    end

    # Linearity of the non-systematic map: E(a) ⊕ E(b) = E(a ⊕ b).
    c = polar_code(16, 8; design_snr_db=2.0)
    for _ in 1:20
        a = rand(0:1, c.k)
        b = rand(0:1, c.k)
        @test polar_encode(c, a) .⊻ polar_encode(c, b) == polar_encode(c, a .⊻ b)
    end

    # Frozen encode vectors (default design 0 dB).
    c84 = polar_code(8, 4)
    @test polar_encode(c84, [1, 0, 1, 1]) == [1, 0, 1, 0, 0, 1, 0, 1]
    @test polar_encode(c84, [1, 0, 1, 1]; systematic=true) == [0, 0, 1, 1, 0, 0, 1, 1]

    # Systematic: the codeword carries the message on info_pos, for a seeded
    # sweep of constructions and messages. The solved pre-image stays zero on
    # the frozen positions, so systematic and non-systematic words span the
    # same linear code.
    Random.seed!(RNG_SEED)
    for n in (2, 4, 8, 16, 32, 64, 256)
        for k in unique!([1, max(1, n ÷ 4), n ÷ 2, n - 1, n])
            for dsnr in (-2.0, 0.0, 2.0)
                code = polar_code(n, k; design_snr_db=dsnr)
                for _ in 1:3
                    msg = rand(0:1, k)
                    x = polar_encode(code, msg; systematic=true)
                    @test x[code.info_pos] == msg
                    llr = 30.0 .* (1 .- 2 .* x)
                    @test polar_sc_decode(code, llr; systematic=true) == msg
                end
            end
        end
    end

    # Same code: exhaustive word sets of the two maps coincide (N = 8, k = 3).
    c83 = polar_code(8, 3)
    words_ns = Set{Vector{Int}}()
    words_sy = Set{Vector{Int}}()
    for m in 0:7
        msg = [(m >> (3 - t)) & 1 for t in 1:3]
        push!(words_ns, polar_encode(c83, msg))
        push!(words_sy, polar_encode(c83, msg; systematic=true))
    end
    @test words_ns == words_sy
    @test length(words_ns) == 8

    # Minimum distance by exhaustive enumeration: d_min is the smallest row
    # weight 2^popcount(i − 1) over the information rows.
    for (n, k, dsnr, dmin) in ((8, 3, 0.0, 4), (8, 4, 0.0, 4), (8, 5, 2.0, 2),
                               (16, 8, 2.0, 4), (16, 4, -1.0, 8))
        code = polar_code(n, k; design_snr_db=dsnr)
        w_min = n
        for m in 1:(2^k - 1)
            msg = [(m >> (k - t)) & 1 for t in 1:k]
            w_min = min(w_min, sum(polar_encode(code, msg)))
        end
        @test w_min == dmin
        @test dmin == 2^minimum(count_ones(i - 1) for i in code.info_pos)
    end

    @test_throws ArgumentError polar_encode(c83, [1, 0])
    @test_throws ArgumentError polar_encode(c83, [1, 0, 2])
    @test_throws ArgumentError polar_encode(c83, [1, 0, -1])
end

@testset "Polar SC roundtrip" begin
    # Exhaustive noiseless roundtrip, N = 8, every k, every message, both maps.
    for k in 1:8
        code = polar_code(8, k)
        for m in 0:(2^k - 1)
            msg = [(m >> (k - t)) & 1 for t in 1:k]
            for sys in (false, true)
                x = polar_encode(code, msg; systematic=sys)
                llr = 30.0 .* (1 .- 2 .* x)
                @test polar_sc_decode(code, llr; systematic=sys) == msg
            end
        end
    end

    # Noiseless roundtrip on larger frames, seeded messages, both maps.
    Random.seed!(RNG_SEED)
    for n in (16, 32, 64, 128, 256)
        for k in unique!([1, n ÷ 2, n - 1])
            code = polar_code(n, k; design_snr_db=2.0)
            for _ in 1:3
                msg = rand(0:1, k)
                for sys in (false, true)
                    x = polar_encode(code, msg; systematic=sys)
                    @test polar_sc_decode(code, 30.0 .* (1 .- 2 .* x);
                                          systematic=sys) == msg
                end
            end
        end
    end

    # Through the actual BPSK + awgn chain at +Inf: no noise, LLR sign is exact.
    c = polar_code(64, 32; design_snr_db=2.0)
    Random.seed!(RNG_SEED)
    msg = rand(0:1, c.k)
    x = polar_encode(c, msg)
    y = awgn(modulate(x, BPSK()), Inf; energy_per_bit=c.n / c.k)
    @test polar_sc_decode(c, demodulate_soft(y, BPSK(), 1e-6)) == msg

    # Tie rule: an all-zero LLR vector decodes to the all-zero message.
    @test polar_sc_decode(c, zeros(c.n)) == zeros(Int, c.k)
    @test polar_sc_decode(c, zeros(c.n); systematic=true) == zeros(Int, c.k)

    @test_throws ArgumentError polar_sc_decode(c, zeros(c.n + 1))
    @test_throws ArgumentError polar_sc_decode(c, [0.1, NaN, fill(0.0, c.n - 2)...])
    @test_throws ArgumentError polar_sc_decode(c, [0.1, Inf, fill(0.0, c.n - 2)...])
    @test_throws ArgumentError polar_sc_decode(c, [0.1, -Inf, fill(0.0, c.n - 2)...])
    @test_throws ArgumentError polar_sc_decode(c, zeros(Int, c.n)) # hard bits
end

@testset "Polar repetition anchors" begin
    # polar_code(2, 1) and polar_code(4, 1) put the single information bit on
    # the best index, which generates the length-n repetition code: x = (d, …, d).
    c2 = polar_code(2, 1)
    c4 = polar_code(4, 1)
    @test c2.info_pos == [2]
    @test c4.info_pos == [4]
    @test polar_encode(c2, [1]) == [1, 1]
    @test polar_encode(c4, [1]) == [1, 1, 1, 1]
    @test polar_encode(c2, [1]; systematic=true) == [1, 1]
    @test polar_encode(c4, [1]; systematic=true) == [1, 1, 1, 1]

    # SC on the repetition code sums the channel LLRs (frozen partial sums are
    # zero), so it is exactly ML: BER = Q(√(2 N R γb)) with R = 1/n, i.e.
    # Q(√(2 γb)) = ber_bpsk_awgn. A repetition code has no AWGN coding gain;
    # this anchors the decoder and the energy bookkeeping to a closed form.
    for (code, snr_db, frames) in ((c2, 2.0, 8000), (c4, 3.0, 12000))
        γ = 10^(snr_db / 10)
        Rinv = code.n / code.k
        σ2 = Rinv / (2 * γ)
        bit_errors = 0
        Random.seed!(RNG_SEED)
        for _ in 1:frames
            msg = rand(0:1, code.k)
            x = polar_encode(code, msg)
            y = awgn(modulate(x, BPSK()), snr_db; energy_per_bit=Rinv)
            hat = polar_sc_decode(code, demodulate_soft(y, BPSK(), σ2))
            bit_errors += count(!=(0), msg .- hat)
        end
        nbits = frames * code.k
        println("Polar repetition n=$(code.n) Eb/N0=$(snr_db) dB " *
                "bits=$(bit_errors)/$(nbits) theory=$(ber_bpsk_awgn(γ))")
        assert_theory_inside(bit_errors, nbits, ber_bpsk_awgn(γ))
    end
end

@testset "Polar AWGN coding gain" begin
    # N = 64, k = 32, designed at 2 dB. Es = 1, Eb = n/k, σ² = Eb/(2 Eb/N0),
    # same call as awgn. 1000 frames per point, seed RNG_SEED, α = 1e-4.
    code = polar_code(64, 32; design_snr_db=2.0)
    nframes = 1000
    Rinv = code.n / code.k
    @test Rinv == 2.0
    intervals = Dict{Float64,NTuple{2,Float64}}()
    Random.seed!(RNG_SEED)
    for snr_db in (2.0, 3.0)
        γ = 10^(snr_db / 10)
        σ2 = Rinv / (2 * γ)
        bit_errors = 0
        block_errors = 0
        for _ in 1:nframes
            bits = rand(0:1, code.k)
            word = polar_encode(code, bits)
            y = awgn(modulate(word, BPSK()), snr_db; energy_per_bit=Rinv)
            hat = polar_sc_decode(code, demodulate_soft(y, BPSK(), σ2))
            errors = count(!=(0), bits .- hat)
            bit_errors += errors
            block_errors += errors > 0
        end
        nbits = nframes * code.k
        p_bit = ber_bpsk_awgn(γ)
        p_block = 1 - (1 - p_bit)^code.k
        lo_b, hi_b = clopper_pearson(bit_errors, nbits, BER_ALPHA)
        lo_B, hi_B = clopper_pearson(block_errors, nframes, BER_ALPHA)
        println("Polar (64, 32) SC Eb/N0=$(snr_db) dB " *
                "bits=$(bit_errors)/$(nbits) BER_CP=($(lo_b), $(hi_b)) " *
                "uncoded_BER=$(p_bit) blocks=$(block_errors)/$(nframes) " *
                "BLER_CP=($(lo_B), $(hi_B)) uncoded_BLER=$(p_block)")
        @test bit_errors >= 15
        @test block_errors >= 4
        @test hi_b < p_bit
        @test hi_B < p_block
        intervals[snr_db] = (lo_b, hi_b)
    end
    # The two operating points separate: more energy, strictly fewer errors.
    @test intervals[3.0][2] < intervals[2.0][1]
end
