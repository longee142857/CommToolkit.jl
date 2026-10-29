@testset "GF(2^m)" begin
    # α = x = 2. α^m equals the low part of the primitive polynomial.
    F2 = GF2m(2)
    @test sprint(show, F2) == "GF(2^2)"
    @test gf_pow(F2, 2, 2) == 3          # α^2 = α + 1
    @test gf_pow(F2, 2, 3) == 1
    F3 = GF2m(3)
    @test gf_pow(F3, 2, 3) == 3          # α^3 = α + 1
    @test gf_pow(F3, 2, 7) == 1
    F4 = GF2m(4)
    @test gf_pow(F4, 2, 4) == 3          # α^4 = α + 1
    @test gf_pow(F4, 2, 15) == 1
    # α^4 = α+1, so (α+1)^2 = α^8.
    @test gf_mul(F4, gf_add(F4, 2, 1), gf_add(F4, 2, 1)) == gf_pow(F4, 2, 8)

    for m in (2, 3, 4)
        F = GF2m(m)
        n = 2^m - 1
        @test gf_pow(F, 2, n) == 1
        seen = Int[]
        for e in 0:(n - 1)
            push!(seen, gf_pow(F, 2, e))
        end
        @test sort(seen) == collect(1:n)
        for a in 0:n
            @test gf_add(F, a, a) == 0
            @test gf_mul(F, a, 0) == 0
            @test gf_mul(F, a, 1) == a
            a == 0 && continue
            @test gf_mul(F, a, gf_inv(F, a)) == 1
            @test gf_div(F, a, a) == 1
            @test gf_pow(F, a, -1) == gf_inv(F, a)
            @test gf_pow(F, a, n) == 1
            b = mod(a * 3 + 1, n + 1)
            @test gf_mul(F, gf_add(F, a, b), gf_add(F, a, b)) ==
                  gf_add(F, gf_mul(F, a, a), gf_mul(F, b, b))
        end
        @test gf_pow(F, 0, 0) == 1
        @test gf_pow(F, 0, 4) == 0
        @test gf_div(F, 0, 1) == 0
    end
    @test_throws ArgumentError GF2m(1)
    @test_throws ArgumentError GF2m(5)
    @test_throws ArgumentError GF2m(8)
    F = GF2m(4)
    @test_throws ArgumentError gf_inv(F, 0)
    @test_throws ArgumentError gf_div(F, 1, 0)
    @test_throws ArgumentError gf_pow(F, 0, -1)
    @test_throws ArgumentError gf_mul(F, -1, 1)
    @test_throws ArgumentError gf_add(F, 16, 1)
end

@testset "Hamming(7, 4) hard decision" begin
    # Lin and Costello column order, data in positions 3, 5, 6, 7.
    @test hamming74_encode([1, 0, 1, 1]) == [0, 1, 1, 0, 0, 1, 1]
    @test hamming74_encode([0, 0, 0, 0]) == [0, 0, 0, 0, 0, 0, 0]
    @test hamming74_decode([0, 1, 1, 0, 0, 1, 1]) == [1, 0, 1, 1]

    code = hamming74_ldpc()
    words = Set{Vector{Int}}()
    ldpc_words = Set{Vector{Int}}()
    for mask in 0:15
        data = [(mask >> j) & 1 for j in 0:3]
        c = hamming74_encode(data)
        @test length(c) == 7
        @test hamming74_decode(c) == data
        @test ldpc_syndrome(code, c) == [0, 0, 0]
        push!(words, c)
        push!(ldpc_words, ldpc_encode(code, data))
        for i in 1:7
            rx = copy(c)
            rx[i] ⊻= 1
            @test hamming74_decode(rx) == data
        end
    end
    @test length(words) == 16
    @test words == ldpc_words

    # Two errors: syndrome points at a third column, so the data changes.
    rx2 = hamming74_encode([0, 0, 0, 0])
    rx2[1] ⊻= 1
    rx2[2] ⊻= 1
    @test hamming74_decode(rx2) == [1, 0, 0, 0]
    @test hamming74_decode(rx2) != [0, 0, 0, 0]

    @test_throws ArgumentError hamming74_encode([0, 1, 1])
    @test_throws ArgumentError hamming74_encode([0, 1, 2, 0])
    @test_throws ArgumentError hamming74_decode([0, 0, 0, 0, 0, 0])
    @test_throws ArgumentError hamming74_decode([0, 0, 0, 0, 0, 0, 2])
end

@testset "BCH(15, 7)" begin
    F = GF2m(4)
    α = 2
    # g(α^j) = 0 for the four consecutive designed roots.
    g = [1, 1, 1, 0, 1, 0, 0, 0, 1]
    for j in 1:4
        @test CommToolkit.Coding._horner(F, g, gf_pow(F, α, j)) == 0
    end

    # Generator, high-degree coefficients first, is a codeword.
    # x^8+x^7+x^6+x^4+1 as a length-15 vector.
    gen_word = [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1]
    @test bch157_encode([0, 0, 0, 0, 0, 0, 0]) == zeros(Int, 15)
    @test bch157_decode(gen_word) == bch157_decode(bch157_encode(bch157_decode(gen_word)))
    cw_g = bch157_encode(bch157_decode(gen_word))
    @test cw_g == gen_word
    for j in 1:4
        @test CommToolkit.Coding._horner(F, gen_word, gf_pow(F, α, j)) == 0
        @test CommToolkit.Coding._horner(F, cw_g, gf_pow(F, α, j)) == 0
    end

    # Binary syndrome identity S2 = S1^2 on a non-codeword.
    scratch = [1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1]
    S1 = CommToolkit.Coding._horner(F, scratch, α)
    S2 = CommToolkit.Coding._horner(F, scratch, gf_pow(F, α, 2))
    @test S2 == gf_mul(F, S1, S1)

    # Every message, every pattern of at most t = 2 errors.
    messages = Vector{Vector{Int}}()
    for mask in 0:127
        push!(messages, [(mask >> j) & 1 for j in 0:6])
    end
    for msg in messages
        cw = bch157_encode(msg)
        @test cw[1:7] == msg
        @test bch157_decode(cw) == msg
        for j in 1:4
            @test CommToolkit.Coding._horner(F, cw, gf_pow(F, α, j)) == 0
        end
    end
    # Weight 1 and 2 on every codeword is 128 * (15 + 105) = 15360 decodes.
    for msg in messages
        cw = bch157_encode(msg)
        for i in 1:15
            rx = copy(cw)
            rx[i] ⊻= 1
            @test bch157_decode(rx) == msg
        end
        for i in 1:15, j in (i + 1):15
            rx = copy(cw)
            rx[i] ⊻= 1
            rx[j] ⊻= 1
            @test bch157_decode(rx) == msg
        end
    end

    # A published-shape vector: the all-zero word with errors in positions 1 and 8
    # (coefficients of x^14 and x^7) comes back to the zero message.
    rx = zeros(Int, 15)
    rx[1] = 1
    rx[8] = 1
    @test bch157_decode(rx) == zeros(Int, 7)

    # Three errors are outside the designed radius. This pattern is not corrected.
    rx3 = zeros(Int, 15)
    rx3[1] = 1
    rx3[2] = 1
    rx3[3] = 1
    @test bch157_decode(rx3) != zeros(Int, 7)

    # Not every weight-3 pattern of the zero word is restored. d = 5, t = 2.
    n_miss = 0
    n3 = 0
    for i in 1:15, j in (i + 1):15, k in (j + 1):15
        rx = zeros(Int, 15)
        rx[i] = 1
        rx[j] = 1
        rx[k] = 1
        n3 += 1
        if bch157_decode(rx) != zeros(Int, 7)
            n_miss += 1
        end
    end
    @test n3 == binomial(15, 3)
    @test n_miss > 0

    @test_throws ArgumentError bch157_encode([0, 0, 0, 0, 0, 0])
    @test_throws ArgumentError bch157_encode([0, 0, 0, 0, 0, 0, 2])
    @test_throws ArgumentError bch157_decode(zeros(Int, 14))
    @test_throws ArgumentError bch157_decode(fill(2, 15))
end

@testset "Reed–Solomon constructors and reference vectors" begin
    C1511 = rs_code(4, 11)
    @test C1511 isa RSCode
    @test sprint(show, C1511) == "RS(15, 11) over GF(2^4)"
    @test (C1511.n, C1511.k, C1511.t) == (15, 11, 2)
    # Lin and Costello, 2nd ed., ch. 7: roots α^1..α^4 of p(x) = x^4+x+1 give
    # g(x) = x^4 + α^13 x^3 + α^6 x^2 + α^3 x + α^10 (α^13 = 13, α^6 = 12,
    # α^3 = 8, α^10 = 7 as integer masks).
    @test C1511.g == [1, 13, 12, 8, 7]
    C159 = rs_code(4, 9)
    @test (C159.n, C159.k, C159.t) == (15, 9, 3)
    # g(x) = x^6 + α^10 x^5 + α^14 x^4 + α^4 x^3 + α^6 x^2 + α^9 x + α^6.
    @test C159.g == [1, 7, 9, 3, 12, 10, 12]
    C73 = rs_code(3, 3)
    @test (C73.n, C73.k, C73.t) == (7, 3, 2)

    # Every designed root α^j is a root of the generator.
    F = GF2m(4)
    for (C, d) in ((C1511, 4), (C159, 6))
        for j in 1:d
            @test CommToolkit.Coding._gf_eval(F, C.g, gf_pow(F, 2, j)) == 0
        end
    end

    # Frozen systematic vectors. c(x) = x^14 + x^13 + ... + 1 has roots
    # α^1..α^4, so eleven 1s encode to fifteen 1s.
    @test rs_encode(C1511, ones(Int, 11)) == ones(Int, 15)
    @test rs_encode(C1511, collect(1:11)) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 10, 14, 6]
    @test rs_encode(C159, ones(Int, 9)) == ones(Int, 15)
    @test rs_encode(C73, [1, 2, 3]) == [1, 2, 3, 0, 0, 1, 3]

    # RS(3, 1) over GF(2^2) is the repetition code: g(x) = x^2 + x + 1.
    C31 = rs_code(2, 1)
    @test C31.g == [1, 1, 1]
    for a in 0:3
        @test rs_encode(C31, [a]) == [a, a, a]
        for i in 1:3, v in 1:3
            rx = [a, a, a]
            rx[i] ⊻= v
            @test rs_decode(C31, rx) == [a]
        end
    end

    @test_throws ArgumentError rs_code(4, 0)
    @test_throws ArgumentError rs_code(4, 15)
    @test_throws ArgumentError rs_code(5, 20)   # no GF(2^5) table
    @test_throws ArgumentError rs_encode(C1511, ones(Int, 10))
    @test_throws ArgumentError rs_encode(C1511, fill(16, 11))
    @test_throws ArgumentError rs_decode(C1511, ones(Int, 14))
    @test_throws ArgumentError rs_decode(C1511, fill(-1, 15))
end

@testset "RS(7, 3) exhaustive" begin
    C = rs_code(3, 3)
    F = GF2m(3)
    messages = Vector{Vector{Int}}()
    for m0 in 0:7, m1 in 0:7, m2 in 0:7
        push!(messages, [m0, m1, m2])
    end

    # Round trip, systematic order, zero syndromes, and the Singleton bound
    # with equality: the minimum weight of the code is d = n - k + 1 = 5.
    min_weight = C.n
    for msg in messages
        cw = rs_encode(C, msg)
        @test cw[1:3] == msg
        @test rs_decode(C, cw) == msg
        for j in 1:4
            @test CommToolkit.Coding._gf_eval(F, cw, gf_pow(F, 2, j)) == 0
        end
        w = count(!=(0), cw)
        w > 0 && (min_weight = min(min_weight, w))
    end
    @test min_weight == 5

    # Every message against every error pattern of weight at most t = 2:
    # position sets and nonzero error values are exhaustive.
    for msg in messages
        cw = rs_encode(C, msg)
        for i in 1:7, v1 in 1:7
            rx = copy(cw)
            rx[i] ⊻= v1
            @test rs_decode(C, rx) == msg
        end
        for i in 1:7, j in (i + 1):7, v1 in 1:7, v2 in 1:7
            rx = copy(cw)
            rx[i] ⊻= v1
            rx[j] ⊻= v2
            @test rs_decode(C, rx) == msg
        end
    end

    # Beyond t. Every weight-3 pattern on the zero word is either returned
    # unchanged or miscorrected to a codeword inside the radius-t ball of the
    # received word. It is never "corrected" to the zero message; the zero
    # message only appears through the unchanged path when all three errors
    # sit in the parity symbols.
    n_unchanged = 0
    n_miscorrect = 0
    for i in 1:7, j in (i + 1):7, k in (j + 1):7, v1 in 1:7, v2 in 1:7, v3 in 1:7
        rx = zeros(Int, 7)
        rx[i] = v1
        rx[j] = v2
        rx[k] = v3
        dec = rs_decode(C, rx)
        if dec == rx[1:3]
            n_unchanged += 1
        else
            n_miscorrect += 1
            @test sum(rs_encode(C, dec) .!= rx) <= C.t
            @test dec != zeros(Int, 3)
        end
    end
    @test n_unchanged + n_miscorrect == binomial(7, 3) * 7^3
    @test n_unchanged > 0
    @test n_miscorrect > 0

    # Frozen decode vectors: one unchanged return, one miscorrection.
    @test rs_decode(C, [1, 1, 1, 0, 0, 0, 0]) == [1, 1, 1]
    dec = rs_decode(C, [1, 3, 0, 2, 0, 0, 0])
    @test dec == [1, 3, 1]
    @test sum(rs_encode(C, dec) .!= [1, 3, 0, 2, 0, 0, 0]) <= C.t
end

@testset "RS(15, k) correction radius" begin
    Random.seed!(RNG_SEED)
    C1511 = rs_code(4, 11)
    messages = [zeros(Int, 11), ones(Int, 11),
                [rand(0:15) for _ in 1:11], [rand(0:15) for _ in 1:11],
                [rand(0:15) for _ in 1:11], [rand(0:15) for _ in 1:11],
                vcat([7], zeros(Int, 10)), vcat(zeros(Int, 10), [13])]
    # Weight at most t = 2, positions and values exhaustive on fixed messages.
    for msg in messages
        cw = rs_encode(C1511, msg)
        @test rs_decode(C1511, cw) == msg
        for i in 1:15, v1 in 1:15
            rx = copy(cw)
            rx[i] ⊻= v1
            @test rs_decode(C1511, rx) == msg
        end
        for i in 1:15, j in (i + 1):15, v1 in 1:15, v2 in 1:15
            rx = copy(cw)
            rx[i] ⊻= v1
            rx[j] ⊻= v2
            @test rs_decode(C1511, rx) == msg
        end
    end

    # t = 3: weight-3 patterns sampled under the fixed seed.
    C159 = rs_code(4, 9)
    Random.seed!(RNG_SEED)
    for _ in 1:20000
        msg = rand(0:15, 9)
        rx = copy(rs_encode(C159, msg))
        for p in sort(randperm(15)[1:3])
            rx[p] ⊻= rand(1:15)
        end
        @test rs_decode(C159, rx) == msg
    end

    # k = n - 1 gives t = 0: the decoder checks the syndrome but corrects
    # nothing, so a single error comes back unchanged.
    C1514 = rs_code(4, 14)
    @test C1514.t == 0
    msg = rand(0:15, 14)
    cw = rs_encode(C1514, msg)
    @test rs_decode(C1514, cw) == msg
    rx = copy(cw)
    rx[3] ⊻= 5
    @test rs_decode(C1514, rx) == rx[1:14]
    @test rs_decode(C1514, rx) != msg
end

@testset "Shortened Reed–Solomon" begin
    Random.seed!(RNG_SEED)
    C106 = rs_shorten(rs_code(4, 11), 5)
    @test sprint(show, C106) == "RS(10, 6) over GF(2^4), shortened by 5"
    @test (C106.n, C106.k, C106.t) == (10, 6, 2)
    @test C106.g == rs_code(4, 11).g
    # Shortening twice composes.
    C2 = rs_shorten(rs_shorten(rs_code(4, 11), 2), 3)
    @test (C2.n, C2.k, C2.t, C2.s) == (10, 6, 2, 5)

    F = GF2m(4)
    @test rs_encode(C106, collect(1:6)) == [1, 2, 3, 4, 5, 6, 13, 2, 4, 14]
    messages = [zeros(Int, 6), ones(Int, 6),
                [rand(0:15) for _ in 1:6], [rand(0:15) for _ in 1:6],
                [rand(0:15) for _ in 1:6], [rand(0:15) for _ in 1:6]]
    # Round trip, syndromes of the padded parent word, and exhaustive weight
    # at most t = 2 on the effective length.
    for msg in messages
        cw = rs_encode(C106, msg)
        @test length(cw) == 10
        @test cw[1:6] == msg
        @test rs_decode(C106, cw) == msg
        padded = vcat(zeros(Int, 5), cw)
        for j in 1:4
            @test CommToolkit.Coding._gf_eval(F, padded, gf_pow(F, 2, j)) == 0
        end
        for i in 1:10, v1 in 1:15
            rx = copy(cw)
            rx[i] ⊻= v1
            @test rs_decode(C106, rx) == msg
        end
        for i in 1:10, j in (i + 1):10, v1 in 1:15, v2 in 1:15
            rx = copy(cw)
            rx[i] ⊻= v1
            rx[j] ⊻= v2
            @test rs_decode(C106, rx) == msg
        end
    end

    # Beyond t the shortened code keeps the same two failure modes.
    n_unchanged = 0
    for _ in 1:2000
        rx = rand(0:15, 10)
        dec = rs_decode(C106, rx)
        if dec == rx[1:6]
            n_unchanged += 1
        else
            @test sum(rs_encode(C106, dec) .!= rx) <= C106.t
        end
    end
    @test n_unchanged > 0

    @test_throws ArgumentError rs_shorten(rs_code(4, 11), -1)
    @test_throws ArgumentError rs_shorten(rs_code(4, 11), 11)
    @test_throws ArgumentError rs_encode(C106, ones(Int, 7))
    @test_throws ArgumentError rs_decode(C106, ones(Int, 11))
end
