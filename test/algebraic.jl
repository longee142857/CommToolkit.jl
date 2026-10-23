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
