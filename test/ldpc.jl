# Teaching LDPC checks.
#
# Constructions: Gallager, IRE Trans. Inf. Theory 8(1):21–28, 1962, Section II;
#   Fan, Proc. 2nd Int. Symp. Turbo Codes, Brest, 2000, pp. 543–546;
#   Lin and Costello, Error Control Coding, 2nd ed., §3.2 (Hamming columns and
#   systematic Gaussian elimination).
# Decoding: Kschischang, Frey, and Loeliger, IEEE Trans. Inf. Theory
#   47(2):498–519, 2001; Fossorier, Mihaljević, and Imai, IEEE Trans. Commun.
#   47(5):673–680, 1999; Chen and Fossorier, IEEE Trans. Commun.
#   50(3):406–414, 2002; Chen and Fossorier, IEEE Commun. Lett. 6(5):208–210, 2002.
# The [7,4] Hamming matrix is the length-7 fixture in Lin and Costello §3.2.
# Noise uses the same σ² = N0/2 as awgn / demodulate_soft.

using SparseArrays

const _LDPC_ALGS = (:sum_product, :min_sum, :normalized_min_sum, :offset_min_sum)

function _has_4_cycle(code::LDPCCode)
    rows = [Int[] for _ in 1:code.m]
    for e in eachindex(code.edge_var)
        push!(rows[code.edge_check[e]], code.edge_var[e])
    end
    for i in 1:code.m
        for j in (i + 1):code.m
            inter = 0
            for a in rows[i]
                for b in rows[j]
                    a == b && (inter += 1)
                end
            end
            inter >= 2 && return true
        end
    end
    return false
end

function _kernel_words(code::LDPCCode)
    code.n <= 16 || throw(ArgumentError("kernel enumeration is only for n ≤ 16"))
    words = Vector{Int}[]
    buf = zeros(Int, code.n)
    for mask in 0:((1 << code.n) - 1)
        for j in 1:code.n
            buf[j] = (mask >> (j - 1)) & 1
        end
        all(iszero, ldpc_syndrome(code, buf)) && push!(words, copy(buf))
    end
    return words
end

function _logaddexp(a::Float64, b::Float64)
    a == -Inf && return b
    b == -Inf && return a
    return a >= b ? a + log1p(exp(b - a)) : b + log1p(exp(a - b))
end

"""Bitwise MAP LLRs by enumerating the kernel. Independent of the encoder."""
function _map_llr(code::LDPCCode, llr::AbstractVector{<:Real})
    n = code.n
    pos = fill(-Inf, n)
    neg = fill(-Inf, n)
    buf = zeros(Int, n)
    for mask in 0:((1 << n) - 1)
        for j in 1:n
            buf[j] = (mask >> (j - 1)) & 1
        end
        all(iszero, ldpc_syndrome(code, buf)) || continue
        logp = 0.0
        for j in 1:n
            buf[j] == 1 && (logp -= float(llr[j]))
        end
        for j in 1:n
            if buf[j] == 0
                pos[j] = _logaddexp(pos[j], logp)
            else
                neg[j] = _logaddexp(neg[j], logp)
            end
        end
    end
    return pos .- neg
end

"""One-check posterior from the tanh product or the independent min formula."""
function _spc_posterior(llr::Vector{Float64}, algorithm::Symbol;
                        alpha::Float64=1.0, beta::Float64=0.0)
    n = length(llr)
    out = Vector{Float64}(undef, n)
    for i in 1:n
        if algorithm === :sum_product
            prod_tanh = 1.0
            for j in 1:n
                j == i && continue
                prod_tanh *= tanh(llr[j] / 2)
            end
            msg = 2 * atanh(clamp(prod_tanh, -1 + 1e-15, 1 - 1e-15))
        else
            mag = Inf
            neg = false
            for j in 1:n
                j == i && continue
                mag = min(mag, abs(llr[j]))
                llr[j] < 0 && (neg = !neg)
            end
            if algorithm === :offset_min_sum
                mag = mag > beta ? mag - beta : 0.0
            end
            if algorithm === :normalized_min_sum
                mag *= alpha
            end
            msg = neg ? -mag : mag
        end
        out[i] = llr[i] + msg
    end
    return out
end

function _assert_code_structure(code::LDPCCode)
    @test code.k >= 1
    @test code.k + code.rank == code.n
    @test code.m >= code.rank
    @test length(code.info_pos) == code.k
    @test length(code.parity_pos) == code.rank
    @test issorted(code.info_pos)
    @test issorted(code.parity_pos)
    @test sort(vcat(code.info_pos, code.parity_pos)) == collect(1:code.n)
    @test allunique(code.info_pos)
    @test allunique(code.parity_pos)
    col = ldpc_column_degrees(code)
    row = ldpc_row_degrees(code)
    @test length(col) == code.n
    @test length(row) == code.m
    @test sum(col) == sum(row) == length(code.edge_var)
    if code.column_weight > 0
        @test all(==(code.column_weight), col)
    end
    if code.row_weight > 0
        @test all(==(code.row_weight), row)
    end
    G = ldpc_generator_matrix(code)
    @test size(G) == (code.k, code.n)
    @test G[:, code.info_pos] == Matrix{Int}(LinearAlgebra.I, code.k, code.k)
    @test G[:, code.parity_pos] == Int.(transpose(code.P))
    HG = Matrix(code.H) * transpose(G)
    @test all(iszero, HG .% 2)
    for i in 1:code.k
        unit = zeros(Int, code.k)
        unit[i] = 1
        word = ldpc_encode(code, unit)
        @test word == vec(G[i, :])
        @test all(iszero, ldpc_syndrome(code, word))
    end
    Random.seed!(11)
    u = rand(0:1, code.k)
    v = rand(0:1, code.k)
    @test ldpc_encode(code, u .⊻ v) ==
          ldpc_encode(code, u) .⊻ ldpc_encode(code, v)
    @test all(iszero, ldpc_syndrome(code, ldpc_encode(code, u)))
    return nothing
end

function _assert_encoder_is_kernel(code::LDPCCode)
    words = _kernel_words(code)
    @test length(words) == 1 << code.k
    image = Set{Vector{Int}}()
    for mask in 0:((1 << code.k) - 1)
        bits = [(mask >> (t - 1)) & 1 for t in 1:code.k]
        push!(image, ldpc_encode(code, bits))
    end
    @test image == Set(words)
    weights = [count(isone, w) for w in words if !all(iszero, w)]
    return minimum(weights)
end

@testset "LDPC" begin
    @testset "Structure and reference matrices" begin
        h = hamming74_ldpc()
        @test h.n == 7 && h.k == 4 && h.m == 3 && h.rank == 3
        @test h.info_pos == [3, 5, 6, 7]
        @test h.parity_pos == [1, 2, 4]
        @test h.column_weight == 0
        @test h.row_weight == 4
        @test ldpc_column_degrees(h) == [1, 1, 2, 1, 2, 2, 3]
        @test ldpc_row_degrees(h) == [4, 4, 4]
        @test Matrix(h.H) == [
            1 0 1 0 1 0 1
            0 1 1 0 0 1 1
            0 0 0 1 1 1 1
        ]
        @test ldpc_encode(h, [1, 0, 0, 0]) == [1, 1, 1, 0, 0, 0, 0]
        @test ldpc_encode(h, [0, 0, 0, 1]) == [1, 1, 0, 1, 0, 0, 1]
        _assert_code_structure(h)
        @test _assert_encoder_is_kernel(h) == 3
        @test _has_4_cycle(h)
        @test sprint(show, h) == "LDPCCode(n=7, k=4, m=3, rank=3)"
        @test sprint(show, MIME"text/plain"(), h) == sprint(show, h)

        arr = array_ldpc(2, 3, 3)
        @test arr.n == 9 && arr.k == 4 && arr.m == 6 && arr.rank == 5
        @test arr.column_weight == 2 && arr.row_weight == 3
        @test Matrix(arr.H) == [
            1 0 0 1 0 0 1 0 0
            0 1 0 0 1 0 0 1 0
            0 0 1 0 0 1 0 0 1
            1 0 0 0 1 0 0 0 1
            0 1 0 0 0 1 1 0 0
            0 0 1 1 0 0 0 1 0
        ]
        @test !_has_4_cycle(arr)
        _assert_code_structure(arr)
        @test _assert_encoder_is_kernel(arr) >= 2
        @test sprint(show, arr) == "LDPCCode(n=9, k=4, m=6, rank=5, regular=(2,3))"
        @test sprint(show, MIME"text/plain"(), arr) == sprint(show, arr)

        g = gallager_ldpc(6, 2, 3)
        g2 = gallager_ldpc(6, 2, 3)
        @test g.H == g2.H
        @test g.n == 6 && g.k == 3 && g.m == 4 && g.rank == 3
        @test g.column_weight == 2 && g.row_weight == 3
        @test Matrix(g.H) == [
            1 1 1 0 0 0
            0 0 0 1 1 1
            1 0 1 0 0 1
            0 1 0 1 1 0
        ]
        @test _has_4_cycle(g)
        _assert_code_structure(g)
        @test _assert_encoder_is_kernel(g) >= 2
        @test sprint(show, g) == "LDPCCode(n=6, k=3, m=4, rank=3, regular=(2,3))"

        Random.seed!(1)
        before = rand()
        Random.seed!(1)
        gallager_ldpc(12, 3, 4)
        @test rand() == before
        @test gallager_ldpc(12, 3, 4).H == gallager_ldpc(12, 3, 4).H

        a7 = array_ldpc(3, 6, 7)
        @test a7.n == 42 && a7.k == 23 && a7.m == 21 && a7.rank == 19
        @test !_has_4_cycle(a7)
        _assert_code_structure(a7)
        a11 = array_ldpc(3, 6, 11)
        @test a11.n == 66 && a11.k == 35 && a11.m == 33 && a11.rank == 31
        @test !_has_4_cycle(a11)
        _assert_code_structure(a11)

        Hz = Int8[1 1 0; 0 0 0; 1 0 1]
        cz = ldpc_from_H(Hz)
        @test cz.m == 3
        @test ldpc_row_degrees(cz) == [2, 0, 2]
        @test all(iszero, ldpc_syndrome(cz, ldpc_encode(cz, zeros(Int, cz.k))))
        _assert_code_structure(cz)
        isolated = ldpc_from_H([1 1 0; 0 0 0])
        @test ldpc_column_degrees(isolated)[3] == 0
        det = ldpc_decode_details(isolated, [1.5, -0.5, 0.25];
                                  algorithm=:sum_product, early_stop=false, max_iter=2)
        @test det.iterations == 2
        @test det.posterior[3] ≈ 0.25

        Hbool = [1 1 1; 1 0 1]
        from_int = ldpc_from_H(Hbool)
        from_bool = ldpc_from_H(Hbool .== 1)
        from_sparse = ldpc_from_H(sparse(Hbool))
        @test from_int.H == from_bool.H == from_sparse.H
        @test from_int.info_pos == from_bool.info_pos == from_sparse.info_pos
        swapped = ldpc_from_H([0 1 1; 1 0 1])
        _assert_code_structure(swapped)
    end

    @testset "Rejected arguments" begin
        @test_throws ArgumentError gallager_ldpc(1, 2, 2)
        @test_throws ArgumentError gallager_ldpc(10, 3, 4)
        @test_throws ArgumentError gallager_ldpc(8, 1, 4)
        @test_throws ArgumentError gallager_ldpc(8, 3, 1)
        @test_throws ArgumentError gallager_ldpc(4, 3, 6)
        @test_throws ArgumentError gallager_ldpc(4, 4, 2)
        @test_throws ArgumentError gallager_ldpc(2, 2, 2)
        @test_throws ArgumentError gallager_ldpc(4098, 2, 2)
        @test_throws ArgumentError array_ldpc(3, 6, 15)
        @test_throws ArgumentError array_ldpc(3, 6, 5)
        @test_throws ArgumentError array_ldpc(1, 3, 5)
        @test_throws ArgumentError array_ldpc(3, 1, 5)
        @test_throws ArgumentError array_ldpc(2, 2, 1)
        @test_throws ArgumentError array_ldpc(2, 2, 4)
        @test_throws ArgumentError array_ldpc(5, 2, 3)
        @test_throws ArgumentError array_ldpc(2, 2, 4099)
        @test_throws ArgumentError array_ldpc(4, 2, 1031)
        @test_throws ArgumentError ldpc_from_H(zeros(Int, 2, 2))
        @test_throws ArgumentError ldpc_from_H(zeros(Int, 0, 2))
        @test_throws ArgumentError ldpc_from_H(zeros(Int, 2, 0))
        @test_throws ArgumentError ldpc_from_H([1 0; 0 1])
        @test_throws ArgumentError ldpc_from_H(reshape([1, 1], 2, 1))
        @test_throws ArgumentError ldpc_from_H([0 2; 1 1])
        @test_throws ArgumentError ldpc_from_H([0.0 1.0; 1.0 0.0])
        @test_throws ArgumentError ldpc_from_H(ones(Int, 4097, 2))
        @test_throws ArgumentError ldpc_from_H(sparse([1, 1], [1, 1], Int8[1, 1], 2, 2))
        @test_throws ArgumentError ldpc_from_H(sparse([1.0 1.0; 0.0 1.0]))

        spc = ldpc_from_H([1 1 1])
        @test_throws ArgumentError ldpc_encode(spc, [0])
        @test_throws ArgumentError ldpc_encode(spc, [0, 2])
        @test_throws ArgumentError ldpc_syndrome(spc, [0, 0])
        @test_throws ArgumentError ldpc_syndrome(spc, [0, 0, 2])
        @test_throws ArgumentError ldpc_decode_details(spc, [0.0, 0.0])
        @test_throws ArgumentError ldpc_decode_details(spc, [NaN, 0.0, 0.0])
        @test_throws ArgumentError ldpc_decode_details(spc, [Inf, 0.0, 0.0])
        @test_throws ArgumentError ldpc_decode_details(spc, [1.0, 1.0, 1.0]; algorithm=:bp)
        @test_throws ArgumentError ldpc_decode_details(spc, [1.0, 1.0, 1.0];
                                                       algorithm=:normalized_min_sum, alpha=0)
        @test_throws ArgumentError ldpc_decode_details(spc, [1.0, 1.0, 1.0];
                                                       algorithm=:normalized_min_sum, alpha=1.1)
        @test_throws ArgumentError ldpc_decode_details(spc, [1.0, 1.0, 1.0];
                                                       algorithm=:normalized_min_sum, alpha=NaN)
        @test_throws ArgumentError ldpc_decode_details(spc, [1.0, 1.0, 1.0];
                                                       algorithm=:offset_min_sum, beta=-0.1)
        @test_throws ArgumentError ldpc_decode_details(spc, [1.0, 1.0, 1.0];
                                                       algorithm=:offset_min_sum, beta=NaN)
        @test_throws ArgumentError ldpc_decode_details(spc, [1.0, 1.0, 1.0]; max_iter=-1)
        @test_throws ArgumentError ldpc_decode_details(spc, [1, 0, 0]; max_iter=3)
        @test_throws ArgumentError ldpc_decode(spc, [1, 0, 0])
        @test_throws ArgumentError ldpc_decode_details(spc, Rational{Int}[1//2, 1//2, 1//2])
    end

    @testset "Check-node updates and tree MAP" begin
        spc4 = ldpc_from_H([1 1 1 1])
        llr = [0.7, -1.2, 2.5, 0.3]
        cases = (
            (:sum_product, 1.0, 0.0),
            (:min_sum, 1.0, 0.0),
            (:normalized_min_sum, 0.75, 0.0),
            (:normalized_min_sum, 0.8, 0.0),
            (:offset_min_sum, 1.0, 0.15),
            (:offset_min_sum, 1.0, 0.5),
        )
        for (algorithm, alpha, beta) in cases
            got = ldpc_decode_details(spc4, llr; algorithm=algorithm, alpha=alpha,
                                      beta=beta, early_stop=false, max_iter=3)
            kind = algorithm === :min_sum ? :min_sum : algorithm
            ref = _spc_posterior(llr, kind; alpha=alpha, beta=beta)
            @test got.posterior ≈ ref atol=1e-10 rtol=1e-10
            @test got.iterations == 3
            # Bitwise APP hard decisions on a parity check need not form a codeword.
            if got.converged
                @test got.codeword == ldpc_encode(spc4, got.bits)
            end
            @test ldpc_decode(spc4, llr; algorithm=algorithm, alpha=alpha, beta=beta,
                              early_stop=false, max_iter=3) == got.bits
        end
        spa_map = ldpc_decode_details(spc4, llr; algorithm=:sum_product,
                                      early_stop=false, max_iter=4)
        @test spa_map.posterior ≈ _map_llr(spc4, llr) atol=1e-8 rtol=1e-8
        clipped = [0.1, 0.1, 0.1, 0.1]
        got = ldpc_decode_details(spc4, clipped; algorithm=:offset_min_sum, beta=0.5,
                                  early_stop=false, max_iter=2)
        @test got.posterior ≈ _spc_posterior(clipped, :offset_min_sum; beta=0.5) atol=1e-12
        plain = ldpc_decode_details(spc4, llr; algorithm=:min_sum, early_stop=false, max_iter=2)
        offset0 = ldpc_decode_details(spc4, llr; algorithm=:offset_min_sum, beta=0.0,
                                      early_stop=false, max_iter=2)
        @test plain.posterior ≈ offset0.posterior atol=1e-12
        spa = ldpc_decode_details(spc4, llr; algorithm=:sum_product, early_stop=false,
                                  max_iter=2, alpha=0.2, beta=0.9)
        spa_default_scale = ldpc_decode_details(spc4, llr; algorithm=:sum_product,
                                                early_stop=false, max_iter=2)
        @test spa.posterior ≈ spa_default_scale.posterior atol=1e-12
        nms = ldpc_decode_details(spc4, llr; early_stop=false, max_iter=2)
        nms_explicit = ldpc_decode_details(spc4, llr; algorithm=:normalized_min_sum,
                                           alpha=0.75, early_stop=false, max_iter=2)
        @test nms.posterior ≈ nms_explicit.posterior atol=1e-12
        @test !isapprox(nms.posterior, spa.posterior; atol=1e-6)
        @test sprint(show, nms) ==
              "LDPCResult(k=$(length(nms.bits)), iterations=$(nms.iterations), converged=$(nms.converged))"
        @test sprint(show, MIME"text/plain"(), nms) == sprint(show, nms)

        deg1 = ldpc_from_H([1 0 0; 1 1 1])
        @test deg1.k == 1 && deg1.info_pos == [3] && deg1.parity_pos == [1, 2]
        @test ldpc_encode(deg1, [1]) == [0, 1, 1]
        @test ldpc_encode(deg1, [0]) == [0, 0, 0]
        pinned = ldpc_decode_details(deg1, [-20.0, -20.0, -20.0]; algorithm=:sum_product,
                                     early_stop=false, max_iter=4)
        @test isinf(pinned.posterior[1]) && pinned.posterior[1] > 0
        @test pinned.posterior[2] ≈ -40 atol=1e-9
        @test pinned.posterior[3] ≈ -40 atol=1e-9
        @test pinned.posterior ≈ _map_llr(deg1, [-20.0, -20.0, -20.0]) atol=1e-8
        @test pinned.bits == [1] && pinned.codeword == [0, 1, 1] && pinned.converged
        for algorithm in _LDPC_ALGS
            other = ldpc_decode_details(deg1, [-20.0, 20.0, 20.0]; algorithm=algorithm,
                                        max_iter=6)
            @test other.codeword == [0, 0, 0]
            @test other.bits == [0]
            @test other.converged
            forced = ldpc_decode_details(deg1, [-20.0, -20.0, -20.0]; algorithm=algorithm,
                                         early_stop=false, max_iter=4)
            @test forced.codeword == [0, 1, 1]
            @test forced.converged
        end

        tree = ldpc_from_H([1 1 1 0 0; 0 0 1 1 1])
        @test tree.k == 3
        _assert_code_structure(tree)
        @test _assert_encoder_is_kernel(tree) >= 2
        @test !_has_4_cycle(tree)
        Random.seed!(22)
        for _ in 1:4
            channel = randn(tree.n)
            app = ldpc_decode_details(tree, channel; algorithm=:sum_product,
                                      early_stop=false, max_iter=16)
            @test app.posterior ≈ _map_llr(tree, channel) atol=1e-8 rtol=1e-8
            # Bitwise MAP can disagree with the block MAP codeword.
            if app.converged
                @test app.codeword == ldpc_encode(tree, app.bits)
            end
        end
        zeros_llr = zeros(tree.n)
        flat = ldpc_decode_details(tree, zeros_llr; algorithm=:sum_product,
                                   early_stop=false, max_iter=8)
        @test flat.posterior ≈ zeros_llr atol=1e-12
        rep = ldpc_from_H([1 1])
        rep_llr = [1.0, -0.5]
        rep_app = ldpc_decode_details(rep, rep_llr; algorithm=:sum_product,
                                      early_stop=false, max_iter=4)
        @test rep_app.posterior ≈ _map_llr(rep, rep_llr) atol=1e-8 rtol=1e-8

        spc = ldpc_from_H([1 1 1])
        @test spc.info_pos == [2, 3] && spc.parity_pos == [1]
        @test ldpc_encode(spc, [1, 0]) == [1, 1, 0]
        @test ldpc_encode(spc, [0, 1]) == [1, 0, 1]
        @test ldpc_encode(spc, [1, 1]) == [0, 1, 1]
        for channel in ([ -0.4, 2.5, 2.5], [2.0, -0.3, 2.0])
            for algorithm in _LDPC_ALGS
                det = ldpc_decode_details(spc, channel; algorithm=algorithm, max_iter=8)
                @test det.codeword == [0, 0, 0]
                @test det.bits == [0, 0]
                @test det.converged
                @test all(iszero, ldpc_syndrome(spc, det.codeword))
            end
        end
        tie = ldpc_decode_details(spc, [0.0, -1.0, -1.0])
        @test tie.iterations == 0
        @test tie.codeword == [0, 1, 1]
        @test tie.converged
        stopped = ldpc_decode_details(spc, [-1.0, 1.0, 1.0]; max_iter=0)
        @test stopped.iterations == 0
        @test !stopped.converged
        @test stopped.codeword == [1, 0, 0]
        @test !all(iszero, ldpc_syndrome(spc, stopped.codeword))
        already = ldpc_decode_details(spc, [2.0, 3.0, 4.0]; early_stop=true, max_iter=5)
        @test already.iterations == 0 && already.converged
        @test already.posterior == [2.0, 3.0, 4.0]
        @test already.bits == [0, 0]
        cw = ldpc_encode(spc, [1, 0])
        kept = ldpc_decode_details(spc, 12.0 .* (1 .- 2 .* cw); algorithm=:sum_product,
                                   early_stop=false, max_iter=4)
        @test kept.iterations == 4 && kept.converged
        @test kept.bits == [1, 0] && kept.codeword == cw
        float32 = ldpc_decode_details(spc, Float32[4, 4, 4])
        @test float32.bits == [0, 0] && float32.converged
    end

    @testset "Known corrections" begin
        # Lin and Costello §3.2. The Tanner graph has 4-cycles, so this is not ML.
        # Column 7 has weight 3. A single flip there is an undetected error: every
        # algorithm stops on the weight-3 word [0,0,1,0,1,1,0], not the all-zero word.
        h = hamming74_ldpc()
        for algorithm in _LDPC_ALGS
            for j in 1:6
                channel = fill(8.0, 7)
                channel[j] = -8.0
                det = ldpc_decode_details(h, channel; algorithm=algorithm, max_iter=20)
                @test det.converged
                @test det.codeword == zeros(Int, 7)
                @test det.bits == zeros(Int, 4)
            end
            channel = fill(8.0, 7)
            channel[7] = -8.0
            det = ldpc_decode_details(h, channel; algorithm=algorithm, max_iter=20)
            @test det.converged
            @test det.iterations == 1
            @test det.codeword == [0, 0, 1, 0, 1, 1, 0]
            @test all(iszero, ldpc_syndrome(h, det.codeword))
        end

        short = array_ldpc(3, 6, 7)
        for algorithm in _LDPC_ALGS
            for j in 1:short.n
                channel = fill(6.0, short.n)
                channel[j] = -6.0
                det = ldpc_decode_details(short, channel; algorithm=algorithm, max_iter=30)
                @test det.converged
                @test det.codeword == zeros(Int, short.n)
            end
        end
    end

    @testset "Noiseless and high-SNR recovery" begin
        codes = (array_ldpc(3, 6, 11), gallager_ldpc(48, 3, 6))
        Random.seed!(33)
        for code in codes
            _assert_code_structure(code)
            for _ in 1:4
                bits = rand(0:1, code.k)
                word = ldpc_encode(code, bits)
                channel = 30.0 .* (1 .- 2 .* word)
                for algorithm in _LDPC_ALGS
                    det = ldpc_decode_details(code, channel; algorithm=algorithm, max_iter=5)
                    @test det.iterations == 0
                    @test det.converged
                    @test det.bits == bits
                    @test det.codeword == word
                    @test all(iszero, ldpc_syndrome(code, det.codeword))
                end
                held = ldpc_decode_details(code, channel; algorithm=:normalized_min_sum,
                                           early_stop=false, max_iter=2)
                @test held.iterations == 2 && held.converged && held.bits == bits
            end
            bits = rand(0:1, code.k)
            word = ldpc_encode(code, bits)
            y = awgn(modulate(word, BPSK()), Inf; energy_per_bit=code.n / code.k)
            llr = demodulate_soft(y, BPSK(), 1e-6)
            det = ldpc_decode(code, llr; max_iter=5)
            @test det == bits
        end
        σ2 = (66 / 35) / (2 * 10^(2.5 / 10))
        @test demodulate_soft([1.0], BPSK(), σ2)[1] ≈ 2 / σ2 rtol=1e-12
    end

    @testset "AWGN coding gain" begin
        # Fan array code, (dv, dc, p) = (3, 6, 11): n = 66, k = 35, no 4-cycles.
        # Es = 1, Eb = n/k, σ² = Eb / (2 Eb/N0) = N0/2, same call as awgn.
        # 1000 frames, seed RNG_SEED, max_iter = 25. Clopper–Pearson α = 1e-4.
        code = array_ldpc(3, 6, 11)
        nframes = 1000
        max_iter = 25
        Rinv = code.n / code.k
        @test Rinv == 66 / 35
        Random.seed!(RNG_SEED)
        for algorithm in (:normalized_min_sum, :sum_product)
            for snr_db in (2.5, 3.5)
                γ = 10^(snr_db / 10)
                σ2 = Rinv / (2 * γ)
                bit_errors = 0
                block_errors = 0
                for _ in 1:nframes
                    bits = rand(0:1, code.k)
                    word = ldpc_encode(code, bits)
                    y = awgn(modulate(word, BPSK()), snr_db; energy_per_bit=Rinv)
                    hat = ldpc_decode(code, demodulate_soft(y, BPSK(), σ2);
                                      algorithm=algorithm, max_iter=max_iter)
                    errors = count(!=(0), bits .- hat)
                    bit_errors += errors
                    block_errors += errors > 0
                end
                nbits = nframes * code.k
                p_bit = ber_bpsk_awgn(γ)
                p_block = 1 - (1 - p_bit)^code.k
                lo_b, hi_b = clopper_pearson(bit_errors, nbits, BER_ALPHA)
                lo_B, hi_B = clopper_pearson(block_errors, nframes, BER_ALPHA)
                println("LDPC array(3,6,11) $(algorithm) Eb/N0=$(snr_db) dB " *
                        "bits=$(bit_errors)/$(nbits) BER_CP=($(lo_b), $(hi_b)) " *
                        "uncoded_BER=$(p_bit) blocks=$(block_errors)/$(nframes) " *
                        "BLER_CP=($(lo_B), $(hi_B)) uncoded_BLER=$(p_block)")
                @test bit_errors >= 15
                @test block_errors >= 4
                @test hi_b < p_bit
                @test hi_B < p_block
            end
        end
    end
end
