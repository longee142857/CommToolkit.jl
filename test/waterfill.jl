# SVD water-filling. Powers meet the total-power constraint. High SNR is
# nearly uniform across positive modes; low SNR sits on the strongest mode.
# Tse and Viswanath, 2005, §7.1.1: p_i = (μ - N0 / λ_i^2)_+.

function _kkt(σ, power, noise, p, level; atol=1e-8)
    @test all(isfinite, p)
    @test all(p .>= -atol)
    @test sum(p) ≈ power atol = atol
    for (s, pi) in zip(σ, p)
        if s == 0
            @test pi ≈ 0 atol = atol
            continue
        end
        ν = noise / abs2(float(s))
        if pi > 1e-8
            @test pi + ν ≈ level atol = atol
        else
            @test ν >= level - atol
        end
    end
end

@testset "SVD water-filling" begin
    # Two positive modes, both above water. Dyadic, so the powers are exact.
    σ2 = [4.0, 1.0]
    p2, μ2 = waterfill(σ2, 1.0, 1.0)
    @test p2 ≈ [0.96875, 0.03125] atol = 1e-12
    @test μ2 ≈ 1.03125 atol = 1e-12
    _kkt(σ2, 1.0, 1.0, p2, μ2)

    # Same channel, weakest mode shut off.
    p2lo, μ2lo = waterfill(σ2, 1.0, 100.0)
    @test p2lo ≈ [1.0, 0.0] atol = 1e-12
    @test μ2lo ≈ 7.25 atol = 1e-12
    _kkt(σ2, 1.0, 100.0, p2lo, μ2lo)

    # Three modes. P = 1, N0 = 1 turns the weakest off and keeps the other two.
    # λ = [4, 1, 0.25], ν = [0.25, 1, 4], μ = 1.125.
    σ3 = [2.0, 1.0, 0.5]
    p3, μ3 = waterfill(σ3, 1.0, 1.0)
    @test p3 ≈ [0.875, 0.125, 0.0] atol = 1e-12
    @test μ3 ≈ 1.125 atol = 1e-12
    _kkt(σ3, 1.0, 1.0, p3, μ3)
    @test issorted(p3; rev=true)

    p3hi, _ = waterfill(σ3, 3.0, 1e-12)
    @test maximum(abs.(p3hi .- 1.0)) < 1e-9
    _kkt(σ3, 3.0, 1e-12, p3hi, waterfill(σ3, 3.0, 1e-12)[2]; atol=1e-8)

    p3lo, _ = waterfill(σ3, 1.0, 10.0)
    @test p3lo[1] ≈ 1.0 atol = 1e-12
    @test p3lo[2] ≈ 0.0 atol = 1e-12
    @test p3lo[3] ≈ 0.0 atol = 1e-12
    # High SNR is the flatter allocation on this same singular-value vector.
    @test maximum(p3lo) - minimum(p3lo) > maximum(p3hi) - minimum(p3hi)

    # Order of the input is the caller's order, not strongest-first.
    p_unsort, μ_unsort = waterfill([1.0, 0.0, 4.0], 1.0, 1.0)
    @test p_unsort ≈ [0.03125, 0.0, 0.96875] atol = 1e-12
    @test μ_unsort ≈ μ2 atol = 1e-12

    # N0 = 0 shares power across every strictly positive mode.
    p0, μ0 = waterfill([2.0, 1.0, 0.0], 3.0, 0.0)
    @test p0 ≈ [1.5, 1.5, 0.0] atol = 1e-12
    @test μ0 ≈ 1.5 atol = 1e-12

    pzero, μzero = waterfill([2.0, 1.0], 0.0, 0.4)
    @test pzero == [0.0, 0.0]
    @test μzero ≈ 0.4 / 4 atol = 1e-12

    pnil, μnil = waterfill([0.0, 0.0], 0.0, 1.0)
    @test pnil == [0.0, 0.0]
    @test μnil == 0.0

    # One mode: the whole budget, level = P + N0/λ^2.
    p1, μ1 = waterfill([2.0], 1.5, 0.3)
    @test p1 ≈ [1.5] atol = 1e-12
    @test μ1 ≈ 1.5 + 0.3 / 4 atol = 1e-12

    # Float32 singular values, and a row-vector view, take the same path.
    @test waterfill(Float32[4, 1], 1, 1)[1] ≈ p2 atol = 1e-6
    @test waterfill(view([1.0, 0.0, 4.0], :), 1.0, 1.0)[1] ≈ p_unsort atol = 1e-12

    @test_throws ArgumentError waterfill(Float64[], 1.0, 1.0)
    @test_throws ArgumentError waterfill([1.0], -0.1, 1.0)
    @test_throws ArgumentError waterfill([1.0], Inf, 1.0)
    @test_throws ArgumentError waterfill([1.0], 1.0, -0.2)
    @test_throws ArgumentError waterfill([1.0], 1.0, NaN)
    @test_throws ArgumentError waterfill([1.0, -0.5], 1.0, 1.0)
    @test_throws ArgumentError waterfill([1.0, NaN], 1.0, 1.0)
    @test_throws ArgumentError waterfill([0.0, 0.0], 1.0, 1.0)

    H = ComplexF64[4 0 0; 0 1 0]
    F, p, level = svd_waterfill(H, 1.0, 1.0)
    V, S = svd_precoder(H)
    @test S ≈ [4.0, 1.0] atol = 1e-12
    @test p ≈ p2 atol = 1e-12
    @test level ≈ μ2 atol = 1e-12
    @test F ≈ V * Diagonal(sqrt.(p)) atol = 1e-10
    @test sum(abs2, F) ≈ 1.0 atol = 1e-10
    @test size(F) == (3, 2)

    Flo, plo, _ = svd_waterfill(H, 1.0, 100.0)
    @test plo ≈ [1.0, 0.0] atol = 1e-10
    @test sum(abs2, Flo[:, 1]) ≈ 1.0 atol = 1e-10
    @test norm(Flo[:, 2]) ≈ 0.0 atol = 1e-10

    Fhi, phi, _ = svd_waterfill(H, 2.0, 1e-12)
    @test maximum(abs.(phi .- 1.0)) < 1e-9
    @test sum(abs2, Fhi) ≈ 2.0 atol = 1e-8

    # Rank-1 rectangular channel: the zero mode stays a zero column.
    Hr = ComplexF64[1 0 0; 0 0 0]
    Fr, pr, _ = svd_waterfill(Hr, 2.0, 0.5)
    @test pr[1] ≈ 2.0 atol = 1e-10
    @test pr[2] ≈ 0.0 atol = 1e-10
    @test sum(abs2, Fr) ≈ 2.0 atol = 1e-10
    @test norm(Fr[:, 2]) ≈ 0.0 atol = 1e-10

    Z, pz, lz = svd_waterfill(zeros(ComplexF64, 2, 2), 0.0, 1.0)
    @test pz == [0.0, 0.0]
    @test lz == 0.0
    @test Z == zeros(ComplexF64, 2, 2)

    @test_throws ArgumentError svd_waterfill(zeros(ComplexF64, 2, 2), 1.0, 1.0)
    @test_throws ArgumentError svd_waterfill(zeros(ComplexF64, 0, 2), 1.0, 1.0)
    @test_throws ArgumentError svd_waterfill(zeros(ComplexF64, 2, 0), 1.0, 1.0)
    @test_throws ArgumentError svd_waterfill(H, -1.0, 1.0)

    Random.seed!(RNG_SEED)
    for (nt, nr, P, N0) in ((2, 2, 1.0, 0.3), (4, 2, 2.5, 0.2), (3, 5, 0.7, 1e-4), (3, 3, 4.0, 10.0))
        Hrand = mimo_channel(nt, nr)
        F, p, level = svd_waterfill(Hrand, P, N0)
        V, S = svd_precoder(Hrand)
        @test F ≈ V * Diagonal(sqrt.(p)) atol = 1e-8
        @test sum(abs2, F) ≈ P atol = 1e-8
        _kkt(S, P, N0, p, level)
        @test issorted(p; rev=true)
    end

    # svd_precoder itself still returns unweighted right singular vectors.
    Hs = ComplexF64[1.2+0.1im 0.3-0.2im; 0.15+0.05im -1.0+0.4im]
    V0, S0 = svd_precoder(Hs)
    @test V0 ≈ svd(Hs).V
    @test S0 ≈ svd(Hs).S
end
