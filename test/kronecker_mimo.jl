# Flat Kronecker MIMO: Loyka exponential correlation, separable moments,
# identity / ρ = 0 collapsing to mimo_channel, and a ZF BER gap.
# Entry moments are the product in the kronecker_channel docstring
# (Kermoal et al., IEEE JSAC 20(6), 2002, Appendix). No new BER closed form.
# Clopper–Pearson uses the same α as reference.jl.

@testset "Exponential correlation" begin
    R = exponential_correlation(3, 0.5)
    @test R ≈ ComplexF64[1 0.5 0.25; 0.5 1 0.5; 0.25 0.5 1]
    @test R ≈ R'

    ρ = -0.4
    Rn = exponential_correlation(3, ρ)
    @test Rn ≈ ComplexF64[1 ρ ρ^2; ρ 1 ρ; ρ^2 ρ 1]
    @test all(iszero, imag.(Rn))

    ρc = 0.5 * cis(0.4)
    Rc = exponential_correlation(3, ρc)
    @test Rc[1, 1] == 1
    @test Rc[1, 2] ≈ ρc
    @test Rc[1, 3] ≈ ρc^2
    @test Rc[2, 1] ≈ conj(ρc)
    @test Rc[3, 1] ≈ conj(ρc)^2
    @test Rc ≈ Rc'
    # |ρ| = 0.5 < 1, so the Toeplitz factor is positive definite.
    @test minimum(eigvals(Hermitian(Rc))) > 0

    # Adjacent spacing 2: lag-1 coefficient is ρ^2, lag-2 is ρ^4.
    Rs = exponential_correlation(3, 0.5; spacing=2)
    @test Rs ≈ ComplexF64[1 0.25 0.0625; 0.25 1 0.25; 0.0625 0.25 1]

    # Negative base, non-integer lag: arg(ρ) = π, so the entry is imaginary.
    Rfrac = exponential_correlation(2, -0.5; spacing=0.5)
    @test Rfrac[1, 2] ≈ sqrt(0.5) * im atol = 1e-12
    @test Rfrac[2, 1] ≈ -sqrt(0.5) * im atol = 1e-12
    @test Rfrac ≈ Rfrac'

    @test exponential_correlation(4, 0) == Matrix{ComplexF64}(I, 4, 4)
    @test exponential_correlation(2, 0; spacing=0) == Matrix{ComplexF64}(I, 2, 2)
    @test exponential_correlation(2, 0.7; spacing=0) ≈ ones(ComplexF64, 2, 2)
    @test only(exponential_correlation(1, -0.3)) == 1
    @test exponential_correlation(3, 1) ≈ ones(ComplexF64, 3, 3)
    @test exponential_correlation(2, 1 + 1e-15) ≈ ones(ComplexF64, 2, 2)

    # |ρ| = 1 is a rank-1 phase ramp, still a correlation matrix.
    φ = 0.4
    Ru = exponential_correlation(3, cis(φ))
    @test rank(Ru, 1e-8) == 1
    @test Ru ≈ Ru'
    @test all(abs.(Ru) .≈ 1)

    Random.seed!(1)
    mark = randn()
    Random.seed!(1)
    @test exponential_correlation(5, 0.2) ≈ exponential_correlation(5, 0.2)
    @test randn() == mark

    @test_throws ArgumentError exponential_correlation(0, 0.2)
    @test_throws ArgumentError exponential_correlation(-1, 0.2)
    @test_throws ArgumentError exponential_correlation(2, 1.01)
    @test_throws ArgumentError exponential_correlation(2, NaN)
    @test_throws ArgumentError exponential_correlation(2, Inf)
    @test_throws ArgumentError exponential_correlation(2, 0.2; spacing=-0.1)
    @test_throws ArgumentError exponential_correlation(2, 0.2; spacing=Inf)
    @test_throws ArgumentError exponential_correlation(2, 0.2; spacing=NaN)
end

"""Hermitian PSD square root, written out so the test does not call the private helper."""
function _spec_sqrt(R::AbstractMatrix)
    F = eigen(Hermitian(Matrix{ComplexF64}(R)))
    return F.vectors * Diagonal(sqrt.(max.(F.values, 0.0))) * F.vectors'
end

@testset "Kronecker special cases" begin
    # Integer call and exact identities are mimo_channel, including the randn order.
    Random.seed!(RNG_SEED)
    A = kronecker_channel(3, 2)
    Random.seed!(RNG_SEED)
    B = mimo_channel(3, 2)
    @test A == B
    Random.seed!(RNG_SEED)
    C = kronecker_channel(3, 2, Matrix{Float64}(I, 3, 3), Matrix{Float64}(I, 2, 2))
    @test C == B
    Random.seed!(RNG_SEED)
    D = kronecker_channel(Matrix{ComplexF64}(I, 3, 3), Matrix{ComplexF64}(I, 2, 2))
    @test D == B

    # ρ = 0 is I, so the same stream.
    Random.seed!(RNG_SEED)
    E = kronecker_channel(exponential_correlation(3, 0), exponential_correlation(2, 0))
    @test E == B

    # Non-trivial factors use that same H_w, then the Hermitian square roots.
    Rt = exponential_correlation(2, 0.45)
    Rr = exponential_correlation(3, 0.2 + 0.1im)
    Random.seed!(11)
    H = kronecker_channel(Rt, Rr)
    Random.seed!(11)
    Hw = (randn(3, 2) .+ im .* randn(3, 2)) ./ sqrt(2)
    @test H ≈ _spec_sqrt(Rr) * Hw * _spec_sqrt(Rt) atol = 1e-10
    @test size(H) == (3, 2)

    # ρ = 1 at the transmitter: both columns are the same vector.
    Random.seed!(12)
    Hc = kronecker_channel(exponential_correlation(2, 1), Matrix{Float64}(I, 3, 3))
    @test Hc[:, 1] ≈ Hc[:, 2] atol = 1e-8
    # ρ = 1 at the receiver: both rows are the same vector. Transmit factor is exact I.
    Random.seed!(13)
    Hr = kronecker_channel(Matrix{Float64}(I, 2, 2), exponential_correlation(2, 1))
    @test Hr[1, :] ≈ Hr[2, :] atol = 1e-8

    # A rounding-negative eigenvalue is clipped. The killed receive mode stays at 0.
    Rclip = ComplexF64[1 0; 0 -1e-12]
    Random.seed!(14)
    Hk = kronecker_channel(Matrix{Float64}(I, 2, 2), Rclip)
    @test maximum(abs, Hk[2, :]) < 1e-6

    # Hermitian part within 1e-8: accepted, and the square root reproduces R.
    Rnear = ComplexF64[1.0 0.4 + 0.2im; 0.4 - 0.2im + 1e-10 1.0]
    S = _spec_sqrt((Rnear .+ Rnear') ./ 2)
    @test S * S ≈ (Rnear .+ Rnear') ./ 2 atol = 1e-8
    Random.seed!(15)
    @test all(isfinite, kronecker_channel(Matrix{Float64}(I, 2, 2), Rnear))

    # 1×1 covariances: E[|H|^2] = R_rx[1,1] * R_tx[1,1] = 6.
    Random.seed!(16)
    p = 0.0
    n1 = 20_000
    for _ in 1:n1
        p += abs2(only(kronecker_channel(fill(3.0, 1, 1), fill(2.0, 1, 1))))
    end
    @test p / n1 ≈ 6 rtol = 0.05

    bad = [1.0 2.0; 0.0 1.0]
    Random.seed!(1)
    mark = randn()
    Random.seed!(1)
    @test_throws ArgumentError kronecker_channel(bad, Matrix{Float64}(I, 2, 2))
    @test randn() == mark
    @test_throws ArgumentError kronecker_channel([1.0 0.9; 0.9 0.5], Matrix{Float64}(I, 2, 2))
    @test_throws ArgumentError kronecker_channel(ComplexF64[1 NaN; NaN 1], Matrix{Float64}(I, 2, 2))
    @test_throws ArgumentError kronecker_channel(ComplexF64[1 Inf; Inf 1], Matrix{Float64}(I, 2, 2))
    @test_throws ArgumentError kronecker_channel([1.0 0.0; 0.0 -1e-6], Matrix{Float64}(I, 2, 2))
    @test_throws ArgumentError kronecker_channel([1.0 0.0 0.0; 0.0 1.0 0.0], Matrix{Float64}(I, 2, 2))
    @test_throws ArgumentError kronecker_channel(zeros(0, 0), zeros(0, 0))
    @test_throws ArgumentError kronecker_channel(0, 1)
    @test_throws ArgumentError kronecker_channel(2, 2, Matrix{Float64}(I, 3, 3), Matrix{Float64}(I, 2, 2))
    @test_throws ArgumentError kronecker_channel(2, 2, Matrix{Float64}(I, 2, 2), Matrix{Float64}(I, 3, 3))
    @test_throws ArgumentError kronecker_channel(0, 2, ones(1, 1), ones(1, 1))
end

@testset "Kronecker correlation moments" begin
    # Complex receive factor, real transmit factor. A conjugated R_r, or a
    # swapped (R_t)[q,s], misses the imaginary part by far more than the tolerance.
    ρr = 0.4 * cis(0.7)
    Rt = exponential_correlation(2, 0.55)
    Rr = exponential_correlation(2, ρr)
    nr, nt = 2, 2
    N = 20_000
    acc = zeros(ComplexF64, nr, nt, nr, nt)
    acc_rx = zeros(ComplexF64, nr, nr)
    acc_tx = zeros(ComplexF64, nt, nt)
    Random.seed!(RNG_SEED)
    for _ in 1:N
        H = kronecker_channel(Rt, Rr)
        acc_rx .+= H * H'
        acc_tx .+= H' * H
        for s in 1:nt, r in 1:nr, q in 1:nt, p in 1:nr
            acc[p, q, r, s] += H[p, q] * conj(H[r, s])
        end
    end
    acc ./= N
    acc_rx ./= N
    acc_tx ./= N
    for s in 1:nt, r in 1:nr, q in 1:nt, p in 1:nr
        @test acc[p, q, r, s] ≈ Rr[p, r] * Rt[s, q] atol = 0.06
    end
    @test acc_rx ≈ tr(Rt) * Rr atol = 0.06
    @test acc_tx ≈ tr(Rr) * Rt atol = 0.06
    # The imaginary part of E[H_11 conj(H_21)] is the imaginary part of ρr, not of conj(ρr).
    @test imag(acc[1, 1, 2, 1]) > 0.1
    @test imag(ρr) > 0.1
end

@testset "Kronecker link and detectors" begin
    x = ComplexF64[1.0, -1.0]
    Rt = exponential_correlation(2, 0.3)
    Rr = exponential_correlation(2, 0.5)
    Random.seed!(4)
    y0, H0 = apply_kronecker(x, Rt, Rr, 0.0)
    @test y0 ≈ H0 * x
    Random.seed!(4)
    y1, H1 = apply_kronecker(x, Rt, Rr, 0.0)
    @test y1 == y0
    @test H1 == H0

    Random.seed!(5)
    yn, Hn = apply_kronecker(x, Rt, Rr, 0.25)
    Random.seed!(5)
    H2 = kronecker_channel(Rt, Rr)
    y2 = apply_mimo(x, H2, 0.25)
    @test Hn == H2
    @test yn == y2

    Random.seed!(1)
    mark = randn()
    Random.seed!(1)
    @test_throws ArgumentError apply_kronecker(x, Rt, Rr, -0.1)
    @test randn() == mark
    @test_throws ArgumentError apply_kronecker([1.0], Rt, Rr, 0.0)
    @test_throws ArgumentError apply_kronecker(x, [1.0 0.2 0.0; 0.2 1.0 0.0; 0.0 0.0 1.0], Rr, 0.0)

    bpsk = ComplexF64[1.0, -1.0]
    Random.seed!(6)
    H = kronecker_channel(exponential_correlation(2, 0.35), exponential_correlation(2, 0.25))
    @test cond(H) < 30
    y = apply_mimo(x, H, 0.0)
    @test ml_detect(y, H, bpsk) ≈ x
    @test zf_detect(y, H) ≈ x atol = 1e-8
    @test mmse_detect(y, H, 0.0) ≈ x atol = 1e-8
    @test sic_detect(y, H, bpsk) ≈ x

    Random.seed!(8)
    Ht = kronecker_channel(exponential_correlation(2, 0.4), exponential_correlation(3, 0.5))
    @test size(Ht) == (3, 2)
    yt = apply_mimo(x, Ht, 0.0)
    @test zf_detect(yt, Ht) ≈ x atol = 1e-8
    @test mmse_detect(yt, Ht, 0.0) ≈ x atol = 1e-8
    @test ml_detect(yt, Ht, bpsk) ≈ x

    # Same complex noise variance. Strong Kronecker correlation is ill-conditioned,
    # so ZF is worse than the i.i.d. draw. Intervals must separate (α = 1e-4).
    function kronecker_zf_errors(ρ::Float64, noise_var::Float64, n_vec::Int)
        R = exponential_correlation(2, ρ)
        k = 0
        for _ in 1:n_vec
            bits = rand(0:1, 2)
            s = ComplexF64.(modulate(bits, BPSK()))
            yb, Hb = apply_kronecker(s, R, R, noise_var)
            z = zf_detect(yb, Hb)
            k += count(!=(0), bits .- demodulate(real.(z), BPSK()))
        end
        return k, 2 * n_vec
    end
    Random.seed!(RNG_SEED)
    # noise_var = 0.05 is 13 dB below unit-energy BPSK on an identity channel.
    # ρ = 0.95 makes the 2×2 correlation condition number (1+ρ)/(1-ρ) = 39.
    k_iid, n_iid = kronecker_zf_errors(0.0, 0.05, 4_000)
    k_cor, n_cor = kronecker_zf_errors(0.95, 0.05, 4_000)
    _, hi_iid = clopper_pearson(k_iid, n_iid, BER_ALPHA)
    lo_cor, hi_cor = clopper_pearson(k_cor, n_cor, BER_ALPHA)
    @test k_iid >= 20
    @test k_cor >= 20
    @test hi_iid < lo_cor
    @test hi_cor > hi_iid

    # On a fixed moderate correlation, raising the noise raises the ZF error rate.
    Random.seed!(RNG_SEED)
    k_hi, n_hi = kronecker_zf_errors(0.5, 1.0, 2_000)
    k_lo, n_lo = kronecker_zf_errors(0.5, 1e-3, 2_000)
    lo_hi, _ = clopper_pearson(k_hi, n_hi, BER_ALPHA)
    _, hi_lo = clopper_pearson(k_lo, n_lo, BER_ALPHA)
    @test k_hi >= 20
    @test hi_lo < lo_hi
end
