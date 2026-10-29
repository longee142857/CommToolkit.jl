# EXIT chart utilities (1.1.0). Deterministic throughout: the transfer curves
# use the module's seeded xorshift64* + Box–Muller stream, and the estimator
# checks below replicate that stream locally, as test/bicm.jl does, because
# Julia's global `randn` is not stable across 1.x versions.
#
# References:
#   ten Brink, IEEE Trans. Commun. 49(10):1727–1737, 2001 (EXIT chart,
#     histogram MI measurement).
#   Hagenauer, Proc. 12th EUSIPCO, Vienna, 2004, pp. 1541–1548 (histogram-free
#     soft MI estimator).
#   ten Brink, Speidel, Yan, IEEE Trans. Commun. 52(4):670–678, 2004
#     (J(σ) approximation and its inverse, eqs. (9)–(11)).
#   Hagenauer, Offer, Papke, IEEE Trans. Inf. Theory 42(2):429–445, 1996
#     (Le = Lpost − Lprior).
#   BIAWGN capacity: numerical integral in test/reference.jl (biawgn_llr_mi),
#   no elementary closed form.

function _exit_xorshift(state::UInt64)
    z = state
    z = xor(z, z >> 12)
    z = xor(z, z << 25)
    z = xor(z, z >> 27)
    return z
end

function _exit_u01(st::Ref{UInt64})
    st[] = _exit_xorshift(st[])
    return ((st[] * 0x2545f4914f6cdd1d) >> 11) * 2.0^-53
end

function _exit_normal(st::Ref{UInt64})
    u1 = max(_exit_u01(st), 1e-300)
    u2 = _exit_u01(st)
    return sqrt(-2 * log(u1)) * cos(2π * u2)
end

# BIAWGN LLR stream in the package convention: bit 0 -> +1, L = 2y/σ².
function _exit_biawgn(σ2::Float64, N::Int, seed::UInt64)
    st = Ref(seed == 0 ? 0x9e3779b97f4a7c15 : seed)
    bits = Vector{Int}(undef, N)
    L = Vector{Float64}(undef, N)
    σ = sqrt(σ2)
    for i in 1:N
        b = _exit_u01(st) >= 0.5 ? 1 : 0
        bits[i] = b
        L[i] = 2 * ((1 - 2b) + σ * _exit_normal(st)) / σ2
    end
    return bits, L
end

@testset "EXIT estimators vs BIAWGN capacity" begin
    # Known case: consistent BIAWGN LLRs. True MI is the BIAWGN capacity,
    # evaluated by quadrature in test/reference.jl. N = 200_000 keeps the
    # Monte Carlo standard error near 1e-3; the observed deviations at seed
    # 0x0c0de1162 are ≤ 0.0014 for :soft and ≤ 0.0008 for :hist.
    N = 200_000
    for db in (-2.0, 0.0, 2.0, 4.0)
        σ2 = 1 / (2 * 10^(db / 10))
        C = biawgn_llr_mi(σ2)
        bits, L = _exit_biawgn(σ2, N, 0x0c0de1162)
        @test abs(exit_mi_soft(L, bits) - C) < 0.01
        @test abs(exit_mi_hist(L, bits) - C) < 0.01
        # J(σ) with σ_L = 2/σ is the same integral; the polynomial
        # approximation deviates ≤ 0.0004 at these SNRs.
        @test abs(exit_j(2 / sqrt(σ2)) - C) < 0.005
    end
    # Degenerate and infinite LLRs.
    @test exit_mi_soft(zeros(4), [0, 1, 0, 1]) == 0.0
    @test exit_mi_hist(zeros(4), [0, 1, 0, 1]) == 0.0
    @test exit_mi_soft([1e9, -1e9], [0, 1]) == 1.0
    @test exit_mi_soft([Inf, -Inf], [0, 1]) == 1.0
    @test exit_mi_soft([-Inf], [0]) == -Inf
    # Single-class input is legal for :soft (all-zero-codeword practice) and
    # rejected for :hist (two conditional densities are needed).
    @test exit_mi_soft([3.0, 2.0], [0, 0]) > 0
    @test_throws ArgumentError exit_mi_hist([3.0, 2.0], [0, 0])
    # Soft and histogram agree on a non-Gaussian consistent stream: BPSK under
    # flat Rayleigh with ideal CSI (LLR = 2 Re(h* y)/σ², not Gaussian).
    st = Ref(UInt64(0x0c0de1162))
    σ2 = 0.5
    Nb = 100_000
    bits = Vector{Int}(undef, Nb)
    L = Vector{Float64}(undef, Nb)
    for i in 1:Nb
        b = _exit_u01(st) >= 0.5 ? 1 : 0
        h = (_exit_normal(st) + im * _exit_normal(st)) / sqrt(2)
        y = h * (1 - 2b) + sqrt(σ2) * (_exit_normal(st) + im * _exit_normal(st))
        bits[i] = b
        L[i] = 2 * real(conj(h) * y) / σ2
    end
    @test abs(exit_mi_soft(L, bits) - exit_mi_hist(L, bits)) < 0.02
end

@testset "EXIT J(sigma) approximation" begin
    @test exit_j(0.0) == 0.0
    @test exit_j(Inf) == 1.0
    @test exit_jinv(0.0) == 0.0
    @test exit_jinv(1.0) == Inf
    @test_throws ArgumentError exit_j(-1.0)
    @test_throws ArgumentError exit_jinv(-0.1)
    @test_throws ArgumentError exit_jinv(1.1)
    grid = 0.0:0.25:6.0
    @test all(diff(exit_j.(grid)) .> 0)
    igrid = 0.01:0.05:0.99
    @test all(diff(exit_jinv.(igrid)) .> 0)
    # Round trips, ten Brink 2004 eqs. (9)–(11): observed |j(jinv(I)) − I| ≤
    # 0.003 and |jinv(j(σ)) − σ| ≤ 0.022 on these grids.
    for I in (0.05, 0.2, 0.3646, 0.5, 0.8, 0.95, 0.999)
        @test abs(exit_j(exit_jinv(I)) - I) < 0.01
    end
    for σ in (0.5, 1.0, 2.0, 3.0, 5.0)
        @test abs(exit_jinv(exit_j(σ)) - σ) < 0.03
    end
end

@testset "EXIT demapper transfer" begin
    # BPSK at Es/N0 = 0 dB: extrinsic = channel LLR, flat in I_A, equal to the
    # BIAWGN capacity at σ² = 1/2.
    σ2 = 0.5
    c = exit_demap_transfer(BPSK(); ia=0.0:0.3:0.9, esn0_db=0.0, frames=8, n_symbols=4000)
    C = biawgn_llr_mi(σ2)
    @test abs(c.ie[1] - C) < 0.02
    @test all(abs.(c.ie .- C) .< 0.03)
    # Gray QPSK: two independent axes of energy 1/2, so the anchor is the
    # capacity at 2σ²; flat in I_A as well.
    q = exit_demap_transfer(QPSK(); ia=[0.0, 0.9], esn0_db=0.0, frames=8, n_symbols=4000)
    @test abs(q.ie[1] - biawgn_llr_mi(2σ2)) < 0.02
    @test abs(q.ie[2] - q.ie[1]) < 0.03
    # 8-PSK at Es/N0 = 6 dB: set partitioning starts below Gray and gains
    # strongly from a priori; Gray is nearly flat. Observed gaps: SP +0.32,
    # Gray +0.013, Gray − SP at I_A = 0 is 0.18.
    g = exit_demap_transfer(PSK(8); ia=[0.0, 0.98], esn0_db=6.0, frames=4, n_symbols=4000)
    s = exit_demap_transfer(PSK(8); ia=[0.0, 0.98], esn0_db=6.0, frames=4, n_symbols=4000,
                            labels=bicm_psk8_id_labels())
    @test s.ie[2] - s.ie[1] > 0.2
    @test g.ie[2] - g.ie[1] < 0.1
    @test s.ie[1] < g.ie[1] - 0.1
    # Smoke: 16-QAM and max-log run and stay in [0, 1].
    m = exit_demap_transfer(QAM(16); ia=[0.0, 0.5], esn0_db=6.0, frames=2,
                            n_symbols=2000, approx=:maxlog)
    @test all(0 .<= m.ie .<= 1)
    # Histogram method through the same entry point, N >> nbins.
    h = exit_demap_transfer(BPSK(); ia=[0.0], esn0_db=0.0, frames=25,
                            n_symbols=8000, method=:hist, nbins=257)
    @test abs(h.ie[1] - C) < 0.02
end

@testset "EXIT decoder transfer" begin
    # Terminated (7, 5) BCJR. High SNR: the channel alone decides.
    hi = exit_decoder_transfer(; decoder=:bcjr, ia=[0.0], ebn0_db=10.0,
                               frames=4, n_info=200)
    @test hi.ie[1] > 0.99
    # 0 dB: T(0) interior, curve rises with I_A. Observed 0.729 and 0.957 at
    # these parameters.
    lo = exit_decoder_transfer(; decoder=:bcjr, ia=[0.0, 0.99], ebn0_db=0.0,
                               frames=16, n_info=200)
    @test 0.5 < lo.ie[1] < 0.9
    @test lo.ie[2] > lo.ie[1] + 0.1
    # Curve at I_A = 0 rises with SNR: observed 0.983 at 3 dB vs 0.729 at 0 dB.
    mid = exit_decoder_transfer(; decoder=:bcjr, ia=[0.0], ebn0_db=3.0,
                                frames=16, n_info=200)
    @test mid.ie[1] > lo.ie[1] + 0.15
    @test all(0 .<= lo.ie .<= 1)
    # Short LDPC smoke: array (3, 6, 11), normalized min-sum, 1 dB. Observed
    # T(0) = 0.84 and T(0.9) = 1.0.
    code = array_ldpc(3, 6, 11)
    lc = exit_decoder_transfer(; decoder=:ldpc, code=code, ia=[0.0, 0.9],
                               ebn0_db=1.0, frames=8, ldpc_iterations=20)
    @test 0.3 < lc.ie[1] < 1.0
    @test lc.ie[2] > lc.ie[1] + 0.05
    @test all(0 .<= lc.ie .<= 1)
    # Consistency teaching point, pinned to the fixed seed: flooding
    # sum-product on a cyclic graph at 2.5 dB converges to a wrong state on at
    # least one of these 16 frames and returns overconfident wrong LLRs. The
    # soft estimator, which assumes consistent LLRs, punishes that; the
    # histogram estimator still sees the information. Observed: hist = 0.96,
    # soft = -120.
    sp_soft = exit_decoder_transfer(; decoder=:ldpc, code=code, ia=[0.0],
                                    ebn0_db=2.5, frames=16, algorithm=:sum_product)
    sp_hist = exit_decoder_transfer(; decoder=:ldpc, code=code, ia=[0.0],
                                    ebn0_db=2.5, frames=16, algorithm=:sum_product,
                                    method=:hist, nbins=101)
    @test sp_hist.ie[1] > 0.9
    @test sp_hist.ie[1] - sp_soft.ie[1] > 0.2
end

@testset "EXIT container and argument checks" begin
    c = EXITCurve([0.0, 0.5, 1.0], [0.2, 0.6, 0.9])
    @test c.ia == [0.0, 0.5, 1.0]
    @test c.ie == [0.2, 0.6, 0.9]
    @test occursin("EXITCurve", sprint(show, c))
    @test_throws ArgumentError EXITCurve([0.1], [0.1, 0.2])
    @test_throws ArgumentError EXITCurve(Float64[], Float64[])
    @test_throws ArgumentError EXITCurve([NaN], [0.1])
    @test_throws ArgumentError exit_mi_soft([0.1], [0, 1])
    @test_throws ArgumentError exit_mi_soft([0.1], [2])
    @test_throws ArgumentError exit_mi_soft([NaN], [0])
    @test_throws ArgumentError exit_mi_hist([1.0], [0]; nbins=2)
    @test_throws ArgumentError exit_mi_hist([1.0], [0]; clip=0.0)
    @test_throws ArgumentError exit_demap_transfer(DBPSK())
    @test_throws ArgumentError exit_demap_transfer(BPSK(); ia=[1.2])
    @test_throws ArgumentError exit_demap_transfer(BPSK(); esn0_db=Inf)
    @test_throws ArgumentError exit_demap_transfer(BPSK(); method=:bogus)
    @test_throws ArgumentError exit_demap_transfer(BPSK(); frames=0)
    @test_throws ArgumentError exit_demap_transfer(BPSK(); sigma_max=0.0)
    @test_throws ArgumentError exit_decoder_transfer(; decoder=:turbo)
    @test_throws ArgumentError exit_decoder_transfer(; decoder=:ldpc)
    @test_throws ArgumentError exit_decoder_transfer(; decoder=:bcjr, code=array_ldpc(3, 6, 11))
    @test_throws ArgumentError exit_decoder_transfer(; decoder=:ldpc,
                                                     code=array_ldpc(3, 6, 11), n_info=30)
    @test_throws ArgumentError exit_decoder_transfer(; decoder=:bcjr, ebn0_db=-Inf)
end
