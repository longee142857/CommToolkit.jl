# Gray M-PSK, DBPSK/DQPSK, flat Rice, and ideal-CSI MRC.
# Closed forms live in reference.jl. Monte Carlo points use the same
# Clopper–Pearson rule as the AWGN and Rayleigh tests (α = 1e-4).

@testset "M-PSK, differential PSK, Rice, MRC formulas" begin
    for db in (0.0, 4.0, 8.0, 12.0)
        γ = 10^(db / 10)
        @test ber_mpsk_gray(2, γ) ≈ ber_bpsk_awgn(γ) rtol = 1e-5 atol = 1e-10
        @test ber_mpsk_gray(4, γ) ≈ ber_bpsk_awgn(γ) rtol = 1e-5 atol = 1e-10
        @test ber_mpsk_approx(4, γ) ≈ ber_bpsk_awgn(γ) rtol = 1e-12
        @test ser_mpsk(2, γ) ≈ ber_bpsk_awgn(γ) rtol = 1e-6
        p4 = ber_bpsk_awgn(γ)
        @test ser_mpsk(4, 2γ) ≈ 2 * p4 - p4^2 rtol = 1e-5 atol = 1e-10
        @test ser_mpsk(8, 3γ) ≈ ser_mpsk_craig(8, 3γ) rtol = 1e-5
        @test ber_bpsk_rice(0.0, γ) ≈ ber_bpsk_rayleigh(γ) rtol = 1e-6
        @test ber_bpsk_mrc_rayleigh(1, γ) ≈ ber_bpsk_rayleigh(γ) rtol = 1e-12
        @test ber_bpsk_rice(1e6, γ) ≈ ber_bpsk_awgn(γ) rtol = 1e-3 atol = 1e-10
    end

    # Nearest-neighbour 8-PSK tracks the sector BER at high SNR and is not exact.
    γ12 = 10^(12 / 10)
    exact8 = ber_mpsk_gray(8, γ12)
    approx8 = ber_mpsk_approx(8, γ12)
    @test 0.85 < approx8 / exact8 < 1.05

    @test ber_mpsk_gray(8, 0.0) ≈ 0.5
    @test ber_dbpsk(0.0) ≈ 0.5
    @test ber_dqpsk(0.0) ≈ 0.5
    @test ber_bpsk_rice(2.0, 0.0) ≈ 0.5
    @test ber_bpsk_mrc_rayleigh(4, 0.0) ≈ 0.5
    @test ber_dbpsk(Inf) == 0.0
    @test ber_dqpsk(Inf) == 0.0
    @test ber_bpsk_mrc_rayleigh(3, Inf) == 0.0

    γ = 10^(6 / 10)
    @test ber_mpsk_gray(16, γ) > ber_mpsk_gray(8, γ) > ber_mpsk_gray(4, γ)
    @test ber_dbpsk(γ) > ber_bpsk_awgn(γ)
    @test ber_dqpsk(γ) > ber_qpsk_awgn(γ)
    @test ber_bpsk_rice(0.0, γ) > ber_bpsk_rice(4.0, γ) > ber_bpsk_rice(1e6, γ)
    @test ber_bpsk_mrc_rayleigh(1, γ) > ber_bpsk_mrc_rayleigh(2, γ) >
          ber_bpsk_mrc_rayleigh(4, γ)
    @test ber_dbpsk(10^(10 / 10)) < ber_dbpsk(γ)
    @test ber_dqpsk(10^(10 / 10)) < ber_dqpsk(γ)

    p = psk_offset_probs(8, 5.0)
    @test sum(p) ≈ 1 rtol = 1e-8
    @test p[1] > 0.5
    for d in 1:3
        @test p[d + 1] ≈ p[8 - d + 1] rtol = 1e-8
    end

    @test marcumq1(0, 0) == 1.0
    @test marcumq1(4, 0) == 1.0
    @test marcumq1(0, 2) ≈ exp(-2) rtol = 1e-12
    Random.seed!(RNG_SEED)
    a, b = 1.5, 2.0
    nq = 200_000
    kq = count(_ -> hypot(a + randn(), randn()) > b, 1:nq)
    @test abs(kq / nq - marcumq1(a, b)) < 0.008

    @test_throws ArgumentError ber_mpsk_gray(3, 1.0)
    @test_throws ArgumentError ber_mpsk_approx(2, 1.0)
    @test_throws ArgumentError ber_bpsk_mrc_rayleigh(0, 1.0)
    @test_throws ArgumentError ber_bpsk_rice(-1.0, 1.0)
end

@testset "Gray M-PSK maps and LLRs" begin
    @test_throws ArgumentError PSK(2)
    @test_throws ArgumentError PSK(3)
    @test_throws ArgumentError PSK(6)
    @test_throws ArgumentError PSK(0)
    @test bits_per_symbol(PSK(8)) == 3
    @test bits_per_symbol(PSK(16)) == 4
    @test bits_per_symbol(PSK(32)) == 5

    pts4 = constellation(PSK(4))
    @test pts4 ≈ ComplexF64[1, im, -im, -1]
    @test demodulate(pts4, PSK(4)) == [0, 0, 0, 1, 1, 0, 1, 1]
    @test constellation(PSK(8))[1] ≈ 1

    for M in (4, 8, 16, 32)
        scheme = PSK(M)
        pts = constellation(scheme)
        k = bits_per_symbol(scheme)
        @test sum(abs2, pts) / M ≈ 1 rtol = 1e-12
        @test length(pts) == M
        labels = Int[]
        for lab in 0:(M - 1), j in 1:k
            push!(labels, (lab >> (k - j)) & 1)
        end
        @test modulate(labels, scheme) ≈ pts
        @test demodulate(pts, scheme) == labels
        for i in 1:M
            d2 = [j == i ? Inf : abs2(pts[i] - pts[j]) for j in 1:M]
            j1 = argmin(d2)
            dmin = d2[j1]
            d2[j1] = Inf
            j2 = argmin(d2)
            @test d2[j2] ≈ dmin rtol = 1e-12
            @test count_ones(xor(i - 1, j1 - 1)) == 1
            @test count_ones(xor(i - 1, j2 - 1)) == 1
        end
    end

    @test_throws ArgumentError modulate([0, 1], PSK(8))
    @test_throws ArgumentError modulate([0, 1, 2], PSK(8))

    for scheme in (PSK(8), PSK(16))
        Random.seed!(RNG_SEED)
        k = bits_per_symbol(scheme)
        bits = rand(0:1, 300 * k)
        @test demodulate(modulate(bits, scheme), scheme) == bits
        pts = constellation(scheme)
        for approx in (:exact, :maxlog)
            L = demodulate_soft(pts, scheme, 1e-4; approx=approx)
            @test all(isfinite, L)
            @test Int.(L .< 0) == demodulate(pts, scheme)
        end
    end

    # Off-constellation 8-PSK: exact and max-log differ, and max-log matches an
    # independent nearest-subset evaluation with N0 = 2σ².
    y8 = [0.45 * cis(π / 16)]
    σ2 = 0.5
    Lex = demodulate_soft(y8, PSK(8), σ2; approx=:exact)
    Lmx = demodulate_soft(y8, PSK(8), σ2; approx=:maxlog)
    @test !(Lex ≈ Lmx)
    pts = constellation(PSK(8))
    k = 3
    N0 = 2σ2
    manual = Float64[]
    for bit in 1:k
        d0, d1 = Inf, Inf
        for label in 0:7
            d = abs2(y8[1] - pts[label + 1])
            if ((label >> (k - bit)) & 1) == 0
                d < d0 && (d0 = d)
            else
                d < d1 && (d1 = d)
            end
        end
        push!(manual, (d1 - d0) / N0)
    end
    @test Lmx ≈ manual rtol = 1e-12
    @test demodulate_soft([0.3], PSK(8), σ2) ≈ demodulate_soft(ComplexF64[0.3], PSK(8), σ2)
    @test_throws ArgumentError demodulate_soft(y8, PSK(8), 0.0)
    @test_throws ArgumentError demodulate_soft(y8, PSK(8), 0.2; approx=:nope)

    # Phase decisions away from the boundaries agree with minimum distance.
    Random.seed!(RNG_SEED)
    M = 8
    θ = (rand(100) .- 0.5) .* (0.6 * 2π / M)
    ys = cis.(θ)
    phase_bits = Int[]
    for y in ys
        ang = angle(y)
        ang < 0 && (ang += 2π)
        m = mod(round(Int, ang / (2π / M)), M)
        g = xor(m, m >> 1)
        for j in 1:3
            push!(phase_bits, (g >> (3 - j)) & 1)
        end
    end
    @test demodulate(ys, PSK(8)) == phase_bits
end

@testset "DBPSK and DQPSK" begin
    @test bits_per_symbol(DBPSK()) == 1
    @test bits_per_symbol(DQPSK()) == 2
    @test constellation(DBPSK()) ≈ ComplexF64[1, -1]
    @test constellation(DQPSK()) ≈ ComplexF64[1, im, -im, -1]

    @test modulate(Int[], DBPSK()) ≈ [1]
    @test demodulate(modulate(Int[], DBPSK()), DBPSK()) == Int[]
    @test modulate([0, 0, 1, 1], DBPSK()) ≈ ComplexF64[1, 1, 1, -1, 1]
    @test modulate([1, 1, 1], DBPSK()) ≈ ComplexF64[1, -1, 1, -1]
    @test modulate([0, 1, 1, 1], DQPSK()) ≈ ComplexF64[1, im, -im]
    @test modulate([1, 0], DQPSK()) ≈ ComplexF64[1, -im]
    @test modulate(Int[], DQPSK()) ≈ [1]

    for scheme in (DBPSK(), DQPSK())
        Random.seed!(RNG_SEED)
        k = bits_per_symbol(scheme)
        bits = rand(0:1, 400 * k)
        syms = modulate(bits, scheme)
        @test length(syms) == length(bits) ÷ k + 1
        @test syms[1] ≈ 1
        @test all(s -> abs2(s) ≈ 1, syms)
        @test demodulate(syms, scheme) == bits
        @test demodulate(syms .* cis(0.37), scheme) == bits
        for approx in (:exact, :maxlog)
            L = demodulate_soft(syms, scheme, 0.05; approx=approx)
            @test length(L) == length(bits)
            @test Int.(L .< 0) == bits
        end
    end

    db = modulate([0, 1, 0, 1, 1, 0], DBPSK())
    @test demodulate_soft(db, DBPSK(), 0.4; approx=:exact) ≈
          demodulate_soft(db, DBPSK(), 0.4; approx=:maxlog)
    @test demodulate([1.0, 1.0, -1.0], DBPSK()) == [0, 1]
    @test demodulate(real.(modulate([0, 1, 1, 0], DBPSK())), DBPSK()) == [0, 1, 1, 0]

    y0 = 0.8 + 0.3im
    y1 = -0.2 + 0.7im
    for σ2 in (2.0, 0.35, 0.01)
        function logi0(x)
            x == 0 && return 0.0
            return x + log(besselix(0, x))
        end
        lp = logi0(abs(y0 + y1) / σ2)
        lm = logi0(abs(y0 - y1) / σ2)
        @test demodulate_soft([y0, y1], DBPSK(), σ2) ≈ [lp - lm] rtol = 1e-10
    end

    σq = 0.4
    z0 = 1.0 + 0.0im
    z1 = cis(π / 6)
    deltas = constellation(DQPSK())
    ll = [let x = abs(z0 + z1 * conj(d)) / σq
              x + log(besselix(0, x))
          end for d in deltas]
    logaddexp(a, b) = (m = max(a, b); m + log(exp(a - m) + exp(b - m)))
    Lref = [
        logaddexp(ll[1], ll[2]) - logaddexp(ll[3], ll[4]),
        logaddexp(ll[1], ll[3]) - logaddexp(ll[2], ll[4]),
    ]
    @test demodulate_soft([z0, z1], DQPSK(), σq; approx=:exact) ≈ Lref rtol = 1e-10
    Lmx = demodulate_soft([z0, z1], DQPSK(), σq; approx=:maxlog)
    @test Lmx ≈ [
        max(ll[1], ll[2]) - max(ll[3], ll[4]),
        max(ll[1], ll[3]) - max(ll[2], ll[4]),
    ]
    @test !(Lmx ≈ Lref)
    @test demodulate_soft([1.0, -1.0], DBPSK(), 0.5) ≈
          demodulate_soft(ComplexF64[1.0, -1.0], DBPSK(), 0.5)

    @test_throws ArgumentError modulate([0], DQPSK())
    @test_throws ArgumentError modulate([0, 2], DBPSK())
    @test_throws ArgumentError demodulate(ComplexF64[], DBPSK())
    @test_throws ArgumentError demodulate(ComplexF64[], DQPSK())
    @test_throws ArgumentError demodulate_soft(ComplexF64[], DBPSK(), 0.5)
    @test_throws ArgumentError demodulate_soft(ComplexF64[], DQPSK(), 0.5)
    @test_throws ArgumentError demodulate_soft([1.0, -1.0], DBPSK(), -1.0)

    # A single phase on the whole block is invisible. Independent phases are not.
    Random.seed!(RNG_SEED)
    bits = rand(0:1, 500)
    syms = modulate(bits, DBPSK())
    @test demodulate(syms .* cis(1.1), DBPSK()) == bits
    spun = syms .* cis.(2π .* rand(length(syms)))
    @test ber(bits, demodulate(spun, DBPSK())) > 0.3
end

@testset "Rice fading and MRC" begin
    Random.seed!(7)
    h0 = rice_fading(6, 0)
    Random.seed!(7)
    href = (randn(6) .+ im .* randn(6)) ./ sqrt(2)
    @test h0 ≈ href

    Random.seed!(1)
    draw = randn()
    Random.seed!(1)
    hinf = rice_fading(5, Inf)
    @test all(hinf .== 1)
    @test randn() == draw

    Random.seed!(1)
    hK = rice_fading(20_000, 3.0)
    @test sum(abs2, hK) / length(hK) ≈ 1 rtol = 0.05
    @test real(sum(hK)) / length(hK) ≈ sqrt(3 / 4) rtol = 0.05
    @test abs(imag(sum(hK)) / length(hK)) < 0.02

    Random.seed!(3)
    syms = modulate(zeros(Int, 8), BPSK())
    rx, h = rice_channel(syms, 4.0, 1; K=2.0)
    Random.seed!(3)
    h2 = rice_fading(8, 2.0)
    rx2 = awgn(ComplexF64.(syms) .* h2, 4.0, 1)
    @test h ≈ h2
    @test rx ≈ rx2

    Random.seed!(4)
    _, hR = rayleigh_channel(ones(ComplexF64, 9), Inf, 1)
    Random.seed!(4)
    _, hRice = rice_channel(ones(ComplexF64, 9), Inf, 1; K=0)
    @test hR ≈ hRice

    Random.seed!(1)
    faded, hf = rice_channel(ones(5), Inf, 1; K=2.0)
    Random.seed!(1)
    zf = rice_channel(ones(5), Inf, 1, true; K=2.0)
    @test zf ≈ faded ./ hf
    @test zf ≈ ones(ComplexF64, 5)

    Random.seed!(1)
    raw, hraw = rice_channel(ones(4), 0.0, 1; K=1.0)
    Random.seed!(1)
    raw_only = rice_channel(ones(4), 0.0, 1, false; K=1.0)
    @test raw_only ≈ raw
    @test !(raw_only ≈ raw ./ hraw)

    Random.seed!(1)
    awgn_only = rice_channel([1.0, -1.0], 5.0, 1, true; K=Inf)
    Random.seed!(1)
    @test awgn_only ≈ awgn(ComplexF64[1.0, -1.0], 5.0, 1)

    rx_e, h_e = rice_channel(Float64[], 0.0, 1; K=1.5)
    @test isempty(rx_e) && isempty(h_e)

    @test_throws ArgumentError rice_fading(-1, 1.0)
    @test_throws ArgumentError rice_fading(1, -0.1)
    @test_throws ArgumentError rice_fading(1, NaN)
    @test_throws UndefKeywordError rice_channel([1.0], 0.0, 1)
    @test_throws ArgumentError rice_channel([1.0], 0.0, 0; K=1.0)

    H = reshape(ComplexF64[3, 0.1], 2, 1)
    x = 0.5 + 0.25im
    @test only(mrc(H .* x, H)) ≈ x
    Y = H .* x
    Y[2] += 10
    den = abs2(H[1]) + abs2(H[2])
    @test only(mrc(Y, H)) ≈ x + conj(H[2]) * 10 / den

    @test mrc([2.0 -0.5; 0.0 0.0], [2.0 1.0; 0.0 0.0]) ≈ ComplexF64[1, -0.5]
    @test isnan(only(mrc(zeros(ComplexF64, 1, 1), zeros(ComplexF64, 1, 1))))
    @test_throws ArgumentError mrc(zeros(ComplexF64, 0, 2), zeros(ComplexF64, 0, 2))
    @test_throws ArgumentError mrc(zeros(ComplexF64, 1, 2), zeros(ComplexF64, 1, 3))

    Random.seed!(9)
    z1, _ = mrc_channel(ones(ComplexF64, 7), 2.0, 1, 1; K=0)
    Random.seed!(9)
    @test z1 ≈ rayleigh_channel(ones(ComplexF64, 7), 2.0, 1, true)

    bits = [0, 1, 1, 0, 1, 0]
    syms_b = modulate(bits, BPSK())
    zinf, hinfm = mrc_channel(syms_b, Inf, 1, 3; K=Inf)
    @test zinf ≈ ComplexF64.(syms_b)
    @test size(hinfm) == (3, 6)
    @test all(hinfm .== 1)
    @test demodulate(zinf, BPSK()) == bits

    zempty, hempty = mrc_channel(ComplexF64[], Inf, 1, 2)
    @test isempty(zempty)
    @test size(hempty) == (2, 0)
    @test_throws ArgumentError mrc_channel([1.0], 0.0, 1, 0)
    @test_throws ArgumentError mrc_channel([1.0], 0.0, 0, 2)
end

function _mc_bits(make_rx, snr_db; target=220, max_bits=800_000, batch=40_000, kbits=1)
    k = 0
    n = 0
    while k < target && n < max_bits
        nb = min(batch, max_bits - n)
        nb -= mod(nb, kbits)
        nb > 0 || break
        bits = rand(0:1, nb)
        rx_bits = make_rx(bits, snr_db)
        k += count(!=(0), bits .- rx_bits)
        n += nb
    end
    return k, n
end

@testset "M-PSK, differential, Rice, MRC versus closed forms" begin
    Random.seed!(RNG_SEED)
    for (M, snr_db) in ((8, 5.0), (8, 8.0), (16, 8.0))
        γ = 10^(snr_db / 10)
        scheme = PSK(M)
        k, n = _mc_bits(snr_db; kbits=bits_per_symbol(scheme)) do bits, snr
            demodulate(awgn(modulate(bits, scheme), snr, bits_per_symbol(scheme)), scheme)
        end
        assert_theory_inside(k, n, ber_mpsk_gray(M, γ))
    end

    Random.seed!(RNG_SEED)
    Mser, snr_ser = 8, 6.0
    γs = bits_per_symbol(PSK(Mser)) * 10^(snr_ser / 10)
    kser = 0
    nser = 0
    target = 220
    while kser < target && nser < 200_000
        ns = 4_000
        bits = rand(0:1, ns * 3)
        rxb = demodulate(awgn(modulate(bits, PSK(Mser)), snr_ser, 3), PSK(Mser))
        for i in 1:ns
            sl = (i - 1) * 3
            kser += rxb[(sl + 1):(sl + 3)] != bits[(sl + 1):(sl + 3)]
        end
        nser += ns
    end
    assert_theory_inside(kser, nser, ser_mpsk_craig(Mser, γs))

    Random.seed!(RNG_SEED)
    for snr_db in (4.0, 8.0)
        γ = 10^(snr_db / 10)
        k, n = _mc_bits(snr_db) do bits, snr
            demodulate(awgn(modulate(bits, DBPSK()), snr, 1), DBPSK())
        end
        assert_theory_inside(k, n, ber_dbpsk(γ))
    end

    Random.seed!(RNG_SEED)
    for snr_db in (4.0, 7.0)
        γ = 10^(snr_db / 10)
        k, n = _mc_bits(snr_db; kbits=2) do bits, snr
            demodulate(awgn(modulate(bits, DQPSK()), snr, 2), DQPSK())
        end
        assert_theory_inside(k, n, ber_dqpsk(γ))
    end

    Random.seed!(RNG_SEED)
    for (K, snr_db) in ((0.0, 8.0), (4.0, 8.0))
        γ = 10^(snr_db / 10)
        k, n = _mc_bits(snr_db) do bits, snr
            demodulate(rice_channel(modulate(bits, BPSK()), snr, 1, true; K=K), BPSK())
        end
        assert_theory_inside(k, n, ber_bpsk_rice(K, γ))
    end

    Random.seed!(RNG_SEED)
    for (L, snr_db) in ((2, 6.0), (4, 2.0))
        γ = 10^(snr_db / 10)
        k, n = _mc_bits(snr_db; batch=20_000) do bits, snr
            z, _ = mrc_channel(modulate(bits, BPSK()), snr, 1, L; K=0)
            demodulate(z, BPSK())
        end
        assert_theory_inside(k, n, ber_bpsk_mrc_rayleigh(L, γ))
    end
end
