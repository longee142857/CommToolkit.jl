# Flat Clarke/Jakes correlation and i.i.d. Nakagami-m.
# J0 is SpecialFunctions.besselj0 (test dependency only). The Nakagami BER
# integral lives in reference.jl. Monte Carlo points use the same
# Clopper–Pearson rule as the Rice tests (α = 1e-4).

using Statistics

@testset "Nakagami-m formulas" begin
    for db in (0.0, 4.0, 8.0, 12.0)
        γ = 10^(db / 10)
        @test ber_bpsk_nakagami(1, γ) ≈ ber_bpsk_rayleigh(γ) rtol = 1e-10
        @test ber_bpsk_nakagami(1, γ) ≈ ber_bpsk_mrc_rayleigh(1, γ) rtol = 1e-10
        for m in 2:5
            @test ber_bpsk_nakagami(m, γ) ≈ ber_bpsk_mrc_rayleigh(m, γ / m) rtol = 1e-8
        end
    end
    @test ber_bpsk_nakagami(0.5, 0.0) ≈ 0.5
    @test ber_bpsk_nakagami(3, 0.0) ≈ 0.5
    @test ber_bpsk_nakagami(2, Inf) == 0.0
    @test ber_bpsk_nakagami(Inf, 1.5) ≈ ber_bpsk_awgn(1.5) rtol = 1e-12
    @test ber_bpsk_nakagami(10_000, 1.0) ≈ ber_bpsk_awgn(1.0) rtol = 1e-3

    γ = 10^(6 / 10)
    @test ber_bpsk_nakagami(0.5, γ) > ber_bpsk_nakagami(1, γ) >
          ber_bpsk_nakagami(2, γ) > ber_bpsk_nakagami(4, γ) >
          ber_bpsk_nakagami(Inf, γ)

    @test_throws ArgumentError ber_bpsk_nakagami(0.4, 1.0)
    @test_throws ArgumentError ber_bpsk_nakagami(1.0, -0.1)
end

@testset "Nakagami envelope and complex coefficients" begin
    Random.seed!(1)
    draw = randn()
    Random.seed!(1)
    @test nakagami_envelope(0, 2.5) == Float64[]
    @test nakagami_fading(0, 0.5) == ComplexF64[]
    @test nakagami_envelope(4, Inf) == ones(4)
    @test randn() == draw

    Random.seed!(1)
    hinf = nakagami_fading(6, Inf)
    @test all(abs2.(hinf) .≈ 1)
    Random.seed!(1)
    @test hinf ≈ cis.(2π .* rand(6))

    Random.seed!(2)
    h1 = nakagami_fading(8, 1)
    Random.seed!(2)
    href = (randn(8) .+ im .* randn(8)) ./ sqrt(2)
    @test h1 ≈ href
    Random.seed!(2)
    @test abs.(h1) ≈ nakagami_envelope(8, 1)

    Random.seed!(3)
    Rhalf = nakagami_envelope(7, 1 // 2)
    Random.seed!(3)
    @test Rhalf ≈ abs.(randn(7))

    Random.seed!(4)
    h2 = nakagami_fading(5, 2)
    Random.seed!(4)
    R2 = nakagami_envelope(5, 2)
    @test h2 ≈ R2 .* cis.(2π .* rand(5))
    # m = 2: four standard normals per sample, in that order, then the phase.
    Random.seed!(4)
    z = randn(4 * 5)
    manual = [sqrt(sum(z[(i - 1) * 4 + j]^2 for j in 1:4) / 4) for i in 1:5]
    @test R2 ≈ manual

    Random.seed!(5)
    a = nakagami_fading(20, 1.3)
    Random.seed!(5)
    @test a == nakagami_fading(20, 1.3)

    for (m, n) in ((0.5, 80_000), (0.7, 80_000), (1.5, 40_000), (1.7, 40_000), (4.0, 40_000))
        Random.seed!(11)
        R = nakagami_envelope(n, m)
        r2 = abs2.(R)
        @test mean(r2) ≈ 1 rtol = 0.03
        @test mean(r2 .^ 2) ≈ (m + 1) / m rtol = 0.05
        @test all(R .>= 0)
    end

    Random.seed!(6)
    syms = modulate(zeros(Int, 8), BPSK())
    rx, h = nakagami_channel(syms, 4.0, 1; m=2.5)
    Random.seed!(6)
    h2 = nakagami_fading(8, 2.5)
    rx2 = awgn(ComplexF64.(syms) .* h2, 4.0, 1)
    @test h ≈ h2
    @test rx ≈ rx2

    Random.seed!(7)
    ray = rayleigh_channel(ones(ComplexF64, 9), 2.0, 1)
    Random.seed!(7)
    rice = rice_channel(ones(ComplexF64, 9), 2.0, 1; K=0)
    Random.seed!(7)
    nak = nakagami_channel(ones(ComplexF64, 9), 2.0, 1; m=1)
    @test ray[2] ≈ rice[2] ≈ nak[2]
    @test ray[1] ≈ rice[1] ≈ nak[1]

    Random.seed!(8)
    faded, hf = nakagami_channel(ones(5), Inf, 1; m=1.5)
    Random.seed!(8)
    zf = nakagami_channel(ones(5), Inf, 1, true; m=1.5)
    @test zf ≈ faded ./ hf
    @test zf ≈ ones(ComplexF64, 5)

    Random.seed!(9)
    raw, hraw = nakagami_channel(ones(4), 0.0, 1; m=0.8)
    Random.seed!(9)
    raw_only = nakagami_channel(ones(4), 0.0, 1, false; m=0.8)
    @test raw_only ≈ raw
    @test !(raw_only ≈ raw ./ hraw)

    Random.seed!(1)
    phase_only = nakagami_channel([1.0, -1.0], 5.0, 1, true; m=Inf)
    Random.seed!(1)
    # Unit amplitude: perfect CSI rotates the awgn noise and leaves the symbols.
    spun = cis.(2π .* rand(2))
    @test phase_only ≈ awgn(ComplexF64[1.0, -1.0] .* spun, 5.0, 1) ./ spun

    rx_e, h_e = nakagami_channel(Float64[], 0.0, 1; m=2)
    @test isempty(rx_e) && isempty(h_e)

    @test_throws ArgumentError nakagami_envelope(-1, 1)
    @test_throws ArgumentError nakagami_envelope(1, 0.4)
    @test_throws ArgumentError nakagami_envelope(1, NaN)
    @test_throws ArgumentError nakagami_fading(1, -Inf)
    @test_throws ArgumentError nakagami_channel([1.0], 0.0, 0; m=1)
    @test_throws UndefKeywordError nakagami_channel([1.0], 0.0, 1)
end

function _mc_nakagami(m, snr_db; target=220, max_bits=800_000, batch=40_000)
    k = 0
    n = 0
    while k < target && n < max_bits
        nb = min(batch, max_bits - n)
        bits = rand(0:1, nb)
        rx = nakagami_channel(modulate(bits, BPSK()), snr_db, 1, true; m=m)
        k += count(!=(0), bits .- demodulate(rx, BPSK()))
        n += nb
    end
    return k, n
end

@testset "Nakagami-m BPSK versus the MGF integral" begin
    Random.seed!(RNG_SEED)
    for (m, snr_db) in ((0.5, 10.0), (1.5, 8.0), (2, 6.0), (1.3, 8.0))
        γ = 10^(snr_db / 10)
        k, n = _mc_nakagami(m, snr_db)
        assert_theory_inside(k, n, ber_bpsk_nakagami(m, γ))
    end
end

function _dft_bin_power(h::AbstractVector{ComplexF64}, k::Int)
    n = length(h)
    s = 0.0 + 0.0im
    rot = cis(-2π * k / n)
    p = 1.0 + 0.0im
    @inbounds for x in h
        s += x * p
        p *= rot
    end
    return abs2(s) / n^2
end

@testset "Jakes fading" begin
    Random.seed!(1)
    draw = randn()
    Random.seed!(1)
    @test jakes_fading(0, 10.0, 1e-3) == ComplexF64[]
    @test randn() == draw

    Random.seed!(2)
    h0 = jakes_fading(6, 0.0, 1.0)
    Random.seed!(2)
    c0 = (randn() + im * randn()) / sqrt(2)
    @test all(h0 .== c0)

    Random.seed!(3)
    h1 = jakes_fading(1, 12.0, 1e-3)
    Random.seed!(3)
    @test only(h1) ≈ (randn() + im * randn()) / sqrt(2)

    # Odd length hits the non-Nyquist bin partition. Power is checked on the long record below.
    Random.seed!(4)
    hodd = jakes_fading(1001, 20.0, 1e-3)
    @test length(hodd) == 1001
    @test all(isfinite, hodd)

    @test_throws ArgumentError jakes_fading(-1, 1.0, 1.0)
    @test_throws ArgumentError jakes_fading(4, -0.1, 1.0)
    @test_throws ArgumentError jakes_fading(4, NaN, 1.0)
    @test_throws ArgumentError jakes_fading(4, 1.0, 0.0)
    @test_throws ArgumentError jakes_fading(4, 1.0, -1.0)
    @test_throws ArgumentError jakes_fading(4, 1.0, Inf)
    @test_throws ArgumentError jakes_fading(4, 0.5, 1.0)
    @test_throws ArgumentError jakes_fading(4, Inf, 1.0)

    # 50 Hz Doppler, 10 kHz sampling: fd*Ts = 0.005. Bin 2000/4096 = 0.488
    # cycles/sample sits outside the Doppler and is a zeroed DFT cell.
    Random.seed!(5)
    hband = jakes_fading(4096, 50.0, 1e-4)
    @test _dft_bin_power(hband, 2000) < 1e-20

    # U-shape: averaged power in the last occupied bin exceeds DC.
    N = 2048
    fd = 0.05
    Ts = 1.0
    k_edge = floor(Int, fd * N)
    dc = 0.0
    edge = 0.0
    Random.seed!(6)
    for _ in 1:24
        h = jakes_fading(N, fd, Ts)
        dc += _dft_bin_power(h, 0)
        edge += _dft_bin_power(h, k_edge)
    end
    @test edge > 4 * dc

    # Sample autocorrelation against J0. Sixteen records of 8192, fd*Ts = 0.02.
    # A separate stress of this average stayed under 0.03 in absolute error.
    Nacf = 8192
    fdacf = 0.02
    lags = [0, 4, 8, 16, 24, 40, 60]
    acc = zeros(length(lags))
    acc_i = zeros(length(lags))
    Random.seed!(RNG_SEED)
    for _ in 1:16
        h = jakes_fading(Nacf, fdacf, 1.0)
        for (i, lag) in pairs(lags)
            s = 0.0 + 0.0im
            @inbounds for t in 1:(Nacf - lag)
                s += h[t + lag] * conj(h[t])
            end
            s /= (Nacf - lag)
            acc[i] += real(s)
            acc_i[i] += imag(s)
        end
    end
    acc ./= 16
    acc_i ./= 16
    theory = [besselj0(2π * fdacf * lag) for lag in lags]
    @test acc ≈ theory atol = 0.05
    @test acc_i ≈ zeros(length(lags)) atol = 0.05
    @test theory[1] ≈ 1
    @test any(<(0), theory)

    # Marginal of a long correlated record: conditional BPSK BER tracks Rayleigh.
    Random.seed!(8)
    hmarg = jakes_fading(80_000, 0.05, 1.0)
    γb = 10^(8 / 10)
    pb = mean(qfunc(sqrt(2 * abs2(hk) * γb)) for hk in hmarg)
    @test pb ≈ ber_bpsk_rayleigh(γb) atol = 0.006
    @test mean(abs2, hmarg) ≈ 1 rtol = 0.04
    @test mean(abs2.(hmarg) .^ 2) ≈ 2 rtol = 0.08

    Random.seed!(9)
    syms = modulate([0, 1, 1, 0, 1, 0], BPSK())
    rx, h = jakes_channel(syms, 3.0, 1, 25.0, 1e-4)
    Random.seed!(9)
    h2 = jakes_fading(6, 25.0, 1e-4)
    rx2 = awgn(ComplexF64.(syms) .* h2, 3.0, 1)
    @test h ≈ h2
    @test rx ≈ rx2

    Random.seed!(10)
    faded, hf = jakes_channel(ones(5), Inf, 1, 0.0, 1.0)
    Random.seed!(10)
    zf = jakes_channel(ones(5), Inf, 1, 0.0, 1.0, true)
    @test all(hf .== hf[1])
    @test zf ≈ ones(ComplexF64, 5)
    @test zf ≈ faded ./ hf

    Random.seed!(11)
    raw, hraw = jakes_channel(ones(4), 0.0, 1, 10.0, 1e-3)
    Random.seed!(11)
    raw_only = jakes_channel(ones(4), 0.0, 1, 10.0, 1e-3, false)
    @test raw_only ≈ raw
    @test !(raw_only ≈ raw ./ hraw)

    empty_rx, empty_h = jakes_channel(Float64[], 0.0, 1, 10.0, 1e-3)
    @test isempty(empty_rx) && isempty(empty_h)
    @test_throws ArgumentError jakes_channel([1.0], 0.0, 0, 1.0, 1.0)
    @test_throws ArgumentError jakes_channel([1.0], 0.0, 1, 1.0, 1.0)
end
