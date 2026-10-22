# Symbol-spaced LMS and the finite-block DFE.
#
# apply_fir / fir_channel build the same convolution the block equalizers invert.
# LMS is supervised and freezes after the training prefix. The DFE feeds back
# hard decisions; a wrong one is allowed to move later symbols.

function _lms_regressor(received, m::Int, n_taps::Int, delay::Int)
    T = eltype(received)
    u = zeros(T, n_taps)
    k = m + delay
    for i in 1:n_taps
        idx = k - i + 1
        if 1 <= idx <= length(received)
            u[i] = received[idx]
        end
    end
    return u
end

function _bpsk_bits(symbols)
    return Int.(real.(symbols) .< 0)
end

@testset "LMS update" begin
    # Hand calculation, real: w <- w + μ e u, y = wᵀu before the update.
    received = [0.4, -0.2, 1.1]
    training = [1.0, -1.0]
    y, trajectory = lms_equalize(received, training, 2, 0.1; tap_trajectory=true)
    @test y isa Vector{Float64}
    @test y ≈ [0.0, -0.008] atol = 1e-12
    @test trajectory[:, 1] ≈ [0.04, 0.0] atol = 1e-12
    @test trajectory[:, 2] ≈ [0.05984, -0.03968] atol = 1e-12

    # One complex step locks the conjugate on the regressor, not on the error.
    yc, tc = lms_equalize(ComplexF64[1 + im], [1.0], 1, 0.25; tap_trajectory=true)
    @test yc ≈ [0.0]
    @test tc[:, 1] ≈ [0.25 * (1 - im)]

    y0, t0 = lms_equalize([0.4], [1.0], 1, 0.1; w0=[0.5], tap_trajectory=true)
    @test y0 ≈ [0.2] atol = 1e-12
    @test t0[:, 1] ≈ [0.532] atol = 1e-12

    @test_throws ArgumentError lms_equalize([1.0], [1.0], 0, 0.1)
    @test_throws ArgumentError lms_equalize([1.0], [1.0], 2, 0.0)
    @test_throws ArgumentError lms_equalize([1.0], [1.0], 2, -0.1)
    @test_throws ArgumentError lms_equalize([1.0], [1.0], 2, Inf)
    @test_throws ArgumentError lms_equalize([1.0], [1.0], 2, NaN)
    @test_throws ArgumentError lms_equalize([1.0], [1.0], 2, 0.1; delay=-1)
    @test_throws ArgumentError lms_equalize([1.0], Float64[], 1, 0.1)
    @test_throws ArgumentError lms_equalize(Float64[], [1.0], 1, 0.1)
    @test_throws ArgumentError lms_equalize([1.0, 2.0], [1.0], 1, 0.1; n_symbols=0)
    @test_throws ArgumentError lms_equalize([1.0], [1.0, 1.0], 1, 0.1; delay=1)
    @test_throws ArgumentError lms_equalize([1.0, 2.0], [1.0], 1, 0.1; w0=[0.0, 0.0])
end

@testset "LMS on a symbol-rate FIR" begin
    Random.seed!(2)
    h = [1.0, 0.5, 0.25]
    n = 2000
    symbols = rand([-1.0, 1.0], n)
    received = apply_fir(symbols, h)
    n_taps = 7
    y, trajectory = lms_equalize(received, symbols, n_taps, 0.05; delay=0, tap_trajectory=true)
    @test length(y) == n
    @test size(trajectory) == (n_taps, n)

    early = sum(abs2, y[1:100] .- symbols[1:100]) / 100
    late = sum(abs2, y[end - 199:end] .- symbols[end - 199:end]) / 200
    @test late < early / 20
    @test sqrt(late) < 0.05

    # Batch least squares for the same regressor the docstring defines.
    U = zeros(n, n_taps)
    for m in 1:n
        U[m, :] = _lms_regressor(received, m, n_taps, 0)
    end
    w_ls = U \ symbols
    w_final = trajectory[:, end]
    @test norm(w_final .- w_ls) < 0.05 * norm(w_ls)
    @test real(dot(w_final, w_ls)) > 0.99 * norm(w_final) * norm(w_ls)

    # Training outputs use the weights from before that update. The next
    # symbol, still inside a longer block, uses the frozen final weights.
    n_train = 1500
    y_pay, traj_pay = lms_equalize(
        received, symbols[1:n_train], n_taps, 0.05;
        n_symbols=n, tap_trajectory=true,
    )
    @test size(traj_pay, 2) == n_train
    u_last = _lms_regressor(received, n_train, n_taps, 0)
    w_before = traj_pay[:, n_train - 1]
    @test y_pay[n_train] ≈ sum(w_before[i] * u_last[i] for i in 1:n_taps) atol = 1e-12
    w_frozen = traj_pay[:, n_train]
    for m in (n_train + 1):n
        u = _lms_regressor(received, m, n_taps, 0)
        @test y_pay[m] ≈ sum(w_frozen[i] * u[i] for i in 1:n_taps) atol = 1e-12
    end
    payload_mse = sum(abs2, y_pay[n_train + 1:end] .- symbols[n_train + 1:end]) / (n - n_train)
    @test payload_mse < 0.05^2
    @test ber(_bpsk_bits(symbols[n_train + 1:end]), demodulate(y_pay[n_train + 1:end], BPSK())) == 0.0

    # Cursor is not the first tap: delay 0 does not settle, delay 2 does.
    Random.seed!(5)
    h_late = [0.2, 1.0, 0.3]
    s_late = rand([-1.0, 1.0], 3000)
    rx_late = apply_fir(s_late, h_late)
    function _tail_rmse(delay)
        y_d = lms_equalize(rx_late, s_late, 9, 0.04; delay=delay)
        err = y_d[end - 299:end] .- s_late[end - 299:end]
        return sqrt(sum(abs2, err) / length(err))
    end
    @test _tail_rmse(0) > 0.5
    @test _tail_rmse(2) < 0.1

    # Complex baseband, symbol-rate FIR, no noise.
    Random.seed!(4)
    hc = [1.0 + 0.2im, 0.3 - 0.1im]
    qbits = rand(0:1, 6000)
    qsyms = modulate(qbits, QPSK())
    rxc = apply_fir(qsyms, hc)
    yc = lms_equalize(rxc, qsyms, 8, 0.08)
    tail_c = yc[end - 299:end] .- qsyms[end - 299:end]
    @test eltype(yc) <: Complex
    @test sqrt(sum(abs2, tail_c) / length(tail_c)) < 1e-2
    @test demodulate(yc[end - 299:end], QPSK()) == qbits[end - 599:end]

    # High SNR through fir_channel. The acquisition prefix is not required to be error-free.
    Random.seed!(11)
    bits = rand(0:1, 2000)
    s_hi = modulate(bits, BPSK())
    h_hi = [1.0, 0.4]
    rx_hi, = fir_channel(s_hi, h_hi, 40.0, 1)
    y_hi = lms_equalize(rx_hi, s_hi, 5, 0.05)
    @test ber(bits[end - 199:end], demodulate(y_hi[end - 199:end], BPSK())) == 0.0
end

@testset "LMS step-size and settled error" begin
    # One fixed record. A larger stable step leaves a larger tail MSE.
    # This is not the Widrow excess-MSE formula.
    h = [1.0, 0.35]
    n = 6000
    symbols = [((i * 17) % 5) < 2 ? 1.0 : -1.0 for i in 1:n]
    clean = apply_fir(symbols, h)
    received = clean .+ [0.025 * sin(0.31 * i) for i in eachindex(clean)]
    function _tail_mse(mu)
        y = lms_equalize(received, symbols, 5, mu; delay=0)
        return sum(abs2, y[end - 999:end] .- symbols[end - 999:end]) / 1000
    end
    mse_large = _tail_mse(0.18)
    mse_small = _tail_mse(0.02)
    early = let
        y = lms_equalize(received, symbols, 5, 0.02; delay=0)
        sum(abs2, y[1:200] .- symbols[1:200]) / 200
    end
    @test mse_large > 1.5 * mse_small
    @test mse_small < early / 5
    @test mse_large < 0.02
end

@testset "Block DFE" begin
    bpsk = [-1.0, 1.0]

    h = [1.0]
    symbols = [1.0, -1.0, 1.0, -1.0]
    received = apply_fir(symbols, h)
    coef = dfe_equalizer(h, length(symbols), 0.0)
    @test coef.noise_var == 0.0
    @test coef.feedforward ≈ I(length(symbols)) atol = 1e-12
    @test coef.feedback ≈ zeros(length(symbols), length(symbols)) atol = 1e-12
    @test zf_dfe(received, h; alphabet=bpsk) ≈ symbols atol = 1e-12

    for (hc, n) in (
        ([1.0, 0.5], 8),
        ([1.0, 0.5, -0.25], 12),
        ([0.3, 1.0, 0.4], 10),
        ([1.0, -1.2, 0.3, 0.2], 15),
        ([0.0, 1.0, 0.2], 9),
    )
        Random.seed!(1)
        s = rand([-1.0, 1.0], n)
        rx = apply_fir(s, hc)
        design = dfe_equalizer(hc, n, 0.0)
        @test size(design.feedforward) == (n, length(rx))
        @test size(design.feedback) == (n, n)
        @test design.feedback ≈ tril(design.feedback, -1) atol = 1e-12
        y = zf_dfe(rx, hc; alphabet=bpsk)
        @test y ≈ s atol = 1e-9
        @test mmse_dfe(rx, hc, 0.0; alphabet=bpsk) ≈ y atol = 1e-12
        @test zf_dfe(rx, hc; alphabet=bpsk, training=s) ≈ y atol = 1e-12
        @test mmse_dfe(rx, hc, 1e-8; alphabet=bpsk) ≈ s atol = 1e-5
    end

    # Symbol-rate FIR glue, including the noiseless fir_channel variance.
    h = [1.0, 0.5, -0.25]
    symbols = [isodd(i) ? 1.0 : -1.0 for i in 1:11]
    rx_inf, noise_var = fir_channel(symbols, h, Inf)
    @test noise_var == 0.0
    @test rx_inf ≈ apply_fir(symbols, h)
    @test zf_dfe(rx_inf, h; alphabet=bpsk) ≈ symbols atol = 1e-10
    @test mmse_dfe(rx_inf, h, noise_var; alphabet=bpsk) ≈ symbols atol = 1e-10

    hc = [1.0 + 0.2im, 0.4 - 0.1im, -0.3 + 0.05im]
    qbits = [0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0]
    qsyms = modulate(qbits, QPSK())
    rxc = apply_fir(qsyms, hc)
    yq = zf_dfe(rxc, hc; alphabet=constellation(QPSK()))
    @test yq ≈ qsyms atol = 1e-8
    @test demodulate(yq, QPSK()) == qbits

    mixed_h = [0.5 + 0.0im, 0.2im]
    mixed_s = [1.0, -1.0, 1.0, -1.0]
    @test zf_dfe(apply_fir(mixed_s, mixed_h), mixed_h; alphabet=bpsk) ≈ mixed_s atol = 1e-10

    # Entries on and above the feedback diagonal are not read.
    coef_mut = dfe_equalizer([1.0, 0.5], 4, 0.0)
    s4 = [1.0, -1.0, 1.0, -1.0]
    rx4 = apply_fir(s4, [1.0, 0.5])
    y4 = dfe_equalize(rx4, coef_mut; alphabet=bpsk)
    coef_mut.feedback[1, 2] = 50.0
    @test dfe_equalize(rx4, coef_mut; alphabet=bpsk) ≈ y4 atol = 1e-12

    # Tie at 0 keeps the earlier alphabet entry, and that decision enters the next symbol.
    tie = DFECoefficients([0.0 0.0; 0.0 0.0], [0.0 0.0; 1.0 0.0], 0.0)
    @test dfe_equalize([0.0, 0.0], tie; alphabet=[1.0, -1.0]) ≈ [0.0, -1.0]
    @test dfe_equalize([0.0, 0.0], tie; alphabet=[-1.0, 1.0]) ≈ [0.0, 1.0]

    # Wrong feedback on symbol 1. The soft value at symbol 1 is unchanged;
    # later hard decisions run the other way on this post-cursor.
    h_prop = [1.0, 0.9]
    s_prop = [isodd(i) ? 1.0 : -1.0 for i in 1:12]
    rx_prop = apply_fir(s_prop, h_prop)
    good = zf_dfe(rx_prop, h_prop; alphabet=bpsk)
    bad = zf_dfe(rx_prop, h_prop; alphabet=bpsk, training=[-1.0])
    @test good ≈ s_prop atol = 1e-10
    @test bad[1] ≈ s_prop[1] atol = 1e-10
    @test abs(bad[2] - s_prop[2]) > 0.5
    @test ber(_bpsk_bits(s_prop), demodulate(good, BPSK())) == 0.0
    @test ber(_bpsk_bits(s_prop), demodulate(bad, BPSK())) > 0.5

    @test_throws ArgumentError dfe_equalizer([0.0, 0.0], 3, 0.0)
    @test_throws ArgumentError dfe_equalizer(ComplexF64[0.0, 0.0], 2, 0.0)
    @test_throws ArgumentError dfe_equalizer(Float64[], 2, 0.0)
    @test_throws ArgumentError dfe_equalizer([1.0], 0, 0.0)
    @test_throws ArgumentError dfe_equalizer([1.0], 2, -0.1)
    @test_throws ArgumentError dfe_equalizer([1.0], 2, Inf)
    @test_throws ArgumentError zf_dfe([1.0], [1.0, 0.5, 0.25]; alphabet=bpsk)
    @test_throws ArgumentError zf_dfe(apply_fir([1.0, -1.0], [1.0]), [1.0]; alphabet=Float64[])
    @test_throws ArgumentError zf_dfe([1.0, -1.0], [1.0]; alphabet=bpsk, training=Float64[])
    @test_throws ArgumentError zf_dfe([1.0], [1.0]; alphabet=bpsk, training=[1.0, -1.0])
    @test_throws ArgumentError dfe_equalize([0.0], tie; alphabet=bpsk)
    bad_fb = DFECoefficients(reshape([1.0, 0.0], 1, 2), zeros(2, 2), 0.0)
    @test_throws ArgumentError dfe_equalize([0.0, 0.0], bad_fb; alphabet=bpsk)
    empty_ff = DFECoefficients(zeros(0, 1), zeros(0, 0), 0.0)
    @test_throws ArgumentError dfe_equalize(Float64[], empty_ff; alphabet=bpsk)
    @test_throws ArgumentError mmse_dfe(apply_fir([1.0, -1.0], [1.0, 0.2]), [1.0, 0.2], -1.0; alphabet=bpsk)
end
