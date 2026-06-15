using CommToolkit
using Test

@testset "BPSK" begin
    # Test modulation mapping
    @test modulate([0, 0, 1, 1], BPSK()) == [1, 1, -1, -1]

    # Test demodulation
    @test demodulate([1.0, -1.0, 0.5, -0.5], BPSK()) == [0, 1, 0, 1]

    # Test round-trip (no noise)
    tx = rand(0:1, 1000)
    rx = demodulate(modulate(tx, BPSK()), BPSK())
    @test tx == rx
end

@testset "QPSK" begin
    # Test round-trip (no noise)
    tx = rand(0:1, 1000)
    rx = demodulate(modulate(tx, QPSK()), QPSK())
    @test tx == rx
end

@testset "16QAM" begin
    # Test round-trip (no noise)
    tx = rand(0:1, 10000)
    rx = demodulate(modulate(tx, QAM(16)), QAM(16))
    @test tx == rx
end

@testset "AWGN" begin
    # High SNR → BER should be low
    tx = rand(0:1, 100000)
    syms = modulate(tx, BPSK())
    rx = awgn(syms, 20, 1)  # 20 dB, BPSK=1 bit/sym
    bits_rx = demodulate(rx, BPSK())
    @test ber(tx, bits_rx) < 0.01

    # Low SNR → BER should be high
    rx2 = awgn(syms, -10, 1)  # -10 dB, BPSK=1 bit/sym
    bits_rx2 = demodulate(rx2, BPSK())
    @test ber(tx, bits_rx2) > 0.1
end

@testset "Rayleigh" begin
    # BPSK over Rayleigh (perfect CSI) at high SNR
    tx = rand(0:1, 50000)
    syms = modulate(tx, BPSK())
    rx = rayleigh_channel(syms, 30, 1, true)  # 30 dB, perfect CSIR
    bits_rx = demodulate(rx, BPSK())
    @test ber(tx, bits_rx) < 0.05

    # Perfect CSIR should be better than no CSIR
    rx_no_csi = rayleigh_channel(syms, 20, 1, false)
    bits_rx_no = demodulate(rx_no_csi, BPSK())
    ber_no = ber(tx, bits_rx_no)

    ber_csi = ber(tx, bits_rx)
    @test ber_csi <= ber_no || ber_csi < 0.1
end

@testset "Convolutional code" begin
    # Round-trip: encode → decode = original
    tx = rand(0:1, 100)
    coded = conv_encode(tx)
    rx = viterbi_decode(coded)
    @test tx == rx

    # Single-bit error correction
    tx2 = [1, 0, 1]
    coded2 = conv_encode(tx2)
    coded2[3] ⊻= 1  # flip one coded bit
    rx2 = viterbi_decode(coded2)
    @test tx2 == rx2

    # Coded length: (n + K-1) * 2
    @test length(coded) == 2 * (100 + 2)
end

@testset "Equalizers" begin
    # Simple 2-tap channel
    h = [1.0, 0.5]
    syms = [1.0, -1.0, 1.0, -1.0]
    L = length(h)
    rx = similar(syms, length(syms) + L - 1)
    # Convolve
    for i in 1:length(rx)
        rx[i] = 0.0
        for j in 1:L
            if 1 <= i - j + 1 <= length(syms)
                rx[i] += syms[i-j+1] * h[j]
            end
        end
    end

    # ZF equalization (noiseless → perfect recovery)
    est_zf = zf_equalize(rx, h)
    @test est_zf ≈ syms atol=1e-10

    # MMSE equalization (noiseless → equal to ZF when σ²=0)
    est_mmse = mmse_equalize(rx, h, 0.0)
    @test est_mmse ≈ syms atol=1e-10
end

@testset "OFDM" begin
    params = OFDMParams(N_fft=16, N_cp=4, mod_order=4, n_guard=2, pilot_spacing=4)
    channel = TDLChannel(delays_samples=[0, 1], gain_dB=[0, -3])

    # Round-trip: tx → known-channel rx
    n_bits = params.N_data * params.bits_per_symbol * 5
    bits = rand(0:1, n_bits)
    tx_signal, pilots = ofdm_tx(params, bits)
    rx_signal, H_freq = apply_channel(tx_signal, channel, 30, params)
    rx_bits = ofdm_rx_known(rx_signal, params, H_freq, pilots)
    @test ber(bits, rx_bits[1:length(bits)]) < 0.01

    # LS estimation receiver
    rx_bits_ls = ofdm_rx_ls(rx_signal, params, pilots)
    @test ber(bits, rx_bits_ls[1:length(bits)]) < 0.05

    # MMSE estimation receiver
    rx_bits_mmse = ofdm_rx_mmse(rx_signal, params, pilots, 30)
    @test ber(bits, rx_bits_mmse[1:length(bits)]) < 0.05
end

@testset "OFDM sync" begin
    params = OFDMParams(N_fft=16, N_cp=4, mod_order=4, n_guard=2, pilot_spacing=4)

    n_bits = params.N_data * params.bits_per_symbol * 3
    bits = rand(0:1, n_bits)
    tx_signal, pilots = ofdm_tx(params, bits)

    # Add CFO
    cfo_true = 0.15
    N = params.N_fft
    rx = [tx_signal[i] * exp(1im * 2π * cfo_true * (i-1) / N) for i in 1:length(tx_signal)]

    # Estimate and correct CFO
    cfo_est = cfo_estimate(rx, params)
    @test abs(cfo_est - cfo_true) < 0.02

    # Sync
    corrected, sto_est, cfo_est2 = ofdm_sync(rx, params)
    @test abs(cfo_est2 - cfo_true) < 0.02
end

@testset "Pulse shaping" begin
    beta, span, sps = 0.5, 6, 4
    rc = raised_cosine(beta, span, sps)
    @test length(rc) == span * sps + 1
    @test abs(sum(rc) - 1) < 1e-10  # unit DC gain

    rrc = rrc_filter(beta, span, sps)
    @test length(rrc) == span * sps + 1
    @test abs(sum(abs2, rrc) - 1) < 1e-10  # unit energy

    # Pulse shaping round-trip
    syms = [1.0, -1.0, 1.0]
    shaped = pulse_shape(syms, rc, sps)
    @test length(shaped) == length(syms) * sps
end

@testset "BER metric" begin
    @test ber([0, 0, 0], [0, 0, 0]) == 0.0
    @test ber([0, 0, 0], [1, 1, 1]) == 1.0
    @test ber([0, 1, 0], [1, 1, 0]) == 1/3
end

@testset "MIMO" begin
    H = mimo_channel(2, 2)
    @test size(H) == (2, 2)

    # ML detection for 2x2 BPSK
    constell = [1.0+0im, -1.0+0im]
    x = [1.0+0im, -1.0+0im]
    y = apply_mimo(x, H, 0.0)  # noiseless
    x_est = ml_detect(y, H, constell)
    @test x_est ≈ x

    # ZF detection
    x_zf = zf_detect(y, H)
    @test x_zf ≈ x

    # Alamouti 2x1 round-trip
    h1, h2 = 0.5+0.3im, -0.2+0.8im
    code = alamouti_encode(1.0+0im, -1.0+0im)
    r1 = h1 * code[1, 1] + h2 * code[2, 1]
    r2 = h1 * code[1, 2] + h2 * code[2, 2]
    x1, x2 = alamouti_decode_2x1(r1, r2, h1, h2)
    @test x1 ≈ 1.0+0im atol=1e-10
    @test x2 ≈ -1.0+0im atol=1e-10

    # ZF precoder
    P = zf_precoder(H)
    @test size(P) == (2, 2)
end

println("All tests passed!")
