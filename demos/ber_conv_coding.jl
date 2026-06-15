import CommToolkit: modulate, demodulate, ber, awgn, BPSK, conv_encode, viterbi_decode, snr_db_to_linear
import Plots

n_bits = 10000
snr_db_range = 0:2:10

ber_uncoded = zeros(length(snr_db_range))
ber_coded = zeros(length(snr_db_range))

for (i, snr) in enumerate(snr_db_range)
    tx = rand(0:1, n_bits)
    snr_lin = snr_db_to_linear(snr)

    # Uncoded BPSK
    syms = modulate(tx, BPSK())
    rx_uncoded = awgn(syms, snr, 1)
    bits_rx = demodulate(rx_uncoded, BPSK())
    ber_uncoded[i] = ber(tx, bits_rx)

    # Convolutional coded BPSK (rate 1/2)
    tx_coded = conv_encode(tx)
    syms_coded = modulate(tx_coded, BPSK())
    # Energy per info bit = 2x energy per coded bit (rate 1/2)
    energy_per_bit = 2.0
    noise_var = energy_per_bit / (2 * snr_lin)
    noise_std = sqrt(noise_var)
    rx_coded = syms_coded + noise_std * randn(length(syms_coded))
    bits_rx_coded = demodulate(rx_coded, BPSK())
    decoded = viterbi_decode(bits_rx_coded)
    n_info = length(decoded)
    ber_coded[i] = ber(tx[1:n_info], decoded)
end

p = Plots.plot(snr_db_range, ber_uncoded,
    yaxis=:log, marker=:circle,
    label="Uncoded", legend=:bottomleft)
Plots.plot!(p, snr_db_range, ber_coded,
    yaxis=:log, marker=:square,
    label="Conv code (r=1/2, K=3)")
Plots.xlabel!("Eb/N0 (dB)")
Plots.ylabel!("Bit Error Rate")
Plots.title!("Convolutional Coding over AWGN")

Plots.savefig(joinpath(@__DIR__, "ber_conv_coding.png"))
Plots.display(p)
println("Conv code BER simulation complete.")
