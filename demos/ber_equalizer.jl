import CommToolkit: modulate, demodulate, ber, BPSK, zf_equalize, mmse_equalize, snr_db_to_linear
import Plots

n_bits = 2000
snr_db_range = 0:5:30

ber_zf = zeros(length(snr_db_range))
ber_mmse = zeros(length(snr_db_range))

# 2-tap multipath channel: main path + echo
h = [1.0, 0.7]

for (i, snr) in enumerate(snr_db_range)
    tx = rand(0:1, n_bits)
    syms = modulate(tx, BPSK())

    # Multipath convolution with h
    L = length(h)
    rx = zeros(length(syms) + L - 1)
    for k in 1:length(rx)
        val = 0.0
        for j in 1:L
            idx = k - j + 1
            if 1 <= idx <= length(syms)
                val += syms[idx] * h[j]
            end
        end
        rx[k] = val
    end

    # Add noise
    noise_var = 1.0 / (2 * snr_db_to_linear(snr))
    rx += sqrt(noise_var) * randn(length(rx))

    # Equalize (may produce slightly different lengths due to boundary)
    est_zf = zf_equalize(rx, h)
    est_mmse = mmse_equalize(rx, h, noise_var)

    bits_zf = demodulate(est_zf, BPSK())
    bits_mmse = demodulate(est_mmse, BPSK())

    n_info = min(length(bits_zf), length(bits_mmse), n_bits)
    ber_zf[i] = ber(tx[1:n_info], bits_zf[1:n_info])
    ber_mmse[i] = ber(tx[1:n_info], bits_mmse[1:n_info])
end

p = Plots.plot(snr_db_range, ber_zf,
    yaxis=:log, marker=:circle,
    label="ZF", legend=:bottomleft)
Plots.plot!(p, snr_db_range, ber_mmse,
    yaxis=:log, marker=:square,
    label="MMSE")
Plots.xlabel!("Eb/N0 (dB)")
Plots.ylabel!("Bit Error Rate")
Plots.title!("BPSK over Multipath: ZF vs MMSE")

Plots.savefig(joinpath(@__DIR__, "ber_equalizer.png"))
Plots.display(p)
println("Equalizer BER simulation complete.")
