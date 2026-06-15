import CommToolkit: mimo_channel, apply_mimo, zf_detect, mmse_detect, ml_detect
import CommToolkit: modulate, demodulate, ber, awgn, BPSK, QPSK, snr_db_to_linear
import Plots

n_bits = 10000
snr_range = 0:2:16
n_tx, n_rx = 2, 2

ber_zf = Float64[]
ber_mmse = Float64[]

for snr in snr_range
    tx = rand(0:1, n_bits)
    syms = modulate(tx, BPSK())
    noise_var = 1.0 / (2 * snr_db_to_linear(snr))

    err_zf = 0
    err_mmse = 0
    total = 0

    for i in 1:2:n_bits
        x = [syms[i], syms[i+1]]
        H = mimo_channel(n_tx, n_rx)
        y = apply_mimo(x, H, noise_var)

        x_zf = zf_detect(y, H)
        x_mmse = mmse_detect(y, H, noise_var)

        b_zf = demodulate([real(x_zf[1]), real(x_zf[2])], BPSK())
        b_mmse = demodulate([real(x_mmse[1]), real(x_mmse[2])], BPSK())

        err_zf += sum(tx[2i-1:2i] .!= b_zf)
        err_mmse += sum(tx[2i-1:2i] .!= b_mmse)
        total += 2
    end

    push!(ber_zf, err_zf / total)
    push!(ber_mmse, err_mmse / total)
end

# SISO AWGN reference
ber_ref = Float64[]
for snr in snr_range
    tx = rand(0:1, 100000)
    rx = awgn(modulate(tx, BPSK()), snr, 1)
    push!(ber_ref, ber(tx, demodulate(rx, BPSK())))
end

p = Plots.plot(snr_range, ber_ref,
    yaxis=:log, marker=:circle, lw=2, ls=:dash, color=:gray,
    label="SISO AWGN", legend=:bottomleft)
Plots.plot!(p, snr_range, ber_zf,
    yaxis=:log, marker=:square, lw=2,
    label="ZF 2x2")
Plots.plot!(p, snr_range, ber_mmse,
    yaxis=:log, marker=:diamond, lw=2,
    label="MMSE 2x2")
Plots.xlabel!("Eb/N0 (dB)")
Plots.ylabel!("Bit Error Rate")
Plots.title!("MIMO 2x2 BPSK: ZF vs MMSE")

Plots.savefig(joinpath(@__DIR__, "ber_mimo.png"))
Plots.display(p)
println("MIMO BER simulation complete.")
