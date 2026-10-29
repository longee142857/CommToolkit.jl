# Teaching CDL (TR 38.901 CDL-A/B/C delay–power profiles): profile summary,
# sample-quantized taps, a static drop on the textbook OFDM chain, and the
# flat-degenerate Rayleigh cross-check. No Plots. The stream is not seeded.
# The CDL here is the teaching subset: delays and powers only, no angles.

using CommToolkit
using Printf

println("== TR 38.901 NLOS CDL profiles (normalized delays, unit RMS delay spread) ==")
for name in (:A, :B, :C)
    p = cdl_profile(name)
    pmax, imax = findmax(exp10.(p.powers_db ./ 10))
    @printf("CDL-%s: %2d clusters, strongest %5.1f dB at normalized delay %.4f\n",
            string(name), length(p.delays), p.powers_db[imax], p.delays[imax])
end

println()
println("== CDL-A at fs = 100 MHz, nominal delay spread 100 ns (TR 38.901 Table 7.7.3-1) ==")
taps = cdl_taps(cdl_profile(:A), 1, 1e-8; delay_spread=100e-9)
println(" 23 clusters merge to ", length(taps.delays), " taps (powers sum to ",
        sum(taps.powers), ")")
println("  tap   delay [samples]   delay [ns]   power [dB]")
for (i, (d, p)) in enumerate(zip(taps.delays, taps.powers))
    @printf("  %3d   %13d   %10.1f   %10.2f\n", i, d, d * 10.0, 10 * log10(p))
end

println()
println("== Static drop through the textbook OFDM chain (known CSI, QPSK) ==")
params = OFDMParams(N_fft=32, N_cp=8, mod_order=4, n_guard=2, pilot_spacing=4)
bits = rand(0:1, params.N_data * params.bits_per_symbol * 40)
tx, pilots = ofdm_tx(params, bits)
drop = cdl_taps(cdl_profile(:A), 1, 40e-9; delay_spread=30e-9)
println("  CDL-A, fs = 25 MHz, DS = 30 ns: taps at samples ", drop.delays,
        " (N_cp = ", params.N_cp, ")")
tdl = TDLChannel(drop.gains[:, 1], drop.delays)
flat = TDLChannel(delays_samples=[0], gain_dB=[0.0])
for snr_db in (10.0, 20.0)
    rx_cdl, H_cdl = apply_channel(tx, tdl, snr_db, params)
    rx_flat, H_flat = apply_channel(tx, flat, snr_db, params)
    e_cdl = count(!=(0), bits .- ofdm_rx_known(rx_cdl, params, H_cdl, pilots)[1:length(bits)])
    e_flat = count(!=(0), bits .- ofdm_rx_known(rx_flat, params, H_flat, pilots)[1:length(bits)])
    @printf("  %4.1f dB: CDL static drop %4d errors,  flat reference %4d errors  (%d bits, simulation only)\n",
            snr_db, e_cdl, e_flat, length(bits))
end

println()
println("== Per-tap Clarke Doppler (fd = 2 MHz at fs = 100 MHz, fd*Ts = 0.02) ==")
tv = cdl_taps(cdl_profile(:A), 4096, 1e-8; delay_spread=100e-9, fd=2e6)
i_strong = argmax(tv.powers)
@printf("  strongest tap (delay %d samples, power %.2f dB): |g| over 6 samples = ",
        tv.delays[i_strong], 10 * log10(tv.powers[i_strong]))
foreach(g -> @printf("%.3f ", abs(g)), tv.gains[i_strong, 1:6])
println()
m = sum(abs2, tv.gains[i_strong, :]) / size(tv.gains, 2)
@printf("  time-averaged power of that tap %.3f vs tabulated %.3f (one record)\n",
        m, tv.powers[i_strong])

println()
println("== Flat degenerate: delay spread far below one sample is flat Rayleigh ==")
function bpsk_rayleigh_formula(γ)
    return 0.5 * (1 - sqrt(γ / (1 + γ)))
end
snr_db = 8.0
γ = 10^(snr_db / 10)
nbits = 40_000
errors = Ref(0)
for _ in 1:nbits
    bit = rand(0:1)
    rx, t = cdl_channel(modulate([bit], BPSK()), snr_db, 1;
                        Ts=1e-8, delay_spread=1e-12)
    errors[] += demodulate([rx[1] / t.gains[1, 1]], BPSK())[1] != bit
end
@printf("  Eb/N0 = %.1f dB, %d bits: sim %.4e,  flat-Rayleigh formula %.4e\n",
        snr_db, nbits, errors[] / nbits, bpsk_rayleigh_formula(γ))
println("  (Only the fully merged case has a closed form; dispersive CDL has no golden curve.)")
