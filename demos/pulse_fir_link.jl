# pulse_shape → sample-rate FIR → matched_filter → downsample → ZF / MMSE.
# The pulse is a single sample, so the symbol-rate channel is the FIR
# decimated on that phase and the noiseless equaliser recovers the symbols.
# No Plots.

using CommToolkit

sps = 2
pulse = [0.0, 1.0]
h_sample = [0.8, 0.0, -0.3]
h_symbol = [0.8, -0.3]
bits = [0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0]
symbols = modulate(bits, QPSK())

waveform = pulse_shape(symbols, pulse, sps)
faded = apply_fir(waveform, h_sample)
filtered = matched_filter(faded, pulse)
sampled = downsample_symbols(filtered, sps; offset=1)

zf = zf_equalize(sampled, h_symbol)
mmse = mmse_equalize(sampled, h_symbol, 0.0)

zf_err = maximum(abs, zf .- symbols)
mmse_err = maximum(abs, mmse .- symbols)
println("symbols = ", length(symbols))
println("max |ZF error| = ", zf_err)
println("max |MMSE error| = ", mmse_err)
println("ZF bits match = ", demodulate(zf, QPSK()) == bits)
println("MMSE bits match = ", demodulate(mmse, QPSK()) == bits)

zf_err < 1e-8 || error("ZF did not recover the symbols")
mmse_err < 1e-8 || error("MMSE did not recover the symbols")
