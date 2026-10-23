# Gardner open-loop timing on a root-raised-cosine burst, a Costas loop on
# QPSK with a small phase and frequency offset, and a Schmidl–Cox preamble
# with an even-integer CFO. No Plots. Not a production receiver.

using CommToolkit

sps = 4
taps = rrc_filter(0.35, 10, sps)
bits = [((i * 17 + 3) % 7) > 3 ? 1 : 0 for i in 1:160]
symbols = modulate(bits, BPSK())
timed = single_carrier_timing(symbols, taps, sps; n_frac=8)
window = 30:(length(symbols) - 30)
idx = timed.symbol_index .+ (window .- 1)
timing_err = maximum(abs, timed.symbols[idx] .- symbols[window])
shifted = single_carrier_timing(symbols, taps, sps; sample_delay=0.4, n_frac=8)
println("Gardner delay at the matched-filter peak (samples, mod sps) = ", timed.delay)
println("Gardner delay after an extra 0.4 sample = ", shifted.delay)
println("symbol 1 index = ", timed.symbol_index)
println("max |timing error| on the interior = ", timing_err)

qbits = [((i * 17 + 3) % 7) > 3 ? 1 : 0 for i in 1:600]
qpsk = modulate(qbits, QPSK())
θ = 0.2
ω = 0.01
spun = qpsk .* exp.(im .* (θ .+ ω .* (0:(length(qpsk) - 1))))
loop = costas_loop(spun, QPSK())
tail = loop.symbols[(end - 40):end]
tail_bits = qbits[(end - 2 * length(tail) + 1):end]
println("Costas frequency, last 80 mean = ", sum(loop.frequency[(end - 79):end]) / 80)
println("Costas tail bit errors = ", count(!=(0), demodulate(tail, QPSK()) .- tail_bits))

pre = schmidl_preamble(64; n_cp=16)
ε = 2.25
N = length(pre.freq)
rx = [pre.time[i] * exp(im * 2π * ε * (i - 1) / N) for i in eachindex(pre.time)]
sync = schmidl_sync(rx, pre.freq; n_cp=pre.n_cp)
println("Schmidl d = ", sync.d, "  eps_f = ", sync.epsilon_f, "  eps_i = ", sync.epsilon_i)
println("total eps = ", sync.epsilon)

timing_err < 0.05 || error("Gardner interior error is too large")
abs(sync.epsilon - ε) < 1e-6 || error("Schmidl–Cox did not recover the CFO")
sync.d == 1 || error("Schmidl–Cox frame index moved")
