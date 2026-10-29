# Short Reed–Solomon: RS(15, 11) over GF(2^4) and its shortening RS(10, 6).
# Narrow-sense primitive teaching codes with hard-decision Berlekamp–Massey
# decoding. Not CCSDS (255, 223) / (255, 239), not a soft or list decoder.
#
# Symbols are integers 0..15 (bitmasks over GF(2^4), primitive polynomial
# x^4+x+1). Codewords are highest degree first, systematic part in front.
# Seed 0x0c0de1162. Run from the repo root:
#   julia --project=. demos/rs_short.jl

using Random
using Printf
using CommToolkit

const SEED = 0x0c0de1162

println("CommToolkit short Reed–Solomon (teaching, not CCSDS).")
println("Seed ", string(SEED, base=16), ". Symbols are GF(2^4) integers, α = x.")
println()

C = rs_code(4, 11)
msg = collect(1:11)
cw = rs_encode(C, msg)
println("RS(15, 11), t = ", C.t, ", g(x) = x^4 + α^13 x^3 + α^6 x^2 + α^3 x + α^10")
println("  msg    = ", msg)
println("  codew  = ", cw)
rx = copy(cw)
rx[2] ⊻= 0x9
rx[13] ⊻= 0x5
println("  rx     = ", rx, "  (two symbol errors)")
println("  decode = ", rs_decode(C, rx), "  (recovered)")

println()
Cs = rs_shorten(C, 5)
smsg = collect(1:6)
scw = rs_encode(Cs, smsg)
println("RS(10, 6) shortened from RS(15, 11), t = ", Cs.t)
println("  msg    = ", smsg)
println("  codew  = ", scw)
srx = copy(scw)
srx[4] ⊻= 0x7
srx[9] ⊻= 0x2
println("  rx     = ", srx, "  (two symbol errors)")
println("  decode = ", rs_decode(Cs, srx), "  (recovered)")

println()
println("Beyond t the decoder has no failure flag: it returns the systematic")
println("symbols unchanged, or miscorrects to a codeword in the radius-t ball.")
println("Descriptive sweep on RS(15, 11): uniform symbol errors at rate p,")
println("4000 frames per point. tail = P(at most t errors), the guaranteed")
println("part of the correct fraction (errors only in parity can add to it).")
println()
@printf("%6s  %8s  %10s  %10s\n", "p", "correct", "unchanged", "tail")
Random.seed!(SEED)
for p in (0.02, 0.05, 0.10, 0.20)
    n_correct = 0
    n_unchanged = 0
    n_frames = 4000
    for _ in 1:n_frames
        m = rand(0:15, 11)
        r = copy(rs_encode(C, m))
        for i in 1:15
            rand() < p && (r[i] ⊻= rand(1:15))
        end
        d = rs_decode(C, r)
        d == m && (n_correct += 1)
        d == r[1:11] && (n_unchanged += 1)
    end
    tail = sum(binomial(15, i) * p^i * (1 - p)^(15 - i) for i in 0:C.t)
    @printf("%6.2f  %8.4f  %10.4f  %10.4f\n", p, n_correct / n_frames,
            n_unchanged / n_frames, tail)
end
