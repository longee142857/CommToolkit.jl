# Teaching Arikan polar code: Bhattacharyya construction, non-systematic and
# systematic encoding, and successive-cancellation (SC) decoding. Short frames
# only: N = 2^n, 2 ≤ N ≤ 256. Not a 3GPP NR polar code — no CRC attachment,
# no CRC-aided SCL list decoding, no rate matching, no NR reliability
# sequence. The information set comes from the Bhattacharyya recursion
# Z⁻ = 2Z − Z², Z⁺ = Z² (Arikan 2009, Prop. 5), exact on the BEC and a
# teaching heuristic on AWGN.
#
# LLRs are positive towards bit 0, the demodulate_soft convention. Seed
# 0x0c0de1162. Run from the repo root:
#   julia --project=. demos/polar_sc.jl

using Random
using Printf
using CommToolkit

const SEED = 0x0c0de1162

# Q(√(2γ)) for the theory columns, without loading SpecialFunctions:
# Abramowitz–Stegun 7.1.26 erf approximation (|error| ≤ 1.5e-7), enough for
# a printed comparison. The tested anchors live in test/reference.jl.
function _demo_q_bpsk(γ::Real)
    x = sqrt(γ) # erf(x) = 2Φ(x√2) − 1, so Q(√(2γ)) = (1 − erf(√γ)) / 2
    t = 1 / (1 + 0.3275911 * x)
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 +
           t * (-1.453152027 + t * 1.061405429))))
    return 0.5 * (1 - (1 - poly * exp(-x * x)))
end

println("CommToolkit teaching Arikan polar SC (not a 3GPP NR polar code).")
println("Seed ", string(SEED, base=16), ". N = 2^n ≤ 256, Bhattacharyya construction, SC decoding.")
println()

code = polar_code(64, 32; design_snr_db=2.0)
println("polar_code(64, 32; design_snr_db=2.0)")
println("  info_pos = ", code.info_pos)
println("  worst info z = ", round(maximum(code.z[code.info_pos]); sigdigits=4),
        ", best frozen z = ", round(minimum(code.z[code.frozen]); sigdigits=4))
println()

Random.seed!(SEED)
msg = rand(0:1, code.k)
x_ns = polar_encode(code, msg)
x_sy = polar_encode(code, msg; systematic=true)
println("Noiseless roundtrip on a random 32-bit message:")
println("  msg              = ", msg)
println("  non-systematic x = ", x_ns)
println("  systematic x     = ", x_sy, "  (carries msg on info_pos)")
println("  SC decode        = ", polar_sc_decode(code, 30.0 .* (1 .- 2 .* x_ns)), "  (non-systematic)")
println("  SC decode        = ", polar_sc_decode(code, 30.0 .* (1 .- 2 .* x_sy); systematic=true), "  (systematic)")
println()

# k = 1 anchor: the (N, 1) polar code is the length-N repetition code, SC
# sums the channel LLRs and is exactly ML, so BER = Q(√(2 γb)) — a
# repetition code has no AWGN coding gain, which makes it a closed-form
# anchor for the decoder and the energy bookkeeping.
println("Repetition anchor (k = 1): SC is ML, theory Q(√(2γb)) = uncoded BPSK BER.")
@printf("%4s  %8s  %8s  %12s  %12s\n", "N", "Eb/N0_dB", "BER_sim", "Q(√(2γb))", "errors/bits")
for (n, snr_db, frames) in ((2, 2.0, 8000), (4, 3.0, 12000))
    rep = polar_code(n, 1)
    γ = 10^(snr_db / 10)
    σ2 = n / (2 * γ)
    bit_errors = 0
    Random.seed!(SEED)
    for _ in 1:frames
        m = rand(0:1, 1)
        y = awgn(modulate(polar_encode(rep, m), BPSK()), snr_db; energy_per_bit=n)
        bit_errors += count(!=(0), m .- polar_sc_decode(rep, demodulate_soft(y, BPSK(), σ2)))
    end
    theory = _demo_q_bpsk(γ)
    @printf("%4d  %8.1f  %8.5f  %12.5f  %12s\n", n, snr_db,
            bit_errors / frames, theory, "$(bit_errors)/$(frames)")
end
println()

# Small BER sweep for the (64, 32) code against uncoded BPSK at the same
# information-bit energy. Es = 1, Eb = n/k, σ² = Eb/(2 Eb/N0).
println("(64, 32) SC versus uncoded BPSK, 1000 frames per point:")
@printf("%8s  %10s  %14s  %12s  %10s\n", "Eb/N0_dB", "BER_sim", "uncoded_BER", "BLER_sim", "frames")
for snr_db in (2.0, 3.0)
    γ = 10^(snr_db / 10)
    σ2 = 2.0 / (2 * γ)
    bit_errors = 0
    block_errors = 0
    Random.seed!(SEED)
    for _ in 1:1000
        bits = rand(0:1, code.k)
        y = awgn(modulate(polar_encode(code, bits), BPSK()), snr_db; energy_per_bit=2.0)
        errors = count(!=(0), bits .- polar_sc_decode(code, demodulate_soft(y, BPSK(), σ2)))
        bit_errors += errors
        block_errors += errors > 0
    end
    @printf("%8.1f  %10.5f  %14.5f  %12.5f  %10s\n", snr_db,
            bit_errors / (1000 * code.k), _demo_q_bpsk(γ),
            block_errors / 1000, "$(block_errors)/1000")
end
println()
println("Boundaries: no CRC, no CA-SCL, no rate matching, no PW/QPP interleaver,")
println("no 3GPP mother code, no NR reliability sequence. NR polar stays outside.")
