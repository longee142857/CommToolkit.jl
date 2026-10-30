# Teaching soft-output Viterbi on the terminated K = 3, (7, 5) code.
# Information-bit extrinsic LLRs are the default. return_coded=true also
# returns coded-bit extrinsic LLRs: a Hagenauer–Hoeher competing-path gap
# with the channel LLR removed. This is not the BCJR, not MAP, and not a
# turbo decoder. At the same Eb/N0 the information-bit error rate or the
# soft mutual information can be worse than bcjr_decode. A metric tie
# (Δ = 0) cancels the prior. A bit that never splits is capped at
# Σ|L| + Σ|La|, not ±∞.
#
# LLRs are positive towards bit 0, the demodulate_soft convention.
# Es = 1 on coded BPSK, Eb/N0 via energy_per_bit = 2 (nominal rate 1/2).
# Seed 0x0c0de1162. Run from the repo root:
#   julia --project=. demos/sova_soft.jl

using Random
using Printf
using CommToolkit

const SEED = 0x0c0de1162

println("CommToolkit teaching SOVA (not MAP, not BCJR, not a turbo decoder).")
println("Seed ", string(SEED, base=16), ". Terminated K=3 (7, 5). Positive LLR favours bit 0.")
println()

bits = [1, 0, 1, 1, 0]
coded = conv_encode(bits)
strong = 30.0 .* (1.0 .- 2.0 .* coded)
info, Lc = sova_decode(strong; return_coded=true)
println("Noiseless roundtrip, |L| = 30:")
println("  information bits     ", bits)
println("  SOVA hard            ", Int.(info .< 0))
println("  Viterbi hard         ", viterbi_decode(strong))
println("  BCJR hard            ", Int.(bcjr_decode(strong) .< 0))
println("  coded bits           ", coded)
println("  coded extrinsic hard ", Int.(Lc .< 0))
println("  |info LLR|           ", round.(abs.(info); digits=1))
println("  |coded extrinsic|    ", round.(abs.(Lc); digits=1), "  (gap minus |L|)")
println()

one = conv_encode([1])
Lone = 12.0 .* (1.0 .- 2.0 .* one)
_, Lc1 = sova_decode(Lone; return_coded=true)
ρ = abs.(Lc1 .+ Lone)
println("Truncation on a 1-bit word, A = 12, cap Σ|L| = ", sum(abs, Lone), ":")
println("  coded reliability ρ  ", ρ)
println("  (5A = 60 is a measured gap; the cap is the bit that never split)")
println()

zinfo, zc = sova_decode(zeros(2 * (4 + 2)); return_coded=true)
println("Tie, channel LLR = 0: info ", zinfo, " coded max |L| ", maximum(abs, zc))
println("  Equal metrics, Δ = 0. Survivor prefers the smaller input bit.")
println()

function sweep(snr_db, frames, ninfo)
    γ = 10^(snr_db / 10)
    σ2 = 1 / γ
    es = eb = 0
    mis = mib = 0.0
    Random.seed!(SEED)
    for _ in 1:frames
        msg = rand(0:1, ninfo)
        word = conv_encode(msg)
        y = awgn(modulate(word, BPSK()), snr_db; energy_per_bit=2.0)
        L = demodulate_soft(y, BPSK(), σ2)
        Ls = sova_decode(L)
        Lb = bcjr_decode(L)
        es += count(!=(0), Int.(Ls .< 0) .- msg)
        eb += count(!=(0), Int.(Lb .< 0) .- msg)
        mis += exit_mi_soft(Ls, msg)
        mib += exit_mi_soft(Lb, msg)
    end
    nbits = frames * ninfo
    return es / nbits, eb / nbits, mis / frames, mib / frames, es, eb, nbits
end

println("Information-bit comparison. Same frames, same noise. SOVA is allowed to lose.")
@printf("%8s  %10s  %10s  %10s  %10s  %s\n",
        "Eb/N0", "BER_SOVA", "BER_BCJR", "MI_SOVA", "MI_BCJR", "errors sova/bcjr")
for (snr, frames) in ((0.0, 48), (2.0, 40), (6.0, 40))
    ps, pb, mis, mib, es, eb, nbits = sweep(snr, frames, 64)
    @printf("%8.1f  %10.5f  %10.5f  %10.4f  %10.4f  %d/%d of %d\n",
            snr, ps, pb, mis, mib, es, eb, nbits)
end
println()
println("0 dB: SOVA BER is higher and soft MI is lower (over-confident path gaps).")
println("2 dB: both beat uncoded BPSK; the ranking of the two point estimates is not the claim.")
println("6 dB: hard decisions match. Magnitudes are still not required to equal BCJR.")
println("Coded extrinsic is L_e(c) = (1-2ĉ)ρ - L(c). Slicing L_e(c) alone is not the ML bit")
println("when ρ < |L(c)|. BICM and exit_decoder_transfer do not consume this vector.")
