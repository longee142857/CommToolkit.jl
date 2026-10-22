# Terminated (7, 5) BPSK on block Rayleigh fading, with and without a
# rectangular interleaver. The fade is constant for `block_length` symbols
# and independent across blocks. This is not a Jakes model.
#
# Es = 1 on coded BPSK, Eb = 2 (nominal rate 1/2), σ² = N0/2.
# Coherent LLR: L = 2 Re(h* y) / σ². Seed 0x0c0de1162.
# Run from the repo root:
#   julia --project=. demos/ber_block_interleave.jl

using Random
using Printf
using CommToolkit

const SEED = 0x0c0de1162
const N_ROWS = 16
const N_COLS = 40
const BLOCK_LENGTH = 16
const N_FRAMES = 30
const EB = 2.0

n_coded = N_ROWS * N_COLS
n_info = n_coded ÷ 2 - 2

function bpsk_csi_llr(y, h, σ2)
    llr = Vector{Float64}(undef, length(y))
    @inbounds for i in eachindex(y)
        llr[i] = 2 * real(conj(h[i]) * y[i]) / σ2
    end
    return llr
end

function run_snr(snr_db::Float64)
    σ2 = EB / (2 * snr_db_to_linear(snr_db))
    err_plain = 0
    err_intl = 0
    frames_plain = 0
    frames_intl = 0
    Random.seed!(SEED)
    for _ in 1:N_FRAMES
        bits = rand(0:1, n_info)
        coded = conv_encode(bits)
        syms = ComplexF64.(modulate(coded, BPSK()))
        syms_i = block_interleave(syms, N_ROWS, N_COLS)
        h = block_rayleigh(n_coded, BLOCK_LENGTH)
        noise = awgn(zeros(ComplexF64, n_coded), snr_db; energy_per_bit=EB)
        dec = viterbi_decode(bpsk_csi_llr(h .* syms .+ noise, h, σ2))
        dec_i = viterbi_decode(block_deinterleave(
            bpsk_csi_llr(h .* syms_i .+ noise, h, σ2),
            N_ROWS, N_COLS,
        ))
        e = count(!=(0), dec .- bits)
        ei = count(!=(0), dec_i .- bits)
        err_plain += e
        err_intl += ei
        frames_plain += e > 0
        frames_intl += ei > 0
    end
    nbits = n_info * N_FRAMES
    return err_plain / nbits, err_intl / nbits,
           frames_plain / N_FRAMES, frames_intl / N_FRAMES
end

println("Block fading + (7, 5), block_length = ", BLOCK_LENGTH,
        ", interleaver ", N_ROWS, "×", N_COLS,
        ", ", N_FRAMES, " frames, seed ", string(SEED, base=16))
println("Eb/N0   BER plain     BER interleaved    BLER plain   BLER interleaved")
for snr_db in (8.0, 12.0)
    bp, bi, fp, fi = run_snr(snr_db)
    @printf("%5.1f   %.3e      %.3e           %.3f        %.3f\n", snr_db, bp, bi, fp, fi)
end
