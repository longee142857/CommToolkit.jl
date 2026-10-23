# BICM versus BICM-ID on a terminated (7, 5) code and 8-PSK set partitioning.
# Not Gray labelling, not 5G NR, not an EXIT chart, not a capacity claim.
#
# Es = 1. Terminated rate R = 400/804. Eb = 1/(R * 3) = 0.67.
# σ² = Eb / (2 Eb/N0), the same split as awgn. Max-log demapper.
# The bit and noise stream is xorshift + Box–Muller, not Julia's global RNG.
# Run from the repo root:
#   julia --project=. demos/ber_bicm_id.jl

using Printf
using CommToolkit

const N_INFO = 400
const N_ROWS = 12
const N_CODED = 2 * (N_INFO + 2)
const EB = 1 / ((N_INFO / N_CODED) * 3)
const LABELS = bicm_psk8_id_labels()
const SEED = 0x0c0de1162

function _xorshift(state::UInt64)
    z = xor(state, state >> 12)
    z = xor(z, z << 25)
    return xor(z, z >> 27)
end

function _bits(k::Int, seed::UInt64)
    st = seed
    bits = Vector{Int}(undef, k)
    for i in 1:k
        st = _xorshift(st)
        bits[i] = Int((st * 0x2545f4914f6cdd1d) >> 63)
    end
    return bits
end

function _normals(n::Int, seed::UInt64)
    st = Ref(seed)
    function uni()
        st[] = _xorshift(st[])
        u = st[] * 0x2545f4914f6cdd1d
        return (u >> 11) * 2.0^-53
    end
    out = Vector{Float64}(undef, n)
    i = 1
    while i <= n
        u1 = max(uni(), 1e-12)
        r = sqrt(-2 * log(u1))
        θ = 2π * uni()
        out[i] = r * cos(θ)
        i += 1
        if i <= n
            out[i] = r * sin(θ)
            i += 1
        end
    end
    return out
end

function run_point(snr_db::Float64, n_frames::Int; fading::Bool=false)
    σ = sqrt(EB / (2 * 10^(snr_db / 10)))
    bit_err = zeros(Int, 2)
    frame_err = zeros(Int, 2)
    for f in 1:n_frames
        bits = _bits(N_INFO, SEED + UInt64(f))
        syms = bicm_map(bits, PSK(8); n_rows=N_ROWS, labels=LABELS)
        nsym = length(syms)
        if fading
            z = _normals(4 * nsym, SEED + 0x10000 + UInt64(f))
            h = (z[1:2:(2nsym)] .+ im .* z[2:2:(2nsym)]) ./ sqrt(2)
            noise = σ .* (z[(2nsym + 1):2:end] .+ im .* z[(2nsym + 2):2:end])
            y = h .* syms .+ noise
            gains = h
        else
            z = _normals(2 * nsym, SEED + 0x10000 + UInt64(f))
            y = syms .+ σ .* (z[1:2:end] .+ im .* z[2:2:end])
            gains = nothing
        end
        for (j, it) in enumerate((1, 4))
            dec = bicm_decode(y, PSK(8), σ^2, N_INFO; iterations=it, n_rows=N_ROWS,
                              labels=LABELS, approx=:maxlog, gains=gains)
            e = count(!=(0), dec.bits .- bits)
            bit_err[j] += e
            e > 0 && (frame_err[j] += 1)
        end
    end
    nbits = n_frames * N_INFO
    channel = fading ? "Rayleigh" : "AWGN    "
    @printf("%s  %4.1f dB   %8.5f  %8.5f   %6.3f  %6.3f\n",
            channel, snr_db, bit_err[1] / nbits, bit_err[2] / nbits,
            frame_err[1] / n_frames, frame_err[2] / n_frames)
    return nothing
end

println("8-PSK set partitioning, (7, 5), 400 info bits, 12×67 interleaver, max-log")
println("Eb = $(EB) (terminated rate, Es = 1). Not Gray. Not NR.")
@printf("%-8s  %6s   %8s  %8s   %6s  %6s\n", "channel", "Eb/N0", "BER 1", "BER 4", "BLER 1", "BLER 4")
run_point(4.0, 30)
run_point(6.0, 40)
run_point(8.0, 30)
run_point(8.0, 20; fading=true)
