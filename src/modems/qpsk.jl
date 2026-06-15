# QPSK modulation
#
# Maps: bit pair (b0, b1) → complex symbol
# Gray coding:
#   00 → +1+j, 01 → +1-j, 10 → -1+j, 11 → -1-j
# Normalized to unit average power

const QPSK_CONSTELLATION = [ 1+1im,  1-1im, -1+1im, -1-1im] ./ sqrt(2)

function modulate(bits::AbstractVector{<:Integer}, ::QPSK)
    @assert all(b -> b ∈ (0, 1), bits) "bits must be 0 or 1"
    @assert length(bits) % 2 == 0 "bit length must be even for QPSK"

    n_sym = length(bits) ÷ 2
    syms = Vector{ComplexF64}(undef, n_sym)
    for i in 1:n_sym
        b0, b1 = bits[2i - 1], bits[2i]
        idx = 2 * b0 + b1 + 1  # gray mapping: 00→1, 01→2, 10→3, 11→4
        # Corrected Gray mapping:
        # 00 (b0=0,b1=0) → 1+1j  (idx=1)
        # 01 (b0=0,b1=1) → 1-1j  (idx=2)
        # 10 (b0=1,b1=0) → -1+1j (idx=3)
        # 11 (b0=1,b1=1) → -1-1j (idx=4)
        syms[i] = QPSK_CONSTELLATION[idx]
    end
    return syms
end

function demodulate(syms::AbstractVector{<:Complex}, ::QPSK)
    n_sym = length(syms)
    bits = Vector{Int}(undef, 2 * n_sym)
    for i in 1:n_sym
        # Hard decision by quadrant
        bits[2i - 1] = Int(real(syms[i]) < 0)
        bits[2i]     = Int(imag(syms[i]) < 0)
    end
    return bits
end
