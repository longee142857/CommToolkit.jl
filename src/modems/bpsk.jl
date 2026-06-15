# BPSK modulation
#
# Maps: bit 0 → +1, bit 1 → -1
# Constellation: {+1, -1}

function modulate(bits::AbstractVector{<:Integer}, ::BPSK)
    @assert all(b -> b ∈ (0, 1), bits) "bits must be 0 or 1"
    return 1 .- 2 .* bits  # 0→+1, 1→-1
end

function demodulate(syms::AbstractVector{<:Real}, ::BPSK)
    return Int.(syms .< 0)  # negative → 1, positive → 0
end

function demodulate(syms::AbstractVector{<:Complex}, ::BPSK)
    return Int.(real(syms) .< 0)
end

# Soft demodulation: return log-likelihood ratio
function demodulate_soft(syms::AbstractVector{<:Real}, ::BPSK, noise_var::Real)
    return 2 .* syms ./ noise_var
end
