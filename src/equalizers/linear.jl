# Linear equalizers: ZF and MMSE
#
# For frequency-selective channels with impulse response h
# Channel matrix H is Toeplitz (convolution matrix)

using LinearAlgebra

# Build convolution matrix from channel impulse response
function _channel_matrix(h::AbstractVector{<:Number}, n::Int)
    L = length(h)
    H = zeros(promote_type(eltype(h), Float64), n + L - 1, n)
    for i in 1:n
        for j in 1:L
            H[i + j - 1, i] = h[j]
        end
    end
    return H
end

# Zero-Forcing equalizer: W = (H^H * H)^{-1} * H^H
function zf_equalizer(h::AbstractVector{<:Number}, n::Int)
    H = _channel_matrix(h, n)
    return (H' * H) \ H'
end

function zf_equalize(syms::AbstractVector{<:Number},
                     h::AbstractVector{<:Number})
    n = length(syms) - length(h) + 1
    @assert n >= 1 "received signal too short for given channel"
    W = zf_equalizer(h, n)
    return W * syms
end

# MMSE equalizer: W = (H^H * H + σ^2 * I)^{-1} * H^H
function mmse_equalizer(h::AbstractVector{<:Number}, n::Int, noise_var::Real)
    H = _channel_matrix(h, n)
    return (H' * H + noise_var * I) \ H'
end

function mmse_equalize(syms::AbstractVector{<:Number},
                       h::AbstractVector{<:Number},
                       noise_var::Real)
    n = length(syms) - length(h) + 1
    @assert n >= 1 "received signal too short for given channel"
    W = mmse_equalizer(h, n, noise_var)
    return W * syms
end
