module Equalization

using LinearAlgebra

export zf_equalizer, zf_equalize, mmse_equalizer, mmse_equalize

@doc raw"""
    _channel_matrix(h, n) -> Matrix

Tall Toeplitz convolution matrix ``H`` of size ``(n + L - 1) \times n`` for an impulse response `h` of length ``L``. Column ``i`` is `h` placed at row ``i``. The model is ``r = H s + \nu``.
"""
function _channel_matrix(h::AbstractVector{<:Number}, n::Int)
    n >= 1 || throw(ArgumentError("n must be ≥ 1"))
    L = length(h)
    L >= 1 || throw(ArgumentError("channel impulse response must be non-empty"))
    H = zeros(promote_type(eltype(h), Float64), n + L - 1, n)
    @inbounds for i in 1:n
        for j in 1:L
            H[i + j - 1, i] = h[j]
        end
    end
    return H
end

@doc raw"""
    zf_equalizer(h::AbstractVector{<:Number}, n::Integer) -> Matrix

Zero-forcing block equalizer for a length-`n` symbol block observed through the convolution matrix ``H`` of `h`,

```math
W_{\mathrm{ZF}} = (H^{\mathrm{H}} H)^{-1} H^{\mathrm{H}}.
```

`n` is the number of unknown symbols, not the length of the received vector. ``H`` has full column rank when `h` is not identically zero, so the inverse exists (Proakis and Salehi, *Digital Communications*, 5th ed., §9.4; Tse and Viswanath, *Fundamentals of Wireless Communication*, Cambridge, 2005, §8.3).

# Returns
A matrix of size ``n \times (n + L - 1)``. Apply it with `W * r`.
"""
function zf_equalizer(h::AbstractVector{<:Number}, n::Integer)
    H = _channel_matrix(h, Int(n))
    return (H' * H) \ H'
end

@doc raw"""
    zf_equalize(received::AbstractVector{<:Number}, h::AbstractVector{<:Number})

Zero-forcing equalisation of one convolution. `received` has length ``n + L - 1`` and `h` has length ``L``; the result has length ``n``.

Noiseless convolution is inverted exactly (up to linear-algebra roundoff) when ``H`` has full column rank. See [`zf_equalizer`](@ref).
"""
function zf_equalize(syms::AbstractVector{<:Number}, h::AbstractVector{<:Number})
    n = length(syms) - length(h) + 1
    n >= 1 || throw(ArgumentError("received signal is too short for the given channel"))
    return zf_equalizer(h, n) * syms
end

@doc raw"""
    mmse_equalizer(h::AbstractVector{<:Number}, n::Integer, noise_var::Real) -> Matrix

Linear minimum-mean-square-error block equalizer

```math
W_{\mathrm{MMSE}} = \left(H^{\mathrm{H}} H + \sigma^2 I\right)^{-1} H^{\mathrm{H}},
```

where ``\sigma^2`` is `noise_var`, the variance of each entry of the additive noise in ``r = H s + \nu``, and the unknown symbols are modelled as unit-variance and uncorrelated (Proakis and Salehi, 5th ed., §9.4; Tse and Viswanath, 2005, §8.3).

`noise_var = 0` reduces to [`zf_equalizer`](@ref). `noise_var` must be nonnegative.
"""
function mmse_equalizer(h::AbstractVector{<:Number}, n::Integer, noise_var::Real)
    noise_var >= 0 || throw(ArgumentError("noise_var must be ≥ 0, got $noise_var"))
    H = _channel_matrix(h, Int(n))
    return (H' * H + float(noise_var) * I) \ H'
end

@doc raw"""
    mmse_equalize(received, h, noise_var::Real)

Apply [`mmse_equalizer`](@ref) to one received block. `received` has length ``n + L - 1``. `noise_var` is the per-sample noise variance in that same vector, not a decibel SNR.
"""
function mmse_equalize(syms::AbstractVector{<:Number}, h::AbstractVector{<:Number}, noise_var::Real)
    n = length(syms) - length(h) + 1
    n >= 1 || throw(ArgumentError("received signal is too short for the given channel"))
    return mmse_equalizer(h, n, noise_var) * syms
end

end # module Equalization
