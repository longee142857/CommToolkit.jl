module Equalization

using LinearAlgebra

export zf_equalizer, zf_equalize, mmse_equalizer, mmse_equalize,
    lms_equalize, DFECoefficients, dfe_equalizer, dfe_equalize, zf_dfe, mmse_dfe

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

Noiseless convolution is inverted exactly (up to linear-algebra roundoff) when ``H`` has full column rank. [`CommToolkit.Channel.apply_fir`](@ref)`(s, h)` is that convolution and has length ``n + L - 1``. See [`zf_equalizer`](@ref).
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

@doc raw"""
    lms_equalize(received, training, n_taps, mu; delay=0, n_symbols=length(training), w0=nothing, tap_trajectory=false)

Symbol-spaced complex LMS equalizer (Proakis and Salehi, *Digital Communications*, 5th ed., §10.1-2). Real inputs stay on the real algorithm: conjugation is the identity.

At symbol index ``m``, the regressor ``u`` holds `n_taps` received samples, newest first,

```math
u_i = r[m + \delta - i + 1], \qquad i = 1,\ldots,N,
```

with ``\delta`` equal to `delay` and ``r[k] = 0`` outside `received`. The filter and the update are

```math
\hat{y} = w^{\mathsf{T}} u, \qquad e = d - \hat{y}, \qquad w \leftarrow w + \mu\, e\, u^{*}.
```

``d`` is the training symbol. ``\hat{y}`` is computed with the weights from **before** this update, and that value is what the function returns. Column ``t`` of the trajectory, when requested, is ``w`` **after** the supervised update of symbol ``t``.

# Arguments
- `received`: symbol-rate observation. [`CommToolkit.Channel.apply_fir`](@ref)`(symbols, h)` is the matching convolution when `symbols` is the transmitted block. A sample-rate waveform that has not been reduced to one sample per symbol is the wrong input.
- `training`: known symbols in transmission order, starting at symbol 1. Must be non-empty. Updates run for symbols `1:min(length(training), n_symbols)` and then stop.
- `n_taps`: length of ``w``. At least 1.
- `mu`: step size. Finite and strictly positive. The caller chooses it. A step larger than about ``2/\lambda_{\max}(R)`` diverges; this function does not estimate ``\lambda_{\max}`` or the received power, and it does not implement normalized LMS.
- `delay`: nonnegative integer ``\delta``. `delay = 0` puts the newest sample on the same index as the symbol, which matches the first tap of `h` in `apply_fir`. A channel whose energy sits on a later tap needs a larger delay. The delay is not searched.
- `n_symbols`: number of estimates to return. The default is `length(training)`, so the default result is only the acquisition output. A larger value keeps filtering symbols after the training prefix with the **frozen** weights from the last supervised update. Requires `n_symbols + delay ≤ length(received)`.
- `w0`: initial weights, length `n_taps`. The default is zeros, so the first outputs are the acquisition transient.
- `tap_trajectory`: when `true`, also return the `n_taps × n_supervised` weight matrix.

# Returns
With `tap_trajectory = false`, a vector of length `n_symbols`. With `tap_trajectory = true`, `(symbols, trajectory)`.

The element type is `promote_type` of `received`, `training`, `w0` when given, and `Float64`.

# Steady-state error
Settled mean-square error grows with ``\mu`` on a fixed record, once the step is small enough to converge. That direction is checked in `test/lms_dfe.jl`. The classical excess-MSE formula (Proakis and Salehi §10.1-4) is not evaluated and is not returned: it needs a small step, the independence assumption, and a known minimum MSE. With no noise, and with enough taps at a delay that covers the channel cursor, the error can fall to roundoff. Nothing here estimates the misadjustment.

# Limits
Supervised only. After the training prefix the weights freeze; there is no decision-directed mode and no blind (CMA) mode. Not RLS, not a Kalman filter, and not a fractionally spaced equalizer. The training-segment outputs are not the block refiltered by the final weights.
"""
function lms_equalize(
    received::AbstractVector{<:Number},
    training::AbstractVector{<:Number},
    n_taps::Integer,
    mu::Real;
    delay::Integer=0,
    n_symbols::Integer=length(training),
    w0::Union{Nothing,AbstractVector{<:Number}}=nothing,
    tap_trajectory::Bool=false,
)
    n_taps = Int(n_taps)
    n_taps >= 1 || throw(ArgumentError("n_taps must be ≥ 1, got $n_taps"))
    delay = Int(delay)
    delay >= 0 || throw(ArgumentError("delay must be ≥ 0, got $delay"))
    μ = Float64(mu)
    (isfinite(μ) && μ > 0) || throw(ArgumentError(
        "mu must be finite and > 0, got $mu",
    ))
    length(training) >= 1 || throw(ArgumentError("training must be non-empty"))
    length(received) >= 1 || throw(ArgumentError("received must be non-empty"))
    n_symbols = Int(n_symbols)
    n_symbols >= 1 || throw(ArgumentError("n_symbols must be ≥ 1, got $n_symbols"))
    if delay > length(received) || n_symbols > length(received) - delay
        throw(ArgumentError(
            "received (length $(length(received))) is too short for n_symbols=$n_symbols and delay=$delay",
        ))
    end

    T = promote_type(eltype(received), eltype(training), Float64)
    if w0 !== nothing
        length(w0) == n_taps || throw(ArgumentError(
            "w0 must have length $n_taps, got $(length(w0))",
        ))
        T = promote_type(T, eltype(w0))
    end
    w = w0 === nothing ? zeros(T, n_taps) : Vector{T}(w0)
    u = Vector{T}(undef, n_taps)
    symbols = Vector{T}(undef, n_symbols)
    n_supervised = min(length(training), n_symbols)
    trajectory = tap_trajectory ? zeros(T, n_taps, n_supervised) : zeros(T, 0, 0)

    r = received
    @inbounds for m in 1:n_symbols
        k = m + delay
        for i in 1:n_taps
            idx = k - i + 1
            u[i] = 1 <= idx <= length(r) ? convert(T, r[idx]) : zero(T)
        end
        yhat = zero(T)
        for i in 1:n_taps
            yhat += w[i] * u[i]
        end
        symbols[m] = yhat
        if m <= n_supervised
            e = convert(T, training[m]) - yhat
            for i in 1:n_taps
                w[i] += μ * e * conj(u[i])
            end
            if tap_trajectory
                trajectory[:, m] = w
            end
        end
    end
    return tap_trajectory ? (symbols, trajectory) : symbols
end

@doc raw"""
    DFECoefficients

Finite-block decision-feedback equalizer for one convolution matrix.

# Fields
- `feedforward`: ``n \times (n + L - 1)`` matrix ``F``. Row ``m`` is the feedforward filter for symbol ``m``.
- `feedback`: ``n \times n`` matrix ``B``, strictly lower triangular. Entry ``B_{mj}`` (``j < m``) multiplies the decision of symbol ``j`` and is subtracted. The diagonal and the strict upper triangle are zero and are not read.
- `noise_var`: the value used to design the rows. `0` is the zero-forcing design. A positive value is the MMSE design. [`dfe_equalize`](@ref) does not read this field; the rows already contain the design.

Build with [`dfe_equalizer`](@ref). Apply with [`dfe_equalize`](@ref).
"""
struct DFECoefficients
    feedforward::Matrix
    feedback::Matrix
    noise_var::Float64
end

function _received_block_length(received::AbstractVector, h::AbstractVector)
    n = length(received) - length(h) + 1
    n >= 1 || throw(ArgumentError("received signal is too short for the given channel"))
    return n
end

function _nearest_symbol(z, alphabet)
    best = first(alphabet)
    best_d = abs2(z - best)
    @inbounds for i in 2:length(alphabet)
        d = abs2(z - alphabet[i])
        if d < best_d
            best = alphabet[i]
            best_d = d
        end
    end
    return best
end

@doc raw"""
    dfe_equalizer(h::AbstractVector{<:Number}, n::Integer, noise_var::Real) -> DFECoefficients

Finite-block decision-feedback equalizer for the same convolution matrix ``H`` as [`zf_equalizer`](@ref) and [`mmse_equalizer`](@ref). `n` is the number of symbols. `h` has length ``L``, so the observation has length ``n + L - 1``.

`noise_var = 0` is zero forcing. `noise_var > 0` is MMSE, with the same per-sample noise variance [`mmse_equalizer`](@ref) uses (real observation: ``N_0/2``; complex observation: ``N_0``). The value must be finite. Symbols are modelled as unit-variance and uncorrelated, as in the linear MMSE equalizer.

# Design
Let ``\Pi`` reverse the ``n`` columns of ``H``, so symbol 1 is the last column. Zero forcing takes the thin QR factorization ``H\Pi = QR``. MMSE takes

```math
\begin{bmatrix} H\Pi \\ \sqrt{\sigma^2}\, I \end{bmatrix} = QR
```

and keeps the top ``n + L - 1`` rows of ``Q`` as ``Q_1``. In either case ``R^{\mathrm{H}} R = \Pi^{\mathsf{T}}(H^{\mathrm{H}} H + \sigma^2 I)\Pi``, the normal matrix of [`mmse_equalizer`](@ref) with columns reversed. Back-substitution on ``z = Q_1^{\mathrm{H}} r``, with decisions standing in for symbols that have already been detected, is the block DFE. Dividing each row by its diagonal pivot folds the cursor gain into `feedforward`; `feedback` holds the strict upper part of ``R`` after that scaling, written in the original symbol order so it is strictly lower triangular.

Detection runs from symbol 1 to symbol ``n``. There is no V-BLAST reordering.

The feedforward row changes with the symbol index because the block is finite. This is not the time-invariant short feedforward/feedback pair of Proakis and Salehi, 5th ed., §9.5. The structure those pages describe — a feedforward filter, a feedback filter, coefficients from a zero-forcing or MSE criterion, and the assumption that past decisions are correct — is the one implemented here. A channel that is not identically zero gives ``H`` full column rank, as for [`zf_equalizer`](@ref). An identically zero `h` is rejected.

# Bias
The MMSE rows are the biased successive estimator, consistent with [`mmse_equalize`](@ref), which also does not remove the MMSE scale. At `noise_var = 0` that scale is 1 whenever the factorization exists.
"""
function dfe_equalizer(h::AbstractVector{<:Number}, n::Integer, noise_var::Real)
    n = Int(n)
    n >= 1 || throw(ArgumentError("n must be ≥ 1, got $n"))
    noise_var >= 0 && isfinite(noise_var) || throw(ArgumentError(
        "noise_var must be finite and ≥ 0, got $noise_var",
    ))
    H = _channel_matrix(h, n)
    Nr, ncol = size(H)
    T = eltype(H)
    # QR detects the last column first. Reversing columns detects symbol 1 first,
    # so feedback refers only to earlier symbols.
    Hperm = H[:, ncol:-1:1]
    σ = sqrt(Float64(noise_var))
    A = σ == 0 ? Hperm : vcat(Hperm, σ * Matrix{T}(I, ncol, ncol))
    F = qr(A)
    R = F.R
    Q1 = Matrix(F.Q[:, 1:ncol])[1:Nr, :]
    feedforward = zeros(T, ncol, Nr)
    feedback = zeros(T, ncol, ncol)
    @inbounds for k in 1:ncol
        ρ = R[k, k]
        iszero(ρ) && throw(ArgumentError(
            "channel matrix does not have full column rank",
        ))
        m = ncol - k + 1
        inv_ρ = inv(ρ)
        for i in 1:Nr
            feedforward[m, i] = conj(Q1[i, k]) * inv_ρ
        end
        for j in (k + 1):ncol
            m_past = ncol - j + 1
            feedback[m, m_past] = R[k, j] * inv_ρ
        end
    end
    return DFECoefficients(feedforward, feedback, Float64(noise_var))
end

@doc raw"""
    dfe_equalize(received, coef::DFECoefficients; alphabet, training=nothing)

Apply a [`DFECoefficients`](@ref) design. `received` has length `size(coef.feedforward, 2)`.

For symbol ``m``,

```math
\tilde{s}_m = F_{m,:}\, r - \sum_{j < m} B_{mj}\, \hat{s}_j.
```

``\hat{s}_j`` is `training[j]` when `j ≤ length(training)`, and otherwise the nearest point of `alphabet` to ``\tilde{s}_j``. Ties keep the earlier alphabet entry. The returned value is ``\tilde{s}``, the soft output before the slicer, not the decision that entered the feedback loop.

`training` is optional ideal feedback for a prefix. It is not checked against `alphabet`. A prefix that covers the whole block is genie-aided feedback: the slicer is never called, and error propagation is out of the problem. On the positions `training` does not cover, a wrong decision is multiplied by that row of `feedback` and added into every later soft output. Whether a later hard decision flips depends on the feedback tap and on which side of a decision boundary the soft value lands; a small tap can leave the decision unchanged. This function does not bound the length of an error burst. There is no Tomlinson–Harashima precoder. Proakis and Salehi, 5th ed., §9.5, state the same correct-decision assumption and the same propagation mechanism.

`alphabet` must be non-empty. The feedback filter does not read `feedback` on or above the diagonal.
"""
function dfe_equalize(
    received::AbstractVector{<:Number},
    coef::DFECoefficients;
    alphabet::AbstractVector{<:Number},
    training::Union{Nothing,AbstractVector{<:Number}}=nothing,
)
    n, Nr = size(coef.feedforward)
    (n >= 1 && Nr >= 1) || throw(ArgumentError("feedforward must be non-empty"))
    size(coef.feedback) == (n, n) || throw(ArgumentError(
        "feedback must be $(n)×$(n), got $(size(coef.feedback))",
    ))
    length(received) == Nr || throw(ArgumentError(
        "received length must be $Nr, got $(length(received))",
    ))
    length(alphabet) >= 1 || throw(ArgumentError("alphabet must be non-empty"))
    n_train = 0
    if training !== nothing
        n_train = length(training)
        1 <= n_train <= n || throw(ArgumentError(
            "training length must be between 1 and $n, got $n_train",
        ))
    end
    T = promote_type(eltype(received), eltype(coef.feedforward), eltype(alphabet), Float64)
    if training !== nothing
        T = promote_type(T, eltype(training))
    end
    ff = coef.feedforward * received
    soft = Vector{T}(undef, n)
    decisions = Vector{T}(undef, n)
    alpha = Vector{T}(alphabet)
    B = coef.feedback
    @inbounds for m in 1:n
        acc = convert(T, ff[m])
        for j in 1:(m - 1)
            acc -= convert(T, B[m, j]) * decisions[j]
        end
        soft[m] = acc
        decisions[m] = m <= n_train ? convert(T, training[m]) : _nearest_symbol(acc, alpha)
    end
    return soft
end

@doc raw"""
    zf_dfe(received, h; alphabet, training=nothing)

Zero-forcing block DFE. `received` is the length-``n + L - 1`` convolution of `n` symbols with `h`, the same layout as [`zf_equalize`](@ref) and as [`CommToolkit.Channel.apply_fir`](@ref).

Equivalent to `dfe_equalize(received, dfe_equalizer(h, n, 0); alphabet, training)`. See [`dfe_equalizer`](@ref) for the factorization and [`dfe_equalize`](@ref) for decisions and error propagation.

With every feedback decision equal to the transmitted symbol, a channel that is not identically zero is inverted up to linear-algebra roundoff. A single wrong decision can move later symbols; the noiseless statement does not survive that.
"""
function zf_dfe(
    received::AbstractVector{<:Number},
    h::AbstractVector{<:Number};
    alphabet::AbstractVector{<:Number},
    training::Union{Nothing,AbstractVector{<:Number}}=nothing,
)
    n = _received_block_length(received, h)
    return dfe_equalize(received, dfe_equalizer(h, n, 0); alphabet, training)
end

@doc raw"""
    mmse_dfe(received, h, noise_var; alphabet, training=nothing)

MMSE block DFE. `noise_var` is the per-sample variance used by [`mmse_equalize`](@ref), not a decibel SNR. `noise_var = 0` is [`zf_dfe`](@ref).

Equivalent to `dfe_equalize(received, dfe_equalizer(h, n, noise_var); alphabet, training)`. The soft output keeps the MMSE bias. Error propagation is the same as [`dfe_equalize`](@ref): ideal feedback holds only on the `training` prefix.
"""
function mmse_dfe(
    received::AbstractVector{<:Number},
    h::AbstractVector{<:Number},
    noise_var::Real;
    alphabet::AbstractVector{<:Number},
    training::Union{Nothing,AbstractVector{<:Number}}=nothing,
)
    n = _received_block_length(received, h)
    return dfe_equalize(received, dfe_equalizer(h, n, noise_var); alphabet, training)
end

end # module Equalization
