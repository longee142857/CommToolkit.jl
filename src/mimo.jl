module MIMO

using LinearAlgebra

export mimo_channel, exponential_correlation, kronecker_channel,
    apply_mimo, apply_kronecker,
    ml_detect, zf_detect, mmse_detect, sic_detect,
    zf_precoder, svd_precoder, waterfill, svd_waterfill,
    alamouti_encode, alamouti_decode_2x1, alamouti_decode_2x2

@doc raw"""
    mimo_channel(n_tx::Integer, n_rx::Integer) -> Matrix{ComplexF64}

i.i.d. Rayleigh matrix of size `(n_rx, n_tx)`. Each entry is ``\mathcal{CN}(0, 1)``, formed as `(randn + j·randn) / √2` (Tse and Viswanath, *Fundamentals of Wireless Communication*, Cambridge, 2005, §7.1).
"""
function mimo_channel(n_tx::Integer, n_rx::Integer)
    n_tx >= 1 && n_rx >= 1 || throw(ArgumentError("n_tx and n_rx must be ≥ 1"))
    return (randn(Int(n_rx), Int(n_tx)) .+ im .* randn(Int(n_rx), Int(n_tx))) ./ sqrt(2)
end

# Reject a Hermitian defect larger than this fraction of the largest entry.
# A negative eigenvalue inside the same absolute tolerance is rounding and is clipped.
const _CORR_HERM_TOL = 1e-8
# |ρ| may overshoot 1 by this much (for example abs(cis(θ))) and is then clipped to 1.
const _CORR_ABS_TOL = 1e-12

@doc raw"""
    exponential_correlation(n::Integer, ρ::Number; spacing::Real=1) -> Matrix{ComplexF64}

Loyka's exponential correlation matrix for `n` antennas on a line (Loyka, IEEE Commun. Lett., 5(9):369–371, 2001).

`ρ` is the correlation of two antennas **one unit of distance** apart, and ``|\rho| \le 1``. `spacing` is the distance between adjacent antennas, in that same unit. Elements `i` and `j` are ``|i-j|\,\texttt{spacing}`` apart, and

```math
R_{ij}
  = |\rho|^{|i-j|\,\texttt{spacing}}
    \exp\bigl(j\,\arg(\rho)\,(j-i)\,\texttt{spacing}\bigr),
```

with ``R_{ii} = 1``. For `spacing = 1` the first row is ``1, \rho, \rho^2, \ldots`` and the matrix is Hermitian Toeplitz. A real nonnegative `ρ`, and a negative real `ρ` whose lags ``|i-j|\,\texttt{spacing}`` are all integers, give a real symmetric matrix ``R_{ij} = \rho^{|i-j|\,\texttt{spacing}}`` (ordinary real power). A negative base with a non-integer lag uses ``\arg(\rho) = \pi``, so that entry is complex.

# Special cases
- `|ρ| = 0` returns `I` for every finite `spacing`, including `0`. No power of zero is formed.
- `spacing = 0` and `|ρ| > 0` returns the all-ones matrix: every pair sits at distance 0.
- `n = 1` returns a `1 × 1` matrix `[1]`, for every accepted `ρ` and `spacing`.
- `|ρ| = 1` is a rank-1 phase ramp. `|ρ|` may exceed 1 by at most `1e-12` and is then clipped to 1; a larger modulus is rejected.

`n ≥ 1`. `spacing` must be finite and `≥ 0`. `ρ` must be finite. Every accepted pair `(ρ, spacing)` produces a positive semidefinite correlation matrix with unit diagonal.

This is a Toeplitz coefficient for [`kronecker_channel`](@ref). It is not ``J_0(2\pi d/\lambda)``, not an array manifold, and not a 3GPP antenna pattern.
"""
function exponential_correlation(n::Integer, ρ::Number; spacing::Real=1)
    n = Int(n)
    n >= 1 || throw(ArgumentError("n must be ≥ 1, got $n"))
    ρc = ComplexF64(ρ)
    isfinite(ρc) || throw(ArgumentError("ρ must be finite, got $ρ"))
    s = Float64(spacing)
    (isfinite(s) && s >= 0) || throw(ArgumentError(
        "spacing must be finite and ≥ 0, got $spacing",
    ))
    mag = abs(ρc)
    mag <= 1 + _CORR_ABS_TOL || throw(ArgumentError(
        "correlation coefficient must satisfy |ρ| ≤ 1, got |ρ| = $mag",
    ))
    mag = min(mag, 1.0)
    R = Matrix{ComplexF64}(I, n, n)
    mag == 0 && return R
    # A real base with an integer lag uses the ordinary real power, so a
    # negative ρ stays real. A non-integer lag keeps arg(ρ) = π.
    ρ_real = imag(ρc) == 0
    ρr = copysign(mag, real(ρc))
    θ = angle(ρc)
    @inbounds for j in 1:n
        for i in 1:n
            i == j && continue
            lag = (j - i) * s
            dist = abs(lag)
            if ρ_real && dist == round(dist)
                R[i, j] = ρr^Int(round(dist))
            else
                amp = mag == 1 ? 1.0 : mag^dist
                R[i, j] = amp * cis(θ * lag)
            end
        end
    end
    return R
end

function _exact_identity(R::AbstractMatrix)
    n, m = size(R)
    n == m || return false
    @inbounds for j in 1:n
        for i in 1:n
            v = R[i, j]
            if i == j
                v == 1 || return false
            else
                v == 0 || return false
            end
        end
    end
    return true
end

"""Hermitian PSD square root. `what` names the argument in the error text."""
function _hermitian_sqrt(R::AbstractMatrix{<:Number}, what::AbstractString)
    n, m = size(R)
    n == m && n >= 1 || throw(ArgumentError(
        "$what must be square with at least one row, got $(size(R))",
    ))
    C = Matrix{ComplexF64}(undef, n, n)
    scale = 1.0
    @inbounds for j in 1:n
        for i in 1:n
            v = ComplexF64(R[i, j])
            isfinite(v) || throw(ArgumentError("$what must be finite"))
            C[i, j] = v
            scale = max(scale, abs(v))
        end
    end
    herm_err = 0.0
    @inbounds for j in 1:n
        for i in 1:n
            herm_err = max(herm_err, abs(C[i, j] - conj(C[j, i])))
        end
    end
    herm_err <= _CORR_HERM_TOL * scale || throw(ArgumentError(
        "$what must be Hermitian; maximum |R - Rᴴ| is $herm_err",
    ))
    F = eigen(Hermitian((C .+ C') ./ 2))
    tol = _CORR_HERM_TOL * scale
    λmin = F.values[1]
    λmin >= -tol || throw(ArgumentError(
        "$what must be positive semidefinite; smallest eigenvalue is $λmin",
    ))
    λ = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        v = F.values[i]
        λ[i] = v < 0 ? 0.0 : v
    end
    return ComplexF64.(F.vectors * Diagonal(sqrt.(λ)) * F.vectors')
end

@doc raw"""
    kronecker_channel(n_tx::Integer, n_rx::Integer) -> Matrix{ComplexF64}
    kronecker_channel(R_tx, R_rx) -> Matrix{ComplexF64}
    kronecker_channel(n_tx::Integer, n_rx::Integer, R_tx, R_rx) -> Matrix{ComplexF64}

One draw of a **flat Kronecker** spatially correlated Rayleigh matrix, size `n_rx × n_tx`,

```math
H = R_r^{1/2}\, H_w\, R_t^{1/2}.
```

`H_w` is the i.i.d. matrix of [`mimo_channel`](@ref): every real part, then every imaginary part, each ``\mathcal{CN}(0, 1)``. ``R_t^{1/2}`` and ``R_r^{1/2}`` are the Hermitian positive-semidefinite square roots, from the eigenvalue decomposition. On a real correlation matrix this is the same product as CommPy's `MIMOFlatChannel` (`sqrtm(R_rx) @ H_w @ sqrtm(R_tx)`). On a complex Hermitian factor the right multiplier is that Hermitian square root itself.

The separable correlation, with `p, r` receive indices and `q, s` transmit indices, is

```math
\mathbb{E}[H_{pq} H_{rs}^{*}] = (R_r)_{pr}\,(R_t)_{sq}.
```

Column-major `vec` (Julia's `vec`) therefore has covariance ``R_t^{\mathsf{T}} \otimes R_r``, because ``\mathrm{vec}(ABC) = (C^{\mathsf{T}} \otimes A)\,\mathrm{vec}(B)`` with ``A = R_r^{1/2}`` and ``C = R_t^{1/2}``. The one-sided moments are ``\mathbb{E}[H H^{\mathrm{H}}] = \mathrm{tr}(R_t)\, R_r`` and ``\mathbb{E}[H^{\mathrm{H}} H] = \mathrm{tr}(R_r)\, R_t``.

Kermoal, Schumacher, Pedersen, Mogensen, and Frederiksen (IEEE J. Sel. Areas Commun., 20(6):1211–1226, 2002) factor the narrowband covariance as a Kronecker product of the two link ends; the Appendix shows that the correlation of two channel coefficients is the product of the antenna correlations. The square-root product above is that factorization. `H_w` is the i.i.d. matrix of Tse and Viswanath, 2005, §7.1.

# Arguments
- `R_tx`: transmit correlation, `n_tx × n_tx`, Hermitian and positive semidefinite. Passed as given: the diagonal is not scaled to 1.
- `R_rx`: receive correlation, `n_rx × n_rx`, same requirements.

A unit diagonal is the usual correlation normalisation. Then ``\mathbb{E}[|H_{ij}|^2] = 1``, matching [`mimo_channel`](@ref), and ``\mathrm{tr}(R_t) = n_t``. A matrix whose largest entry of ``R - R^{\mathrm{H}}`` is at most `1e-8` times ``\max(1, \|R\|_\infty)`` is replaced by its Hermitian part before the square root. An eigenvalue below ``-10^{-8}`` times that scale is rejected. A smaller negative eigenvalue is clipped to 0, so a singular factor such as the all-ones matrix (`ρ = 1`) is accepted.

# Calls
- `kronecker_channel(n_tx, n_rx)` is [`mimo_channel`](@ref). Both factors are identities.
- Exact identity matrices take that same path: the draw is `H_w`, not a product by `I`, and the global random stream matches `mimo_channel`. Any other factor consumes the same `H_w` and then multiplies. Validation of `R_tx` and `R_rx` happens before the draw, so a rejected matrix does not advance the stream.
- The four-argument method checks that `R_tx` is `n_tx × n_tx` and `R_rx` is `n_rx × n_rx`, then calls the two-matrix method.

# What this is not
One flat matrix per call. Distinct calls are independent. There is no Doppler spectrum and no delay tap. This is not the Weichselberger model (no coupling matrix between the two ends) and not a 3GPP TR 38.901 CDL or UMi geometry. [`exponential_correlation`](@ref) supplies a Loyka factor; it is not ``J_0``.

Pass `H` to [`apply_mimo`](@ref), [`ml_detect`](@ref), [`zf_detect`](@ref), [`mmse_detect`](@ref), and [`sic_detect`](@ref). Those functions condition on the realized `H`. [`mmse_detect`](@ref) still assumes uncorrelated unit-energy symbols and spatially white noise; it does not whiten `R_t` or `R_r`. [`apply_kronecker`](@ref) draws `H` and then calls `apply_mimo`.
"""
function kronecker_channel(R_tx::AbstractMatrix{<:Number}, R_rx::AbstractMatrix{<:Number})
    n_tx = size(R_tx, 1)
    n_rx = size(R_rx, 1)
    n_tx >= 1 && n_rx >= 1 || throw(ArgumentError(
        "R_tx and R_rx must each have at least one row, got $(size(R_tx)) and $(size(R_rx))",
    ))
    # Square roots before the draw: a rejected factor must not consume randn.
    tx_sqrt = _exact_identity(R_tx) ? nothing : _hermitian_sqrt(R_tx, "R_tx")
    rx_sqrt = _exact_identity(R_rx) ? nothing : _hermitian_sqrt(R_rx, "R_rx")
    Hw = (randn(n_rx, n_tx) .+ im .* randn(n_rx, n_tx)) ./ sqrt(2)
    if rx_sqrt !== nothing
        Hw = rx_sqrt * Hw
    end
    if tx_sqrt !== nothing
        Hw = Hw * tx_sqrt
    end
    return Hw
end

function kronecker_channel(n_tx::Integer, n_rx::Integer, R_tx::AbstractMatrix{<:Number},
                           R_rx::AbstractMatrix{<:Number})
    nt = Int(n_tx)
    nr = Int(n_rx)
    nt >= 1 && nr >= 1 || throw(ArgumentError("n_tx and n_rx must be ≥ 1"))
    size(R_tx) == (nt, nt) || throw(ArgumentError(
        "R_tx must be $(nt)×$(nt), got $(size(R_tx))",
    ))
    size(R_rx) == (nr, nr) || throw(ArgumentError(
        "R_rx must be $(nr)×$(nr), got $(size(R_rx))",
    ))
    return kronecker_channel(R_tx, R_rx)
end

function kronecker_channel(n_tx::Integer, n_rx::Integer)
    return mimo_channel(n_tx, n_rx)
end

@doc raw"""
    apply_mimo(x, H, noise_var::Real) -> Vector{ComplexF64}

``y = Hx + n`` with ``n \sim \mathcal{CN}(0, \texttt{noise\_var}\, I)``.

# Arguments
- `x`: length-`n_tx` symbol vector. Real entries are promoted to complex.
- `H`: `n_rx × n_tx` matrix.
- `noise_var`: total complex noise variance ``\mathbb{E}[|n_i|^2]``. Each real dimension has variance `noise_var / 2`. `0` disables the noise. For unit-energy BPSK and an ``E_b/N_0`` of ``\rho``, the matching variance on an identity channel is ``N_0 = 1/\rho``, not ``N_0/2``.

# References
Tse and Viswanath, 2005, §7.1.
"""
function apply_mimo(x::AbstractVector{<:Number}, H::AbstractMatrix{<:Number}, noise_var::Real)
    noise_var >= 0 || throw(ArgumentError("noise_var must be ≥ 0, got $noise_var"))
    size(H, 2) == length(x) || throw(ArgumentError(
        "H has $(size(H, 2)) columns but x has length $(length(x))",
    ))
    y = ComplexF64.(H) * ComplexF64.(x)
    if noise_var > 0
        σ = sqrt(float(noise_var) / 2)
        y = y .+ σ .* randn(length(y)) .+ im .* σ .* randn(length(y))
    end
    return y
end

@doc raw"""
    apply_kronecker(x, R_tx, R_rx, noise_var::Real) -> (y, H)

Draw one flat Kronecker channel and apply it.

`H` is [`kronecker_channel`](@ref)`(R_tx, R_rx)`. `y` is [`apply_mimo`](@ref)`(x, H, noise_var)`. `length(x)` must equal the number of rows of `R_tx`, which is `n_tx`. `R_rx` sets `n_rx`.

The draw order on the global stream is the `H_w` of `kronecker_channel` (all real parts, then all imaginary parts), then the noise of `apply_mimo` (all real parts, then all imaginary parts) when `noise_var > 0`. `noise_var = 0` draws `H` and does not draw noise. A negative `noise_var`, or an `x` whose length disagrees with `R_tx`, is rejected before either draw.

This is the same `H` the existing detectors take. It does not change [`mimo_channel`](@ref) or [`apply_mimo`](@ref).
"""
function apply_kronecker(
    x::AbstractVector{<:Number},
    R_tx::AbstractMatrix{<:Number},
    R_rx::AbstractMatrix{<:Number},
    noise_var::Real,
)
    noise_var >= 0 || throw(ArgumentError("noise_var must be ≥ 0, got $noise_var"))
    length(x) == size(R_tx, 1) || throw(ArgumentError(
        "R_tx is $(size(R_tx, 1))×$(size(R_tx, 2)) but x has length $(length(x))",
    ))
    H = kronecker_channel(R_tx, R_rx)
    return apply_mimo(x, H, noise_var), H
end

@doc raw"""
    ml_detect(y, H, constellation) -> Vector{ComplexF64}

Exhaustive maximum-likelihood detection

```math
\hat x = \arg\min_{x \in \mathcal{X}^{n_t}} \|y - Hx\|^2
```

(Tse and Viswanath, 2005, §7.1). The search is ``|\mathcal{X}|^{n_t}``. Ties keep the first vector in lexicographic order of the constellation index. Intended for small ``n_t`` (two antennas, BPSK or QPSK).
"""
function ml_detect(y::AbstractVector{<:Number}, H::AbstractMatrix{<:Number},
                   constellation::AbstractVector{<:Number})
    n_tx = size(H, 2)
    n_tx >= 1 || throw(ArgumentError("H must have at least one column"))
    isempty(constellation) && throw(ArgumentError("constellation must be non-empty"))
    y = ComplexF64.(y)
    Hc = ComplexF64.(H)
    alphabet = ComplexF64.(constellation)
    length(y) == size(Hc, 1) || throw(ArgumentError("length(y) must equal the number of rows of H"))
    best_x = zeros(ComplexF64, n_tx)
    best_dist = Inf
    prefix = ComplexF64[]
    function search(depth::Int)
        if depth > n_tx
            dist = sum(abs2, y .- Hc * prefix)
            if dist < best_dist
                best_dist = dist
                best_x .= prefix
            end
            return
        end
        for s in alphabet
            push!(prefix, s)
            search(depth + 1)
            pop!(prefix)
        end
        return nothing
    end
    search(1)
    return best_x
end

@doc raw"""
    zf_detect(y, H) -> Vector{ComplexF64}

Zero-forcing detection ``\hat x = H^{+} y`` via `LinearAlgebra.pinv` (Tse and Viswanath, 2005, §7.2). `H` must have full column rank for the estimate to be finite and useful; a fat or square full-rank matrix is the usual case (``n_r \ge n_t``).
"""
function zf_detect(y::AbstractVector{<:Number}, H::AbstractMatrix{<:Number})
    length(y) == size(H, 1) || throw(ArgumentError("length(y) must equal the number of rows of H"))
    return pinv(ComplexF64.(H)) * ComplexF64.(y)
end

@doc raw"""
    mmse_detect(y, H, noise_var::Real) -> Vector{ComplexF64}

Linear MMSE detection for uncorrelated unit-energy symbols,

```math
\hat x = \left(H^{\mathrm{H}} H + \sigma^2 I\right)^{-1} H^{\mathrm{H}} y,
```

where ``\sigma^2`` is `noise_var`, the same complex noise variance passed to [`apply_mimo`](@ref) (Tse and Viswanath, 2005, §7.2). `noise_var = 0` is the left pseudo-inverse and agrees with ZF when ``H`` has full column rank.
"""
function mmse_detect(y::AbstractVector{<:Number}, H::AbstractMatrix{<:Number}, noise_var::Real)
    noise_var >= 0 || throw(ArgumentError("noise_var must be ≥ 0, got $noise_var"))
    length(y) == size(H, 1) || throw(ArgumentError("length(y) must equal the number of rows of H"))
    Hc = ComplexF64.(H)
    n_tx = size(Hc, 2)
    return (Hc' * Hc + float(noise_var) * I(n_tx)) \ (Hc' * ComplexF64.(y))
end

@doc raw"""
    sic_detect(y, H, constellation) -> Vector{ComplexF64}

Ordered successive interference cancellation with a ZF filter at each stage (V-BLAST nulling and cancellation; Wolniansky, Foschini, Golden, and Valenzuela, Proc. URSI ISSSE, 1998; the ordering rule is also in Tse and Viswanath, 2005, §7.3).

At each stage the remaining stream with the smallest ZF-row norm is sliced to the nearest point of `constellation`, its contribution is subtracted, and its column is removed. The soft statistic is the corresponding entry of ``H^{+} y`` (no conjugate on the filter row). There is no MMSE nulling and no optimal ordering beyond that row-norm rule.
"""
function sic_detect(y::AbstractVector{<:Number}, H::AbstractMatrix{<:Number},
                    constellation::AbstractVector{<:Number})
    n_tx = size(H, 2)
    n_tx >= 1 || throw(ArgumentError("H must have at least one column"))
    isempty(constellation) && throw(ArgumentError("constellation must be non-empty"))
    alphabet = ComplexF64.(constellation)
    x_est = zeros(ComplexF64, n_tx)
    y_rem = ComplexF64.(y)
    H_rem = ComplexF64.(H)
    active = collect(1:n_tx)
    for _ in 1:n_tx
        W = pinv(H_rem)
        row_norms = [norm(@view W[i, :]) for i in 1:size(W, 1)]
        best = argmin(row_norms)
        # (W*y)[i] = sum_j W[i,j] y[j]. The row must not be conjugated:
        # `row' * y` is the adjoint product and is wrong for a complex ZF filter.
        z = (W * y_rem)[best]
        _, idx = findmin(abs2.(z .- alphabet))
        x_est[active[best]] = alphabet[idx]
        y_rem = y_rem .- H_rem[:, best] .* alphabet[idx]
        keep = setdiff(1:size(H_rem, 2), best)
        H_rem = H_rem[:, keep]
        deleteat!(active, best)
    end
    return x_est
end

@doc raw"""
    zf_precoder(H) -> Matrix{ComplexF64}

Zero-forcing precoder for ``n_t \ge n_r``,

```math
P = H^{\mathrm{H}} (H H^{\mathrm{H}})^{-1},
```

so ``H P = I_{n_r}`` when ``H`` has full row rank (Tse and Viswanath, 2005, §7.4). The columns are not normalised to a power constraint: ``\|P\|_F`` grows as ``H`` becomes ill-conditioned. `H` is `n_rx × n_tx`.
"""
function zf_precoder(H::AbstractMatrix{<:Number})
    size(H, 2) >= size(H, 1) || throw(ArgumentError(
        "ZF precoder requires n_tx ≥ n_rx, got $(size(H, 2)) × $(size(H, 1)) (columns × rows)",
    ))
    Hc = ComplexF64.(H)
    return Hc' / (Hc * Hc')
end

@doc raw"""
    svd_precoder(H) -> (V, S)

Right-singular-vector precoder of `H`.

# Returns
- `V`: `n_tx × r` matrix whose columns are the right singular vectors (`LinearAlgebra.svd(H).V`).
- `S`: the singular values, largest first.

With the thin SVD ``H = U \operatorname{diag}(S) V^{\mathrm{H}}``, the noiseless link `H * V` equals `U * Diagonal(S)` (Tse and Viswanath, 2005, §7.1). No water-filling powers are applied. [`waterfill`](@ref) allocates a total power across `S`, and [`svd_waterfill`](@ref) folds that allocation into the precoder.

# Notes
Releases before 0.2 returned `Vᴴ` instead of `V`, because the third value of `svd`'s iterator was treated as `Vᵀ`. That matrix is not a precoder: it has the wrong shape whenever `n_tx ≠ rank`.
"""
function svd_precoder(H::AbstractMatrix{<:Number})
    F = svd(ComplexF64.(H))
    return F.V, F.S
end

@doc raw"""
    waterfill(singular_values, power::Real, noise_var::Real) -> (p, level)

Total-power water-filling on the singular values of one MIMO channel (Tse and Viswanath, *Fundamentals of Wireless Communication*, Cambridge, 2005, §7.1.1, eq. (7.10)–(7.11); the parallel-channel solution is §5.3.3).

# Arguments
- `singular_values`: nonnegative singular values ``\lambda_i``, in any order. These are the ``S`` returned by [`svd_precoder`](@ref), not the eigenvalues of ``H^{\mathrm{H}} H``. A zero singular value is a mode that cannot carry power.
- `power`: total transmit power ``P \ge 0``, finite. This is ``\sum_i p_i``, the same constraint as ``\mathbb{E}[\|x\|^2] = P`` in §7.1.1.
- `noise_var`: complex noise variance ``N_0 = \mathbb{E}[|n_i|^2]`` on each receive antenna, finite and ``\ge 0``. This is the `noise_var` of [`apply_mimo`](@ref), not ``N_0/2``.

# Returns
- `p`: powers, one per input singular value, in the same order. ``p_i = \max(0,\, \mu - N_0/\lambda_i^2)`` and ``p_i = 0`` when ``\lambda_i = 0``.
- `level`: the water level ``\mu`` in eq. (7.11). It is not the Lagrange multiplier of §5.3.3. When `power` is 0 and some ``\lambda_i > 0``, ``\mu = N_0 / \lambda_{\max}^2``, which leaves every mode dry. When every singular value is 0 and `power` is 0, `level` is `0.0`.

`power > 0` with every singular value equal to 0 throws: the constraint cannot be met on the range of ``H``. ``N_0 = 0`` shares `power` uniformly across the strictly positive singular values.

At high SNR (``N_0`` negligible beside ``P \lambda_{\min}^2``) every positive mode receives nearly `power` divided by the number of positive modes. At low SNR the allocation collapses onto the strongest mode (Tse and Viswanath, 2005, §7.1.1). The search is the water-level sweep: modes ordered from strongest to weakest, and the largest active set whose level still clears that set's noisiest mode.
"""
function waterfill(singular_values::AbstractVector{<:Real}, power::Real, noise_var::Real)
    isempty(singular_values) && throw(ArgumentError("singular_values must be non-empty"))
    P = float(power)
    N0 = float(noise_var)
    isfinite(P) && P >= 0 || throw(ArgumentError(
        "power must be finite and ≥ 0, got $power",
    ))
    isfinite(N0) && N0 >= 0 || throw(ArgumentError(
        "noise_var must be finite and ≥ 0, got $noise_var",
    ))
    n = length(singular_values)
    σ = Vector{Float64}(undef, n)
    i = 0
    for s0 in singular_values
        i += 1
        s = float(s0)
        isfinite(s) && s >= 0 || throw(ArgumentError(
            "singular values must be finite and ≥ 0, got $s0",
        ))
        σ[i] = s
    end
    p = zeros(Float64, n)
    positive = Int[]
    λmax = 0.0
    for j in 1:n
        λ = abs2(σ[j])
        if λ > 0
            push!(positive, j)
            λmax = max(λmax, λ)
        end
    end
    if isempty(positive)
        P == 0 && return p, 0.0
        throw(ArgumentError("power is positive but every singular value is zero"))
    end
    if P == 0
        return p, N0 / λmax
    end
    # Strongest eigenmode first, so the noise-over-gain ν is nondecreasing.
    sort!(positive; by=j -> σ[j], rev=true)
    m = length(positive)
    ν = Vector{Float64}(undef, m)
    @inbounds for j in 1:m
        ν[j] = N0 == 0 ? 0.0 : N0 / abs2(σ[positive[j]])
    end
    # k = 1 always clears: μ = P + ν_1 > ν_1. Larger sets are tried first.
    k_star = 1
    level = P + ν[1]
    for k in m:-1:2
        μ = (P + sum(@view ν[1:k])) / k
        if μ > ν[k]
            k_star = k
            level = μ
            break
        end
    end
    @inbounds for j in 1:k_star
        p[positive[j]] = level - ν[j]
    end
    return p, level
end

@doc raw"""
    svd_waterfill(H, power::Real, noise_var::Real) -> (F, p, level)

SVD precoder with the total-power water-fill of [`waterfill`](@ref).

# Arguments
- `H`: `n_rx × n_tx` channel, at least ``1 \times 1``.
- `power`, `noise_var`: forwarded to [`waterfill`](@ref). `noise_var` is the complex per-antenna variance used by [`apply_mimo`](@ref).

# Returns
- `F`: `n_tx × r` precoder, ``r = \min(n_{rx}, n_{tx})``,

```math
F = V \operatorname{diag}\bigl(\sqrt{p_1}, \ldots, \sqrt{p_r}\bigr),
```

where `V` and the singular values are those of [`svd_precoder`](@ref). Columns with ``p_i = 0`` are kept, as zeros, so column `i` is still singular mode `i` (largest singular value first). For unit-energy uncorrelated symbols ``s``, ``x = F s`` uses transmit power ``\|F\|_F^2 = \sum_i p_i``.
- `p`: the power allocation, length `r`.
- `level`: the water level ``\mu``.

[`svd_precoder`](@ref) is unchanged and still applies no powers. This is single-link water-filling, not iterative multiuser water-filling.
"""
function svd_waterfill(H::AbstractMatrix{<:Number}, power::Real, noise_var::Real)
    size(H, 1) >= 1 && size(H, 2) >= 1 || throw(ArgumentError(
        "H must have at least one row and one column, got $(size(H))",
    ))
    V, S = svd_precoder(H)
    p, level = waterfill(S, power, noise_var)
    F = V * Diagonal(sqrt.(p))
    return F, p, level
end

@doc raw"""
    alamouti_encode(x1::Number, x2::Number) -> Matrix{ComplexF64}

Alamouti space-time block code (Alamouti, IEEE J. Sel. Areas Commun., 16(8):1451–1458, 1998). The returned `2 × 2` matrix has transmitters in the rows and time slots in the columns,

```math
\begin{bmatrix}
x_1 & -x_2^{*} \\
x_2 & x_1^{*}
\end{bmatrix}.
```
"""
function alamouti_encode(x1::Number, x2::Number)
    return ComplexF64[x1  -conj(x2);
                      x2   conj(x1)]
end

@doc raw"""
    alamouti_decode_2x1(r1, r2, h1, h2) -> (x1_hat, x2_hat)

Decode one Alamouti block with a single receive antenna and perfect CSI. `r1`, `r2` are the received samples in the two time slots; `h1`, `h2` are the coefficients from the two transmitters.

```math
\hat x_1 = \frac{h_1^{*} r_1 + h_2 r_2^{*}}{|h_1|^2 + |h_2|^2}, \qquad
\hat x_2 = \frac{h_2^{*} r_1 - h_1 r_2^{*}}{|h_1|^2 + |h_2|^2}.
```

The denominator is the array gain of the code (Alamouti, 1998, eq. (11)–(12) of the complex-baseband form). A zero channel returns `(0, 0)`.
"""
function alamouti_decode_2x1(r1::Number, r2::Number, h1::Number, h2::Number)
    h1 = ComplexF64(h1)
    h2 = ComplexF64(h2)
    r1 = ComplexF64(r1)
    r2 = ComplexF64(r2)
    denom = abs2(h1) + abs2(h2)
    if denom == 0
        return ComplexF64(0), ComplexF64(0)
    end
    x1 = (conj(h1) * r1 + h2 * conj(r2)) / denom
    x2 = (conj(h2) * r1 - h1 * conj(r2)) / denom
    return x1, x2
end

@doc raw"""
    alamouti_decode_2x2(r1::AbstractVector, r2::AbstractVector, H::AbstractMatrix) -> (x1_hat, x2_hat)

Decode one Alamouti block with two receive antennas. `H` is `2 × 2` with `H[rx, tx]`. `r1` and `r2` are the length-2 received vectors in the two time slots.

The combiner adds the two single-antenna Alamouti statistics and normalises by ``\sum_{i,j} |H_{ij}|^2`` (Alamouti, 1998, §III). A zero channel returns `(0, 0)`.
"""
function alamouti_decode_2x2(r1::AbstractVector{<:Number}, r2::AbstractVector{<:Number},
                             H::AbstractMatrix{<:Number})
    size(H) == (2, 2) || throw(ArgumentError("H must be 2×2, got $(size(H))"))
    length(r1) == 2 && length(r2) == 2 || throw(ArgumentError("r1 and r2 must have length 2"))
    h11, h21 = ComplexF64(H[1, 1]), ComplexF64(H[2, 1])
    h12, h22 = ComplexF64(H[1, 2]), ComplexF64(H[2, 2])
    r11, r21 = ComplexF64(r1[1]), ComplexF64(r1[2])
    r12, r22 = ComplexF64(r2[1]), ComplexF64(r2[2])
    denom = abs2(h11) + abs2(h21) + abs2(h12) + abs2(h22)
    if denom == 0
        return ComplexF64(0), ComplexF64(0)
    end
    x1 = (conj(h11) * r11 + h12 * conj(r12) + conj(h21) * r21 + h22 * conj(r22)) / denom
    x2 = (conj(h12) * r11 - h11 * conj(r12) + conj(h22) * r21 - h21 * conj(r22)) / denom
    return x1, x2
end

end # module MIMO
