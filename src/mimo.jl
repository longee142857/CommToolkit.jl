module MIMO

using LinearAlgebra

export mimo_channel, apply_mimo,
    ml_detect, zf_detect, mmse_detect, sic_detect,
    zf_precoder, svd_precoder,
    alamouti_encode, alamouti_decode_2x1, alamouti_decode_2x2

@doc raw"""
    mimo_channel(n_tx::Integer, n_rx::Integer) -> Matrix{ComplexF64}

i.i.d. Rayleigh matrix of size `(n_rx, n_tx)`. Each entry is ``\mathcal{CN}(0, 1)``, formed as `(randn + j·randn) / √2` (Tse and Viswanath, *Fundamentals of Wireless Communication*, Cambridge, 2005, §7.1).
"""
function mimo_channel(n_tx::Integer, n_rx::Integer)
    n_tx >= 1 && n_rx >= 1 || throw(ArgumentError("n_tx and n_rx must be ≥ 1"))
    return (randn(Int(n_rx), Int(n_tx)) .+ im .* randn(Int(n_rx), Int(n_tx))) ./ sqrt(2)
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

With the thin SVD ``H = U \operatorname{diag}(S) V^{\mathrm{H}}``, the noiseless link `H * V` equals `U * Diagonal(S)` (Tse and Viswanath, 2005, §7.1). No water-filling powers are applied; `S` is returned so a caller can compute them.

# Notes
Releases before 0.2 returned `Vᴴ` instead of `V`, because the third value of `svd`'s iterator was treated as `Vᵀ`. That matrix is not a precoder: it has the wrong shape whenever `n_tx ≠ rank`.
"""
function svd_precoder(H::AbstractMatrix{<:Number})
    F = svd(ComplexF64.(H))
    return F.V, F.S
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
