# MIMO detectors: ML, ZF, MMSE, SIC
#
# y = H * x + n, where x is modulation symbol vector

using LinearAlgebra

# ML detector: exhaustive search over all possible symbol vectors
function ml_detect(y::Vector{ComplexF64}, H::Matrix{ComplexF64},
                   constellation::Vector{ComplexF64})
    n_tx = size(H, 2)
    best_x = zeros(ComplexF64, n_tx)
    best_dist = typemax(Float64)

    # Recursive search over all combinations
    function search(prefix::Vector{ComplexF64}, depth::Int)
        if depth > n_tx
            x = H * prefix
            dist = sum(abs2, y - x)
            if dist < best_dist
                best_dist = dist
                best_x .= prefix
            end
            return
        end
        for s in constellation
            push!(prefix, s)
            search(prefix, depth + 1)
            pop!(prefix)
        end
    end

    search(ComplexF64[], 1)
    best_x
end

# ZF detector: x_hat = pinv(H) * y
function zf_detect(y::Vector{ComplexF64}, H::Matrix{ComplexF64})
    pinv(H) * y
end

# MMSE detector: x_hat = (H'*H + σ²*I)^{-1} * H' * y
function mmse_detect(y::Vector{ComplexF64}, H::Matrix{ComplexF64}, noise_var::Real)
    n_tx = size(H, 2)
    (H' * H + noise_var * I) \ (H' * y)
end

# SIC (Successive Interference Cancellation) with ZF ordering
# Orders streams by SNR, detects strongest, cancels, repeats
function sic_detect(y::Vector{ComplexF64}, H::Matrix{ComplexF64},
                   constellation::Vector{ComplexF64})
    n_tx = size(H, 2)
    x_est = zeros(ComplexF64, n_tx)
    y_rem = copy(y)
    H_rem = copy(H)
    active = collect(1:n_tx)

    for _ in 1:n_tx
        # Compute ZF filter and pick best stream (max SNR = min row norm)
        W = pinv(H_rem)
        row_norms = [norm(W[i, :]) for i in 1:size(W, 1)]
        best = argmin(row_norms)

        # Detect best stream
        z = W[best, :]' * y_rem
        _, idx = findmin(abs2.(z .- constellation))
        x_est[active[best]] = constellation[idx]

        # Cancel interference
        y_rem .-= H_rem[:, best] * constellation[idx]
        H_rem = H_rem[:, setdiff(1:size(H_rem, 2), best)]
        deleteat!(active, best)
    end

    x_est
end
