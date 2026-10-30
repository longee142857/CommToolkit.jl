module Coding

export conv_encode, viterbi_decode, bcjr_decode, sova_decode, distance_spectrum
export puncture_mask, puncture, depuncture
export turbo_interleaver, turbo_encode, turbo_decode
export GF2m, gf_add, gf_mul, gf_div, gf_inv, gf_pow
export hamming74_encode, hamming74_decode
export bch157_encode, bch157_decode
export RSCode, rs_code, rs_shorten, rs_encode, rs_decode
export LDPCCode, LDPCResult
export gallager_ldpc, array_ldpc, hamming74_ldpc, ldpc_from_H
export ldpc_encode, ldpc_decode, ldpc_decode_details, ldpc_syndrome
export ldpc_generator_matrix, ldpc_column_degrees, ldpc_row_degrees
export PolarCode, polar_code, polar_encode, polar_sc_decode

@doc raw"""
Rate-1/2 feedforward convolutional code with constraint length ``K = 3`` (memory ``m = 2``, four states).

Generators, octal `(7, 5)`:

```math
g_1(D) = 1 + D + D^2, \qquad g_2(D) = 1 + D^2.
```

This is the standard textbook code of free distance ``d_{\mathrm{free}} = 5`` (Proakis and Salehi, *Digital Communications*, 5th ed., §8.2; Lin and Costello, *Error Control Coding*, 2nd ed., Pearson, 2004, §11.1, where the same code is called a ``(2,1,2)`` code). Output bit order is ``(g_1, g_2)`` on each branch. Swapping the two generators does not change the distance spectrum.
"""
const CONV_G1 = 0b111
const CONV_G2 = 0b101
const CONV_K = 3
const CONV_N_STATES = 4

function _trellis()
    next_state = zeros(Int, CONV_N_STATES, 2)
    output = zeros(Int, CONV_N_STATES, 2)
    for s in 0:(CONV_N_STATES - 1)
        for b in 0:1
            reg = (s << 1) | b
            g1 = count_ones(reg & CONV_G1) & 1
            g2 = count_ones(reg & CONV_G2) & 1
            next_state[s + 1, b + 1] = reg & (CONV_N_STATES - 1)
            output[s + 1, b + 1] = (g1 << 1) | g2
        end
    end
    return next_state, output
end

const CONV_NEXT, CONV_OUT = _trellis()

@doc raw"""
    conv_encode(bits::AbstractVector{<:Integer}) -> Vector{Int}

Encode information bits with the rate-1/2, ``K = 3``, generators-`(7, 5)` code, then append ``K - 1 = 2`` zero tail bits so the encoder returns to state 0.

# Arguments
- `bits`: entries in ``\{0, 1\}``.

# Returns
A vector of length ``2(\texttt{length(bits)} + 2)``. Each information or tail bit produces the pair ``(g_1, g_2)``.

# References
Proakis and Salehi, 5th ed., §8.1–§8.2; Lin and Costello, 2nd ed., §11.1. The code is terminated, so the actual rate on a block of ``n`` information bits is ``n / (2n + 4)``.
"""
function conv_encode(bits::AbstractVector{<:Integer})
    all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("bits must be 0 or 1"))
    n = length(bits)
    coded = Vector{Int}(undef, 2 * (n + CONV_K - 1))
    state = 0
    idx = 1
    @inbounds for b in bits
        out = CONV_OUT[state + 1, b + 1]
        coded[idx] = (out >> 1) & 1
        coded[idx + 1] = out & 1
        state = CONV_NEXT[state + 1, b + 1]
        idx += 2
    end
    @inbounds for _ in 1:(CONV_K - 1)
        out = CONV_OUT[state + 1, 1] # input 0
        coded[idx] = (out >> 1) & 1
        coded[idx + 1] = out & 1
        state = CONV_NEXT[state + 1, 1]
        idx += 2
    end
    return coded
end

@doc raw"""
    viterbi_decode(coded::AbstractVector{<:Integer}) -> Vector{Int}
    viterbi_decode(llr::AbstractVector{<:AbstractFloat}) -> Vector{Int}

Viterbi decoding of a tail-terminated block from [`conv_encode`](@ref).

# Arguments
- `coded`: hard channel bits in ``\{0, 1\}``, even length, including the two tail bits. Branch metric is the Hamming distance on each pair.
- `llr`: finite channel LLRs of those coded bits, in the same order,

```math
L_j = \log\frac{p(y_j\mid c_j=0)}{p(y_j\mid c_j=1)}.
```

  `AbstractFloat` inputs are LLRs. They are not probabilities in ``[0,1]`` and not hard bits stored as floating-point 0/1. Integer vectors take the hard-decision method.

# Returns
The maximum-likelihood information bits. The ``K - 1`` tail bits are removed. Neither method is constrained to end in state 0; with a correct termination and few channel errors the best survivor ends at state 0 anyway. [`bcjr_decode`](@ref) and [`sova_decode`](@ref) do enforce termination.

# Branch metric
The soft decoder minimises, on the coded pair of a branch,

```math
\lambda(c_0,c_1) = -\bigl(L_0(1-2c_0)+L_1(1-2c_1)\bigr).
```

For BPSK in AWGN, ``L_j = 2y_j/\sigma^2`` (`demodulate_soft`), so ``\lambda`` is an affine function, with positive slope, of the squared Euclidean distance. A common positive scale does not change the survivor. The same correlation metric is the soft-decision Viterbi metric of Lin and Costello, 2nd ed., §12.1, and of Hagenauer and Hoeher, Proc. IEEE GLOBECOM, 1989, pp. 1680–1686. Substituting ``L_j = 1-2\hat{c}_j`` reproduces hard Hamming decoding, including ties. ``L_j = 0`` adds nothing to ``\lambda``. That is the erasure [`depuncture`](@ref) writes for a punctured mother bit; the integer method has no erasure symbol.

# Notes
``d_{\mathrm{free}} = 5``, so on a BSC every pattern of at most ``t = \lfloor(d_{\mathrm{free}}-1)/2\rfloor = 2`` coded-bit errors is correctable (Lin and Costello, 2nd ed., §12.1; Proakis and Salehi, 5th ed., §8.2).

The hard-decision bit-error union bound on a BSC with crossover ``p`` is

```math
P_b < \sum_{d=d_{\mathrm{free}}}^{\infty} c_d P_2(d),
```

where ``c_d`` is the total information weight of the simple error events of Hamming weight ``d`` ([`distance_spectrum`](@ref)) and ``P_2(d)`` is the pairwise error probability of a weight-``d`` path, including the factor ``1/2`` on the middle binomial term when ``d`` is even (Proakis and Salehi, 5th ed., §8.2; Lin and Costello, 2nd ed., §12.3).

With soft decisions and antipodal modulation the pairwise error probability at information-bit ``E_b/N_0`` is ``Q(\sqrt{2\,d\,R\,E_b/N_0})`` with ``R = 1/2`` (Proakis and Salehi, 5th ed., §8.2).
"""
function viterbi_decode(coded::AbstractVector{<:Integer})
    n_coded = length(coded)
    iseven(n_coded) || throw(ArgumentError("coded length must be even, got $n_coded"))
    n_sym = n_coded ÷ 2
    n_bits = n_sym - (CONV_K - 1)
    n_bits >= 0 || throw(ArgumentError("coded sequence is shorter than the termination tail"))

    n_states = CONV_N_STATES
    INF = typemax(Int) ÷ 2
    pm = fill(INF, n_states)
    new_pm = fill(INF, n_states)
    pm[1] = 0
    tb_state = zeros(Int, n_sym, n_states)
    tb_bit = zeros(Int, n_sym, n_states)

    @inbounds for t in 1:n_sym
        rx0 = coded[2t - 1]
        rx1 = coded[2t]
        fill!(new_pm, INF)
        for s in 0:(n_states - 1)
            pm[s + 1] == INF && continue
            for b in 0:1
                out = CONV_OUT[s + 1, b + 1]
                bm = (rx0 ⊻ ((out >> 1) & 1)) + (rx1 ⊻ (out & 1))
                ns = CONV_NEXT[s + 1, b + 1]
                cand = pm[s + 1] + bm
                if cand < new_pm[ns + 1]
                    new_pm[ns + 1] = cand
                    tb_state[t, ns + 1] = s
                    tb_bit[t, ns + 1] = b
                end
            end
        end
        pm, new_pm = new_pm, pm
    end

    _, best = findmin(pm)
    state = best - 1
    bits = Vector{Int}(undef, n_bits)
    @inbounds for t in n_sym:-1:1
        if t <= n_bits
            bits[t] = tb_bit[t, state + 1]
        end
        state = tb_state[t, state + 1]
    end
    return bits
end

function viterbi_decode(llr::AbstractVector{<:AbstractFloat})
    all(isfinite, llr) || throw(ArgumentError("LLRs must be finite"))
    n_coded = length(llr)
    iseven(n_coded) || throw(ArgumentError("LLR length must be even, got $n_coded"))
    n_sym = n_coded ÷ 2
    n_bits = n_sym - (CONV_K - 1)
    n_bits >= 0 || throw(ArgumentError("coded sequence is shorter than the termination tail"))

    n_states = CONV_N_STATES
    pm = fill(Inf, n_states)
    new_pm = similar(pm)
    pm[1] = 0.0
    tb_state = zeros(Int, n_sym, n_states)
    tb_bit = zeros(Int, n_sym, n_states)

    @inbounds for t in 1:n_sym
        L0 = Float64(llr[2t - 1])
        L1 = Float64(llr[2t])
        fill!(new_pm, Inf)
        for s in 0:(n_states - 1)
            isfinite(pm[s + 1]) || continue
            ps = pm[s + 1]
            for b in 0:1
                out = CONV_OUT[s + 1, b + 1]
                c0 = (out >> 1) & 1
                c1 = out & 1
                # Negative correlation with the BPSK image of the hypothesised bits.
                bm = -((1 - 2 * c0) * L0 + (1 - 2 * c1) * L1)
                ns = CONV_NEXT[s + 1, b + 1]
                cand = ps + bm
                if cand < new_pm[ns + 1]
                    new_pm[ns + 1] = cand
                    tb_state[t, ns + 1] = s
                    tb_bit[t, ns + 1] = b
                end
            end
        end
        shift = Inf
        for v in new_pm
            v < shift && (shift = v)
        end
        isfinite(shift) && (new_pm .-= shift)
        pm, new_pm = new_pm, pm
    end

    _, best = findmin(pm)
    state = best - 1
    bits = Vector{Int}(undef, n_bits)
    @inbounds for t in n_sym:-1:1
        if t <= n_bits
            bits[t] = tb_bit[t, state + 1]
        end
        state = tb_state[t, state + 1]
    end
    return bits
end

@doc raw"""
    bcjr_decode(llr::AbstractVector{<:AbstractFloat}; apriori=nothing) -> Vector{Float64}

Log-domain BCJR (MAP) decoding of a tail-terminated `(7, 5)` block. Returns extrinsic LLRs of the information bits.

# Arguments
- `llr`: finite channel LLRs of the coded bits from [`conv_encode`](@ref), same convention as the soft method of [`viterbi_decode`](@ref). Even length, at least the two tail pairs.
- `apriori`: optional finite a-priori LLRs of the information bits, ``L_a(u)=\log P(u=0)/P(u=1)``, length ``\texttt{length(llr)}/2 - 2``. `nothing` means ``L_a = 0``.

# Returns
Extrinsic information-bit LLRs, tail removed,

```math
L_e(u_k)=\log\frac{P(u_k=0\mid\mathbf{L})}{P(u_k=1\mid\mathbf{L})}-L_a(u_k).
```

A positive value favours bit 0. The code is not systematic: the information bit is not a channel observation, so no systematic channel LLR is subtracted. With ``L_a = 0``, ``L_e`` is the a-posteriori LLR. Hard decisions are ``L_e .< 0``.

# Recursion
Bahl, Cocke, Jelinek, and Raviv, IEEE Trans. Inf. Theory, 20(2):284–287, 1974, evaluated in the log domain. Branch transition ``(s',s)`` labelled by input ``u`` and coded bits ``(c_0,c_1)`` uses

```math
\log\gamma = \tfrac12 L_a(u)\,(1-2u) + \tfrac12\sum_{j=0}^{1} L(c_j)\,(1-2c_j).
```

Terms independent of the branch cancel in the normalised forward and backward recursions. ``\alpha_0`` is concentrated on state 0. The encoder's two tail bits are fixed at 0 and ``\beta`` at the end of the block is concentrated on state 0, which is the terminated-code boundary condition. ``\alpha`` and ``\beta`` are normalised by their maximum at each step. A channel LLR of 0, the value [`depuncture`](@ref) writes for a punctured bit, drops that bit from ``\gamma``.

The split into a priori and extrinsic LLRs is the L-value decomposition of Hagenauer, Offer, and Papke, IEEE Trans. Inf. Theory, 42(2):429–445, 1996. This function does not iterate with an outer code. [`sova_decode`](@ref) is the soft-output Viterbi approximation of the same terminated code. [`turbo_decode`](@ref) iterates a different, recursive systematic trellis; it does not call this function.
"""
function bcjr_decode(llr::AbstractVector{<:AbstractFloat}; apriori::Union{Nothing,AbstractVector{<:Real}}=nothing)
    all(isfinite, llr) || throw(ArgumentError("LLRs must be finite"))
    n_coded = length(llr)
    iseven(n_coded) || throw(ArgumentError("LLR length must be even, got $n_coded"))
    n_sym = n_coded ÷ 2
    n_bits = n_sym - (CONV_K - 1)
    n_bits >= 0 || throw(ArgumentError("coded sequence is shorter than the termination tail"))
    if apriori === nothing
        La = zeros(Float64, n_bits)
    else
        length(apriori) == n_bits || throw(ArgumentError(
            "apriori length must be $n_bits, got $(length(apriori))",
        ))
        all(isfinite, apriori) || throw(ArgumentError("apriori LLRs must be finite"))
        La = Float64.(apriori)
    end

    n_states = CONV_N_STATES
    logα = fill(-Inf, n_sym + 1, n_states)
    logβ = fill(-Inf, n_sym + 1, n_states)
    logα[1, 1] = 0.0
    logβ[n_sym + 1, 1] = 0.0

    @inbounds for t in 1:n_sym
        L0 = Float64(llr[2t - 1])
        L1 = Float64(llr[2t])
        tail = t > n_bits
        for s in 0:(n_states - 1)
            isfinite(logα[t, s + 1]) || continue
            for b in 0:1
                tail && b == 1 && continue
                ns = CONV_NEXT[s + 1, b + 1]
                γ = _branch_logγ(s, b, L0, L1, tail ? 0.0 : La[t])
                logα[t + 1, ns + 1] = _logaddexp(logα[t + 1, ns + 1], logα[t, s + 1] + γ)
            end
        end
        _normalise_log!(@view logα[t + 1, :])
    end

    @inbounds for t in n_sym:-1:1
        L0 = Float64(llr[2t - 1])
        L1 = Float64(llr[2t])
        tail = t > n_bits
        for s in 0:(n_states - 1)
            for b in 0:1
                tail && b == 1 && continue
                ns = CONV_NEXT[s + 1, b + 1]
                γ = _branch_logγ(s, b, L0, L1, tail ? 0.0 : La[t])
                logβ[t, s + 1] = _logaddexp(logβ[t, s + 1], logβ[t + 1, ns + 1] + γ)
            end
        end
        _normalise_log!(@view logβ[t, :])
    end

    Le = Vector{Float64}(undef, n_bits)
    acc0 = Vector{Float64}(undef, n_states * 2)
    acc1 = similar(acc0)
    @inbounds for t in 1:n_bits
        L0 = Float64(llr[2t - 1])
        L1 = Float64(llr[2t])
        n0 = 0
        n1 = 0
        for s in 0:(n_states - 1)
            for b in 0:1
                ns = CONV_NEXT[s + 1, b + 1]
                v = logα[t, s + 1] + _branch_logγ(s, b, L0, L1, La[t]) + logβ[t + 1, ns + 1]
                if b == 0
                    n0 += 1
                    acc0[n0] = v
                else
                    n1 += 1
                    acc1[n1] = v
                end
            end
        end
        Le[t] = _logsumexp(acc0, n0) - _logsumexp(acc1, n1) - La[t]
    end
    return Le
end

function _branch_logγ(s::Int, b::Int, L0::Float64, L1::Float64, La::Float64)
    out = CONV_OUT[s + 1, b + 1]
    c0 = (out >> 1) & 1
    c1 = out & 1
    return 0.5 * ((1 - 2 * c0) * L0 + (1 - 2 * c1) * L1 + (1 - 2 * b) * La)
end

function _logaddexp(a::Float64, b::Float64)
    if a == Inf || b == Inf
        return Inf
    elseif a == -Inf
        return b
    elseif b == -Inf
        return a
    end
    return a >= b ? a + log1p(exp(b - a)) : b + log1p(exp(a - b))
end

function _logsumexp(xs::Vector{Float64}, n::Int)
    m = -Inf
    @inbounds for i in 1:n
        xs[i] > m && (m = xs[i])
    end
    m == -Inf && return -Inf
    m == Inf && return Inf
    s = 0.0
    @inbounds for i in 1:n
        s += exp(xs[i] - m)
    end
    return m + log(s)
end

function _normalise_log!(v::AbstractVector{Float64})
    m = maximum(v)
    isfinite(m) && (v .-= m)
    return v
end

@doc raw"""
    sova_decode(llr::AbstractVector{<:AbstractFloat}; apriori=nothing) -> Vector{Float64}

Soft-output Viterbi (SOVA) decoding of a tail-terminated `(7, 5)` block. Returns extrinsic LLRs of the information bits.

# Arguments
- `llr`: finite channel LLRs of the coded bits from [`conv_encode`](@ref), same convention as [`bcjr_decode`](@ref). Even length, at least the two tail pairs. A `0` is an erasure, including a punctured bit restored by [`depuncture`](@ref).
- `apriori`: optional finite a-priori LLRs of the information bits, ``L_a(u)=\log P(u=0)/P(u=1)``, length ``\texttt{length(llr)}/2 - 2``. `nothing` means ``L_a = 0``.

# Returns
Extrinsic information-bit LLRs, tail removed,

```math
L_e(u_k)=(1-2\hat{u}_k)\,\rho_k-L_a(u_k).
```

``\hat{u}`` is the maximum-likelihood information sequence on the terminated trellis. ``\rho_k\ge 0`` is the Hagenauer–Hoeher reliability: the minimum metric gap ``\Delta`` between the survivor and the runner-up at merges along that path where the two paths disagree in bit ``k``. The branch metric is ``-\log\gamma`` with the same ``\log\gamma`` as [`bcjr_decode`](@ref), so a positive scale matches the BCJR LLR unit. A positive ``L_e`` favours bit 0. Hard decisions are ``L_e .< 0`` when ``L_a = 0``.

The two tail inputs are fixed at 0 and the survivor is the best path that ends in state 0, the same boundary condition as [`bcjr_decode`](@ref). [`viterbi_decode`](@ref) does not force that boundary. On a terminated word and a decisive channel the three hard decisions agree.

# Approximation
SOVA is not the BCJR. ``\rho_k`` is a competing-path gap, not a sum over paths, and it is often over-confident relative to the exact a-posteriori LLR (Hagenauer and Hoeher, Proc. IEEE GLOBECOM, 1989, pp. 1680–1686). Signs agree with [`bcjr_decode`](@ref) when the channel is noiseless or at high SNR and the ML and MAP bits coincide. Magnitudes are not required to match. If no competing path along the survivor ever disagrees in bit ``k``, ``\rho_k`` is replaced by ``\sum |L|+\sum |L_a|``, which is an upper bound on a path-metric gap in these units. The update window is the whole block, so time and memory grow as the square of the block length. This is a short-block decoder, not a sliding-window chip.

The code is not systematic, so no channel LLR is subtracted beyond the ``L_a`` term. This function does not iterate with an outer code.
"""
function sova_decode(llr::AbstractVector{<:AbstractFloat}; apriori::Union{Nothing,AbstractVector{<:Real}}=nothing)
    all(isfinite, llr) || throw(ArgumentError("LLRs must be finite"))
    n_coded = length(llr)
    iseven(n_coded) || throw(ArgumentError("LLR length must be even, got $n_coded"))
    n_sym = n_coded ÷ 2
    n_bits = n_sym - (CONV_K - 1)
    n_bits >= 0 || throw(ArgumentError("coded sequence is shorter than the termination tail"))
    if apriori === nothing
        La = zeros(Float64, n_bits)
    else
        length(apriori) == n_bits || throw(ArgumentError(
            "apriori length must be $n_bits, got $(length(apriori))",
        ))
        all(isfinite, apriori) || throw(ArgumentError("apriori LLRs must be finite"))
        La = Float64.(apriori)
    end

    n_states = CONV_N_STATES
    tb_state = zeros(Int, n_sym, n_states)
    tb_bit = zeros(Int, n_sym, n_states)
    tb_delta = fill(Inf, n_sym, n_states)
    tb_comp = fill(-1, n_sym, n_states)
    tb_comp_bit = zeros(Int, n_sym, n_states)
    pm = fill(Inf, n_states)
    pm[1] = 0.0
    best_m = Vector{Float64}(undef, n_states)
    second_m = Vector{Float64}(undef, n_states)
    best_s = Vector{Int}(undef, n_states)
    second_s = Vector{Int}(undef, n_states)
    best_b = Vector{Int}(undef, n_states)
    second_b = Vector{Int}(undef, n_states)

    @inbounds for t in 1:n_sym
        L0 = Float64(llr[2t - 1])
        L1 = Float64(llr[2t])
        tail = t > n_bits
        La_t = tail ? 0.0 : La[t]
        fill!(best_m, Inf)
        fill!(second_m, Inf)
        fill!(best_s, -1)
        fill!(second_s, -1)
        for s in 0:(n_states - 1)
            isfinite(pm[s + 1]) || continue
            ps = pm[s + 1]
            for b in 0:1
                tail && b == 1 && continue
                met = ps - _branch_logγ(s, b, L0, L1, La_t)
                ns = CONV_NEXT[s + 1, b + 1] + 1
                better = best_s[ns] < 0 || met < best_m[ns] ||
                    (met == best_m[ns] && (b < best_b[ns] || (b == best_b[ns] && s < best_s[ns])))
                if better
                    second_m[ns] = best_m[ns]
                    second_s[ns] = best_s[ns]
                    second_b[ns] = best_b[ns]
                    best_m[ns] = met
                    best_s[ns] = s
                    best_b[ns] = b
                elseif second_s[ns] < 0 || met < second_m[ns] ||
                       (met == second_m[ns] && (b < second_b[ns] || (b == second_b[ns] && s < second_s[ns])))
                    second_m[ns] = met
                    second_s[ns] = s
                    second_b[ns] = b
                end
            end
        end
        shift = Inf
        for v in best_m
            v < shift && (shift = v)
        end
        for ns in 1:n_states
            if best_s[ns] < 0
                pm[ns] = Inf
                tb_comp[t, ns] = -1
                tb_delta[t, ns] = Inf
                continue
            end
            pm[ns] = best_m[ns] - shift
            tb_state[t, ns] = best_s[ns]
            tb_bit[t, ns] = best_b[ns]
            if second_s[ns] >= 0
                tb_delta[t, ns] = second_m[ns] - best_m[ns]
                tb_comp[t, ns] = second_s[ns]
                tb_comp_bit[t, ns] = second_b[ns]
            else
                tb_delta[t, ns] = Inf
                tb_comp[t, ns] = -1
            end
        end
    end
    isfinite(pm[1]) || error("internal: SOVA found no path that ends in state 0")

    state_after = zeros(Int, n_sym + 1)
    @inbounds for t in n_sym:-1:1
        state_after[t] = tb_state[t, state_after[t + 1] + 1]
    end
    ml_bit = Vector{Int}(undef, n_sym)
    @inbounds for t in 1:n_sym
        ml_bit[t] = tb_bit[t, state_after[t + 1] + 1]
    end

    ρ = fill(Inf, n_bits)
    @inbounds for t in 1:n_sym
        s_end = state_after[t + 1] + 1
        cs = tb_comp[t, s_end]
        cs < 0 && continue
        Δ = tb_delta[t, s_end]
        if t <= n_bits && tb_comp_bit[t, s_end] != ml_bit[t]
            ρ[t] = min(ρ[t], Δ)
        end
        cj = cs
        for j in (t - 1):-1:1
            bj = tb_bit[j, cj + 1]
            if j <= n_bits && bj != ml_bit[j]
                ρ[j] = min(ρ[j], Δ)
            end
            cj = tb_state[j, cj + 1]
        end
    end

    bound = 0.0
    for x in llr
        bound += abs(Float64(x))
    end
    for x in La
        bound += abs(x)
    end
    Le = Vector{Float64}(undef, n_bits)
    @inbounds for t in 1:n_bits
        r = ρ[t]
        isfinite(r) || (r = bound)
        Le[t] = (1 - 2 * ml_bit[t]) * r - La[t]
    end
    return Le
end

function _binary_bits(bits::AbstractVector{<:Integer})
    out = Vector{Int}(undef, length(bits))
    i = 0
    for b in bits
        (b == 0 || b == 1) || throw(ArgumentError("bits must be 0 or 1"))
        i += 1
        @inbounds out[i] = Int(b)
    end
    return out
end

@doc raw"""
    distance_spectrum(dmax::Integer=16) -> (ad, cd)

Information-weight spectrum of the simple error events of the `(7, 5)`, ``K = 3`` code, enumerated on its trellis.

# Arguments
- `dmax`: largest output Hamming weight to retain. Must be at least 1.

# Returns
- `ad[d]`: number of simple paths of output weight ``d`` that leave state 0 and return to it.
- `cd[d]`: total information-bit weight of those paths (the coefficient in the bit-error union bound).

Both vectors are 1-based and of length `dmax`. Weights above `dmax` are dropped. For this code ``d_{\mathrm{free}} = 5`` and, through weight 16,

```math
a_d = 2^{d-5}, \qquad c_d = (d-4)\,2^{d-5} \quad (d \ge 5).
```

Those coefficients are the classical spectrum of the ``(2,1,2)`` code (Lin and Costello, 2nd ed., §11.1 and the transfer-function examples in §11.2; Proakis and Salehi, 5th ed., §8.2). The enumeration, not the closed form, is what this function returns. The weights are for the unpunctured mother code. Punctured spectra are not enumerated.
"""
function distance_spectrum(dmax::Integer=16)
    dmax = Int(dmax)
    dmax >= 1 || throw(ArgumentError("dmax must be ≥ 1"))
    A = zeros(Int, CONV_N_STATES, dmax + 1)
    B = zeros(Int, CONV_N_STATES, dmax + 1)
    # First branch that leaves the all-zero path: input 1 from state 0.
    ns0 = CONV_NEXT[1, 2]
    out0 = CONV_OUT[1, 2]
    ow0 = (out0 & 1) + ((out0 >> 1) & 1)
    if ow0 <= dmax
        A[ns0 + 1, ow0 + 1] = 1
        B[ns0 + 1, ow0 + 1] = 1
    end
    ad = zeros(Int, dmax)
    cd = zeros(Int, dmax)
    # No zero-output cycle exists off state 0; 2*dmax steps covers every simple event of weight ≤ dmax.
    for _ in 1:(2 * dmax)
        An = zeros(Int, CONV_N_STATES, dmax + 1)
        Bn = zeros(Int, CONV_N_STATES, dmax + 1)
        for s in 0:(CONV_N_STATES - 1)
            for w in 0:dmax
                paths = A[s + 1, w + 1]
                paths == 0 && continue
                in_weight = B[s + 1, w + 1]
                for b in 0:1
                    ns = CONV_NEXT[s + 1, b + 1]
                    out = CONV_OUT[s + 1, b + 1]
                    nw = w + (out & 1) + ((out >> 1) & 1)
                    nw > dmax && continue
                    if ns == 0
                        ad[nw] += paths
                        cd[nw] += in_weight + b * paths
                    else
                        An[ns + 1, nw + 1] += paths
                        Bn[ns + 1, nw + 1] += in_weight + b * paths
                    end
                end
            end
        end
        A, B = An, Bn
    end
    return ad, cd
end

# Standard periodic keep-masks on the (g1, g2) stream. A 1 keeps the bit.
# 2/3 and 3/4 are the Yasuda / IEEE 802.11 masks, not a dfree search for (7, 5).
const _PUNCTURE_12 = (1, 1)
const _PUNCTURE_23 = (1, 1, 1, 0)
const _PUNCTURE_34 = (1, 1, 1, 0, 0, 1)

@doc raw"""
    puncture_mask(rate::Rational) -> Vector{Int}

Periodic keep-mask taking the rate-1/2 `(g_1, g_2)` stream to `rate`. A 1 keeps the coded bit. [`conv_encode`](@ref) writes that order, tail included.

- `1//2` → `[1, 1]`, nothing deleted.
- `2//3` → `[1, 1, 1, 0]`, column-wise ``(g_1: 1\, 1),\ (g_2: 1\, 0)``.
- `3//4` → `[1, 1, 1, 0, 0, 1]`, column-wise ``(g_1: 1\, 1\, 0),\ (g_2: 1\, 0\, 1)``.

These are the puncturing patterns of Yasuda, Kashiki, and Hirata, IEEE Trans. Commun. 32(3):315–319, 1984, and the masks IEEE 802.11 applies to a rate-1/2 convolutional mother code. They are not a free-distance search for the `(7, 5)` generators. Any other rate, including a floating-point ratio, raises `ArgumentError`. Each call returns a fresh vector.

The nominal rate is ``(1/2) \cdot P/W`` with mask length ``P`` and Hamming weight ``W``. On a terminated block the tail and a partial last period remain, so the true rate is the number of information bits divided by the number of kept bits.
"""
function puncture_mask(rate::Rational)
    if rate == 1//2
        return Int[_PUNCTURE_12...]
    elseif rate == 2//3
        return Int[_PUNCTURE_23...]
    elseif rate == 3//4
        return Int[_PUNCTURE_34...]
    else
        throw(ArgumentError("puncture mask is defined for 1//2, 2//3, and 3//4, got $rate"))
    end
end

function puncture_mask(rate::Real)
    throw(ArgumentError("puncture rate must be a Rational such as 2//3, got $(repr(rate))"))
end

function _puncture_mask(mask::AbstractVector{<:Integer})
    isempty(mask) && throw(ArgumentError("puncture mask must be non-empty"))
    out = Vector{Int}(undef, length(mask))
    kept = 0
    i = 0
    for b in mask
        (b == 0 || b == 1) || throw(ArgumentError("puncture mask entries must be 0 or 1, got $b"))
        i += 1
        out[i] = Int(b)
        kept += Int(b)
    end
    kept > 0 || throw(ArgumentError("puncture mask must keep at least one position"))
    return out
end

function _nkept(n::Int, mask::Vector{Int})
    period = length(mask)
    nfull, rem = divrem(n, period)
    kept_period = 0
    @inbounds for b in mask
        kept_period += b
    end
    kept = nfull * kept_period
    @inbounds for i in 1:rem
        kept += mask[i]
    end
    return kept
end

@doc raw"""
    puncture(x::AbstractVector, mask::AbstractVector{<:Integer}) -> Vector
    puncture(x::AbstractVector, rate::Rational) -> Vector

Delete every entry of `x` at which the periodic keep-mask is 0. `rate` is passed through [`puncture_mask`](@ref).

`mask` entries are 0 or 1, and at least one entry is 1. The mask repeats from the first element of `x`. A trailing partial period uses the leading entries of `mask`. `x` is not modified. The result has the element type of `x`, so the same call punctures hard bits, symbols, or LLRs.
"""
function puncture(x::AbstractVector, mask::AbstractVector{<:Integer})
    m = _puncture_mask(mask)
    period = length(m)
    out = Vector{eltype(x)}(undef, _nkept(length(x), m))
    j = 0
    i = 0
    for v in x
        i += 1
        if m[mod1(i, period)] == 1
            j += 1
            @inbounds out[j] = v
        end
    end
    return out
end

function puncture(x::AbstractVector, rate::Rational)
    return puncture(x, puncture_mask(rate))
end

@doc raw"""
    depuncture(y, mask, n_mother::Integer; fill=zero(eltype(y))) -> Vector
    depuncture(y, rate::Rational, n_mother::Integer; fill=zero(eltype(y))) -> Vector

Write `y` back into a vector of length `n_mother`, placing `fill` where the periodic keep-mask is 0. `length(y)` must be the number of entries [`puncture`](@ref) would keep from a vector of length `n_mother`. `rate` is passed through [`puncture_mask`](@ref).

`n_mother` has to be supplied. Two mother lengths can share one punctured length when a partial period ends on deleted positions, so the punctured vector alone does not determine the mother length.

For channel LLRs in the sign convention of [`viterbi_decode`](@ref), `fill = 0` is an erasure: the soft branch metric and the BCJR branch weight ignore that bit. Pass the restored vector to the soft method of `viterbi_decode` or to [`bcjr_decode`](@ref). The integer method of `viterbi_decode` treats every entry as a hard bit and has no erasure. Hard decisions on the kept bits are connected by

```math
L = 1 - 2\hat{c}
```

on those bits and `fill = 0` on the deleted ones, then calling the soft decoder. Filling hard zeros and calling the integer decoder does not implement erasure decoding.
"""
function depuncture(
    y::AbstractVector,
    mask::AbstractVector{<:Integer},
    n_mother::Integer;
    fill=zero(eltype(y)),
)
    m = _puncture_mask(mask)
    n = Int(n_mother)
    n >= 0 || throw(ArgumentError("mother length must be ≥ 0, got $n"))
    expected = _nkept(n, m)
    ny = length(y)
    ny == expected || throw(ArgumentError(
        "punctured length $ny does not match mother length $n (mask keeps $expected)",
    ))
    T = promote_type(eltype(y), typeof(fill))
    out = Vector{T}(undef, n)
    fill!(out, convert(T, fill))
    src = eachindex(y)
    j = 0
    period = length(m)
    @inbounds for i in 1:n
        if m[mod1(i, period)] == 1
            j += 1
            out[i] = y[src[j]]
        end
    end
    return out
end

function depuncture(
    y::AbstractVector,
    rate::Rational,
    n_mother::Integer;
    fill=zero(eltype(y)),
)
    return depuncture(y, puncture_mask(rate), n_mother; fill=fill)
end

include("ldpc.jl")
include("turbo.jl")
include("algebraic.jl")
include("polar.jl")

end # module Coding
