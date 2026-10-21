module Coding

export conv_encode, viterbi_decode, bcjr_decode, distance_spectrum

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
The maximum-likelihood information bits. The ``K - 1`` tail bits are removed. Neither method is constrained to end in state 0; with a correct termination and few channel errors the best survivor ends at state 0 anyway. [`bcjr_decode`](@ref) does enforce termination.

# Branch metric
The soft decoder minimises, on the coded pair of a branch,

```math
\lambda(c_0,c_1) = -\bigl(L_0(1-2c_0)+L_1(1-2c_1)\bigr).
```

For BPSK in AWGN, ``L_j = 2y_j/\sigma^2`` (`demodulate_soft`), so ``\lambda`` is an affine function, with positive slope, of the squared Euclidean distance. A common positive scale does not change the survivor. The same correlation metric is the soft-decision Viterbi metric of Lin and Costello, 2nd ed., §12.1, and of Hagenauer and Hoeher, Proc. IEEE GLOBECOM, 1989, pp. 1680–1686. Substituting ``L_j = 1-2\hat{c}_j`` reproduces hard Hamming decoding, including ties.

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

Terms independent of the branch cancel in the normalised forward and backward recursions. ``\alpha_0`` is concentrated on state 0. The encoder's two tail bits are fixed at 0 and ``\beta`` at the end of the block is concentrated on state 0, which is the terminated-code boundary condition. ``\alpha`` and ``\beta`` are normalised by their maximum at each step.

The split into a priori and extrinsic LLRs is the L-value decomposition of Hagenauer, Offer, and Papke, IEEE Trans. Inf. Theory, 42(2):429–445, 1996. This function does not iterate with an outer code.
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

Those coefficients are the classical spectrum of the ``(2,1,2)`` code (Lin and Costello, 2nd ed., §11.1 and the transfer-function examples in §11.2; Proakis and Salehi, 5th ed., §8.2). The enumeration, not the closed form, is what this function returns.
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

end # module Coding
