module Coding

export conv_encode, viterbi_decode, distance_spectrum

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

Hard-decision Viterbi decoding of a tail-terminated block from [`conv_encode`](@ref).

# Arguments
- `coded`: bits in ``\{0, 1\}``, even length, produced by the `(7, 5)` encoder including its two tail bits. Branch metric is the Hamming distance on each pair of coded bits.

# Returns
The maximum-likelihood information bits. The ``K - 1`` tail bits are removed. The decoder is not constrained to end in state 0; with a correct termination and few channel errors the best survivor ends at state 0 anyway.

# Notes
This is hard-decision decoding on a binary symmetric channel. It does not accept soft metrics. ``d_{\mathrm{free}} = 5``, so every pattern of at most ``t = \lfloor (d_{\mathrm{free}} - 1)/2 \rfloor = 2`` coded-bit errors is correctable (Lin and Costello, 2nd ed., §12.1; Proakis and Salehi, 5th ed., §8.2).

The bit-error union bound on a BSC with crossover ``p`` is

```math
P_b < \sum_{d = d_{\mathrm{free}}}^{\infty} c_d P_2(d),
```

where ``c_d`` is the total information weight of the simple error events of Hamming weight ``d`` ([`distance_spectrum`](@ref)) and ``P_2(d)`` is the pairwise error probability of a weight-``d`` path, including the factor ``1/2`` on the middle binomial term when ``d`` is even (Proakis and Salehi, 5th ed., §8.2; Lin and Costello, 2nd ed., §12.3).
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
