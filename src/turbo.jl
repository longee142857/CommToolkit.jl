# Short parallel-concatenated turbo code. Not a 3GPP / LTE turbo encoder:
# memory-2 RSC components, a rectangular or stride interleaver, and a hard cap
# on the information length and the iteration count.

const _TURBO_M = 2
const _TURBO_MAX_K = 1024
const _TURBO_MAX_ITER = 32

function _rsc_tables()
    next = zeros(Int, 4, 2)
    par = zeros(Int, 4, 2)
    for s in 0:3
        d = s & 1
        d2 = (s >> 1) & 1
        for u in 0:1
            a = xor(u, d, d2)
            next[s + 1, u + 1] = a | (d << 1)
            par[s + 1, u + 1] = xor(a, d2)
        end
    end
    return next, par
end

const _RSC_NEXT, _RSC_PAR = _rsc_tables()

function _rsc_tail_input(state::Int)
    return (state & 1) ⊻ ((state >> 1) & 1)
end

function _rsc_encode(bits::Vector{Int})
    k = length(bits)
    parity = Vector{Int}(undef, k + _TURBO_M)
    tail = Vector{Int}(undef, _TURBO_M)
    state = 0
    @inbounds for i in 1:k
        u = bits[i]
        parity[i] = _RSC_PAR[state + 1, u + 1]
        state = _RSC_NEXT[state + 1, u + 1]
    end
    @inbounds for j in 1:_TURBO_M
        u = _rsc_tail_input(state)
        tail[j] = u
        parity[k + j] = _RSC_PAR[state + 1, u + 1]
        state = _RSC_NEXT[state + 1, u + 1]
    end
    state == 0 || error("internal: RSC tail did not return to state 0")
    return parity, tail
end

function _turbo_perm(perm::AbstractVector{<:Integer}, k::Int)
    length(perm) == k || throw(ArgumentError(
        "interleaver length must be $k, got $(length(perm))",
    ))
    seen = falses(k)
    out = Vector{Int}(undef, k)
    i = 0
    for p in perm
        pi = Int(p)
        (1 <= pi <= k) || throw(ArgumentError("interleaver entry $pi is outside 1:$k"))
        seen[pi] && throw(ArgumentError("interleaver repeats index $pi"))
        seen[pi] = true
        i += 1
        @inbounds out[i] = pi
    end
    return out
end

@doc raw"""
    turbo_interleaver(k::Integer) -> Vector{Int}

Length-`k` interleaver for [`turbo_encode`](@ref). Entries are a permutation of `1:k`: encoder 2 reads information bit `perm[i]` at its time `i`.

`k` is an integer from 2 to 1024.

- If `k` has a factor `r` with `2 ≤ r ≤ k/2` and `k/r ≥ 2`, the permutation is the rectangular block interleaver (`block_interleave`): write `1:k` by rows and read by columns. The row count is the factor closest to a square, and the smaller factor wins a tie. It is the same map as `block_interleave(collect(1:k), n_rows, n_cols)`, not a second interleaver family.
- Otherwise (`k` prime, and the length-2 case) the permutation is the stride ``\pi(i)=\bigl(s(i-1)\bmod k\bigr)+1`` with `s` the integer closest to ``\sqrt{k}`` from below that is coprime to `k`. Length 2 returns `[2, 1]`.

This is not a 3GPP quadratic permutation polynomial, and it is not required to be a derangement. The caller can pass any other permutation of `1:k` to `turbo_encode` and `turbo_decode`.
"""
function turbo_interleaver(k::Integer)
    n = Int(k)
    (2 <= n <= _TURBO_MAX_K) || throw(ArgumentError(
        "turbo interleaver length must be from 2 to $_TURBO_MAX_K, got $n",
    ))
    if n == 2
        return [2, 1]
    end
    best_r = 0
    best_gap = typemax(Int)
    for r in 2:(n ÷ 2)
        n % r == 0 || continue
        c = n ÷ r
        c >= 2 || continue
        gap = abs(r - c)
        if gap < best_gap || (gap == best_gap && (best_r == 0 || r < best_r))
            best_r = r
            best_gap = gap
        end
    end
    if best_r != 0
        nr = best_r
        nc = n ÷ best_r
        perm = Vector{Int}(undef, n)
        @inbounds for c in 1:nc
            for r in 1:nr
                perm[(c - 1) * nr + r] = (r - 1) * nc + c
            end
        end
        return perm
    end
    s = isqrt(n)
    s = s < 1 ? 1 : (s > n - 1 ? n - 1 : s)
    while s > 1 && gcd(s, n) != 1
        s -= 1
    end
    if s <= 1
        s = n - 1
        while s > 1 && gcd(s, n) != 1
            s -= 1
        end
    end
    return [mod((i - 1) * s, n) + 1 for i in 1:n]
end

@doc raw"""
    turbo_encode(bits::AbstractVector{<:Integer}, perm::AbstractVector{<:Integer}) -> Vector{Int}

Parallel concatenation of two tail-terminated recursive systematic convolutional encoders.

Each component has memory ``m = 2`` and generators ``(1,\ 5/7)`` in octal: feedback ``1+D+D^2``, parity ``1+D^2``. That is the recursive systematic form of the feedforward `(7, 5)` code used by [`conv_encode`](@ref). It is not that feedforward encoder, and it is not the memory-4 LTE / NR turbo component.

`bits` has length ``k\in[2,1024]``, entries in ``\{0,1\}``. `perm` is a permutation of `1:k` in the sense of [`turbo_interleaver`](@ref): the second encoder reads `bits[perm]`.

Register convention, so the map can be reproduced without this source file. State bits are `(D, D^2)` with the least significant bit equal to the most recent value shifted into the register. The register input is ``a = u + D + D^2`` over GF(2). The parity bit is ``a + D^2``. The systematic bit is ``u``. Termination drives ``a = 0`` for two steps; those two input bits are the tail, and they are not part of `bits` and not part of `perm`.

# Returns
One vector of length ``3k+8 = 3k+4m``, laid out as

1. `k` systematic information bits,
2. `k+m` parity bits of encoder 1, including its tail parity,
3. `m` systematic tail bits of encoder 1,
4. `k+m` parity bits of encoder 2, including its tail parity,
5. `m` systematic tail bits of encoder 2.

The information bits are transmitted once. Nominal rate ``k/(3k+8)`` approaches ``1/3`` and is not exactly ``1/3``. Encoder 2's systematic symbols are the interleaved information bits already stored in part 1; they are not repeated.

# References
Berrou, Glavieux, and Thitimajshima, Proc. IEEE ICC, 1993, pp. 1064–1070, for the parallel concatenation. The polynomials, the interleaver, the tail, and the block length here are a teaching choice, not that paper's memory-4 code and not 3GPP TS 36.212.
"""
function turbo_encode(bits::AbstractVector{<:Integer}, perm::AbstractVector{<:Integer})
    info = _binary_bits(bits)
    k = length(info)
    (2 <= k <= _TURBO_MAX_K) || throw(ArgumentError(
        "turbo information length must be from 2 to $_TURBO_MAX_K, got $k",
    ))
    π = _turbo_perm(perm, k)
    p1, t1 = _rsc_encode(info)
    p2, t2 = _rsc_encode(info[π])
    return vcat(info, p1, t1, p2, t2)
end

function _turbo_split(llr::AbstractVector{<:AbstractFloat}, k::Int)
    m = _TURBO_M
    need = 3k + 4m
    n = length(llr)
    n == need || throw(ArgumentError("turbo LLR length must be $need, got $n"))
    all(isfinite, llr) || throw(ArgumentError("LLRs must be finite"))
    Ls = Float64.(@view llr[1:k])
    Lp1 = Float64.(@view llr[(k + 1):(2k + m)])
    Lt1 = Float64.(@view llr[(2k + m + 1):(2k + 2m)])
    Lp2 = Float64.(@view llr[(2k + 2m + 1):(3k + 3m)])
    Lt2 = Float64.(@view llr[(3k + 3m + 1):(3k + 4m)])
    return Ls, Lp1, Lt1, Lp2, Lt2
end

function _rsc_bcjr(sys::Vector{Float64}, par::Vector{Float64}, La::Vector{Float64})
    k = length(La)
    m = _TURBO_M
    n_sym = k + m
    n_states = 4
    logα = fill(-Inf, n_sym + 1, n_states)
    logβ = fill(-Inf, n_sym + 1, n_states)
    logα[1, 1] = 0.0
    logβ[n_sym + 1, 1] = 0.0

    @inbounds for t in 1:n_sym
        tail = t > k
        for s in 0:(n_states - 1)
            isfinite(logα[t, s + 1]) || continue
            u0 = tail ? _rsc_tail_input(s) : 0
            u1 = tail ? u0 : 1
            for u in (tail ? (u0,) : (u0, u1))
                ns = _RSC_NEXT[s + 1, u + 1]
                pbit = _RSC_PAR[s + 1, u + 1]
                prior = tail ? 0.0 : La[t]
                γ = 0.5 * ((1 - 2 * u) * (sys[t] + prior) + (1 - 2 * pbit) * par[t])
                logα[t + 1, ns + 1] = _logaddexp(logα[t + 1, ns + 1], logα[t, s + 1] + γ)
            end
        end
        _normalise_log!(@view logα[t + 1, :])
    end

    @inbounds for t in n_sym:-1:1
        tail = t > k
        for s in 0:(n_states - 1)
            u0 = tail ? _rsc_tail_input(s) : 0
            u1 = tail ? u0 : 1
            for u in (tail ? (u0,) : (u0, u1))
                ns = _RSC_NEXT[s + 1, u + 1]
                pbit = _RSC_PAR[s + 1, u + 1]
                prior = tail ? 0.0 : La[t]
                γ = 0.5 * ((1 - 2 * u) * (sys[t] + prior) + (1 - 2 * pbit) * par[t])
                logβ[t, s + 1] = _logaddexp(logβ[t, s + 1], logβ[t + 1, ns + 1] + γ)
            end
        end
        _normalise_log!(@view logβ[t, :])
    end

    Le = Vector{Float64}(undef, k)
    Lapp = Vector{Float64}(undef, k)
    acc0 = Vector{Float64}(undef, n_states * 2)
    acc1 = similar(acc0)
    @inbounds for t in 1:k
        n0 = 0
        n1 = 0
        for s in 0:(n_states - 1)
            αs = logα[t, s + 1]
            for u in 0:1
                ns = _RSC_NEXT[s + 1, u + 1]
                pbit = _RSC_PAR[s + 1, u + 1]
                γ = 0.5 * ((1 - 2 * u) * (sys[t] + La[t]) + (1 - 2 * pbit) * par[t])
                v = αs + γ + logβ[t + 1, ns + 1]
                if u == 0
                    n0 += 1
                    acc0[n0] = v
                else
                    n1 += 1
                    acc1[n1] = v
                end
            end
        end
        Lapp[t] = _logsumexp(acc0, n0) - _logsumexp(acc1, n1)
        Le[t] = Lapp[t] - La[t] - sys[t]
    end
    return Le, Lapp
end

@doc raw"""
    turbo_decode(llr::AbstractVector{<:AbstractFloat}, perm::AbstractVector{<:Integer}; iterations::Integer=6) -> Vector{Float64}

Iterative log-BCJR decoding of a block from [`turbo_encode`](@ref).

# Arguments
- `llr`: finite channel LLRs in the layout returned by `turbo_encode`, ``L=\log P(c=0\mid y)/P(c=1\mid y)``. For BPSK with ``E_s=1`` this is ``2y/\sigma^2``, the same convention as [`bcjr_decode`](@ref).
- `perm`: the interleaver that was passed to `turbo_encode`.
- `iterations`: how many times decoder 2's extrinsic LLRs are fed back as decoder 1's a priori LLRs. An integer from 1 to 32. The default 6 is a teaching choice, not a standard iteration count (3GPP does not define one for this code, and this code is not the 3GPP encoder).

# Returns
A-posteriori LLRs of the `k` information bits, tail removed. A positive value favours bit 0. Hard decisions are `.< 0`.

# Schedule
Each iteration runs the encoder-1 log-BCJR, interleaves its extrinsic LLRs into the encoder-2 log-BCJR, and deinterleaves encoder 2's extrinsic LLRs back to encoder 1. The returned vector is a further encoder-1 a-posteriori pass that consumes the latest feedback. Extrinsic LLRs follow the systematic split of Hagenauer, Offer, and Papke, IEEE Trans. Inf. Theory 42(2):429–445, 1996:

```math
L_e = L_{\mathrm{app}} - L_a - L_{\mathrm{systematic}}.
```

The component decoder is this log-BCJR on the RSC trellis, not [`sova_decode`](@ref) and not [`bcjr_decode`](@ref). Both of those decode the feedforward `(7, 5)` code. One iteration is already a complete exchange; it is not maximum-likelihood decoding of the concatenated code. More iterations are not monotone on every frame.

`k` is inferred from `perm` and must lie in ``[2,1024]``. The LLR vector must have length ``3k+8``.
"""
function turbo_decode(
    llr::AbstractVector{<:AbstractFloat},
    perm::AbstractVector{<:Integer};
    iterations::Integer=6,
)
    k = length(perm)
    (2 <= k <= _TURBO_MAX_K) || throw(ArgumentError(
        "turbo information length must be from 2 to $_TURBO_MAX_K, got $k",
    ))
    it = Int(iterations)
    (1 <= it <= _TURBO_MAX_ITER) || throw(ArgumentError(
        "turbo iterations must be from 1 to $_TURBO_MAX_ITER, got $it",
    ))
    π = _turbo_perm(perm, k)
    Ls, Lp1, Lt1, Lp2, Lt2 = _turbo_split(llr, k)
    sys1 = vcat(Ls, Lt1)
    sys2 = vcat(Ls[π], Lt2)
    La = zeros(Float64, k)
    Le2_as_La = Vector{Float64}(undef, k)
    for _ in 1:it
        Le1, _ = _rsc_bcjr(sys1, Lp1, La)
        Le2, _ = _rsc_bcjr(sys2, Lp2, Le1[π])
        @inbounds for i in eachindex(π)
            Le2_as_La[π[i]] = Le2[i]
        end
        La, Le2_as_La = Le2_as_La, La
    end
    _, Lapp = _rsc_bcjr(sys1, Lp1, La)
    return Lapp
end
