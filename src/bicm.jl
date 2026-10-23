"""
    BICM

Teaching bit-interleaved coded modulation with iterative demapping.

Gray [`demodulate_soft`](@ref) is the equal-prior demapper. This module adds an optional rectangular bit interleaver and, for a convolutional log-MAP decoder or a short LDPC decoder, feeds coded-bit extrinsic LLRs back into the demapper. It is not 5G NR BICM: no NR rate matching, no NR LDPC base graph, no layer mapping, no scrambler, and no EXIT chart. It is not wired into the MCS / HARQ loop.

`:sova` and `:turbo` run one demapper pass and then the existing information-bit decoder. Those two public decoders do not return coded-bit extrinsic LLRs, so they are not iterated with the demapper.
"""
module BICM

using ..Modulation: ModulationScheme, BPSK, QPSK, PSK, QAM, DBPSK, DQPSK,
    modulate, demodulate_soft, constellation, bits_per_symbol
using ..Coding: conv_encode, sova_decode,
    LDPCCode, ldpc_encode, ldpc_decode_details,
    turbo_encode, turbo_decode
using ..Interleaving: block_interleave, block_deinterleave

export BICMResult, bicm_labels, bicm_psk8_id_labels, bicm_map, bicm_demap, bicm_decode

const _MAX_ITER = 32
const _G1 = 0b111
const _G2 = 0b101
const _N_STATES = 4
const _TAIL = 2

function _trellis()
    next_state = zeros(Int, _N_STATES, 2)
    output = zeros(Int, _N_STATES, 2)
    for s in 0:(_N_STATES - 1)
        for b in 0:1
            reg = (s << 1) | b
            g1 = count_ones(reg & _G1) & 1
            g2 = count_ones(reg & _G2) & 1
            next_state[s + 1, b + 1] = reg & (_N_STATES - 1)
            output[s + 1, b + 1] = (g1 << 1) | g2
        end
    end
    return next_state, output
end

const _NEXT, _OUT = _trellis()

@doc raw"""
    BICMResult

Information-bit decision from [`bicm_decode`](@ref).

# Fields
- `bits`: hard decisions, length equal to the information block. Bit 1 is an LLR strictly below 0. A tie at 0 is bit 0.
- `llr`: a-posteriori information-bit LLRs from the last decoder pass, ``\log P(u=0\mid\mathbf{y})/P(u=1\mid\mathbf{y})``. A positive value favours bit 0. For `:bcjr` this matches [`CommToolkit.Coding.bcjr_decode`](@ref) on the deinterleaved demapper LLRs when `iterations == 1` and the information-bit prior is zero. For `:ldpc` the entries are `posterior[info_pos]` from [`ldpc_decode_details`](@ref). For `:sova` and `:turbo` they are the vectors those decoders return.
- `iterations`: demapper passes actually run. Equal to the requested `iterations`.
"""
struct BICMResult
    bits::Vector{Int}
    llr::Vector{Float64}
    iterations::Int
end

function Base.show(io::IO, result::BICMResult)
    print(io, "BICMResult(k=", length(result.bits), ", iterations=", result.iterations, ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", result::BICMResult)
    return show(io, result)
end

function _reject_differential(scheme)
    if scheme isa DBPSK || scheme isa DQPSK
        throw(ArgumentError(
            "BICM demapping is coherent and memoryless; $(typeof(scheme)) is differential",
        ))
    end
    return nothing
end

function _scheme_M(scheme::ModulationScheme)
    _reject_differential(scheme)
    return 1 << bits_per_symbol(scheme)
end

@doc raw"""
    bicm_labels(scheme::ModulationScheme) -> Vector{Int}

Gray bit labels of [`constellation`](@ref)`(scheme)`, in that vector's index order.

Entry `i` is the integer `i - 1`. [`modulate`](@ref) and [`demodulate_soft`](@ref) already store the Gray point for those bits, MSB first. `DBPSK` and `DQPSK` are rejected: this receiver does not demap a differential phase.
"""
function bicm_labels(scheme::ModulationScheme)
    M = _scheme_M(scheme)
    return collect(0:(M - 1))
end

@doc raw"""
    bicm_psk8_id_labels() -> Vector{Int}

Ungerboeck set-partitioning labels for [`PSK`](@ref)`(8)`, in [`constellation`](@ref) index order.

The constellation point at angle ``2\pi m/8`` (index ``m \oplus (m \gg 1)`` inside `constellation(PSK(8))`) receives

| ``m`` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|------:|--:|--:|--:|--:|--:|--:|--:|--:|
| label | 0 | 4 | 2 | 6 | 1 | 5 | 3 | 7 |

The most-significant bit equals ``m \bmod 2``, so each value of that bit is a QPSK subset of minimum Euclidean distance ``\sqrt{2}``. Given that bit, the middle bit selects an antipodal pair. This is the set partition of Ungerboeck, IEEE Trans. Inf. Theory, 28(1):55–67, 1982. It is not Gray labelling, not a 3GPP label, and not used by [`modulate`](@ref).

Gray 8-PSK leaves little extrinsic information for a second demapper pass (Li, Chindapol, and Ritcey, IEEE Trans. Commun., 50(8):1250–1257, 2002). [`bicm_decode`](@ref) still accepts Gray labels; the iterative BER comparison in the tests uses this map because the gain is visible on a short terminated `(7, 5)` block. The labels are not claimed to be EXIT-optimal.
"""
function bicm_psk8_id_labels()
    angle_labels = (0, 4, 2, 6, 1, 5, 3, 7)
    labels = Vector{Int}(undef, 8)
    @inbounds for m in 0:7
        labels[xor(m, m >> 1) + 1] = angle_labels[m + 1]
    end
    return labels
end

function _checked_labels(labels, M::Int)
    labels === nothing && return nothing
    length(labels) == M || throw(ArgumentError(
        "labels length must be $M, got $(length(labels))",
    ))
    labs = Vector{Int}(undef, M)
    seen = falses(M)
    i = 0
    for lab in labels
        v = Int(lab)
        (0 <= v < M) || throw(ArgumentError("label $v is outside 0:$(M - 1)"))
        seen[v + 1] && throw(ArgumentError("label $v is repeated"))
        seen[v + 1] = true
        i += 1
        @inbounds labs[i] = v
    end
    return labs
end

function _is_identity_labels(labels, M::Int)
    labels === nothing && return true
    @inbounds for i in 1:M
        labels[i] == i - 1 || return false
    end
    return true
end

function _checked_rows(n_rows::Integer, n::Int)
    nr = Int(n_rows)
    nr >= 0 || throw(ArgumentError("n_rows must be ≥ 0, got $nr"))
    if nr == 0 || n == 0
        return 0
    end
    n % nr == 0 || throw(ArgumentError(
        "coded length $n is not divisible by n_rows $nr",
    ))
    return nr
end

function _interleave(x::AbstractVector, n_rows::Int)
    n_rows == 0 && return x
    return block_interleave(x, n_rows)
end

function _deinterleave(x::AbstractVector, n_rows::Int)
    n_rows == 0 && return x
    return block_deinterleave(x, n_rows)
end

function _bits_per(scheme)
    _reject_differential(scheme)
    return bits_per_symbol(scheme)
end

function _require_multiple(n_coded::Int, k::Int)
    n_coded % k == 0 || throw(ArgumentError(
        "coded length $n_coded is not a multiple of $k bits per symbol",
    ))
    return n_coded ÷ k
end

@doc raw"""
    bicm_map(info, scheme::ModulationScheme; encoder::Symbol=:conv, n_rows::Integer=0,
             labels=nothing, code=nothing, perm=nothing)

Encode, optionally bit-interleave, and modulate.

# Arguments
- `info`: information bits in ``\{0,1\}``.
- `scheme`: [`BPSK`](@ref), [`QPSK`](@ref), [`PSK`](@ref), or square [`QAM`](@ref). Differential schemes are rejected.
- `encoder`: `:conv` calls [`conv_encode`](@ref). `:ldpc` calls [`ldpc_encode`](@ref) and requires `code`. `:turbo` calls [`turbo_encode`](@ref) and requires `perm`.
- `n_rows`: passed to [`block_interleave`](@ref) on the coded bits. `0` leaves them in encoder order.
- `labels`: a permutation of `0:M-1` giving the integer bit label of each [`constellation`](@ref) point, or `nothing` for Gray ([`bicm_labels`](@ref)). [`bicm_psk8_id_labels`](@ref) is the 8-PSK set-partitioning map.

# Returns
Channel symbols in modulator order. Gray `BPSK` with `labels === nothing` is the real vector from [`modulate`](@ref). Every other path returns `Vector{ComplexF64}`.

The coded length, including a convolutional or turbo tail, must be a multiple of [`bits_per_symbol`](@ref)`(scheme)`. The terminated convolutional rate on `n` information bits is ``n/(2n+4)``, not exactly ``1/2``.
"""
function bicm_map(
    info::AbstractVector{<:Integer},
    scheme::ModulationScheme;
    encoder::Symbol=:conv,
    n_rows::Integer=0,
    labels=nothing,
    code=nothing,
    perm=nothing,
)
    k = _bits_per(scheme)
    coded = _encode(info, encoder, code, perm)
    nr = _checked_rows(n_rows, length(coded))
    _require_multiple(length(coded), k)
    bits = _interleave(coded, nr)
    return _modulate_bits(bits, scheme, labels)
end

function _encode(info, encoder::Symbol, code, perm)
    if encoder === :conv
        return conv_encode(info)
    elseif encoder === :ldpc
        code isa LDPCCode || throw(ArgumentError("encoder=:ldpc requires an LDPCCode"))
        return ldpc_encode(code, info)
    elseif encoder === :turbo
        perm === nothing && throw(ArgumentError("encoder=:turbo requires perm"))
        return turbo_encode(info, perm)
    else
        throw(ArgumentError(
            "encoder must be :conv, :ldpc, or :turbo, got $(repr(encoder))",
        ))
    end
end

function _modulate_bits(bits::AbstractVector{<:Integer}, scheme::ModulationScheme, labels)
    M = _scheme_M(scheme)
    labs = _checked_labels(labels, M)
    if _is_identity_labels(labs, M)
        return modulate(bits, scheme)
    end
    k = bits_per_symbol(scheme)
    nbits = length(bits)
    nbits % k == 0 || throw(ArgumentError(
        "bit length $nbits is not a multiple of $k",
    ))
    all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("bits must be 0 or 1"))
    inv = Vector{Int}(undef, M)
    @inbounds for i in 1:M
        inv[labs[i] + 1] = i
    end
    pts = constellation(scheme)
    nsym = nbits ÷ k
    out = Vector{ComplexF64}(undef, nsym)
    @inbounds for s in 1:nsym
        v = 0
        base = (s - 1) * k
        for j in 1:k
            v = (v << 1) | bits[base + j]
        end
        out[s] = pts[inv[v + 1]]
    end
    scheme isa BPSK && return real.(out)
    return out
end

function _check_noise(noise_var::Real)
    σ2 = float(noise_var)
    (isfinite(σ2) && σ2 > 0) || throw(ArgumentError(
        "noise_var must be finite and positive, got $noise_var",
    ))
    return σ2
end

function _check_approx(approx::Symbol)
    (approx === :exact || approx === :maxlog) || throw(ArgumentError(
        "approx must be :exact or :maxlog, got $(repr(approx))",
    ))
    return nothing
end

function _as_complex_symbols(symbols, scheme)
    if scheme isa BPSK || scheme isa PSK
        return symbols isa AbstractVector{<:Complex} ? symbols : ComplexF64.(symbols)
    end
    symbols isa AbstractVector{<:Complex} || throw(ArgumentError(
        "soft demodulation of $(typeof(scheme)) requires complex symbols",
    ))
    return symbols
end

function _as_gains(gains, nsym::Int)
    gains === nothing && return nothing
    length(gains) == nsym || throw(ArgumentError(
        "gains length must be $nsym, got $(length(gains))",
    ))
    all(isfinite, gains) || throw(ArgumentError("gains must be finite"))
    return ComplexF64.(gains)
end

function _as_apriori(apriori, nbits::Int)
    apriori === nothing && return nothing
    length(apriori) == nbits || throw(ArgumentError(
        "apriori length must be $nbits, got $(length(apriori))",
    ))
    all(isfinite, apriori) || throw(ArgumentError("apriori LLRs must be finite"))
    return Float64.(apriori)
end

function _logsumexp_bit!(metrics::Vector{Float64}, bitmat::Matrix{Int}, j::Int, bitval::Int)
    m = -Inf
    M = size(bitmat, 1)
    @inbounds for i in 1:M
        bitmat[i, j] == bitval || continue
        metrics[i] > m && (m = metrics[i])
    end
    m == -Inf && return -Inf
    m == Inf && return Inf
    s = 0.0
    @inbounds for i in 1:M
        bitmat[i, j] == bitval || continue
        s += exp(metrics[i] - m)
    end
    return m + log(s)
end

function _max_bit(metrics::Vector{Float64}, bitmat::Matrix{Int}, j::Int, bitval::Int)
    best = -Inf
    M = size(bitmat, 1)
    @inbounds for i in 1:M
        bitmat[i, j] == bitval || continue
        metrics[i] > best && (best = metrics[i])
    end
    return best
end

function _bitmat(labels::Vector{Int}, k::Int)
    M = length(labels)
    bitmat = Matrix{Int}(undef, M, k)
    @inbounds for i in 1:M
        lab = labels[i]
        for j in 1:k
            bitmat[i, j] = (lab >> (k - j)) & 1
        end
    end
    return bitmat
end

"""
Demap one block. `La === nothing` means a zero prior. Metrics use ``-|y-hs|^2/N_0`` with ``N_0 = 2`` `noise_var`, plus the bit prior of every bit except the one whose extrinsic LLR is being formed.
"""
function _demap_general(symbols, scheme, noise_var::Float64, La::Union{Nothing,Vector{Float64}},
                        approx::Symbol, labels::Vector{Int}, gains::Union{Nothing,Vector{ComplexF64}})
    k = bits_per_symbol(scheme)
    M = length(labels)
    pts = constellation(scheme)
    bitmat = _bitmat(labels, k)
    invN0 = inv(2 * noise_var)
    nsym = length(symbols)
    out = Vector{Float64}(undef, nsym * k)
    channel = Vector{Float64}(undef, M)
    metric = Vector{Float64}(undef, M)
    prior = Vector{Float64}(undef, M)
    @inbounds for s in 1:nsym
        y = ComplexF64(symbols[s])
        h = gains === nothing ? ComplexF64(1) : gains[s]
        for i in 1:M
            channel[i] = -abs2(y - h * pts[i]) * invN0
            prior[i] = 0.0
        end
        base = (s - 1) * k
        if La !== nothing
            for i in 1:M
                acc = 0.0
                for j in 1:k
                    b = bitmat[i, j]
                    acc += 0.5 * (1 - 2 * b) * La[base + j]
                end
                prior[i] = acc
            end
        end
        for j in 1:k
            for i in 1:M
                b = bitmat[i, j]
                bit_prior = La === nothing ? 0.0 : 0.5 * (1 - 2 * b) * La[base + j]
                metric[i] = channel[i] + prior[i] - bit_prior
            end
            if approx === :maxlog
                out[base + j] = _max_bit(metric, bitmat, j, 0) - _max_bit(metric, bitmat, j, 1)
            else
                out[base + j] = _logsumexp_bit!(metric, bitmat, j, 0) -
                                _logsumexp_bit!(metric, bitmat, j, 1)
            end
        end
    end
    return out
end

@doc raw"""
    bicm_demap(symbols, scheme::ModulationScheme, noise_var::Real; apriori=nothing,
               approx::Symbol=:exact, labels=nothing, gains=nothing) -> Vector{Float64}

Extrinsic bit LLRs in modulator order.

# Arguments
- `symbols`: one coherent observation per symbol. Real samples are accepted for [`BPSK`](@ref) and [`PSK`](@ref), as in [`demodulate_soft`](@ref). `QPSK` and `QAM` require complex samples.
- `scheme`: the constellation. `DBPSK` and `DQPSK` are rejected.
- `noise_var`: ``\sigma^2 = N_0/2`` of one real noise component, the same number [`demodulate_soft`](@ref) and [`CommToolkit.Channel.awgn`](@ref) use. Finite and positive.
- `apriori`: bit LLRs in modulator order, ``L_a=\log P(b=0)/P(b=1)``, length ``\texttt{length(symbols)}\times k``. `nothing` is a zero prior. The returned value does not include `apriori`: it is extrinsic.
- `approx`: `:exact` marginalises the constellation. `:maxlog` keeps the best point of each bit subset (Robertson, Villebrun, and Hoeher, Proc. IEEE ICC, 1995, pp. 1009–1013).
- `labels`: see [`bicm_map`](@ref). `nothing` is Gray.
- `gains`: optional per-symbol coefficients. The metric is ``|y - h s|^2/N_0``. `nothing` means ``h = 1``. This is ideal CSI on ``y = h x + n`` with a common noise variance. It is not the zero-forced output of `rayleigh_channel(..., true)`, whose noise variance is ``\sigma^2/|h|^2`` and is no longer a single `noise_var`.

# Returns
Length ``\texttt{length(symbols)}\times k``. A positive value favours bit 0.

With Gray labels, a zero prior, and `gains === nothing`, the vector is [`demodulate_soft`](@ref) of the same arguments, including `:exact` and `:maxlog`.

The bit prior enters as ``\tfrac12(1-2b)L_a``. The term that depends only on the bit being computed is left out, which is the extrinsic split of ten Brink, Speidel, and Yan, Proc. IEEE GLOBECOM, 1998, pp. 579–584, and of Hagenauer, Offer, and Papke, IEEE Trans. Inf. Theory, 42(2):429–445, 1996. For BPSK and for Gray QPSK the extrinsic LLR does not depend on `apriori`: each bit is already a sufficient statistic. Iteration needs a labelling or a product constellation in which bits that share a symbol are coupled. Square Gray QAM couples the bits that share an axis and does not couple in-phase bits with quadrature bits.
"""
function bicm_demap(
    symbols::AbstractVector{<:Number},
    scheme::ModulationScheme,
    noise_var::Real;
    apriori=nothing,
    approx::Symbol=:exact,
    labels=nothing,
    gains=nothing,
)
    _reject_differential(scheme)
    σ2 = _check_noise(noise_var)
    _check_approx(approx)
    all(isfinite, symbols) || throw(ArgumentError("symbols must be finite"))
    M = _scheme_M(scheme)
    labs = _checked_labels(labels, M)
    k = bits_per_symbol(scheme)
    nsym = length(symbols)
    nbits = nsym * k
    La = _as_apriori(apriori, nbits)
    H = _as_gains(gains, nsym)
    if _is_identity_labels(labs, M) && La === nothing && H === nothing
        return demodulate_soft(symbols, scheme, σ2; approx=approx)
    end
    ys = _as_complex_symbols(symbols, scheme)
    use = labs === nothing ? collect(0:(M - 1)) : labs
    return _demap_general(ys, scheme, σ2, La, approx, use, H)
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

"""
Coded-bit extrinsic log-MAP on the terminated (7, 5) trellis.

Information-bit APP uses the same forward–backward recursion, the same branch metric, and the same termination as `bcjr_decode` with a zero information-bit prior. Coded-bit extrinsic LLRs subtract the channel LLR. Tail inputs are fixed at 0.
"""
function _conv_siso(llr::Vector{Float64})
    n_coded = length(llr)
    iseven(n_coded) || throw(ArgumentError("LLR length must be even, got $n_coded"))
    n_sym = n_coded ÷ 2
    n_bits = n_sym - _TAIL
    n_bits >= 0 || throw(ArgumentError("coded sequence is shorter than the termination tail"))
    all(isfinite, llr) || throw(ArgumentError("LLRs must be finite"))

    logα = fill(-Inf, n_sym + 1, _N_STATES)
    logβ = fill(-Inf, n_sym + 1, _N_STATES)
    logα[1, 1] = 0.0
    logβ[n_sym + 1, 1] = 0.0

    @inbounds for t in 1:n_sym
        L0 = llr[2t - 1]
        L1 = llr[2t]
        tail = t > n_bits
        for s in 0:(_N_STATES - 1)
            isfinite(logα[t, s + 1]) || continue
            for b in 0:1
                tail && b == 1 && continue
                out = _OUT[s + 1, b + 1]
                c0 = (out >> 1) & 1
                c1 = out & 1
                γ = 0.5 * ((1 - 2 * c0) * L0 + (1 - 2 * c1) * L1)
                ns = _NEXT[s + 1, b + 1]
                logα[t + 1, ns + 1] = _logaddexp(logα[t + 1, ns + 1], logα[t, s + 1] + γ)
            end
        end
        _normalise_log!(@view logα[t + 1, :])
    end

    @inbounds for t in n_sym:-1:1
        L0 = llr[2t - 1]
        L1 = llr[2t]
        tail = t > n_bits
        for s in 0:(_N_STATES - 1)
            for b in 0:1
                tail && b == 1 && continue
                out = _OUT[s + 1, b + 1]
                c0 = (out >> 1) & 1
                c1 = out & 1
                γ = 0.5 * ((1 - 2 * c0) * L0 + (1 - 2 * c1) * L1)
                ns = _NEXT[s + 1, b + 1]
                logβ[t, s + 1] = _logaddexp(logβ[t, s + 1], logβ[t + 1, ns + 1] + γ)
            end
        end
        _normalise_log!(@view logβ[t, :])
    end

    info = Vector{Float64}(undef, n_bits)
    coded_ext = Vector{Float64}(undef, n_coded)
    u0 = Vector{Float64}(undef, _N_STATES * 2)
    u1 = similar(u0)
    c00 = similar(u0)
    c01 = similar(u0)
    c10 = similar(u0)
    c11 = similar(u0)
    @inbounds for t in 1:n_sym
        L0 = llr[2t - 1]
        L1 = llr[2t]
        tail = t > n_bits
        n_u0 = 0
        n_u1 = 0
        n_c00 = 0
        n_c01 = 0
        n_c10 = 0
        n_c11 = 0
        for s in 0:(_N_STATES - 1)
            αs = logα[t, s + 1]
            for b in 0:1
                tail && b == 1 && continue
                out = _OUT[s + 1, b + 1]
                c0 = (out >> 1) & 1
                c1 = out & 1
                γ = 0.5 * ((1 - 2 * c0) * L0 + (1 - 2 * c1) * L1)
                ns = _NEXT[s + 1, b + 1]
                v = αs + γ + logβ[t + 1, ns + 1]
                if b == 0
                    n_u0 += 1
                    u0[n_u0] = v
                else
                    n_u1 += 1
                    u1[n_u1] = v
                end
                if c0 == 0
                    n_c00 += 1
                    c00[n_c00] = v
                else
                    n_c01 += 1
                    c01[n_c01] = v
                end
                if c1 == 0
                    n_c10 += 1
                    c10[n_c10] = v
                else
                    n_c11 += 1
                    c11[n_c11] = v
                end
            end
        end
        if !tail
            info[t] = _logsumexp(u0, n_u0) - _logsumexp(u1, n_u1)
        end
        coded_ext[2t - 1] = _logsumexp(c00, n_c00) - _logsumexp(c01, n_c01) - L0
        coded_ext[2t] = _logsumexp(c10, n_c10) - _logsumexp(c11, n_c11) - L1
    end
    return info, coded_ext
end

function _clip_feedback(x::Vector{Float64}, clip::Float64)
    y = Vector{Float64}(undef, length(x))
    @inbounds for i in eachindex(x)
        v = x[i]
        isnan(v) && throw(ArgumentError("decoder extrinsic LLR is NaN"))
        y[i] = v > clip ? clip : (v < -clip ? -clip : v)
    end
    return y
end

function _hard_bits(llr::Vector{Float64})
    bits = Vector{Int}(undef, length(llr))
    @inbounds for i in eachindex(llr)
        bits[i] = llr[i] < 0 ? 1 : 0
    end
    return bits
end

function _check_iterations(iterations::Integer, decoder::Symbol)
    it = Int(iterations)
    (1 <= it <= _MAX_ITER) || throw(ArgumentError(
        "iterations must be from 1 to $_MAX_ITER, got $it",
    ))
    if (decoder === :sova || decoder === :turbo) && it != 1
        throw(ArgumentError(
            "decoder=$(repr(decoder)) returns information-bit LLRs only; " *
            "coded-bit extrinsic feedback is not available, so iterations must be 1",
        ))
    end
    return it
end

function _expected_coded_length(decoder::Symbol, n_info::Int, code, perm)
    if decoder === :bcjr || decoder === :sova
        n_info >= 0 || throw(ArgumentError("n_info must be ≥ 0, got $n_info"))
        return 2 * (n_info + _TAIL)
    elseif decoder === :ldpc
        code isa LDPCCode || throw(ArgumentError("decoder=:ldpc requires an LDPCCode"))
        n_info == code.k || throw(ArgumentError(
            "n_info must be $(code.k) for this LDPC code, got $n_info",
        ))
        return code.n
    elseif decoder === :turbo
        perm === nothing && throw(ArgumentError("decoder=:turbo requires perm"))
        n_info >= 2 || throw(ArgumentError("turbo n_info must be ≥ 2, got $n_info"))
        length(perm) == n_info || throw(ArgumentError(
            "perm length must be $n_info, got $(length(perm))",
        ))
        return 3 * n_info + 8
    else
        throw(ArgumentError(
            "decoder must be :bcjr, :sova, :ldpc, or :turbo, got $(repr(decoder))",
        ))
    end
end

@doc raw"""
    bicm_decode(symbols, scheme::ModulationScheme, noise_var::Real, n_info::Integer;
                decoder::Symbol=:bcjr, iterations::Integer=4, n_rows::Integer=0,
                approx::Symbol=:exact, labels=nothing, gains=nothing,
                code=nothing, ldpc_iterations::Integer=8,
                algorithm::Symbol=:normalized_min_sum,
                perm=nothing, turbo_iterations::Integer=4,
                llr_clip::Real=50) -> BICMResult

Iterative demapping and decoding of one block from [`bicm_map`](@ref).

# Arguments
- `symbols`, `scheme`, `noise_var`, `approx`, `labels`, `gains`: the same contract as [`bicm_demap`](@ref). `n_rows` and `labels` must be the values used to map the block.
- `n_info`: number of information bits. For `:bcjr` and `:sova` the coded length is ``2(n_{\mathrm{info}}+2)``. For `:ldpc` it is `code.n` and `n_info` must equal `code.k`. For `:turbo` it is ``3 n_{\mathrm{info}}+8``.
- `decoder`:
  - `:bcjr` — coded-bit log-MAP on the terminated `(7, 5)` trellis. The branch metric and the termination match [`CommToolkit.Coding.bcjr_decode`](@ref). Information-bit LLRs at `iterations == 1` match `bcjr_decode` on the deinterleaved demapper output. Coded-bit extrinsic LLRs, clipped and interleaved, are the next demapper prior.
  - `:ldpc` — [`ldpc_decode_details`](@ref) with `early_stop=false` and `max_iter=ldpc_iterations`. Extrinsic feedback is `posterior - L_channel`, where `L_channel` is the demapper extrinsic and not the previous prior. `algorithm` is forwarded.
  - `:sova` — one pass: deinterleaved demapper LLRs go to [`sova_decode`](@ref). `iterations` must be 1.
  - `:turbo` — one pass: those LLRs go to [`turbo_decode`](@ref) with `perm` and `turbo_iterations`. `iterations` must be 1.
- `iterations`: demapper–decoder passes, from 1 to 32. The default 4 is a teaching choice. One pass is BICM without feedback.
- `llr_clip`: positive finite cap on the magnitude of coded-bit extrinsic LLRs before they re-enter the demapper. The default 50 keeps the exact demapper's exponentials in a comfortable range. Information-bit LLRs returned in the result are not clipped.

# Schedule
The demapper writes an extrinsic LLR for every modulated bit. [`block_deinterleave`](@ref) restores encoder order when `n_rows > 0`. The decoder treats that vector as its channel observation. Feedback to the demapper is interleaved again with [`block_interleave`](@ref). The information-bit prior inside the convolutional decoder stays zero: the prior lives on the coded bits, which are the bits the modulator used.

Hard decisions are `llr .< 0`. This does not approach coded-modulation capacity, and Gray labelling is not expected to move much between pass 1 and pass 4. The set-partitioning 8-PSK labels are the map on which a short AWGN test separates those two passes.

Not implemented: puncturing inside the loop, EXIT charts, a SOVA or turbo coded-bit SISO, differential BICM, and the MCS table.
"""
function bicm_decode(
    symbols::AbstractVector{<:Number},
    scheme::ModulationScheme,
    noise_var::Real,
    n_info::Integer;
    decoder::Symbol=:bcjr,
    iterations::Integer=4,
    n_rows::Integer=0,
    approx::Symbol=:exact,
    labels=nothing,
    gains=nothing,
    code=nothing,
    ldpc_iterations::Integer=8,
    algorithm::Symbol=:normalized_min_sum,
    perm=nothing,
    turbo_iterations::Integer=4,
    llr_clip::Real=50,
)
    k = _bits_per(scheme)
    n_info = Int(n_info)
    it = _check_iterations(iterations, decoder)
    clip = float(llr_clip)
    (isfinite(clip) && clip > 0) || throw(ArgumentError(
        "llr_clip must be finite and positive, got $llr_clip",
    ))
    n_coded = length(symbols) * k
    expect = _expected_coded_length(decoder, n_info, code, perm)
    n_coded == expect || throw(ArgumentError(
        "symbol block has $n_coded coded bits; decoder=$(repr(decoder)) expects $expect",
    ))
    nr = _checked_rows(n_rows, n_coded)
    _require_multiple(n_coded, k)

    apriori = nothing
    info = zeros(Float64, n_info)
    for pass in 1:it
        ext_mod = bicm_demap(symbols, scheme, noise_var;
                             apriori=apriori, approx=approx, labels=labels, gains=gains)
        ext_enc = _deinterleave(ext_mod, nr)
        feedback = nothing
        if decoder === :bcjr
            info, coded_ext = _conv_siso(ext_enc)
            feedback = coded_ext
        elseif decoder === :sova
            info = sova_decode(ext_enc)
        elseif decoder === :ldpc
            detail = ldpc_decode_details(code, ext_enc;
                                         algorithm=algorithm,
                                         max_iter=ldpc_iterations,
                                         early_stop=false)
            info = detail.posterior[code.info_pos]
            feedback = detail.posterior .- ext_enc
        elseif decoder === :turbo
            info = turbo_decode(ext_enc, perm; iterations=turbo_iterations)
        end
        if feedback !== nothing && pass < it
            apriori = _interleave(_clip_feedback(feedback, clip), nr)
        end
    end
    return BICMResult(_hard_bits(info), info, it)
end

end # module BICM
