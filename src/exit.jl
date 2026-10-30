"""
    EXIT

Teaching EXIT-chart utilities: mutual-information measurement of soft bit
information and transfer characteristics of the existing demapper and decoders.

The measured quantity is the mutual information ``I(B;L)`` between a bit
``B\\in\\{0,1\\}`` and its soft value ``L``, in bits per symbol. All LLRs use the
package convention ``L=\\log P(b=0)/P(b=1)``: a positive value favours bit 0,
and extrinsic values are posterior minus prior (Hagenauer, Offer, and Papke,
IEEE Trans. Inf. Theory, 42(2):429–445, 1996). Two Monte Carlo estimators are
provided: the histogram estimator of ten Brink (IEEE Trans. Commun.,
49(10):1727–1737, 2001) and the histogram-free soft estimator of Hagenauer
(Proc. 12th EUSIPCO, Vienna, 2004, pp. 1541–1548). The Gaussian-symmetric
``J(\\sigma)`` curve and its inverse (ten Brink, Speidel, and Yan, IEEE Trans.
Commun., 52(4):670–678, 2004) are an optional, documented approximation — not
the default estimator — used here mainly to synthesise a-priori LLRs.

[`exit_demap_transfer`](@ref) measures the extrinsic transfer characteristic
``I_E = T(I_A)`` of [`CommToolkit.BICM.bicm_demap`](@ref) on an AWGN channel.
[`exit_decoder_transfer`](@ref) measures the information-bit transfer
characteristic of the terminated ``(7, 5)`` log-BCJR decoder
([`CommToolkit.Coding.bcjr_decode`](@ref)) or of a short flooding LDPC decoder
([`CommToolkit.Coding.ldpc_decode_details`](@ref)); those decoders are not
rewritten here. Synthetic a-priori LLRs are consistent Gaussian samples
``L_a = (\\sigma_A^2/2)(1-2b) + \\sigma_A n``, the standard EXIT measurement
assumption; its nature is stated on each function.

This is a teaching measurement surface. It is not NR LDPC EXIT production
tooling: there is no density evolution, no degree-distribution optimisation,
no 3GPP code or rate matching, and no trajectory tracking of a running
iterative receiver. [`plot_exit_chart`](@ref) draws the two measured curves
after `using Plots`. Decoder curves call `bcjr_decode` and
`ldpc_decode_details` only. As of 1.2.2, `sova_decode(...; return_coded=true)`
can return a teaching coded-bit extrinsic (a competing-path gap, not MAP);
this module does not measure it.
"""
module EXIT

using ..Metrics: snr_db_to_linear
using ..Modulation: ModulationScheme, BPSK, QPSK, PSK, QAM, DBPSK, DQPSK,
    modulate, constellation, bits_per_symbol
using ..Coding: conv_encode, bcjr_decode, LDPCCode, ldpc_encode, ldpc_decode_details
using ..BICM: bicm_demap

export EXITCurve, exit_mi_soft, exit_mi_hist, exit_j, exit_jinv,
    exit_demap_transfer, exit_decoder_transfer, plot_exit_chart

# Fixed teaching seed: the same constant the test suite uses. The stream below
# is a local xorshift64* + Box–Muller generator, identical in construction to
# the one in test/bicm.jl, so a given `seed` reproduces the same samples on
# every Julia version (Julia's global `randn` stream is not version-stable).
const _DEFAULT_SEED = 0x0c0de1162
const _MULT = 0x2545f4914f6cdd1d

mutable struct _ExitRng
    state::UInt64
    spare::Float64
    has_spare::Bool
end

function _ExitRng(seed::Integer)
    s = UInt64(seed)
    s == 0 && (s = 0x9e3779b97f4a7c15) # xorshift degenerates at state 0
    return _ExitRng(s, 0.0, false)
end

function _next_u64(rng::_ExitRng)
    z = rng.state
    z = xor(z, z >> 12)
    z = xor(z, z << 25)
    z = xor(z, z >> 27)
    rng.state = z
    return z * _MULT
end

_next_u01(rng::_ExitRng) = (_next_u64(rng) >> 11) * 2.0^-53

_next_bit(rng::_ExitRng) = Int(_next_u64(rng) >> 63)

function _next_normal(rng::_ExitRng)
    if rng.has_spare
        rng.has_spare = false
        return rng.spare
    end
    u1 = max(_next_u01(rng), 1e-300)
    u2 = _next_u01(rng)
    r = sqrt(-2.0 * log(u1))
    θ = 2.0π * u2
    rng.spare = r * sin(θ)
    rng.has_spare = true
    return r * cos(θ)
end

_bits(rng::_ExitRng, n::Int) = [_next_bit(rng) for _ in 1:n]

# log2(1 + exp(x)), stable for large |x| and for x = ±Inf.
function _log2_1pexp(x::Float64)
    sp = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))
    return sp / log(2.0)
end

function _check_llr_bits(llr, bits)
    n = length(llr)
    n == length(bits) || throw(ArgumentError(
        "llr and bits must have the same length, got $n and $(length(bits))",
    ))
    n >= 1 || throw(ArgumentError("at least one sample is required"))
    all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("bits must be 0 or 1"))
    any(isnan, llr) && throw(ArgumentError("LLRs must not be NaN"))
    return n
end

@doc raw"""
    exit_mi_soft(llr, bits) -> Float64

Monte Carlo mutual information between bits and their soft LLRs, without a histogram:

```math
\hat{I} = 1 - \frac{1}{N}\sum_{j=1}^{N}
\log_2\!\bigl(1 + e^{-(1-2b_j)\,L_j}\bigr).
```

# Arguments
- `llr`: soft values in the package convention, positive towards bit 0. `NaN` is rejected. ``\pm\infty`` is accepted: a consistent infinite LLR contributes a full noiseless sample, and a contradicting one drives the estimate to ``-\infty``, which is what an infinitely confident wrong LLR means.
- `bits`: the true bits in ``\{0,1\}``. They are known to the measurement, exactly as in every published EXIT measurement on a simulated link.

# Assumptions
The estimator converges to the true ``I(B;L)`` when the bits are equiprobable
and ``L`` is the *consistent* (true a-posteriori) LLR of the bit — this is the
identity ``I(B;L) = 1 - \mathbb{E}[\log_2(1+e^{-(1-2B)L})]`` evaluated on
samples. For inconsistent soft values (for example min-sum magnitudes, or
clipped LLRs) the same formula is still the standard EXIT measurement, but it
is then an estimate whose bias relative to ``I(B;L)`` is not corrected; use
[`exit_mi_hist`](@ref) to cross-check with the distribution-free definition.
A single bit class is accepted: measuring on the all-zero codeword is standard
practice under the channel symmetry of the decoders in this package.

# References
Hagenauer, Proc. 12th EUSIPCO, Vienna, 2004, pp. 1541–1548, eq. for the
histogram-free measurement of ``I_E``; the underlying identity goes back to the
L-value decomposition of Hagenauer, Offer, and Papke, IEEE Trans. Inf. Theory,
42(2):429–445, 1996.
"""
function exit_mi_soft(llr::AbstractVector{<:Real}, bits::AbstractVector{<:Integer})
    n = _check_llr_bits(llr, bits)
    acc = 0.0
    @inbounds for j in 1:n
        s = 1 - 2 * bits[j]
        acc += _log2_1pexp(-s * Float64(llr[j]))
    end
    return 1.0 - acc / n
end

@doc raw"""
    exit_mi_hist(llr, bits; nbins::Integer=257, clip::Real=50.0) -> Float64

Histogram mutual information between bits and their soft LLRs, the original EXIT measurement:

```math
\hat{I} = \frac12\sum_{b=0}^{1}\sum_{i}
\hat p_i(b)\,\log_2\frac{2\,\hat p_i(b)}{\hat p_i(0)+\hat p_i(1)},
```

where ``\hat p_i(b)`` is the empirical probability that an LLR of a bit with
value ``b`` falls in bin ``i``. With equiprobable bits this is the discrete
form of ten Brink's integral ``I_E = \frac12\sum_c\int p(\xi\mid c)\,
\mathrm{ld}\,\frac{2p(\xi\mid c)}{p(\xi\mid{+1})+p(\xi\mid{-1})}\,d\xi``.

# Arguments
- `llr`, `bits`: as in [`exit_mi_soft`](@ref). Both bit values must be present:
  the per-class densities are undefined otherwise (the all-zero-codeword
  shortcut belongs to [`exit_mi_soft`](@ref)).
- `nbins`: number of bins on ``[-\texttt{clip}, \texttt{clip}]``, at least 4.
- `clip`: finite positive range limit. Samples outside the range, including
  ``\pm\infty``, accumulate in the edge bins; the information they carry
  beyond `clip` is discarded. `NaN` is rejected.

# Notes
No consistency assumption is made: the estimator uses only the two conditional
histograms. Binning is a bias source — too few bins merge discriminative
detail (bias downward), too many bins amplify sampling noise (bias upward) —
so treat the result as an approximation and check stability against `nbins`.
For consistent LLRs and many samples it agrees with [`exit_mi_soft`](@ref),
and both agree with the BIAWGN capacity on a BPSK/AWGN LLR stream.

# References
ten Brink, IEEE Trans. Commun., 49(10):1727–1737, 2001.
"""
function exit_mi_hist(llr::AbstractVector{<:Real}, bits::AbstractVector{<:Integer};
                      nbins::Integer=257, clip::Real=50.0)
    n = _check_llr_bits(llr, bits)
    nb = Int(nbins)
    nb >= 4 || throw(ArgumentError("nbins must be ≥ 4, got $nb"))
    c = float(clip)
    (isfinite(c) && c > 0) || throw(ArgumentError(
        "clip must be finite and positive, got $clip",
    ))
    c0 = zeros(Int, nb)
    c1 = zeros(Int, nb)
    w = 2c / nb
    @inbounds for j in 1:n
        l = Float64(llr[j])
        i = clamp(floor(Int, (l + c) / w) + 1, 1, nb)
        if bits[j] == 0
            c0[i] += 1
        else
            c1[i] += 1
        end
    end
    n0 = sum(c0)
    n1 = sum(c1)
    (n0 > 0 && n1 > 0) || throw(ArgumentError(
        "both bit values must be present for the histogram estimator",
    ))
    acc = 0.0
    @inbounds for i in 1:nb
        p0 = c0[i] / n0
        p1 = c1[i] / n1
        den = p0 + p1
        den == 0 && continue
        p0 > 0 && (acc += p0 * log2(2p0 / den))
        p1 > 0 && (acc += p1 * log2(2p1 / den))
    end
    return acc / 2
end

# ten Brink, Speidel, Yan 2004, eq. (9)–(11) constants.
const _J_H1 = 0.3073
const _J_H2 = 0.8935
const _J_H3 = 1.1064
const _JINV_A1 = 1.09542
const _JINV_B1 = 0.214217
const _JINV_C1 = 2.33737
const _JINV_A2 = 0.706692
const _JINV_B2 = 0.386013
const _JINV_C2 = -1.75017
const _JINV_SPLIT = 0.3646

@doc raw"""
    exit_j(sigma::Real) -> Float64

Mutual information of a consistent symmetric-Gaussian LLR, ``L \mid b \sim \mathcal{N}\bigl((1-2b)\sigma^2/2,\ \sigma^2\bigr)``:

```math
J(\sigma) = 1 - \mathbb{E}\!\left[\log_2\!\bigl(1 + e^{-L}\bigr)\,\middle|\, b = 0\right],
\qquad \sigma \ge 0.
```

Evaluated with the polynomial approximation of ten Brink, Speidel, and Yan,
IEEE Trans. Commun., 52(4):670–678, 2004, eq. (9):

```math
J(\sigma) \approx \bigl(1 - 2^{-H_1 \sigma^{2H_2}}\bigr)^{H_3},
\quad H_1 = 0.3073,\ H_2 = 0.8935,\ H_3 = 1.1064.
```

This is the Gaussian-symmetric shortcut, **not** a measurement: real demapper
or decoder extrinsic LLRs are only approximately symmetric Gaussian, so EXIT
work should measure with [`exit_mi_soft`](@ref) or [`exit_mi_hist`](@ref) and
use `exit_j` / [`exit_jinv`](@ref) for synthesis and cross-checks. The
approximation error against the exact integral is on the order of ``10^{-2}``
bit in the worst case; the BIAWGN validation in `test/exit.jl` pins it at the
tested SNRs. `exit_j(0) == 0` and `exit_j(Inf) == 1`.
"""
function exit_j(sigma::Real)
    σ = float(sigma)
    σ >= 0 || throw(ArgumentError("sigma must be ≥ 0, got $sigma"))
    σ == 0 && return 0.0
    isinf(σ) && return 1.0
    t = exp((-_J_H1) * σ^(2 * _J_H2) * log(2.0))
    return (1 - t)^_J_H3
end

@doc raw"""
    exit_jinv(I::Real) -> Float64

Inverse of [`exit_j`](@ref): the ``\sigma`` of a consistent symmetric-Gaussian
LLR with mutual information ``I``. Piecewise polynomial approximation of ten
Brink, Speidel, and Yan, IEEE Trans. Commun., 52(4):670–678, 2004, eq. (11):

```math
J^{-1}(I) \approx
\begin{cases}
a_1 I^2 + b_1 I + c_1\sqrt{I}, & 0 \le I \le 0.3646,\[0.4em]
-a_2\ln\!\bigl[b_2(1-I)\bigr] - c_2 I, & 0.3646 < I < 1,
\end{cases}
```

with ``a_1 = 1.09542``, ``b_1 = 0.214217``, ``c_1 = 2.33737``,
``a_2 = 0.706692``, ``b_2 = 0.386013``, ``c_2 = -1.75017``. `I` must lie in
``[0, 1]``. `exit_jinv(0) == 0` and `exit_jinv(1) == Inf`. The same
approximation caveat as [`exit_j`](@ref) applies: this is the
Gaussian-symmetric model, not a measurement.
"""
function exit_jinv(I::Real)
    x = float(I)
    (0 <= x <= 1) || throw(ArgumentError("I must lie in [0, 1], got $I"))
    x == 0 && return 0.0
    x == 1 && return Inf
    if x <= _JINV_SPLIT
        return _JINV_A1 * x^2 + _JINV_B1 * x + _JINV_C1 * sqrt(x)
    end
    return -_JINV_A2 * log(_JINV_B2 * (1 - x)) - _JINV_C2 * x
end

@doc raw"""
    EXITCurve

One measured extrinsic transfer characteristic: mutual-information pairs
``(I_A, I_E)``.

# Fields
- `ia::Vector{Float64}`: a-priori mutual information, the curve input.
- `ie::Vector{Float64}`: extrinsic mutual information measured at that input.

Produced by [`exit_demap_transfer`](@ref) and [`exit_decoder_transfer`](@ref).
Plot with [`plot_exit_chart`](@ref) after `using Plots`.
"""
struct EXITCurve
    ia::Vector{Float64}
    ie::Vector{Float64}
    function EXITCurve(ia::AbstractVector{<:Real}, ie::AbstractVector{<:Real})
        length(ia) == length(ie) || throw(ArgumentError(
            "ia and ie must have the same length, got $(length(ia)) and $(length(ie))",
        ))
        length(ia) >= 1 || throw(ArgumentError("an EXITCurve needs at least one point"))
        all(isfinite, ia) || throw(ArgumentError("ia must be finite"))
        return new(Float64.(ia), Float64.(ie))
    end
end

function Base.show(io::IO, curve::EXITCurve)
    print(io, "EXITCurve(", length(curve.ia), " points, ie ",
          round.(extrema(curve.ie); digits=4), ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", curve::EXITCurve)
    return show(io, curve)
end

function _check_ia_grid(ia)
    grid = Float64.(collect(ia))
    length(grid) >= 1 || throw(ArgumentError("ia grid must not be empty"))
    all(x -> 0 <= x <= 1, grid) || throw(ArgumentError(
        "ia values must lie in [0, 1], got extrema $(extrema(grid))",
    ))
    return grid
end

function _check_pos_int(x::Integer, name::String)
    v = Int(x)
    v >= 1 || throw(ArgumentError("$name must be ≥ 1, got $v"))
    return v
end

function _check_sigma_max(sigma_max::Real)
    s = float(sigma_max)
    (isfinite(s) && s > 0) || throw(ArgumentError(
        "sigma_max must be finite and positive, got $sigma_max",
    ))
    return s
end

function _check_method(method::Symbol)
    (method === :soft || method === :hist) || throw(ArgumentError(
        "method must be :soft or :hist, got $(repr(method))",
    ))
    return nothing
end

_mi(method::Symbol, llr, bits; nbins, clip) =
    method === :soft ? exit_mi_soft(llr, bits) :
                       exit_mi_hist(llr, bits; nbins=nbins, clip=clip)

# Consistent symmetric-Gaussian a-priori LLR of one bit (EXIT synthesis).
_prior_llr(rng::_ExitRng, bit::Int, σA::Float64) =
    0.5 * σA^2 * (1 - 2 * bit) + σA * _next_normal(rng)

# Channel LLR of one coded bit on BPSK over AWGN, symbol energy 1,
# per-component noise variance σ²: L = 2y/σ² with y = (1 - 2b) + σn.
_channel_llr(rng::_ExitRng, bit::Int, σ2::Float64) =
    2 * ((1 - 2 * bit) + sqrt(σ2) * _next_normal(rng)) / σ2

function _reject_differential_exit(scheme)
    if scheme isa DBPSK || scheme isa DQPSK
        throw(ArgumentError(
            "EXIT demapper measurement is coherent; $(typeof(scheme)) is differential",
        ))
    end
    return nothing
end

function _checked_labels_exit(labels, M::Int)
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

# Modulate raw bits with an explicit label permutation, the same convention as
# bicm_map: labels[i] is the MSB-first integer label of constellation point i.
# Duplicated from BICM internals on purpose: EXIT depends only on the public
# BICM surface (bicm_demap), and bicm_map always encodes first.
function _map_labeled(bits::Vector{Int}, scheme::ModulationScheme, labels)
    _reject_differential_exit(scheme)
    k = bits_per_symbol(scheme)
    M = 1 << k
    labs = _checked_labels_exit(labels, M)
    labs === nothing && return modulate(bits, scheme)
    length(bits) % k == 0 || throw(ArgumentError(
        "bit length $(length(bits)) is not a multiple of $k",
    ))
    inv = Vector{Int}(undef, M)
    @inbounds for i in 1:M
        inv[labs[i] + 1] = i
    end
    pts = constellation(scheme)
    nsym = length(bits) ÷ k
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

function _awgn_exit(rng::_ExitRng, symbols, σ::Float64)
    if eltype(symbols) <: Real
        out = Vector{Float64}(undef, length(symbols))
        @inbounds for i in eachindex(symbols)
            out[i] = symbols[i] + σ * _next_normal(rng)
        end
    else
        out = Vector{ComplexF64}(undef, length(symbols))
        @inbounds for i in eachindex(symbols)
            out[i] = symbols[i] + σ * (_next_normal(rng) + im * _next_normal(rng))
        end
    end
    return out
end

@doc raw"""
    exit_demap_transfer(scheme::ModulationScheme; ia=0.0:0.1:1.0, esn0_db::Real=0.0,
                        frames::Integer=8, n_symbols::Integer=2000, labels=nothing,
                        approx::Symbol=:exact, method::Symbol=:soft,
                        nbins::Integer=257, clip::Real=50.0,
                        sigma_max::Real=64.0, seed::Integer=0x0c0de1162) -> EXITCurve

Extrinsic transfer characteristic of the a-priori demapper
[`CommToolkit.BICM.bicm_demap`](@ref) on AWGN: for each ``I_A`` in `ia`,
synthesise consistent Gaussian priors, demap, and measure the mutual
information of the extrinsic output against the true bits.

# Arguments
- `scheme`: [`BPSK`](@ref), [`QPSK`](@ref), [`PSK`](@ref), or square [`QAM`](@ref). Differential schemes are rejected, as in `bicm_demap`.
- `ia`: a-priori mutual-information grid, values in ``[0, 1]``.
- `esn0_db`: ``E_s/N_0`` in dB of the coded symbols, unit symbol energy. The per-component noise variance is ``\sigma^2 = 1/(2\,E_s/N_0)``, the `awgn` convention. This is **not** information-bit ``E_b/N_0``: with a code of rate ``R`` on ``\log_2 M`` bits per symbol, ``E_b/N_0 = E_s/N_0 / (R\log_2 M)``.
- `frames`: independent blocks measured per grid point.
- `n_symbols`: symbols per frame; `frames × n_symbols × log2(M)` LLR samples feed each point.
- `labels`: bit-label permutation as in [`CommToolkit.BICM.bicm_map`](@ref), e.g. [`CommToolkit.BICM.bicm_psk8_id_labels`](@ref); `nothing` is Gray.
- `approx`: forwarded to `bicm_demap`, `:exact` or `:maxlog`.
- `method`: MI estimator, `:soft` ([`exit_mi_soft`](@ref), default) or `:hist` ([`exit_mi_hist`](@ref) with `nbins` and `clip`).
- `sigma_max`: the a-priori ``\sigma_A`` is ``\min(J^{-1}(I_A), \texttt{sigma\_max})``, so `I_A = 1` is realised as a very reliable but finite Gaussian prior rather than ``\pm\infty`` LLRs. Documented Gaussian-symmetric synthesis: the prior *is* consistent symmetric Gaussian by construction — that is the assumption under which decoder/demapper transfer curves are universally measured — while the measured extrinsic output is *not* assumed Gaussian.
- `seed`: seed of the local xorshift64* + Box–Muller stream. Fixed seeds reproduce exactly across Julia versions; the global RNG is not touched.

# Returns
[`EXITCurve`](@ref) with `ia` as requested and `ie` measured.

For BPSK and Gray QPSK the extrinsic LLR does not depend on the prior, so the
curve is flat at the BIAWGN mutual information — the validation in
`test/exit.jl` checks the ``I_A = 0`` point against the BIAWGN capacity
integral. Set-partitioning 8-PSK is the teaching example where
``T(I_A)`` rises with ``I_A``. Reference for demapper transfer
characteristics: ten Brink, AEU Int. J. Electron. Commun., 54(6):389–398,
2000.
"""
function exit_demap_transfer(
    scheme::ModulationScheme;
    ia=0.0:0.1:1.0,
    esn0_db::Real=0.0,
    frames::Integer=8,
    n_symbols::Integer=2000,
    labels=nothing,
    approx::Symbol=:exact,
    method::Symbol=:soft,
    nbins::Integer=257,
    clip::Real=50.0,
    sigma_max::Real=64.0,
    seed::Integer=_DEFAULT_SEED,
)
    _reject_differential_exit(scheme)
    k = bits_per_symbol(scheme)
    M = 1 << k
    labs = _checked_labels_exit(labels, M)
    grid = _check_ia_grid(ia)
    nfr = _check_pos_int(frames, "frames")
    nsym = _check_pos_int(n_symbols, "n_symbols")
    isfinite(esn0_db) || throw(ArgumentError("esn0_db must be finite, got $esn0_db"))
    _check_method(method)
    smax = _check_sigma_max(sigma_max)
    σ2 = inv(2 * snr_db_to_linear(esn0_db))
    σ = sqrt(σ2)
    rng = _ExitRng(seed)
    nbits = nsym * k
    ie = Vector{Float64}(undef, length(grid))
    for (idx, Ia) in enumerate(grid)
        σA = min(exit_jinv(Ia), smax)
        le = sizehint!(Float64[], nfr * nbits)
        bt = sizehint!(Int[], nfr * nbits)
        for _ in 1:nfr
            bits = _bits(rng, nbits)
            symbols = _map_labeled(bits, scheme, labs)
            rx = _awgn_exit(rng, symbols, σ)
            La = [_prior_llr(rng, bits[j], σA) for j in 1:nbits]
            ext = bicm_demap(rx, scheme, σ2; apriori=La, approx=approx, labels=labs)
            append!(le, ext)
            append!(bt, bits)
        end
        ie[idx] = _mi(method, le, bt; nbins=nbins, clip=clip)
    end
    return EXITCurve(grid, ie)
end

@doc raw"""
    exit_decoder_transfer(; decoder::Symbol=:bcjr, ia=0.0:0.1:1.0,
                          ebn0_db::Real=1.0, frames::Integer=8,
                          n_info::Union{Nothing,Integer}=nothing, code=nothing,
                          algorithm::Symbol=:normalized_min_sum,
                          ldpc_iterations::Integer=20, method::Symbol=:soft,
                          nbins::Integer=257, clip::Real=50.0,
                          sigma_max::Real=64.0, seed::Integer=0x0c0de1162) -> EXITCurve

Information-bit extrinsic transfer characteristic of an existing decoder on
coded BPSK over AWGN. The decoders are called, not rewritten.

# Arguments
- `decoder`:
  - `:bcjr` — the terminated ``(7, 5)`` rate-1/2 code: [`conv_encode`](@ref CommToolkit.Coding.conv_encode) plus [`bcjr_decode`](@ref CommToolkit.Coding.bcjr_decode). `bcjr_decode` already returns extrinsic information-bit LLRs, ``L_e = L_{\mathrm{post}} - L_a``. The actual rate ``R = n/(2n+4)`` including the two tail bits sets the coded-bit noise.
  - `:ldpc` — a short [`LDPCCode`](@ref CommToolkit.Coding.LDPCCode), required as `code`. The decoder input is channel LLR plus synthetic prior (the decoder has no separate prior entry; see [`ldpc_decode_details`](@ref CommToolkit.Coding.ldpc_decode_details)). The measured extrinsic on information positions is ``L_e = L_{\mathrm{post}} - L_a`` with `early_stop = false`, `max_iter = ldpc_iterations`, and `algorithm` forwarded.

  `:sova` is not accepted. `sova_decode(...; return_coded=true)` (1.2.2) can return a teaching coded-bit extrinsic, but that vector is a competing-path gap rather than the MAP split this curve measures, and it is not wired in.
- `ia`, `method`, `nbins`, `clip`, `sigma_max`, `seed`: as in [`exit_demap_transfer`](@ref).
- `ebn0_db`: information-bit ``E_b/N_0`` in dB. Coded bits are BPSK of energy 1, so the per-component noise variance is ``\sigma^2 = (1/R)/(2\,E_b/N_0)``, the convention documented on `ldpc_decode_details`.
- `frames`: frames measured per grid point.
- `n_info`: information bits per frame. For `:bcjr` the default is 384. For `:ldpc` it must equal `code.k` and defaults to it.

# Returns
[`EXITCurve`](@ref) with `ia` as requested and `ie` measured on information bits.

Anchors checked in `test/exit.jl`: at high ``E_b/N_0`` the ``(7, 5)`` curve
starts near ``I_E = 1`` even at ``I_A = 0``; the curve rises with ``I_A``; and
``I_E \in [0, 1]`` at every point. Reference for decoder transfer
characteristics: ten Brink, IEEE Trans. Commun., 49(10):1727–1737, 2001.
"""
function exit_decoder_transfer(;
                               decoder::Symbol=:bcjr,
                               ia=0.0:0.1:1.0,
                               ebn0_db::Real=1.0,
                               frames::Integer=8,
                               n_info::Union{Nothing,Integer}=nothing,
                               code=nothing,
                               algorithm::Symbol=:normalized_min_sum,
                               ldpc_iterations::Integer=20,
                               method::Symbol=:soft,
                               nbins::Integer=257,
                               clip::Real=50.0,
                               sigma_max::Real=64.0,
                               seed::Integer=_DEFAULT_SEED)
    grid = _check_ia_grid(ia)
    nfr = _check_pos_int(frames, "frames")
    isfinite(ebn0_db) || throw(ArgumentError("ebn0_db must be finite, got $ebn0_db"))
    _check_method(method)
    smax = _check_sigma_max(sigma_max)
    γb = snr_db_to_linear(ebn0_db)
    if decoder === :bcjr
        code === nothing || throw(ArgumentError("decoder=:bcjr takes no LDPCCode"))
        ninf = n_info === nothing ? 384 : _check_pos_int(n_info, "n_info")
        R = ninf / (2 * (ninf + 2))
    elseif decoder === :ldpc
        code isa LDPCCode || throw(ArgumentError("decoder=:ldpc requires an LDPCCode as code"))
        ninf = code.k
        n_info === nothing || Int(n_info) == ninf || throw(ArgumentError(
            "n_info must be $ninf for this LDPC code, got $n_info",
        ))
        R = code.k / code.n
    else
        throw(ArgumentError("decoder must be :bcjr or :ldpc, got $(repr(decoder))"))
    end
    ldpc_it = Int(ldpc_iterations)
    ldpc_it >= 0 || throw(ArgumentError("ldpc_iterations must be ≥ 0, got $ldpc_it"))
    σ2 = inv(R) / (2 * γb)
    rng = _ExitRng(seed)
    ie = Vector{Float64}(undef, length(grid))
    for (idx, Ia) in enumerate(grid)
        σA = min(exit_jinv(Ia), smax)
        le = sizehint!(Float64[], nfr * ninf)
        bt = sizehint!(Int[], nfr * ninf)
        for _ in 1:nfr
            info = _bits(rng, ninf)
            coded = decoder === :bcjr ? conv_encode(info) : ldpc_encode(code, info)
            Lch = [_channel_llr(rng, coded[j], σ2) for j in eachindex(coded)]
            if decoder === :bcjr
                La = [_prior_llr(rng, info[j], σA) for j in 1:ninf]
                ext = bcjr_decode(Lch; apriori=La)
                append!(le, ext)
            else
                La = [_prior_llr(rng, coded[j], σA) for j in eachindex(coded)]
                detail = ldpc_decode_details(code, Lch .+ La;
                                             algorithm=algorithm,
                                             max_iter=ldpc_it,
                                             early_stop=false)
                ip = code.info_pos
                append!(le, detail.posterior[ip] .- La[ip])
            end
            append!(bt, info)
        end
        ie[idx] = _mi(method, le, bt; nbins=nbins, clip=clip)
    end
    return EXITCurve(grid, ie)
end

function _plots_extension_exit()
    ext = Base.get_extension(parentmodule(@__MODULE__), :CommToolkitPlotsExt)
    if ext === nothing
        error("plotting requires Plots.jl. Load it with `using Plots` so the CommToolkitPlotsExt extension is activated. `using CommToolkit` alone does not load Plots.")
    end
    return ext
end

@doc raw"""
    plot_exit_chart(demapper::EXITCurve, decoder::EXITCurve; title::String="EXIT chart")
    plot_exit_chart(curve::EXITCurve; title::String="EXIT chart")

The classical EXIT chart: the demapper curve drawn as ``I_E`` versus ``I_A``
and the decoder curve with swapped axes, ``I_A`` versus ``I_E``, so a decoding
trajectory would zig-zag between the two lines. The single-curve method draws
one transfer characteristic as measured. This function only plots measured
[`EXITCurve`](@ref) data; it does not simulate a trajectory.

Requires Plots.jl. Loading Plots activates the package extension
`CommToolkitPlotsExt`; `using CommToolkit` by itself does not import Plots.
Returns the `Plots.Plot`. Reference: ten Brink, IEEE Trans. Commun.,
49(10):1727–1737, 2001.
"""
function plot_exit_chart(args...; kwargs...)
    return _plots_extension_exit().plot_exit_chart(args...; kwargs...)
end

end # module EXIT
