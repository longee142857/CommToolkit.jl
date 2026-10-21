module Modulation

export ModulationScheme, BPSK, QPSK, QAM, modulate, demodulate, demodulate_soft,
    constellation, bits_per_symbol

@doc raw"""
    abstract type ModulationScheme end

Marker for a digital constellation. Concrete schemes are [`BPSK`](@ref), [`QPSK`](@ref), and square [`QAM`](@ref).

Symbols are normalised to unit average energy ``E_s = 1``. With ``k`` bits per symbol, ``E_b = 1/k``.
"""
abstract type ModulationScheme end

@doc raw"""
    BPSK()

Binary phase-shift keying. The bit map is ``0 \mapsto +1``, ``1 \mapsto -1``, so ``E_s = E_b = 1``.

In AWGN the bit error probability is ``Q(\sqrt{2 E_b/N_0})``
(Proakis and Salehi, *Digital Communications*, 5th ed., McGraw-Hill, 2008, §5.2).
"""
struct BPSK <: ModulationScheme end

# Gray QPSK, unit average energy. Index 2*b0 + b1 + 1:
# 00 → (+1+j)/√2, 01 → (+1−j)/√2, 10 → (−1+j)/√2, 11 → (−1−j)/√2.
const _QPSK = ComplexF64[1 + 1im, 1 - 1im, -1 + 1im, -1 - 1im] ./ sqrt(2)

@doc raw"""
    QPSK()

Gray-labelled quadrature phase-shift keying with unit average symbol energy.

Bit pairs ``(b_0, b_1)``, ``b_0`` first, map to

| bits | symbol |
|------|--------|
| 00 | ``(+1+j)/\sqrt{2}`` |
| 01 | ``(+1-j)/\sqrt{2}`` |
| 10 | ``(-1+j)/\sqrt{2}`` |
| 11 | ``(-1-j)/\sqrt{2}`` |

Adjacent quadrants differ by one bit. Each component is an independent BPSK decision of amplitude ``1/\sqrt{2}``, so the bit error probability in AWGN equals the BPSK formula ``Q(\sqrt{2 E_b/N_0})`` with ``E_b = E_s/2`` (Proakis and Salehi, 5th ed., §5.2).
"""
struct QPSK <: ModulationScheme end

@doc raw"""
    QAM(m::Integer)

Square Gray-labelled ``M``-QAM, ``M \in \{4, 16, 64, 256, \ldots\}``.

Each axis is a Gray-labelled ``\sqrt{M}``-PAM alphabet ``\{-(\sqrt{M}-1), \ldots, \sqrt{M}-1\}`` scaled by ``\sqrt{3 / (2(M-1))}``, which gives ``E_s = 1``. The two most-significant bits of each ``\log_2 M``-bit block select the in-phase Gray label; the remaining bits select the quadrature label (MSB first).

The exact average bit error rate of this labelling in AWGN is the square-QAM formula of Cho and Yoon, IEEE Trans. Commun., 50(7):1074–1080, 2002. A common high-SNR approximation is

```math
P_b \approx \frac{4}{\log_2 M}\left(1 - \frac{1}{\sqrt{M}}\right)
Q\left(\sqrt{\frac{3 \log_2 M}{M-1}\,\frac{E_b}{N_0}}\right)
```

(Proakis and Salehi, 5th ed., §5.2). ``M = 4`` is Gray 4-QAM, distinct in type from [`QPSK`](@ref) but with the same geometry.
"""
struct QAM <: ModulationScheme
    m::Int
    function QAM(m::Integer)
        M = Int(m)
        if M < 4 || !ispow2(M) || isqrt(M)^2 != M
            throw(ArgumentError("M must be a square power of two ≥ 4 (4, 16, 64, …), got $m"))
        end
        return new(M)
    end
end

@doc raw"""
    bits_per_symbol(scheme::ModulationScheme) -> Int

Number of bits carried by one symbol: 1 for BPSK, 2 for QPSK, and ``\log_2 M`` for `QAM(M)`.
"""
bits_per_symbol(::BPSK) = 1
bits_per_symbol(::QPSK) = 2
bits_per_symbol(scheme::QAM) = Int(log2(scheme.m))

function _qam_constellation(M::Int)
    n = isqrt(M)
    scale = sqrt(2 * (M - 1) / 3)
    pam = [-(n - 1):2:(n - 1);] ./ scale
    gray = [xor(i, i >> 1) for i in 0:(n - 1)]
    constell = Vector{ComplexF64}(undef, M)
    for i in 1:n, j in 1:n
        constell[gray[i] * n + gray[j] + 1] = pam[i] + pam[j] * im
    end
    return constell
end

@doc raw"""
    constellation(scheme::ModulationScheme) -> Vector{ComplexF64}

Ideal constellation points in the modulator's index order, with average energy 1.

For `QAM`, index ``i`` (1-based) is the symbol for the ``\log_2 M`` bits of ``i - 1`` written MSB first.
"""
constellation(::BPSK) = ComplexF64[1, -1]
constellation(::QPSK) = copy(_QPSK)
constellation(scheme::QAM) = _qam_constellation(scheme.m)

function _check_bits(bits)
    all(b -> b == 0 || b == 1, bits) || throw(ArgumentError("bits must be 0 or 1"))
    return nothing
end

@doc raw"""
    modulate(bits::AbstractVector{<:Integer}, scheme::ModulationScheme)

Map bits to unit-energy symbols.

# Arguments
- `bits`: entries in ``\{0, 1\}``. Length must be a multiple of [`bits_per_symbol`](@ref)`(scheme)`.
- `scheme`: [`BPSK`](@ref), [`QPSK`](@ref), or [`QAM`](@ref).

# Returns
- `BPSK`: real vector, ``0 \mapsto +1``, ``1 \mapsto -1``.
- `QPSK` and `QAM`: `Vector{ComplexF64}` using the Gray maps documented on those types.

# References
Proakis and Salehi, *Digital Communications*, 5th ed., §5.2 (constellation normalisation and Gray labelling).
"""
function modulate end

function modulate(bits::AbstractVector{<:Integer}, ::BPSK)
    _check_bits(bits)
    return 1 .- 2 .* bits
end

function modulate(bits::AbstractVector{<:Integer}, ::QPSK)
    _check_bits(bits)
    length(bits) % 2 == 0 || throw(ArgumentError("bit length must be even for QPSK"))
    n_sym = length(bits) ÷ 2
    syms = Vector{ComplexF64}(undef, n_sym)
    @inbounds for i in 1:n_sym
        b0 = bits[2i - 1]
        b1 = bits[2i]
        syms[i] = _QPSK[2 * b0 + b1 + 1]
    end
    return syms
end

function modulate(bits::AbstractVector{<:Integer}, scheme::QAM)
    _check_bits(bits)
    M = scheme.m
    k = bits_per_symbol(scheme)
    length(bits) % k == 0 ||
        throw(ArgumentError("bit length must be a multiple of log2(M) = $k"))
    constell = _qam_constellation(M)
    n_sym = length(bits) ÷ k
    syms = Vector{ComplexF64}(undef, n_sym)
    @inbounds for i in 1:n_sym
        idx = 0
        base = (i - 1) * k
        for j in 1:k
            idx = (idx << 1) | bits[base + j]
        end
        syms[i] = constell[idx + 1]
    end
    return syms
end

@doc raw"""
    demodulate(symbols, scheme::ModulationScheme) -> Vector{Int}

Hard decisions, inverse to [`modulate`](@ref) in the noiseless case.

# Arguments
- `symbols`: real or complex samples for BPSK (only the real part is tested). Complex samples for QPSK and QAM.
- `scheme`: the modulator that produced the symbols.

# Returns
Bits in ``\{0, 1\}``, MSB first for QPSK and QAM. QPSK decides each component by its sign. QAM uses minimum-distance decoding, which is maximum likelihood in isotropic AWGN (Proakis and Salehi, 5th ed., §4.3 and §5.2).
"""
function demodulate end

function demodulate(syms::AbstractVector{<:Real}, ::BPSK)
    return Int.(syms .< 0)
end

function demodulate(syms::AbstractVector{<:Complex}, ::BPSK)
    return Int.(real.(syms) .< 0)
end

function demodulate(syms::AbstractVector{<:Complex}, ::QPSK)
    n_sym = length(syms)
    bits = Vector{Int}(undef, 2 * n_sym)
    @inbounds for i in 1:n_sym
        bits[2i - 1] = Int(real(syms[i]) < 0)
        bits[2i] = Int(imag(syms[i]) < 0)
    end
    return bits
end

function demodulate(syms::AbstractVector{<:Complex}, scheme::QAM)
    M = scheme.m
    k = bits_per_symbol(scheme)
    constell = _qam_constellation(M)
    bits = Vector{Int}(undef, length(syms) * k)
    @inbounds for i in eachindex(syms)
        _, min_idx = findmin(abs2.(syms[i] .- constell))
        label = min_idx - 1
        base = (i - 1) * k
        for j in 1:k
            bits[base + j] = (label >> (k - j)) & 1
        end
    end
    return bits
end

@doc raw"""
    demodulate_soft(symbols, scheme::ModulationScheme, noise_var::Real; approx::Symbol=:exact) -> Vector{Float64}

Bit log-likelihood ratios of Gray-labelled [`BPSK`](@ref), [`QPSK`](@ref), and square [`QAM`](@ref) observed in AWGN.

# Arguments
- `symbols`: real or complex samples for BPSK. Complex samples for QPSK and QAM. BPSK uses only the real part.
- `scheme`: the modulator whose labelling [`demodulate`](@ref) inverts.
- `noise_var`: variance ``\sigma^2 = N_0/2`` of one real noise component. Finite and positive. This is the variance added by `awgn`:

```math
\sigma^2 = \frac{E_b}{2\,E_b/N_0}, \qquad N_0 = 2\sigma^2.
```

  The positional call `awgn(symbols, snr_db, bits_per_symbol(scheme))` has unit symbol energy, so ``E_b = 1/\texttt{bits\_per\_symbol}`` and ``\sigma^2 = 1/(2\,\texttt{bits\_per\_symbol}\,10^{\texttt{snr\_db}/10})``. For a complex symbol the total noise power is ``N_0``.
- `approx`: `:exact` enumerates the constellation. `:maxlog` keeps the nearest point of each bit subset.

# Returns
Length ``\texttt{length(symbols)}\times k`` vector, ``k = `` [`bits_per_symbol`](@ref)`(scheme)`, in the bit order of [`demodulate`](@ref). The entry for bit ``b`` is

```math
L(b) = \log\frac{p(y\mid b=0)}{p(y\mid b=1)}.
```

A positive value favours bit 0. Symbol priors are equal, so these are channel LLRs (the bit prior ratio is 1).

# Formula
Write ``\mathcal{S}_i^{(0)}`` and ``\mathcal{S}_i^{(1)}`` for the points whose label has 0 or 1 in bit position ``i`` (MSB first). With total noise power ``N_0``,

```math
p(y\mid s)\propto\exp\bigl(-|y-s|^2/N_0\bigr).
```

Exact marginalisation:

```math
L_i(y)=\log\sum_{s\in\mathcal{S}_i^{(0)}}\exp\bigl(-|y-s|^2/N_0\bigr)
-\log\sum_{s\in\mathcal{S}_i^{(1)}}\exp\bigl(-|y-s|^2/N_0\bigr).
```

Max-log approximation (each sum replaced by its largest term):

```math
L_i(y)\approx\frac{1}{N_0}\Bigl(\min_{s\in\mathcal{S}_i^{(1)}}|y-s|^2-\min_{s\in\mathcal{S}_i^{(0)}}|y-s|^2\Bigr).
```

For BPSK, ``0\mapsto +1``, ``1\mapsto -1``, both expressions equal ``2\Re(y)/\sigma^2``. For Gray QPSK the two bits reduce to that expression on the real and imaginary parts with amplitude ``1/\sqrt{2}``, and exact and max-log agree. On square QAM the in-phase bits depend only on ``\Re(y)`` and the quadrature bits only on ``\Im(y)``.

# References
Proakis and Salehi, *Digital Communications*, 5th ed., §4.2 (binary AWGN likelihood ratio) and §5.2 (constellations and Gray labelling). The max-log replacement of the log-sum is the approximation in Robertson, Villebrun, and Hoeher, Proc. IEEE ICC, 1995, pp. 1009–1013.
"""
function demodulate_soft(syms::AbstractVector{<:Real}, scheme::ModulationScheme, noise_var::Real;
                          approx::Symbol=:exact)
    scheme isa BPSK || throw(ArgumentError(
        "soft demodulation of $(typeof(scheme)) requires complex symbols",
    ))
    return demodulate_soft(ComplexF64.(syms), scheme, noise_var; approx=approx)
end

function demodulate_soft(syms::AbstractVector{<:Complex}, scheme::ModulationScheme, noise_var::Real;
                          approx::Symbol=:exact)
    _check_soft_args(noise_var, approx)
    pts = constellation(scheme)
    k = bits_per_symbol(scheme)
    M = length(pts)
    idx0 = [Int[] for _ in 1:k]
    idx1 = [Int[] for _ in 1:k]
    for label in 0:(M - 1)
        for j in 1:k
            bucket = ((label >> (k - j)) & 1) == 0 ? idx0[j] : idx1[j]
            push!(bucket, label + 1)
        end
    end
    invN0 = inv(2 * float(noise_var))
    out = Vector{Float64}(undef, length(syms) * k)
    scratch = Vector{Float64}(undef, M)
    @inbounds for (i, y) in enumerate(syms)
        yc = ComplexF64(y)
        base = (i - 1) * k
        if approx === :maxlog
            for m in 1:M
                scratch[m] = abs2(yc - pts[m])
            end
            for j in 1:k
                d0 = Inf
                d1 = Inf
                for m in idx0[j]
                    scratch[m] < d0 && (d0 = scratch[m])
                end
                for m in idx1[j]
                    scratch[m] < d1 && (d1 = scratch[m])
                end
                out[base + j] = (d1 - d0) * invN0
            end
        else
            for m in 1:M
                scratch[m] = -abs2(yc - pts[m]) * invN0
            end
            for j in 1:k
                out[base + j] = _logsumexp_at(scratch, idx0[j]) - _logsumexp_at(scratch, idx1[j])
            end
        end
    end
    return out
end

function _check_soft_args(noise_var::Real, approx::Symbol)
    (isfinite(noise_var) && noise_var > 0) || throw(ArgumentError(
        "noise_var must be finite and positive, got $noise_var",
    ))
    (approx === :exact || approx === :maxlog) || throw(ArgumentError(
        "approx must be :exact or :maxlog, got $(repr(approx))",
    ))
    return nothing
end

function _logsumexp_at(metrics::Vector{Float64}, idx::Vector{Int})
    m = -Inf
    @inbounds for i in idx
        metrics[i] > m && (m = metrics[i])
    end
    m == -Inf && return -Inf
    m == Inf && return Inf
    s = 0.0
    @inbounds for i in idx
        s += exp(metrics[i] - m)
    end
    return m + log(s)
end

end # module Modulation
