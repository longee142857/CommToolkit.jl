module Modulation

export ModulationScheme, BPSK, QPSK, PSK, DBPSK, DQPSK, QAM, modulate, demodulate,
    demodulate_soft, constellation, bits_per_symbol

@doc raw"""
    abstract type ModulationScheme end

Marker for a digital constellation or a differential phase scheme. Concrete schemes are [`BPSK`](@ref), [`QPSK`](@ref), Gray [`PSK`](@ref), differential [`DBPSK`](@ref) and [`DQPSK`](@ref), and square [`QAM`](@ref).

Symbols are normalised to unit average energy ``E_s = 1``. With ``k`` bits per symbol, ``E_b = 1/k``. Differential modulators prepend one reference symbol that carries no information bits; ``k`` still counts only the information bits on each data symbol.
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
    PSK(m::Integer)

Coherent Gray ``M``-PSK, ``M`` a power of two and ``M \ge 4``, with unit average symbol energy.

Natural index ``m \in \{0,\ldots,M-1\}`` sits at

```math
s = \exp(j 2\pi m / M),
```

so index 0 is on the positive real axis. The modulator stores that point under the binary-reflected Gray label ``m \oplus (m \gg 1)``. Bit blocks are MSB first, the same index order as [`QAM`](@ref). Adjacent phases differ by one bit. ``M = 4`` is Gray 4-PSK on the axes. It has the same geometry as [`QPSK`](@ref) up to a rotation of ``\pi/4``, and it is not the same bit map.

Hard decisions are minimum-distance, which is maximum likelihood in isotropic AWGN and, for equal-energy PSK, nearest-phase. Soft bit LLRs use the coherent AWGN metric in [`demodulate_soft`](@ref).

A nearest-neighbour approximation of the Gray bit error rate in AWGN is

```math
P_b \approx \frac{1}{\log_2 M}
\operatorname{erfc}\!\left(\sqrt{\log_2 M \cdot E_b/N_0}\,\sin\frac{\pi}{M}\right).
```

For ``M = 4`` this is exactly ``Q(\sqrt{2 E_b/N_0})``. For ``M > 4`` it is not exact. The exact average is the mean Gray Hamming distance over the phase-error sectors (Proakis and Salehi, 5th ed., §5.2 and the phase distribution in §4.2). That integral is the reference used by the tests, not a function of this package.
"""
struct PSK <: ModulationScheme
    m::Int
    function PSK(m::Integer)
        M = Int(m)
        if M < 4 || !ispow2(M)
            throw(ArgumentError("M must be a power of two ≥ 4 (4, 8, 16, …), got $m"))
        end
        return new(M)
    end
end

# Phase increments for Gray DQPSK, index 2*b0 + b1 + 1.
# 00 → 0, 01 → +π/2, 10 → −π/2, 11 → π.
const _DQPSK_DELTA = ComplexF64[1, im, -im, -1]

@doc raw"""
    DBPSK()

Differential binary PSK. [`modulate`](@ref) prepends the reference symbol ``1`` and then encodes each bit as a phase increment of ``0`` (bit 0) or ``\pi`` (bit 1):

```math
s_0 = 1, \qquad s_i = s_{i-1}\,(1 - 2 b_i), \qquad i = 1,\ldots,N.
```

Every symbol, including the reference, has energy 1. The block has length ``N + 1``. [`demodulate`](@ref) decides

```math
\hat b_i = 0 \quad\text{when}\quad \Re(y_i y_{i-1}^*) \ge 0,
```

and 1 otherwise. A common phase rotation of the whole block cancels. This is two-symbol differential detection. It is not coherent BPSK, and it is not a noncoherent orthogonal waveform.

In circular AWGN the bit error probability of each detected bit is

```math
P_b = \tfrac12 \exp(-E_b/N_0),
```

with ``E_b = E_s = 1`` (Proakis and Salehi, *Digital Communications*, 5th ed., §4.5). The formula is the error rate of a pair. It does not amortise the reference symbol as an information-bit energy penalty; that symbol is a preamble, so a block of ``N`` information bits uses ``N + 1`` channel uses.

Symbol-by-symbol independent fading destroys the phase increment. A phase that is constant across the block does not.
"""
struct DBPSK <: ModulationScheme end

@doc raw"""
    DQPSK()

Gray differential QPSK. [`modulate`](@ref) prepends the reference ``1`` and maps each dibit, MSB first, to a phase increment

| bits | increment |
|------|-----------|
| 00 | ``0`` |
| 01 | ``+\pi/2`` |
| 10 | ``-\pi/2`` |
| 11 | ``\pi`` |

Adjacent increments differ by one bit. Data symbols lie in ``\{1, j, -1, -j\}`` when the reference is ``1``, each with energy 1. ``E_b = 1/2``. The block length is ``N_{\mathrm{dibit}} + 1``. This is not ``\pi/4``-DQPSK: increments are multiples of ``\pi/2``, not of ``\pi/4``.

[`demodulate`](@ref) takes the quadrant of ``y_i y_{i-1}^*``. A common phase rotation of the whole block cancels. Soft metrics, when requested, are the two-symbol noncoherent likelihood in [`demodulate_soft`](@ref).

With Gray labelling the AWGN bit error probability is

```math
P_b = Q_1(a, b) - \tfrac12 I_0(ab)\,\exp\!\left(-(a^2+b^2)/2\right),
```

```math
a = \sqrt{2\gamma_b\bigl(1 - \sqrt{1/2}\bigr)}, \qquad
b = \sqrt{2\gamma_b\bigl(1 + \sqrt{1/2}\bigr)}, \qquad
\gamma_b = E_b/N_0,
```

where ``Q_1`` is the first-order Marcum Q-function (Proakis and Salehi, 5th ed., §4.5). The same preamble remark as on [`DBPSK`](@ref) applies.
"""
struct DQPSK <: ModulationScheme end

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

Number of information bits carried by one data symbol: 1 for BPSK and DBPSK, 2 for QPSK and DQPSK, and ``\log_2 M`` for `PSK(M)` and `QAM(M)`.

Differential schemes prepend a reference symbol. That symbol is not included in this count.
"""
bits_per_symbol(::BPSK) = 1
bits_per_symbol(::QPSK) = 2
bits_per_symbol(scheme::PSK) = Int(log2(scheme.m))
bits_per_symbol(::DBPSK) = 1
bits_per_symbol(::DQPSK) = 2
bits_per_symbol(scheme::QAM) = Int(log2(scheme.m))

function _psk_constellation(M::Int)
    pts = Vector{ComplexF64}(undef, M)
    @inbounds for m in 0:(M - 1)
        pts[xor(m, m >> 1) + 1] = cis(2π * m / M)
    end
    return pts
end

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

For `QAM` and `PSK`, index ``i`` (1-based) is the symbol for the ``\log_2 M`` bits of ``i - 1`` written MSB first.

For [`DBPSK`](@ref) and [`DQPSK`](@ref) the vector is the set of phase increments, not the absolute points of a memoryless constellation. Absolute DBPSK/DQPSK symbols depend on the preceding symbol.
"""
constellation(::BPSK) = ComplexF64[1, -1]
constellation(::QPSK) = copy(_QPSK)
constellation(scheme::PSK) = _psk_constellation(scheme.m)
constellation(::DBPSK) = ComplexF64[1, -1]
constellation(::DQPSK) = copy(_DQPSK_DELTA)
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
- `scheme`: [`BPSK`](@ref), [`QPSK`](@ref), [`PSK`](@ref), [`DBPSK`](@ref), [`DQPSK`](@ref), or [`QAM`](@ref).

# Returns
- `BPSK`: real vector, ``0 \mapsto +1``, ``1 \mapsto -1``.
- `QPSK`, `PSK`, and `QAM`: `Vector{ComplexF64}` using the Gray maps documented on those types. Length `length(bits) / k`.
- `DBPSK` and `DQPSK`: `Vector{ComplexF64}` of length `length(bits) / k + 1`. Entry 1 is the reference ``1 + 0j``. Each later entry is the previous symbol times the phase increment of that bit or dibit.

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

function modulate(bits::AbstractVector{<:Integer}, scheme::PSK)
    _check_bits(bits)
    k = bits_per_symbol(scheme)
    length(bits) % k == 0 ||
        throw(ArgumentError("bit length must be a multiple of log2(M) = $k"))
    constell = _psk_constellation(scheme.m)
    return _symbols_from_labels(bits, constell, k)
end

function modulate(bits::AbstractVector{<:Integer}, ::DBPSK)
    _check_bits(bits)
    n = length(bits)
    syms = Vector{ComplexF64}(undef, n + 1)
    s = 1.0 + 0.0im
    syms[1] = s
    @inbounds for i in 1:n
        s *= bits[i] == 0 ? 1.0 : -1.0
        syms[i + 1] = s
    end
    return syms
end

function modulate(bits::AbstractVector{<:Integer}, ::DQPSK)
    _check_bits(bits)
    length(bits) % 2 == 0 || throw(ArgumentError("bit length must be even for DQPSK"))
    n = length(bits) ÷ 2
    syms = Vector{ComplexF64}(undef, n + 1)
    s = 1.0 + 0.0im
    syms[1] = s
    @inbounds for i in 1:n
        s *= _DQPSK_DELTA[2 * bits[2i - 1] + bits[2i] + 1]
        syms[i + 1] = s
    end
    return syms
end

function _symbols_from_labels(bits, constell, k::Int)
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
- `symbols`: real or complex samples for BPSK (only the real part is tested). Complex samples for QPSK, PSK, and QAM. For DBPSK and DQPSK, the block produced by [`modulate`](@ref), including the leading reference; a real vector is treated as a purely real complex vector.
- `scheme`: the modulator that produced the symbols.

# Returns
Bits in ``\{0, 1\}``, MSB first for QPSK, PSK, DQPSK, and QAM.

QPSK decides each component by its sign. PSK and QAM use minimum-distance decoding, which is maximum likelihood in isotropic AWGN (Proakis and Salehi, 5th ed., §4.3 and §5.2). DBPSK and DQPSK use the differential rules on those types: the reference contributes no output bits, and a common phase factor on the block does not change the decisions. An empty differential block (no reference) is rejected.
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

function demodulate(syms::AbstractVector{<:Complex}, scheme::PSK)
    return _demodulate_nearest(syms, _psk_constellation(scheme.m), bits_per_symbol(scheme))
end

function demodulate(syms::AbstractVector{<:Real}, scheme::Union{PSK,DBPSK,DQPSK})
    return demodulate(ComplexF64.(syms), scheme)
end

function demodulate(syms::AbstractVector{<:Complex}, ::DBPSK)
    n = length(syms)
    n >= 1 || throw(ArgumentError("DBPSK demodulation needs the leading reference symbol"))
    bits = Vector{Int}(undef, n - 1)
    @inbounds for i in 1:(n - 1)
        bits[i] = Int(real(syms[i + 1] * conj(syms[i])) < 0)
    end
    return bits
end

function demodulate(syms::AbstractVector{<:Complex}, ::DQPSK)
    n = length(syms)
    n >= 1 || throw(ArgumentError("DQPSK demodulation needs the leading reference symbol"))
    bits = Vector{Int}(undef, 2 * (n - 1))
    @inbounds for i in 1:(n - 1)
        θ = angle(syms[i + 1] * conj(syms[i]))
        if -π / 4 <= θ < π / 4
            b0, b1 = 0, 0
        elseif π / 4 <= θ < 3π / 4
            b0, b1 = 0, 1
        elseif -3π / 4 <= θ < -π / 4
            b0, b1 = 1, 0
        else
            b0, b1 = 1, 1
        end
        bits[2i - 1] = b0
        bits[2i] = b1
    end
    return bits
end

function _demodulate_nearest(syms::AbstractVector{<:Complex}, constell, k::Int)
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

function demodulate(syms::AbstractVector{<:Complex}, scheme::QAM)
    return _demodulate_nearest(syms, _qam_constellation(scheme.m), bits_per_symbol(scheme))
end

@doc raw"""
    demodulate_soft(symbols, scheme::ModulationScheme, noise_var::Real; approx::Symbol=:exact) -> Vector{Float64}

Bit log-likelihood ratios of Gray-labelled [`BPSK`](@ref), [`QPSK`](@ref), [`PSK`](@ref), and square [`QAM`](@ref) in coherent AWGN, and two-symbol noncoherent LLRs for [`DBPSK`](@ref) and [`DQPSK`](@ref).

# Arguments
- `symbols`: real or complex samples for BPSK. Complex samples for QPSK, PSK, and QAM. BPSK uses only the real part. A real vector passed with PSK, DBPSK, or DQPSK is treated as purely real. Differential blocks include the leading reference from [`modulate`](@ref).
- `scheme`: the modulator whose labelling [`demodulate`](@ref) inverts.
- `noise_var`: variance ``\sigma^2 = N_0/2`` of one real noise component. Finite and positive. This is the variance added by `awgn`:

```math
\sigma^2 = \frac{E_b}{2\,E_b/N_0}, \qquad N_0 = 2\sigma^2.
```

  The positional call `awgn(symbols, snr_db, bits_per_symbol(scheme))` has unit symbol energy, so ``E_b = 1/\texttt{bits\_per\_symbol}`` and ``\sigma^2 = 1/(2\,\texttt{bits\_per\_symbol}\,10^{\texttt{snr\_db}/10})``. For a complex symbol the total noise power is ``N_0``.
- `approx`: `:exact` enumerates the constellation. `:maxlog` keeps the nearest point of each bit subset.

# Returns
For a coherent scheme, length ``\texttt{length(symbols)}\times k``, ``k = `` [`bits_per_symbol`](@ref)`(scheme)`. For DBPSK or DQPSK, length ``(\texttt{length(symbols)} - 1) \times k``: the reference produces no LLR. Order matches [`demodulate`](@ref). The entry for bit ``b`` is

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

For BPSK, ``0\mapsto +1``, ``1\mapsto -1``, both expressions equal ``2\Re(y)/\sigma^2``. For Gray QPSK the two bits reduce to that expression on the real and imaginary parts with amplitude ``1/\sqrt{2}``, and exact and max-log agree. On square QAM the in-phase bits depend only on ``\Re(y)`` and the quadrature bits only on ``\Im(y)``. On ``M``-PSK with ``M > 4`` the bits do not separate across I and Q, and exact and max-log differ.

# Differential schemes
DBPSK and DQPSK do not use the coherent sum above. With an absolute phase that is constant on the pair ``(y_{i-1}, y_i)`` and is not estimated, a phase increment ``\Delta`` has likelihood

```math
p(\Delta)\propto I_0\bigl(|y_{i-1} + y_i e^{-j\Delta}| / \sigma^2\bigr),
```

``\sigma^2 = `` `noise_var`, ``I_0`` the order-zero modified Bessel function of the first kind. Bit LLRs marginalise the two (DBPSK) or four (DQPSK) increments with equal priors, by the same log-sum-exp (`:exact`) or nearest-term (`:maxlog`) replacement. DBPSK has one increment per bit value, so the two approximations agree. A positive value still favours bit 0.

The metric uses two symbols only. It is not a coherent constellation LLR and it is not a multi-symbol window. Per-symbol independent fading is outside the model.

# References
Proakis and Salehi, *Digital Communications*, 5th ed., §4.2 (binary AWGN likelihood ratio), §4.5 (differential detection), and §5.2 (constellations and Gray labelling). The max-log replacement of the log-sum is the approximation in Robertson, Villebrun, and Hoeher, Proc. IEEE ICC, 1995, pp. 1009–1013.
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

function _logaddexp(a::Float64, b::Float64)
    m = max(a, b)
    isinf(m) && return m
    return m + log(exp(a - m) + exp(b - m))
end

"""log I₀(x) for x ≥ 0. Series below 12, Debye asymptotic above. No SpecialFunctions dependency."""
function _log_besseli0(x::Float64)
    x < 0 && (x = -x)
    x == 0 && return 0.0
    isfinite(x) || return x > 0 ? Inf : NaN
    if x <= 12
        t2 = (0.5 * x)^2
        term = 1.0
        s = 1.0
        for k in 1:48
            term *= t2 / (k * k)
            s += term
            term < eps(s) && break
        end
        return log(s)
    end
    # I₀(x) ∼ e^x / √(2πx) · (1 + 1/(8x) + 9/(128x²) + 225/(3072x³) + 11025/(98304x⁴))
    u = inv(x)
    series = 1 + u * (1 / 8 + u * (9 / 128 + u * (225 / 3072 + u * (11025 / 98304))))
    return x - 0.5 * log(2π * x) + log(series)
end

function _diff_loglik(y0::ComplexF64, y1::ComplexF64, delta::ComplexF64, σ2::Float64)
    return _log_besseli0(abs(y0 + y1 * conj(delta)) / σ2)
end

function demodulate_soft(syms::AbstractVector{<:Real}, scheme::Union{PSK,DBPSK,DQPSK}, noise_var::Real;
                          approx::Symbol=:exact)
    return demodulate_soft(ComplexF64.(syms), scheme, noise_var; approx=approx)
end

function demodulate_soft(syms::AbstractVector{<:Complex}, ::DBPSK, noise_var::Real;
                          approx::Symbol=:exact)
    _check_soft_args(noise_var, approx)
    n = length(syms)
    n >= 1 || throw(ArgumentError("DBPSK soft demodulation needs the leading reference symbol"))
    σ2 = float(noise_var)
    out = Vector{Float64}(undef, n - 1)
    @inbounds for i in 1:(n - 1)
        y0 = ComplexF64(syms[i])
        y1 = ComplexF64(syms[i + 1])
        # One phase per bit value, so exact and max-log are the same.
        out[i] = _diff_loglik(y0, y1, 1.0 + 0.0im, σ2) - _diff_loglik(y0, y1, -1.0 + 0.0im, σ2)
    end
    return out
end

function demodulate_soft(syms::AbstractVector{<:Complex}, ::DQPSK, noise_var::Real;
                          approx::Symbol=:exact)
    _check_soft_args(noise_var, approx)
    n = length(syms)
    n >= 1 || throw(ArgumentError("DQPSK soft demodulation needs the leading reference symbol"))
    σ2 = float(noise_var)
    out = Vector{Float64}(undef, 2 * (n - 1))
    @inbounds for i in 1:(n - 1)
        y0 = ComplexF64(syms[i])
        y1 = ComplexF64(syms[i + 1])
        ll1 = _diff_loglik(y0, y1, _DQPSK_DELTA[1], σ2)
        ll2 = _diff_loglik(y0, y1, _DQPSK_DELTA[2], σ2)
        ll3 = _diff_loglik(y0, y1, _DQPSK_DELTA[3], σ2)
        ll4 = _diff_loglik(y0, y1, _DQPSK_DELTA[4], σ2)
        # Index 2*b0+b1+1: 00, 01, 10, 11.
        if approx === :maxlog
            out[2i - 1] = max(ll1, ll2) - max(ll3, ll4)
            out[2i] = max(ll1, ll3) - max(ll2, ll4)
        else
            out[2i - 1] = _logaddexp(ll1, ll2) - _logaddexp(ll3, ll4)
            out[2i] = _logaddexp(ll1, ll3) - _logaddexp(ll2, ll4)
        end
    end
    return out
end

end # module Modulation
