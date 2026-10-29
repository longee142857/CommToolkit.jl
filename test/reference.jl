# Closed forms used by the Monte Carlo tests.
#
# BPSK / Gray QPSK in AWGN: Proakis and Salehi, Digital Communications, 5th ed., §5.2.
# Square Gray M-QAM: Cho and Yoon, IEEE Trans. Commun. 50(7):1074–1080, 2002, eq. (16)
#   (average of the per-bit PAM error rates; erfc argument matches their normalisation).
# Coherent BPSK in flat Rayleigh fading, perfect CSI: Goldsmith, Wireless Communications,
#   Cambridge, 2005, §6.3; also Proakis and Salehi, 5th ed., §13.3.
# Gray M-PSK: phase-tail sector probabilities, Proakis and Salehi, 5th ed., §4.2 and §5.2.
#   The Craig form of the symbol error rate is the cross-check. The erfc expression is
#   the nearest-neighbour approximation and is exact only for M = 4.
# DBPSK: Pb = (1/2) exp(−Eb/N0), Proakis and Salehi, 5th ed., §4.5.
# Gray DQPSK: Marcum-Q bit error rate, Proakis and Salehi, 5th ed., §4.5.
# Coherent BPSK in flat Rice, perfect CSI: Craig integral of the Rician MGF,
#   Proakis and Salehi, 5th ed., §13.3. K is linear.
# Rayleigh MRC, coherent BPSK, L i.i.d. branches: Proakis and Salehi, 5th ed., §13.4.
# Nakagami-m coherent BPSK, perfect CSI: Craig integral of the Gamma MGF.
#   Integer m equals MRC with L = m and per-branch Eb/N0 = γb/m.
# Hard-decision union bound: Proakis and Salehi, 5th ed., §8.2; Lin and Costello,
#   Error Control Coding, 2nd ed., §12.3. c_d comes from Coding.distance_spectrum.
# Soft-decision pairwise error on antipodal AWGN: Proakis and Salehi, 5th ed., §8.2,
#   P2(d) = Q(sqrt(2 d R Eb/N0)) with R = 1/2.
# Clopper–Pearson interval: Clopper and Pearson, Biometrika 26(4):404–413, 1934.
# BIAWGN consistent-LLR mutual information (1.1.0, EXIT tests): the identity
#   I = 1 − E[log2(1 + exp(−L))] on L = 2y/σ², y ~ N(1, σ²), Hagenauer,
#   Proc. 12th EUSIPCO, Vienna, 2004, pp. 1541–1548; histogram form ten Brink,
#   IEEE Trans. Commun. 49(10):1727–1737, 2001. No elementary closed form:
#   numerical Simpson integration, the exact value that exit_j(2/σ) approximates.

using SpecialFunctions

const BER_ALPHA = 1e-4
const RNG_SEED = 0x0c0de1162

qfunc(x::Real) = erfc(x / sqrt(2)) / 2

ber_bpsk_awgn(γb::Real) = qfunc(sqrt(2 * float(γb)))
ber_qpsk_awgn(γb::Real) = ber_bpsk_awgn(γb)

"""Exact average BER of square Gray M-QAM at linear Eb/N0 = γb (Cho and Yoon, 2002)."""
function ber_square_mqam(M::Int, γb::Real)
    L = isqrt(M)
    L^2 == M && ispow2(M) || throw(ArgumentError("M must be a square power of two"))
    m = round(Int, log2(L))
    k = 2m
    γb = float(γb)
    acc = 0.0
    for j in 1:m
        sj = 0.0
        pj = 1 << (j - 1)
        imax = (L - (L >> j)) - 1
        for i in 0:imax
            sgn = iseven(fld(i * pj, L)) ? 1.0 : -1.0
            mag = pj - floor(Int, i * pj / L + 0.5)
            arg = (2i + 1) * sqrt(3 * k * γb / (2 * (M - 1)))
            sj += sgn * mag * erfc(arg)
        end
        acc += sj
    end
    return acc / (L * m)
end

ber_bpsk_rayleigh(γb::Real) = 0.5 * (1 - sqrt(float(γb) / (1 + float(γb))))

"""BIAWGN capacity at per-component noise variance σ², unit symbol energy.
The mutual information of the consistent LLR L = 2y/σ², y ~ N(1, σ²), and the
exact integral that the Gaussian-symmetric approximation `exit_j(2/√σ²)`
approximates. Numerical integration, not a closed form."""
function biawgn_llr_mi(σ2::Real)
    σ2 = float(σ2)
    σ2 > 0 || throw(ArgumentError("noise variance must be positive"))
    σ = sqrt(σ2)
    invln2 = inv(log(2))
    function integrand(y)
        x = -2y / σ2
        sp = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))
        return (1 - sp * invln2) * exp(-abs2(y - 1) / (2σ2)) / (σ * sqrt(2π))
    end
    return _simpson(integrand, 1 - 12σ, 1 + 12σ, 8192)
end

function _simpson(f, a::Real, b::Real, n::Int)
    a = float(a)
    b = float(b)
    b == a && return 0.0
    n += n % 2
    h = (b - a) / n
    s = f(a) + f(b)
    for i in 1:(n - 1)
        s += (isodd(i) ? 4 : 2) * f(a + i * h)
    end
    return s * h / 3
end

# P(phase of √γs + CN(0,1) noise exceeds φ), φ ∈ [−π, π].
# Equivalent to (1/2π) ∫_0^{π−φ} exp(−γs sin²φ / sin²θ) dθ for φ ≥ 0,
# which is the tail behind Craig's M-PSK symbol-error integral.
function _phase_tail(φ::Real, γs::Real)
    φ = float(φ)
    γs = float(γs)
    if φ <= -π
        return 1.0
    elseif φ >= π
        return 0.0
    elseif φ < 0
        return 1.0 - _phase_tail(-φ, γs)
    end
    upper = π - φ
    upper <= 0 && return 0.0
    s2 = sin(φ)^2
    integrand(θ) = begin
        sθ = sin(θ)
        if sθ == 0
            return (γs == 0 || s2 == 0) ? 1.0 : 0.0
        end
        return exp(-γs * s2 / (sθ * sθ))
    end
    return _simpson(integrand, 0.0, upper, 2048) / (2π)
end

function _offset_prob(d::Int, M::Int, γs::Real)
    lo = (d - 0.5) * 2π / M
    hi = (d + 0.5) * 2π / M
    while lo > π
        lo -= 2π
        hi -= 2π
    end
    while lo <= -π
        lo += 2π
        hi += 2π
    end
    if hi <= π
        return _phase_tail(lo, γs) - _phase_tail(hi, γs)
    end
    return (_phase_tail(lo, γs) - _phase_tail(π, γs)) +
           (_phase_tail(-π, γs) - _phase_tail(hi - 2π, γs))
end

"""P(nearest-phase decision is natural index tx+d), d = 0..M−1. Rows sum to 1."""
function psk_offset_probs(M::Integer, γs::Real)
    M = Int(M)
    M >= 2 || throw(ArgumentError("M must be ≥ 2"))
    return [_offset_prob(d, M, γs) for d in 0:(M - 1)]
end

"""Exact SER of coherent M-PSK from the phase sectors. γs is Es/N0."""
function ser_mpsk(M::Integer, γs::Real)
    return 1 - psk_offset_probs(M, γs)[1]
end

"""Craig SER of coherent M-PSK. γs is Es/N0. Proakis and Salehi, 5th ed., §5.2."""
function ser_mpsk_craig(M::Integer, γs::Real)
    M = Int(M)
    γs = float(γs)
    ψ = π / M
    upper = (M - 1) * π / M
    integrand(θ) = begin
        sθ = sin(θ)
        sθ == 0 && return 0.0
        return exp(-γs * sin(ψ)^2 / (sθ * sθ))
    end
    return _simpson(integrand, 0.0, upper, 2048) / π
end

"""Exact average BER of binary-reflected Gray M-PSK. γb is Eb/N0, Es = log2(M) Eb."""
function ber_mpsk_gray(M::Integer, γb::Real)
    M = Int(M)
    (M >= 2 && ispow2(M)) || throw(ArgumentError("M must be a power of two ≥ 2"))
    γb = float(γb)
    γb >= 0 || throw(ArgumentError("γb must be ≥ 0"))
    k = trailing_zeros(M)
    γs = k * γb
    γs == 0 && return 0.5
    p = psk_offset_probs(M, γs)
    gray = [xor(m, m >> 1) for m in 0:(M - 1)]
    acc = 0.0
    for tx in 0:(M - 1), d in 0:(M - 1)
        rx = mod(tx + d, M)
        acc += count_ones(xor(gray[tx + 1], gray[rx + 1])) * p[d + 1]
    end
    return acc / (M * k)
end

"""Nearest-neighbour Gray M-PSK BER. Exact for M = 4; an approximation for M > 4."""
function ber_mpsk_approx(M::Integer, γb::Real)
    M = Int(M)
    (M >= 4 && ispow2(M)) || throw(ArgumentError("M must be a power of two ≥ 4"))
    k = float(trailing_zeros(M))
    return erfc(sqrt(k * float(γb)) * sin(π / M)) / k
end

ber_dbpsk(γb::Real) = 0.5 * exp(-float(γb))

"""First-order Marcum Q, Q₁(a, b) = ∫_b^∞ x exp(−(x−a)²/2) I₀(ax) e^{−ax} dx, via besselix."""
function marcumq1(a::Real, b::Real)
    a = float(a)
    b = float(b)
    (a >= 0 && b >= 0) || throw(ArgumentError("Marcum Q arguments must be ≥ 0"))
    b == 0 && return 1.0
    a == 0 && return exp(-0.5 * b * b)
    lo = b
    hi = max(a, b) + 12.0
    integrand(x) = x * exp(-0.5 * (x - a)^2) * besselix(0, a * x)
    return _simpson(integrand, lo, hi, 2000)
end

"""Gray DQPSK bit error rate. γb is Eb/N0 with Es = 2 Eb. Proakis and Salehi, 5th ed., §4.5."""
function ber_dqpsk(γb::Real)
    γb = float(γb)
    γb >= 0 || throw(ArgumentError("γb must be ≥ 0"))
    γb == 0 && return 0.5
    isinf(γb) && return 0.0
    s = sqrt(0.5)
    a = sqrt(2 * γb * (1 - s))
    b = sqrt(2 * γb * (1 + s))
    # I₀(ab) exp(−(a²+b²)/2) = besselix(0, ab) exp(−(a−b)²/2)
    bessel_term = besselix(0, a * b) * exp(-0.5 * (a - b)^2)
    return marcumq1(a, b) - 0.5 * bessel_term
end

"""Coherent BPSK BER in flat Rice fading, perfect CSI. K is linear; γb is average Eb/N0."""
function ber_bpsk_rice(K::Real, γb::Real)
    K = float(K)
    γb = float(γb)
    (K >= 0 && γb >= 0) || throw(ArgumentError("K and γb must be ≥ 0"))
    γb == 0 && return 0.5
    isinf(γb) && return 0.0
    isinf(K) && return ber_bpsk_awgn(γb)
    integrand(θ) = begin
        s2 = sin(θ)^2
        den = (1 + K) * s2 + γb
        den == 0 && return 0.0
        return ((1 + K) * s2 / den) * exp(-K * γb / den)
    end
    return _simpson(integrand, 0.0, π / 2, 2048) / π
end

"""Coherent BPSK BER, L-branch i.i.d. Rayleigh MRC. γb is the average Eb/N0 per branch."""
function ber_bpsk_mrc_rayleigh(L::Integer, γb::Real)
    L = Int(L)
    L >= 1 || throw(ArgumentError("L must be ≥ 1"))
    γb = float(γb)
    γb >= 0 || throw(ArgumentError("γb must be ≥ 0"))
    isinf(γb) && return 0.0
    μ = sqrt(γb / (1 + γb))
    α = (1 - μ) / 2
    β = (1 + μ) / 2
    s = 0.0
    for k in 0:(L - 1)
        s += float(binomial(L - 1 + k, k)) * β^k
    end
    return α^L * s
end

"""Coherent BPSK BER in i.i.d. Nakagami-m fading, perfect CSI.

`m ≥ 1/2` and `γb` is the average Eb/N0. Instantaneous SNR is Gamma(shape = m, scale = γb/m), with Ω = E[|h|²] = 1. Craig form of `Q` against that MGF. `m = 1` is Rayleigh. Positive integer `m` equals `ber_bpsk_mrc_rayleigh(m, γb/m)`. `m = Inf` is AWGN.
"""
function ber_bpsk_nakagami(m::Real, γb::Real)
    m = float(m)
    γb = float(γb)
    (m >= 0.5 && γb >= 0) || throw(ArgumentError("m must be ≥ 1/2 and γb must be ≥ 0"))
    γb == 0 && return 0.5
    isinf(γb) && return 0.0
    isinf(m) && return ber_bpsk_awgn(γb)
    integrand(θ) = begin
        s2 = sin(θ)^2
        den = m * s2 + γb
        den == 0 && return 0.0
        return ((m * s2) / den)^m
    end
    return _simpson(integrand, 0.0, π / 2, 2048) / π
end

# SpecialFunctions 2.9 returns `(x, 1 - x)` from `beta_inc_inv`; earlier 2.x releases
# returned the quantile alone. Keep the first element either way.
function _beta_quantile(a::Real, b::Real, p::Real)
    result = beta_inc_inv(float(a), float(b), float(p))
    return result isa Tuple ? float(first(result)) : float(result)
end

"""Two-sided Clopper–Pearson interval for k successes in n Bernoulli trials."""
function clopper_pearson(k::Int, n::Int, α::Float64=BER_ALPHA)
    0 <= k <= n || throw(ArgumentError("need 0 ≤ k ≤ n"))
    0 < α < 1 || throw(ArgumentError("α must lie in (0, 1)"))
    lower = k == 0 ? 0.0 : _beta_quantile(k, n - k + 1, α / 2)
    upper = k == n ? 1.0 : _beta_quantile(k + 1, n - k, 1 - α / 2)
    return lower, upper
end

"""Pairwise error probability of a weight-d path on a BSC with crossover p."""
function pairwise_bsc(d::Int, p::Float64)
    0 < p < 0.5 || throw(ArgumentError("BSC crossover must lie in (0, 1/2)"))
    s = 0.0
    if iseven(d)
        k0 = d ÷ 2
        s += exp(-log(2) + _logbinom(d, k0) + k0 * log(p) + (d - k0) * log(1 - p))
        start = k0 + 1
    else
        start = d ÷ 2 + 1
    end
    for k in start:d
        s += exp(_logbinom(d, k) + k * log(p) + (d - k) * log(1 - p))
    end
    return s
end

function _logbinom(n::Int, k::Int)
    return log(binomial(n, k))
end

"""Truncated soft-decision union bound at linear information-bit Eb/N0 = γb.

Pairwise term `Q(√(2 d R γb))` with `R = 1/2`, i.e. `Q(√(d γb))`, times `c_d`.
Proakis and Salehi, 5th ed., §8.2. The rate assumes `Es = 1` on coded BPSK and
`awgn(..., energy_per_bit = 2)`.
"""
function conv_soft_union_bound(γb::Real; dmax::Int=24)
    _, cd = distance_spectrum(dmax)
    γb = float(γb)
    s = 0.0
    for d in 1:dmax
        cd[d] == 0 && continue
        s += cd[d] * qfunc(sqrt(d * γb))
    end
    return s
end

"""Truncated hard-decision union bound using the trellis spectrum up to `dmax`."""
function conv_union_bound(p::Float64; dmax::Int=24)
    _, cd = distance_spectrum(dmax)
    s = 0.0
    for d in 1:dmax
        cd[d] == 0 && continue
        s += cd[d] * pairwise_bsc(d, p)
    end
    return s
end

function assert_theory_inside(k::Int, n::Int, p_theory::Float64;
                              min_errors::Int=200, α::Float64=BER_ALPHA)
    lo, hi = clopper_pearson(k, n, α)
    if k < min_errors
        @error "too few errors for an informative interval" k n min_errors p_theory
    end
    if !(lo <= p_theory <= hi)
        @error "theory BER outside Clopper–Pearson interval" k n p_hat=k / n p_theory lo hi α
    end
    @test k >= min_errors
    @test lo <= p_theory <= hi
    return lo, hi
end
