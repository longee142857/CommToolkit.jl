# Closed forms used by the Monte Carlo tests.
#
# BPSK / Gray QPSK in AWGN: Proakis and Salehi, Digital Communications, 5th ed., §5.2.
# Square Gray M-QAM: Cho and Yoon, IEEE Trans. Commun. 50(7):1074–1080, 2002, eq. (16)
#   (average of the per-bit PAM error rates; erfc argument matches their normalisation).
# Coherent BPSK in flat Rayleigh fading, perfect CSI: Goldsmith, Wireless Communications,
#   Cambridge, 2005, §6.3; also Proakis and Salehi, 5th ed., §13.3.
# Hard-decision union bound: Proakis and Salehi, 5th ed., §8.2; Lin and Costello,
#   Error Control Coding, 2nd ed., §12.3. c_d comes from Coding.distance_spectrum.
# Soft-decision pairwise error on antipodal AWGN: Proakis and Salehi, 5th ed., §8.2,
#   P2(d) = Q(sqrt(2 d R Eb/N0)) with R = 1/2.
# Clopper–Pearson interval: Clopper and Pearson, Biometrika 26(4):404–413, 1934.

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
