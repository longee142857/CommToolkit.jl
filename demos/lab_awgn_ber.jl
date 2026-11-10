# Uncoded BPSK and Gray QPSK on AWGN. No plots.
#
# Run from the repository root:
#   julia --project=. demos/lab_awgn_ber.jl
#
# Noise is drawn from a local RNG, not from awgn (that function uses the
# global stream). Es = 1. With k bits per symbol, Eb = Es/k and
# N0 = Eb / (Eb/N0). Each real dimension has variance σ² = N0/2.
# BPSK is real, so its total noise power is σ². QPSK is complex, so its
# total noise power is N0.
#
# The terminated (7, 5) code is not in this script. conv_encode turns n
# information bits into 2(n+2) coded bits, so the effective rate is
# n/(2(n+2)), not 1/2.
#
# Pass rule, fixed before looking at a sample:
#   Eb/N0 ∈ {0, 2, 4} dB, 20_000 bits, both schemes
#   two-sided Wilson interval at α = 1e-3
#   theory lies inside the interval, and the error count is at least 30
# A count of 0 is printed with its upper bound and does not pass.

using Random
using Printf
using CommToolkit

const LAB_ALPHA = 1e-3
const LAB_NBITS = 20_000
const LAB_SNRS = (0.0, 2.0, 4.0)
const LAB_MIN_ERRORS = 30
const LAB_SEED = 0x0c0de1162

# Hastings rational approximation of erfc, relative error below 1.5e-7
# (Abramowitz and Stegun, 7.1.26). The Wilson interval below is orders of
# magnitude wider than that, at every point this script runs.
function _erfc(x::Float64)
    x < 0 && return 2 - _erfc(-x)
    t = inv(1 + 0.3275911 * x)
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 +
           t * (-1.453152027 + t * 1.061405429))))
    return poly * exp(-x * x)
end

"Bit error rate Q(√(2γ)) for uncoded BPSK and Gray QPSK at linear Eb/N0 = γ."
function ber_uncoded_awgn(γ::Float64)
    return 0.5 * _erfc(sqrt(γ))
end

# Solve erfc(z/√2) = α, i.e. the standard normal z with two-sided tail α.
# Bisection: erfc is decreasing, and a Newton step from a cold start overshoots the tail.
function _z_twosided(α::Float64)
    lo = 0.0
    hi = 6.0
    for _ in 1:60
        mid = (lo + hi) / 2
        if _erfc(mid) > α
            lo = mid
        else
            hi = mid
        end
    end
    return ((lo + hi) / 2) * sqrt(2)
end

"""Two-sided Wilson score interval for k errors in n bits."""
function wilson_interval(k::Int, n::Int, α::Float64)
    z = _z_twosided(α)
    phat = k / n
    z2 = z * z
    den = 1 + z2 / n
    centre = phat + z2 / (2n)
    half = z * sqrt(phat * (1 - phat) / n + z2 / (4 * n * n))
    return (centre - half) / den, (centre + half) / den
end

function _awgn_local(rng::AbstractRNG, syms::AbstractVector, ebn0_db::Float64, k::Int)
    ρ = exp10(ebn0_db / 10)
    Eb = 1 / k
    N0 = Eb / ρ
    σ = sqrt(N0 / 2)
    n = length(syms)
    if eltype(syms) <: Complex
        return syms .+ σ .* (randn(rng, n) .+ im .* randn(rng, n))
    else
        return syms .+ σ .* randn(rng, n)
    end
end

function _count_errors(rng::AbstractRNG, scheme, ebn0_db::Float64, nbits::Int)
    k = bits_per_symbol(scheme)
    nbits % k == 0 || throw(ArgumentError("nbits must be a multiple of bits per symbol"))
    bits = rand(rng, 0:1, nbits)
    rx = _awgn_local(rng, modulate(bits, scheme), ebn0_db, k)
    got = demodulate(rx, scheme)
    errors = count(!=(0), bits .- got)
    return errors, nbits
end

function _check_formula!()
    # Published values, not a second call of the same expression.
    # Q(1) = erfc(1/√2)/2, and BER at Eb/N0 = 0 dB is erfc(1)/2.
    abs(0.5 * _erfc(1 / sqrt(2)) - 0.15865525393145707) < 1e-7 ||
        error("erfc approximation missed Q(1)")
    abs(ber_uncoded_awgn(1.0) - 0.07864960352514258) < 1e-7 ||
        error("uncoded BER at 0 dB missed erfc(1)/2")
    z = _z_twosided(LAB_ALPHA)
    abs(z - 3.2905267314919255) < 1e-4 ||
        error("normal quantile for α=$LAB_ALPHA is $z, expected about 3.2905")
    lo, hi = wilson_interval(0, 1000, LAB_ALPHA)
    abs(lo) < 1e-15 && 0 < hi < 1 ||
        error("Wilson interval at 0 errors must be [0, hi) with hi in (0, 1)")
    return nothing
end

function main()
    _check_formula!()
    rng = Xoshiro(LAB_SEED)
    z = _z_twosided(LAB_ALPHA)
    println("Uncoded AWGN. Es = 1. Local Xoshiro seed ", string(LAB_SEED, base=16), ".")
    println("BPSK: Eb = Es. QPSK: Eb = Es/2. Per-dimension variance σ² = N0/2, N0 = Eb/(Eb/N0).")
    println("Wilson interval, α = $LAB_ALPHA, z = $(round(z; digits=4)). Pass requires k ≥ $LAB_MIN_ERRORS and theory inside the interval.")
    @printf("%-6s %8s %8s %8s %12s %12s %12s %12s  %s\n",
            "scheme", "Eb/N0", "n", "errors", "BER", "theory", "lo", "hi", "note")
    ok = true
    for scheme in (BPSK(), QPSK())
        name = scheme isa BPSK ? "BPSK" : "QPSK"
        for snr in LAB_SNRS
            errors, n = _count_errors(rng, scheme, snr, LAB_NBITS)
            ber_hat = errors / n
            theory = ber_uncoded_awgn(exp10(snr / 10))
            lo, hi = wilson_interval(errors, n, LAB_ALPHA)
            note = if errors == 0
                "k=0: hi is not a claim that BER is 0"
            elseif !(lo <= theory <= hi) || errors < LAB_MIN_ERRORS
                "FAIL"
            else
                "pass"
            end
            if note != "pass"
                ok = false
            end
            @printf("%-6s %8.1f %8d %8d %12.6g %12.6g %12.6g %12.6g  %s\n",
                    name, snr, n, errors, ber_hat, theory, lo, hi, note)
        end
    end
    println(ok ? "criterion: pass" : "criterion: fail")
    return ok ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
