# Monte Carlo markers against the closed forms used by test/reference.jl.
#
# From the repository root:
#   julia --project=test -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd()))'
#   GKSwstype=nul julia --project=test demos/validation_figure.jl
# Writes docs/src/assets/ber_validation.png.
# `test/Project.toml` names this unregistered package, so the develop step is required
# before a bare `Pkg.instantiate()` in that environment. `Pkg.test` does not need it.

using Random
using SpecialFunctions
using Plots
using CommToolkit

const SEED = 0x0c0de1162

qfunc(x::Real) = erfc(x / sqrt(2)) / 2
ber_bpsk_awgn(γb::Real) = qfunc(sqrt(2 * float(γb)))
ber_qpsk_awgn(γb::Real) = ber_bpsk_awgn(γb)
ber_bpsk_rayleigh(γb::Real) = 0.5 * (1 - sqrt(float(γb) / (1 + float(γb))))

# Cho and Yoon, IEEE Trans. Commun. 50(7):1074–1080, 2002, eq. (16).
function ber_square_mqam(M::Int, γb::Real)
    L = isqrt(M)
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

function mc_awgn(scheme, snr_db::Real, n_bits::Int)
    bits = rand(0:1, n_bits)
    k = bits_per_symbol(scheme)
    rx = awgn(modulate(bits, scheme), snr_db, k)
    return ber(bits, demodulate(rx, scheme))
end

function mc_rayleigh(snr_db::Real, n_bits::Int)
    bits = rand(0:1, n_bits)
    rx = rayleigh_channel(modulate(bits, BPSK()), snr_db, 1, true)
    return ber(bits, demodulate(rx, BPSK()))
end

function curve(snrs, theory, simulate)
    γ = exp10.(snrs ./ 10)
    return theory.(γ), [simulate(s) for s in snrs]
end

function main()
    Random.seed!(SEED)
    snr_psk = 0:2:8
    snr_qam = 0:2:12
    snr_ray = 0:4:20

    bpsk_th, bpsk_mc = curve(snr_psk, ber_bpsk_awgn, s -> mc_awgn(BPSK(), s, 200_000))
    qpsk_th, qpsk_mc = curve(snr_psk, ber_qpsk_awgn, s -> mc_awgn(QPSK(), s, 200_000))
    qam_th, qam_mc = curve(snr_qam, γ -> ber_square_mqam(16, γ), s -> mc_awgn(QAM(16), s, 200_000))
    ray_th, ray_mc = curve(snr_ray, ber_bpsk_rayleigh, s -> mc_rayleigh(s, 100_000))

    plt = plot(;
        yaxis = :log10,
        xlabel = "Eb/N0 (dB)",
        ylabel = "Bit error rate",
        legend = :bottomleft,
        size = (760, 480),
        framestyle = :box,
    )
    plot!(plt, snr_psk, bpsk_th; label = "BPSK AWGN, Q(√(2γb))", color = :black, lw = 1.5)
    scatter!(plt, snr_psk, bpsk_mc; label = "BPSK AWGN, Monte Carlo", color = :black, marker = :circle)
    plot!(plt, snr_psk, qpsk_th; label = "Gray QPSK AWGN", color = :gray, lw = 1.5, ls = :dash)
    scatter!(plt, snr_psk, qpsk_mc; label = "Gray QPSK, Monte Carlo", color = :gray, marker = :diamond)
    plot!(plt, snr_qam, qam_th; label = "16-QAM AWGN, Cho–Yoon", color = 1, lw = 1.5)
    scatter!(plt, snr_qam, qam_mc; label = "16-QAM, Monte Carlo", color = 1, marker = :square)
    plot!(plt, snr_ray, ray_th; label = "BPSK Rayleigh, perfect CSI", color = 2, lw = 1.5)
    scatter!(plt, snr_ray, ray_mc; label = "BPSK Rayleigh, Monte Carlo", color = 2, marker = :utriangle)

    out = joinpath(@__DIR__, "..", "docs", "src", "assets", "ber_validation.png")
    mkpath(dirname(out))
    savefig(plt, out)
    println("wrote ", abspath(out))
    return out
end

main()
