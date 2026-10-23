# Clarke/Jakes lag correlation against J0, and coherent BPSK on Nakagami-m.
# J0 is the power series. The Nakagami curve is the Craig integral of the Gamma MGF.
# No Plots. The stream is not seeded.

using CommToolkit
using Printf

function bessel_j0(x::Real)
    x = float(x)
    t2 = (0.5 * x)^2
    term = 1.0
    s = 1.0
    for k in 1:40
        term *= -t2 / (k * k)
        s += term
        abs(term) < 1e-15 * max(1.0, abs(s)) && break
    end
    return s
end

function simpson(f, a, b, n)
    n += n % 2
    h = (b - a) / n
    s = f(a) + f(b)
    for i in 1:(n - 1)
        s += (isodd(i) ? 4 : 2) * f(a + i * h)
    end
    return s * h / 3
end

function ber_nakagami(m, γ)
    # Craig form. m = Inf is Q(√(2γ)), the same limit as the Gamma MGF.
    f = if isinf(m)
        θ -> begin
            s = sin(θ)
            s == 0 && return 0.0
            return exp(-γ / (s * s))
        end
    else
        θ -> ((m * sin(θ)^2) / (m * sin(θ)^2 + γ))^m
    end
    return simpson(f, 0.0, π / 2, 2048) / π
end

N = 8192
fd = 0.02
h = jakes_fading(N, fd, 1.0)
println("jakes_fading(N=$N, fd=$fd Hz, Ts=1 s), one record")
for lag in (0, 8, 16, 40)
    s = sum(h[t + lag] * conj(h[t]) for t in 1:(N - lag)) / (N - lag)
    j0 = bessel_j0(2π * fd * lag)
    @printf("  lag %-3d   Re R %.4f   J0 %.4f\n", lag, real(s), j0)
end

snr_db = 8.0
γ = 10^(snr_db / 10)
nbits = 80_000
println("Nakagami-m coherent BPSK, Eb/N0 = ", snr_db, " dB, bits = ", nbits)
for m in (0.5, 1.0, 2.0, Inf)
    bits = rand(0:1, nbits)
    rx = nakagami_channel(modulate(bits, BPSK()), snr_db, 1, true; m=m)
    sim = count(!=(0), bits .- demodulate(rx, BPSK())) / nbits
    @printf("  m = %-6s  sim %.4e   formula %.4e\n", string(m), sim, ber_nakagami(m, γ))
end
