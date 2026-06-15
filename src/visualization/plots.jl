# Visualization utilities: constellation plots, eye diagrams

import Plots

function _get_ideal_constellation(scheme)
    if scheme isa BPSK
        return [1.0+0im, -1.0+0im]
    elseif scheme isa QPSK
        return [1+1im, 1-1im, -1+1im, -1-1im] ./ sqrt(2)
    elseif scheme isa QAM
        m = scheme.m
        n = Int(sqrt(m))
        pam = [-(n-1):2:(n-1)...] ./ sqrt(2 * (m - 1) / 3)
        return [pam[i] + pam[j]*im for i in 1:n, j in 1:n][:]
    end
end

# Plot constellation: received symbols vs ideal
function plot_constellation(rx_syms::AbstractVector{<:Complex}, scheme;
                            title::String="Constellation Diagram")
    ideal = _get_ideal_constellation(scheme)

    p = Plots.plot(real(rx_syms), imag(rx_syms),
        seriestype=:scatter, markersize=3, color=:blue,
        label="Received", legend=:topright)
    Plots.scatter!(p, real(ideal), imag(ideal),
        markersize=8, color=:red, markerstrokewidth=2,
        label="Ideal", alpha=0.6)
    Plots.plot!(p, xlabel="In-Phase", ylabel="Quadrature",
        title=title, aspect_ratio=:equal,
        grid=true, framestyle=:box)
    p
end

# Eye diagram: overlay segments of a signal
function plot_eye_diagram(signal::AbstractVector{<:Real}, sps::Int;
                          n_segments::Int=100, title::String="Eye Diagram")
    seg_len = 2 * sps
    n_total = length(signal)
    n_usable = n_total ÷ seg_len * seg_len
    n_plot = min(n_segments, n_usable ÷ seg_len)

    p = Plots.plot(xlabel="Time (symbols)", ylabel="Amplitude",
        title=title, grid=true, framestyle=:box, legend=false)

    for i in 1:n_plot
        idx = (i-1) * seg_len + 1
        seg = signal[idx:idx+seg_len-1]
        t = (0:seg_len-1) ./ sps
        Plots.plot!(p, t, seg, color=:blue, alpha=0.3, lw=0.5)
    end
    p
end

# Plot pulse shape impulse response and frequency response
function plot_pulse_response(taps::AbstractVector{<:Real};
                              title::String="Pulse Shape Response")

    p1 = Plots.plot(1:length(taps), taps,
        lw=2, color=:blue,
        xlabel="Tap index", ylabel="Amplitude",
        title="$title (Time)", grid=true, legend=false)

    nfft = 1024
    H = fft(taps, nfft)
    freq = (0:nfft-1) ./ nfft
    p2 = Plots.plot(freq[1:nfft÷2], abs.(H[1:nfft÷2]),
        lw=2, color=:red,
        xlabel="Normalized Frequency", ylabel="|H(f)|",
        title="$title (Frequency)", grid=true, legend=false)

    Plots.plot(p1, p2, layout=(1, 2), size=(800, 350))
end
