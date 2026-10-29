module CommToolkitPlotsExt

using CommToolkit
using Plots
import FFTW

function plot_constellation(rx_syms::AbstractVector{<:Complex}, scheme; title::String="Constellation Diagram")
    ideal = constellation(scheme)
    plt = plot(real.(rx_syms), imag.(rx_syms);
        seriestype=:scatter, markersize=3, color=:blue,
        label="Received", legend=:topright)
    scatter!(plt, real.(ideal), imag.(ideal);
        markersize=8, color=:red, markerstrokewidth=2,
        label="Ideal", alpha=0.6)
    plot!(plt; xlabel="In-Phase", ylabel="Quadrature",
        title=title, aspect_ratio=:equal,
        grid=true, framestyle=:box)
    return plt
end

function plot_eye_diagram(signal::AbstractVector{<:Real}, sps::Int;
                          n_segments::Int=100, title::String="Eye Diagram")
    sps >= 1 || throw(ArgumentError("sps must be ≥ 1"))
    seg_len = 2 * sps
    n_total = length(signal)
    n_usable = (n_total ÷ seg_len) * seg_len
    n_plot = min(n_segments, n_usable ÷ seg_len)
    plt = plot(; xlabel="Time (symbols)", ylabel="Amplitude",
        title=title, grid=true, framestyle=:box, legend=false)
    for i in 1:n_plot
        idx = (i - 1) * seg_len + 1
        seg = signal[idx:(idx + seg_len - 1)]
        t = (0:(seg_len - 1)) ./ sps
        plot!(plt, t, seg; color=:blue, alpha=0.3, lw=0.5)
    end
    return plt
end

function plot_pulse_response(taps::AbstractVector{<:Real}; title::String="Pulse Shape Response")
    p1 = plot(1:length(taps), taps;
        lw=2, color=:blue,
        xlabel="Tap index", ylabel="Amplitude",
        title="$title (Time)", grid=true, legend=false)
    nfft = 1024
    # `fft(x, n)` is not a padded transform in AbstractFFTs; the integer is a dimension.
    buf = zeros(ComplexF64, nfft)
    ncopy = min(nfft, length(taps))
    buf[1:ncopy] .= taps[1:ncopy]
    H = FFTW.fft(buf)
    freq = (0:(nfft - 1)) ./ nfft
    half = 1:(nfft ÷ 2)
    p2 = plot(freq[half], abs.(H[half]);
        lw=2, color=:red,
        xlabel="Normalized Frequency", ylabel="|H(f)|",
        title="$title (Frequency)", grid=true, legend=false)
    return plot(p1, p2; layout=(1, 2), size=(800, 350))
end

function plot_exit_chart(curve::EXITCurve; title::String="EXIT chart")
    plt = plot(curve.ia, curve.ie;
        lw=2, color=:blue, marker=:circle, markersize=4,
        label="measured T(I_A)", legend=:bottomright)
    plot!(plt; xlabel="I_A", ylabel="I_E",
        xlims=(0, 1), ylims=(0, 1), aspect_ratio=:equal,
        title=title, grid=true, framestyle=:box)
    return plt
end

function plot_exit_chart(demapper::EXITCurve, decoder::EXITCurve;
                         title::String="EXIT chart")
    plt = plot(demapper.ia, demapper.ie;
        lw=2, color=:blue, marker=:circle, markersize=4,
        label="demapper I_E = T(I_A)", legend=:bottomright)
    plot!(plt, decoder.ie, decoder.ia;
        lw=2, color=:red, marker=:square, markersize=4,
        label="decoder I_A = T(I_E), swapped axes")
    plot!(plt; xlabel="I_A demapper / I_E decoder", ylabel="I_E demapper / I_A decoder",
        xlims=(0, 1), ylims=(0, 1), aspect_ratio=:equal,
        title=title, grid=true, framestyle=:box)
    return plt
end

end # module CommToolkitPlotsExt
