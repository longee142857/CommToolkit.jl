module Visualization

import ..Modulation

export raised_cosine, rrc_filter, pulse_shape,
    plot_constellation, plot_eye_diagram, plot_pulse_response

function _sinc(x::Real)
    return x == 0 ? 1.0 : sin(π * x) / (π * x)
end

function _plots_extension()
    ext = Base.get_extension(parentmodule(@__MODULE__), :CommToolkitPlotsExt)
    if ext === nothing
        error("plotting requires Plots.jl. Load it with `using Plots` so the CommToolkitPlotsExt extension is activated. `using CommToolkit` alone does not load Plots.")
    end
    return ext
end

@doc raw"""
    raised_cosine(beta::Real, span::Integer, sps::Integer) -> Vector{Float64}

Samples of the unit-symbol-period raised-cosine pulse (Proakis and Salehi, *Digital Communications*, 5th ed., §9.2),

```math
h(\tau) = \operatorname{sinc}(\tau)\,
\frac{\cos(\pi \beta \tau)}{1 - (2\beta\tau)^2},
\qquad \tau = t/T,
```

with the usual removable values ``h(0) = 1`` and ``h(\pm 1/(2\beta)) = (\pi/4)\operatorname{sinc}(1/(2\beta))``. `span` is the one-sided length in symbols (the support runs from `-span/2` to `span/2`). `sps` is the number of samples per symbol. `beta = 0` returns the sinc pulse. The taps are scaled to unit DC gain (`sum(taps) = 1`).

`beta` must lie in ``[0, 1]``.
"""
function raised_cosine(beta::Real, span::Integer, sps::Integer)
    0 <= beta <= 1 || throw(ArgumentError("beta must lie in [0, 1], got $beta"))
    span >= 1 || throw(ArgumentError("span must be ≥ 1"))
    sps >= 1 || throw(ArgumentError("sps must be ≥ 1"))
    β = float(beta)
    t = range(-span / 2, span / 2; length=span * sps + 1)
    taps = Vector{Float64}(undef, length(t))
    for (i, τ) in enumerate(t)
        if β == 0 || τ == 0
            taps[i] = _sinc(τ)
        elseif abs(abs(τ) - 1 / (2 * β)) < 1e-12
            taps[i] = (π / 4) * _sinc(1 / (2 * β))
        else
            taps[i] = _sinc(τ) * cos(π * β * τ) / (1 - (2 * β * τ)^2)
        end
    end
    return taps ./ sum(taps)
end

@doc raw"""
    rrc_filter(beta::Real, span::Integer, sps::Integer) -> Vector{Float64}

Samples of the root-raised-cosine pulse with roll-off `beta`, support `span` symbols, and `sps` samples per symbol (Proakis and Salehi, 5th ed., §9.2). The continuous-time shape, before truncation, is

```math
g(\tau) =
\begin{cases}
\dfrac{1}{\sqrt{T}}\bigl(1 - \beta + 4\beta/\pi\bigr), & \tau = 0, \\[1em]
\dfrac{\beta}{\sqrt{2T}}\left[(1+2/\pi)\sin\dfrac{\pi}{4\beta} + (1-2/\pi)\cos\dfrac{\pi}{4\beta}\right],
  & |\tau| = 1/(4\beta), \\[1em]
\dfrac{1}{\sqrt{T}}
  \dfrac{\sin\bigl(\pi\tau(1-\beta)\bigr) + 4\beta\tau\cos\bigl(\pi\tau(1+\beta)\bigr)}
        {\pi\tau\bigl(1 - (4\beta\tau)^2\bigr)},
  & \text{otherwise,}
\end{cases}
```

with ``\tau = t/T``. `beta = 0` is the sinc pulse. Taps are then scaled to unit energy, ``\sum g_n^2 = 1``. A cascade of two such filters, sampled at symbol instants about the peak, is approximately Nyquist; the residual ISI is the truncation error of the finite `span`.
"""
function rrc_filter(beta::Real, span::Integer, sps::Integer)
    0 <= beta <= 1 || throw(ArgumentError("beta must lie in [0, 1], got $beta"))
    span >= 1 || throw(ArgumentError("span must be ≥ 1"))
    sps >= 1 || throw(ArgumentError("sps must be ≥ 1"))
    β = float(beta)
    t = range(-span / 2, span / 2; length=span * sps + 1)
    taps = Vector{Float64}(undef, length(t))
    for (i, τ) in enumerate(t)
        if β == 0
            taps[i] = _sinc(τ)
        elseif τ == 0
            taps[i] = (1 - β + 4 * β / π) / sqrt(sps)
        elseif abs(abs(τ) - 1 / (4 * β)) < 1e-12
            a = β / sqrt(2)
            b = (1 + 2 / π) * sin(π / (4 * β))
            c = (1 - 2 / π) * cos(π / (4 * β))
            taps[i] = a * (b + c) / sqrt(sps)
        else
            num = sin(π * τ * (1 - β)) + 4 * β * τ * cos(π * τ * (1 + β))
            den = π * τ * (1 - (4 * β * τ)^2)
            taps[i] = num / den / sqrt(sps)
        end
    end
    return taps ./ sqrt(sum(abs2, taps))
end

@doc raw"""
    pulse_shape(symbols, taps::AbstractVector{<:Real}, sps::Integer)

Upsample `symbols` by inserting `sps - 1` zeros after every symbol, convolve with `taps`, and **truncate the result to `length(symbols) * sps` samples** (the filter tail of length `length(taps) - 1` is discarded).

The truncation keeps the start of the convolution, including the filter's group delay, and drops the tail. It is the waveform used by the eye-diagram helper. It does not apply a matched filter or remove the group delay; that is [`CommToolkit.Pulse.matched_filter`](@ref) followed by [`CommToolkit.Pulse.downsample_symbols`](@ref). Because the tail is already gone, a later matched filter is not the full cascade on the last symbols of the burst.
"""
function pulse_shape(symbols::AbstractVector{<:Number}, taps::AbstractVector{<:Real}, sps::Integer)
    sps >= 1 || throw(ArgumentError("sps must be ≥ 1"))
    up = zeros(promote_type(eltype(symbols), eltype(taps), Float64), length(symbols) * sps)
    up[1:sps:end] .= symbols
    n, m = length(up), length(taps)
    full = zeros(eltype(up), n + m - 1)
    @inbounds for i in 1:n, j in 1:m
        full[i + j - 1] += up[i] * taps[j]
    end
    return full[1:n]
end

@doc raw"""
    plot_constellation(rx_syms::AbstractVector{<:Complex}, scheme; title::String="Constellation Diagram")

Scatter `rx_syms` against [`Modulation.constellation`](@ref)`(scheme)`.

Requires Plots.jl. Loading Plots activates the package extension `CommToolkitPlotsExt`; `using CommToolkit` by itself does not import Plots. Returns the `Plots.Plot`.
"""
function plot_constellation(rx_syms::AbstractVector{<:Complex}, scheme; title::String="Constellation Diagram")
    return _plots_extension().plot_constellation(rx_syms, scheme; title=title)
end

@doc raw"""
    plot_eye_diagram(signal::AbstractVector{<:Real}, sps::Integer; n_segments::Integer=100, title::String="Eye Diagram")

Overlay `n_segments` successive blocks of `2 * sps` samples. `sps` is the number of samples per symbol. Requires Plots.jl; see [`plot_constellation`](@ref).
"""
function plot_eye_diagram(signal::AbstractVector{<:Real}, sps::Integer;
                          n_segments::Integer=100, title::String="Eye Diagram")
    return _plots_extension().plot_eye_diagram(signal, Int(sps); n_segments=Int(n_segments), title=title)
end

@doc raw"""
    plot_pulse_response(taps::AbstractVector{<:Real}; title::String="Pulse Shape Response")

Time-domain taps and the magnitude of a length-1024 FFT on one figure. Requires Plots.jl; see [`plot_constellation`](@ref). The FFT is unnormalised (`FFTW.fft`).
"""
function plot_pulse_response(taps::AbstractVector{<:Real}; title::String="Pulse Shape Response")
    return _plots_extension().plot_pulse_response(taps; title=title)
end

end # module Visualization
