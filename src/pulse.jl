"""
    Pulse

Matched filter for a transmit pulse, and symbol-rate sampling at a caller-chosen phase. This module does not choose the sampling phase or the carrier. Gardner timing, a Costas loop, and Schmidl–Cox live in `Sync` and call these two functions; they do not change them.
"""
module Pulse

import ..Channel

export matched_filter, downsample_symbols

@doc raw"""
    matched_filter(samples, taps) -> Vector

Filter `samples` with the matched filter of `taps`: full linear convolution with `conj.(reverse(taps))`, the same sum as [`Channel.apply_fir`](@ref).

The output length is `length(samples) + length(taps) - 1`. The tail is kept. No group delay is removed and no boundary samples are dropped. For a real symmetric pulse ([`CommToolkit.Visualization.raised_cosine`](@ref), [`CommToolkit.Visualization.rrc_filter`](@ref)) the matched taps equal `taps`.

An impulse at index 1 produces those matched taps and nothing else. On a symmetric pulse the largest sample of that impulse response is the centre tap, 1-based index `(length(taps) + 1) ÷ 2`, not the last sample.

The deterministic autocorrelation is this filter applied to the taps themselves. Its zero-lag sample is 1-based index `length(taps)`,

```math
\bigl(g * \tilde g\bigr)[\texttt{length}(g)] = \sum_n |g_n|^2,
```

with ``\tilde g = \operatorname{conj}(\operatorname{reverse}(g))``. [`downsample_symbols`](@ref) does not move to that index on its own.

# Cascade with `pulse_shape`
[`CommToolkit.Visualization.pulse_shape`](@ref) places symbol ``k`` (1-based) at upsampled index ``(k-1)\,\texttt{sps}+1``, convolves with `taps`, and truncates to `length(symbols) * sps` samples. The sampling phase that hits the zero-lag autocorrelation of symbol 1 is

```math
\texttt{offset} = \texttt{length(taps)} - 1,
```

for any `taps`, provided `length(symbols) * sps >= length(taps)` so the transmit truncation has not already cut that first pulse off. `offset` counts samples skipped before the first kept sample; see [`downsample_symbols`](@ref). The sample equals the symbol only when the autocorrelation is Nyquist and scaled to one (a unit-energy root-raised-cosine pair is approximately that; a single raised-cosine filter is already Nyquist and should not be matched again). Symbols within about half a filter span of the **end** of the burst are missing the tail `pulse_shape` discarded, so they are not the full matched-filter response. This function does not search for the peak.
"""
function matched_filter(samples::AbstractVector{<:Number}, taps::AbstractVector{<:Number})
    return Channel.apply_fir(samples, conj.(reverse(taps)))
end

@doc raw"""
    downsample_symbols(samples, sps::Integer; offset::Integer=0) -> Vector

Keep every `sps`-th sample. The first kept sample is 1-based index `offset + 1`, then `offset + 1 + sps`, and so on while the index lies inside `samples`.

# Arguments
- `sps`: samples per symbol. Must be ≥ 1. The function does not check that `samples` was built with [`CommToolkit.Visualization.pulse_shape`](@ref) at this `sps`, and `length(samples)` need not be a multiple of `sps`.
- `offset`: samples skipped before the first kept sample. This is the sampling phase. It must be ≥ 0 and strictly less than `length(samples)`. It is **not** reduced modulo `sps`: `offset ≥ sps` also drops whole leading symbols. The default `0` keeps indices `1, 1+sps, 1+2sps, ...`, which is the autocorrelation peak only when that peak is already at index 1.

There is no timing-error detector and no interpolation between samples. The result is a newly allocated vector of the same element type as `samples`.
"""
function downsample_symbols(samples::AbstractVector{<:Number}, sps::Integer; offset::Integer=0)
    sps = Int(sps)
    offset = Int(offset)
    sps >= 1 || throw(ArgumentError("sps must be ≥ 1, got $sps"))
    offset >= 0 || throw(ArgumentError("offset must be ≥ 0, got $offset"))
    n = length(samples)
    offset < n || throw(ArgumentError(
        "offset $offset is past the end of a length-$n vector",
    ))
    idx = (offset + 1):sps:n
    out = Vector{eltype(samples)}(undef, length(idx))
    @inbounds for (k, i) in enumerate(idx)
        out[k] = samples[i]
    end
    return out
end

end # module Pulse
