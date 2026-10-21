module Metrics

export ber, snr_db_to_linear

@doc raw"""
    snr_db_to_linear(snr_db::Real) -> Float64

Convert a ratio expressed in decibels to a linear ratio,

```math
10^{\mathrm{snr\_db}/10}.
```

`snr_db` may be `+Inf`, which returns `Inf`.
"""
function snr_db_to_linear(snr_db::Real)
    return exp10(float(snr_db) / 10)
end

@doc raw"""
    ber(tx_bits::AbstractVector{<:Integer}, rx_bits::AbstractVector{<:Integer}) -> Float64

Empirical bit error rate ``k/n``.

# Arguments
- `tx_bits`, `rx_bits`: integer vectors of equal, nonzero length. Compared elementwise; no type of coding is assumed.

# Returns
The fraction of positions at which the two vectors differ.
"""
function ber(tx_bits::AbstractVector{<:Integer}, rx_bits::AbstractVector{<:Integer})
    n = length(tx_bits)
    length(rx_bits) == n || throw(ArgumentError(
        "bit vectors must have the same length, got $n and $(length(rx_bits))",
    ))
    n > 0 || throw(ArgumentError("bit vectors must be non-empty"))
    k = 0
    @inbounds for i in eachindex(tx_bits)
        k += tx_bits[i] != rx_bits[i]
    end
    return k / n
end

end # module Metrics
