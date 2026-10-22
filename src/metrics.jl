module Metrics

export ber, ser, fer, snr_db_to_linear

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

@doc raw"""
    ser(tx, rx) -> Float64

Empirical symbol error rate ``k/n``.

# Arguments
- `tx`, `rx`: vectors of equal, nonzero length, compared elementwise with `!=`. One entry is one symbol (a constellation point, a symbol index, or any other label).

# Returns
The fraction of positions at which the two vectors differ.

Noisy samples are not sliced, and bits are not grouped into multi-bit symbols. On integer bit vectors the value equals [`ber`](@ref).
"""
function ser(tx::AbstractVector, rx::AbstractVector)
    n = length(tx)
    length(rx) == n || throw(ArgumentError(
        "symbol vectors must have the same length, got $n and $(length(rx))",
    ))
    n > 0 || throw(ArgumentError("symbol vectors must be non-empty"))
    k = 0
    for (a, b) in zip(tx, rx)
        k += a != b
    end
    return k / n
end

@doc raw"""
    fer(tx, rx, block_length::Integer) -> Float64

Empirical frame error rate. `tx` and `rx` are partitioned into contiguous blocks of `block_length`. A block is in error when any of its entries differs.

# Arguments
- `tx`, `rx`: vectors of equal, nonzero length. The length must be a multiple of `block_length`.
- `block_length`: positive number of entries per block.

# Returns
The fraction of blocks that contain at least one mismatch under `!=`, the same comparison as [`ber`](@ref) and [`ser`](@ref).
"""
function fer(tx::AbstractVector, rx::AbstractVector, block_length::Integer)
    n = length(tx)
    length(rx) == n || throw(ArgumentError(
        "frame vectors must have the same length, got $n and $(length(rx))",
    ))
    n > 0 || throw(ArgumentError("frame vectors must be non-empty"))
    bl = Int(block_length)
    bl >= 1 || throw(ArgumentError("block_length must be ≥ 1, got $bl"))
    n % bl == 0 || throw(ArgumentError(
        "length $n is not a multiple of block_length $bl",
    ))
    bad = 0
    run = 0
    frame_bad = false
    for (a, b) in zip(tx, rx)
        frame_bad |= a != b
        run += 1
        if run == bl
            bad += frame_bad
            frame_bad = false
            run = 0
        end
    end
    return bad / (n ÷ bl)
end

end # module Metrics
