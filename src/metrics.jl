module Metrics

export ber, ser, fer, evm, snr_db_to_linear

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

@doc raw"""
    evm(reference, received; unit=:percent) -> Float64

RMS error-vector magnitude of `received` against `reference`.

```math
\mathrm{EVM}_{\mathrm{rms}}
  = \sqrt{ \frac{\sum_i |r_i - s_i|^2}{\sum_i |s_i|^2} }.
```

`reference` holds the ideal symbols ``s_i`` and `received` the observed symbols ``r_i``. One entry is one symbol. The two vectors must have the same nonzero length.

# Units
- `unit = :percent` (default): ``100 \cdot \mathrm{EVM}_{\mathrm{rms}}``. This is a percentage of the RMS reference amplitude, not a percentage of the reference power. An RMS ratio of ``0.1`` is `10.0`, i.e. 10%.
- `unit = :db`: ``20 \log_{10}(\mathrm{EVM}_{\mathrm{rms}})``, which is the same number as ``10 \log_{10}\) of the error power over the reference power. 10% is ``-20``. A perfect match returns `-Inf`.
- `unit = :linear`: the dimensionless ratio ``\mathrm{EVM}_{\mathrm{rms}}``.

# Notes
The sums are over the vectors as passed in. There is no constellation template, DC removal, gain fit, or phase rotation: a common phase offset is an error. Scaling both vectors by the same nonzero constant does not change the result. The ratio is undefined when the reference energy is zero, and when any symbol is non-finite.
"""
function evm(reference::AbstractVector{<:Number}, received::AbstractVector{<:Number};
             unit::Symbol=:percent)
    n = length(reference)
    length(received) == n || throw(ArgumentError(
        "symbol vectors must have the same length, got $n and $(length(received))",
    ))
    n > 0 || throw(ArgumentError("symbol vectors must be non-empty"))
    unit === :percent || unit === :db || unit === :linear || throw(ArgumentError(
        "unit must be :percent, :db, or :linear, got $(repr(unit))",
    ))
    err = 0.0
    refp = 0.0
    for (s, r) in zip(reference, received)
        es = abs2(float(s))
        er = abs2(float(r) - float(s))
        if !isfinite(es) || !isfinite(er)
            throw(ArgumentError("reference and received symbols must be finite"))
        end
        refp += es
        err += er
    end
    refp > 0 || throw(ArgumentError("reference energy is zero; EVM is undefined"))
    ratio = sqrt(err / refp)
    if unit === :linear
        return ratio
    elseif unit === :percent
        return 100 * ratio
    else
        return 20 * log10(ratio)
    end
end

end # module Metrics
