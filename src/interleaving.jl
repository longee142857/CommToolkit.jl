"""
    Interleaving

Rectangular block interleaver: write by rows, read by columns. The same permutation applies to bits, symbols, and LLRs.
"""
module Interleaving

export block_interleave, block_deinterleave

function _block_dims(n_rows::Integer, n_cols::Integer, n::Int)
    nr = Int(n_rows)
    nc = Int(n_cols)
    nr >= 1 || throw(ArgumentError("n_rows must be ≥ 1, got $nr"))
    nc >= 1 || throw(ArgumentError("n_cols must be ≥ 1, got $nc"))
    nr * nc == n || throw(ArgumentError(
        "length $n is not n_rows * n_cols = $(nr * nc)",
    ))
    return nr, nc
end

function _checked_rows(n_rows::Integer, n::Int)
    nr = Int(n_rows)
    nr >= 1 || throw(ArgumentError("n_rows must be ≥ 1, got $nr"))
    n % nr == 0 || throw(ArgumentError("length $n is not divisible by n_rows $nr"))
    return nr
end

@doc raw"""
    block_interleave(x::AbstractVector, n_rows::Integer, n_cols::Integer) -> Vector
    block_interleave(x::AbstractVector, n_rows::Integer) -> Vector

Write `x` into an `n_rows` × `n_cols` matrix by rows and read it by columns.

The two-argument method uses `n_cols = length(x) ÷ n_rows`. `length(x)` must equal `n_rows * n_cols`, and both dimensions must be at least 1. Linear index `i` (from the start of `x`) at matrix position `(r, c)`, with `r` and `c` 1-based, moves from `(r - 1) n_{\mathrm{cols}} + c` to `(c - 1) n_{\mathrm{rows}} + r`.

[`block_deinterleave`](@ref) is the inverse for the same dimensions. The alphabet is not inspected: a bit vector, a symbol vector, and an LLR vector of equal length take the same permutation.
"""
function block_interleave(x::AbstractVector, n_rows::Integer, n_cols::Integer)
    nr, nc = _block_dims(n_rows, n_cols, length(x))
    out = Vector{eltype(x)}(undef, nr * nc)
    base = firstindex(x)
    @inbounds for c in 1:nc
        for r in 1:nr
            out[(c - 1) * nr + r] = x[base + (r - 1) * nc + (c - 1)]
        end
    end
    return out
end

function block_interleave(x::AbstractVector, n_rows::Integer)
    nr = _checked_rows(n_rows, length(x))
    return block_interleave(x, nr, length(x) ÷ nr)
end

@doc raw"""
    block_deinterleave(x::AbstractVector, n_rows::Integer, n_cols::Integer) -> Vector
    block_deinterleave(x::AbstractVector, n_rows::Integer) -> Vector

Inverse of [`block_interleave`](@ref) at the same `n_rows` and `n_cols`.

The two-argument method divides `length(x)` by the original `n_rows`. Pass the row count that was used to interleave, not the swapped column count of the read-out matrix.
"""
function block_deinterleave(x::AbstractVector, n_rows::Integer, n_cols::Integer)
    nr, nc = _block_dims(n_rows, n_cols, length(x))
    return block_interleave(x, nc, nr)
end

function block_deinterleave(x::AbstractVector, n_rows::Integer)
    nr = _checked_rows(n_rows, length(x))
    return block_interleave(x, length(x) ÷ nr, nr)
end

end # module Interleaving
