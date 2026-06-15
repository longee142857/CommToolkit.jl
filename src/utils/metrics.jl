# BER calculation

function ber(tx_bits::AbstractVector{<:Integer}, rx_bits::AbstractVector{<:Integer})
    @assert length(tx_bits) == length(rx_bits) "bit vectors must have same length"
    n_errors = sum(tx_bits .!= rx_bits)
    return n_errors / length(tx_bits)
end
