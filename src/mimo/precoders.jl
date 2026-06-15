# MIMO precoders: ZF and SVD-based

using LinearAlgebra

# ZF precoder: P = H' * inv(H * H')
# Assumes n_tx >= n_rx
function zf_precoder(H::Matrix{ComplexF64})
    H' * inv(H * H')
end

# SVD-based precoder: water-filling inspired
# P = V, treating H = U * S * V'
# Returns precoder P and the diagonal gains S
function svd_precoder(H::Matrix{ComplexF64})
    U, S, Vt = svd(H)
    V = Vt'
    V, S
end
