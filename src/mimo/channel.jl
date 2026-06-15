# MIMO channel: i.i.d. Rayleigh flat fading
#
# y = H * x + n
# H: (n_rx x n_tx) complex Gaussian matrix, each entry ~ CN(0,1)

function mimo_channel(n_tx::Int, n_rx::Int)
    (randn(n_rx, n_tx) + im * randn(n_rx, n_tx)) / sqrt(2)
end

# Apply MIMO channel: y = H*x + n
# noise_var is the variance per real dimension
function apply_mimo(x::Vector{ComplexF64}, H::Matrix{ComplexF64}, noise_var::Real)
    y = H * x
    noise = sqrt(noise_var / 2) * (randn(length(y)) + im * randn(length(y)))
    y + noise
end
