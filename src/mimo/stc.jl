# Alamouti space-time block code (2x1 and 2x2)
#
# Tx1: [x1, -conj(x2)]
# Tx2: [x2,  conj(x1)]

# Alamouti encoder: two symbols → two time slots
function alamouti_encode(x1::ComplexF64, x2::ComplexF64)
    [x1  -conj(x2);
     x2   conj(x1)]
end

# Alamouti 2x1 decoder (single receive antenna)
# h1, h2: channel coefficients from Tx1, Tx2 to Rx
# r1, r2: received signals at time 1, 2
function alamouti_decode_2x1(r1::ComplexF64, r2::ComplexF64,
                              h1::ComplexF64, h2::ComplexF64)
    denom = abs2(h1) + abs2(h2)
    if denom == 0
        return 0.0 + 0.0im, 0.0 + 0.0im
    end
    x1_est = (conj(h1) * r1 + h2 * conj(r2)) / denom
    x2_est = (conj(h2) * r1 - h1 * conj(r2)) / denom
    x1_est, x2_est
end

# Alamouti 2x2 decoder (two receive antennas)
function alamouti_decode_2x2(r1::Vector{ComplexF64}, r2::Vector{ComplexF64},
                              H::Matrix{ComplexF64})
    h11, h21 = H[1, 1], H[2, 1]
    h12, h22 = H[1, 2], H[2, 2]
    r11, r21 = r1[1], r1[2]
    r12, r22 = r2[1], r2[2]

    denom = abs2(h11) + abs2(h21) + abs2(h12) + abs2(h22)
    if denom == 0
        return 0.0 + 0.0im, 0.0 + 0.0im
    end

    x1_est = (conj(h11)*r11 + h12*conj(r12) + conj(h21)*r21 + h22*conj(r22)) / denom
    x2_est = (conj(h12)*r11 - h11*conj(r12) + conj(h22)*r21 - h21*conj(r22)) / denom
    x1_est, x2_est
end
