# QAM modulation
#
# Generic M-QAM with Gray coding
# Currently supports: 16QAM (m=16)

function _qam_constellation(m::Int)
    @assert m >= 4 && ispow2(m) "M must be a power of 2, got $m"
    n = Int(sqrt(m))
    @assert n^2 == m "M must be a perfect square for square QAM"

    # Generate n-PAM on each axis: -(n-1):2:n-1, then normalize
    pam = [-(n-1):2:(n-1)...] ./ sqrt(2 * (m - 1) / 3)

    # Gray coding for PAM levels
    gray = zeros(Int, n)
    for i in 1:n
        gray[i] = xor(i - 1, (i - 1) >> 1)
    end

    # Build constellation
    constell = Vector{ComplexF64}(undef, m)
    for i in 1:n, j in 1:n
        idx = gray[i] * n + gray[j] + 1
        constell[idx] = pam[i] + pam[j] * im
    end
    return constell
end

function modulate(bits::AbstractVector{<:Integer}, scheme::QAM)
    @assert all(b -> b ∈ (0, 1), bits) "bits must be 0 or 1"
    m = scheme.m
    k = Int(log2(m))
    @assert length(bits) % k == 0 "bit length must be multiple of log2(M)"

    constell = _qam_constellation(m)
    n_sym = length(bits) ÷ k
    syms = Vector{ComplexF64}(undef, n_sym)

    for i in 1:n_sym
        idx = 1
        for j in 1:k
            idx += bits[(i-1)*k + j] << (k - j)
        end
        syms[i] = constell[idx]
    end
    return syms
end

function demodulate(syms::AbstractVector{<:Complex}, scheme::QAM)
    m = scheme.m
    k = Int(log2(m))
    n_sym = length(syms)
    bits = Vector{Int}(undef, n_sym * k)

    constell = _qam_constellation(m)

    for i in 1:n_sym
        # Minimum distance detection
        _, min_idx = findmin(abs2.(syms[i] .- constell))
        # Convert to bits (reverse Gray)
        for j in 1:k
            bits[(i-1)*k + j] = ((min_idx - 1) >> (k - j)) & 1
        end
    end
    return bits
end
