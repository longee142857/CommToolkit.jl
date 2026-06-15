# Rayleigh fading channel models
#
# Flat Rayleigh fading: h ~ CN(0, 1), i.i.d. per symbol
# The received signal is: y = h * x + n

function rayleigh_channel(syms::AbstractVector, snr_db::Real, bits_per_sym::Int)
    n = length(syms)
    # Convert to complex if real
    syms_c = ComplexF64.(syms)

    # Generate fading coefficients (unit power)
    h = (randn(n) + im * randn(n)) / sqrt(2)

    # Apply fading
    faded = syms_c .* h

    # Add AWGN
    energy_per_bit = 1.0 / bits_per_sym
    snr_lin = snr_db_to_linear(snr_db)
    noise_var = energy_per_bit / (2 * snr_lin)
    noise_std = sqrt(noise_var)
    noise = noise_std * (randn(n) + im * randn(n))

    return faded .+ noise, h
end

# Known-CSIR version: return equalized symbols
function rayleigh_channel(syms::AbstractVector, snr_db::Real, bits_per_sym::Int, perfect_csi::Bool)
    rx, h = rayleigh_channel(syms, snr_db, bits_per_sym)
    if perfect_csi
        return rx ./ h
    end
    return rx
end
