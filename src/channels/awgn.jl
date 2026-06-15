# AWGN channel model
#
# Adds complex or real white Gaussian noise
# SNR defined as Eb/N0 (energy per bit to noise power spectral density)

function awgn(syms::AbstractVector{T}, snr_db::Real; energy_per_bit::Real=1.0) where T
    snr_lin = snr_db_to_linear(snr_db)
    noise_var = energy_per_bit / (2 * snr_lin)
    noise_std = sqrt(noise_var)

    if T <: Complex
        noise = noise_std * (randn(length(syms)) + im * randn(length(syms)))
    else
        noise = noise_std * randn(length(syms))
    end
    return syms .+ noise
end

function awgn(syms::AbstractVector{T}, snr_db::Real, bits_per_sym::Int) where T
    energy_per_bit = 1.0 / bits_per_sym
    return awgn(syms, snr_db; energy_per_bit = energy_per_bit)
end
