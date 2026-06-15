module CommToolkit

export
    # Modulation types
    BPSK, QPSK, QAM,
    # Modulation functions
    modulate, demodulate,
    # Channel models
    awgn, rayleigh_channel,
    # Coding
    conv_encode, viterbi_decode,
    # Equalizers
    zf_equalizer, mmse_equalizer, zf_equalize, mmse_equalize,
    # OFDM
    OFDMParams, TDLChannel,
    generate_ofdm_symbol, ofdm_tx, apply_multipath, apply_channel,
    extract_ofdm_symbols, ofdm_rx_known, ofdm_rx_ls, ofdm_rx_mmse,
    interpolate_channel,
    # OFDM Sync
    sto_estimate, cfo_estimate, cfo_correct, ofdm_sync,
    # Metrics
    ber,
    # Utilities
    snr_db_to_linear

# Modulation types
abstract type ModulationScheme end
struct BPSK <: ModulationScheme end
struct QPSK <: ModulationScheme end
struct QAM <: ModulationScheme
    m::Int  # constellation size, must be power of 4
end

# Utility functions
snr_db_to_linear(snr_db::Real) = 10^(snr_db / 10)

# Channel models
include("channels/awgn.jl")
include("channels/fading.jl")

# Modems
include("modems/bpsk.jl")
include("modems/qpsk.jl")
include("modems/qam.jl")

# Coding
include("coding/conv.jl")

# Equalizers
include("equalizers/linear.jl")

# OFDM
include("ofdm/params.jl")
include("ofdm/transmitter.jl")
include("ofdm/channel.jl")
include("ofdm/receiver.jl")
include("ofdm/sync.jl")

# Metrics
include("utils/metrics.jl")

end # module
