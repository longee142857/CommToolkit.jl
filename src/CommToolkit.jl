"""
    CommToolkit

Baseband communication-link simulations: Gray PSK/QAM, AWGN and flat Rayleigh channels, a terminated ``K = 3`` convolutional code, linear ZF/MMSE equalisation, a textbook CP-OFDM chain, small MIMO detectors and Alamouti coding, and raised-cosine pulses.

The public names are re-exported from the submodules `Modulation`, `Channel`, `Coding`, `Equalization`, `OFDM`, `MIMO`, `Visualization`, and `Metrics`. Calling `using CommToolkit` does not load Plots; constellation, eye, and pulse figures are provided by the package extension `CommToolkitPlotsExt` after `using Plots`.

This package is a simulation library for checking algorithms against textbook formulas. It does not implement 3GPP NR, LTE, or an O-RAN fronthaul, and it is not a production physical-layer stack.
"""
module CommToolkit

include("metrics.jl")
include("modulation.jl")
include("channel.jl")
include("coding.jl")
include("equalization.jl")
include("ofdm.jl")
include("mimo.jl")
include("visualization.jl")

using .Metrics: ber, snr_db_to_linear
using .Modulation: ModulationScheme, BPSK, QPSK, QAM, modulate, demodulate,
    demodulate_soft, constellation, bits_per_symbol
using .Channel: awgn, rayleigh_channel
using .Coding: conv_encode, viterbi_decode, distance_spectrum
using .Equalization: zf_equalizer, zf_equalize, mmse_equalizer, mmse_equalize
using .OFDM: OFDMParams, TDLChannel, modulation,
    generate_ofdm_symbol, ofdm_tx,
    apply_multipath, apply_channel,
    extract_ofdm_symbols, ofdm_rx_known, ofdm_rx_ls, ofdm_rx_mmse,
    interpolate_channel,
    sto_estimate, cfo_estimate, cfo_correct, ofdm_sync
using .MIMO: mimo_channel, apply_mimo,
    ml_detect, zf_detect, mmse_detect, sic_detect,
    zf_precoder, svd_precoder,
    alamouti_encode, alamouti_decode_2x1, alamouti_decode_2x2
using .Visualization: raised_cosine, rrc_filter, pulse_shape,
    plot_constellation, plot_eye_diagram, plot_pulse_response

export
    ModulationScheme, BPSK, QPSK, QAM,
    modulate, demodulate, demodulate_soft, constellation, bits_per_symbol,
    awgn, rayleigh_channel,
    conv_encode, viterbi_decode, distance_spectrum,
    zf_equalizer, mmse_equalizer, zf_equalize, mmse_equalize,
    OFDMParams, TDLChannel,
    generate_ofdm_symbol, ofdm_tx, apply_multipath, apply_channel,
    extract_ofdm_symbols, ofdm_rx_known, ofdm_rx_ls, ofdm_rx_mmse,
    interpolate_channel,
    sto_estimate, cfo_estimate, cfo_correct, ofdm_sync,
    mimo_channel, apply_mimo,
    ml_detect, zf_detect, mmse_detect, sic_detect,
    zf_precoder, svd_precoder,
    alamouti_encode, alamouti_decode_2x1, alamouti_decode_2x2,
    raised_cosine, rrc_filter, pulse_shape,
    plot_constellation, plot_eye_diagram, plot_pulse_response,
    ber, snr_db_to_linear

end # module CommToolkit
