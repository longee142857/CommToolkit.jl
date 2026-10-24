"""
    CommToolkit

Baseband communication-link simulations: Gray PSK/QAM (including ``M``-PSK) with hard decisions and AWGN bit LLRs, DBPSK/DQPSK differential detection, AWGN, symbol-wise flat Rayleigh and Rice, block-constant Rayleigh, flat Clarke/Jakes time correlation, i.i.d. Nakagami-m, ideal-CSI receive MRC, a deterministic single-carrier FIR, a terminated ``K = 3`` convolutional code with hard-decision Viterbi, soft-input Viterbi, BCJR, and SOVA, puncturing of that mother code to nominal rates 2/3 and 3/4, a short parallel-concatenated turbo code on the recursive systematic form of those generators, syndrome decoding of the Hamming ``[7,4]`` code, a binary BCH (15, 7) decoder over GF(16), a rectangular block interleaver, short regular LDPC codes with systematic encoding and flooding sum-product / min-sum decoding, a teaching CRC / MCS / Chase-HARQ loop, linear ZF/MMSE equalisation, trained symbol-spaced LMS, a finite-block ZF/MMSE decision-feedback equalizer, a textbook CP-OFDM chain, small MIMO detectors, SVD water-filling, and Alamouti coding, an additive LFSR scrambler, raised-cosine pulses, a matched filter, symbol-rate downsampling, teaching-grade Gardner, Costas, and Schmidl–Cox synchronisers, and a teaching BICM-ID receiver.

The public names are re-exported from the submodules `Modulation`, `Channel`, `Coding`, `Interleaving`, `BICM`, `Link`, `Equalization`, `OFDM`, `MIMO`, `Visualization`, `Pulse`, `Scrambling`, `Sync`, and `Metrics`. Those exported names are the stability surface. `modulation` is imported from `OFDM` and is intentionally not exported; call it as `CommToolkit.modulation`. Calling `using CommToolkit` does not load Plots; constellation, eye, and pulse figures are provided by the package extension `CommToolkitPlotsExt` after `using Plots`.

This package is a simulation library for checking algorithms against textbook formulas. It does not implement 3GPP NR or LTE PHY, an O-RAN fronthaul, or a geometric 3GPP TR 38.901 clustered delay line, and it is not a production physical-layer stack. The LDPC codes are short teaching constructions, not NR base graphs. The turbo code is a short teaching construction: two memory-2 recursive systematic encoders, a rectangular or stride interleaver, at most 1024 information bits, and at most 32 log-BCJR iterations. It is not a 3GPP turbo code. The BCH decoder is only the binary length-15, dimension-7, designed-distance-5 code, not a Reed–Solomon library and not a general BCH constructor. Hamming hard decisions are syndrome decoding of the ``[7,4]`` code; the LDPC fixture on the same parity-check matrix is not that decoder. The MCS table and HARQ loop are a link-level teaching reference, not an NR scheduler. Convolutional puncturing is the classical periodic mask on the ``K = 3`` mother code, not NR rate matching. Block fading holds one coefficient fixed inside a block. `jakes_fading` is a separate flat Clarke spectrum, synthesised by an inverse DFT of integrated bins, not a classical Jakes sum of sinusoids and not a frequency-selective tap model. Nakagami-m coefficients are i.i.d. The FIR is a tap vector supplied by the caller, not a fading process. The matched filter does not estimate timing or carrier phase. `Sync` adds a Gardner detector, a Costas loop for BPSK and QPSK, and a Schmidl–Cox preamble estimator; those are teaching loops with a stated capture range, not a production synchroniser, and they do not change `ofdm_sync`. LMS adapts on a training sequence and then freezes; it is not blind and it is not RLS. The block DFE feeds back sliced decisions on the same convolution matrix as the linear equalizers, and a wrong decision can propagate. `evm` is the RMS error-vector magnitude of a received block against the caller's reference symbols, in percent or decibels, with no gain or phase fit. The scrambler XORs a block with an LFSR sequence and restarts from the same seed to undo it; it is not a self-synchronizing scrambler, not a Gold code, and it is not wired into the MCS loop. `svd_waterfill` puts a single total-power water-fill on the singular modes; `svd_precoder` still returns the unweighted right singular vectors. SOVA, the turbo decoder, and the BICM-ID loop are not wired into the MCS table. `BICM` interleaves coded bits, maps them onto the existing constellations, and can feed coded-bit extrinsic LLRs from a convolutional log-MAP decoder or a short LDPC decoder back into the demapper. It is not 5G NR BICM and it has no EXIT chart. `:sova` and `:turbo` are a single demapper pass in that loop: those decoders return information-bit LLRs, not coded-bit extrinsic values. Gray labelling is the default; an optional 8-PSK set-partitioning map is what makes a short-block iterative gain visible.
"""
module CommToolkit

include("metrics.jl")
include("scrambling.jl")
include("modulation.jl")
include("channel.jl")
include("pulse.jl")
include("coding.jl")
include("interleaving.jl")
include("bicm.jl")
include("link.jl")
include("equalization.jl")
include("ofdm.jl")
include("mimo.jl")
include("visualization.jl")
include("sync.jl")

using .Metrics: ber, ser, fer, evm, snr_db_to_linear
using .Scrambling: prbs_polynomial, lfsr_sequence, scramble, descramble
using .Modulation: ModulationScheme, BPSK, QPSK, PSK, DBPSK, DQPSK, QAM, modulate, demodulate,
    demodulate_soft, constellation, bits_per_symbol
using .Channel: awgn, rayleigh_channel, block_rayleigh, block_rayleigh_channel,
    jakes_fading, jakes_channel,
    rice_fading, rice_channel,
    nakagami_envelope, nakagami_fading, nakagami_channel,
    mrc, mrc_channel,
    apply_fir, fir_channel
using .Pulse: matched_filter, downsample_symbols
using .Coding: conv_encode, viterbi_decode, bcjr_decode, sova_decode, distance_spectrum,
    puncture_mask, puncture, depuncture,
    turbo_interleaver, turbo_encode, turbo_decode,
    GF2m, gf_add, gf_mul, gf_div, gf_inv, gf_pow,
    hamming74_encode, hamming74_decode, bch157_encode, bch157_decode,
    LDPCCode, LDPCResult,
    gallager_ldpc, array_ldpc, hamming74_ldpc, ldpc_from_H,
    ldpc_encode, ldpc_decode, ldpc_decode_details, ldpc_syndrome,
    ldpc_generator_matrix, ldpc_column_degrees, ldpc_row_degrees
using .Interleaving: block_interleave, block_deinterleave
using .BICM: BICMResult, bicm_labels, bicm_psk8_id_labels, bicm_map, bicm_demap, bicm_decode
using .Link: CRC16_CCITT_POLY, crc16_ccitt, crc_append, crc_check,
    MCS, mcs_convolutional, mcs_ldpc, teaching_mcs_table, spectral_efficiency,
    ThresholdAMC, OuterLoopAMC, select_mcs, observe_acknack!, teaching_amc_thresholds,
    estimate_esn0_db, LinkReport, simulate_link
using .Equalization: zf_equalizer, zf_equalize, mmse_equalizer, mmse_equalize,
    lms_equalize, DFECoefficients, dfe_equalizer, dfe_equalize, zf_dfe, mmse_dfe
using .OFDM: OFDMParams, TDLChannel, modulation,
    generate_ofdm_symbol, ofdm_tx,
    apply_multipath, apply_channel,
    extract_ofdm_symbols, ofdm_rx_known, ofdm_rx_ls, ofdm_rx_mmse,
    interpolate_channel,
    sto_estimate, cfo_estimate, cfo_correct, ofdm_sync
using .MIMO: mimo_channel, apply_mimo,
    ml_detect, zf_detect, mmse_detect, sic_detect,
    zf_precoder, svd_precoder, waterfill, svd_waterfill,
    alamouti_encode, alamouti_decode_2x1, alamouti_decode_2x2
using .Visualization: raised_cosine, rrc_filter, pulse_shape,
    plot_constellation, plot_eye_diagram, plot_pulse_response
using .Sync: fractional_delay,
    gardner_ted, gardner_loop, gardner_open_loop, single_carrier_timing,
    costas_detector, costas_loop,
    schmidl_preamble, schmidl_cox, integer_cfo_estimate, schmidl_sync

export
    ModulationScheme, BPSK, QPSK, PSK, DBPSK, DQPSK, QAM,
    modulate, demodulate, demodulate_soft, constellation, bits_per_symbol,
    awgn, rayleigh_channel, block_rayleigh, block_rayleigh_channel,
    jakes_fading, jakes_channel,
    rice_fading, rice_channel,
    nakagami_envelope, nakagami_fading, nakagami_channel,
    mrc, mrc_channel,
    apply_fir, fir_channel,
    matched_filter, downsample_symbols,
    conv_encode, viterbi_decode, bcjr_decode, sova_decode, distance_spectrum,
    puncture_mask, puncture, depuncture,
    turbo_interleaver, turbo_encode, turbo_decode,
    GF2m, gf_add, gf_mul, gf_div, gf_inv, gf_pow,
    hamming74_encode, hamming74_decode, bch157_encode, bch157_decode,
    block_interleave, block_deinterleave,
    BICMResult, bicm_labels, bicm_psk8_id_labels, bicm_map, bicm_demap, bicm_decode,
    LDPCCode, LDPCResult,
    gallager_ldpc, array_ldpc, hamming74_ldpc, ldpc_from_H,
    ldpc_encode, ldpc_decode, ldpc_decode_details, ldpc_syndrome,
    ldpc_generator_matrix, ldpc_column_degrees, ldpc_row_degrees,
    CRC16_CCITT_POLY, crc16_ccitt, crc_append, crc_check,
    MCS, mcs_convolutional, mcs_ldpc, teaching_mcs_table, spectral_efficiency,
    ThresholdAMC, OuterLoopAMC, select_mcs, observe_acknack!, teaching_amc_thresholds,
    estimate_esn0_db, LinkReport, simulate_link,
    zf_equalizer, mmse_equalizer, zf_equalize, mmse_equalize,
    lms_equalize, DFECoefficients, dfe_equalizer, dfe_equalize, zf_dfe, mmse_dfe,
    OFDMParams, TDLChannel,
    generate_ofdm_symbol, ofdm_tx, apply_multipath, apply_channel,
    extract_ofdm_symbols, ofdm_rx_known, ofdm_rx_ls, ofdm_rx_mmse,
    interpolate_channel,
    sto_estimate, cfo_estimate, cfo_correct, ofdm_sync,
    mimo_channel, apply_mimo,
    ml_detect, zf_detect, mmse_detect, sic_detect,
    zf_precoder, svd_precoder, waterfill, svd_waterfill,
    alamouti_encode, alamouti_decode_2x1, alamouti_decode_2x2,
    raised_cosine, rrc_filter, pulse_shape,
    plot_constellation, plot_eye_diagram, plot_pulse_response,
    prbs_polynomial, lfsr_sequence, scramble, descramble,
    fractional_delay,
    gardner_ted, gardner_loop, gardner_open_loop, single_carrier_timing,
    costas_detector, costas_loop,
    schmidl_preamble, schmidl_cox, integer_cfo_estimate, schmidl_sync,
    ber, ser, fer, evm, snr_db_to_linear

end # module CommToolkit
