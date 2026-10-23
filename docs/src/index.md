# CommToolkit.jl API

API reference for the exported interface. This page is not a tutorial.

The package simulates textbook baseband links so that Monte Carlo error rates can be checked against closed-form formulas. It does not implement 3GPP NR or LTE, it does not implement an O-RAN fronthaul, and it is not a production physical-layer stack. LDPC support is a short teaching reference (Gallager and array constructions), not an NR base-graph implementation. The MCS table, CRC, and Chase combining loop are a link-level teaching reference, not an NR scheduler or HARQ entity. `PSK` is coherent Gray ``M``-PSK. `DBPSK` and `DQPSK` are two-symbol differential detectors, not ``\pi/4``-DQPSK. `rice_channel` takes a linear Rice factor. `mrc` is ideal-CSI receive combining, not a MIMO detector.

Names are defined in the submodules `Modulation`, `Channel`, `Coding`, `Interleaving`, `Link`, `Equalization`, `OFDM`, `MIMO`, `Visualization`, `Pulse`, `Scrambling`, and `Metrics`, and re-exported here. Plotting functions need Plots.jl; `using CommToolkit` does not load it.

```@docs
CommToolkit
```

## Modulation

```@docs
ModulationScheme
BPSK
QPSK
PSK
DBPSK
DQPSK
QAM
bits_per_symbol
constellation
modulate
demodulate
demodulate_soft
```

## Channel and metrics

```@docs
snr_db_to_linear
awgn
rayleigh_channel
block_rayleigh
block_rayleigh_channel
rice_fading
rice_channel
mrc
mrc_channel
apply_fir
fir_channel
ber
ser
fer
evm
```

## Coding

```@docs
conv_encode
viterbi_decode
bcjr_decode
distance_spectrum
puncture_mask
puncture
depuncture
LDPCCode
LDPCResult
gallager_ldpc
array_ldpc
hamming74_ldpc
ldpc_from_H
ldpc_encode
ldpc_decode
ldpc_decode_details
ldpc_syndrome
ldpc_generator_matrix
ldpc_column_degrees
ldpc_row_degrees
```

## Scrambling

Additive LFSR. `descramble` restarts the same polynomial and seed; it is not a self-synchronizing scrambler, not a Gold code, and not an NR scrambler. Named degrees are connection polynomials. The sequence is the raw shift-register output, not an inverted ITU-T O.150 test pattern.

```@docs
CommToolkit.Scrambling
prbs_polynomial
lfsr_sequence
scramble
descramble
```

## Interleaving

Rectangular block interleaver. Not wired into the MCS table.

```@docs
CommToolkit.Interleaving
block_interleave
block_deinterleave
```

## Link

Teaching CRC, MCS table, Chase combining, and threshold or outer-loop MCS selection. Not an NR protocol.

```@docs
CommToolkit.Link
CRC16_CCITT_POLY
crc16_ccitt
crc_append
crc_check
MCS
spectral_efficiency
mcs_convolutional
mcs_ldpc
teaching_mcs_table
ThresholdAMC
OuterLoopAMC
select_mcs
observe_acknack!
teaching_amc_thresholds
estimate_esn0_db
LinkReport
simulate_link
```

## Equalization

```@docs
zf_equalizer
zf_equalize
mmse_equalizer
mmse_equalize
lms_equalize
DFECoefficients
dfe_equalizer
dfe_equalize
zf_dfe
mmse_dfe
```

## OFDM

`modulation` is available as `CommToolkit.modulation` and is not part of the exported name list.

```@docs
OFDMParams
TDLChannel
CommToolkit.modulation
generate_ofdm_symbol
ofdm_tx
apply_multipath
apply_channel
extract_ofdm_symbols
ofdm_rx_known
ofdm_rx_ls
ofdm_rx_mmse
interpolate_channel
sto_estimate
cfo_estimate
cfo_correct
ofdm_sync
```

## MIMO

```@docs
mimo_channel
apply_mimo
ml_detect
zf_detect
mmse_detect
sic_detect
zf_precoder
svd_precoder
waterfill
svd_waterfill
alamouti_encode
alamouti_decode_2x1
alamouti_decode_2x2
```

## Pulses and figures

`pulse_shape` truncates the transmit convolution. `matched_filter` keeps the full matched-filter tail and does not choose the sampling phase; pass that phase as `offset` to `downsample_symbols`. These functions do not estimate timing.

```@docs
raised_cosine
rrc_filter
pulse_shape
CommToolkit.Pulse
matched_filter
downsample_symbols
plot_constellation
plot_eye_diagram
plot_pulse_response
```
