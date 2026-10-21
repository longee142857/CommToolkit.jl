# CommToolkit.jl API

API reference for the exported interface. This page is not a tutorial.

The package simulates textbook baseband links so that Monte Carlo error rates can be checked against closed-form formulas. It does not implement 3GPP NR or LTE, it does not implement an O-RAN fronthaul, and it is not a production physical-layer stack. LDPC support is a short teaching reference (Gallager and array constructions), not an NR base-graph implementation. The MCS table, CRC, and Chase combining loop are a link-level teaching reference, not an NR scheduler or HARQ entity.

Names are defined in the submodules `Modulation`, `Channel`, `Coding`, `Link`, `Equalization`, `OFDM`, `MIMO`, `Visualization`, and `Metrics`, and re-exported here. Plotting functions need Plots.jl; `using CommToolkit` does not load it.

```@docs
CommToolkit
```

## Modulation

```@docs
ModulationScheme
BPSK
QPSK
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
ber
```

## Coding

```@docs
conv_encode
viterbi_decode
bcjr_decode
distance_spectrum
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
alamouti_encode
alamouti_decode_2x1
alamouti_decode_2x2
```

## Pulses and figures

```@docs
raised_cosine
rrc_filter
pulse_shape
plot_constellation
plot_eye_diagram
plot_pulse_response
```
