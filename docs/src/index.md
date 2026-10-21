# CommToolkit.jl API

API reference for the exported interface. This page is not a tutorial.

The package simulates textbook baseband links so that Monte Carlo error rates can be checked against closed-form formulas. It does not implement 3GPP NR or LTE, it does not implement an O-RAN fronthaul, and it is not a production physical-layer stack.

Names are defined in the submodules `Modulation`, `Channel`, `Coding`, `Equalization`, `OFDM`, `MIMO`, `Visualization`, and `Metrics`, and re-exported here. Plotting functions need Plots.jl; `using CommToolkit` does not load it.

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
