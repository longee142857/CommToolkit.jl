# CommToolkit.jl API

API reference for the exported interface of version 1.1.1. This page is not a tutorial. Those exported names are the stable surface; breaking changes require a major version bump. Versions 1.1.x are additive over 1.0.0: 1.1.0 adds short Reed–Solomon and 1.1.1 adds the `EXIT` submodule; neither changes an existing signature.

The stability surface, the validation matrix, and the scope page are written in Chinese, matching the README: [边界](boundaries.md), [公开接口](api.md), [验证矩阵](validation.md), [已知限制](limitations.md).

The package simulates textbook baseband links so that Monte Carlo error rates can be checked against closed-form formulas. It does not implement 3GPP NR or LTE PHY, it does not implement an O-RAN fronthaul, it does not implement a geometric 3GPP TR 38.901 clustered delay line, and it is not a production physical-layer stack. LDPC support is a short teaching reference (Gallager and array constructions), not an NR base-graph implementation. The turbo code is a short memory-2 teaching construction, not a 3GPP turbo code. `bch157_encode` is the binary BCH code of length 15 and dimension 7, not a general BCH constructor. `rs_code` adds short narrow-sense primitive Reed–Solomon codes over the same small fields (lengths 3, 7, 15, and shortened forms) with a hard-decision Berlekamp–Massey decoder; they are not the CCSDS long codes and not a soft algebraic decoder. The MCS table, CRC, and Chase combining loop are a link-level teaching reference, not an NR scheduler or HARQ entity. Convolutional MCS rows may use the periodic puncture masks; retransmissions still repeat one word. `receiver=:bicm` is a teaching BICM path on those convolutional rows, not NR BICM. `PSK` is coherent Gray ``M``-PSK. `DBPSK` and `DQPSK` are two-symbol differential detectors, not ``\pi/4``-DQPSK. `rice_channel` takes a linear Rice factor. `jakes_fading` is a flat Clarke spectrum on a DFT grid, not a classical Jakes sum of sinusoids. `nakagami_channel` takes the Nakagami parameter `m`. `mrc` is ideal-CSI receive combining, not a MIMO detector. `kronecker_channel` is one flat Kronecker matrix, not a Weichselberger model and not a TR 38.901 CDL or UMi.

Names are defined in the submodules `Modulation`, `Channel`, `Coding`, `Interleaving`, `BICM`, `EXIT`, `Link`, `Equalization`, `OFDM`, `MIMO`, `Visualization`, `Pulse`, `Scrambling`, `Sync`, and `Metrics`, and re-exported here. Plotting functions need Plots.jl; `using CommToolkit` does not load it. Gardner, Costas, and Schmidl–Cox in `Sync` are teaching loops with a stated capture range. They do not replace `ofdm_sync`.

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
jakes_fading
jakes_channel
rice_fading
rice_channel
nakagami_envelope
nakagami_fading
nakagami_channel
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

SOVA, the short turbo code, and the algebraic decoders are teaching implementations. They are not a 3GPP turbo code, not NR LDPC, and not a general BCH library; the Reed–Solomon codes of `rs_code` are the short narrow-sense primitive ones, not CCSDS long codes. `turbo_decode` is not part of the MCS / HARQ loop.

```@docs
conv_encode
viterbi_decode
bcjr_decode
sova_decode
distance_spectrum
puncture_mask
puncture
depuncture
turbo_interleaver
turbo_encode
turbo_decode
GF2m
gf_add
gf_mul
gf_div
gf_inv
gf_pow
hamming74_encode
hamming74_decode
bch157_encode
bch157_decode
RSCode
rs_code
rs_shorten
rs_encode
rs_decode
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

## BICM

Teaching bit-interleaved coded modulation. Not 5G NR BICM. The receiver loop itself draws no EXIT chart; soft-information mutual information and transfer characteristics are measured by the `EXIT` submodule below. With Gray labels and no bit prior, `bicm_demap` is `demodulate_soft`. `bicm_psk8_id_labels` is Ungerboeck set partitioning on 8-PSK; that is the map used when an iterative gain has to be visible. `:sova` and `:turbo` are one demapper pass plus the existing information-bit decoder. Convolutional `rate` of `2//3` or `3//4` punctures the mother code before the interleaver. The MCS loop calls this receiver through `simulate_link(...; receiver=:bicm)`; this module does not own the MCS table.

```@docs
CommToolkit.BICM
BICMResult
bicm_labels
bicm_psk8_id_labels
bicm_map
bicm_demap
bicm_decode
```

## EXIT charts

Teaching EXIT measurement, added in 1.1.1: mutual information of soft bit values from histograms (ten Brink 2001) or the histogram-free soft estimator (Hagenauer, EUSIPCO 2004), the Gaussian-symmetric `J(σ)` approximation and its inverse (ten Brink, Speidel, and Yan 2004) as a documented option, and transfer characteristics measured on the existing demapper and decoders without rewriting them. Not NR LDPC EXIT production tooling: no density evolution, no degree optimisation, no 3GPP codes, and no trajectory tracking. The chart figure needs Plots.jl.

```@docs
CommToolkit.EXIT
EXITCurve
exit_mi_soft
exit_mi_hist
exit_j
exit_jinv
exit_demap_transfer
exit_decoder_transfer
plot_exit_chart
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

Rectangular block interleaver. The MCS table does not select an interleaver. `receiver=:bicm` may pass `bicm_n_rows` into this same permutation after puncturing.

```@docs
CommToolkit.Interleaving
block_interleave
block_deinterleave
```

## Link

Teaching CRC, MCS tables, Chase combining, and threshold or outer-loop MCS selection. Convolutional rows can take puncture rates 1/2, 2/3, and 3/4. `receiver=:bicm` is Chase combining into a convolutional BICM decoder. Not an NR protocol and not incremental redundancy.

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
teaching_punctured_mcs_table
ThresholdAMC
OuterLoopAMC
select_mcs
observe_acknack!
teaching_amc_thresholds
teaching_punctured_amc_thresholds
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

`mimo_channel` is i.i.d. Rayleigh. `kronecker_channel` is one flat draw ``H = R_r^{1/2} H_w R_t^{1/2}``. `exponential_correlation` is Loyka's Toeplitz factor. Neither is a Weichselberger coupling model nor a 3GPP TR 38.901 CDL or UMi. Detectors take the realized `H` and are unchanged. `apply_kronecker` draws `H` and then calls `apply_mimo`.

```@docs
mimo_channel
exponential_correlation
kronecker_channel
apply_mimo
apply_kronecker
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

`pulse_shape` truncates the transmit convolution. `matched_filter` keeps the full matched-filter tail and does not choose the sampling phase; pass that phase as `offset` to `downsample_symbols`. The sampling phase and the carrier are estimated in `Sync`, not here.

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

## Synchronisation

Teaching Gardner timing, a BPSK/QPSK Costas loop, and a Schmidl–Cox preamble. Not a production synchroniser. `ofdm_sync` is unchanged: van de Beek, fractional CFO in ``(-1/2, 1/2)``, no integer CFO, buffer not shifted by the timing index.

```@docs
CommToolkit.Sync
fractional_delay
gardner_ted
gardner_loop
gardner_open_loop
single_carrier_timing
costas_detector
costas_loop
schmidl_preamble
schmidl_cox
integer_cfo_estimate
schmidl_sync
```
