# 公开接口

## v1.0 承诺

从 1.0.0 起，下面「稳定导出」列出的名字是稳定面。`using CommToolkit` 带入的符号与这份清单相同，也与各子模块的 `export` 相同，并都出现在 [API Reference](index.md) 的 `@docs` 块里。删除、改名，或改变位置参数和已有返回值的含义，是破坏性变更，必须升主版本号，并写进 `CHANGELOG`。1.x 可以增加名字；增加本身不升主版本。

1.0.0 不增加算法，也不改已有函数的位置参数和返回值：清单停在 0.9.3 已经导出的名字上（0.9.1 把当时的导出写成契约；0.9.2 增加 `teaching_punctured_mcs_table` 与 `teaching_punctured_amc_thresholds`，把删余码率和 BICM 接进 MCS / Chase；0.9.3 增加 `exponential_correlation`、`kronecker_channel` 与 `apply_kronecker`，即平坦 Kronecker 空间相关）。1.1.0 按承诺只增不破：增加 `RSCode`、`rs_code`、`rs_shorten`、`rs_encode`、`rs_decode` 五个名字，即短 Reed–Solomon；已有名字一个不动。1.1.1 按同一加法规则增加 `EXIT` 一组名字（`EXITCurve`、`exit_mi_soft`、`exit_mi_hist`、`exit_j`、`exit_jinv`、`exit_demap_transfer`、`exit_decoder_transfer`、`plot_exit_chart`），即 EXIT 图教学互信息测量；已有名字的位置参数、返回值与语义都不变。1.2.0 按同一加法规则增加 `XApp` 一组名字（`XAppUE`、`UEKPI`、`XAppAction`、`FixedXAppPolicy`、`RuleXAppPolicy`、`xapp_actions`、`observe_kpi!`、`jain_index`、`XAppReport`、`simulate_xapp`），即教学多用户 KPI→动作迷你 xApp；已有名字一个不动，`simulate_link` 一个字符都没有改。这是契约页，不是教程。

未导出的名字不是契约。这包括下划线开头的内部函数，以及子模块里的常数（例如 `Coding.CONV_G1`、`Coding._rsc_encode`）。子模块可以写成 `CommToolkit.Modulation`，但 `using CommToolkit` 不把子模块名字本身带进当前作用域。

定义位置：

| 子模块 | 源文件 |
|---|---|
| `Metrics` | `src/metrics.jl` |
| `Scrambling` | `src/scrambling.jl` |
| `Modulation` | `src/modulation.jl` |
| `Channel` | `src/channel.jl`，并包含 `src/cdl.jl` |
| `Pulse` | `src/pulse.jl` |
| `Coding` | `src/coding.jl`，并包含 `src/turbo.jl`、`src/algebraic.jl`、`src/ldpc.jl` |
| `Interleaving` | `src/interleaving.jl` |
| `BICM` | `src/bicm.jl` |
| `EXIT` | `src/exit.jl` |
| `Link` | `src/link.jl` |
| `XApp` | `src/xapp.jl` |
| `Equalization` | `src/equalization.jl` |
| `OFDM` | `src/ofdm.jl` |
| `MIMO` | `src/mimo.jl` |
| `Visualization` | `src/visualization.jl`。星座、眼图、脉冲响应在扩展 `ext/CommToolkitPlotsExt.jl`，`using Plots` 之后才加载 |
| `Sync` | `src/sync.jl` |

1.1.2 是 1.x 内的增加：`simulate_link` 增加关键字 `harq`（`:chase` 默认，与 1.0.0 逐比特一致；`:ir` 是删余卷积行上的教学增量冗余），稳定清单不增删名字，已有位置参数和返回值不变。1.1.3 按同一加法规则增加 6 个名字（`CDLProfile`、`cdl_profile`、`CDLTaps`、`cdl_taps`、`apply_cdl`、`cdl_channel`），即 TR 38.901 CDL-A/B/C 延迟–功率剖面上的教学簇延迟线；已有名字一个不动。1.2.0 把路线 C 收进范围：`XApp` 子模块的教学多用户 KPI→动作迷你 xApp，增加 10 个名字，只增不破；它是教学闭环，不是 O-RAN。

`test/contract.jl` 要求 `pkgversion(CommToolkit) == v"1.2.0"`，并核对三件事。`names(CommToolkit)` 比本节清单多一个名字，即模块自己的 `CommToolkit`；去掉它之后清单与 `names` 相同。各子模块的 `names` 同样含有该子模块自己的名字；去掉这些模块名之后，子模块导出的并集比清单多且只多 `modulation`。每个清单里的名字，以及 `CommToolkit` 和 `CommToolkit.modulation`，都出现在 `docs/src/index.md` 的 `@docs` 块里。

## 稳定导出

### Metrics

- `ber`
- `ser`
- `fer`
- `evm`
- `snr_db_to_linear`

### Scrambling

- `prbs_polynomial`
- `lfsr_sequence`
- `scramble`
- `descramble`

### Modulation

- `ModulationScheme`
- `BPSK`
- `QPSK`
- `PSK`
- `DBPSK`
- `DQPSK`
- `QAM`
- `modulate`
- `demodulate`
- `demodulate_soft`
- `constellation`
- `bits_per_symbol`

### Channel

- `awgn`
- `rayleigh_channel`
- `block_rayleigh`
- `block_rayleigh_channel`
- `jakes_fading`
- `jakes_channel`
- `rice_fading`
- `rice_channel`
- `nakagami_envelope`
- `nakagami_fading`
- `nakagami_channel`
- `mrc`
- `mrc_channel`
- `apply_fir`
- `fir_channel`
- `CDLProfile`
- `cdl_profile`
- `CDLTaps`
- `cdl_taps`
- `apply_cdl`
- `cdl_channel`

### Pulse

- `matched_filter`
- `downsample_symbols`

### Coding

- `conv_encode`
- `viterbi_decode`
- `bcjr_decode`
- `sova_decode`
- `distance_spectrum`
- `puncture_mask`
- `puncture`
- `depuncture`
- `turbo_interleaver`
- `turbo_encode`
- `turbo_decode`
- `GF2m`
- `gf_add`
- `gf_mul`
- `gf_div`
- `gf_inv`
- `gf_pow`
- `hamming74_encode`
- `hamming74_decode`
- `bch157_encode`
- `bch157_decode`
- `RSCode`
- `rs_code`
- `rs_shorten`
- `rs_encode`
- `rs_decode`
- `LDPCCode`
- `LDPCResult`
- `gallager_ldpc`
- `array_ldpc`
- `hamming74_ldpc`
- `ldpc_from_H`
- `ldpc_encode`
- `ldpc_decode`
- `ldpc_decode_details`
- `ldpc_syndrome`
- `ldpc_generator_matrix`
- `ldpc_column_degrees`
- `ldpc_row_degrees`

### Interleaving

- `block_interleave`
- `block_deinterleave`

### BICM

- `BICMResult`
- `bicm_labels`
- `bicm_psk8_id_labels`
- `bicm_map`
- `bicm_demap`
- `bicm_decode`

### EXIT

- `EXITCurve`
- `exit_mi_soft`
- `exit_mi_hist`
- `exit_j`
- `exit_jinv`
- `exit_demap_transfer`
- `exit_decoder_transfer`
- `plot_exit_chart`

### Link

- `CRC16_CCITT_POLY`
- `crc16_ccitt`
- `crc_append`
- `crc_check`
- `MCS`
- `mcs_convolutional`
- `mcs_ldpc`
- `teaching_mcs_table`
- `teaching_punctured_mcs_table`
- `spectral_efficiency`
- `ThresholdAMC`
- `OuterLoopAMC`
- `select_mcs`
- `observe_acknack!`
- `teaching_amc_thresholds`
- `teaching_punctured_amc_thresholds`
- `estimate_esn0_db`
- `LinkReport`
- `simulate_link`

### XApp

- `XAppUE`
- `UEKPI`
- `XAppAction`
- `FixedXAppPolicy`
- `RuleXAppPolicy`
- `xapp_actions`
- `observe_kpi!`
- `jain_index`
- `XAppReport`
- `simulate_xapp`

### Equalization

- `zf_equalizer`
- `zf_equalize`
- `mmse_equalizer`
- `mmse_equalize`
- `lms_equalize`
- `DFECoefficients`
- `dfe_equalizer`
- `dfe_equalize`
- `zf_dfe`
- `mmse_dfe`

### OFDM

- `OFDMParams`
- `TDLChannel`
- `generate_ofdm_symbol`
- `ofdm_tx`
- `apply_multipath`
- `apply_channel`
- `extract_ofdm_symbols`
- `ofdm_rx_known`
- `ofdm_rx_ls`
- `ofdm_rx_mmse`
- `interpolate_channel`
- `sto_estimate`
- `cfo_estimate`
- `cfo_correct`
- `ofdm_sync`

### MIMO

- `mimo_channel`
- `exponential_correlation`
- `kronecker_channel`
- `apply_mimo`
- `apply_kronecker`
- `ml_detect`
- `zf_detect`
- `mmse_detect`
- `sic_detect`
- `zf_precoder`
- `svd_precoder`
- `waterfill`
- `svd_waterfill`
- `alamouti_encode`
- `alamouti_decode_2x1`
- `alamouti_decode_2x2`

### Visualization

- `raised_cosine`
- `rrc_filter`
- `pulse_shape`
- `plot_constellation`
- `plot_eye_diagram`
- `plot_pulse_response`

### Sync

- `fractional_delay`
- `gardner_ted`
- `gardner_loop`
- `gardner_open_loop`
- `single_carrier_timing`
- `costas_detector`
- `costas_loop`
- `schmidl_preamble`
- `schmidl_cox`
- `integer_cfo_estimate`
- `schmidl_sync`

## 限定名

`CommToolkit.modulation` 根据 `OFDMParams` 选择 `QPSK()` 或 `QAM`。它从 `OFDM` 引入父模块，**不在** `export` 列表里，`using CommToolkit` 不会带入 `modulation`。文档字符串在 [API Reference](index.md)。这是仅有的导出差：子模块导出的并集比上面的稳定面多这一个名字。
