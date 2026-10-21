# 变更记录

## 0.5.0

链路级教学闭环：CRC-16、离散 MCS、Chase 合并、按 ``E_s/N_0`` 切换的自适应 MCS。不是 3GPP NR/LTE 调度器，没有 CQI 表、速率匹配、HARQ 进程、时序或多用户。

- CRC-16：生成多项式 ``x^{16}+x^{12}+x^5+1``（CCITT / ITU-T V.41）。系统循环码约定：初值 0、不反射、异或输出 0，与 CRC-16/XMODEM 的参数一致。ASCII `"123456789"` 的余数是 `0x31C3`。
- MCS：五档教学参数（BPSK/QPSK/16-QAM/64-QAM 与终止卷积码，外加 QPSK 上的 Fan 阵列 LDPC `array_ldpc(3, 8, 11)`）。频谱效率是每符号的净荷比特，CRC 和卷积分组的尾比特算开销。
- 链路：净荷 → CRC → 编码 → 调制 → AWGN → 软解调 → 译码 → CRC。``E_s = 1``，``\sigma^2 = N_0/2``。Chase 合并把同一编码块的比特 LLR 相加，``max_tx`` 次失败后丢弃。
- 自适应：门限策略用理想 ``E_s/N_0``，门限来自固定 MCS 的吞吐包络。外环按 ACK/NACK 移动一个 dB 偏移。可选上一块的二阶矩 ``E_s/N_0``。没有功率控制，没有增量冗余。

## 0.4.0

短 LDPC 的教学 / 参考实现。不是 3GPP NR / LTE 基图，不含速率匹配、CRC 附加或 HARQ。

- 码：Gallager 正则块构造（确定性 SplitMix64 列置换）、Fan 阵列码、``[7,4]`` Hamming 夹具，以及 `ldpc_from_H`。校验矩阵为稀疏 0/1。
- 编码：构造时一次 GF(2) 高斯消元，系统形式。行数和列数不超过 4096。
- 译码：泛洪 log-domain sum-product、min-sum、归一化 min-sum（默认，``\alpha = 0.75``）、偏移 min-sum（``\beta = 0.15``）。输入为信道 LLR，正 LLR 表示比特 0；可用迭代上限，并以综合征为 0 早停。
- 噪声与 v0.3 的 `awgn` / `demodulate_soft` 一致：``\sigma^2 = N_0/2``。编码 BPSK 取 ``E_s = 1``、``E_b = n/k``。

## 0.3.0

AWGN 逐比特软判决 LLR（精确边际与 max-log），软输入硬输出 Viterbi，终止卷积码上的 log-BCJR。

## 0.2.0

规范 Julia 包布局，以及与闭式公式对照的蒙特卡洛误码率验证。
