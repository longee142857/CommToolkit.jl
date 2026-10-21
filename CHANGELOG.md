# 变更记录

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
