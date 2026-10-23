# 变更记录

## 0.8.0

RMS EVM、加性 LFSR 扰码，以及 SVD 预编码上的总功率注水。不是 Turbo，没有 Jakes，没有 Gardner，没有 BICM，也没有接到 MCS / HARQ 上。

- 指标：`evm`。参考符号对接收符号的 RMS 比。默认 `unit=:percent`，是 ``100`` 乘幅度比，不是误差功率占参考功率的百分比。`unit=:db` 是 ``20\log_{10}``（幅度）也就是 ``10\log_{10}``（功率）；10% 对应 ``-20``，完全对齐时是 `-Inf`。`unit=:linear` 是无量纲幅度比。不做幅度、相位或直流拟合。
- 扰码：`prbs_polynomial`、`lfsr_sequence`、`scramble`、`descramble`，放在 `Scrambling`。Fibonacci LFSR，比特 0 先出；`seed=nothing` 是全 1。经典三项式次数 7、9、11、15、20、23、29、31。次数 7 是 ``x^7+x^6+1``，不是 802.11 的 ``x^7+x^4+1``。次数 9 及以上取 ITU-T O.150 的连接多项式；次数 20 是 ``x^{20}+x^3+1``，不是另一条 ``x^{20}+x^{17}+1``。输出不取反、不塞 1，不是 O.150 测试图样。同一多项式和种子再异或一次即解扰，多项式不必本原。不是自同步扰码，不是 Gold 码。
- MIMO：`waterfill` 在奇异值（不是 ``H^{\mathsf H}H`` 的特征值）上做水位线。`svd_waterfill` 把 ``\sqrt{p_i}`` 乘进 `svd_precoder` 的右奇异向量；功率为 0 的模式保留成零列。`noise_var` 与 `apply_mimo` 一样是 ``\mathbb{E}[|n_i|^2]``，不是 ``N_0/2``。`svd_precoder` 的返回值不变，仍然不分配功率。
- 版本：已有函数的位置参数和返回值不变，所以不是破坏性变更。这一刀新开了扰码模块，并给指标和 MIMO 各加了一组接口，不是 0.7.1 上的小补丁，按 0.x 记为 0.8.0，而不是 0.7.2，也不升到 1.0。

## 0.7.1

符号率 LMS，以及和块 ZF / MMSE 同一卷积矩阵上的有限长块 DFE。不是 RLS，不是盲均衡，没有接到 MCS / HARQ 上。

- LMS：`lms_equalize`。复基带约定是 ``\hat{y} = w^{\mathsf{T}} u``，``w \leftarrow w + \mu e u^*``；实信号上共轭是恒等。步长、抽头数、时延和训练序列由调用者给定。`tap_trajectory=true` 时另外返回每次监督更新之后的抽头。训练段之后抽头冻结。
- DFE：`dfe_equalizer` 给出每个符号一行的前馈，以及严格下三角的反馈。`noise_var = 0` 是 ZF，正的有限噪声方差是 MMSE，Gram 矩阵与 `mmse_equalizer` 相同。`zf_dfe` / `mmse_dfe` 把反馈判决取成 `alphabet` 上的最近邻。`training` 覆盖到的位置改用已知符号，不经过切片。
- 版本：已有均衡器的位置参数和返回值不变。这一刀补在 0.7.0 的块线性均衡旁边，不新开调制或信道族，所以记为补丁 0.7.1，而不是 0.8.0，也不升到 1.0。

## 0.7.0

Gray M-PSK、DBPSK/DQPSK、平坦 Rice，以及理想 CSI 的 L 支路 MRC。不是 ``\pi/4``-DQPSK，没有 Jakes，没有接到 MCS / HARQ 上。

- 调制：`PSK(M)`，``M \ge 4`` 且为 2 的幂。二进制反射 Gray，符号 0 在正实轴，单位能量。硬判决是最小距离；软 LLR 走原来的相干 AWGN 边际（精确与 max-log）。`PSK(4)` 不是 `QPSK()` 的比特映射。
- 差分：`DBPSK`、`DQPSK`。`modulate` 在块首放参考符号 ``1``，输出长度比数据符号多 1。硬判决看 ``y_i y_{i-1}^*``。软 LLR 是两符号、公共相位未知时的 ``I_0`` 似然。DQPSK 的相位增量是 ``0, \pm\pi/2, \pi``，不是 ``\pi/4``。
- 信道：`rice_fading` / `rice_channel`。`K` 是线性 ``P_{\mathrm{LOS}}/P_{\mathrm{scatter}}``，不是 dB；直射相位为 0；``\mathbb{E}[|h|^2] = 1``。噪声仍是 `awgn` 的 ``E_b/N_0``。`K = 0` 与瑞利同分布，`K = Inf` 是 ``h = 1``。
- 分集：`mrc` 对 ``L \times N`` 的已知 ``h`` 做最大比合并并按 ``\sum|h|^2`` 归一化。`mrc_channel` 生成 L 条独立 Rice/瑞利支路（默认 `K = 0`）再合并。每条支路的平均 ``E_b/N_0`` 是调用参数，不是合并之后的信噪比。
- 对照：Gray M-PSK 扇区积分（Craig 符号错误率交叉核对）、DBPSK 指数式、DQPSK 的 Marcum Q、Rice 上相干 BPSK 的 Craig 积分、瑞利 MRC 的 Proakis 求和式。实现在 `test/reference.jl`，不从包里导出。
- 版本：已有函数的位置参数和返回值不变，所以不是破坏性变更。这一刀是一组新的调制和信道，不是 0.6.1 上的小补丁，按 0.x 记为 0.7.0，而不是 0.6.2，也不升到 1.0。

## 0.6.1

单载波接收侧匹配滤波、按每符号样本数抽样，以及确定性 FIR 多径。不是定时恢复，不是 Jakes，也没有接到 OFDM 的抽头延迟线上。

- 脉冲：`matched_filter` 对 `conj.(reverse(taps))` 做全长线性卷积，尾部保留，不去掉群时延。`downsample_symbols` 从 `offset + 1` 起每 `sps` 个样本留一个；`offset` 不按 `sps` 取模。
- 信道：`apply_fir` 与 ZF / MMSE 用的是同一套全长卷积，无噪声时 ZF 可以恢复符号。`fir_channel` 在卷积后按 `awgn` 的 ``E_b/N_0`` 加噪，并返回该噪声向量每一项的方差（实数为 ``N_0/2``，复数为 ``N_0``）。
- 版本：0.6.0 已经包含交织、删余和 `ser` / `fer`。这一刀只增加接口，不改已有调用，记为补丁版本 0.6.1。

## 0.6.0

矩形块交织、卷积母码的标准删余、符号 / 帧错误率，以及块恒定的瑞利衰落。不是 3GPP 速率匹配，没有 Jakes 时间相关，删余和交织都没有接到 MCS / HARQ 闭环上。

- 交织：`block_interleave` / `block_deinterleave`。按行写入、按列读出。对比特、符号和 LLR 是同一种置换。
- 删余：母码仍是终止的 ``K = 3``、(7, 5)。`puncture_mask` 给出 1/2、2/3、3/4。2/3 的图案是 ``[1, 1, 1, 0]``，3/4 是 ``[1, 1, 1, 0, 0, 1]``，按 ``(g_1, g_2)`` 的顺序周期作用，含尾比特。`depuncture` 把删余位填成调用者指定的值；软译码路径填 LLR 0。硬判决 Viterbi 不认识擦除。
- 指标：`ser` 与 `ber` 一样按元素比较；`fer` 按给定块长统计至少含一个错误的块。
- 信道：`block_rayleigh` 生成一块内恒定的 ``\mathcal{CN}(0,1)`` 系数，`block_rayleigh_channel` 按与 `rayleigh_channel` 相同的 ``E_b/N_0`` 约定加噪。块之间独立。

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
