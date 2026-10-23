# CommToolkit.jl

基带通信链路仿真库（v0.8.5），用来把蒙特卡洛误码率和教材中的闭式公式对上，并在 AWGN 上跑一条教学级的 MCS / HARQ 闭环。作者刘恒宇（longee142857）。

本包**不是** 3GPP NR / LTE 物理层，**不是** O-RAN 前传实现，也**不是**生产级协议栈。OFDM 是教科书循环前缀模型，没有资源栅格、DMRS 或空口时隙结构。LDPC 是短码教学 / 参考实现，不是 NR BG1/BG2。Turbo 是记忆 2 的短帧教学码，不是 3GPP Turbo。BCH 只有二进制 (15, 7) 这一种，不是 Reed–Solomon。MCS 表、CRC 和 Chase 合并是链路级教学闭环，不是 NR CQI/MCS 表，也不是 HARQ 实体或调度器。卷积删余是 ``K = 3`` 母码上的周期图案，不是 NR 速率匹配。块衰落在一块内系数恒定。Clarke/Jakes 时间相关是另一条平坦过程，不是频率选择性抽头，也不是 Jakes 1974 的正弦和。

## 功能

公开接口从子模块扁平再导出，原有调用方式保持不变（例如 `modulate(bits, BPSK())`、`awgn(syms, snr_db, bits_per_symbol)`）。没有重命名，因此没有弃用别名。

- **调制**（`Modulation`）：BPSK、Gray QPSK、Gray M-PSK（`PSK(M)`，``M \ge 4`` 且为 2 的幂，符号 0 在正实轴）、方形 Gray M-QAM（`M ≥ 4` 且为平方、2 的幂）。单位符号能量。硬判决解调；相干 AWGN 下逐比特 LLR（精确边际与 max-log），噪声方差与 `awgn` 的 ``\sigma^2 = N_0/2`` 相同。另有 DBPSK / DQPSK：`modulate` 在块首放参考符号 ``1``，硬判决是两符号差分检测，软 LLR 是未知公共相位下的 ``I_0`` 边际。不是 ``\pi/4``-DQPSK。
- **信道**（`Channel`）：AWGN，`snr_db` 是 ``E_b/N_0``（dB）；平坦瑞利衰落，可选理想信道补偿。平坦 Rice：`rice_fading` / `rice_channel`，`K` 是线性 ``P_{\mathrm{LOS}}/P_{\mathrm{scatter}}``，不是 dB，``\mathbb{E}[|h|^2] = 1``，直射相位为 0，噪声仍走 `awgn`。Clarke/Jakes 时间相关：`jakes_fading` / `jakes_channel`，最大多普勒 `fd`（Hz）和采样间隔 `Ts`（秒），``\mathbb{E}[|h|^2] = 1``，谱是 Clarke 的 U 形，实现是 DFT 格子上对 ``S(f)`` 积分后再做 `ifft`，不是正弦和。Nakagami-m：`nakagami_envelope` / `nakagami_fading` / `nakagami_channel`，`m ≥ 1/2`，``\mathbb{E}[R^2] = 1``；`m = 1` 与瑞利同一套 `randn` 顺序；`m = Inf` 是单位幅度加均匀相位。理想 CSI 的 L 支路 MRC：`mrc` 合并已知的 ``L \times N`` 观测，`mrc_channel` 生成独立 Rice/瑞利支路再合并。另有块衰落：同一块内一个 ``h \sim \mathcal{CN}(0,1)``，块之间独立。单载波 FIR：`apply_fir` 做全长线性卷积，`fir_channel` 在卷积之后按同一 ``E_b/N_0`` 约定加 AWGN，并返回 `mmse_equalize` 用的噪声方差。抽头由调用者给定，不随时间变化。
- **编码**（`Coding`）：码率 1/2、约束长度 ``K = 3``、八进制生成多项式 (7, 5) 的卷积码，尾比特终止。硬判决 Viterbi、以 LLR 为分支度量的软输入硬输出 Viterbi、终止码上的 log-BCJR 外信息、同一条终止网格上的 SOVA（Hagenauer–Hoeher 竞争路径可靠度），以及网格距离谱。母码可用标准周期图案删余到名义码率 2/3 和 3/4；译码前把删余位填回 LLR 0。短帧 Turbo：两路记忆 2 的递归系统卷积 (1, 5/7)，信息比特只发一次，块交织或步长交织，log-BCJR 迭代，信息长度不超过 1024，迭代不超过 32。代数硬判决：``\mathrm{GF}(2^m)``（``m=2,3,4``）、Hamming(7,4) 伴随式译码（与 `hamming74_ldpc` 同一组码字，但是硬判决）、二进制 BCH(15,7)（设计距离 5，纠正 2 个错误）。另有短 LDPC 的教学实现：Gallager 正则构造与 Fan 阵列码（稀疏校验矩阵）、构造时一次 GF(2) 高斯消元得到的系统编码、泛洪 log-domain sum-product 与 min-sum（默认归一化 min-sum，``\alpha = 0.75``），综合征为 0 时早停。信道 LLR 与 `demodulate_soft` 相同：正 LLR 表示比特 0。不是 NR 基图，不是 3GPP Turbo，没有 BICM-ID。
- **交织**（`Interleaving`）：矩形块交织，按行写入、按列读出，对比特、符号和 LLR 是同一种置换。
- **链路**（`Link`）：CRC-16（CCITT 多项式 ``x^{16}+x^{12}+x^5+1``，系统码，初值 0）、五档教学 MCS、Chase 合并（比特 LLR 相加，有限次重传）和按 ``E_s/N_0`` 切换的 MCS。仿真信噪比是 ``E_s/N_0``，``E_s = 1``，``\sigma^2 = N_0/2``，与 `awgn` / `demodulate_soft` 相同。返回 BLER、BER、吞吐（送达净荷比特 / 信道符号，含重传）和频谱效率。没有速率匹配、增量冗余、多进程 HARQ、调度或功率控制。
- **均衡**（`Equalization`）：频率选择性信道上的线性块 ZF / MMSE。无噪声时反演 `apply_fir` 的全长卷积。同一条卷积矩阵上的有限长块 DFE（`dfe_equalizer`，`noise_var = 0` 为 ZF，正的噪声方差为 MMSE；`zf_dfe` / `mmse_dfe` 做最近邻反馈）。符号率 LMS（`lms_equalize`）：训练序列监督，步长和抽头数由调用者给定，可返回抽头轨迹；训练段之后抽头冻结。不是 RLS，不是盲均衡。
- **OFDM**（`OFDM`）：循环前缀 OFDM 收发、抽头延迟线、已知 CSI / LS / 逐导频收缩后线性插值、van de Beek 符号定时与小数倍频偏。van de Beek 的捕获范围是 ``(-1/2, 1/2)`` 个子载波间隔，不估计整数倍频偏。Schmidl–Cox 前导在 `Sync` 里，不改 `ofdm_sync`。
- **MIMO**（`MIMO`）：i.i.d. 瑞利信道，ML、ZF、MMSE、按行范数排序的 SIC，ZF 与 SVD 预编码，SVD 上的总功率注水（`waterfill` / `svd_waterfill`），Alamouti ``2\times1`` 与 ``2\times2``。`svd_precoder` 本身仍然不分配功率。
- **扰码**（`Scrambling`）：加性 LFSR。`scramble` 与 `descramble` 使用同一条从种子重起的序列。`prbs_polynomial` 给出次数 7、9、11、15、20、23、29、31 的经典三项式。不是自同步扰码，不是 Gold 码，没有接到链路闭环上。
- **脉冲**（`Visualization`、`Pulse`）：升余弦与根升余弦，`pulse_shape`。接收侧 `matched_filter` 用发射脉冲的共轭反转做全长卷积，`downsample_symbols` 按每符号样本数抽一个相位。星座图、眼图、脉冲频谱由包扩展 `CommToolkitPlotsExt` 提供，只有 `using Plots` 之后才加载。
- **同步**（`Sync`）：Gardner 定时误差（偶数 `sps`）、线性内插的一阶环和开环过零，以及 `pulse_shape` → `matched_filter` → 分数延迟 → `downsample_symbols` 的胶水。Costas 环只覆盖 `BPSK()` 和 `QPSK()`。Schmidl–Cox 重复半符号前导做帧检测和 ``(-1, 1)`` 内的小数 CFO，偶数倍整数 CFO 另搜。不是生产级同步。
- **指标**（`Metrics`）：误码率、符号错误率、按块计的帧错误率、RMS 误差矢量幅度（默认百分数，也可 dB 或线性比）、dB 与线性信噪比换算。

`CommToolkit.modulation` 根据 `OFDMParams` 选择 QPSK 或 QAM，不在 `export` 列表里。

## 安装

需要 Julia 1.10 或更新版本。

```julia
using Pkg
Pkg.add(url="https://github.com/longee142857/CommToolkit.jl")
using CommToolkit
```

`using CommToolkit` 不加载 Plots。作图时再执行 `using Plots`。

```julia
bits = rand(0:1, 10_000)
rx = awgn(modulate(bits, BPSK()), 8, 1)  # 8 dB Eb/N0
ber(bits, demodulate(rx, BPSK()))
```

## 验证

下图是固定种子 `0x0c0de1162` 的蒙特卡洛点与闭式曲线。实线是理论值，标记是仿真。生成脚本是 `demos/validation_figure.jl`。统计检验在 `test/`：双侧 Clopper–Pearson 区间（``\alpha = 10^{-4}``，并要求足够的错误计数），理论值必须落在区间内。公式与出处见该脚本和 `test/reference.jl`。

- AWGN 下 BPSK 与 Gray QPSK：``Q(\sqrt{2 E_b/N_0})``。Proakis 与 Salehi，《Digital Communications》第 5 版，§5.2。
- 方形 Gray M-QAM 的精确平均误比特率：Cho 与 Yoon，IEEE Trans. Commun. 50(7):1074–1080，2002，式 (16)。图中 16-QAM 用的是这个精确式，不是高信噪比近似。
- 理想 CSI 下平坦瑞利相干 BPSK：``\tfrac12\bigl(1 - \sqrt{\bar\gamma/(1+\bar\gamma)}\bigr)``，``\mathbb{E}[|h|^2] = 1``。Goldsmith，《Wireless Communications》，2005，§6.3；Proakis 与 Salehi 第 5 版 §13.3。
- Gray M-PSK：相位误差扇区上的 Gray 汉明距离平均（Proakis 与 Salehi §4.2、§5.2）。符号错误率同时用 Craig 积分交叉核对。最近邻近似 ``(1/\log_2 M)\operatorname{erfc}(\sqrt{\log_2 M\cdot E_b/N_0}\sin(\pi/M))`` 在 ``M = 4`` 时等于 QPSK，``M > 4`` 不是精确式。蒙特卡洛对照的是扇区积分。
- DBPSK：``\tfrac12\exp(-E_b/N_0)``。Gray DQPSK：一阶 Marcum Q 比特错误率（Proakis 与 Salehi §4.5）。差分检测的参考符号是前导，公式不把它折进信息比特能量。
- 理想 CSI 相干 BPSK、平坦 Rice：Craig 形式的 ``Q`` 对 Rician 矩生成函数积分。`K = 0` 退回瑞利闭式，`K \to \infty` 退回 AWGN。`K` 是线性比。
- L 支路 i.i.d. 瑞利、理想 CSI、相干 BPSK 的 MRC：Proakis 与 Salehi §13.4 的求和式。``\bar\gamma_b`` 是每条支路的平均 ``E_b/N_0``。`L = 1` 退回单支路瑞利。
- Clarke/Jakes：滞后相关对照 ``J_0(2\pi f_d \tau)``。有限 DFT 格子和有限样本都留有误差。理想 CSI 相干 BPSK 的平均仍用瑞利闭式，因为边际是 ``\mathcal{CN}(0,1)``。
- Nakagami-m、理想 CSI、相干 BPSK：Craig 形式对 Gamma 矩生成函数积分。`m = 1` 退回瑞利。正整数 `m` 等于 `m` 支路瑞利 MRC，每条支路的平均是 ``\bar\gamma_b/m``。`m \to \infty` 在理想 CSI 下退回 AWGN。

![BPSK、QPSK、16-QAM 在 AWGN 下以及 BPSK 在平坦瑞利下的蒙特卡洛误码率与闭式曲线](docs/src/assets/ber_validation.png)

上图仍是 v0.6 之前的四条曲线。Gray M-PSK、DBPSK、DQPSK、Rice 上的相干 BPSK 和瑞利 MRC 用同一套 Clopper–Pearson 区间（``\alpha = 10^{-4}``）写在 `test/psk_rice_mrc.jl`，公式在 `test/reference.jl`。可复现的打印在 `demos/ber_psk_rice_mrc.jl`，不依赖 Plots。Jakes 的滞后相关和 Nakagami-m 的相干 BPSK 在 `test/jakes_nakagami.jl`，打印在 `demos/ber_jakes_nakagami.jl`，也不依赖 Plots。

卷积码在测试里与硬判决联合界比较（Proakis 与 Salehi §8.2；Lin 与 Costello，《Error Control Coding》第 2 版 §12.3）。同一 ``E_b/N_0`` 下，软输入 Viterbi 与 BCJR 硬判决的 Clopper–Pearson 区间低于硬判决 Viterbi，也低于未编码 BPSK，并落在软判决联合界 ``\sum c_d Q(\sqrt{2 d R E_b/N_0})`` 的同一数量级（Proakis 与 Salehi §8.2；Bahl、Cocke、Jelinek 与 Raviv，IEEE Trans. Inf. Theory 20(2):284–287，1974）。BPSK / Gray QPSK 的 LLR 与闭式条件高斯分布对照（Proakis 与 Salehi §4.2、§5.2）。OFDM / MIMO 测试覆盖无噪声完美恢复，以及高低信噪比下互不重叠的 Clopper–Pearson 区间。单载波 FIR 与匹配滤波在 `test/pulse_fir.jl`：无噪声时 `zf_equalize` 反演 `apply_fir`；冲激脉冲走 `pulse_shape`、样本率 FIR、`matched_filter`、`downsample_symbols` 之后，ZF 与 `noise_var = 0` 的 MMSE 恢复符号。无多径、无噪声时，根升余弦级联在符号点上贴近原符号。可复现的打印在 `demos/pulse_fir_link.jl`，不依赖 Plots。

`test/lms_dfe.jl` 把符号率 `apply_fir` / `fir_channel` 接到 LMS 和块 DFE。无噪声时 LMS 的抽头朝这条回归量的批最小二乘解走，末段误差低于开头；同一条确定性噪声记录上，更大的稳定步长留下更大的末段均方误差，不对照 Widrow 超额 MSE 闭式。ZF-DFE 在判决正确时把全列秩 FIR 恢复到舍入误差；第一个反馈符号被换成错误值时，后面的硬判决可以成串翻掉。

`test/sync.jl` 不跑误码率曲线。Gardner 在根升余弦上核对 S 曲线的符号和开环过零（相对真抽样相位的圆周距离），一阶环在偏晚和偏早两种初值下把抽样时刻拉向自相关峰。Costas 在小相偏和小频偏下核对鉴相器符号、相位和频率状态的移动方向，以及锁定后的 QPSK 硬判决。Schmidl–Cox 在无噪声前导上核对两半相等、小数 CFO、偶数倍整数 CFO，以及 ``\varepsilon = 0.7`` 时 van de Beek 仍报成接近 ``-0.3``。可复现的打印在 `demos/sync_timing.jl`，不依赖 Plots。

LDPC 使用同一套 `awgn` / `demodulate_soft` 标定：BPSK 符号能量 ``E_s = 1``，信息比特能量 ``E_b = n/k``，``\sigma^2 = N_0/2``。`test/ldpc.jl` 在 Fan 阵列码 `array_ldpc(3, 6, 11)`（``n = 66``，``k = 35``，无 4 环）上，于 2.5 dB 与 3.5 dB、各 1000 帧、种子 `0x0c0de1162`，对归一化 min-sum 与 sum-product 检查信息比特 BER 和 BLER 的 Clopper–Pearson 上界（``\alpha = 10^{-4}``）低于未编码 BPSK。小码与 ``[7,4]`` Hamming 码做穷举对照（Lin 与 Costello §3.2）；树上的 sum-product 与逐比特 MAP 一致。Hamming 校验图含 4 环，不是最大似然译码。

链路测试在 `test/link.jl`。CRC-16 对照 CRC-16/XMODEM 的 `"123456789" → 0x31C3`，并检查单比特翻转和 ``g(x)`` 自身作为漏检图样。无噪声时每一档 MCS 的吞吐等于该档频谱效率，Chase 不再重传。``E_s/N_0 = 5\,\mathrm{dB}`` 的 16-QAM 卷积码上，3 次 Chase 合并相对单次传输的 CRC 失败率区间（Clopper–Pearson，``\alpha = 10^{-4}``，100 个块）互不重叠，合并后的吞吐更高。同一 ``\alpha`` 下，5 dB 的门限自适应（选中 QPSK 阵列 LDPC）吞吐下界高于固定 BPSK，也高于固定 64-QAM；16 dB 的自适应（选中 64-QAM）下界高于固定 BPSK。门限本身来自 Julia 1.10.12、种子 `0x0c0de1162`、每点 80 块的固定 MCS 吞吐包络，见 `teaching_amc_thresholds`。可复现的打印表在 `demos/link_mcs_harq.jl`，不依赖 Plots。

本地运行：

```julia
using Pkg
Pkg.test("CommToolkit")
```

## 源码布局

```
src/CommToolkit.jl    模块入口与再导出
src/modulation.jl     Modulation
src/channel.jl        Channel（含 Rice、Clarke/Jakes、Nakagami-m、理想 CSI MRC、单载波 FIR）
src/pulse.jl          Pulse：匹配滤波与按 sps 抽样
src/coding.jl         Coding（包含 ldpc.jl、turbo.jl、algebraic.jl；卷积与 SOVA 在本文件）
src/turbo.jl          短帧 Turbo：记忆 2 RSC、块/步长交织、log-BCJR 迭代
src/algebraic.jl      GF(2^m)、Hamming(7,4) 伴随式硬判决、BCH(15,7)
src/interleaving.jl   Interleaving：矩形块交织
src/ldpc.jl           短 LDPC：Gallager、阵列码、系统编码、泛洪译码
src/link.jl           Link：CRC-16、教学 MCS、Chase 合并、自适应选择
src/equalization.jl   Equalization
src/ofdm.jl           OFDM
src/mimo.jl           MIMO
src/scrambling.jl     Scrambling：加性 LFSR
src/visualization.jl  Visualization（作图函数在未加载 Plots 时抛错）
src/sync.jl           Sync：Gardner、Costas、Schmidl–Cox
src/metrics.jl        Metrics
ext/CommToolkitPlotsExt.jl
test/                 闭式对照与单元测试
docs/                 Documenter API 参考（没有教程）
demos/                可单独运行的脚本，不在测试范围内
```

API 参考由 Documenter 从导出函数的 docstring 生成，页面只含接口说明。

## 许可

MIT。Copyright (c) 2026 刘恒宇 (longee142857)。见 `LICENSE`。

## 已知限制

- 卷积码母码仍是 (7, 5)、``K = 3``。删余图案只有两张：1/2→2/3 的 ``[1, 1, 1, 0]`` 和 1/2→3/4 的 ``[1, 1, 1, 0, 0, 1]``，按 ``(g_1, g_2)`` 码流顺序周期作用，尾比特一并删余。这是 Yasuda 等人与 IEEE 802.11 用在 1/2 母码上的图案，不是针对 (7, 5) 重新搜索的自由距表。没有其它码率。``depuncture`` 在软 LLR 上把删余位填 0（擦除）。硬判决 Viterbi 只接受 0/1，不把填进去的数当成擦除；删余后的硬判决要先换成 ``1-2\hat{c}`` 再填 0，送给软 Viterbi、SOVA 或 BCJR。软输入 Viterbi 输出硬比特，并且和硬判决译码一样不把最终状态钉在 0。BCJR 和 SOVA 按终止码把尾比特固定为 0、最终状态固定在 0。BCJR 返回信息比特外信息（后验 LLR 减去先验）。SOVA 的外信息是 ``(1-2\hat{u})\rho-L_a``，``\rho`` 是幸存路径上竞争路径度量差的最小值，不是路径求和；高信噪比下符号与 BCJR 一致，幅度不要求相等，而且常常偏大。整块更新，块长很大时是平方复杂度。码是非系统的，没有信道系统位可减。联合界仍按未删余、未终止母码的第一事件界计算，不覆盖删余码。
- Turbo 只做短帧教学码。分量码是记忆 2 的递归系统卷积 (1, 5/7)，不是 `conv_encode` 那条前向码，也不是 3GPP 的 (13, 15) 或 (15, 13)。信息长度 2 到 1024。码字长度是 ``3k+8``，不是精确 1/3：两路各有 2 个不进交织器的尾比特。交织器在长度可分解时等于矩形块交织，质数长度是步长置换，索引 1 可以是不动点；不是 QPP，也不保证是错位排列。迭代 1 到 32 次，默认 6，用的是 RSC 网格上的 log-BCJR，不是 SOVA。一次迭代不是级联码的最大似然译码，增加迭代也不保证每一帧的错误都下降。没有打孔，没有 3GPP 速率匹配，没有接到 MCS / HARQ 上。
- 有限域只实现 ``m=2,3,4``，各一个固定本原多项式。元素是整数比特掩码，不是独立的元素类型。没有 ``\mathrm{GF}(2^8)``，没有一般的 Conway 搜索。Hamming(7,4) 硬判决只纠正 1 个错误；两个错误一定被改到另一个码字，函数仍返回 4 个比特，没有失败标志。信息位固定在第 3、5、6、7 位，和 `ldpc_encode` 的主元顺序不一定相同，但 16 个码字是同一集合。`hamming74_ldpc` 上的泛洪 sum-product 不是这个硬判决译码器。BCH 只有 (15, 7)、设计距离 5 这一种，生成多项式 ``x^8+x^7+x^6+x^4+1``，比特从高次到低次，系统位在前。Peterson–Chien 在定位失败时把收到的 7 个系统位原样返回。3 个错误可能译错，也可能被发现后不改。不是 Reed–Solomon，没有一般的 ``n,k,t`` 构造器，也没有软判决代数译码。
- 交织只有矩形块交织（按行写入、按列读出）。没有卷积交织，也没有和删余图案或 MCS 绑在一起的交织器。逐符号独立衰落上，块交织不改变每个符号的边际信噪比。
- 平坦瑞利 `rayleigh_channel` 仍是逐符号独立的 ``h``。块衰落在一块内 ``h`` 恒定、块之间独立，``\mathbb{E}[|h|^2] = 1``。时间相关平坦瑞利在 `jakes_fading`：Clarke 谱，`n` 点 DFT 上对 ``S(f)`` 做格子积分，不是 Jakes 1974 的有限正弦和，也不是 IIR。`fd` 的单位是 Hz，`Ts` 的单位是秒，`fd * Ts < 1/2`。集合自相关是格子功率的逆 DFT，只在 `n * fd * Ts` 比较大时接近 ``J_0``。序列以 `n` 为周期。一条平坦系数，不是抽头延迟线，没有接到 OFDM 的 `TDLChannel`。``rx ./ h`` 之后各符号噪声方差是 ``\sigma^2/|h|^2``，不能再拿一个标量方差做软解调。相干 BPSK 的 LLR 是 ``2\Re(h^* y)/\sigma^2``。理想 CSI 下的平均误码率与逐符号瑞利相同；相关把错误聚成大约 ``1/(f_d T_s)`` 个样本的突发，不能把一次实现上的比特错误当成独立伯努利试验来套置信区间。
- M-PSK 只接受不小于 4 的 2 的幂。`PSK(4)` 的星座在坐标轴上，Gray 序号是二进制反射码，不是 `QPSK()` 的 ``\pi/4`` 正交映射。没有 ``\pi/4``-DQPSK，差分只做到 DBPSK 和 DQPSK。`modulate` 多输出一个参考符号 ``1``，``N`` 个信息比特占 ``N/k + 1`` 个信道使用；闭式误码率是每一对符号的错误率，不把这个前导折成信息比特能量损失。软 LLR 只用相邻两个符号、公共相位未知时的 ``I_0`` 似然，不做更长的多符号窗。逐符号独立衰落会毁掉差分相位；整块乘同一个相位则不会。
- Rice 的 `K` 是线性 ``P_{\mathrm{LOS}}/P_{\mathrm{scatter}}``，不是 dB。直射分量是正实数 ``\sqrt{K/(K+1)}``，``\mathbb{E}[|h|^2] = 1``。`K = 0` 与 `rayleigh_channel` 的系数同分布、同样的 `randn` 顺序；`K = Inf` 使 ``h = 1``，且不消耗随机数。没有时间相关。`rice_channel` 在衰落之后调用 `awgn`，`snr_db` 是平均 ``E_b/N_0``。
- Nakagami-m 是逐符号独立的。`m ≥ 1/2`，``\mathbb{E}[R^2] = 1``。`m = 1` 的复包络与 `rayleigh_channel` 同一套 `randn` 顺序。`2m` 在 `Float64` 中为整数时，幅度平方是 `2m` 个标准正态平方的平均（`m = 1` 不走这条，以免打乱与瑞利对齐的随机数）；其它有限 `m` 用 Marsaglia–Tsang gamma 拒绝抽样，接受与否决定消耗多少个 `randn` / `rand`。`m ≠ 1` 的相位是独立的 ``\mathrm{Unif}[0, 2\pi)``，不是高斯 I/Q。`m = Inf` 幅度恒为 1、相位均匀，理想 CSI 抵消相位之后是 AWGN；这不是 Rice `K = Inf` 的 ``h = 1``。没有多普勒。
- MRC 是接收分集：同一个符号进入 L 条独立支路，每条支路的平均 ``E_b/N_0`` 就是调用时的 `snr_db`，符号能量不在支路之间拆分。`mrc` 假定 ``h`` 已知且就是生成 ``y`` 的那组系数，返回按 ``\sum|h|^2`` 归一化的样本。硬判决可以直接解这条样本；软判决不能再喂一个标量 ``\sigma^2``。全零列得到 `NaN`。不是 Alamouti，不在 `MIMO` 的检测器里，没有信道估计误差，也没有天线相关。Rice 与 MRC 可以叠在 `mrc_channel(...; K)` 上，但测试里的闭式只覆盖瑞利 MRC 和单支路 Rice BPSK。
- LDPC 只覆盖短的教学构造。Gallager 的列置换由本实现的 SplitMix64 规则确定，不是 Gallager 1962 系综里的某一个唯一矩阵，也不保证围长；没有 PEG，没有非正则优化，也没有按密度进化重估 ``\alpha``、``\beta``（固定为 0.75 与 0.15）。编码是构造时的一次稠密高斯消元，边长不超过 4096，不是 Richardson–Urbanke 的线性时间编码器。译码是泛洪调度，不是分层调度。有环时 sum-product 不等于 MAP。``[7,4]`` Hamming 码含 4 环，全零字上第 7 位（列重 3）的单比特翻转会收敛到另一个码字。信道硬判决已经是码字且早停打开时，`posterior` 就是信道 LLR，迭代次数为 0。冗余校验行留在 Tanner 图里。没有 NR BG1/BG2、LDPC 速率匹配、LDPC 打孔、缩短，也没有把解调器外信息再送回译码器。CRC 和 Chase 合并在 `Link` 里，不在 LDPC 译码器内部。
- 链路闭环是单链路 AWGN。MCS 表只有五档，频谱效率把 16 位 CRC 和卷积尾比特算作开销，所以低于母码率乘 ``\log_2 M``。门限是这条表、这个块长上的教学开关，不是 NR CQI，也不会随信道重新估计。3 dB 和 11 dB 附近两条吞吐曲线靠得近，换随机数流可能把交点挪约 1 dB。Chase 合并重复同一编码块；闭环里没有增量冗余，卷积删余也没有接到 Chase 上。重传次数有限，失败的块丢弃。外环只移动一个 SNR 偏移，不控制功率，也没有 HARQ 进程号或时序。二阶矩 ``E_s/N_0`` 估计在短块上有偏。CRC-16 检不出 ``g(x)`` 的倍式。没有多用户，也没有调度器。
- 软解调假定星座等先验、符号独立。LLR 是信道比特边际，没有把译码器外信息反馈进解调器。`approx = :maxlog` 只保留每个比特子集里的最近星座点。衰落系数要由调用者先补偿；本函数不读取 CSI。
- OFDM 的 “MMSE” 是逐导频幅度收缩后再做线性插值，不是 Edfors 的 LMMSE。多径上线性插值有与信噪比无关的误差平底。`ofdm_sync` 校正频偏，但不按定时估计平移缓存。`cfo_estimate` 只在 ``(-1/2, 1/2)`` 个子载波间隔内无模糊：``0.7`` 会报成接近 ``-0.3``，函数不估计整数倍。这是 van de Beek 的循环前缀相关，0.8.4 没有改它的公式或返回值。
- `pulse_shape` 仍把卷积截成 `length(symbols) * sps`，滤波器尾部被丢掉。`matched_filter` 是 `conj.(reverse(taps))` 的全长线性卷积：不截断、不搜索定时、不做分数延迟内插。`downsample_symbols` 的 `offset` 是调用者给出的前导样本数（抽样相位），不按 `sps` 取模。根升余弦在符号点上只是截断后的近似奈奎斯特；突发末尾大约半个滤波器跨度的符号缺少被丢掉的尾部。抽样相位要对上自相关峰值时，`offset = length(taps) - 1`（符号 1 放在上采样后的下标 1，且 `length(symbols) * sps >= length(taps)`）。`Pulse` 自己仍然不估计定时或载波。
- 同步不是生产级。Gardner 要求偶数 `sps ≥ 2`，中点必须落在样点上；奇数过采样直接拒绝，而不是内插到半符号。平均误差在采样偏晚时为正。半个符号附近还有一个负斜率的零点，开环只留正斜率的那个。没有符号钟频率状态，捕获的是一个符号周期里的抽样相位，不是整数个符号的延迟；`single_carrier_timing` 的 `symbol_index` 用的是 `length(taps)`。连续相同的符号使 `late - early ≈ 0`，开环会因为找不到过零而报错。分数延迟是两点线性内插，不是带限 Farrow。环路增益按单位能量根升余弦上的 BPSK 来试，QPSK 的 S 曲线更平，同一个 `alpha` 收得更慢。
- Costas 只定义了 `BPSK()` 和 `QPSK()`。`PSK(4)` 在坐标轴上，鉴相器会拒绝。BPSK 对 ``\pi`` 旋转不敏感（比特全翻，误差仍可以是 0），QPSK 对 ``\pi/2`` 不敏感。线性区大约是 BPSK 的 ``|\theta| < \pi/2``、QPSK 的 ``|\theta| < \pi/4``。默认 ``\alpha``、``\beta`` 用来锁大约 ``10^{-2}`` 弧度/符号的频偏；更大的频偏不在这组增益的捕获范围里。没有锁定检测，没有 AGC，也不跟踪定时。输入是每符号一个样点。
- Schmidl–Cox 必须用 `schmidl_preamble` 那种两半相同的前导，`ofdm_tx` 不会加。小数 CFO 的无模糊范围是 ``(-1, 1)``，再加 ``2z`` 的偶数倍整数它看不见。整数搜索默认只到 ``|2z| \le 8``，而且是平坦信道上的相干相关；FFT 窗起点不对、或者小数部分还没消掉，峰值就会偏。循环前缀把定时度量变成平台。无噪声时最早的最大值就是前导的第一个样点，`schmidl_sync` 把 FFT 窗放在 `d + n_cp`。噪声若把最大值挪到平台后部，这个窗会滑出符号。代码没有做平台平均。`schmidl_sync` 和 `ofdm_sync` 一样不按定时估计裁剪波形。
- FIR 抽头是确定性的，不是随机抽头过程，也不是 Jakes 谱，并且没有接到 OFDM 的 `TDLChannel` 上。`fir_channel` 的噪声按调用者给出的 ``E_b`` 标定（位置参数假定 ``E_s = 1``），不按滤波后的样本功率重标。返回的 `noise_var` 是加在**这条卷积结果**上的噪声方差：实序列为 ``N_0/2``，复序列为 ``N_0``。它只在均衡器的输入就是这个 `rx` 时交给 `mmse_equalize` 或 `mmse_dfe`；匹配滤波之后噪声被着色，这个标量不再是抽样点上的方差。`+Inf` 不加噪声，也不消耗全局随机数。ZF / MMSE 要的是抽样之后的符号率冲激响应。中间有脉冲成形时，它不等于样本率上的 FIR 抽头，除非该脉冲在抽样相位上是冲激。
- LMS 是符号率、监督训练的复基带算法；实信号上共轭退化为恒等。输出是 ``w^{\mathsf{T}} u``（最新样本在 ``u`` 的第一项），更新是 ``w \leftarrow w + \mu e u^*``。步长由调用者选择。稳定范围大约是 ``0 < \mu < 2/\lambda_{\max}``，函数不估计这个上界，也不做归一化 LMS。抽头初值默认是 0。`delay` 是期望符号相对最新接收样本的滞后；能量不在第一抽头上时，`delay = 0` 对不齐。训练序列用完之后抽头冻结，没有判决引导，也没有盲模式。返回的训练段输出用的是该步更新前的抽头，不是用最终抽头把整块重滤一遍。稳态误调不作为返回值，也不对照 Widrow / Proakis §10.1-4 的超额 MSE 公式。无噪声且抽头和时延足够时，误差可以收到舍入；有噪声时更大的稳定步长留下更大的超额误差。不是 RLS，不是 Kalman，不是分数间隔均衡。
- 块 DFE 和 `zf_equalizer` / `mmse_equalizer` 用同一条卷积矩阵。前馈来自列反转后的 QR；MMSE 把 ``\sqrt{\sigma^2} I`` 接在矩阵下面。反馈是上三角因子里严格在游标之后的项，检测顺序从第一个符号到最后一个，不按 V-BLAST 重排。每个符号一行，因为块是有限长的；不是时不变的短抽头 streaming DFE。反馈送进的是 `alphabet` 上的最近邻，平局取更靠前的点。`training` 覆盖到的符号改用已知符号，这一段不经过切片，整块都给训练就是理想反馈。没被训练盖住的位置一旦判错，错误乘上该行反馈抽头加进后面的软输出。小抽头不一定翻得过判决边界，大抽头可以连着翻一段；函数不限制这个突发有多长，也没有 Tomlinson–Harashima 预编码。无噪声且每个反馈判决都正确时，ZF-DFE 把非全零 FIR 恢复到线性代数舍入。MMSE 软输出带有和 `mmse_equalize` 相同的偏差，不另做去偏。`noise_var` 必须有限。
- EVM 是调用者给出的参考符号和接收符号之间的 RMS 幅度比。默认百分数：``100\sqrt{\sum|e|^2/\sum|s|^2}``。dB 是 ``20\log_{10}`` 这个比，10% 为 ``-20\,\mathrm{dB}``。不做幅度、相位或直流拟合，也不换成星座模板功率。
- 扰码是加性 LFSR：同一多项式、同一种子再异或一次即解扰。每次调用从种子重起，没有跨块保留的相位。次数不超过 62。经典多项式是三项式；次数 7 是 ``x^7+x^6+1``，不是 802.11 的 ``x^7+x^4+1``。次数 9 及以上取 ITU-T O.150 的连接多项式；次数 20 是 ``x^{20}+x^3+1``，不是 ``x^{20}+x^{17}+1``。序列不按 O.150 取反，也不塞入额外的 1。不是自同步乘性扰码，不是 Gold 码，没有接到 MCS / HARQ。
- ML 检测是穷举。`svd_precoder` 只返回右奇异向量，不分配功率。`waterfill` / `svd_waterfill` 是单用户总功率注水（Tse 与 Viswanath §7.1.1）：低信噪比集中到最强特征模式，高信噪比在正奇异值上接近均匀。`noise_var` 是 ``\mathbb{E}[|n_i|^2]``，不是 ``N_0/2``。功率为 0 的模式保留成零列。不是多用户迭代注水，没有速率约束。ZF 预编码矩阵未做功率归一化。
- 导频和抽头相位用固定线性同余发生器，不走全局随机数流。
