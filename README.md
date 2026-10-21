# CommToolkit.jl

基带通信链路仿真库（v0.3.0），用来把蒙特卡洛误码率和教材中的闭式公式对上。作者刘恒宇（longee142857）。

本包**不是** 3GPP NR / LTE 物理层，**不是** O-RAN 前传实现，也**不是**生产级协议栈。OFDM 是教科书循环前缀模型，没有资源栅格、DMRS 或空口时隙结构。

## 功能

公开接口从子模块扁平再导出，原有调用方式保持不变（例如 `modulate(bits, BPSK())`、`awgn(syms, snr_db, bits_per_symbol)`）。没有重命名，因此没有弃用别名。

- **调制**（`Modulation`）：BPSK、Gray QPSK、方形 Gray M-QAM（`M ≥ 4` 且为平方、2 的幂）。单位符号能量。硬判决解调；AWGN 下逐比特 LLR（精确边际与 max-log），噪声方差与 `awgn` 的 ``\sigma^2 = N_0/2`` 相同。星座点与每符号比特数。
- **信道**（`Channel`）：AWGN，`snr_db` 是 ``E_b/N_0``（dB）；平坦瑞利衰落，可选理想信道补偿。
- **编码**（`Coding`）：码率 1/2、约束长度 ``K = 3``、八进制生成多项式 (7, 5) 的卷积码，尾比特终止。硬判决 Viterbi、以 LLR 为分支度量的软输入硬输出 Viterbi、终止码上的 log-BCJR 外信息，以及网格距离谱。
- **均衡**（`Equalization`）：频率选择性信道上的线性 ZF / MMSE。
- **OFDM**（`OFDM`）：循环前缀 OFDM 收发、抽头延迟线、已知 CSI / LS / 逐导频收缩后线性插值、van de Beek 符号定时与小数倍频偏。
- **MIMO**（`MIMO`）：i.i.d. 瑞利信道，ML、ZF、MMSE、按行范数排序的 SIC，ZF 与 SVD 预编码，Alamouti ``2\times1`` 与 ``2\times2``。
- **脉冲**（`Visualization`）：升余弦与根升余弦。星座图、眼图、脉冲频谱由包扩展 `CommToolkitPlotsExt` 提供，只有 `using Plots` 之后才加载。
- **指标**（`Metrics`）：误码率、dB 与线性信噪比换算。

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

![BPSK、QPSK、16-QAM 在 AWGN 下以及 BPSK 在平坦瑞利下的蒙特卡洛误码率与闭式曲线](docs/src/assets/ber_validation.png)

卷积码在测试里与硬判决联合界比较（Proakis 与 Salehi §8.2；Lin 与 Costello，《Error Control Coding》第 2 版 §12.3）。同一 ``E_b/N_0`` 下，软输入 Viterbi 与 BCJR 硬判决的 Clopper–Pearson 区间低于硬判决 Viterbi，也低于未编码 BPSK，并落在软判决联合界 ``\sum c_d Q(\sqrt{2 d R E_b/N_0})`` 的同一数量级（Proakis 与 Salehi §8.2；Bahl、Cocke、Jelinek 与 Raviv，IEEE Trans. Inf. Theory 20(2):284–287，1974）。BPSK / Gray QPSK 的 LLR 与闭式条件高斯分布对照（Proakis 与 Salehi §4.2、§5.2）。OFDM / MIMO 测试覆盖无噪声完美恢复，以及高低信噪比下互不重叠的 Clopper–Pearson 区间。

本地运行：

```julia
using Pkg
Pkg.test("CommToolkit")
```

## 源码布局

```
src/CommToolkit.jl    模块入口与再导出
src/modulation.jl     Modulation
src/channel.jl        Channel
src/coding.jl         Coding
src/equalization.jl   Equalization
src/ofdm.jl           OFDM
src/mimo.jl           MIMO
src/visualization.jl  Visualization（作图函数在未加载 Plots 时抛错）
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

- 卷积码仍是未删余的 (7, 5)、``K = 3``。软输入 Viterbi 输出硬比特，并且和硬判决译码一样不把最终状态钉在 0。BCJR 按终止码把尾比特固定为 0、最终状态固定在 0，返回信息比特外信息（后验 LLR 减去先验）；码是非系统的，没有信道系统位可减。没有 SOVA、交织、删余、Turbo / LDPC 迭代。联合界按未终止码的第一事件界计算。
- 软解调假定星座等先验、符号独立。LLR 是信道比特边际，没有把译码器外信息反馈进解调器。`approx = :maxlog` 只保留每个比特子集里的最近星座点。衰落系数要由调用者先补偿；本函数不读取 CSI。
- OFDM 的 “MMSE” 是逐导频幅度收缩后再做线性插值，不是 Edfors 的 LMMSE。多径上线性插值有与信噪比无关的误差平底。`ofdm_sync` 校正频偏，但不按定时估计平移缓存。
- `pulse_shape` 把卷积结果截成 `length(symbols) * sps`，滤波器尾部被丢掉。
- ML 检测是穷举。SVD 预编码不做注水。ZF 预编码矩阵未做功率归一化。
- 导频和抽头相位用固定线性同余发生器，不走全局随机数流。
