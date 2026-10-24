# CommToolkit.jl

基带通信链路仿真库（v0.9.1），用来把蒙特卡洛误码率和教材中的闭式公式对上，并在 AWGN 上跑一条教学级的 MCS / HARQ 闭环。作者刘恒宇（longee142857）。

公开名字是 [稳定面](docs/src/api.md)：`using CommToolkit` 带入的符号，这一版不改名。本包不是什么，见 [边界](docs/src/boundaries.md)。对照表见 [验证矩阵](docs/src/validation.md)。参数级限制见 [已知限制](docs/src/limitations.md)。

这一版不升到 1.0.0。0.9.1 只把上述三页和这份导航固定下来，不改已有函数的位置参数和返回值。

## 功能

公开接口从子模块扁平再导出。例如 `modulate(bits, BPSK())`、`awgn(syms, snr_db, bits_per_symbol)`。没有重命名，因此没有弃用别名。`CommToolkit.modulation` 不在 `export` 列表里。

- **调制**（`Modulation`）：BPSK、Gray QPSK、Gray M-PSK、方形 Gray M-QAM，单位符号能量。硬判决，以及相干 AWGN 下的逐比特 LLR。DBPSK / DQPSK 在块首放参考符号，不是 ``\pi/4``-DQPSK。
- **信道**（`Channel`）：AWGN 的 `snr_db` 是 ``E_b/N_0``。平坦瑞利、块衰落、平坦 Rice（`K` 为线性比）、Clarke/Jakes 时间相关、Nakagami-m、理想 CSI 的 MRC，以及调用者给定抽头的单载波 FIR。
- **编码**（`Coding`）：终止的 (7, 5) 卷积码，硬 / 软 Viterbi、log-BCJR、SOVA，以及到 2/3、3/4 的周期删余。短帧 Turbo（记忆 2 的 RSC）。``\mathrm{GF}(2^m)``（``m = 2,3,4``）、Hamming(7,4)、二进制 BCH(15,7)。短 LDPC：Gallager、Fan 阵列码、系统编码、泛洪 sum-product 与 min-sum。
- **交织**（`Interleaving`）：矩形块交织，按行写入、按列读出。
- **BICM**（`BICM`）：比特交织编码调制的教学接收机。默认 Gray。`bicm_psk8_id_labels` 是 8-PSK 集合分割。没有接到 MCS / HARQ。
- **链路**（`Link`）：CRC-16、五档教学 MCS、Chase 合并、按 ``E_s/N_0`` 切换的 MCS。``E_s = 1``，``\sigma^2 = N_0/2``。
- **均衡**（`Equalization`）：线性块 ZF / MMSE、符号率 LMS、有限长块 DFE。
- **OFDM**（`OFDM`）：循环前缀收发、静态抽头延迟线、已知 CSI / LS / 逐导频收缩后线性插值、van de Beek 定时与小数倍频偏。
- **MIMO**（`MIMO`）：i.i.d. 瑞利上的 ML、ZF、MMSE、SIC，ZF 与 SVD 预编码，总功率注水，Alamouti ``2\times1`` 与 ``2\times2``。
- **扰码**（`Scrambling`）：加性 LFSR。`scramble` 与 `descramble` 从同一颗种子重起。
- **脉冲**（`Pulse`、`Visualization`）：升余弦与根升余弦，匹配滤波，按每符号样本数抽样。星座图、眼图和脉冲频谱在 `using Plots` 之后由 `CommToolkitPlotsExt` 提供。
- **同步**（`Sync`）：Gardner 定时、BPSK/QPSK 的 Costas 环、Schmidl–Cox 前导。不改 `ofdm_sync`。
- **指标**（`Metrics`）：误码率、符号错误率、帧错误率、RMS EVM，以及 dB 与线性信噪比换算。

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

固定种子 `0x0c0de1162` 的蒙特卡洛点与闭式曲线。实线是理论值，标记是仿真。生成脚本是 `demos/validation_figure.jl`。统计检验在 `test/`：双侧 Clopper–Pearson 区间（``\alpha = 10^{-4}``，并要求足够的错误计数），理论值必须落在区间内。

公式、出处、测试文件和演示的对应关系在 [验证矩阵](docs/src/validation.md)。实现在 `test/reference.jl`。下面这张图仍是 v0.6 之前的四条曲线；Gray M-PSK、差分 PSK、Rice、MRC、Jakes 和 Nakagami-m 的对照在同一页。

![BPSK、QPSK、16-QAM 在 AWGN 下以及 BPSK 在平坦瑞利下的蒙特卡洛误码率与闭式曲线](docs/src/assets/ber_validation.png)

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
src/bicm.jl           BICM：比特交织编码调制与迭代解调
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
docs/                 公开接口、验证矩阵、边界、已知限制，以及 Documenter API 参考
demos/                可单独运行的脚本，不在测试范围内
```

API 参考由 Documenter 从导出函数的 docstring 生成，页面只含接口说明。

## 许可

MIT。Copyright (c) 2026 刘恒宇 (longee142857)。见 `LICENSE`。

## 已知限制

完整条目在 [docs/src/limitations.md](docs/src/limitations.md)，与 v0.9.0 的技术边界相同。这里只留检索用的几条：

- 卷积母码仍是 (7, 5)。删余只有 2/3 和 3/4 两张图案。联合界不覆盖删余码。
- Turbo 是记忆 2 的短帧，不是 3GPP。BCH 只有 (15, 7)，不是 Reed–Solomon。
- LDPC 是短教学码，不是 NR 基图。BICM-ID 没有 EXIT 图，也没有接到 MCS / HARQ。
- Jakes 是平坦 Clarke 谱的 DFT 实现，不是正弦和，也不是抽头延迟线。`TDLChannel` 不是 TR 38.901 的 CDL。
- OFDM 的 “MMSE” 是导频收缩后的线性插值。`ofdm_sync` 只在 ``(-1/2, 1/2)`` 个子载波间隔内无模糊。
- Gardner、Costas、Schmidl–Cox 有写明的捕获范围，不是生产级接收机，也不是调度器。
- 链路闭环是单链路 AWGN 上的五档教学 MCS 和有限次 Chase 合并。没有多用户调度。
