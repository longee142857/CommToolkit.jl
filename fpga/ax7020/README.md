# AX7020 教学链路

黑金 ALINX **AX7020**（Zynq-7020，器件 **XC7Z020**，常见封装 CLG400）上的 **PL 教学**示例。和 Julia 包 `CommToolkit` v1.0 里的调制、信道、卷积硬判决、短 OFDM、平坦 MIMO，以及链路闭环里的 CRC / 离散 MCS / Chase **概念对齐**。定点噪声是 12 个离散均匀数的和，不是 `randn`。**不要把这里的误码率或 BLER 当成和 Julia 浮点仿真 bit-exact，也不要当成已经和闭式公式做过区间验收。** 本目录没有附带定点对照的 Julia 脚本；没做那一步，就不声称数值一致。

刀6 把这条 FPGA 教学路线收到 v1 主链路的骨架为止。主链路在这里的意思是：硬判决调制、AWGN、卷积母码的硬 Viterbi、N=16 CP-OFDM、2×2 平坦 ZF/MMSE，再加上 CRC-16、三档 MCS 和硬 Chase 的闭环计数。六刀都可以在 Icarus 里仿真。这不是把 v1 的每条算法搬上 PL。

| 刀 | 内容 | 状态 |
| --- | --- | --- |
| 1 | BPSK + 实 AWGN + BER | 已在 main |
| 2 | QPSK，以及方形 Gray 16-QAM、64-QAM 的硬映射 / 硬判决 | 已在 main |
| 3 | K=3 (7, 5) 卷积母码 + 硬判决 Viterbi，接到 BPSK | 已在 main |
| 4 | N=16 的 CP-OFDM：IFFT/FFT、循环前缀、QPSK 子载波 | 已在 main |
| 5 | 2×2 平坦 MIMO，理想 CSI 下的 ZF 与 MMSE，QPSK | 已在 main |
| 6 | 三档 MCS、CRC-16、硬 Chase、BLER/吞吐计数 | 本目录，路线收口 |
| — | 软判决全链、删余、LDPC、Turbo、Kronecker、38.901、EXIT、路线 C 等 | 刻意不做，见文末对照表 |

刀1 的链路是：PRBS 比特 → BPSK 映射 → 实 AWGN → 硬判决 → 差错计数。刀2 在旁边加复星座，噪声源和差错计数的做法不变。刀3 在 BPSK 前面加上 (7, 5) 编码，后面加上硬 Viterbi，误码改数信息比特。刀4 是 16 点 CP-OFDM：QPSK 子载波、IFFT、循环前缀、可选时域噪声、去前缀、FFT、硬判决。刀5 是固定的 2×2 平坦信道，理想 CSI，ZF 和 MMSE 硬判决。刀6 把已有模块收成一条短帧闭环：16 比特净荷、CRC-16、按 Es/N0 选三档里的一档、硬 Chase 重传，直到 CRC 通过或达到 `max_tx`。

## 边界

- 教学级定点，不是量产调制解调器。没有 AGC、成形滤波、定时或载波恢复。
- 不是 NR，不是 LTE，不是 O-RAN。没有资源栅格、前传、eCPRI 或 xApp。
- 和 Julia `CommToolkit` 是概念对照，不是数值契约。数据 LFSR 的移位规则可以和 `lfsr_sequence` 对上（见下表）；噪声和 BER 不对上。
- 刀1 是 **BPSK + 实 AWGN + BER**。刀2 是 **QPSK** 和方形 Gray **16-QAM / 64-QAM** 硬判决。刀3 是终止的 **K=3 (7, 5) 卷积母码**和**硬判决 Viterbi**。刀4 是 **N=16、CP=4** 的 CP-OFDM 骨架。刀5 是 **2×2 平坦 MIMO** 的 ZF 与 MMSE。刀6 是 **三档 MCS + CRC-16 + 硬 Chase**，统计 ACK/失败和净荷比特 / 信道使用。刻意不做的部分列在文末「和 v1.0 的对照」：软判决全链、删余、LDPC、Turbo、Kronecker、TR 38.901、EXIT、路线 C 的多用户 xApp，以及 NR 调度和增量冗余。
- 引脚是 AX7020 手册里的常见占位。**没有在这块板上实测过**，综合和时序也没有在本仓库里跑过。

GitHub Actions 仍只测 Julia 包（`src/`、`ext/`）。这里不进那个覆盖率，CI 也不跑 Vivado 或 Icarus。

## 和 Julia 怎么对

| 这里 | Julia | 对齐到哪 |
| --- | --- | --- |
| `bpsk_mapper`：`0 → +AMP`，`1 → −AMP` | `modulate(bits, BPSK())`，即 `1 - 2*bits` | 极性相同。Julia 是单位能量浮点 ±1；这里是整数 ±1024。 |
| `bpsk_slicer`：样本 `< 0` 判 1，否则判 0 | `demodulate` 对 BPSK 用 `syms .< 0` | 判决规则相同，包括 sample = 0 判成 0。 |
| `awgn_clt` 的刻度 | `awgn(syms, snr_db, 1)` | 都把 `snr_db` 当成 \(E_b/N_0\)（dB）。实噪声方差的目标是 \(\sigma^2 = E_b / (2 E_b/N_0)\)。Julia 用 `randn`；这里用下面的整数近似。 |
| `n_err / n_bit` | `ber(tx, rx) = k/n` | 定义相同：逐比特比较，计差错占比。 |
| 数据 `lfsr_fib`，\(x^{15}+x^{14}+1\)，种子全 1 | `lfsr_sequence(prbs_polynomial(15); seed=nothing)` | **比特序列的移位规则相同**。测试台锁了前 32 位。这只说明扰码那一节的 Fibonacci 约定。 |
| 噪声用的 \(x^{31}+x^{28}+1\) | `prbs_polynomial(31)` 的连接多项式 | 只借用多项式。12 路并行、每次取 8 位，不是 `lfsr_sequence` 的逐比特输出。 |
| `qpsk_mapper` / `qpsk_slicer` | `QPSK()`，`modulate` / `demodulate` | 比特顺序和符号判决相同：b0 看实部，b1 看虚部，负数判 1。幅度是整数 724，不是 \((\pm 1 \pm j)/\sqrt{2}\)。不是 `QAM(4)`，也不是 `PSK(4)`。 |
| `qam_mapper` / `qam_slicer` | `QAM(16)`、`QAM(64)` | MSB 先入。前一半比特是 I 路 Gray，后一半是 Q 路。电平取最近整数单位。门限上的平局归右格，和 `findmin` 的平局规则不一定相同。 |
| `conv75_encode` | `conv_encode` | 生成多项式 \(1+D+D^2\)、\(1+D^2\)，输出顺序 \((g_1,g_2)\)，再补两个 0 尾比特。测试台锁了固定向量。 |
| `viterbi_hard75` | 整数 `viterbi_decode` | 汉明支路度量，严格小于才替换，从最小度量状态回溯，平局取最小状态号。去掉尾比特。无噪声和给定硬比特可以对上。不是 LLR 那条，没有删除符号。 |
| `fft16` 的 IFFT / FFT | `FFTW.ifft` / `fft` | IFFT 每级右移 1 位，四级合起来除以 16。FFT 不除以 N。旋转因子是 1024 倍的余弦，乘完再右移 10 位。无噪声时 `fft(ifft(X))` 回到 X 附近几位，不是逐位相同。 |
| CP 长度 4，加在 16 点前面 | `generate_ofdm_symbol` 的前缀 | 前缀是时域块的最后 4 个样本。接收丢掉前 4 个再做 FFT。长度不是 Julia 默认的 64/16，只是同一个 1/4。 |
| 16 个 QPSK 子载波，含 DC | `OFDMParams` 的数据/导频/保护 | **不是**那张默认图。没有保护带，没有导频，没有把 DC 置零。 |
| `mimo2x2_detect` | `zf_detect` / `mmse_detect` | 2×2，\(G=H^{\mathrm{H}}H+\sigma^2 I\)，再 \(G^{-1}H^{\mathrm{H}}y\)。`noise_var=0` 是 ZF。除法朝 0 截断。不是 `pinv` 的逐位结果。 |
| 链路里的固定 H | `mimo_channel` | **不是** i.i.d. 瑞利抽样，也不对 H 做平均。默认是实矩阵 `[[1024,256],[128,1024]]`，1024 表示 1。 |
| `crc16_ccitt` | `crc16_ccitt` / `crc_append` / `crc_check` | 多项式 0x1021，初值 0，高位先入，不反射，xor-out 0。`"123456789"` 的余数是 `0x31C3`。整字余数 0 才算通过。 |
| `link_mcs_harq` 的三档 | `teaching_mcs_table` | 只借用「离散表 + 按 Es/N0 选一档」。档数、净荷长度和门限都不是那张五档表，也不是 `teaching_amc_thresholds()`。 |
| 硬 Chase：比特 0 记 +1，比特 1 记 −1，和为负判 1 | `simulate_link` 的 Chase | 都是重传同一份符号，CRC 过了就停，不是增量冗余。Julia 把软比特 LLR 相加；这里是硬判决的 ±1 累计。平局（和为 0）判成比特 0。 |

想自己做定点对照时，按下一节的 AMP、均值 1530、算术右移（负数朝 \(-\infty\) 截断）和 README 里的种子，在 Julia 里重写**整数**模型，再比 `n_err`。浮点 `awgn` 对不上这些比特，这是预期的。

## 定点模型

符号幅度 `AMP = 1024`（有符号 16 位）。BPSK 的 \(E_s = E_b\)，判决门限在 0。若噪声真是高斯、标准差为

\[
\sigma = \frac{\mathrm{AMP}}{\sqrt{2 \cdot 10^{\gamma/10}}}, \quad \gamma = E_b/N_0\ \text{(dB)},
\]

则差错概率是 \(Q(\sqrt{2 \cdot 10^{\gamma/10}})\)，和教材里单位能量 BPSK 同一个自变量。刻度表用的是 `round(σ)`，只为这个 AMP 算过。`awgn_clt` 的参数不是 1024 时，仿真和综合都会因为实例化一个不存在的模块而失败，避免静默用错信噪比。

噪声不是高斯采样：

1. 12 路 31 级 Fibonacci LFSR，多项式 \(x^{31}+x^{28}+1\)，种子各不相同（写在 `rtl/awgn_clt.v` 的 `seed_of`）。
2. 每路取状态低 8 位，当作 \(\{0,\ldots,255\}\) 上的均匀整数。
3. 12 个数相加，减去均值 \(12 \times 127.5 = 1530\)。和的标准差是 \(\sqrt{65535} \approx 255.998\)。
4. `noise = floor(centered * MULT / 256)`。`MULT = round(σ)`，右移 8 位正好把标准差从约 256 收到 σ。负数用算术右移，朝 \(-\infty\) 截断，不四舍五入。

| `snr_index` | \(E_b/N_0\) | MULT |
| --- | ---: | ---: |
| 0 | 0 dB | 724 |
| 1 | 2 dB | 575 |
| 2 | 4 dB | 457 |
| 3 | 6 dB | 363 |
| 4 | 8 dB | 288 |
| 5 | 10 dB | 229 |
| 6、7 | 无噪声 | 0 |

`bits_per_symbol = 1` 用上表。QPSK（k = 2）、16-QAM（k = 4）、64-QAM（k = 6）在 `awgn_clt` 里另有刻度，公式是 \(\sigma = 1024 / \sqrt{2 k \cdot 10^{\gamma/10}}\)。刀1 的 BPSK 链路把 k 接成 1、种子组保持 0，上面的 BPSK 计数不因此改动。

和的取值被限制在 \([-1530, 1530]\)，大约 \(\pm 6\sigma\) 之后没有尾。0 dB 到 8 dB 的判决门限还落在这个范围里，所以能看到差错随信噪比下降。10 dB 的 BPSK 理论误码率已经在 \(10^{-6}\) 量级，短帧里经常是 0，测试台只要求它不超过 8 dB 的差错数。

最大 `|centered| * 724 / 256` 约 4328，加上 AMP 1024 仍在有符号 16 位以内。不要把 AMP 加大到和噪声相加会翻符号位。

数据比特是另一颗 LFSR：\(x^{15}+x^{14}+1\)，种子 `15'h7FFF`（全 1）。每个时钟一个符号。

## 目录

```
fpga/ax7020/
  README.md                 本说明
  rtl/lfsr_fib.v            Fibonacci LFSR
  rtl/bpsk_mapper.v         BPSK 映射
  rtl/bpsk_slicer.v         硬判决
  rtl/awgn_clt.v            12 路 LFSR 的均匀和噪声
  rtl/bpsk_awgn_ber.v       刀1 链路：计数 n_err、n_bit
  rtl/ax7020_bpsk_top.v     板级薄封装，仍然只接 BPSK
  rtl/qpsk_mapper.v         QPSK 硬映射
  rtl/qpsk_slicer.v         QPSK 硬判决
  rtl/qam_gray.v            16/64-QAM 硬映射和硬判决
  rtl/qpsk_qam_awgn_ber.v   刀2 链路
  rtl/conv75_trellis.vh     (7, 5) 的下一状态和输出，被编码器和译码器 include
  rtl/conv75_encode.v       卷积编码
  rtl/viterbi_hard75.v      硬判决 Viterbi
  rtl/conv_bpsk_awgn_ber.v  刀3 链路：编码、BPSK、AWGN、译码、BER
  rtl/fft16.v               16 点基-2 FFT / IFFT
  rtl/ofdm16_qpsk_ber.v     刀4 链路：QPSK、IFFT、CP、FFT、BER
  rtl/mimo2x2_detect.v      2×2 ZF/MMSE
  rtl/mimo2x2_qpsk_ber.v    刀5 链路：固定 H、QPSK、AWGN、BER
  rtl/crc16_ccitt.v         CRC-16/XMODEM，按比特移位
  rtl/link_mcs_harq.v       刀6：三档 MCS、硬 Chase、BLER/吞吐计数
  tb/tb_bpsk_unit.v         刀1：映射、判决、LFSR、噪声刻度
  tb/tb_bpsk_awgn_ber.v     刀1：BER 对 Eb/N0
  tb/tb_ax7020_top.v        顶层复位和按键的短冒烟
  tb/tb_qpsk_qam_unit.v     刀2：星座坐标、往返、噪声刻度
  tb/tb_qpsk_qam_ber.v      刀2：BER 对 Eb/N0
  tb/tb_conv75_unit.v       刀3：编码和硬 Viterbi 向量
  tb/tb_conv_bpsk_ber.v     刀3：编码 BPSK 的信息比特 BER
  tb/tb_fft16.v             刀4：直流冲激和 QPSK 的 IFFT/FFT 往返
  tb/tb_ofdm16_ber.v        刀4：CP-OFDM 的 BER
  tb/tb_mimo2x2_unit.v      刀5：ZF/MMSE 的定点向量
  tb/tb_mimo2x2_ber.v       刀5：固定信道上的 BER
  tb/tb_crc16.v             刀6：0x31C3、追加后余数 0、单比特翻转
  tb/tb_link_mcs.v          刀6：无噪声三档、自动选档、硬 Chase
  constr/ax7020_bpsk.xdc    占位引脚，需按板本改。刀2 及以后没有新管脚
  sim/run_iverilog.sh       本地 Icarus 仿真
```

源文件是 UTF-8。标识符都是 ASCII，中文只在注释里。

## 仿真

工具是 [Icarus Verilog](http://iverilog.icarus.com/) 的 `iverilog` 和 `vvp`。不需要 Vivado，也不需要板子。版本不锁；脚本按 Verilog-2012（`-g2012`）编译。下面的计数表是 **Icarus Verilog 12.0** 的一次输出。

```bash
# Debian / Ubuntu 尚未安装时
sudo apt-get install iverilog

cd fpga/ax7020
bash sim/run_iverilog.sh
```

脚本依次跑十三个测试台（刀1 三个，刀2 两个，刀3 两个，刀4 两个，刀5 两个，刀6 两个）。编译带 `-I rtl`，因为编码器和译码器 include 同目录的 `conv75_trellis.vh`。每个测试台打印一行 `PASS`。刀1 的 `NBITS = 16384`。刀2 的 `NBITS = 16380`。刀3 每帧 64 个信息比特、32 帧。刀4 是 32 个 OFDM 符号、共 1024 比特。刀5 是 256 个 2×2 向量、共 1024 比特。刀6 默认 8 个传输块，每块 16 个净荷比特。种子写在 RTL 里。同一版本的 Icarus 应得到同一张表。

通过条件：

- 无噪声（`snr_index = 7`）差错为 0，且数据 LFSR 前 32 位与 Fibonacci 规则一致。
- 0 dB 到 8 dB 差错数严格下降。
- 10 dB 的差错数不超过 8 dB。
- 刀1：0 dB 的差错比例落在 5% 到 12%（教材 \(Q\) 约 7.9%）。4 dB 落在 0.4% 到 4%（教材约 1.25%）。
- 刀2：QPSK 的 0 dB 落在 4% 到 15%，4 dB 落在 0.4% 到 4%。16-QAM 的 0 dB 落在 6% 到 25%。64-QAM 的 0 dB 落在 8% 到 35%。三种都要求 0 dB 到 8 dB 严格下降。
- 刀3：无噪声信息比特差错为 0。0 dB、2 dB、4 dB 的差错严格下降，6 dB 和 8 dB 不超过前一档。0 dB 落在 2% 到 40%。这不是联合界验收。
- 刀4：无噪声差错为 0。0 dB 到 4 dB 严格下降，6 dB 和 8 dB 不超过前一档。0 dB 落在 2% 到 20%。不是 `apply_channel` 的样本 SNR 验收。
- 刀5：ZF 和 MMSE 无噪声都是 0。各自 0 dB 到 4 dB 严格下降，6 dB 和 8 dB 不超过前一档。0 dB 落在 2% 到 25%。不要求 MMSE 每一档都低于 ZF。
- 刀6：无噪声三档都是 8 个 ACK、0 个失败，信道使用分别是 256、544、128。`mcs_force = 3` 时，`snr_index` 0 和 1 选档 0，2 和 3 选档 1，4 以及无噪声的 6、7 选档 2。0 dB 上 `max_tx = 3` 的失败数必须少于 `max_tx = 1`，并且信道使用更多。这是这一套种子上的教学计数，不是 `simulate_link` 的 Clopper–Pearson 验收。

这是为了抓住「噪声没加上」或「刻度错一个数量级」。**不是** Julia `test/` 里那种 Clopper–Pearson 对闭式的验收，也不是和 `demos/ber_bpsk.jl` 比比特。

波形（可选，文件在 `sim/out/`，已被 gitignore）：

```bash
vvp sim/out/tb_bpsk_awgn_ber.vvp +vcd
gtkwave sim/out/bpsk_awgn_ber.vcd
```

`sim/out/` 里的 `.vvp`、`.log`、`.vcd` 不要提交。

刀1（BPSK）这张表是当前 RTL、上述种子、`NBITS = 16384`、Icarus Verilog 12.0 的一次输出。给 `awgn_clt` 加上 k 之后，k = 1 的这列数字没有变。换种子、改 `NBITS` 或换一套舍入，数字会变。闭式列是教材参照，**不是**这张表的验收门限，也不是 Julia 蒙特卡洛的结果。这一次从 0 dB 到 8 dB 单调下降，并且和 \(Q\) 同量级；10 dB 的期望差错大约 0.06 个，这一帧是 0。

| \(E_b/N_0\) | 差错 / 16384 | 这次的 BER | 教材 \(Q(\sqrt{2 E_b/N_0})\) |
| --- | ---: | ---: | ---: |
| 无噪声 | 0 | 0 | 0 |
| 0 dB | 1289 | \(7.87 \times 10^{-2}\) | \(7.86 \times 10^{-2}\) |
| 2 dB | 592 | \(3.61 \times 10^{-2}\) | \(3.75 \times 10^{-2}\) |
| 4 dB | 211 | \(1.29 \times 10^{-2}\) | \(1.25 \times 10^{-2}\) |
| 6 dB | 37 | \(2.26 \times 10^{-3}\) | \(2.39 \times 10^{-3}\) |
| 8 dB | 3 | \(1.83 \times 10^{-4}\) | \(1.91 \times 10^{-4}\) |
| 10 dB | 0 | 0 | \(3.87 \times 10^{-6}\) |

## 刀2：QPSK 与方形 Gray QAM

每个时钟一个复符号。数据 LFSR 仍是 \(x^{15}+x^{14}+1\)、全 1 种子；一次取出 k 个比特，最早的比特是标号的 MSB。I、Q 各接一个 `awgn_clt`：I 用种子组 0，Q 用种子组 1，k 分别是 2、4、6。差错按比特计，和 `ber` 一样是 \(k_{\mathrm{err}}/n\)。

整数格点（单位能量的近似，不是浮点坐标）：

| 制式 | 每维电平 | 说明 |
| --- | --- | --- |
| QPSK | \(\pm 724\) | \(724 = \mathrm{round}(1024/\sqrt{2})\)。00 在第一象限。这是 `QPSK()`，不是 `QAM(4)`，也不是轴上的 `PSK(4)`。 |
| 16-QAM | \(\pm 324,\ \pm 972\) | \(324 = \mathrm{round}(1024/\sqrt{10})\)。Gray：00→−3、01→−1、11→+1、10→+3，再乘 324。 |
| 64-QAM | \(\pm 158\) 到 \(\pm 1106\)，步长 316 | \(158 = \mathrm{round}(1024/\sqrt{42})\)。3 比特二进制反射 Gray。 |

硬判决看每维和门限的位置。门限在相邻电平正中；样本正好等于门限时用 `<`，归到右边那一格。Julia `demodulate` 对 QAM 用 `findmin`，距离完全相等时取星座表里更小的下标。只在这个平局上可能不一致。无噪声的星座点不在门限上，测试台要求无噪声差错为 0。

QPSK 的教材比特差错率和 BPSK 相同，都是 \(Q(\sqrt{2 E_b/N_0})\)，因为每维是幅度 \(1/\sqrt{2}\) 的 BPSK，而 \(E_b = E_s/2\)。16-QAM / 64-QAM 的精确平均比特差错率是 Cho and Yoon（2002）那条式子，Julia 的 `QAM` 文档写了高信噪比近似。下表不拿它们做门限。

Icarus Verilog 12.0、`NBITS = 16380`、当前种子的一次输出。0 dB 到 8 dB 都是严格下降。QPSK 在 10 dB 这一帧是 0，期望差错大约 0.06 个。

| \(E_b/N_0\) | QPSK 差错 | QPSK BER | 16-QAM 差错 | 16-QAM BER | 64-QAM 差错 | 64-QAM BER |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 无噪声 | 0 | 0 | 0 | 0 | 0 | 0 |
| 0 dB | 1283 | \(7.83 \times 10^{-2}\) | 2353 | \(1.44 \times 10^{-1}\) | 3345 | \(2.04 \times 10^{-1}\) |
| 2 dB | 581 | \(3.55 \times 10^{-2}\) | 1656 | \(1.01 \times 10^{-1}\) | 2616 | \(1.60 \times 10^{-1}\) |
| 4 dB | 209 | \(1.28 \times 10^{-2}\) | 981 | \(5.99 \times 10^{-2}\) | 1983 | \(1.21 \times 10^{-1}\) |
| 6 dB | 44 | \(2.69 \times 10^{-3}\) | 451 | \(2.75 \times 10^{-2}\) | 1429 | \(8.72 \times 10^{-2}\) |
| 8 dB | 1 | \(6.11 \times 10^{-5}\) | 158 | \(9.65 \times 10^{-3}\) | 900 | \(5.49 \times 10^{-2}\) |
| 10 dB | 0 | 0 | 31 | \(1.89 \times 10^{-3}\) | 424 | \(2.59 \times 10^{-2}\) |

板上的 `ax7020_bpsk_top` 仍然只跑 BPSK。刀2 及以后的模块可以综合，但没有接到按键或 LED，也没有新的 `.xdc`。要上板得自己改顶层。本仓库没有为它们跑过 Vivado，也没有比特流。

## 刀3：(7, 5) 卷积码和硬判决 Viterbi

母码对齐 Julia 的 `conv_encode`：约束长度 \(K=3\)，生成多项式 \(g_1=1+D+D^2\)、\(g_2=1+D^2\)（八进制 7 和 5），每步先出 \(g_1\) 再出 \(g_2\)。信息比特之后补两个 0，编码器回到状态 0。块上的真实码率是 \(n/(2n+4)\)，不是恰好 \(1/2\)。

硬 Viterbi 对齐整数版 `viterbi_decode`：支路度量是两比特汉明距离；候选度量严格更小才替换；回溯从路径度量最小的状态出发，平局取状态号最小的那个；不强制结束在状态 0；输出丢掉两个尾比特。测试台锁了 `10110010` 加尾的编码字、其中 1 个和 2 个编码比特的翻转，以及自由距事件 `100`。这些是同一套网格上的离散规则，所以硬比特可以对上。加上 AWGN 之后的误码计数不是和 Julia `randn` 或联合界的数值契约。

**不做删余。** 2/3、3/4 不在这条 FPGA 路线里。Julia 的硬判决接口没有删除符号，删余是在 LLR 里把丢掉的母比特写成 0，那是软判决路径。这里没有软度量，也没有 LLR。

链路是：LFSR 信息比特 → 编码 → 每个编码比特做 BPSK → 两路独立的均匀和噪声 → 硬判决 → Viterbi → 只对信息比特计 BER。`snr_index` 是信息比特 \(E_b/N_0\)。名义码率 \(1/2\)、编码符号 \(E_s=1\)，所以 \(E_b=2\)，和测试里 `awgn(..., energy_per_bit=2)` 的噪声公式一样：\(\sigma = 1024/\sqrt{10^{\gamma/10}}\)。尾比特造成的速率损失没有再折进 \(\sigma\)。这张刻度是 `awgn_clt` 的 `SCALE_MODE=1`。刀1 的 `k=1` 表不动，未编码 BPSK 和这条链路的噪声不是同一套数，不能把两条 BER 曲线当成同一条 \(\sigma\)。

Icarus Verilog 12.0，每帧 64 个信息比特、32 帧（2048 比特），当前种子的一次输出。0 dB 到 4 dB 下降；6 dB 和 8 dB 这一次是 0。2048 比特太短，不能用来声称已经低于未编码 BPSK 或低于联合界。

| \(E_b/N_0\) | 信息比特差错 / 2048 | BER |
| --- | ---: | ---: |
| 无噪声 | 0 | 0 |
| 0 dB | 428 | \(2.09 \times 10^{-1}\) |
| 2 dB | 173 | \(8.45 \times 10^{-2}\) |
| 4 dB | 17 | \(8.30 \times 10^{-3}\) |
| 6 dB | 0 | 0 |
| 8 dB | 0 | 0 |

## 刀4：N=16 的 CP-OFDM

长度写死为 **N = 16**，循环前缀 **4** 个样本。Julia 默认是 64 和 16，这里只沿用「前缀是 FFT 长度的四分之一」这一点，方便在 Icarus 里把整条链路跑完。旋转因子只有 0、±1、±j 和 \(\pi/8\) 的余弦，用 `round(1024 * cos)` 存成整数。

发射：每个子载波一个 QPSK（刀2 的映射，幅度 724）→ IFFT → 把时域最后 4 点接到前面。16 个格子全是数据，**包括 DC**。没有保护带，没有导频，没有 LCG 导频序列。所以这不是 `OFDMParams` 那张默认子载波图，也不是 `ofdm_tx` 的整段行为。

接收：假定突发从样本 0 对齐，丢掉前缀，做 FFT，对每个格子硬判决。平坦、`H = 1`，没有除法均衡。

和 Julia OFDM 的差距，这里故意留着：

- 没有 `ofdm_sync`，没有 van de Beek 的定时，没有小数倍 CFO，更没有整数倍 CFO。
- 没有 `TDLChannel`，没有多径，因此前缀只在无失真回路里被加上再丢掉，没有用来吸收时延扩展。
- 没有 `ofdm_rx_ls` / `ofdm_rx_mmse`，也没有把已知 `H_freq` 除掉的 `ofdm_rx_known`。理想平坦信道下 FFT 结果直接判决。
- Julia `apply_channel` 的 `snr_db` 是时域样本功率比上噪声功率，不是 \(E_b/N_0\)。刀4 的 `snr_index` 仍是想让每个子载波看到和刀2 QPSK 一样的 \(E_b/N_0\)：时域噪声取 `awgn_clt` 的 k=2 刻度，再算术右移 2 位。未归一化的 16 点 FFT 把标准差乘 4，右移先除掉。定点舍入之后这只是近似，不是那条样本 SNR 的定义。

Icarus Verilog 12.0、32 个符号、1024 比特的一次输出。直流冲激的 IFFT/FFT 在 `tb_fft16` 里是精确的 64 和 1024；QPSK 回路无噪声差错为 0，但格子上的数值和 724 差几个 LSB，所以不是 bit-exact。

| \(E_b/N_0\) | 差错 / 1024 | BER |
| --- | ---: | ---: |
| 无噪声 | 0 | 0 |
| 0 dB | 77 | \(7.52 \times 10^{-2}\) |
| 2 dB | 36 | \(3.52 \times 10^{-2}\) |
| 4 dB | 9 | \(8.79 \times 10^{-3}\) |
| 6 dB | 2 | \(1.95 \times 10^{-3}\) |
| 8 dB | 0 | 0 |

## 刀5：2×2 平坦 MIMO 的 ZF 与 MMSE

两根发射、两根接收，一个时隙一个 QPSK 向量。信道在整段仿真里固定，检测器知道这份 H，也就是理想 CSI。模型是 \(y = Hx + n\)，和 `apply_mimo` 同一件事。默认

\[
H = \begin{bmatrix} 1024 & 256 \\ 128 & 1024 \end{bmatrix},
\]

实数，1024 表示 1。两列不正交，所以有干扰。这不是 `mimo_channel` 的 i.i.d. 瑞利，也不对 H 做遍历平均。换 H 要把参数传进 `mimo2x2_qpsk_ber`。

检测和 Julia 一样：

\[
\hat x = (H^{\mathrm{H}} H + \sigma^2 I)^{-1} H^{\mathrm{H}} y.
\]

`MODE=0` 把 \(\sigma^2\) 置 0，就是 ZF。`MODE=1` 是 MMSE。硬判决只看实部和虚部的符号，0 判成比特 0，和 `qpsk_slicer` 一样。

数值边界要单独看，不能当成 `pinv`：

- H 的实部、虚部用有符号 16 位，1024 表示 1，绝对值大约别超过 8192，否则 Gram 的平方和会溢出。
- 2×2 的行列式用有符号 64 位。`det ≤ 1024` 当成奇异，输出 0。默认这组 H 的行列式大约 \(10^{12}\)，离这条线很远。
- 除法是组合逻辑的 64 位有符号除法，商朝 0 截断，再饱和到 16 位。没有流水，也没有时序报告。
- 接收噪声的标准差是刀2 QPSK 刻度 \(\sigma\) 再乘 1024，好跟 H 的「1.0」对齐。MMSE 加到 Gram 对角上的却是 \(2\sigma^2\)，不是 y 上那个已经乘过 1024 的样本方差。差一个 \(1024^2\)。若把 y 的真实方差代进去，行列式装不进 64 位，商会被截成 0。所以这条 MMSE 比「真实 y 方差」更靠近 ZF，低信噪比下仍然通常少几个差错，但不是 Julia 浮点 MMSE 的逐位实现。

不做的部分：`ml_detect` 的穷搜、`sic_detect`、Alamouti、`zf_precoder` / `svd_precoder`、注水，以及 `kronecker_channel`。相关信道不在这条 FPGA 路线里。

正交矩阵 \(\begin{bmatrix}724&724\\724&-724\end{bmatrix}\) 上，无噪声 ZF 在测试台里回到精确的 \(\pm 724\)。默认耦合矩阵上，无噪声 ZF 也回到 \(\pm 724\)；MMSE 在 \(\sigma^2=524288\)（0 dB 那一档）把幅度收到 409 和 \(-416\)，符号还在。

Icarus Verilog 12.0、256 个向量、1024 比特、当前种子的一次输出：

| \(E_b/N_0\) | ZF 差错 | ZF BER | MMSE 差错 | MMSE BER |
| --- | ---: | ---: | ---: | ---: |
| 无噪声 | 0 | 0 | 0 | 0 |
| 0 dB | 84 | \(8.20 \times 10^{-2}\) | 79 | \(7.71 \times 10^{-2}\) |
| 2 dB | 39 | \(3.81 \times 10^{-2}\) | 37 | \(3.61 \times 10^{-2}\) |
| 4 dB | 12 | \(1.17 \times 10^{-2}\) | 8 | \(7.81 \times 10^{-3}\) |
| 6 dB | 1 | \(9.77 \times 10^{-4}\) | 1 | \(9.77 \times 10^{-4}\) |
| 8 dB | 1 | \(9.77 \times 10^{-4}\) | 1 | \(9.77 \times 10^{-4}\) |

## 刀6：三档 MCS 和硬 Chase

这一刀把前面的映射、噪声、卷积和硬 Viterbi 收成一条可仿真的闭环。它不是 NR 调度器，没有 HARQ 进程号，没有增量冗余，也没有 BICM。`link_mcs_harq` 按可综合 Verilog 写成一个状态机，把已有模块串起来；`crc16_ccitt` 是独立的小移位寄存器。本仓库没有为它们跑过 Vivado，也没有接到 `ax7020_bpsk_top` 的按键或 LED。板上顶层仍然只跑刀1 的 BPSK。

每个传输块 16 个净荷比特，来自和刀1 相同的数据 LFSR（块与块之间不复位，一次 `start` 里连续往下走）。CRC-16/XMODEM 附加 16 位，高位先送。三档：

| 档 | 内容 | 一次发送的符号 | 净荷 / 符号 |
| --- | --- | ---: | ---: |
| 0 | 未编码 BPSK | 32 | 16/32 = 0.5 |
| 1 | (7, 5) 母码 + 两个 0 尾比特，再 BPSK | 68 | 16/68 |
| 2 | 未编码 QPSK | 16 | 16/16 = 1 |

频谱效率只数净荷。CRC 和卷积尾比特是开销。这不是 `teaching_mcs_table()` 的五行（净荷 48 或 41，含 QPSK 卷积、阵列 LDPC、16-QAM 和 64-QAM），也不是净荷 54 的删余表。

`snr_index` 在这一刀里是 \(E_s/N_0\)，不是刀3 的信息比特 \(E_b/N_0\)。两路 `awgn_clt` 都把 `bits_per_symbol` 接成 1，用刀1 的 k=1 刻度：BPSK 幅度 1024，QPSK 幅度 724，符号能量都按 \(1024^2\) 对齐，\(\sigma^2\) 的目标和 `awgn(..., energy_per_bit=1)` 一样。`snr_index ≥ 6` 仍是无噪声。BPSK 只用 I 路噪声；Q 路发生器在 BPSK 的发送拍上也会往前走，序列是确定的，只是那一拍的样本没被加到符号上。

选档：`mcs_force` 为 0、1、2 时强制该档；为 3 时看 `snr_index`。0 和 1（0 dB、2 dB）选档 0，2 和 3（4 dB、6 dB）选档 1，其余（含无噪声）选档 2。这些边界是这三档自己的，**不是** `teaching_amc_thresholds()` 的 \([-\infty, 2, 3, 7, 11]\) dB，也没有外环偏移。

重传的是同一份符号。硬 Chase 对每个信道比特做有符号累加：硬判决 0 加 +1，硬判决 1 加 −1。合并后再判决，和为负才判 1，和为 0 判 0。卷积档把合并后的 68 个硬比特送进硬 Viterbi，再对 32 个译码比特做 CRC。不是 LLR 相加。CRC 通过就停，并给 `n_good` 加 16；否则若已发送次数还小于 `max_tx` 就再发一次，累加器不清零。`max_tx` 含第一次发送，0 被收成 1，最大 7（累加器是 4 位有符号）。失败的块丢掉，不进入 `n_good`。

计数：

- `n_ack`、`n_fail`：CRC 通过和最终失败的块数。BLER 用 `n_fail / NBLOCKS`。
- `n_uses`：发出的信道符号数，重传也算。BPSK 一个符号一个比特，QPSK 一个符号两比特。
- `n_good`：通过 CRC 的净荷比特。吞吐用 `n_good / n_uses`。

默认 `NBLOCKS = 8`。Icarus Verilog 12.0、当前种子的一次输出：

| 条件 | ACK | 失败 | 信道使用 | 净荷比特 | BLER | 吞吐 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 无噪声，档 0，`max_tx=1` | 8 | 0 | 256 | 128 | 0 | 0.5 |
| 无噪声，档 1，`max_tx=1` | 8 | 0 | 544 | 128 | 0 | 16/68 |
| 无噪声，档 2，`max_tx=1` | 8 | 0 | 128 | 128 | 0 | 1 |
| 0 dB，档 0，`max_tx=1` | 1 | 7 | 256 | 16 | 7/8 | 16/256 |
| 0 dB，档 0，`max_tx=3` | 2 | 6 | 768 | 32 | 6/8 | 32/768 |
| 0 dB，档 1，`max_tx=1` | 6 | 2 | 544 | 96 | 2/8 | 96/544 |
| 0 dB，档 1，`max_tx=3` | 8 | 0 | 952 | 128 | 0 | 128/952 |

8 个块太短，不能拿来和 Julia 的 80 块吞吐扫频比，也不能说已经贴近联合界。0 dB 的未编码 BPSK 上，三次硬合并只多救回 1 块，吞吐还因为重传下降。卷积档在同一种子上从 2 块失败变成 0，信道使用从 544 增到 952。这只说明这条定点闭环里合并和早停在起作用。

## 和 v1.0 的对照

下表是这条 FPGA 路线相对 Julia **v1.0.0** 公开面的收口说明。左列是包里已经有的能力，右列是本目录做到哪。对齐一律是概念和定点骨架，不是数值契约。

| Julia v1.0 | 本目录 | 结论 |
| --- | --- | --- |
| BPSK 硬映射 / 硬判决，实 AWGN，BER | 刀1 | 已覆盖。整数 ±1024，噪声是均匀和。 |
| Gray QPSK，方形 Gray 16-QAM / 64-QAM 硬判决 | 刀2 | 已覆盖。不是 `PSK(4)`，门限平局规则可能和 `findmin` 不同。 |
| 终止 (7, 5) 母码，硬 Viterbi | 刀3 | 已覆盖。无噪声硬比特可以对上固定向量。 |
| CRC-16/XMODEM，`crc_append` / `crc_check` | 刀6 | 已覆盖。`"123456789"` 余数 `0x31C3`。 |
| 离散 MCS、按 \(E_s/N_0\) 选档、Chase、BLER 和吞吐 | 刀6 | 教学骨架。三档、净荷 16、硬 ±1 合并、本地门限。不是五档表，不是软 LLR Chase，不是外环。 |
| CP-OFDM：IFFT/FFT、循环前缀 | 刀4 | 骨架。N=16、CP=4、全数据子载波。没有同步，没有信道估计。 |
| 平坦 MIMO 的 ZF / MMSE，理想 CSI | 刀5 | 已覆盖 2×2 固定实 H。不是 i.i.d. 瑞利，MMSE 的噪声方差刻度和浮点不完全同一套。 |
| Fibonacci LFSR 扰码的移位规则 | 数据 `lfsr_fib` | 15 级多项式和全 1 种子的比特序列可以对上。噪声 LFSR 只借用 31 级多项式。 |

刻意不做，也不在后续刀次里补：

| 留下不做的 | 原因 |
| --- | --- |
| 软判决全链：逐比特 LLR、软 Viterbi、log-BCJR、SOVA、`demodulate_soft` | 刀3 和刀6 都停在硬度量。删余要往丢掉的母比特上写 LLR 0，所以一并不做。 |
| 删余 2/3、3/4，以及 `teaching_punctured_mcs_table` | 没有删除符号，也没有那张净荷 54 的表。 |
| BICM 和 BICM-ID | 没有比特交织，没有迭代解调。 |
| LDPC（Gallager、Fan 阵列、NR 基图） | 五档表里的阵列码行不在 FPGA 上。 |
| Turbo | 短帧 RSC 不在这条路线。 |
| Kronecker、Loyka 指数相关、Alamouti、ML、SIC、SVD、注水 | 刀5 只做固定 H 上的 ZF/MMSE。 |
| TR 38.901 的 CDL、几何信道；`TDLChannel` 的多径 | 噪声只有 AWGN。OFDM 前缀没有用来吸收时延。 |
| EXIT 图 | 没有软信息的互信息估计。 |
| 路线 C 的多用户 xApp、NR 调度器、资源栅格、O-RAN / eCPRI | 单链路教学计数器，没有用户间调度。 |
| 增量冗余 | 重传的是同一份符号。 |
| 外环 AMC | 门限不随 ACK/NACK 移动。 |
| OFDM 的 STO/CFO、LS/MMSE 信道估计 | 刀4 假定对齐，H=1。 |
| 线性均衡、LMS、DFE；Gardner、Costas、Schmidl–Cox | 没有定时或载波恢复。 |
| BCH、Hamming、Reed–Solomon | 帧校验只用 CRC-16。 |
| Rice、Jakes、Nakagami-m、MRC、脉冲成形、眼图 | 不在 PL 教学目录里。 |

**完成声明。** FPGA 教学复刻到 v1 主链路，到刀6 为止，在 `fpga/ax7020` 收口。已经落地的是上表「已覆盖」和「教学骨架」那几行，并且十三个 Icarus 测试台可以跑完。边界也写在同一张表里：软判决全链、LDPC、Turbo、Kronecker、TR 38.901、EXIT、路线 C 不在这条 FPGA 路线中。没有比特流，没有板上实测，计数也不是和 Julia 浮点的逐块对照。

## 上板

**没有实测。** 下面是把现有 RTL 放进 Vivado 的步骤，不是已经生成并下载过的工程。不要把比特流、IP 缓存或工程二进制提交进仓；`fpga/ax7020/.gitignore` 排除了常见的 `.bit`、`.xpr`、`.cache`、`.runs`、`.gen` 等。

工具链：Xilinx Vivado，**推荐 2020.x 或 2022.x** 这种还能覆盖 Zynq-7000 的常见版本。不锁小版本。器件库里有 XC7Z020 的更新版本也可以。仿真不依赖 Vivado。

器件按芯片丝印选。AX7020 常见是 **XC7Z020-2CLG400I**（工程里往往写成 `xc7z020clg400-2`）。以板上的丝印为准。

1. 新建 RTL 工程，不要建 Block Design。本顶层不用 PS、AXI、DDR 或 UART。
2. 加入 `rtl/` 下全部 `.v`，以及 `conv75_trellis.vh`（和 `.v` 放在同一目录，或把该目录加进 include 路径）。不要加入 `tb/`。现成顶层只用 BPSK。
3. 顶层设为 `ax7020_bpsk_top`（只演示 BPSK）。
4. 加入 `constr/ax7020_bpsk.xdc`。先对照原理图改引脚和 `IOSTANDARD`。
5. 综合、实现、生成比特流。本仓库没跑过这一步，所以没有时序报告可以引用。50 MHz 直接用 PL 脚上的时钟，没有 Clocking Wizard。
6. 用板卡自己的 JTAG 流程下载。PL 时钟是 **U18 的 50 MHz**，不是 PS 那颗 33.333 MHz（E7）。若下载后心跳灯不闪，先查时钟脚和 BANK 电压。

按键和灯（低电平有效，仍以原理图为准）：

| 信号 | 手册常见引脚 | 作用 |
| --- | --- | --- |
| `clk` | U18 | PL 50 MHz |
| `rst_n` | N15（KEY1） | 复位 |
| `key_n` | N16（KEY2） | 在 0/2/4/6/8/10 dB 之间循环，LED2..LED4 显示索引的低三位 |
| `led[0]` | M14（LED1） | 心跳，只说明时钟在走 |
| `led[3:1]` | J16、K16、M15 | `snr_index` |

板上的帧长默认是 \(2^{23}\) 个符号（50 MHz 下大约 0.17 s 一帧），然后自动重开。差错计数在模块内部的 `n_err`、`n_bit`，**没有接到灯或串口**。要看数，自己加 ILA 或 UART；本仓库不附带 ILA IP。灯不能代替 BER 曲线。趋势以仿真打印为准。

## 许可

与仓库根目录相同，MIT。Copyright (c) 2026 刘恒宇 (longee142857)。
