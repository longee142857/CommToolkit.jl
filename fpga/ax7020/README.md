# AX7020 教学链路

黑金 ALINX **AX7020**（Zynq-7020，器件 **XC7Z020**，常见封装 CLG400）上的 **PL 教学**示例。和 Julia 包 `CommToolkit` 里的 `Modulation` / `Channel` / `Metrics` / `Scrambling` **概念对齐**。定点噪声是 12 个离散均匀数的和，不是 `randn`。**不要把这里的误码率当成和 Julia 浮点仿真 bit-exact，也不要当成已经和闭式公式做过区间验收。** 本目录没有附带定点对照的 Julia 脚本；没做那一步，就不声称数值一致。

| 刀 | 内容 | 状态 |
| --- | --- | --- |
| 1 | BPSK + 实 AWGN + BER | 已在 main |
| 2 | QPSK，以及方形 Gray 16-QAM、64-QAM 的硬映射 / 硬判决 | 当前目录 |
| 以后 | OFDM、MIMO、信道码、软判决 | 不做 |

刀1 的链路是：PRBS 比特 → BPSK 映射 → 实 AWGN → 硬判决 → 差错计数。刀2 在旁边加复星座，噪声源和差错计数的做法不变。

## 边界

- 教学级定点，不是量产调制解调器。没有 AGC、成形滤波、定时或载波恢复。
- 不是 NR，不是 LTE，不是 O-RAN。没有资源栅格、前传、eCPRI 或 xApp。
- 和 Julia `CommToolkit` 是概念对照，不是数值契约。数据 LFSR 的移位规则可以和 `lfsr_sequence` 对上（见下表）；噪声和 BER 不对上。
- 刀1 是 **BPSK + 实 AWGN + BER**。刀2 加上 **QPSK** 和方形 Gray **16-QAM / 64-QAM** 的硬判决。仍然不做 OFDM、MIMO、信道编码、交织、均衡、同步或软信息。
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
  tb/tb_bpsk_unit.v         刀1：映射、判决、LFSR、噪声刻度
  tb/tb_bpsk_awgn_ber.v     刀1：BER 对 Eb/N0
  tb/tb_ax7020_top.v        顶层复位和按键的短冒烟
  tb/tb_qpsk_qam_unit.v     刀2：星座坐标、往返、噪声刻度
  tb/tb_qpsk_qam_ber.v      刀2：BER 对 Eb/N0
  constr/ax7020_bpsk.xdc    占位引脚，需按板本改。刀2 没有新管脚
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

脚本依次跑五个测试台（刀1 三个，刀2 两个）。每个测试台打印一行 `PASS`。趋势测试台会打印无噪声帧和 0、2、4、6、8、10 dB 的差错计数。刀1 的 `NBITS = 16384`。刀2 的 `NBITS = 16380`，这样 2、4、6 比特每符号都能整除。种子写在 RTL 里。同一版本的 Icarus 应得到同一张表。

通过条件：

- 无噪声（`snr_index = 7`）差错为 0，且数据 LFSR 前 32 位与 Fibonacci 规则一致。
- 0 dB 到 8 dB 差错数严格下降。
- 10 dB 的差错数不超过 8 dB。
- 刀1：0 dB 的差错比例落在 5% 到 12%（教材 \(Q\) 约 7.9%）。4 dB 落在 0.4% 到 4%（教材约 1.25%）。
- 刀2：QPSK 的 0 dB 落在 4% 到 15%，4 dB 落在 0.4% 到 4%。16-QAM 的 0 dB 落在 6% 到 25%。64-QAM 的 0 dB 落在 8% 到 35%。三种都要求 0 dB 到 8 dB 严格下降。

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

板上的 `ax7020_bpsk_top` 仍然只跑 BPSK。刀2 的模块可以综合，但没有接到按键或 LED，也没有新的 `.xdc`。要上板得自己改顶层。本仓库没有为刀2 跑过 Vivado，也没有比特流。

## 上板

**没有实测。** 下面是把现有 RTL 放进 Vivado 的步骤，不是已经生成并下载过的工程。不要把比特流、IP 缓存或工程二进制提交进仓；`fpga/ax7020/.gitignore` 排除了常见的 `.bit`、`.xpr`、`.cache`、`.runs`、`.gen` 等。

工具链：Xilinx Vivado，**推荐 2020.x 或 2022.x** 这种还能覆盖 Zynq-7000 的常见版本。不锁小版本。器件库里有 XC7Z020 的更新版本也可以。仿真不依赖 Vivado。

器件按芯片丝印选。AX7020 常见是 **XC7Z020-2CLG400I**（工程里往往写成 `xc7z020clg400-2`）。以板上的丝印为准。

1. 新建 RTL 工程，不要建 Block Design。本顶层不用 PS、AXI、DDR 或 UART。
2. 加入 `rtl/` 下全部 `.v`。不要加入 `tb/`。刀2 的映射模块可以加进去，但现成顶层不用它们。
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
