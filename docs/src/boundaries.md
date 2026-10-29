# 不是什么

本包把教材基带算法和闭式或参考实现对上。下面是类别上的边界。参数级的限制在 [已知限制](limitations.md)。验证到哪一步在 [验证矩阵](validation.md)。

1.0.0 不扩大这些类别。1.1.0 收进了短 Reed–Solomon：`rs_code` 的窄义本原码，``\mathrm{GF}(2^m)``、``m \le 4``，长度 3、7、15 及其缩短，硬判决有界距离译码，不是 CCSDS 长码，也不是软代数译码。1.1.1 把 EXIT 图收进范围：直方图与蒙特卡洛软信息互信息测量、写明的 ``J(\sigma)`` 高斯对称近似，以及在现有解映射器和译码器上测的传递曲线。这是教学 EXIT，不是 NR LDPC 的 EXIT 生产工具，也没有密度进化。1.1.2 把教学增量冗余收进链路范围：`simulate_link(...; harq=:ir)` 在删余卷积行上按互补图案重传、按母码位软合并，默认 `harq=:chase` 不变；它不是 NR 速率匹配，也不是冗余版本 RV0–RV3。1.1.3 把教学簇延迟线收进范围：TR 38.901 CDL-A/B/C 表的归一化延迟与功率两列、式 (7.7-1) 的延迟缩放、样本量化与合并、每抽头一个瑞利增益、可选逐抽头 Clarke 多普勒。有意留在范围之外的仍是完整几何 TR 38.901 CDL（离开角 / 到达角、簇角扩展、XPR、每簇 20 射线、天线几何、空间一致性、LOS 的 CDL-D/E、MIMO 全矩阵）、路线 C 的多用户 xApp，以及 NR LDPC。教学 MCS / HARQ 上的删余和 BICM 路径，以及平坦 Kronecker MIMO，留在范围内。

## 不是 3GPP NR / LTE 物理层

OFDM 是教科书循环前缀模型。没有资源栅格，没有 DMRS 或 PTRS，没有 PSS / SSS，没有空口时隙。导频是固定线性同余序列，不是 Gold 码，也不是 Zadoff–Chu。

LDPC 是短的 Gallager 正则构造和 Fan 阵列码，不是 NR BG1/BG2。Turbo 是记忆 2 的递归系统卷积 (1, 5/7)，信息长度不超过 1024，不是 3GPP 的 (13, 15) 或 (15, 13)，也不是 LTE 的记忆 4 分量码。卷积删余是 ``K = 3`` 母码上的两张周期图案，不是 NR 速率匹配。BCH 只有二进制 (15, 7)。Reed–Solomon 只有 `rs_code` 的短码：``\mathrm{GF}(2^m)``、``m \le 4``，长度 3、7、15 及其缩短，Berlekamp–Massey 硬判决。不是 CCSDS 的 (255, 223) 或 (255, 239)，不是软代数译码，也不是 NR 信道码。

BICM 用现有星座和现有译码器做教学迭代。没有 NR 速率匹配、层映射、扰码或 NR LDPC 基图；EXIT 图不在迭代环内，是 `EXIT` 模块在环外的互信息测量。`simulate_link(...; receiver=:bicm)` 可以在卷积 MCS 行上调用它。LDPC 和 Turbo 的 BICM 不进 MCS 表。

扰码是加性 LFSR，不是 NR 扰码。CRC 是 CRC-16/XMODEM 那条多项式，不是 NR 的 CRC24。

## 不是 O-RAN

没有前传映射，没有 eCPRI，没有开放前传的功能切分，也没有路线 C 的多用户 xApp。`OFDMParams` 的文档写明：子载波图不是 O-RAN 前传映射。

## 不是完整几何 38.901 CDL

1.1.3 收进的 `CDLProfile` / `cdl_profile` / `CDLTaps` / `cdl_taps` / `apply_cdl` / `cdl_channel` 是教学簇延迟线：只有 CDL-A/B/C 表的归一化延迟与功率两列、式 (7.7-1) 的延迟缩放、按样本量化并合并、每抽头一个瑞利增益（复用 `jakes_fading` 的 Clarke 谱），外加教学便利的 `n_taps` 抽稀。它不是完整 3GPP 信道实现，也不声称 NR 合规：没有离开角或到达角，没有簇角扩展与 XPR，没有每簇 20 条射线的系数生成（§7.5 第 10–11 步），没有天线几何与方向图，没有空间一致性，没有 LOS 的 CDL-D/E，也没有 MIMO 空间相关全矩阵。逐抽头多普勒是 Clarke 谱，不是 38.901 的射线求和多普勒。

`TDLChannel` 是调用者给定的静态抽头：时延用样本数，增益用 dB，相位来自固定线性同余。抽头功率不归一。文档字符串写明它不是 3GPP TR 38.901 的 TDL-A/B/C/D/E，也不是衰落过程。静态 CDL drop 可以包成 `TDLChannel` 接进 OFDM 链，但 `TDLChannel` 本身仍然不是 CDL：没有归一化延迟剖面，没有瑞利簇增益，也没有多普勒。

`jakes_fading` 是另一条平坦过程。谱是 Clarke 的 U 形，实现是 DFT 格子上对 ``S(f)`` 积分后再做 `ifft`，不是 Jakes 1974 的正弦和，也不是频率选择性抽头，没有接到 `TDLChannel`。

单载波 `apply_fir` 的抽头由调用者给定，不随时间变化，也不是 38.901 的抽头剖面。`mimo_channel` 仍是 i.i.d. 瑞利矩阵。`kronecker_channel` 是平坦 Kronecker 空间相关：相关矩阵由调用者给定，或由 Loyka 指数模型 `exponential_correlation` 生成。它不是 Weichselberger 的两端联合耦合，也不是 TR 38.901 的 CDL 或 UMi。没有离开角或到达角，没有簇，也没有天线方向图。

## 不是生产级同步或调度器

`Sync` 里的 Gardner、Costas 和 Schmidl–Cox 是给定捕获范围的教学环。没有锁定检测，没有 AGC，没有 PSS。Gardner 要求偶数每符号样本数，不估计整数个符号的延迟。Costas 只覆盖 `BPSK()` 和 `QPSK()`。Schmidl–Cox 必须用它自己的重复半符号前导，`ofdm_tx` 不会插入。

`ofdm_sync` 仍是 van de Beek 的循环前缀相关。捕获范围是 ``(-1/2, 1/2)`` 个子载波间隔，不估计整数倍频偏，也不按定时估计平移缓存。Schmidl–Cox 不修改这三个函数。

`Link` 是单链路 AWGN 上的教学 MCS、CRC-16 和有限次 HARQ。原来的表仍是五档。另一张表是净荷 54 比特的删余卷积码，名义码率 1/2、2/3、3/4，默认重传重复同一删余字（Chase），`harq=:ir` 时偶数次重传改发补图案（教学增量冗余，不是 NR 速率匹配，也不是 RV0–RV3）。没有 HARQ 进程号，没有时序，没有多用户，没有调度器，没有功率控制。门限不是 NR CQI。SOVA 和 Turbo 没有接到这两张表上。`receiver=:bicm` 只接受卷积行，默认接收机仍是软输入 Viterbi。
