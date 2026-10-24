# 不是什么

本包把教材基带算法和闭式或参考实现对上。下面是类别上的边界。参数级的限制在 [已知限制](limitations.md)。验证到哪一步在 [验证矩阵](validation.md)。

## 不是 3GPP NR / LTE 物理层

OFDM 是教科书循环前缀模型。没有资源栅格，没有 DMRS 或 PTRS，没有 PSS / SSS，没有空口时隙。导频是固定线性同余序列，不是 Gold 码，也不是 Zadoff–Chu。

LDPC 是短的 Gallager 正则构造和 Fan 阵列码，不是 NR BG1/BG2。Turbo 是记忆 2 的递归系统卷积 (1, 5/7)，信息长度不超过 1024，不是 3GPP 的 (13, 15) 或 (15, 13)，也不是 LTE 的记忆 4 分量码。卷积删余是 ``K = 3`` 母码上的两张周期图案，不是 NR 速率匹配。BCH 只有二进制 (15, 7)，不是 Reed–Solomon。

BICM 用现有星座和现有译码器做教学迭代。没有 NR 速率匹配、层映射、扰码或 NR LDPC 基图，也没有 EXIT 图。`simulate_link(...; receiver=:bicm)` 可以在卷积 MCS 行上调用它。LDPC 和 Turbo 的 BICM 不进 MCS 表。

扰码是加性 LFSR，不是 NR 扰码。CRC 是 CRC-16/XMODEM 那条多项式，不是 NR 的 CRC24。

## 不是 O-RAN

没有前传映射，没有 eCPRI，没有开放前传的功能切分。`OFDMParams` 的文档写明：子载波图不是 O-RAN 前传映射。

## 不是几何 38.901 CDL

`TDLChannel` 是调用者给定的静态抽头：时延用样本数，增益用 dB，相位来自固定线性同余。抽头功率不归一。文档字符串写明它不是 3GPP TR 38.901 的 TDL-A/B/C/D/E，也不是衰落过程。它同样不是 TR 38.901 的簇延迟线（CDL）：没有离开角或到达角，没有天线几何，没有空间簇。

`jakes_fading` 是另一条平坦过程。谱是 Clarke 的 U 形，实现是 DFT 格子上对 ``S(f)`` 积分后再做 `ifft`，不是 Jakes 1974 的正弦和，也不是频率选择性抽头，没有接到 `TDLChannel`。

单载波 `apply_fir` 的抽头由调用者给定，不随时间变化，也不是 38.901 的抽头剖面。`mimo_channel` 仍是 i.i.d. 瑞利矩阵。`kronecker_channel` 是平坦 Kronecker 空间相关：相关矩阵由调用者给定，或由 Loyka 指数模型 `exponential_correlation` 生成。它不是 Weichselberger 的两端联合耦合，也不是 TR 38.901 的 CDL 或 UMi。没有离开角或到达角，没有簇，也没有天线方向图。

## 不是生产级同步或调度器

`Sync` 里的 Gardner、Costas 和 Schmidl–Cox 是给定捕获范围的教学环。没有锁定检测，没有 AGC，没有 PSS。Gardner 要求偶数每符号样本数，不估计整数个符号的延迟。Costas 只覆盖 `BPSK()` 和 `QPSK()`。Schmidl–Cox 必须用它自己的重复半符号前导，`ofdm_tx` 不会插入。

`ofdm_sync` 仍是 van de Beek 的循环前缀相关。捕获范围是 ``(-1/2, 1/2)`` 个子载波间隔，不估计整数倍频偏，也不按定时估计平移缓存。Schmidl–Cox 不修改这三个函数。

`Link` 是单链路 AWGN 上的教学 MCS、CRC-16 和有限次 Chase 合并。原来的表仍是五档。另一张表是净荷 54 比特的删余卷积码，名义码率 1/2、2/3、3/4，重传重复同一删余字。没有 HARQ 进程号，没有时序，没有多用户，没有调度器，没有功率控制，没有增量冗余。门限不是 NR CQI。SOVA 和 Turbo 没有接到这两张表上。`receiver=:bicm` 只接受卷积行，默认接收机仍是软输入 Viterbi。
