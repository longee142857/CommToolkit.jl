# 变更记录

## 0.9.2

把卷积删余和一条 BICM 接收机接到已有的 MCS / Chase 闭环上。不是 NR CQI，不是 NR 速率匹配，不是增量冗余，也没有 EXIT 图。`teaching_mcs_table` 的五档、门限和默认 `receiver=:viterbi` 都不变。已有函数的位置参数和返回值不变。这一版坐在 0.9.1 的公开接口契约上：稳定清单增加两个导出名字，不改名、不删除，也不改已有位置参数和返回值。它把现成的删余和 BICM 接进链路，不新开调制或信道族，所以不升 0.10，也不升 1.0，记为 0.9.2。

- 删余 MCS：`mcs_convolutional(...; rate=1//2)` 也可取 `2//3`、`3//4`，图案就是 `puncture_mask`。母码长度必须是图案周期的整数倍，尾比特算在里面。`n_coded` 是删余后真正调制的比特数，`n_mother` 是删余前的长度。频谱效率仍是净荷比特 / 发送符号。LDPC 行没有删余码率。
- 删余表：`teaching_punctured_mcs_table`。净荷 54 比特，母码 144 比特。五档按频谱效率递增：BPSK 1/2（``54/144``）、QPSK 1/2（``54/72``）、QPSK 2/3（``54/54``）、QPSK 3/4（``54/48``）、16-QAM 1/2（``54/36``）。门限 `teaching_punctured_amc_thresholds` 是 ``-\infty, 2, 4, 6, 8`` dB，来自 Julia 1.10.12、种子 `0x0c0de1162`、每点 80 块、无 HARQ 的固定 MCS 吞吐包络。5 dB 上码率 2/3 和 3/4 差大约 0.01 bit/use；同一条种子在 Julia 1.13.2 上这一档会换边。测试用 0 dB、3 dB 和 10 dB，不拿交点做判据。这组门限不能套到原来的五档表上。
- Chase：重传的是同一个删余字。Viterbi 路径把发送比特的 LLR 相加，删掉的母码位保持 LLR 0，再进软 Viterbi。不是把删掉的比特在下一次重传里补上。
- BICM 路径：`simulate_link(...; receiver=:bicm)`。只接受卷积行。顺序是编码、删余、`bicm_n_rows` 块交织、调制。译码是 `bicm_decode` 的 `:bcjr`。多次传输把同一符号平均，噪声方差除以已合并的次数，这是重复 AWGN 符号的充分统计，不是比特度量 Chase。`bicm_iterations` 默认 1（没有外信息反馈），不是 `bicm_decode` 的默认 4。`bicm_labels=nothing` 是 Gray；8-PSK 集合分割要调用者传入 `bicm_psk8_id_labels`。``E_s/N_0 = +\infty`` 时 BICM 用 ``\sigma^2 = 10^{-6}``，因为解调器不接受 0。Viterbi 接收机不能顺便带上这些关键字。
- `bicm_map` / `bicm_decode` 增加关键字 `rate`。默认 `1//2` 不删余，原有调用的码长不变。`2//3` 和 `3//4` 只给卷积的 `:conv` / `:bcjr` / `:sova`。外信息先删余再交织，送回的是实际调制过的比特。
- 对照：`test/link.jl`。3 dB 自适应（QPSK 1/2）相对固定 BPSK 和固定 16-QAM，10 dB 自适应（16-QAM）相对固定 BPSK，吞吐用 Clopper–Pearson（``\alpha = 10^{-4}``，80 块）分开。2 dB 的 QPSK 3/4，3 次 Chase 相对单次的 CRC 失败率区间分开，吞吐更高。8-PSK 集合分割、12 行交织、max-log、``E_s/N_0 = 6\,\mathrm{dB}``、80 块：链路内 4 次迭代相对 1 次，净荷比特错误率和 CRC 失败率的区间都分开。随机数是全局 `Random.seed!`，与原来的链路测试相同。不宣称逼近容量，不宣称 NR。打印表在 `demos/link_puncture_bicm.jl`。
- 契约：`docs/src/api.md` 增加 `teaching_punctured_mcs_table` 和 `teaching_punctured_amc_thresholds`。边界、已知限制和验证矩阵改成与这条链路一致。SOVA 和 Turbo 仍然不进 MCS 表。

## 0.9.1

工程冻结这一刀：公开接口当成契约写下来，验证对照收成一张表，README 改成导航。不是新算法，不改已有函数的位置参数和返回值，也不升到 1.0.0。0.9.0 已经把 BICM-ID 收进 0.9，这一版按补丁记为 0.9.1，而不是 0.10.0。

- 稳定面：`docs/src/api.md` 列出 `using CommToolkit` 带入的全部名字，按子模块分组。`test/contract.jl` 要求这份清单等于 `names(CommToolkit)`，并且各子模块 `export` 的并集只多一个名字。核对结果：没有漏导出，也没有改名。`CommToolkit.modulation` 仍不导出，和 0.9.0 的 README 相同。
- 验证：`docs/src/validation.md`。每一行是模块、量、闭式或参照、测试文件、演示。公式只引用 `test/reference.jl` 和测试文件头里已经写明的文献，不新造表达式。`test/` 下每个 `*.jl` 都在这一页出现。
- 边界：`docs/src/boundaries.md`。不是 3GPP NR / LTE 物理层，不是 O-RAN，不是 TR 38.901 的几何 CDL，不是生产级同步或调度器。`TDLChannel` 的文档字符串补上同一句：静态抽头，不是 38.901 的 TDL-A/B/C/D/E，也不是 CDL。
- 导航：README 保留安装、四条曲线的图和源码布局。长公式改到验证页，参数级限制改到 `docs/src/limitations.md`，正文与 0.9.0 的「已知限制」相同。
- 文档站：`docs/make.jl` 把边界、公开接口、验证、已知限制放进 Documenter 侧栏。API 参考仍由 docstring 生成。

## 0.9.0

比特交织编码调制的迭代解调 / 译码（BICM-ID）。不是 5G NR BICM，没有 EXIT 图，没有接到 MCS / HARQ，也没有动 Jakes、同步或路线 C。P2 的这一块到此收口，所以记为 0.9.0，而不是 0.8.6，也不升到 1.0。已有函数的位置参数和返回值不变。

- 映射：`bicm_map`。`:conv` 走 `conv_encode`，`:ldpc` 走 `ldpc_encode`，`:turbo` 走 `turbo_encode`。可选 `block_interleave`（`n_rows = 0` 表示不交织）。默认标签是现有 Gray 星座的下标，MSB 在前。`bicm_psk8_id_labels` 把 8-PSK 在角度 ``2\pi m/8`` 上标成 ``(0,4,2,6,1,5,3,7)``：最高位是 ``m \bmod 2``，每个子集是最小距离 ``\sqrt{2}`` 的 QPSK。这是 Ungerboeck 的集合分割，不是 Gray，也不是 3GPP 标签。
- 解映射：`bicm_demap`。Gray、零先验、`gains === nothing` 时直接返回 `demodulate_soft`。有先验时，比特 ``i`` 的外信息不再加进自己的度量，精确边际或 max-log 用的仍是 ``|y-hs|^2/N_0``，``N_0 = 2\sigma^2``。`gains` 是逐符号理想 CSI。不做 DBPSK / DQPSK。BPSK 和 Gray QPSK 的外信息与先验无关。
- 迭代：`bicm_decode`，1 到 32 次，默认 4。`:bcjr` 在终止的 (7, 5) 上做编码比特 log-MAP，分支度量和尾比特边界与 `bcjr_decode` 相同；信息比特 LLR 在一次迭代、零信息先验时与 `bcjr_decode` 一致。编码比特外信息截到 `llr_clip`（默认 50）后再交织，作为下一次解调的先验。`:ldpc` 调用 `ldpc_decode_details`，`early_stop=false`，外信息是后验减去解调器送来的信道 LLR。`:sova` 和 `:turbo` 只能 `iterations = 1`：公开接口没有编码比特外信息。
- 对照：`test/bicm.jl`。AWGN、``E_b/N_0 = 6\,\mathrm{dB}``、80 帧、每帧 400 个信息比特、12×67 块交织、max-log、集合分割 8-PSK。4 次迭代相对 1 次，信息比特 BER 和 BLER 的 Clopper–Pearson 区间（``\alpha = 10^{-4}``）互不重叠。同一条种子下，理想 CSI 的逐符号瑞利在 8 dB、40 帧，也分开。随机数是 xorshift，不是 Julia 全局 `randn`。不宣称逼近容量。Gray 标签不拿来做这条增益曲线。

## 0.8.5

平坦 Clarke/Jakes 时间相关瑞利，以及逐符号独立的 Nakagami-m。不是 Jakes 1974 的正弦和，不是频率选择性抽头，没有接到 OFDM、MIMO 或 MCS / HARQ 上。

- Jakes：`jakes_fading(n, fd, Ts)`。`fd` 是最大多普勒（Hz），`Ts` 是采样间隔（秒），`fd * Ts` 必须小于 1/2。谱是 Clarke 的 U 形 ``S(f) = 1/(\pi f_d \sqrt{1-(f/f_d)^2})``，``|f| < f_d``。实现是 `n` 点 DFT：每个格子的功率等于 ``S`` 在该格子上的积分（原函数 ``\arcsin(f/f_d)/\pi``），非零格子上放独立 ``\mathcal{CN}(0, P_k)``，时间序列是 ``n \cdot \mathrm{ifft}``。集合自相关是这些格子功率的逆 DFT，格子变细时趋于 ``J_0(2\pi f_d \tau)``。序列以 `n` 为周期。`fd = 0` 不跑 FFT，整段是同一个 ``\mathcal{CN}(0,1)``。
- 信道：`jakes_channel` 在 ``h x`` 之后按 `awgn` 的 ``E_b/N_0`` 加噪。每个样本的边际仍是 ``\mathcal{CN}(0,1)``，所以理想 CSI 相干 BPSK 的平均误码率与逐符号瑞利相同。时间相关改变的是突发，不是这个平均。
- Nakagami-m：`nakagami_envelope` 生成 ``\mathbb{E}[R^2] = 1`` 的幅度，`m ≥ 1/2`。`nakagami_fading` 是复包络。`m = 1` 与 `rayleigh_channel` 使用同一套 `randn` 顺序。`2m` 在 `Float64` 里是整数时，幅度平方是 `2m` 个标准正态平方的平均；其余 `m` 用 Marsaglia–Tsang 的 gamma 抽样，相位独立均匀。`m = Inf` 是单位幅度加均匀相位，不是 Rice `K = Inf` 的 ``h = 1``。
- 信道：`nakagami_channel`，`m` 是必填关键字。理想 CSI 相干 BPSK 对照 Gamma 矩生成函数的 Craig 积分。整数 `m` 等于 `m` 条 i.i.d. 瑞利支路的 MRC，每条支路的平均 ``E_b/N_0`` 是 ``\bar\gamma_b/m``。
- 版本：已有函数的位置参数和返回值不变。0.8.4 已经是 Gardner、Costas 与 Schmidl–Cox。这一刀是新的信道族，按 0.x 记为 0.8.5，而不是再占 0.8.4，也不升到 1.0。

## 0.8.4

单载波定时和载波的教学环，以及 OFDM 的 Schmidl–Cox 前导。不是生产级同步。没有接到 MCS / HARQ，没有改 van de Beek 的 `ofdm_sync`，也没有动 SOVA、Turbo、BCH、Jakes、NR 或 O-RAN。0.8.0 的 EVM、LFSR 扰码和 SVD 注水，以及 0.8.3 的 SOVA、短 Turbo 和 BCH(15,7) 保持原样。

- Gardner：`gardner_ted` 要求偶数 `sps ≥ 2`，和 `pulse_shape` / `downsample_symbols` 用同一个每符号样本数。符号点是 1-based 下标 `phase, phase+sps, ...`，中点在当前符号点前 `sps÷2` 个样本。误差是 ``\Re(y_{\mathrm{mid}}^*(y_{\mathrm{late}}-y_{\mathrm{early}}))``。采样点偏晚（落在自相关峰之后）时平均误差为正，所以环路用 ``\mu \leftarrow \mu - \alpha e``。`gardner_loop` 是线性内插的一阶环，没有单独的符号钟频率状态。`gardner_open_loop` 在一个符号周期里找正斜率过零。整数个符号的延迟不由 Gardner 给出。`single_carrier_timing` 把 `pulse_shape`、`matched_filter`、分数延迟和 `downsample_symbols` 串起来，并用 `length(taps)` 算出输入符号 1 落在哪一个抽样点。
- Costas：只接受 `BPSK()` 和 `QPSK()`。`QPSK()` 在 ``\pi/4``，不是轴上的 `PSK(4)`。小残余相位下鉴相器与相位同号。BPSK 的锁定周期是 ``\pi``，QPSK 是 ``\pi/2``。默认增益 ``\alpha = 0.06``、``\beta = 0.0015`` 只覆盖大约 ``10^{-2}`` 弧度/符号这一量级的频偏。没有锁定检测，没有 AGC，输入必须已经是每符号一个样点。
- Schmidl–Cox：`schmidl_preamble` 在 0-based 偶数子载波上放确定性 ``\pm 1``，两半时域相同。`schmidl_cox` 的小数 CFO 在 ``(-1, 1)`` 个子载波间隔内无模糊，对再加上 ``2z`` 的偶数倍整数频偏不敏感。`integer_cfo_estimate` 在调用者给定的 ``z`` 范围内做已知前导的相干相关，默认 ``z_{\max} = 4``（即 ``|m| \le 8`` 的偶数）。它假定平坦信道，并且 FFT 窗已经落在无 ISI 的有用符号上。`schmidl_sync` 先消小数再搜整数，返回的波形不按 `d` 平移。`ofdm_tx` 不会插入这个前导。
- 和 `ofdm_sync` 的关系：那三个函数（`sto_estimate`、`cfo_estimate`、`ofdm_sync`）的计算公式和返回值都没有改。它们仍是循环前缀相关，捕获范围 ``(-1/2, 1/2)``，不估计整数 CFO，也不按定时下标搬移缓存。Schmidl–Cox 是旁边的另一套估计器。循环前缀会使 Schmidl 度量变成平台；无噪声时取最早的最大值。代码不做 Schmidl 与 Cox 原文里的平台平均，所以噪声把峰值挪到平台末端时，`d + n_cp` 会偏出符号。
- 版本：main 已是 0.8.3（SOVA、短 Turbo、BCH）。这一刀只加同步接口，已有函数的位置参数和返回值不变，记为 0.8.4。

## 0.8.3

软输出 Viterbi、短帧 Turbo，以及 Hamming(7,4) / 二进制 BCH(15,7) 的硬判决。不是 3GPP Turbo，不是 Reed–Solomon，没有接到 MCS / HARQ 上，也没有 BICM-ID。

- SOVA：`sova_decode`。终止的 ``K = 3``、(7, 5) 上，Hagenauer–Hoeher 的竞争路径可靠度，返回信息比特外信息 LLR（正值偏向 0）。分支度量是 BCJR 同一个 ``-\log\gamma``。尾比特固定为 0、终点固定在状态 0，边界和 `bcjr_decode` 相同；硬判决 `viterbi_decode` 不钉终点。无先验时硬判决是这条终止网格上的最大似然路径。幅度是近似，不要求等于 BCJR；无噪声和高信噪比下符号与 BCJR 一致。整块回溯，时间和存储随块长平方增长。某比特沿幸存路径从未与竞争路径分叉时，可靠度用 ``\sum|L|+\sum|L_a|`` 封顶。
- Turbo：`turbo_interleaver`、`turbo_encode`、`turbo_decode`。两路记忆 2 的递归系统卷积，八进制 (1, 5/7)，即 (7, 5) 的递归系统形式，不是前向 `conv_encode`，也不是 LTE 的记忆 4 分量码。信息比特只发一次；两路奇偶各带 2 个尾比特，尾比特不进交织器。码字长度 ``3k+8``，``2 \le k \le 1024``。``k`` 有不小于 2 的因子时，交织器就是按行写入、按列读出的块交织；质数长度用接近 ``\sqrt{k}`` 的步长。不是 3GPP QPP。译码是 log-BCJR 交换外信息，不是 SOVA，也不是前向码上的 `bcjr_decode`。`iterations` 为 1 到 32，默认 6。返回信息比特后验 LLR。名义码率 ``k/(3k+8)``，不是精确的 1/3。
- 有限域：`GF2m`，只做 ``m = 2, 3, 4``。本原多项式 ``x^2+x+1``、``x^3+x+1``、``x^4+x+1``。元素是比特掩码，bit 0 是常数项。加法是异或。
- Hamming：`hamming74_encode` / `hamming74_decode`。与 `hamming74_ldpc` 同一张校验矩阵、同一组 16 个码字，但是独立的伴随式硬判决，纠正 1 个错误。信息位在码字第 3、5、6、7 位。两个错误会落到另一个码字上。`hamming74_encode([1, 0, 1, 1])` 得到 `[0, 1, 1, 0, 0, 1, 1]`。
- BCH：`bch157_encode` / `bch157_decode`。GF(16) 上的二进制本原 BCH，``n=15``，``k=7``，设计距离 5，纠正 2 个错误。生成多项式 ``x^8+x^7+x^6+x^4+1``。比特从 ``x^{14}`` 排到 ``x^0``，系统位在前。Peterson 直接解加 Chien 搜索。定位多项式的根个数对不上时，返回收到的系统位，不翻转。3 个及以上错误不保证纠正。
- 版本：已有函数的位置参数和返回值不变。0.8.0 已经是 EVM、LFSR 扰码和 SVD 注水。0.8.1 与 0.8.2 留给并行的同步和 Jakes 两刀。这一刀是新的译码族，按 0.x 记为 0.8.3，而不是 0.8.1，也不升到 1.0。

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
