# 验证矩阵

蒙特卡洛对照用双侧 Clopper–Pearson 区间（Clopper 与 Pearson，Biometrika 26(4):404–413，1934）。``\alpha = 10^{-4}``，并要求足够的错误计数，理论值要落在区间内。种子是 `0x0c0de1162`，除非某个测试写明用了别的种子。公式实现在 `test/reference.jl`，不从包里导出。本页不新写公式。

编码和交织之后，比特错误和帧错误不是独立伯努利试验。区间仍按独立计数使用，和现有测试一样。

`demos/` 不在 `Pkg.test` 里。`demos/ber_16qam.jl` 画的是高信噪比近似，不是 Cho–Yoon 精确式。`demos/ber_comparison.jl` 只画仿真点，没有闭式曲线。`demos/validation_figure.jl` 的四条曲线是 v0.6 之前的图：AWGN 下的 BPSK、QPSK、16-QAM，以及平坦瑞利下的 BPSK。`demos/link_puncture_bicm.jl` 打印删余表的固定行、自适应、Chase 与 `harq=:ir` 的对照，以及 8-PSK 的 1 次与 4 次迭代。

统计规则和公式出处也写在 `test/reference.jl` 文件头。下表把每个测试文件都点到：`test/algebraic.jl`、`test/bicm.jl`、`test/contract.jl`、`test/evm.jl`、`test/exit.jl`、`test/interleave.jl`、`test/jakes_nakagami.jl`、`test/kronecker_mimo.jl`、`test/ldpc.jl`、`test/link.jl`、`test/lms_dfe.jl`、`test/psk_rice_mrc.jl`、`test/pulse_fir.jl`、`test/reference.jl`、`test/runtests.jl`、`test/scrambling.jl`、`test/sova.jl`、`test/sync.jl`、`test/turbo.jl`、`test/waterfill.jl`。

| 模块 | 量 | 闭式或参照 | 测试 | 演示 |
|---|---|---|---|---|
| Modulation，Metrics | AWGN 相干 BPSK 与 Gray QPSK 误比特率 | ``Q(\sqrt{2 E_b/N_0})``。Proakis 与 Salehi，《Digital Communications》第 5 版，§5.2。`ber_bpsk_awgn` / `ber_qpsk_awgn` | `test/runtests.jl`，公式在 `test/reference.jl` | `demos/validation_figure.jl`。`demos/ber_bpsk.jl` 用 ``\tfrac12\mathrm{erfc}(\sqrt{E_b/N_0})``，与这条 Q 函数相同。`demos/ber_qpsk.jl` 同样 |
| Modulation | 方形 Gray M-QAM 平均误比特率 | Cho 与 Yoon，IEEE Trans. Commun. 50(7):1074–1080，2002，式 (16)。`ber_square_mqam`。不是高信噪比近似 | `test/runtests.jl` | `demos/validation_figure.jl`（16-QAM 用精确式） |
| Channel | 理想 CSI、平坦瑞利、相干 BPSK | ``\tfrac12\bigl(1 - \sqrt{\bar\gamma/(1+\bar\gamma)}\bigr)``，``\mathbb{E}[\lvert h\rvert^2] = 1``。Goldsmith，《Wireless Communications》，2005，§6.3；Proakis 与 Salehi §13.3。`ber_bpsk_rayleigh` | `test/runtests.jl` | `demos/validation_figure.jl`；`demos/ber_bpsk_rayleigh.jl` |
| Modulation | Gray M-PSK 误比特率 | 相位误差扇区上的 Gray 汉明距离平均，`ber_mpsk_gray`。Proakis 与 Salehi §4.2、§5.2。最近邻 ``(1/\log_2 M)\operatorname{erfc}(\sqrt{\log_2 M\cdot E_b/N_0}\sin(\pi/M))`` 只在 ``M = 4`` 时等于 QPSK | `test/psk_rice_mrc.jl` | 无。`demos/ber_psk_rice_mrc.jl` 不打印这条积分 |
| Modulation | M-PSK 符号错误率 | Craig 积分 `ser_mpsk_craig`。Proakis 与 Salehi §5.2。与扇区概率交叉核对 | `test/psk_rice_mrc.jl` | 无 |
| Modulation | DBPSK 误比特率 | ``\tfrac12\exp(-E_b/N_0)``。Proakis 与 Salehi §4.5。参考符号不折进信息比特能量 | `test/psk_rice_mrc.jl` | `demos/ber_psk_rice_mrc.jl` |
| Modulation | Gray DQPSK 误比特率 | 一阶 Marcum Q，`ber_dqpsk`。Proakis 与 Salehi §4.5 | `test/psk_rice_mrc.jl` | 无。该演示不打印 Marcum Q |
| Channel | 理想 CSI、平坦 Rice、相干 BPSK | Craig 形式的 ``Q`` 对 Rician 矩生成函数积分，`ber_bpsk_rice`。Proakis 与 Salehi §13.3。`K` 是线性比。`K = 0` 退回瑞利，`K \to \infty` 退回 AWGN | `test/psk_rice_mrc.jl` | `demos/ber_psk_rice_mrc.jl` 只打印 `K = 0` |
| Channel | L 支路 i.i.d. 瑞利、理想 CSI、相干 BPSK 的 MRC | Proakis 与 Salehi §13.4 的求和式，`ber_bpsk_mrc_rayleigh`。``\bar\gamma_b`` 是每条支路的平均 ``E_b/N_0``。`L = 1` 退回单支路瑞利 | `test/psk_rice_mrc.jl` | `demos/ber_psk_rice_mrc.jl`（`L = 2`） |
| Channel | Clarke/Jakes 滞后相关 | 对照 ``J_0(2\pi f_d \tau)``。有限 DFT 格子和有限样本都留有误差。理想 CSI 相干 BPSK 的平均仍用瑞利闭式，因为边际是 ``\mathcal{CN}(0,1)`` | `test/jakes_nakagami.jl` | `demos/ber_jakes_nakagami.jl` |
| Channel | Nakagami-m、理想 CSI、相干 BPSK | Craig 形式对 Gamma 矩生成函数积分，`ber_bpsk_nakagami`。`m = 1` 退回瑞利。正整数 `m` 等于 `m` 支路瑞利 MRC，每条支路的平均是 ``\bar\gamma_b/m``。`m \to \infty` 退回 AWGN | `test/jakes_nakagami.jl` | `demos/ber_jakes_nakagami.jl` |
| Coding | (7, 5) 距离谱 | Lin 与 Costello，《Error Control Coding》第 2 版，(2,1,2) 转移函数系数。`distance_spectrum` | `test/runtests.jl` | 无 |
| Coding | 硬判决 Viterbi 的信息比特误码 | 截断联合界。Proakis 与 Salehi §8.2；Lin 与 Costello §12.3。区间上界不超过该界，并低于未编码 BPSK | `test/runtests.jl` | `demos/ber_conv_coding.jl` 只画仿真曲线，不画联合界 |
| Coding | 软输入 Viterbi 与 log-BCJR | 同一 ``E_b/N_0`` 下，Clopper–Pearson 区间低于硬判决，也低于未编码 BPSK。点估计低于软判决联合界 ``\sum c_d Q(\sqrt{2 d R E_b/N_0})``（``R = 1/2``），并大于该界的十分之一。Proakis 与 Salehi §8.2。BCJR：Bahl、Cocke、Jelinek 与 Raviv，IEEE Trans. Inf. Theory 20(2):284–287，1974 | `test/runtests.jl` | 无 |
| Modulation | BPSK 与 Gray QPSK 的比特 LLR | 闭式条件高斯。Proakis 与 Salehi §4.2、§5.2。``\sigma^2 = N_0/2``，与 `awgn` 相同 | `test/runtests.jl` | 无 |
| Coding | SOVA 相对 BCJR | 无闭式幅度。无噪声和高信噪比下符号与 `bcjr_decode` 一致，幅度不要求相等。Hagenauer 与 Hoeher，Proc. IEEE GLOBECOM，1989（`sova_decode` 的文档字符串） | `test/sova.jl` | 无 |
| Coding | 删余到名义码率 2/3 与 3/4 | 无新的自由距闭式。图案是 Yasuda 等人与 IEEE 802.11 用在 1/2 母码上的掩码（`puncture_mask` 的文档字符串）。4 dB 上 3/4 的区间高于母码；2/3 只要求点估计落在二者之间 | `test/interleave.jl` | 无 |
| Interleaving | 矩形块交织 | 无闭式误码。置换可逆。逐符号独立衰落上，块交织不改变每个符号的边际信噪比 | `test/interleave.jl` | `demos/ber_block_interleave.jl` |
| Channel | 块恒定瑞利 | 一块内一个 ``h \sim \mathcal{CN}(0,1)``，块之间独立。没有单独的闭式误码率测试 | `test/interleave.jl` | `demos/ber_block_interleave.jl` |
| Coding | 短帧 Turbo | 无闭式误码。手算码字、终止回到状态 0、短块上 log-BCJR 与穷举后验对照。不是 3GPP BLER 表。外信息分解见 Hagenauer、Offer 与 Papke，IEEE Trans. Inf. Theory 42(2):429–445，1996（`turbo_decode` 的文档字符串） | `test/turbo.jl` | 无 |
| Coding | ``\mathrm{GF}(2^m)``、Hamming(7,4)、BCH(15,7) | 域运算对照固定本原多项式。Hamming 伴随式用 Lin 与 Costello §3.2 的列顺序。BCH 设计距离 5，纠正 2 个错误。没有误码率曲线 | `test/algebraic.jl` | 无 |
| Coding | 短 Reed–Solomon | 生成根 ``\alpha^1, \ldots, \alpha^{n-k}`` 逐一核对；Singleton 界取等 ``d = n - k + 1``，RS(7,3) 穷举全部 512 个码字验证最小重量 5。``t`` 以内按位置和取值穷举：RS(7,3) 全部消息、RS(15,11) 与缩短 RS(10,6) 的固定消息集；RS(15,9) 的重量 3 用种子 `0x0c0de1162` 抽 20000 条。``t`` 以外只分两类：原样返回（错误全落在校验位时恰好就是发送消息），或误纠到接收字 ``t`` 球内的另一码字。缩短码里落在被缩短位置的错误定位根按失败处理。冻结的生成多项式与编解码向量。没有误码率曲线，不做蒙特卡洛统计判据。Reed 与 Solomon，J. SIAM 8(2):300–304，1960；Massey，IEEE Trans. Inf. Theory 15(1):122–127，1969；Chien，IEEE Trans. Inf. Theory 10(4):357–363，1964；Forney，IEEE Trans. Inf. Theory 11(4):549–557，1965；Lin 与 Costello 第 7 章 | `test/algebraic.jl` | `demos/rs_short.jl` |
| Coding | 短 LDPC | 树上的 sum-product 与逐比特 MAP 一致。`array_ldpc(3, 6, 11)`（``n = 66``，``k = 35``）在 2.5 dB 与 3.5 dB、各 1000 帧、种子 `0x0c0de1162`，归一化 min-sum 与 sum-product 的信息比特 BER 和 BLER 上界低于未编码 BPSK。构造与译码的文献在 `test/ldpc.jl` 文件头：Gallager 1962，Fan 2000，Lin 与 Costello §3.2，Kschischang、Frey 与 Loeliger 2001，Fossorier、Mihaljević 与 Imai 1999，Chen 与 Fossorier 2002。不是密度进化阈值。Hamming 校验图含 4 环，不是最大似然 | `test/ldpc.jl` | 无 |
| BICM | 8-PSK 集合分割上的迭代增益 | 无容量闭式。终止 (7, 5)，``E_s = 1``，``R = 400/804``，``E_b = 1/(3R) = 0.67``，max-log。AWGN、``E_b/N_0 = 6\,\mathrm{dB}``、80 帧：4 次迭代相对 1 次，信息比特 BER 和 BLER 的区间互不重叠。理想 CSI 逐符号瑞利在 8 dB、40 帧同样分开。标签是 Ungerboeck，IEEE Trans. Inf. Theory 28(1):55–67，1982，不是 Gray | `test/bicm.jl` | `demos/ber_bicm_id.jl` |
| EXIT | 软信息互信息与解映射 / 译码传递曲线 | BIAWGN 一致 LLR 的互信息对照数值积分（`biawgn_llr_mi`，无初等闭式）：``N = 200000`` 时 `:soft` 与 `:hist` 估计的偏差都在 0.01 以内，``J(2/\sigma)`` 多项式近似在 0.005 以内。恒等式见 Hagenauer，Proc. 12th EUSIPCO，2004，pp. 1541–1548；直方图形式见 ten Brink，IEEE Trans. Commun. 49(10):1727–1737，2001；``J(\sigma)`` 与逆近似见 ten Brink、Speidel 与 Yan，IEEE Trans. Commun. 52(4):670–678，2004，式 (9)–(11)；``L_e = L_{\mathrm{post}} - L_a`` 见 Hagenauer、Offer 与 Papke，IEEE Trans. Inf. Theory 42(2):429–445，1996。解调锚点：BPSK 与 Gray QPSK 在 ``I_A = 0`` 等于该积分，曲线平坦；6 dB 的 8-PSK 集合分割 ``T(0.98) - T(0) > 0.2``，Gray 接近平坦且起点更高。译码锚点：``E_b/N_0 = 10\,\mathrm{dB}`` 时 (7, 5) BCJR 的 ``T(0) > 0.99``；0 dB 时 ``T(0.99)`` 明显高于 ``T(0)``；阵列 LDPC 烟测。一致性反面教材固定在种子 ``0x0c0de1162`` 里：循环图上泛洪 sum-product 在 2.5 dB 给出过度自信的错误 LLR，`:soft` 远低于 `:hist` | `test/exit.jl` | `demos/exit_chart.jl` |
| Link | CRC-16 | CRC-16/XMODEM：ASCII `"123456789"` 的余数是 `0x31C3`。``g(x)`` 的倍式漏检。单比特翻转被检出 | `test/link.jl` | 无 |
| Link | 五档 MCS、Chase 合并、门限自适应 | 无闭式 BLER。无噪声时每一档吞吐等于该档频谱效率。``E_s/N_0 = 5\,\mathrm{dB}`` 的 16-QAM 卷积码上，3 次 Chase 相对单次传输的 CRC 失败率区间（100 个块）互不重叠，合并后吞吐更高。5 dB 的自适应（QPSK 阵列 LDPC）吞吐下界高于固定 BPSK，也高于固定 64-QAM；16 dB 的自适应（64-QAM）下界高于固定 BPSK。门限来自 `teaching_amc_thresholds` | `test/link.jl` | `demos/link_mcs_harq.jl` |
| Link | 删余卷积 MCS、同一删余字的 Chase、链路内 BICM | 无新的闭式 BLER。种子 `0x0c0de1162`，每点 80 块，Clopper–Pearson ``\alpha = 10^{-4}``。3 dB 自适应（QPSK 1/2）吞吐下界高于固定 BPSK，也高于固定 16-QAM；10 dB 自适应（16-QAM）下界高于固定 BPSK。2 dB 的 QPSK 3/4，3 次 Chase 相对单次的 CRC 失败率区间分开。8-PSK 集合分割、``E_s/N_0 = 6\,\mathrm{dB}``：链路内 4 次迭代相对 1 次，净荷比特错误率和 CRC 失败率的区间分开。5 dB 上码率 2/3 与 3/4 几乎打平，不拿来做判据。门限是 `teaching_punctured_amc_thresholds`，不能套到五档表 | `test/link.jl`；删余 BICM 的无噪声往返在 `test/bicm.jl` | `demos/link_puncture_bicm.jl` |
| Link | `harq=:ir` 教学增量冗余（互补删余图案、按母码位软合并） | 无闭式 BLER。种子 `0x0c0de1162`，QPSK 3/4，400 块，Clopper–Pearson ``\alpha = 10^{-4}``。``E_s/N_0 = 1\,\mathrm{dB}``、`max_tx=3`：IR 与 Chase 的 CRC 失败率区间互相重叠，各自点估计落在对方区间内（49/400 对 51/400），IR 信道使用更少、吞吐更高，即 IR 在该点与 Chase 打平且更省符号。2 dB、`max_tx=2`：Chase 的 NACK 点估计明显更低（96/400 对 136/400，此样本量下区间仍重叠），IR 吞吐更高。`harq=:chase` 与省略该关键字逐比特一致；`max_tx=1` 与码率 1/2 行上 `:ir` 退化为 `:chase`。参照 Hagenauer，IEEE Trans. Commun. 36(4):389–400，1988（速率兼容删余）与 Chase 1985 | `test/link.jl` | `demos/link_puncture_bicm.jl` |
| Equalization | 块 ZF / MMSE | 无噪声时反演 `apply_fir` 的全长卷积。文档字符串引用 Proakis 与 Salehi §9.4 与 Tse 与 Viswanath，《Fundamentals of Wireless Communication》，2005，§8.3。测试不对照误码闭式 | `test/runtests.jl`；`test/pulse_fir.jl` | `demos/ber_equalizer.jl`；`demos/pulse_fir_link.jl` |
| Equalization | 符号率 LMS | 无噪声时抽头朝这条回归量的批最小二乘解走，末段误差低于开头。同一条确定性噪声记录上，更大的稳定步长留下更大的末段均方误差。不对照 Widrow 超额 MSE 闭式 | `test/lms_dfe.jl` | 无 |
| Equalization | 有限长块 DFE | 判决正确时，ZF-DFE 把全列秩 FIR 恢复到舍入误差。第一个反馈符号被换成错误值时，后面的硬判决可以成串翻掉。无闭式 | `test/lms_dfe.jl` | 无 |
| Pulse，Visualization | 升余弦与根升余弦，匹配滤波后抽样 | 截断脉冲的和与能量。形状见 Proakis 与 Salehi §9.2（`raised_cosine` 的文档字符串）。无多径、无噪声时，根升余弦级联在符号点上贴近原符号。冲激脉冲走 `pulse_shape`、样本率 FIR、`matched_filter`、`downsample_symbols` 之后，ZF 与 `noise_var = 0` 的 MMSE 恢复符号 | `test/runtests.jl`；`test/pulse_fir.jl` | `demos/pulse_fir_link.jl`；`demos/viz_demo.jl` |
| OFDM | 循环前缀收发、已知 CSI / LS / 逐导频收缩 | 无闭式误码。无噪声时完美恢复。4 dB 与 18 dB 的 Clopper–Pearson 区间互不重叠。定时和小数倍频偏是 van de Beek、Sandell 与 Börjesson，IEEE Trans. Signal Process. 45(7):1800–1805，1997，捕获范围 ``(-1/2, 1/2)`` | `test/runtests.jl` | `demos/ofdm_ber.jl`（仿真曲线，无闭式） |
| MIMO | i.i.d. 瑞利上的 ZF 硬判决 | 无闭式误码曲线。无噪声时恢复符号。噪声方差 1 与 ``10^{-3}`` 的区间互不重叠。信道条目是 ``\mathcal{CN}(0,1)``，Tse 与 Viswanath §7.1（`mimo_channel` 的文档字符串） | `test/runtests.jl` | `demos/ber_mimo.jl`（仿真曲线，无闭式） |
| MIMO | 平坦 Kronecker 空间相关 | 无新的误码闭式。``H = R_r^{1/2} H_w R_t^{1/2}``，``\mathbb{E}[H_{pq} H_{rs}^{*}] = (R_r)_{pr}(R_t)_{sq}``。单位阵以及 Loyka ``\rho = 0`` 与 `mimo_channel` 逐元素相同。``N = 20000`` 的样本矩对照该乘积，容差 0.06。``\rho = 1`` 的发射相关使各列相等。无噪声时 ML / ZF / MMSE / SIC 恢复符号。``\texttt{noise\_var} = 0.05`` 时，``\rho = 0.95`` 的 ZF 误码区间高于 ``\rho = 0``。定义在 `kronecker_channel` 的文档字符串：Kermoal 等，IEEE JSAC 20(6)，2002，附录的乘积结构；指数矩阵是 Loyka，IEEE Commun. Lett. 5(9)，2001 | `test/kronecker_mimo.jl` | `demos/kronecker_mimo.jl`（打印相关矩阵、样本矩和两条 ZF 误码，无图） |
| MIMO | Alamouti 与 SVD 预编码 | 无噪声时恢复。`svd_precoder` 返回右奇异向量，不分配功率 | `test/runtests.jl` | 无 |
| MIMO | 总功率注水 | Tse 与 Viswanath §7.1.1：``p_i = (\mu - N_0/\lambda_i^2)_+``。功率和等于预算，正功率模式满足同一水位。低信噪比集中到最强模式，高信噪比在正奇异值上接近均匀 | `test/waterfill.jl` | 无 |
| Scrambling | 加性 LFSR | 对照文档字符串里的线性递推，不是另一份移位实现。具名多项式的周期是 ``2^n-1``。序列不按 ITU-T O.150 取反 | `test/scrambling.jl` | 无 |
| Sync | Gardner、Costas、Schmidl–Cox | 不跑误码率曲线。Gardner 核对 S 曲线符号和开环过零（Gardner，IEEE Trans. Commun. 34(5):423–429，1986）。Costas 核对鉴相器符号和锁定后的 QPSK 硬判决。Schmidl–Cox 在无噪声前导上核对两半相等、小数 CFO 和偶数倍整数 CFO（Schmidl 与 Cox，IEEE Trans. Commun. 45(12):1613–1621，1997）。``\varepsilon = 0.7`` 时 van de Beek 仍报成接近 ``-0.3`` | `test/sync.jl` | `demos/sync_timing.jl` |
| Metrics | RMS EVM | 默认百分数 ``100\sqrt{\sum\lvert e\rvert^2/\sum\lvert s\rvert^2}``。dB 是这个比的 ``20\log_{10}``，10% 为 ``-20\,\mathrm{dB}``。不做幅度、相位或直流拟合 | `test/evm.jl` | 无 |
| Metrics | 误码率、符号错误率、帧错误率 | 按元素或按块计数。没有信道闭式 | `test/runtests.jl`；`test/interleave.jl` | 无 |
| Visualization | 作图扩展 | 未 `using Plots` 时 `plot_constellation` 抛错；加载后返回 `Plots.Plot` | `test/runtests.jl` | `demos/viz_demo.jl` |
| （包） | 导出清单 | 不是信道公式。`names(CommToolkit)` 与 [公开接口](api.md) 一致，子模块导出只多一个不导出的 `modulation` | `test/contract.jl` | 无 |

## 闭式摘录

下列句子里，Kronecker 那一条是 0.9.3 加上的，定义在 `kronecker_channel` 的文档字符串，不在 `test/reference.jl`。1.0.0 不新增闭式，也不改这张表。1.1.0 只在矩阵里增加短 Reed–Solomon 一行，它没有误码闭式，本节的公式清单不变。1.1.1 增加一条数值积分对照（不是闭式）：BIAWGN 一致 LLR 的互信息 ``I = 1 - \mathbb{E}[\log_2(1+e^{-L})]``，``L = 2y/\sigma^2``，``y \sim \mathcal{N}(1, \sigma^2)``，Simpson 积分在 `test/reference.jl` 的 `biawgn_llr_mi`，它是 ``J(2/\sigma)`` 的精确值（Hagenauer，Proc. 12th EUSIPCO，2004）。1.1.2 增加 `harq=:ir` 教学增量冗余一行：没有闭式 BLER，是种子 `0x0c0de1162` 下 IR 对 Chase 的两个蒙特卡洛对照点，公式清单不变。其余与 v0.9.0 的 README 验证段相同，公式仍以 `test/reference.jl` 为准。

- AWGN 下 BPSK 与 Gray QPSK：``Q(\sqrt{2 E_b/N_0})``。Proakis 与 Salehi，第 5 版，§5.2。
- 方形 Gray M-QAM 的精确平均误比特率：Cho 与 Yoon，2002，式 (16)。图中 16-QAM 用的是这个精确式，不是高信噪比近似。
- 理想 CSI 下平坦瑞利相干 BPSK：``\tfrac12\bigl(1 - \sqrt{\bar\gamma/(1+\bar\gamma)}\bigr)``，``\mathbb{E}[\lvert h\rvert^2] = 1``。Goldsmith，2005，§6.3；Proakis 与 Salehi §13.3。
- Gray M-PSK：相位误差扇区上的 Gray 汉明距离平均（Proakis 与 Salehi §4.2、§5.2）。符号错误率同时用 Craig 积分交叉核对。最近邻近似在 ``M = 4`` 时等于 QPSK，``M > 4`` 不是精确式。蒙特卡洛对照的是扇区积分。
- DBPSK：``\tfrac12\exp(-E_b/N_0)``。Gray DQPSK：一阶 Marcum Q 比特错误率（Proakis 与 Salehi §4.5）。差分检测的参考符号是前导，公式不把它折进信息比特能量。
- 理想 CSI 相干 BPSK、平坦 Rice：Craig 形式的 ``Q`` 对 Rician 矩生成函数积分。`K = 0` 退回瑞利闭式，`K \to \infty` 退回 AWGN。`K` 是线性比。
- L 支路 i.i.d. 瑞利、理想 CSI、相干 BPSK 的 MRC：Proakis 与 Salehi §13.4 的求和式。``\bar\gamma_b`` 是每条支路的平均 ``E_b/N_0``。
- Clarke/Jakes：滞后相关对照 ``J_0(2\pi f_d \tau)``。理想 CSI 相干 BPSK 的平均仍用瑞利闭式。
- Nakagami-m、理想 CSI、相干 BPSK：Craig 形式对 Gamma 矩生成函数积分。`m = 1` 退回瑞利。正整数 `m` 等于 `m` 支路瑞利 MRC，每条支路的平均是 ``\bar\gamma_b/m``。`m \to \infty` 在理想 CSI 下退回 AWGN。
- 平坦 Kronecker：没有新的误码闭式。样本矩对照 ``\mathbb{E}[H_{pq} H_{rs}^{*}] = (R_r)_{pr}(R_t)_{sq}``，定义在 `kronecker_channel` 的文档字符串。单位阵与 ``\rho = 0`` 退回 i.i.d. 瑞利。

![BPSK、QPSK、16-QAM 在 AWGN 下以及 BPSK 在平坦瑞利下的蒙特卡洛误码率与闭式曲线](assets/ber_validation.png)

上图由 `demos/validation_figure.jl` 生成，种子 `0x0c0de1162`。
