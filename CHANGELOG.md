# 变更记录

## 1.2.1

教学 Arikan 极化码：短帧构造、系统/非系统编码与串行抵消（SC）译码。坐在 1.2.0（短 RS、EXIT、IR-HARQ、教学 CDL、教学迷你 xApp）之上，是 1.x 内的增加：稳定导出清单增加 4 个名字（`PolarCode`、`polar_code`、`polar_encode`、`polar_sc_decode`），已有名字的位置参数、返回值与语义都不变，`simulate_link`、xApp、CDL、EXIT、RS 一个字符都没有改，v1.0.0 的全部默认行为不变。这是教学 Arikan polar SC，不是 3GPP NR polar：没有 CRC 级联（也就没有 CRC-aided 译码），没有 CA-SCL 或任何列表译码，没有速率匹配（``N = 2^n`` 原样发送，无删余 / 缩短 / 重复），没有 PW 或 QPP 交织，没有 3GPP 母码，也没有 38.212 的 NR 可靠度序列。有意留在范围之外的仍是完整几何 TR 38.901 CDL（角度、射线、天线几何、空间一致性、MIMO 全矩阵）、NR LDPC、NR polar、真 O-RAN。

- 构造：`polar_code(n, k; design_snr_db=0.0)` 返回 `PolarCode`（字段 `n`、`k`、`info_pos`、`frozen`、`design_snr_db`、`z`）。码长 ``N = 2^n``、``2 \le N \le 256``（教学上限），``1 \le k \le N``。信息集由 Bhattacharyya 参数递归选出：根信道取二元输入 AWGN 在设计 ``E_s/N_0`` 下的 ``Z = e^{-E_s/N_0}``，按 ``Z^- = 2Z - Z^2``、``Z^+ = Z^2`` 展开到 ``N`` 个比特信道（Arikan，IEEE Trans. Inf. Theory 55(7):3051–3073，2009，Prop. 5），取 ``Z`` 最小的 ``k`` 个索引（1-based 升序，稳定排序）。递归对二元擦除信道精确，对 AWGN 是启发式可靠度排序——不是密度进化，更不是 NR 可靠度序列；`design_snr_db` 只是设计点，信息集不随实际信道自适应。
- 编码：`polar_encode(code, bits; systematic=false)`。非系统编码是 ``x = u F^{\otimes n}`` 的自然序蝴蝶变换（``F = [1 0; 1 1]``），与 Arikan 论文的 ``B_N F^{\otimes n}`` 只差码字坐标的比特反转置换，码集合相同；冻结位填 0。`systematic=true` 走 Arikan 2011（IEEE Trans. Inf. Theory 57(5):3195–3200）思想的三角回代：``F^{\otimes n}`` 单位对角下三角、信息位主子式必可逆，保证 ``x_A = \texttt{bits}``。两种映射的码字集合相同，系统编码不是近似。
- 译码：`polar_sc_decode(code, llr; systematic=false)`。逐位串行抵消，f 运算用精确 box-plus（稳定形式，复用 LDPC 的内部函数），g 运算用 ``(1-2\hat u)\alpha_i + \alpha_{i+m}``，冻结位强制判 0，末位 LLR 平局判 0。LLR 约定与全包一致（正值表示比特 0，与 `demodulate_soft` 相同）；只接受有限浮点，整数硬比特被拒绝。复杂度每帧 ``O(N \log N)``。不是最大似然：早判错的比特沿格图传播；没有 CRC 级联，也没有 SCL / SCF / BP 译码器。
- 诚实的对照（`test/polar.jl`）：Bhattacharyya 递归逐步手算核对；``N \in \{2, 4, 8, 16, 32\}`` 多个设计点的信息集冻结为固定向量；蝴蝶变换的对合性、线性性与 ``N = 2, 4`` 手算；``N = 8, k = 3`` 穷举全部 256 条消息验证两种编码的码字集合相同、系统性逐点成立；最小距离穷举等于 ``2^{\min \operatorname{popcount}(i-1)}``；``N = 8`` 全部 ``k`` 的无噪声 SC 穷举往返。两条统计锚点都在种子 `0x0c0de1162`、Clopper–Pearson ``\alpha = 10^{-4}`` 下：``k = 1`` 时 ``(N, 1)`` 极化码就是重复码、SC 退化为 LLR 求和，``N = 2``（2 dB，8000 帧）与 ``N = 4``（3 dB，12000 帧）的误码率回到 `ber_bpsk_awgn` 的 ``Q(\sqrt{2\gamma_b})``；``N = 64, k = 32``（`design_snr_db = 2.0`，每点 1000 帧）在 2.0 与 3.0 dB 上的信息比特 BER 与 BLER（``1-(1-p_b)^{32}``）区间上界低于未编码 BPSK，且两点区间互不重叠。演示 `demos/polar_sc.jl` 打印构造信息集、无噪声往返、重复码锚点闭式行与小 BER 扫描，不加载 Plots。
- 文档与契约：`docs/src/api.md` 的 Coding 稳定导出加 4 个名字、源文件行注明 `src/polar.jl`；`docs/src/index.md` 的 `@docs` 块同步；边界、已知限制、验证矩阵与 README 同步，「教学 Arikan polar SC，非 NR polar」写进边界与 README；`test/runtests.jl` 与 `test/contract.jl` 的版本断言升到 `v"1.2.1"`。没有触碰 `fpga/ax7020`，也没有触碰 1.1.0 的 RS 面、1.1.1 的 EXIT、1.1.2 的 `harq` 关键字、1.1.3 的 `src/cdl.jl`、1.2.0 的 `src/xapp.jl`。

## 1.2.0

路线 C：教学多用户 KPI→动作迷你 xApp。坐在 1.1.3（短 RS、EXIT、IR-HARQ、教学 CDL）之上，是 1.x 内的增加：稳定导出清单增加 10 个名字，已有名字的位置参数、返回值与语义都不变，`simulate_link` 一个字符都没有改，v1.0.0 的全部默认行为不变。这是教学迷你 xApp，不是 O-RAN：没有 E2 接口，没有 near-RT RIC，没有 E2SM，没有前传 / eCPRI 或功能切分，没有 NR 调度器，没有 NR CQI，也没有 O-RAN 合规声称。有意留在范围之外的仍是完整几何 TR 38.901 CDL（角度、射线、天线几何、空间一致性、MIMO 全矩阵）、NR LDPC、真 O-RAN。

- 模型：`XAppUE(name, esn0_db)` 描述一个 UE（名字加基准 ``E_s/N_0``）。`simulate_xapp(ues, policy, n_rounds; blocks_per_round=30, table=teaching_punctured_mcs_table(), max_tx=1, harq=:chase, approx=:exact, channel=:awgn, power_budget=1.0)` 每轮让策略读累计 `UEKPI`（块数、ACK、交付块、交付比特、信道使用，以及由此算出的 BLER 与吞吐），给每个 UE 指派一个 `XAppAction`——MCS 行、调度份额、功率偏移——再把该轮块预算按份额以最大余数法分块（并列取下标小者），各 UE 在 `基准 + 偏移`（`:rayleigh` 时再加该轮衰落）上跑已有的 `simulate_link`。多用户就是独立单链路加两个共享约束：每轮块预算和份额加权平均功率预算 ``\sum_u \texttt{share}_u 10^{\texttt{offset}_u/10} \le \texttt{power\_budget}``；没有 UE 间干扰模型，没有小区几何。`harq=:ir`、`max_tx`、`approx` 原样转发给 `simulate_link`，检查语义不变。
- 策略：`FixedXAppPolicy(; mcs=1, shares=nothing, power_offset_db=0.0)` 是固定分配基线，不读 KPI，不做功率控制，超预算即 `ArgumentError`。`RuleXAppPolicy(n_ues; target_bler=0.1, offsets_db=[-3, 0, 3], thresholds_db=teaching_punctured_amc_thresholds())` 是教学规则表：份额 ``\propto 1/(\text{已交付比特}+1)``（Tse 与 Viswanath §6.3 的比例公平风味，不是 3GPP 调度器）；`observe_kpi!` 按累计 BLER 是否超过 `target_bler` 或低于 `target_bler/10` 每次走一档并夹到梯子两端；MCS 用 `select_mcs` 读基准 ``E_s/N_0`` 加投影后偏移（静态 AWGN 上是 genie，`:rayleigh` 下衰落不进决策）。预算投影在超限时把在用最高档降一档，最低档仍不可行才报错。KPI 从第 1 轮累计，没有滑动窗，这是写明的教学简化。
- 报告：`XAppReport` 汇总每 UE 的 `UEKPI`、系统吞吐（bit/use，符号率 1、奈奎斯特带宽 1 Hz 时即聚合频谱效率）、系统 BLER，以及按每 UE 交付净荷比特算的 `jain_index`（Jain、Chiu 与 Hawe，DEC-TR-301，1984；全零返回 `NaN`）。`:rayleigh` 是理想 CSI 的轮粒度块衰落，每轮每 UE 一个 ``\mathcal{CN}(0,1)``，dB 夹到不低于 −60，与 `rayleigh_channel(..., true)` 同一条件等价，不是新 PHY；随机流顺序写在 `simulate_xapp` 的文档字符串里，同一颗种子两种策略看到同一条流。
- 诚实的对照（`test/xapp.jl`）：构造与参数校验、Jain 定义值、规则表在合成 KPI 上的份额排序 / 走档夹持 / 预算投影与 MCS 跟随、无噪声闭环的精确频谱效率与 Jain 值、`:awgn` / `:rayleigh` 同种子逐字段复现，外加一组种子固定的基线对照——3 个 UE（0、3、6 dB），20 轮 × 30 块，固定 QPSK-1/2 均匀份额对规则表，两臂各 600 块且各臂 NACK ≥ 30，交付块数的 Clopper–Pearson 区间（``\alpha = 10^{-4}``）分开，系统吞吐与 Jain 是同一流上的点估计对照。计数不锁死。演示 `demos/xapp_multiuser.jl` 打印对照表，不加载 Plots。
- 文档与契约：`docs/src/api.md` 增加 `XApp` 一节与 `src/xapp.jl` 一行；`docs/src/index.md` 增加 `@docs` 块；边界、已知限制、验证矩阵与 README 同步，路线 C 收入范围；`docs/make.jl` 的 `modules` 加 `CommToolkit.XApp`；`test/runtests.jl` 与 `test/contract.jl` 的版本断言升到 `v"1.2.0"`，子模块清单加 `:XApp`。没有触碰 `fpga/ax7020`，也没有触碰 1.1.0 的 RS 面、1.1.1 的 EXIT、1.1.2 的 `harq` 关键字、1.1.3 的 `src/cdl.jl`。

## 1.1.3

教学簇延迟线（Clustered Delay Line，CDL）：TR 38.901 CDL-A/B/C 延迟–功率剖面上的可复现多径，供 OFDM 与链路 demo 使用。坐在 1.1.2（短 RS、EXIT、IR-HARQ）之上，是 1.x 内的增加：稳定导出清单增加 6 个名字，已有名字的位置参数、返回值与语义都不变，v1.0.0 的全部默认行为不变。这是教学 CDL，不是完整 3GPP 信道实现，也不声称 NR 合规：没有角度、没有每簇 20 射线、没有天线几何、没有空间一致性、没有 MIMO 全矩阵。有意留在范围之外的仍是完整几何 TR 38.901 CDL、路线 C 的多用户 xApp、NR LDPC。

- 剖面：`CDLProfile(delays, powers_db; name=:custom)` 与 `cdl_profile(:A | :B | :C)`。表值取自 3GPP TR 38.901 V16.1.0 表 7.7.1-1/2/3 的归一化延迟与功率两列（23、23、24 簇，NLOS）；归一化延迟的功率加权 RMS 延迟扩展为 1（§7.7.3），测试逐值冻结。LOS 的 CDL-D/E 拒绝。离开角 / 到达角、簇角扩展、XPR 不进剖面。
- 抽头：`cdl_taps(profile, n, Ts; delay_spread=100e-9, fd=0.0, n_taps=nothing)` 返回 `CDLTaps(delays, gains, powers)`。延迟按式 (7.7-1) 线性缩放（`delay_spread` 单位秒，默认 100 ns 是表 7.7.3-1 的 nominal 一档；`delay_spread = 1.0` 时剖面延迟直接读成秒），按 `round` 量化到样本，落进同一样本的簇按功率相加合并。抽头功率归一到 ``\sum p_i = 1``。每抽头一条独立的 `jakes_fading(n, fd, Ts)`：`fd = 0` 是单次静态 drop（一列增益，每抽头一对 `randn`）；`fd > 0` 是逐抽头 Clarke 谱（`fd * Ts < 1/2`），不是 38.901 的射线求和多普勒。`n_taps` 抽稀保留最强抽头并重新归一，是教学便利，不是标准内容。同一 `Random.seed!` 复现同一组抽头；随机流顺序写明在文档字符串里。
- 卷积与信道：`apply_cdl(signal, taps)` 是按输入样本索引增益的时变线性卷积，输出全长 ``n + d_{\max}``；静态单列 drop 与 `apply_fir`、`TDLChannel` + `apply_multipath` 逐样本一致。`cdl_channel(symbols, snr_db, bits_per_symbol; Ts, profile, ...)` 在卷积后按 `awgn` 的 ``E_b/N_0`` 约定加噪（含卷积尾部，与 `fir_channel` 相同），`Ts` 是必填关键字；功率归一使 `snr_db` 仍是平均 ``E_b/N_0``。
- OFDM 对接：静态 drop 用 `TDLChannel(taps.gains[:, 1], taps.delays)` 直接进 `apply_channel` / `ofdm_rx_known`。时变增益没有对应 `TDLChannel` 接口。
- 诚实的对照（`test/cdl.jl`）：表值与单位 RMS 延迟扩展冻结；量化合并向量冻结；`fd = 0` 抽取与 `jakes_fading(1, 0, Ts)` 逐值一致；`fd > 0` 单条记录的时间平均功率对照剖面功率（容差写明）；静态 drop 与 FIR / TDL 逐样本一致；内部样本输出功率 ``\mathbb{E}|y|^2 = 1``；唯一闭式锚点是全部簇并到 0 号抽头的平坦退化——合成系数 ``\mathcal{CN}(0,1)``，理想 CSI 相干 BPSK 回到 `ber_bpsk_rayleigh`（种子 `0x0c0de1162`，逐符号独立抽取用 Clopper–Pearson ``\alpha = 10^{-4}``，相关记录用条件平均对照，与 Jakes 边际测试同一写法）。色散情形没有黄金 BER 曲线；OFDM 对照是循环前缀内的无噪声完美恢复烟测。演示 `demos/cdl_ofdm.jl` 打印剖面、量化抽头、静态 drop 的 OFDM 烟测和平坦退化对照，不加载 Plots。
- 文档与契约：`docs/src/api.md` 的 Channel 一节增加 6 个名字，`Channel` 源文件行注明包含 `src/cdl.jl`；`docs/src/index.md` 的 `@docs`、边界、已知限制、验证矩阵与 README 同步；`test/runtests.jl` 与 `test/contract.jl` 的版本断言升到 `v"1.1.3"`。没有触碰 `fpga/ax7020`，也没有触碰 1.1.0 的 RS 面、1.1.1 的 EXIT、1.1.2 的 `harq` 关键字。

## 1.1.2

`simulate_link` 增加关键字 `harq`：`:chase`（默认，与 1.0.0 逐比特一致）或 `:ir`（教学增量冗余）。坐在 1.1.1（短 RS、EXIT）之上，是 1.x 内的增加：稳定导出清单不增删名字，已有位置参数和返回值不变，v1.0.0 的全部默认行为不变。不是 NR 速率匹配，不是 RV0–RV3，也不是 NR HARQ 实体。1.1.0 与 1.1.1 已被并行的短 RS 与 EXIT 两刀占用，这一刀记为补丁 1.1.2，不升主版本。

- 增量冗余：`simulate_link(...; harq=:ir)`，只接 `:viterbi` 接收机与卷积行。删余行上第 1 次发送 `puncture_mask` 图案 A，第 2 次发补图案（A 删掉的母码比特），更长的循环按 A、B 交替。比特 LLR 按母码位相加（擦除位仍是 LLR 0，复用 `depuncture`），所以第 2 次之后译码器见过每个母码比特一次，等效码率从删余码率降到母码 1/2。补字比 A 字短，`:ir` 重传花的信道使用更少，`n_channel_uses` 按实际发送计数。CRC 仍决定 ACK，`max_tx` 语义不变。码率 1/2 的行没有可补的比特，`:ir` 退化为 `:chase`；`max_tx=1` 只发图案 A，与 `:chase` 一致。LDPC 行、`receiver=:bicm`、补字长度不整除每符号比特数，都在仿真开始前报 `ArgumentError`。
- 诚实的对照（静态 AWGN，``K = 3`` 母码，种子 `0x0c0de1162`，QPSK 3/4，400 块，Clopper–Pearson ``\alpha = 10^{-4}``）：1 dB、`max_tx=3` 时 IR 与 Chase 的 CRC 失败率区间互相重叠（49/400 对 51/400），IR 少花信道使用、吞吐更高；2 dB、`max_tx=2` 时 Chase 的 NACK 点估计明显更低（96/400 对 136/400），因为 Chase 把删余字的能量翻倍，胜过把每个母码比特各看一次。这是小码在静态信道上的教科书式折中，不宣称衰落分集增益。参照 Hagenauer，IEEE Trans. Commun. 36(4):389–400，1988；Chase，IEEE Trans. Commun. 33(5):385–393，1985。
- 测试：`test/link.jl` 增加「HARQ mode keyword and IR fallbacks」与「Incremental redundancy versus Chase」。前者锁定 `:chase` 默认回归（显式 `harq=:chase` 与省略该关键字的报告逐字段一致）、`max_tx=1` 与码率 1/2 的退化、无噪声形状、自适应策略与 `:ir` 的组合，以及关键字参数检查。后者是上面的两个对照点。
- 文档与契约：`docs/src/api.md` 注明 1.1.2 只加关键字、不增删导出名字；`docs/src/limitations.md`、`docs/src/boundaries.md`、`docs/src/validation.md`、`docs/src/index.md` 与 README 改成与这条链路一致。`demos/link_puncture_bicm.jl` 增加 IR×3 列与 `max_tx=2` 的诚实角落。`test/runtests.jl` 与 `test/contract.jl` 的版本断言升到 `v"1.1.2"`。没有触碰 `fpga/ax7020`，也没有触碰 1.1.0 的 `src/algebraic.jl` RS 面与 1.1.1 的 `src/exit.jl`。

## 1.1.1

EXIT 图教学测量（软信息互信息）。这是 1.x 的加法：稳定面新增 `EXIT` 子模块的 8 个名字，1.0.0 与 1.1.0 已有名字的位置参数、返回值与语义都不变，`bicm_demap`、`bcjr_decode`、`ldpc_decode_details` 一个字符都没有改。教学 EXIT / 直方图互信息，不是 NR LDPC 的 EXIT 生产工具：没有密度进化，没有次数分布优化，没有 3GPP 码或速率匹配，也没有迭代轨迹跟踪。有意留在范围之外的仍是 TR 38.901 的 CDL、路线 C 的多用户 xApp、NR LDPC。

- 互信息估计：`exit_mi_soft(llr, bits)` 是 Hagenauer（Proc. 12th EUSIPCO，Vienna，2004，pp. 1541–1548）的无直方图软估计 ``1 - \frac1N\sum\log_2(1+e^{-(1-2b)L})``，假定等概比特与一致（真后验）LLR；单类输入合法（全零字惯例），``\pm\infty`` 按定义处理。`exit_mi_hist(llr, bits; nbins, clip)` 是 ten Brink（IEEE Trans. Commun. 49(10):1727–1737，2001）的直方图估计，只做两类条件直方图，不需要一致性，但有分箱偏差，两类比特都必须出现。LLR 约定与全包一致：bit0 为正，外信息 = 后验 − 先验（Hagenauer、Offer 与 Papke，IEEE Trans. Inf. Theory 42(2):429–445，1996）。
- 高斯对称捷径：`exit_j(σ)` 与 `exit_jinv(I)` 是 ten Brink、Speidel 与 Yan（IEEE Trans. Commun. 52(4):670–678，2004，式 (9)–(11)）的多项式近似。这是写明的可选捷径，不是测量，也不是精确积分；文档字符串写清了近似量级。
- 传递曲线：`exit_demap_transfer(scheme; ia, esn0_db, labels, approx, ...)` 在单位能量 AWGN 上测 `bicm_demap` 的 ``I_E = T(I_A)``，先验是合成的一致高斯 LLR（EXIT 测量的标准假设），``I_A = 1`` 用 `sigma_max` 截成有限。`exit_decoder_transfer(; decoder=:bcjr | :ldpc, ebn0_db, ...)` 测信息比特外信息的传递特性：`:bcjr` 走 `conv_encode` + `bcjr_decode`（返回的已是外信息），`:ldpc` 走 `ldpc_decode_details`（`early_stop=false`，外信息 = 后验 − 先验，先验折进信道 LLR 入口）。SOVA 与 Turbo 不接：公开接口没有编码比特外信息。曲线装在 `EXITCurve(ia, ie)` 里。
- 随机流：传递曲线用模块内的 xorshift64* + Box–Muller 局部流（与 `test/bicm.jl` 的构造相同），`seed` 固定时跨 Julia 版本可复现，不消耗全局随机数。默认种子 `0x0c0de1162`。
- 作图：`plot_exit_chart`（`using Plots` 后由 `CommToolkitPlotsExt` 提供）画解调曲线与交换坐标轴的译码曲线；只画测到的曲线，不模拟轨迹。
- 对照：`test/exit.jl`。BIAWGN 一致 LLR 的互信息对照数值积分（`test/reference.jl` 的 `biawgn_llr_mi`，Simpson，无初等闭式）：``N = 200000`` 时两种估计偏差 ≤ 0.01，``J(2/\sigma)`` 近似偏差 ≤ 0.005。解调锚点：BPSK 与 Gray QPSK 在 ``I_A = 0`` 等于该积分、曲线平坦；6 dB 的 8-PSK 集合分割 ``T(0.98) - T(0) > 0.2``，Gray 接近平坦。译码锚点：10 dB 时 (7, 5) BCJR 的 ``T(0) > 0.99``；0 dB 时 ``T(0.99)`` 明显高于 ``T(0)``；阵列 LDPC 烟测。一致性反面教材固定在种子里：循环图上泛洪 sum-product 在 2.5 dB 的过度自信错误 LLR 使 `:soft` 远低于 `:hist`。打印对照在 `demos/exit_chart.jl`，不加载 Plots。
- 契约：`docs/src/api.md` 增加 EXIT 一节与 `src/exit.jl` 一行；`test/contract.jl` 的子模块清单加 `:EXIT`，版本钉到 1.1.1。边界、已知限制、验证矩阵、README 与之一致。

## 1.1.0

短 Reed–Solomon。按 `docs/src/api.md` 的承诺只增名字：不改名、不删除，已有函数的位置参数和返回值不变，按 SemVer 记次版本。这是 1.0.0 有意留在外面的那一块；EXIT 图、TR 38.901 的 CDL、路线 C 的多用户 xApp、NR LDPC 仍在范围外。

- 构造：`rs_code(m, k)`。``\mathrm{GF}(2^m)`` 上的窄义本原 RS(``n, k``)，``n = 2^m - 1``，``m \in \{2,3,4\}``，即 ``n \in \{3, 7, 15\}``。生成多项式 ``g(x) = \prod_{j=1}^{n-k} (x + \alpha^j)``，``\alpha = x``，与各域的固定本原多项式一致。``d = n - k + 1``（MDS），``t = \lfloor (n-k)/2 \rfloor``。`rs_shorten(C, s)` 缩短成 ``(n-s, k-s)``：同一生成多项式，``t`` 不变，可复合。
- 接口：五个新导出 `RSCode`、`rs_code`、`rs_shorten`、`rs_encode`、`rs_decode`。系统编码，码元从高次到低次，与 `bch157_encode` 同序。译码是 Berlekamp–Massey、Chien 搜索和 Forney 公式的硬判决有界距离译码，失败约定与 `bch157_decode` 相同：没有失败标志；定位器次数超 ``t``、根数对不上、根落在被缩短的位置、或纠完伴随式仍非零时，原样返回 ``k`` 个系统符号；超过 ``t`` 可能误纠到接收字 ``t`` 球内的另一码字。不是纠删译码，不是 Guruswami–Sudan 列表译码，不是 Kötter–Vardy 软代数译码，不是 CCSDS 的 (255, 223) 或 (255, 239)，仍没有 ``\mathrm{GF}(2^8)``。
- 对照：`test/algebraic.jl` 新增四节。生成根逐一核对；RS(15,11) 与 RS(15,9) 生成多项式冻结向量（Lin 与 Costello 第 7 章例，``x^4+\alpha^{13}x^3+\alpha^6x^2+\alpha^3x+\alpha^{10}`` 与 ``x^6+\alpha^{10}x^5+\alpha^{14}x^4+\alpha^4x^3+\alpha^6x^2+\alpha^9x+\alpha^6``）；RS(7,3) 全部 512 条消息往返、最小重量 5 穷举（Singleton 取等）；``t`` 以内按位置和取值穷举（RS(7,3) 全部消息，RS(15,11) 与缩短 RS(10,6) 固定消息集，RS(15,9) 重量 3 按种子 `0x0c0de1162` 抽 20000 条）；``t`` 以外分类为原样返回或 ``t`` 球内误纠；冻结编解码向量。无误码率闭式，不加蒙特卡洛统计判据。
- 契约：`docs/src/api.md` 增加五个名字，版本 1.1.0。`test/runtests.jl` 与 `test/contract.jl` 的版本号同步。边界、已知限制、验证矩阵与 README 同步。

## 1.0.0

公开接口按 SemVer 冻成 1.0。这一版不增加算法，不改已有函数的位置参数和返回值，也不改 `docs/src/api.md` 里已经列出的名字。`using CommToolkit` 带入的名字就是稳定面。删除、改名，或改变位置参数和已有返回值的含义，必须升主版本号，并写进本文件。1.x 可以增加名字；增加本身不升主版本。

收进的是已经落地的三块，不在这一版重做：

- 0.9.1 工程契约：`docs/src/api.md` 的导出清单、`test/contract.jl`、验证矩阵、边界页。README 是导航。
- 0.9.2 删余码率与 BICM 接入 MCS / Chase。五档表、默认 `receiver=:viterbi` 和删余表都不变。重传仍是同一个删余字。
- 0.9.3 平坦 Kronecker 空间相关。`mimo_channel` 与 ML / ZF / MMSE / SIC 的位置参数、返回值不变。不是 Weichselberger，也不是 TR 38.901 的 CDL 或 UMi。

范围内的是教学 MCS / HARQ（删余和 BICM 路径）以及平坦 Kronecker MIMO。有意留在 1.0 之外：短 Reed–Solomon（代数码只有 BCH(15,7)）、EXIT 图、TR 38.901 的 CDL、路线 C 的多用户 xApp、NR LDPC。不是 NR，不是 O-RAN，也不是生产协议栈。

- 版本：`Project.toml` 为 `1.0.0`。`test/runtests.jl` 与 `test/contract.jl` 都要求 `pkgversion(CommToolkit) == v"1.0.0"`。Documenter 侧栏和 `objects.inv` 清单版本都取 `pkgversion`，即 `1.0.0`。

## 0.9.3

平坦 Kronecker 空间相关 MIMO。不是 Weichselberger 的两端联合耦合，不是 TR 38.901 的 CDL 或 UMi，没有多普勒，也没有接到 OFDM 或 MCS / HARQ 上。`mimo_channel` 和 ML / ZF / MMSE / SIC 的位置参数、返回值都不变。0.9.2 已经把删余码率和 BICM 接进 MCS / Chase，这一版不动那条路径。记为 0.9.3，不升到 1.0.0，也不记成 0.10.0。

- 信道：`kronecker_channel`。``H = R_r^{1/2} H_w R_t^{1/2}``。``H_w`` 与 `mimo_channel` 是同一套 ``\mathcal{CN}(0,1)``：先全部实部，再全部虚部。``R^{1/2}`` 是厄米半正定平方根。可分离相关是 ``\mathbb{E}[H_{pq} H_{rs}^{*}] = (R_r)_{pr}(R_t)_{sq}``（``p,r`` 为接收天线，``q,s`` 为发射天线）。按列堆叠时协方差是 ``R_t^{\mathsf{T}} \otimes R_r``。这是 Kermoal、Schumacher、Pedersen、Mogensen 与 Frederiksen（IEEE JSAC，20(6):1211–1226，2002）附录里「两个信道系数的相关等于两端天线相关之积」的平方根生成。实相关矩阵上，乘积与 CommPy `MIMOFlatChannel` 的 `sqrtm(R_rx) @ H_w @ sqrtm(R_tx)` 相同；复厄米因子的右乘是厄米平方根本身，不是它的非共轭转置。相关矩阵不重新归一。对角为 1 时，每个元素的平均功率仍是 1。``R - R^{\mathsf{H}}`` 大于 ``\max(1,\|R\|_\infty)`` 的 `1e-8` 倍会被拒绝；更小的负特征值夹到 0。两个整数参数的调用就是 `mimo_channel`。精确单位阵不乘 ``I``，随机数流与 i.i.d. 相同。校验失败发生在抽取之前。
- 指数相关：`exponential_correlation(n, ρ; spacing=1)`。Loyka，IEEE Commun. Lett.，5(9):369–371，2001。``\rho`` 是单位距离上的相关系数，``|\rho| \le 1``。`spacing` 是相邻阵元间距，用同一套距离单位。``R_{ij} = |\rho|^{|i-j|\,\mathrm{spacing}} \exp(j \arg(\rho)\,(j-i)\,\mathrm{spacing})``，对角为 1。``\rho = 0`` 在任何 `spacing` 下都是单位阵。`spacing = 0` 且 ``|\rho| > 0`` 是全 1 矩阵。不是 ``J_0(2\pi d/\lambda)``。
- 接入：`apply_kronecker(x, R_tx, R_rx, noise_var)` 先抽 ``H``，再调用 `apply_mimo`，返回 `(y, H)`。ML / ZF / MMSE / SIC 仍直接吃这个 ``H``。`mmse_detect` 仍假设符号不相关、噪声在接收端空间白，不按 ``R_t`` 或 ``R_r`` 白化。
- 对照：`test/kronecker_mimo.jl`。单位阵和 ``\rho = 0`` 与 `mimo_channel` 逐比特相同。``N = 20000`` 的样本矩对照上面的乘积，绝对误差不超过 0.06。``\rho = 1`` 时发射相关使两列相等、接收相关使两行相等。无噪声时 ML / ZF / MMSE / SIC 恢复符号。``\texttt{noise\_var} = 0.05``、各 8000 个比特，``\rho = 0.95`` 的 ZF 误码率 Clopper–Pearson 区间（``\alpha = 10^{-4}``）整个高于 ``\rho = 0``。没有新的误码闭式。演示 `demos/kronecker_mimo.jl` 不加载 Plots。

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
