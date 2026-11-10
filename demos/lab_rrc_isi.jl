# 不画图的根升余弦实验。没有同步环、均衡器、AWGN 误码，也不新增导出。
#
#   julia --project=. demos/lab_rrc_isi.jl
#
# pulse_shape 会把卷积截成 length(symbols)*sps，尾部丢掉。那是眼图辅助的
# 接口，不是这条链路。这里用 apply_fir 保留全长，再 matched_filter，按已知
# 下标抽样。升余弦只打印参照，不要求有限根升余弦级联逐点等于同跨度升余弦。
#
# Plot-free root-raised-cosine lab. No synchronizer, equalizer, AWGN BER,
# or new export. Full convolution keeps the finite-FIR tail; pulse_shape
# remains the eye-diagram helper and is not used here. Raised cosine is a
# reference printout, not a pointwise target for the finite RRC cascade.
#
# 有限支撑留下的符号间隔残余是截断，按抽头测量。γ_n 只包一次求和的舍入，
# 尺度是 Σ|a_i b_i|（Higham，Accuracy and Stability of Numerical Algorithms，
# 第 2 版，(3.5)；u = eps(Float64)/2，γ_n = n u / (1 - n u)）。它不含
# 三角函数、除法、奇点分支、归一化，也不含截断。
#
# 眼图：实 BPSK 标签 ±1。g0 是级联符号间隔抽头的主项，S 是其余 |g_k| 之和。
# 最坏类间分离下界是 2(g0-S)，允许为负，不截成 0。g0 没有收成 1，所以不写
# 2(1-S)。测量是这一条序列、这一个相位上去掉非符号样点之后的
# min(正类) - max(负类)，不是连续眼宽，也不是对所有序列的保证。
# 错开一个样本只对照该相位的预测，不是同步。

using Random
using Printf
using CommToolkit

const LAB_SEED = 0x0c0de1162
const LAB_NSYM = 32

function gamma_n(n::Int)
    n == 0 && return 0.0
    u = eps(Float64) / 2
    nu = n * u
    nu < 1 || error("gamma_n undefined for n = $n")
    return nu / (1 - nu)
end

"Upper bound on γ_n · Σ|a_i b_i| from the Float64 sum of the absolute products."
function round_allowance(abs_products::Vector{Float64})
    n = length(abs_products)
    n == 0 && return 0.0
    γ = gamma_n(n)
    return γ * sum(abs_products) / (1 - γ)
end

function overlap_abs(x::Vector{Float64}, h::Vector{Float64}, k::Int)
    prods = Float64[]
    for i in eachindex(x)
        j = k - i + 1
        if 1 <= j <= length(h)
            push!(prods, abs(x[i] * h[j]))
        end
    end
    return prods
end

function sample_round_bounds(x::Vector{Float64}, h::Vector{Float64})
    return [round_allowance(overlap_abs(x, h, k)) for k in 1:(length(x) + length(h) - 1)]
end

"Sum of incoming absolute bounds. Inflate the positive sum so the Float64 total stays an upper bound."
function propagate_bounds(bounds::Vector{Float64}, h::Vector{Float64}, p::Int)
    acc = 0.0
    n = 0
    for i in eachindex(bounds)
        j = p - i + 1
        if 1 <= j <= length(h)
            acc += bounds[i] * abs(h[j])
            n += 1
        end
    end
    n == 0 && return 0.0
    return acc / (1 - gamma_n(n))
end

function upsample_symbols(symbols::Vector{Float64}, sps::Int)
    x = zeros(length(symbols) * sps)
    x[1:sps:end] .= symbols
    return x
end

function symbol_model(symbols::Vector{Float64}, r::Vector{Float64}, bound_r::Vector{Float64}, L::Int, sps::Int, k::Int)
    acc = 0.0
    prods = Float64[]
    prop = 0.0
    nprop = 0
    for m in eachindex(symbols)
        q = L + (k - m) * sps
        if 1 <= q <= length(r)
            acc += symbols[m] * r[q]
            push!(prods, abs(symbols[m] * r[q]))
            prop += bound_r[q]
            nprop += 1
        end
    end
    prop_allow = nprop == 0 ? 0.0 : prop / (1 - gamma_n(nprop))
    return acc, round_allowance(prods) + prop_allow
end

function sidelobes(r::Vector{Float64}, bound_r::Vector{Float64}, L::Int, sps::Int)
    acc = 0.0
    energy = 0.0
    prop = 0.0
    n = 0
    m = 1
    while true
        hit = false
        for q in (L + m * sps, L - m * sps)
            if 1 <= q <= length(r)
                acc += abs(r[q])
                energy += r[q]^2
                prop += bound_r[q]
                n += 1
                hit = true
            end
        end
        hit || break
        m += 1
    end
    if n == 0
        return 0.0, 0.0, 0.0
    end
    γ = gamma_n(n)
    allow = prop / (1 - γ) + γ * acc / (1 - γ)
    return acc, energy, allow
end

function rc_symbol_abs(beta::Float64, span::Int, sps::Int)
    h = raised_cosine(beta, span, sps)
    L = length(h)
    if iseven(L)
        return "L=$(L) 偶数，网格上没有 τ=0，不把空的主项写成 0"
    end
    c = (L + 1) ÷ 2
    acc = 0.0
    m = 1
    while true
        hit = false
        for i in (c + m * sps, c - m * sps)
            if 1 <= i <= L
                acc += abs(h[i])
                hit = true
            end
        end
        hit || break
        m += 1
    end
    return @sprintf("%.6e", acc)
end

function impulse_buffer(p0::Int)
    buf = zeros(p0)
    buf[p0] = 1.0
    return buf
end

# An impulse at index p0 has its autocorrelation zero-lag at p0+L-1.
# downsample's offset is that index minus one: samples skipped, not L-1 in general.
function impulse_aligned(g::Vector{Float64}, p0::Int, sps::Int)
    L = length(g)
    rev = reverse(g)
    buf = impulse_buffer(p0)
    tx = apply_fir(buf, g)
    rx = matched_filter(tx, g)
    r = matched_filter(g, g)
    idx = p0 + L - 1
    bound_tx = sample_round_bounds(buf, g)
    bound_r = round_allowance(overlap_abs(g, rev, L))
    stage = round_allowance(overlap_abs(tx, rev, idx)) + propagate_bounds(bound_tx, rev, idx)
    structural = 1 <= idx <= length(rx) &&
                 downsample_symbols(rx, sps; offset=idx - 1)[1] == rx[idx]
    return structural, idx, idx - 1, abs(rx[idx] - r[L]) <= stage + bound_r
end

function run_case(symbols::Vector{Float64}, beta::Float64, span::Int, sps::Int)
    g = rrc_filter(beta, span, sps)
    L = length(g)
    rev = reverse(g)
    x = upsample_symbols(symbols, sps)
    tx = apply_fir(x, g)
    rx = matched_filter(tx, g)
    r = matched_filter(g, g)
    bound_tx = sample_round_bounds(x, g)
    bound_r = sample_round_bounds(g, rev)
    nsym = length(symbols)
    g0 = r[L]
    S, isi_energy, allow_S = sidelobes(r, bound_r, L, sps)
    allow_g0 = bound_r[L]
    lb = 2 * (g0 - S)
    kept = downsample_symbols(rx, sps; offset=L - 1)
    ok = true
    ok &= length(x) == nsym * sps
    ok &= length(tx) == length(x) + L - 1
    ok &= length(rx) == length(x) + 2 * (L - 1)
    ok &= length(kept) == nsym + span
    ok &= g0 > 0
    max_allow = 0.0
    pos = Float64[]
    neg = Float64[]
    for k in 1:nsym
        p = L + (k - 1) * sps
        model, model_allow = symbol_model(symbols, r, bound_r, L, sps, k)
        stage2 = round_allowance(overlap_abs(tx, rev, p))
        prop = propagate_bounds(bound_tx, rev, p)
        stage_allow = stage2 + prop
        total = stage_allow + model_allow
        max_allow = max(max_allow, total)
        sample = rx[p]
        ok &= kept[k] == sample
        ok &= abs(sample - model) <= total
        ok &= symbols[k] * sample + total + allow_g0 + allow_S >= g0 - S
        if symbols[k] > 0
            push!(pos, sample)
        else
            push!(neg, sample)
        end
    end
    sep = minimum(pos) - maximum(neg)
    ok &= sep + 2 * max_allow + 2 * allow_g0 + 2 * allow_S >= lb

    # One sample later than the zero-lag index. This is a known mismatch, not a loop.
    δ = 1
    phase_pos = Float64[]
    phase_neg = Float64[]
    phase_ok = true
    for k in 1:nsym
        p = L + δ + (k - 1) * sps
        p <= length(rx) || (phase_ok = false; break)
        acc = 0.0
        prods = Float64[]
        prop_r = 0.0
        nprop = 0
        for m in 1:nsym
            q = L + δ + (k - m) * sps
            if 1 <= q <= length(r)
                acc += symbols[m] * r[q]
                push!(prods, abs(symbols[m] * r[q]))
                prop_r += bound_r[q]
                nprop += 1
            end
        end
        prop_allow = nprop == 0 ? 0.0 : prop_r / (1 - gamma_n(nprop))
        model_allow = round_allowance(prods) + prop_allow
        stage_allow = round_allowance(overlap_abs(tx, rev, p)) + propagate_bounds(bound_tx, rev, p)
        sample = rx[p]
        phase_ok &= abs(sample - acc) <= stage_allow + model_allow
        if symbols[k] > 0
            push!(phase_pos, sample)
        else
            push!(phase_neg, sample)
        end
    end
    phase_sep = phase_ok ? minimum(phase_pos) - maximum(phase_neg) : NaN
    ok &= phase_ok

    # Index 1: the skip count equals the two-filter delay L-1.
    # Index 4: the same definition moves, so the offset argument is not L-1.
    s1, idx1, off1, v1 = impulse_aligned(g, 1, sps)
    s4, idx4, off4, v4 = impulse_aligned(g, 4, sps)
    ok &= s1 && v1 && off1 == L - 1 && idx1 == L
    ok &= s4 && v4 && off4 == idx4 - 1 && off4 != L - 1

    return (; ok, beta, span, sps, L, g0, S, lb, sep, isi_energy, max_allow,
            phase_sep, up=length(x), tx=length(tx), rx=length(rx),
            first=L, n_valid=nsym, n_down=length(kept),
            rc=rc_symbol_abs(beta, span, sps), half=(L - 1) / 2)
end

function main()
    rng = Xoshiro(LAB_SEED)
    symbols = Float64[rand(rng) < 0.5 ? -1.0 : 1.0 for _ in 1:LAB_NSYM]
    npos = count(>(0), symbols)
    nneg = count(<(0), symbols)
    println("RRC lab. Local Xoshiro seed 0x", string(LAB_SEED, base=16, pad=9), ". Real BPSK ±1. n = ", LAB_NSYM, ".")
    println("Positive labels ", npos, ", negative labels ", nneg, ". One sequence, one phase.")
    println("γ_n bounds only accumulated rounding of one sum, scaled by Σ|a_i b_i|. It excludes the pulse formula, normalization, and truncation.")
    println("Eye lower bound is 2(g0-S), unclipped. 2(1-S) is not used. Wrong phase is a sampling mismatch, not a synchronizer.")
    ok = npos > 0 && nneg > 0
    cases = (
        (0.35, 6, 4),
        (0.35, 10, 4),
        (0.5, 1, 1),
        (0.0, 6, 4),
        (1.0, 3, 1),
    )
    rows = NamedTuple[]
    for (beta, span, sps) in cases
        row = run_case(symbols, beta, span, sps)
        push!(rows, row)
        ok &= row.ok
        half_note = isinteger(row.half) ? "on a sample" : "half-sample; constructor still accepts it"
        @printf("beta=%.2f span=%d sps=%d L=%d single_delay=%s (%s)\n",
                row.beta, row.span, row.sps, row.L, string(row.half), half_note)
        @printf("  upsample=%d tx=%d rx=%d first_index=%d valid=%d downsample=%d\n",
                row.up, row.tx, row.rx, row.first, row.n_valid, row.n_down)
        @printf("  g0=%.16e |g0-1|=%.3e S=%.6e lb=%.6e sep=%.6e isi_energy=%.6e allow=%.3e\n",
                row.g0, abs(row.g0 - 1), row.S, row.lb, row.sep, row.isi_energy, row.max_allow)
        @printf("  wrong_phase_sep=%.6e  rc_symbol_abs=%s  %s\n",
                row.phase_sep, row.rc, row.ok ? "pass" : "FAIL")
    end
    println("span 6 and span 10 at β=0.35, sps=4 are listed side by side. The pair does not show that residual ISI is monotone in span.")
    println("Raised-cosine numbers are a reference. The finite RRC cascade is not required to match the same-span raised cosine.")
    println(ok ? "criterion: pass" : "criterion: fail")
    return ok ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
