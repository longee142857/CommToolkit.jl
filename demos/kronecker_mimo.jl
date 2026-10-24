# 平坦 Kronecker MIMO，无 Plots。
# 打印 Loyka 指数相关、样本矩，以及同一噪声下 ρ = 0 与 ρ = 0.9 的 ZF 误码。
# 模型是 H = R_r^{1/2} H_w R_t^{1/2}，不是 Weichselberger，也不是 TR 38.901 CDL/UMi。
# 运行：julia --project=. demos/kronecker_mimo.jl

using Random
using CommToolkit

Random.seed!(0x0c0de1162)

ρ = 0.6
R = exponential_correlation(3, ρ)
println("Loyka R, n = 3, ρ = $ρ, spacing = 1")
show(stdout, "text/plain", R)
println()

function moment_error(Rt, Rr, N)
    acc = zeros(ComplexF64, 2, 2, 2, 2)
    for _ in 1:N
        H = kronecker_channel(Rt, Rr)
        for s in 1:2, r in 1:2, q in 1:2, p in 1:2
            acc[p, q, r, s] += H[p, q] * conj(H[r, s])
        end
    end
    acc ./= N
    err = 0.0
    for s in 1:2, r in 1:2, q in 1:2, p in 1:2
        err = max(err, abs(acc[p, q, r, s] - Rr[p, r] * Rt[s, q]))
    end
    return err
end

Rt = exponential_correlation(2, 0.5)
Rr = exponential_correlation(2, 0.4 * cis(0.7))
N = 4000
err = moment_error(Rt, Rr, N)
println("样本矩最大绝对误差（N = $(N)，理论值 (R_r)[p,r] * (R_t)[s,q]）: $err")

function zf_ber(ρ, noise_var, n_vec)
    Rc = exponential_correlation(2, ρ)
    k = 0
    for _ in 1:n_vec
        bits = rand(0:1, 2)
        s = ComplexF64.(modulate(bits, BPSK()))
        y, H = apply_kronecker(s, Rc, Rc, noise_var)
        z = zf_detect(y, H)
        k += count(!=(0), bits .- demodulate(real.(z), BPSK()))
    end
    return k / (2 * n_vec)
end

noise_var = 0.05
n_vec = 2000
println("2×2 BPSK ZF，noise_var = $(noise_var)，$(2 * n_vec) 比特")
println("  ρ = 0   BER = ", zf_ber(0.0, noise_var, n_vec))
println("  ρ = 0.9 BER = ", zf_ber(0.9, noise_var, n_vec))
println("i.i.d. 路径仍是 mimo_channel；检测器不按 R_t、R_r 白化。")
