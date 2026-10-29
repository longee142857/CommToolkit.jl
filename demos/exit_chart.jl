# EXIT chart utilities: mutual-information estimators checked against the
# BIAWGN capacity integral, and demapper / decoder transfer curves printed as
# tables. Teaching measurement, not NR LDPC EXIT production tooling.
# The transfer curves use the module's seeded stream; this demo prints with
# the default seed 0x0c0de1162 and does not load Plots.
# Run from the repo root:
#   julia --project=. demos/exit_chart.jl

using Printf
using CommToolkit

println("== BIAWGN consistent LLRs: estimators vs capacity integral ==")
println("   (exit_demap_transfer on BPSK at I_A = 0 measures the same quantity)")
@printf("%6s  %9s  %9s  %9s  %9s\n", "Es/N0", "capacity", "demap soft", "demap hist", "J(2/σ)")
for db in (-2.0, 0.0, 2.0, 4.0)
    σ2 = 1 / (2 * 10^(db / 10))
    # capacity: quadrature of 1 - E[log2(1 + exp(-L))], L = 2y/σ², y ~ N(1, σ²)
    σ = sqrt(σ2)
    f(y) = begin
        x = -2y / σ2
        sp = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))
        (1 - sp / log(2)) * exp(-(y - 1)^2 / (2σ2)) / (σ * sqrt(2π))
    end
    a, b = 1 - 12σ, 1 + 12σ
    n = 8192
    h = (b - a) / n
    s = f(a) + f(b)
    for i in 1:(n - 1)
        s += (isodd(i) ? 4 : 2) * f(a + i * h)
    end
    cap = s * h / 3
    soft = exit_demap_transfer(BPSK(); ia=[0.0], esn0_db=db, frames=25, n_symbols=8000)
    hist = exit_demap_transfer(BPSK(); ia=[0.0], esn0_db=db, frames=25, n_symbols=8000,
                               method=:hist)
    @printf("%5.1f   %9.5f  %9.5f  %9.5f  %9.5f\n",
            db, cap, soft.ie[1], hist.ie[1], exit_j(2 / σ))
end

println()
println("== Demapper transfer I_E = T(I_A), Es/N0 = 6 dB ==")
ia = 0.0:0.2:1.0
gray = exit_demap_transfer(PSK(8); ia=ia, esn0_db=6.0, frames=4, n_symbols=4000)
sp = exit_demap_transfer(PSK(8); ia=ia, esn0_db=6.0, frames=4, n_symbols=4000,
                         labels=bicm_psk8_id_labels())
@printf("%6s  %12s  %12s\n", "I_A", "Gray 8-PSK", "set-part 8-PSK")
for i in eachindex(ia)
    @printf("%5.2f   %11.4f   %11.4f\n", ia[i], gray.ie[i], sp.ie[i])
end

println()
println("== Decoder transfer of the terminated (7, 5) log-BCJR, 384 info bits ==")
iabc = [0.0, 0.3, 0.6, 0.9, 0.99]
@printf("%8s  %8s  %8s  %8s\n", "I_A", "0 dB", "1 dB", "3 dB")
rows = [exit_decoder_transfer(; decoder=:bcjr, ia=iabc, ebn0_db=db, frames=8).ie
        for db in (0.0, 1.0, 3.0)]
for i in eachindex(iabc)
    @printf("%7.2f   %7.4f   %7.4f   %7.4f\n", iabc[i], rows[1][i], rows[2][i], rows[3][i])
end

println()
println("== Consistency caveat: flooding sum-product on a cyclic LDPC graph ==")
println("   overconfident wrong LLRs punish the soft estimator, not the histogram")
code = array_ldpc(3, 6, 11)
for db in (1.5, 2.5, 3.5)
    s = exit_decoder_transfer(; decoder=:ldpc, code=code, ia=[0.0], ebn0_db=db,
                              frames=16, algorithm=:sum_product)
    h = exit_decoder_transfer(; decoder=:ldpc, code=code, ia=[0.0], ebn0_db=db,
                              frames=16, algorithm=:sum_product, method=:hist, nbins=101)
    @printf("%5.1f dB: soft = %9.4f   hist = %7.4f\n", db, s.ie[1], h.ie[1])
end
println()
println("Plot the two curves with plot_exit_chart(demap, decoder) after `using Plots`.")
