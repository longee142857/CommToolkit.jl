# BICM-ID. The Monte Carlo stream is xorshift + Box–Muller, not Julia's global
# RNG: 1.10 and later 1.x do not share `randn`, and both CI jobs include this file.

function _bicm_xorshift(state::UInt64)
    z = state
    z = xor(z, z >> 12)
    z = xor(z, z << 25)
    z = xor(z, z >> 27)
    return z
end

function _bicm_bits(k::Int, seed::Integer)
    st = UInt64(seed)
    bits = Vector{Int}(undef, k)
    for i in 1:k
        st = _bicm_xorshift(st)
        bits[i] = Int((st * 0x2545f4914f6cdd1d) >> 63)
    end
    return bits
end

function _bicm_normals(n::Int, seed::Integer)
    st = Ref(UInt64(seed))
    function uni()
        st[] = _bicm_xorshift(st[])
        u = st[] * 0x2545f4914f6cdd1d
        return (u >> 11) * 2.0^-53
    end
    out = Vector{Float64}(undef, n)
    i = 1
    while i <= n
        u1 = max(uni(), 1e-12)
        u2 = uni()
        r = sqrt(-2 * log(u1))
        θ = 2π * u2
        out[i] = r * cos(θ)
        i += 1
        if i <= n
            out[i] = r * sin(θ)
            i += 1
        end
    end
    return out
end

function _bicm_logadd(a::Float64, b::Float64)
    a == -Inf && return b
    b == -Inf && return a
    return a >= b ? a + log1p(exp(b - a)) : b + log1p(exp(a - b))
end

@testset "BICM labels and Gray demap" begin
    @test bicm_labels(BPSK()) == [0, 1]
    @test bicm_labels(QPSK()) == 0:3
    @test bicm_labels(PSK(8)) == 0:7
    @test bicm_labels(QAM(16)) == 0:15
    @test_throws ArgumentError bicm_labels(DBPSK())
    @test_throws ArgumentError bicm_labels(DQPSK())

    # Angle-order set partition: MSB = m mod 2, each subset a QPSK.
    angle = [bicm_psk8_id_labels()[xor(m, m >> 1) + 1] for m in 0:7]
    @test angle == [0, 4, 2, 6, 1, 5, 3, 7]
    @test all(m -> ((angle[m + 1] >> 2) & 1) == (m % 2), 0:7)

    y16 = ComplexF64[0.2 + 0.5im, -0.7 - 0.1im, 0.4 + 0.0im]
    for scheme in (BPSK(), QPSK(), PSK(8), QAM(4), QAM(16), QAM(64))
        k = bits_per_symbol(scheme)
        y = scheme isa BPSK ? real.(y16) : y16
        for approx in (:exact, :maxlog)
            @test bicm_demap(y, scheme, 0.25; approx=approx) ==
                  demodulate_soft(y, scheme, 0.25; approx=approx)
            @test bicm_demap(y, scheme, 0.25; approx=approx, labels=bicm_labels(scheme)) ==
                  demodulate_soft(y, scheme, 0.25; approx=approx)
        end
        La0 = zeros(length(y) * k)
        @test bicm_demap(y, scheme, 0.3; apriori=La0) ≈
              demodulate_soft(y, scheme, 0.3) atol = 1e-8
    end
    @test bicm_demap(ComplexF64[], QPSK(), 0.4) == Float64[]
    @test bicm_demap(Float64[], BPSK(), 0.4; approx=:maxlog) == Float64[]

    # Constant gain is zero-forcing with noise variance σ²/|h|².
    h = 0.7 + 0.2im
    σ2 = 0.3
    y = [0.15 - 0.4im, -0.55 + 0.25im, 0.05 + 0.8im]
    for scheme in (QPSK(), QAM(16), PSK(8))
        for approx in (:exact, :maxlog)
            got = bicm_demap(h .* y, scheme, σ2; gains=fill(h, length(y)), approx=approx)
            ref = demodulate_soft(y, scheme, σ2 / abs2(h); approx=approx)
            @test got ≈ ref atol = 1e-8
        end
    end

    # BPSK coherent LLR: L = 2 Re(h* y) / σ². Imaginary noise cancels when h = 1.
    yb = ComplexF64[0.3 - 0.8im]
    hb = 0.4 - 0.5im
    σb = 0.2
    @test bicm_demap(yb, BPSK(), σb; gains=[hb]) ≈ [2 * real(conj(hb) * yb[1]) / σb]
    @test bicm_demap(real.(yb), BPSK(), σb) ≈ demodulate_soft(real.(yb), BPSK(), σb)

    # Gray 16-QAM: an in-phase prior does not move the quadrature LLRs, and a
    # bit's own prior is not part of its extrinsic LLR.
    yq = ComplexF64[0.2 + 0.5im]
    L0 = bicm_demap(yq, QAM(16), 0.3)
    La = [2.5, 0.0, 0.0, 0.0]
    L1 = bicm_demap(yq, QAM(16), 0.3; apriori=La)
    @test L1[1] ≈ L0[1] atol = 1e-8
    @test !(L1[2] ≈ L0[2])
    @test L1[3:4] ≈ L0[3:4] atol = 1e-8
    Lq = bicm_demap(yq, QAM(16), 0.3; apriori=[0.0, 0.0, -1.5, 0.0], approx=:maxlog)
    Lq0 = bicm_demap(yq, QAM(16), 0.3; approx=:maxlog)
    @test Lq[1:2] ≈ Lq0[1:2] atol = 1e-8
    @test Lq[3] ≈ Lq0[3] atol = 1e-8

    # Independent enumeration of one 8-PSK symbol, both approximations.
    labs = bicm_psk8_id_labels()
    pts = constellation(PSK(8))
    y8 = 0.15 - 0.42im
    h8 = 0.8 + 0.3im
    La8 = [0.4, -1.1, 0.7]
    σ8 = 0.35
    invN0 = 1 / (2 * σ8)
    for approx in (:exact, :maxlog)
        ref = zeros(3)
        for j in 1:3
            a0 = -Inf
            a1 = -Inf
            for m in 1:8
                b = [(labs[m] >> (3 - t)) & 1 for t in 1:3]
                prior = 0.0
                for t in 1:3
                    t == j && continue
                    prior += 0.5 * (1 - 2 * b[t]) * La8[t]
                end
                met = -abs2(y8 - h8 * pts[m]) * invN0 + prior
                if approx === :maxlog
                    if b[j] == 0
                        a0 = max(a0, met)
                    else
                        a1 = max(a1, met)
                    end
                else
                    if b[j] == 0
                        a0 = _bicm_logadd(a0, met)
                    else
                        a1 = _bicm_logadd(a1, met)
                    end
                end
            end
            ref[j] = a0 - a1
        end
        @test bicm_demap([y8], PSK(8), σ8; apriori=La8, labels=labs, gains=[h8], approx=approx) ≈
              ref atol = 1e-8
    end

    @test_throws ArgumentError bicm_demap([1.0], BPSK(), 0.0)
    @test_throws ArgumentError bicm_demap([1.0], BPSK(), -0.2)
    @test_throws ArgumentError bicm_demap([1.0], BPSK(), Inf)
    @test_throws ArgumentError bicm_demap([NaN], BPSK(), 0.2)
    @test_throws ArgumentError bicm_demap([0.1], QPSK(), 0.2)
    @test_throws ArgumentError bicm_demap([0.1], QAM(16), 0.2)
    @test_throws ArgumentError bicm_demap([1.0 + 0im], DBPSK(), 0.2)
    @test_throws ArgumentError bicm_demap([1.0 + 0im], QPSK(), 0.2; approx=:approx)
    @test_throws ArgumentError bicm_demap([1.0 + 0im], QPSK(), 0.2; apriori=[0.0])
    @test_throws ArgumentError bicm_demap([1.0 + 0im], QPSK(), 0.2; apriori=[0.0, NaN])
    @test_throws ArgumentError bicm_demap([1.0 + 0im], QPSK(), 0.2; gains=[1.0, 1.0])
    @test_throws ArgumentError bicm_demap([1.0 + 0im], QPSK(), 0.2; gains=[NaN])
    @test_throws ArgumentError bicm_demap([1.0 + 0im], PSK(8), 0.2; labels=[0, 1])
    @test_throws ArgumentError bicm_demap([1.0 + 0im], PSK(8), 0.2; labels=[0, 0, 1, 2, 3, 4, 5, 6])
    @test_throws ArgumentError bicm_demap([1.0 + 0im], PSK(8), 0.2; labels=[0, 1, 2, 3, 4, 5, 6, 8])
end

@testset "BICM convolutional SISO" begin
    # Coded-bit extrinsic and information APP against an enumeration of the
    # terminated (7, 5) code. n = 4 information bits, 16 codewords.
    llr = _bicm_normals(2 * (4 + 2), 0x11)
    info, coded = CommToolkit.BICM._conv_siso(llr)
    n_bits = 4
    acc_u0 = fill(-Inf, n_bits)
    acc_u1 = fill(-Inf, n_bits)
    acc_c0 = fill(-Inf, length(llr))
    acc_c1 = fill(-Inf, length(llr))
    for mask in 0:(2^n_bits - 1)
        bits = [(mask >> (k - 1)) & 1 for k in 1:n_bits]
        word = conv_encode(bits)
        logp = 0.0
        for i in eachindex(word)
            logp += 0.5 * (1 - 2 * word[i]) * llr[i]
        end
        for k in 1:n_bits
            if bits[k] == 0
                acc_u0[k] = _bicm_logadd(acc_u0[k], logp)
            else
                acc_u1[k] = _bicm_logadd(acc_u1[k], logp)
            end
        end
        for i in eachindex(word)
            if word[i] == 0
                acc_c0[i] = _bicm_logadd(acc_c0[i], logp)
            else
                acc_c1[i] = _bicm_logadd(acc_c1[i], logp)
            end
        end
    end
    @test info ≈ acc_u0 .- acc_u1 atol = 1e-8
    @test coded ≈ (acc_c0 .- acc_c1) .- llr atol = 1e-8
    @test info ≈ bcjr_decode(llr) atol = 1e-8

    bits = _bicm_bits(14, 0x22)
    # 2*(14+2) = 32 coded bits: 16-QAM, 8 symbols, 8 rows.
    syms = bicm_map(bits, QAM(16); n_rows=8)
    z = _bicm_normals(2 * length(syms), 0x23)
    y = syms .+ 0.45 .* (z[1:2:end] .+ im .* z[2:2:end])
    σ2 = 0.45^2
    r = bicm_decode(y, QAM(16), σ2, 14; iterations=1, n_rows=8)
    L = block_deinterleave(demodulate_soft(y, QAM(16), σ2), 8)
    @test r.llr ≈ bcjr_decode(L) atol = 1e-8
    @test r.bits == Int.(bcjr_decode(L) .< 0)
    @test r.iterations == 1
    @test sprint(show, r) == "BICMResult(k=14, iterations=1)"

    rs = bicm_decode(y, QAM(16), σ2, 14; decoder=:sova, iterations=1, n_rows=8)
    @test rs.llr ≈ sova_decode(L) atol = 1e-8
    @test rs.bits == Int.(sova_decode(L) .< 0)

    # BPSK and Gray QPSK: extrinsic demapper LLRs ignore the bit prior, so a
    # second pass repeats the first.
    for scheme in (BPSK(), QPSK())
        nb = 20
        info_bits = _bicm_bits(nb, scheme isa BPSK ? 0x31 : 0x32)
        tx = bicm_map(info_bits, scheme)
        nz = scheme isa BPSK ? _bicm_normals(length(tx), 0x33) : _bicm_normals(2 * length(tx), 0x34)
        yy = scheme isa BPSK ? tx .+ 0.55 .* nz :
             tx .+ 0.55 .* (nz[1:2:end] .+ im .* nz[2:2:end])
        a = bicm_decode(yy, scheme, 0.55^2, nb; iterations=1)
        b = bicm_decode(yy, scheme, 0.55^2, nb; iterations=3, llr_clip=20)
        @test a.llr ≈ b.llr atol = 1e-6
    end

    # Noiseless set-partitioning 8-PSK, interleaved, several passes.
    sp = bicm_psk8_id_labels()
    info8 = [1, 0, 1, 1, 0, 0, 1]
    tx8 = bicm_map(info8, PSK(8); n_rows=6, labels=sp)
    @test length(tx8) == 6
    clean = bicm_decode(tx8, PSK(8), 1e-6, length(info8); iterations=3, n_rows=6, labels=sp)
    @test clean.bits == info8
    @test all(clean.llr .* (1 .- 2 .* info8) .> 0)

    # Same block in noise: feedback moves the information LLRs, and the clip binds.
    z8 = _bicm_normals(2 * length(tx8), 0x41)
    y8n = tx8 .+ 0.55 .* (z8[1:2:end] .+ im .* z8[2:2:end])
    p1 = bicm_decode(y8n, PSK(8), 0.55^2, length(info8); iterations=1, n_rows=6, labels=sp)
    p4 = bicm_decode(y8n, PSK(8), 0.55^2, length(info8); iterations=4, n_rows=6, labels=sp)
    pclip = bicm_decode(y8n, PSK(8), 0.55^2, length(info8); iterations=4, n_rows=6,
                        labels=sp, llr_clip=0.05)
    @test !(p1.llr ≈ p4.llr)
    @test !(p4.llr ≈ pclip.llr)

    # Ideal-CSI gain on a noiseless Gray QPSK block. 2*(10+2) = 24 coded bits.
    info_q = _bicm_bits(10, 0x51)
    tq = bicm_map(info_q, QPSK(); n_rows=4)
    @test length(tq) == 12
    hq = _bicm_normals(2 * length(tq), 0x52)
    h = hq[1:2:end] .+ im .* hq[2:2:end]
    h .+= 1.0
    yq = h .* tq
    rq = bicm_decode(yq, QPSK(), 1e-5, 10; iterations=1, n_rows=4, gains=h)
    @test rq.bits == info_q

    empty_tx = bicm_map(Int[], BPSK())
    re = bicm_decode(empty_tx, BPSK(), 0.3, 0; iterations=1)
    @test re.bits == Int[]
    @test re.llr == Float64[]

    @test_throws ArgumentError bicm_decode(tx8, PSK(8), 1e-3, 7; decoder=:sova, iterations=2, labels=sp)
    @test_throws ArgumentError bicm_decode(tx8, PSK(8), 1e-3, 7; iterations=0)
    @test_throws ArgumentError bicm_decode(tx8, PSK(8), 1e-3, 7; iterations=33)
    @test_throws ArgumentError bicm_decode(tx8, PSK(8), 1e-3, 7; llr_clip=0)
    @test_throws ArgumentError bicm_decode(tx8, PSK(8), 1e-3, 3; n_rows=6, labels=sp)
    @test_throws ArgumentError bicm_decode(tx8, PSK(8), 1e-3, 7; n_rows=-1, labels=sp)
    @test_throws ArgumentError bicm_decode(tx8, PSK(8), 1e-3, 7; decoder=:unknown)
    @test_throws ArgumentError bicm_map([0, 1], DBPSK())
    @test_throws ArgumentError bicm_map([0, 1, 0], PSK(8))
    @test_throws ArgumentError bicm_map([0, 1], BPSK(); encoder=:unknown)
    @test_throws ArgumentError bicm_map([0, 1], BPSK(); n_rows=-2)
end

@testset "BICM LDPC and turbo docks" begin
    code = array_ldpc(2, 4, 5)
    @test code.n == 20
    info = _bicm_bits(code.k, 0x61)
    syms = bicm_map(info, QPSK(); encoder=:ldpc, code=code, n_rows=4)
    @test length(syms) == 10
    quiet = bicm_decode(syms, QPSK(), 1e-6, code.k; decoder=:ldpc, code=code,
                        n_rows=4, iterations=2, ldpc_iterations=8, algorithm=:sum_product)
    @test quiet.bits == info

    z = _bicm_normals(2 * length(syms), 0x62)
    y = syms .+ 0.35 .* (z[1:2:end] .+ im .* z[2:2:end])
    σ2 = 0.35^2
    one = bicm_decode(y, QPSK(), σ2, code.k; decoder=:ldpc, code=code, n_rows=4,
                      iterations=1, ldpc_iterations=6, algorithm=:sum_product)
    L = block_deinterleave(demodulate_soft(y, QPSK(), σ2), 4)
    detail = ldpc_decode_details(code, L; algorithm=:sum_product, max_iter=6, early_stop=false)
    @test one.llr ≈ detail.posterior[code.info_pos] atol = 1e-8
    @test one.bits == detail.bits
    two = bicm_decode(y, QPSK(), σ2, code.k; decoder=:ldpc, code=code, n_rows=4,
                      iterations=2, ldpc_iterations=6, algorithm=:normalized_min_sum)
    @test length(two.llr) == code.k
    @test all(isfinite, two.llr)
    @test two.iterations == 2

    # 8-PSK needs a multiple of 3. n = 12.
    g = gallager_ldpc(12, 2, 3)
    @test g.n % 3 == 0
    ig = _bicm_bits(g.k, 0x63)
    sg = bicm_map(ig, PSK(8); encoder=:ldpc, code=g, labels=bicm_psk8_id_labels())
    rg = bicm_decode(sg, PSK(8), 1e-6, g.k; decoder=:ldpc, code=g,
                     labels=bicm_psk8_id_labels(), iterations=2, ldpc_iterations=10)
    @test rg.bits == ig

    perm = turbo_interleaver(8)
    itb = _bicm_bits(8, 0x71)
    tb = bicm_map(itb, BPSK(); encoder=:turbo, perm=perm, n_rows=4)
    @test length(tb) == 3 * 8 + 8
    rt = bicm_decode(tb, BPSK(), 1e-5, 8; decoder=:turbo, perm=perm, n_rows=4,
                     iterations=1, turbo_iterations=3)
    @test rt.bits == itb
    zn = _bicm_normals(length(tb), 0x72)
    yt = tb .+ 0.4 .* zn
    got = bicm_decode(yt, BPSK(), 0.16, 8; decoder=:turbo, perm=perm, n_rows=4,
                      iterations=1, turbo_iterations=2)
    Lt = block_deinterleave(demodulate_soft(yt, BPSK(), 0.16), 4)
    @test got.llr ≈ turbo_decode(Lt, perm; iterations=2) atol = 1e-8

    @test_throws ArgumentError bicm_decode(yt, BPSK(), 0.16, 8; decoder=:turbo, perm=perm, iterations=2)
    @test_throws ArgumentError bicm_decode(yt, BPSK(), 0.16, 8; decoder=:turbo, iterations=1)
    @test_throws ArgumentError bicm_map(itb, BPSK(); encoder=:turbo)
    @test_throws ArgumentError bicm_map(itb, BPSK(); encoder=:ldpc)
    @test_throws ArgumentError bicm_decode(syms, QPSK(), 0.2, code.k; decoder=:ldpc)
    @test_throws ArgumentError bicm_decode(syms, QPSK(), 0.2, code.k + 1; decoder=:ldpc, code=code)
end

@testset "BICM-ID gain" begin
    # Terminated (7, 5), 8-PSK set partitioning, 12×67 block interleaver.
    # Es = 1, R = 400/804, Eb = 1/(R * 3) = 0.67.
    # σ² = Eb / (2 Eb/N0). Max-log demapper. Seed 0x0c0de1162, per frame.
    # Not Gray, not NR, not a capacity claim. Clopper–Pearson α = 1e-4 treats
    # bit and frame errors as i.i.d., the same counting interval as the other tests.
    labels = bicm_psk8_id_labels()
    n_info = 400
    n_rows = 12
    n_coded = 2 * (n_info + 2)
    @test n_coded == 804
    Eb = 1 / ((n_info / n_coded) * 3)
    @test Eb ≈ 0.67 atol = 1e-12

    function count_passes(snr_db::Float64, n_frames::Int; fading::Bool)
        σ = sqrt(Eb / (2 * 10^(snr_db / 10)))
        σ2 = σ^2
        bit_err = zeros(Int, 2)
        frame_err = zeros(Int, 2)
        for f in 1:n_frames
            bits = _bicm_bits(n_info, 0x0c0de1162 + UInt64(f))
            syms = bicm_map(bits, PSK(8); n_rows=n_rows, labels=labels)
            nsym = length(syms)
            if fading
                z = _bicm_normals(4 * nsym, 0x0c0de1162 + 0x10000 + UInt64(f))
                h = (z[1:2:(2nsym)] .+ im .* z[2:2:(2nsym)]) ./ sqrt(2)
                noise = σ .* (z[(2nsym + 1):2:end] .+ im .* z[(2nsym + 2):2:end])
                y = h .* syms .+ noise
                gains = h
            else
                z = _bicm_normals(2 * nsym, 0x0c0de1162 + 0x10000 + UInt64(f))
                noise = σ .* (z[1:2:end] .+ im .* z[2:2:end])
                y = syms .+ noise
                gains = nothing
            end
            for (j, it) in enumerate((1, 4))
                dec = bicm_decode(y, PSK(8), σ2, n_info; iterations=it, n_rows=n_rows,
                                  labels=labels, approx=:maxlog, gains=gains)
                e = count(!=(0), dec.bits .- bits)
                bit_err[j] += e
                e > 0 && (frame_err[j] += 1)
            end
        end
        return bit_err, frame_err
    end

    bits_awgn, frames_awgn = count_passes(6.0, 80; fading=false)
    nbits = 80 * n_info
    lo1, _ = clopper_pearson(bits_awgn[1], nbits, BER_ALPHA)
    _, hi4 = clopper_pearson(bits_awgn[2], nbits, BER_ALPHA)
    flo1, _ = clopper_pearson(frames_awgn[1], 80, BER_ALPHA)
    _, fhi4 = clopper_pearson(frames_awgn[2], 80, BER_ALPHA)
    @test bits_awgn[1] >= 500
    @test bits_awgn[2] >= 20
    @test frames_awgn[2] >= 5
    @test hi4 < lo1
    @test fhi4 < flo1

    # Ideal-CSI symbol Rayleigh, same Eb/N0 definition, h ~ CN(0, 1) known to the demapper.
    bits_f, frames_f = count_passes(8.0, 40; fading=true)
    nbf = 40 * n_info
    lof, _ = clopper_pearson(bits_f[1], nbf, BER_ALPHA)
    _, hif = clopper_pearson(bits_f[2], nbf, BER_ALPHA)
    flof, _ = clopper_pearson(frames_f[1], 40, BER_ALPHA)
    _, fhif = clopper_pearson(frames_f[2], 40, BER_ALPHA)
    @test bits_f[1] >= 100
    @test bits_f[2] >= 20
    @test hif < lof
    @test fhif < flof
end
