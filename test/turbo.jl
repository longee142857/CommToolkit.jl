# Deterministic bits and N(0, 1) samples. Julia's global RNG is not stable
# across 1.10 and later 1.x, and this file is included by both CI jobs.
function _xorshift(state::UInt64)
    z = state
    z = xor(z, z >> 12)
    z = xor(z, z << 25)
    z = xor(z, z >> 27)
    return z
end

function _det_bits(k::Int, seed::UInt64)
    st = seed
    bits = Vector{Int}(undef, k)
    for i in 1:k
        st = _xorshift(st)
        bits[i] = Int((st * 0x2545f4914f6cdd1d) >> 63)
    end
    return bits
end

function _std_normals(n::Int, seed::UInt64)
    st = Ref(seed)
    function uni()
        st[] = _xorshift(st[])
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

function _rsc_brute(info_llr::Vector{Float64}, par::Vector{Float64}, tail::Vector{Float64}, La::Vector{Float64})
    k = length(La)
    m = 2
    acc0 = fill(-Inf, k)
    acc1 = fill(-Inf, k)
    for mask in 0:(2^k - 1)
        info = [(mask >> (j - 1)) & 1 for j in 1:k]
        parity, tbits = CommToolkit.Coding._rsc_encode(info)
        logp = 0.0
        for i in 1:k
            logp += 0.5 * (1 - 2 * info[i]) * info_llr[i]
            logp += 0.5 * (1 - 2 * info[i]) * La[i]
        end
        for i in 1:(k + m)
            logp += 0.5 * (1 - 2 * parity[i]) * par[i]
        end
        for j in 1:m
            logp += 0.5 * (1 - 2 * tbits[j]) * tail[j]
        end
        for i in 1:k
            if info[i] == 0
                acc0[i] = acc0[i] == -Inf ? logp : max(acc0[i], logp) + log1p(exp(-abs(acc0[i] - logp)))
            else
                acc1[i] = acc1[i] == -Inf ? logp : max(acc1[i], logp) + log1p(exp(-abs(acc1[i] - logp)))
            end
        end
    end
    app = acc0 .- acc1
    return app .- La .- info_llr, app
end

@testset "Turbo interleaver and encoder" begin
    @test turbo_interleaver(2) == [2, 1]
    for k in (4, 6, 8, 9, 15, 16, 25, 32)
        perm = turbo_interleaver(k)
        @test sort(perm) == collect(1:k)
        @test perm != collect(1:k)
    end
    # Square-ish factor: 15 = 3 × 5, rows = 3.
    @test turbo_interleaver(15) == block_interleave(collect(1:15), 3, 5)
    @test turbo_interleaver(8) == block_interleave(collect(1:8), 2, 4)
    # Prime length: stride, not the rectangular map (7 is prime).
    p7 = turbo_interleaver(7)
    @test sort(p7) == collect(1:7)
    @test p7 == [mod((i - 1) * 2, 7) + 1 for i in 1:7]
    @test p7 != collect(1:7)

    # Hand-checked (1, 5/7) termination. perm swaps the two information bits.
    bits = [1, 0]
    perm = [2, 1]
    cw = turbo_encode(bits, perm)
    @test cw == [1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1]
    @test length(cw) == 3 * 2 + 8
    @test turbo_encode([0, 0], [2, 1]) == zeros(Int, 14)

    # State returns to 0 for every information word of length 6, and the
    # second encoder really sees the interleaved bits.
    perm6 = turbo_interleaver(6)
    for mask in 0:63
        info = [(mask >> j) & 1 for j in 0:5]
        p1, t1 = CommToolkit.Coding._rsc_encode(info)
        p2, t2 = CommToolkit.Coding._rsc_encode(info[perm6])
        @test turbo_encode(info, perm6) == vcat(info, p1, t1, p2, t2)
    end

    @test_throws ArgumentError turbo_interleaver(1)
    @test_throws ArgumentError turbo_interleaver(1025)
    @test_throws ArgumentError turbo_encode([0, 1, 2], [1, 2, 3])
    @test_throws ArgumentError turbo_encode([0], [1])
    @test_throws ArgumentError turbo_encode([0, 1], [1, 1])
    @test_throws ArgumentError turbo_encode([0, 1], [0, 1])
    @test_throws ArgumentError turbo_encode([0, 1], [1, 2, 3])
end

@testset "Turbo component BCJR" begin
    # k = 5 is 32 codewords. Log-BCJR extrinsic and APP match the enumeration.
    k = 5
    sys = [0.2, -0.5, 1.1, 0.0, 0.4, 0.3, -0.2]
    par = [0.6, -0.1, 0.8, 0.2, -0.4, 0.5, 0.1]
    La = [0.3, -0.7, 0.0, 0.2, -0.1]
    Le, app = CommToolkit.Coding._rsc_bcjr(sys, par, La)
    e_ref, a_ref = _rsc_brute(sys[1:k], par, sys[(k + 1):end], La)
    @test Le ≈ e_ref atol = 1e-6 rtol = 1e-6
    @test app ≈ a_ref atol = 1e-6 rtol = 1e-6

    La0 = zeros(k)
    Le0, app0 = CommToolkit.Coding._rsc_bcjr(sys, par, La0)
    e0, a0 = _rsc_brute(sys[1:k], par, sys[(k + 1):end], La0)
    @test Le0 ≈ e0 atol = 1e-6
    @test app0 ≈ a0 atol = 1e-6
    # No channel, no prior: the terminated code is balanced at every information bit.
    @test CommToolkit.Coding._rsc_bcjr(zeros(k + 2), zeros(k + 2), zeros(k))[2] ≈ zeros(k) atol = 1e-8
end

@testset "Turbo decode" begin
    perm = [2, 1]
    bits = [1, 0]
    cw = turbo_encode(bits, perm)
    llr = 30.0 .* (1.0 .- 2.0 .* cw)
    for it in (1, 2, 4)
        app = turbo_decode(llr, perm; iterations=it)
        @test length(app) == 2
        @test all(isfinite, app)
        @test Int.(app .< 0) == bits
        @test all(abs.(app) .> 5)
    end

    # One flipped channel sign on an otherwise decisive word.
    llr1 = copy(llr)
    llr1[4] *= -1
    @test Int.(turbo_decode(llr1, perm; iterations=2) .< 0) == bits

    for k in (2, 7, 8, 15)
        π = turbo_interleaver(k)
        info = _det_bits(k, 0x0c0de1162)
        word = turbo_encode(info, π)
        hard = turbo_decode(25.0 .* (1.0 .- 2.0 .* word), π; iterations=2)
        @test Int.(hard .< 0) == info
    end

    # Same deterministic noise record, two iteration counts. Rate is k/(3k+8),
    # Es = 1, so Eb = (3k+8)/k and σ² = Eb / (2 Eb/N0). On Julia 1.10.10 this
    # point (k = 64, 24 frames, 0.7 dB) was 79 errors at 1 iteration and 33 at
    # 6. The thresholds sit below that gap. Not a 3GPP BLER table, and not a
    # claim that every extra iteration helps on every frame.
    k = 64
    perm = turbo_interleaver(k)
    n = 3k + 8
    eb = n / k
    snr_db = 0.7
    γ = 10^(snr_db / 10)
    σ2 = eb / (2 * γ)
    nframes = 24
    e1 = 0
    e6 = 0
    for f in 1:nframes
        info = _det_bits(k, 0x1000 + UInt64(f))
        word = turbo_encode(info, perm)
        y = (1.0 .- 2.0 .* word) .+ sqrt(σ2) .* _std_normals(n, 0x9000 + UInt64(f))
        ch = (2.0 / σ2) .* y
        d1 = Int.(turbo_decode(ch, perm; iterations=1) .< 0)
        d6 = Int.(turbo_decode(ch, perm; iterations=6) .< 0)
        e1 += count(!=(0), d1 .- info)
        e6 += count(!=(0), d6 .- info)
    end
    @test e1 >= 40
    @test e1 >= e6 + 20

    @test_throws ArgumentError turbo_decode(llr, perm; iterations=0)
    @test_throws ArgumentError turbo_decode(llr, perm; iterations=33)
    @test_throws ArgumentError turbo_decode(llr[1:end-1], perm; iterations=1)
    @test_throws ArgumentError turbo_decode(fill(NaN, length(llr)), perm; iterations=1)
    @test_throws ArgumentError turbo_decode(llr, [1, 1]; iterations=1)
end
