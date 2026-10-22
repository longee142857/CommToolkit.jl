# Block interleaving, convolutional puncturing, symbol/frame error rates,
# and block-constant Rayleigh fading.
#
# Puncturing masks are the Yasuda / IEEE 802.11 keep-patterns on the (7, 5)
# stream. Eb/N0 uses the realised information-bit energy
# (kept coded symbols) / (information bits), Es = 1, σ² = N0/2.
# Clopper–Pearson comparisons use BER_ALPHA from reference.jl.
# The 2/3 rate sits close to the mother code, so only 3/4 is required to
# clear the mother-code interval; 2/3 is a point-estimate ordering.

@testset "Symbol and frame error rates" begin
    bits_tx = [0, 1, 0, 1]
    bits_rx = [0, 1, 1, 1]
    @test ser(bits_tx, bits_rx) == ber(bits_tx, bits_rx) == 0.25
    @test ser([1.0, -1.0, 1.0], [1.0, 1.0, 1.0]) == 1 / 3
    tx_q = modulate([0, 0, 0, 1, 1, 0, 1, 1], QPSK())
    rx_q = modulate([0, 0, 1, 1, 1, 0, 1, 1], QPSK())
    @test ser(tx_q, rx_q) == 0.25
    @test ser([1.0 + 0im], [1.0 + 0im]) == 0.0
    @test_throws ArgumentError ser([1, 0], [1])
    @test_throws ArgumentError ser(Int[], Int[])

    @test fer([0, 0, 1, 0], [0, 0, 1, 1], 2) == 0.5
    @test fer([0, 0, 0, 0], [0, 0, 0, 0], 4) == 0.0
    @test fer([0, 0, 1], [0, 0, 0], 3) == 1.0
    @test fer([0.0 + 0im, 1.0 + 0im], [0.0 + 0im, 0.0 + 0im], 2) == 1.0
    @test_throws ArgumentError fer([0, 1, 0], [0, 1, 0], 2)
    @test_throws ArgumentError fer([0, 1], [0, 1], 0)
    @test_throws ArgumentError fer([0], [1, 0], 1)
    @test_throws ArgumentError fer(Int[], Int[], 1)
end

@testset "Block interleaver" begin
    x = [1, 2, 3, 4, 5, 6]
    @test block_interleave(x, 2, 3) == [1, 4, 2, 5, 3, 6]
    @test block_deinterleave([1, 4, 2, 5, 3, 6], 2, 3) == x
    @test block_interleave(x, 2) == [1, 4, 2, 5, 3, 6]
    @test block_deinterleave(block_interleave(x, 2), 2) == x
    @test block_interleave([1, 2, 3, 4], 1, 4) == [1, 2, 3, 4]
    @test block_interleave([1, 2, 3, 4], 4, 1) == [1, 2, 3, 4]
    @test block_interleave(view(x, :), 2, 3) == [1, 4, 2, 5, 3, 6]

    Random.seed!(RNG_SEED)
    for (n_rows, n_cols) in ((1, 7), (7, 1), (3, 5), (8, 8), (5, 13))
        for src in (
            rand(0:5, n_rows * n_cols),
            rand(n_rows * n_cols),
            randn(ComplexF64, n_rows * n_cols),
        )
            y = block_interleave(src, n_rows, n_cols)
            @test length(y) == length(src)
            @test block_deinterleave(y, n_rows, n_cols) == src
            @test block_deinterleave(y, n_rows) == src
            @test eltype(y) == eltype(src)
        end
    end

    @test_throws ArgumentError block_interleave([1, 2, 3], 2)
    @test_throws ArgumentError block_interleave([1, 2, 3, 4], 0)
    @test_throws ArgumentError block_interleave([1, 2, 3, 4], 2, 3)
    @test_throws ArgumentError block_interleave(Int[], 1)
    @test_throws ArgumentError block_deinterleave([1, 2, 3, 4], 3)
    @test_throws ArgumentError block_deinterleave([1, 2], 1, 1)
    @test_throws ArgumentError block_deinterleave([1, 2, 3, 4], -1, 4)
end

@testset "Convolutional puncturing" begin
    @test puncture_mask(1//2) == [1, 1]
    @test puncture_mask(2//3) == [1, 1, 1, 0]
    @test puncture_mask(3//4) == [1, 1, 1, 0, 0, 1]
    for (rate, mask) in ((1//2, puncture_mask(1//2)), (2//3, puncture_mask(2//3)), (3//4, puncture_mask(3//4)))
        period = length(mask)
        weight = count(isone, mask)
        @test (1//2) * period // weight == rate
    end
    stolen = puncture_mask(2//3)
    stolen[1] = 0
    @test puncture_mask(2//3) == [1, 1, 1, 0]

    @test_throws ArgumentError puncture_mask(4//5)
    @test_throws ArgumentError puncture_mask(2 / 3)
    @test_throws ArgumentError puncture_mask(1)
    @test_throws ArgumentError puncture(Int[1, 0, 1], Int[])
    @test_throws ArgumentError puncture(Int[1, 0, 1], [1, 2])
    @test_throws ArgumentError puncture(Int[1, 0, 1], [0, 0])
    @test_throws ArgumentError depuncture(Int[1], [1, 0], -1)
    @test_throws ArgumentError depuncture(Int[1, 0], [1, 0], 2)

    @test puncture(Int[], [1, 0]) == Int[]
    @test depuncture(Int[], [1, 1, 1, 0], 0) == Int[]
    @test puncture([1, 0, 1, 1, 0, 0], [1, 0]) == [1, 1, 0]
    @test depuncture([1, 2, 3], [1, 1, 1, 0], 3) == [1, 2, 3]
    @test depuncture([1, 2, 3], [1, 1, 1, 0], 4; fill=-1) == [1, 2, 3, -1]
    restored = depuncture(Int[1, 0, 1], [1, 0], 5; fill=0.0)
    @test restored isa Vector{Float64}
    @test restored == [1.0, 0.0, 0.0, 0.0, 1.0]

    Random.seed!(RNG_SEED)
    mother = rand(0:1, 9)
    kept = puncture(mother, [1, 1, 0])
    @test puncture(depuncture(kept, [1, 1, 0], 9; fill=0), [1, 1, 0]) == kept
    @test depuncture(kept, [1, 1, 0], 9; fill=0)[3] == 0

    # Free-distance event. The 3/4 mask deletes a 1, so integer fill is not
    # the mother word; LLR 0 at that position still decodes.
    event_bits = [1, 0, 0]
    event = conv_encode(event_bits)
    @test event == [1, 1, 1, 0, 1, 1, 0, 0, 0, 0]
    punched = puncture(event, 3//4)
    @test depuncture(punched, 3//4, length(event); fill=0) != event
    erased = depuncture(Float64.(1 .- 2 .* punched), 3//4, length(event); fill=0.0)
    @test viterbi_decode(erased) == event_bits
    @test Int.(bcjr_decode(erased) .< 0) == event_bits

    Random.seed!(RNG_SEED)
    for n in 0:24
        bits = rand(0:1, n)
        coded = conv_encode(bits)
        @test puncture(coded, 1//2) == coded
        @test depuncture(coded, 1//2, length(coded)) == coded
        for rate in (2//3, 3//4)
            punched_bits = puncture(coded, rate)
            @test length(punched_bits) == length(puncture(coded, puncture_mask(rate)))
            for scale in (1.0, 25.0)
                llr = depuncture(scale .* (1 .- 2 .* Float64.(punched_bits)), rate, length(coded); fill=0.0)
                @test viterbi_decode(llr) == bits
                post = bcjr_decode(llr)
                hard = isempty(post) ? Int[] : Int.(post .< 0)
                @test hard == bits
            end
        end
    end
end

function _punctured_error_counts(rate, snr_db, ninfo, nframes)
    k_vit = 0
    k_map = 0
    for _ in 1:nframes
        bits = rand(0:1, ninfo)
        coded = conv_encode(bits)
        punched = puncture(coded, rate)
        energy = length(punched) / length(bits)
        rx = awgn(modulate(punched, BPSK()), snr_db; energy_per_bit=energy)
        σ2 = energy / (2 * snr_db_to_linear(snr_db))
        llr = depuncture(demodulate_soft(rx, BPSK(), σ2), rate, length(coded); fill=0.0)
        k_vit += count(!=(0), viterbi_decode(llr) .- bits)
        k_map += count(!=(0), Int.(bcjr_decode(llr) .< 0) .- bits)
    end
    return k_vit, k_map
end

@testset "Punctured convolutional BER direction" begin
    # 4 dB, 256 information bits, 120 frames. Seed RNG_SEED, rates in this order.
    # Nominal keep fractions make Eb ≈ 2, 1.5 and 4/3; the tail is inside `energy`.
    Random.seed!(RNG_SEED)
    ninfo = 256
    nframes = 120
    nbits = ninfo * nframes
    snr_db = 4.0
    vit = Dict{Rational{Int},Int}()
    mapa = Dict{Rational{Int},Int}()
    for rate in (1//2, 2//3, 3//4)
        vit[rate], mapa[rate] = _punctured_error_counts(rate, snr_db, ninfo, nframes)
    end
    for k in (vit, mapa)
        @test k[1//2] >= 5
        @test k[2//3] >= 20
        @test k[3//4] >= 40
        @test k[1//2] < k[2//3] < k[3//4]
        _, hi_mother = clopper_pearson(k[1//2], nbits, BER_ALPHA)
        lo_34, _ = clopper_pearson(k[3//4], nbits, BER_ALPHA)
        @test hi_mother < lo_34
    end
end

@testset "Block Rayleigh" begin
    Random.seed!(RNG_SEED)
    h = block_rayleigh(10, 4)
    @test length(h) == 10
    @test all(h[1:4] .== h[1])
    @test all(h[5:8] .== h[5])
    @test h[9] == h[10]
    @test h[1] != h[5]
    @test h[5] != h[9]
    @test length(unique(block_rayleigh(8, 1))) == 8
    @test length(unique(block_rayleigh(5, 9))) == 1
    @test isempty(block_rayleigh(0, 3))

    Random.seed!(RNG_SEED)
    power = block_rayleigh(20_000, 5)
    @test sum(abs2, power) / length(power) ≈ 1.0 atol = 0.05

    syms = [1.0, -1.0, 0.5, -0.25]
    rx, hf = block_rayleigh_channel(syms, Inf, 1, 2)
    @test rx ./ hf ≈ ComplexF64.(syms) atol = 1e-12
    @test hf[1] == hf[2] && hf[3] == hf[4]
    rx_c, hc = block_rayleigh_channel(ComplexF64.(syms), Inf, 1, 3)
    @test rx_c ./ hc ≈ ComplexF64.(syms) atol = 1e-12
    empty_rx, empty_h = block_rayleigh_channel(ComplexF64[], Inf, 1, 4)
    @test isempty(empty_rx) && isempty(empty_h)

    Random.seed!(RNG_SEED)
    probe = ones(ComplexF64, 6)
    rx1, h1 = block_rayleigh_channel(probe, 0.0, 1, 3)
    Random.seed!(RNG_SEED)
    h2 = block_rayleigh(6, 3)
    rx2 = awgn(h2 .* probe, 0.0, 1)
    @test h1 == h2
    @test rx1 == rx2
    @test all(h1[1:3] .== h1[1])
    @test all(h1[4:6] .== h1[4])

    # Coherent BPSK LLR matches demodulate_soft on the zero-forced sample.
    h0 = 0.3 - 0.4im
    y0 = h0 * (-1.0 + 0im) + (0.1 + 0.2im)
    σ2 = 0.25
    manual = 2 * real(conj(h0) * y0) / σ2
    @test manual ≈ only(demodulate_soft([y0 / h0], BPSK(), σ2 / abs2(h0)))

    @test_throws ArgumentError block_rayleigh(-1, 2)
    @test_throws ArgumentError block_rayleigh(4, 0)
    @test_throws ArgumentError block_rayleigh_channel([1.0, -1.0], 0.0, 0, 1)
    @test_throws ArgumentError block_rayleigh_channel([1.0], 0.0, 1, -2)
end

@testset "Interleaving on block fading" begin
    # Paired comparison: one h and one noise draw, applied in channel order.
    # Coded length 640 = 16 × 40, block_length 16 equals the interleaver depth,
    # so one fade block is one column and deinterleaves to a stride of 40.
    # Eb = 2, Es = 1, nominal rate 1/2. Eb/N0 = 10 dB, 40 frames, seed RNG_SEED.
    Random.seed!(RNG_SEED)
    n_rows, n_cols = 16, 40
    block_length = n_rows
    n_coded = n_rows * n_cols
    n_info = n_coded ÷ 2 - 2
    nframes = 40
    snr_db = 10.0
    Eb = 2.0
    σ2 = Eb / (2 * snr_db_to_linear(snr_db))
    plain_tx = Int[]
    plain_rx = Int[]
    intl_rx = Int[]
    for _ in 1:nframes
        bits = rand(0:1, n_info)
        coded = conv_encode(bits)
        length(coded) == n_coded || error("coded length $(length(coded))")
        syms = ComplexF64.(modulate(coded, BPSK()))
        syms_i = block_interleave(syms, n_rows, n_cols)
        h = block_rayleigh(n_coded, block_length)
        noise = awgn(zeros(ComplexF64, n_coded), snr_db; energy_per_bit=Eb)
        llr = Vector{Float64}(undef, n_coded)
        llr_i = Vector{Float64}(undef, n_coded)
        y = h .* syms .+ noise
        y_i = h .* syms_i .+ noise
        @inbounds for i in 1:n_coded
            llr[i] = 2 * real(conj(h[i]) * y[i]) / σ2
            llr_i[i] = 2 * real(conj(h[i]) * y_i[i]) / σ2
        end
        append!(plain_tx, bits)
        append!(plain_rx, viterbi_decode(llr))
        append!(intl_rx, viterbi_decode(block_deinterleave(llr_i, n_rows, n_cols)))
    end
    nbits = length(plain_tx)
    k_plain = count(!=(0), plain_tx .- plain_rx)
    k_intl = count(!=(0), plain_tx .- intl_rx)
    lo_plain, _ = clopper_pearson(k_plain, nbits, BER_ALPHA)
    _, hi_intl = clopper_pearson(k_intl, nbits, BER_ALPHA)
    @test k_plain >= 50
    @test hi_intl < lo_plain
    @test ber(plain_tx, intl_rx) < ber(plain_tx, plain_rx)

    n_frame_err_plain = 0
    n_frame_err_intl = 0
    for f in 1:nframes
        idx = (f - 1) * n_info .+ (1:n_info)
        n_frame_err_plain += any(!=(0), @view(plain_tx[idx]) .- @view(plain_rx[idx]))
        n_frame_err_intl += any(!=(0), @view(plain_tx[idx]) .- @view(intl_rx[idx]))
    end
    @test fer(plain_tx, plain_rx, n_info) == n_frame_err_plain / nframes
    @test fer(plain_tx, intl_rx, n_info) == n_frame_err_intl / nframes
    lo_fp, _ = clopper_pearson(n_frame_err_plain, nframes, BER_ALPHA)
    _, hi_fi = clopper_pearson(n_frame_err_intl, nframes, BER_ALPHA)
    @test n_frame_err_plain >= 10
    @test hi_fi < lo_fp
    @test ser(plain_tx, intl_rx) == ber(plain_tx, intl_rx)
end
