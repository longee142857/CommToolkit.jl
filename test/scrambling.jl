# Additive LFSR. The oracle below is the linear recurrence in the docstring,
# not a second copy of the shift. Named polynomials are the classical
# trinomials; the bits are not the inverted ITU-T O.150 test pattern.

function _recurrence(poly::Integer, nbits::Int, seed::Integer)
    n = ndigits(poly, base=2) - 1
    s = Vector{Int}(undef, max(nbits, n))
    for k in 0:(n - 1)
        s[k + 1] = (Int(seed) >> k) & 1
    end
    for idx in (n + 1):nbits
        t = idx - n - 1
        acc = 0
        for k in 0:(n - 1)
            if ((Int(poly) >> k) & 1) == 1
                acc ⊻= s[t + k + 1]
            end
        end
        s[idx] = acc
    end
    return s[1:nbits]
end

function _maximal_period(poly::Integer)
    n = ndigits(poly, base=2) - 1
    period = (1 << n) - 1
    seq = lfsr_sequence(poly, period + n)
    seq[(period + 1):(period + n)] == seq[1:n] || return false
    seen = Set{UInt64}()
    sizehint!(seen, period)
    acc = UInt64(0)
    for i in 1:n
        acc |= UInt64(seq[i]) << (i - 1)
    end
    push!(seen, acc)
    for i in 1:(period - 1)
        acc = (acc >> 1) | (UInt64(seq[i + n]) << (n - 1))
        push!(seen, acc)
    end
    return length(seen) == period && !in(zero(UInt64), seen)
end

@testset "LFSR scrambling" begin
    @test prbs_polynomial(7) == (1 << 7) | (1 << 6) | 1
    @test prbs_polynomial(9) == (1 << 9) | (1 << 5) | 1
    @test prbs_polynomial(11) == (1 << 11) | (1 << 9) | 1
    @test prbs_polynomial(15) == (1 << 15) | (1 << 14) | 1
    @test prbs_polynomial(20) == (1 << 20) | (1 << 3) | 1
    @test prbs_polynomial(23) == (1 << 23) | (1 << 18) | 1
    @test prbs_polynomial(29) == (1 << 29) | (1 << 27) | 1
    @test prbs_polynomial(31) == (1 << 31) | (1 << 28) | 1

    # All-ones seed: first n bits are 1. A trinomial then feeds back a 0.
    for d in (7, 9, 11, 15, 20, 23, 29, 31)
        poly = prbs_polynomial(d)
        @test ndigits(poly, base=2) - 1 == d
        seq = lfsr_sequence(poly, d + 1)
        @test seq[1:d] == ones(Int, d)
        @test seq[d + 1] == 0
        bits = [1, 0, 1, 1, 0, 0, 1, 0, 1]
        @test scramble(bits, poly) isa Vector{Int}
        @test descramble(scramble(bits, poly), poly) == bits
        @test scramble(scramble(bits, poly; seed=d), poly; seed=d) == bits
    end

    prbs7 = prbs_polynomial(7)
    prefix7 = [1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1]
    @test lfsr_sequence(prbs7, length(prefix7)) == prefix7
    @test lfsr_sequence(prbs7, 80) == _recurrence(prbs7, 80, (1 << 7) - 1)
    @test lfsr_sequence(prbs7, 40; seed=0b0010110) == _recurrence(prbs7, 40, 0b0010110)
    period7 = lfsr_sequence(prbs7, 127)
    @test count(==(1), period7) == 64
    @test count(==(0), period7) == 63
    @test _maximal_period(prbs7)
    @test _maximal_period(prbs_polynomial(9))
    @test _maximal_period(prbs_polynomial(11))
    @test _maximal_period(prbs_polynomial(15))

    # x^4 + x + 1, seed all ones. Hand recurrence, period 15.
    poly4 = (1 << 4) | (1 << 1) | 1
    prefix4 = [1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0]
    @test lfsr_sequence(poly4, 15) == prefix4
    @test lfsr_sequence(poly4, 15 + 4)[16:19] == prefix4[1:4]
    @test lfsr_sequence(poly4, 30) == _recurrence(poly4, 30, (1 << 4) - 1)
    @test _maximal_period(poly4)

    # Not primitive, still an involution. x^4 + 1.
    poly_np = (1 << 4) | 1
    bits = [1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0]
    @test descramble(scramble(bits, poly_np; seed=0b1010), poly_np; seed=0b1010) == bits
    @test scramble(scramble(bits, poly_np), poly_np) == bits

    # Different seeds are different phases. The wrong seed does not invert.
    long = [isodd(i) ? 1 : 0 for i in 1:64]
    a = scramble(long, prbs7; seed=1)
    b = scramble(long, prbs7; seed=3)
    @test a != b
    @test descramble(a, prbs7; seed=1) == long
    @test descramble(a, prbs7) != long
    @test descramble(a, prbs_polynomial(9); seed=1) != long

    held = [1, 0, 1, 0, 1]
    snapshot = copy(held)
    @test scramble(held, prbs7) != held
    @test held == snapshot
    @test scramble(Bool[true, false, true, true], prbs7) ==
          scramble([1, 0, 1, 1], prbs7)
    @test scramble(Int[], prbs7) == Int[]
    @test lfsr_sequence(prbs7, 0) == Int[]
    @test lfsr_sequence(prbs7, 5) == lfsr_sequence(prbs7, 5; seed=(1 << 7) - 1)

    # Degree 62 still shifts inside UInt64.
    poly62 = (Int(1) << 62) | 1
    @test descramble(scramble(long, poly62), poly62) == long
    @test lfsr_sequence(poly62, 8; seed=big(5)) == lfsr_sequence(poly62, 8; seed=5)

    @test_throws ArgumentError prbs_polynomial(8)
    @test_throws ArgumentError prbs_polynomial(0)
    @test_throws ArgumentError prbs_polynomial(-3)
    @test_throws ArgumentError lfsr_sequence(0, 4)
    @test_throws ArgumentError lfsr_sequence(-4, 4)
    @test_throws ArgumentError lfsr_sequence(1, 4)
    @test_throws ArgumentError lfsr_sequence(2, 4)
    @test_throws ArgumentError lfsr_sequence((BigInt(1) << 63) | 1, 1)
    @test_throws ArgumentError lfsr_sequence(2, 0)
    @test_throws ArgumentError lfsr_sequence(prbs7, -1)
    @test_throws ArgumentError lfsr_sequence(prbs7, 4; seed=0)
    @test_throws ArgumentError lfsr_sequence(prbs7, 4; seed=-2)
    @test_throws ArgumentError lfsr_sequence(prbs7, 4; seed=1 << 7)
    @test_throws ArgumentError lfsr_sequence(prbs7, 4; seed=big(10)^20)
    @test_throws ArgumentError scramble([0, 1, 2], prbs7)
    @test_throws ArgumentError scramble([0, -1], prbs7)
end
