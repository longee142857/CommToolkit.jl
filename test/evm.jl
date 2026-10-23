# RMS EVM. Percent is 100 times the amplitude ratio, not the power ratio.
# Decibels are 20*log10 of that ratio (equivalently 10*log10 of the power ratio).

@testset "EVM" begin
    ref = [1.0, -1.0, 1.0, -1.0]
    @test evm(ref, ref) == 0.0
    @test evm(ref, ref; unit=:linear) == 0.0
    @test evm(ref, ref; unit=:percent) == 0.0
    @test evm(ref, ref; unit=:db) == -Inf
    @test evm([1, -1], [1, -1]) == 0.0

    # |1.1 s - s| / |s| = 0.1. 10% is -20 dB, not -10 dB.
    rx = 1.1 .* ref
    @test evm(ref, rx) ≈ 10.0 atol = 1e-12
    @test evm(ref, rx; unit=:linear) ≈ 0.1 atol = 1e-12
    @test evm(ref, rx; unit=:db) ≈ -20.0 atol = 1e-12
    @test evm(2 .* ref, 2 .* rx) ≈ evm(ref, rx) atol = 1e-12

    # A common phase is not fitted out.
    spun = ComplexF64[-1, -1]
    @test evm(ComplexF64[1, 1], spun; unit=:linear) ≈ 2.0 atol = 1e-12
    @test evm(ComplexF64[1, 1], spun) ≈ 200.0 atol = 1e-12

    # One quadrature error. Reference energy 1, error energy 0.01.
    @test evm(ComplexF64[1], ComplexF64[1 + 0.1im]; unit=:linear) ≈ 0.1 atol = 1e-12

    # A zero symbol does not by itself make the ratio undefined.
    @test evm([0.0, 1.0], [0.0, 1.0]) == 0.0
    @test evm([0.0, 1.0], [0.5, 1.0]; unit=:linear) ≈ 0.5 atol = 1e-12

    # Two unit references, one sign flip: error energy 4, reference energy 2.
    @test evm([1.0, 1.0], [1.0, -1.0]; unit=:linear) ≈ sqrt(2) atol = 1e-12

    @test_throws ArgumentError evm([1.0, -1.0], [1.0])
    @test_throws ArgumentError evm(Float64[], Float64[])
    @test_throws ArgumentError evm([0.0, 0.0], [0.0, 1.0])
    @test_throws ArgumentError evm([0.0, 0.0], [0.0, 0.0])
    @test_throws ArgumentError evm([1.0], [NaN])
    @test_throws ArgumentError evm([Inf], [1.0])
    @test_throws ArgumentError evm([1.0], [1.0]; unit=:power)
    @test_throws ArgumentError evm([1.0], [1.0]; unit=:DB)
end
