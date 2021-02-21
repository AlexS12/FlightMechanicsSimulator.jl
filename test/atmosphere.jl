using Test
using FlightMechanicsSimulator


@testset "ISA1978" begin
    # Test sea level
    T, p, ρ, a = atmosphere_isa(0.0)
    @test isapprox(T, 288.15)
    @test isapprox(p, 101325)
    @test isapprox(ρ, 1.225)
    @test isapprox(a, 340.29398, rtol=1e-7)

    # Test 0-11 Km
    h = [0.0, 50.0, 550.0, 6500.0, 10000.0, 11000.0]

    rv = atmosphere_isa.(h)

    T = [x[1] for x in rv]
    p = [x[2] for x in rv]
    ρ = [x[3] for x in rv]
    a = [x[4] for x in rv]

    @test isapprox(T, [288.150, 287.825, 284.575, 245.900, 223.150, 216.650])
    @test isapprox(ρ, [1.2250, 1.2191, 1.1616, 0.62384, 0.41271, 0.36392], rtol=1e-4)
    @test isapprox(a, [340.29, 340.10, 338.18, 314.36, 299.46, 295.07], rtol=1e-5)
end
