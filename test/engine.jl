using Test
using FlightMechanicsSimulator


@testset "engine.jl" begin

    for thtl in LinRange(0, 1, 20)
        rv1 = FlightMechanicsSimulator.Fortran.tgear(thtl)
        rv2 = FlightMechanicsSimulator.tgear(thtl)
        @test isapprox(rv1, rv2, atol = 1.0e-15)
    end

    for p3 in LinRange(0, 100, 50), p1 in LinRange(0, 100, 50)
        rv1 = FlightMechanicsSimulator.Fortran.pdot(p3, p1)
        rv2 = FlightMechanicsSimulator.pdot(p3, p1)
        @test isapprox(rv1, rv2, atol = 1.0e-15)
    end

    for dp in LinRange(0, 100, 50)
        rv1 = FlightMechanicsSimulator.Fortran.rtau(dp)
        rv2 = FlightMechanicsSimulator.rtau(dp)
        @test isapprox(rv1, rv2, atol = 1.0e-15)
    end

    for pow in LinRange(0, 100, 50),
        alt in LinRange(0, 50000, 500),
        rmach in LinRange(0, 1, 25)

        rv1 = FlightMechanicsSimulator.Fortran.thrust(pow, alt, rmach)
        rv2 = FlightMechanicsSimulator.thrust(pow, alt, rmach)
        @test isapprox(rv1, rv2, atol = 1.0e-10)
    end

end