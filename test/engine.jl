using Test
using CSV
using DataFrames
using FlightMechanicsSimulator


@testset "engine.jl" begin

    df = DataFrame!(CSV.File("data/tgear.csv"))
    for case in eachrow(df)
        rv1 = F16.tgear(case.thtl)
        @test isapprox(rv1, case.tgear, atol = 1.0e-15)
    end

    df = DataFrame!(CSV.File("data/pdot.csv"))
    for case in eachrow(df)
        rv1 = F16.pdot(case.p3, case.p1)
        @test isapprox(rv1, case.pdot, atol = 1.0e-15)
    end

    df = DataFrame!(CSV.File("data/rtau.csv"))
    for case in eachrow(df)
        rv1 = F16.rtau(case.dp)
        @test isapprox(rv1, case.rtau, atol = 1.0e-15)
    end

    df = DataFrame!(CSV.File("data/thrust.csv"))
    for case in eachrow(df)
        rv1 = F16.thrust(case.pow, case.alt, case.rmach)
        @test isapprox(rv1, case.thrust, atol = 1.0e-10)
    end

end
