using Test
using CSV
using DataFrames
using FlightMechanicsSimulator


@testset "aero.jl" begin

    α_test = LinRange(-10, 45, 20)
    β_test = LinRange(-30, 30, 20)

    de_test = LinRange(-25, 25, 20)
    da_test = LinRange(-21.5, 21.5, 20)
    dr_test = LinRange(-30, 30, 20)

    df = DataFrame!(CSV.File("data/damp.csv"))
    for case in eachrow(df)
        rv1 = FlightMechanicsSimulator.damp(case.alpha)
        @test isapprox(rv1, Array(case[2:end]), atol = 1.0e-14)
    end

    df = DataFrame!(CSV.File("data/cx.csv"))
    for case in eachrow(df)
        rv1 = FlightMechanicsSimulator.CX(case.alpha, case.de)
        @test isapprox(rv1, case.cx, atol = 1.0e-15)
    end

    df = DataFrame!(CSV.File("data/cy.csv"))
    for case in eachrow(df)
        rv1 = FlightMechanicsSimulator.CY(case.beta, case.da, case.dr)
        @test isapprox(rv1, case.cy, atol = 1.0e-15)
    end

    df = DataFrame!(CSV.File("data/cz.csv"))
    for case in eachrow(df)
        rv1 = FlightMechanicsSimulator.CZ(case.alpha, case.beta, case.de)
        @test isapprox(rv1, case.cz, atol = 1.0e-15)
    end

    df = DataFrame!(CSV.File("data/cm.csv"))
    for case in eachrow(df)
        rv1 = FlightMechanicsSimulator.CM(case.alpha, case.de)
        @test isapprox(rv1, case.cm, atol = 1.0e-15)
    end

    df = DataFrame!(CSV.File("data/aero_coeffs.csv"))
    for case in eachrow(df)
        rv1 = FlightMechanicsSimulator.CL(case.alpha, case.beta)
        @test isapprox(rv1, case.cl, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.CN(case.alpha, case.beta)
        @test isapprox(rv1, case.cn, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.DLDA(case.alpha, case.beta)
        @test isapprox(rv1, case.dlda, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.DLDR(case.alpha, case.beta)
        @test isapprox(rv1, case.dldr, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.DNDA(case.alpha, case.beta)
        @test isapprox(rv1, case.dnda, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.DNDR(case.alpha, case.beta)
        @test isapprox(rv1, case.dndr, atol = 1.0e-15)
    end
end
