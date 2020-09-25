using Test
using FlightMechanicsSimulator


@testset "aero.jl" begin

    α_test = LinRange(-10, 45, 20)
    β_test = LinRange(-30, 30, 20)

    de_test = LinRange(-25, 25, 20)
    da_test = LinRange(-21.5, 21.5, 20)
    dr_test = LinRange(-30, 30, 20)

    for α in α_test
        rv1 = FlightMechanicsSimulator.Fortran.damp(α)
        rv2 = FlightMechanicsSimulator.damp(α)
        @test isapprox(rv1, rv2, atol = 1.0e-14)
    end

    for α in α_test, de in de_test
        rv1 = FlightMechanicsSimulator.Fortran.CX(α, de)
        rv2 = FlightMechanicsSimulator.CX(α, de)
        @test isapprox(rv1, rv2, atol = 1.0e-15)
    end


    for β in β_test, da in da_test, dr in dr_test
        rv1 = FlightMechanicsSimulator.Fortran.CY(β, da, dr)
        rv2 = FlightMechanicsSimulator.CY(β, da, dr)
        @test isapprox(rv1, rv2, atol = 1.0e-15)
    end


    for α in α_test, β in β_test, de in de_test
        rv1 = FlightMechanicsSimulator.Fortran.CZ(α, β, de)
        rv2 = FlightMechanicsSimulator.CZ(α, β, de)
        @test isapprox(rv1, rv2, atol = 1.0e-15)
    end


    for α in α_test, de in de_test
        rv1 = FlightMechanicsSimulator.Fortran.CM(α, de)
        rv2 = FlightMechanicsSimulator.CM(α, de)
        @test isapprox(rv1, rv2, atol = 1.0e-15)
    end


    for α in α_test, β in de_test
        rv1 = FlightMechanicsSimulator.Fortran.CL(α, β)
        rv2 = FlightMechanicsSimulator.CL(α, β)
        @test isapprox(rv1, rv2, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.Fortran.CN(α, β)
        rv2 = FlightMechanicsSimulator.CN(α, β)
        @test isapprox(rv1, rv2, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.Fortran.DLDA(α, β)
        rv2 = FlightMechanicsSimulator.DLDA(α, β)
        @test isapprox(rv1, rv2, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.Fortran.DLDR(α, β)
        rv2 = FlightMechanicsSimulator.DLDR(α, β)
        @test isapprox(rv1, rv2, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.Fortran.DNDA(α, β)
        rv2 = FlightMechanicsSimulator.DNDA(α, β)
        @test isapprox(rv1, rv2, atol = 1.0e-15)

        rv1 = FlightMechanicsSimulator.Fortran.DNDR(α, β)
        rv2 = FlightMechanicsSimulator.DNDR(α, β)
        @test isapprox(rv1, rv2, atol = 1.0e-15)
    end

end
