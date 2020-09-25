using Test
using FlightMechanicsSimulator


@testset "adc.jl" begin

    for vt = 0.0:10.0:1000.0, alt = 0.0:1000.0:50000.0
        mach_1, qbar_1 = FlightMechanicsSimulator.Fortran.adc(vt, alt)
        mach_2, qbar_2 = FlightMechanicsSimulator.adc(vt, alt)
        @test isapprox(mach_1, mach_2, atol = 1.0e-16)
        @test isapprox(qbar_1, qbar_2)
    end

end
