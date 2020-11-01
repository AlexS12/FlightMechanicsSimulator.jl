using Test
using CSV
using DataFrames
using FlightMechanicsSimulator


df = DataFrame!(CSV.File("data/adc.csv"))
@testset "adc.jl" begin
    for case in eachrow(df)
        mach_1, qbar_1 = FlightMechanicsSimulator.F16.adc(case.vt, case.alt)
        @test isapprox(mach_1, case.mach)
        @test isapprox(qbar_1, case.qbar)
    end
end
