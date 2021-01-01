using Test
using CSV
using DataFrames
using FlightMechanicsSimulator


df = DataFrame!(CSV.File("data/adc.csv"))
@testset "adc.jl" begin
for case in eachrow(df)
        T, ρ, a, p = F16.atmosphere(case.alt)
        mach_1, qbar_1 = F16.adc(case.vt, T, ρ, a, p)
        @test isapprox(mach_1, case.mach)
        @test isapprox(qbar_1, case.qbar)
    end
end
