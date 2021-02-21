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


# The atmosphere in Stevens does not provide exact results when tested against official tables
@testset "ISA1978-Morelli" begin
    # Morelli, Eugene A., and Vladislav Klein. Aircraft system identification: Theory and
    # practice. Williamsburg, VA: Sunflyte Enterprises, 2016. Appendix C, p. 574.
    # Test sea level
    T, ρ, a, p = F16.atmosphere(0.0)  # 0 ft
    @test_broken isapprox(T, 518.67, atol=0.005)  # Rankine
    @test_broken isapprox(p, 2116.2, atol=0.05)  # psf
    @test_broken isapprox(ρ, 0.0023769, atol=5e-8)  # slug/ft³
    @test_broken isapprox(a, 1116.44, atol=0.005)  # ft/s
end
