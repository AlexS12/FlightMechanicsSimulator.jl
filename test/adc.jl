using Test
using CSV
using DataFrames
using FlightMechanicsSimulator


df = DataFrame(CSV.File("data/adc.csv"))
@testset "adc.jl" begin
for case in eachrow(df)
        T, ρ, a, p = F16.atmosphere(case.alt)
        mach_1, qbar_1 = F16.adc(case.vt, T, ρ, a, p)
        @test isapprox(mach_1, case.mach)
        @test isapprox(qbar_1, case.qbar)
    end
end

@testset "atmosphere_f16" begin
    df = DataFrame(CSV.File("data/adc.csv"))
    # df has duplicate altitudes with different vt
    df = unique(df, :alt)
    # Discontinuity alt (35000 ft) produces tests not passing because stevens atmosphere
    # is not continous at 35000 ft, so if conversion FT2M produces slightly different altitude,
    # a wrong atmosphere layer will be chosen
    df = df[(!isapprox).(df[:, :alt], 35000), :]
    for case in eachrow(df)
        T1, ρ1, a1, p1 = F16.atmosphere_f16(case.alt * FT2M)
        T2, ρ2, a2, p2 = F16.atmosphere(case.alt)
        @test isapprox(T1*KEL2RANK, T2)
        @test isapprox(ρ1*KGM32SLUGFT3, ρ2)
        @test isapprox(a1*M2FT, a2)
        @test isapprox(p1*PA2PSF, p2)
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
    # Test in tropopause
    T, ρ, a, p = F16.atmosphere(50000.0)  # 0 ft
    @test_broken isapprox(T, 389.97, atol=0.005)  # Rankine
    @test_broken isapprox(p, 243.6, atol=0.05)  # psf
    @test_broken isapprox(ρ, 0.0003639, atol=5e-8)  # slug/ft³
    @test_broken isapprox(a, 968.08, atol=0.005)  # ft/s
end
