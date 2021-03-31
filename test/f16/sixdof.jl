using Test
using CSV
using DataFrames
using StaticArrays

using FlightMechanicsSimulator
using FlightMechanicsUtils


df = DataFrame(CSV.File("data/sixdof.csv"))
# Test file was wrong before and longitudinal load factor
# needs to be corrected
df[!, "o3"] = df[!, "o3"] / (FlightMechanicsSimulator.F16Stevens.MASS * KG2LB / F16Stevens.GD)

for case in eachrow(df)
    x = Array(case[["x$ii" for ii in 1:13]])

    x[1] = x[1] * FT2M
    x[12] = x[12] * FT2M

    controls = Array(case[["c$ii" for ii in 1:4]])

    xd1, outputs1 =
        f(
            case.time,
            x,
            SixDOFAeroEuler(SVector{12}(x[1:12])),
            controls,
            F16(F16Stevens.MASS, F16Stevens.INERTIA, case.xcg),
            F16StevensAtmosphere(x[12]),
            LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M),
        )

    xd1[1] *= M2FT
    xd1[10] *= M2FT
    xd1[11] *= M2FT
    xd1[12] *= M2FT

    outputs1[4] *= PA2PSF

    @test isapprox(
        xd1,
        Array(case[["x$(ii)_dot" for ii in 1:13]]),
        nans = true,
        atol = 1e-10,
    )

    @test isapprox(
        outputs1,
        Array(case[["o$(ii)" for ii in 1:7]]),
        nans = true,
        atol = 1e-10,
    )
end
