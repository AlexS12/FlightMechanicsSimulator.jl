using Test
using CSV
using DataFrames
using FlightMechanicsSimulator


t = 0.0;

vt_test = [50.0, 75.0]  # ft/s
α_test = deg2rad.([1.0, 5.0])
β_test = deg2rad.([-5.0, 5.0])
ϕ_test = deg2rad.([-30.0, 25.0])
θ_test = deg2rad.([-15.0, 25.0])
ψ_test = deg2rad.([45.0, 175.0])
p_test = deg2rad.([-15.0, 30.0])
q_test = deg2rad.([-5.0, 10.0])
r_test = deg2rad.([-20.0, 30.0])
norh_ft = 0.0  # ft
east_ft = 0.0  # ft
alt_test = [5000.0, 45000.0]  # ft
pow_test = [10., 80.]  # %

de_test = [-25.0, 20.0]  # deg
da_test = [-15.0, 10.0]  # deg
dr_test = [-25.0, 20.0]  # deg
thtl_test = [0.2, 0.6]

xcg_test = [0.35, 0.25]

df = DataFrame!(CSV.File("data/sixdof.csv"))

for case in eachrow(df)
    x = Array(case[["x$ii" for ii in 1:13]])

    controls = Array(case[["c$ii" for ii in 1:4]])

    xd1, outputs1 =
        FlightMechanicsSimulator.f(
            case.time,
            x,
            case.xcg,
            controls,
        )

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
