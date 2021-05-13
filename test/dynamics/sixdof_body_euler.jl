# SixDOFAeroEuler is tested in f16 with the validation data.
# This is a check of sixdof body euler results against sixdof aero.

using Test

using FlightMechanicsUtils
using OrdinaryDiffEq

using FlightMechanicsSimulator


# Check conversion
dss_base = SixDOFAeroEuler([
    250 * FT2M,
    0.2392628,
    0.0005061803,
    1.366289,
    0.05000808,
    0.2340769,
    -0.01499617,
    0.2933811,
    0.06084932,
    0.0 * FT2M,
    0.0 * FT2M,
    0.0 * FT2M,
    64.12363,
])

ac = F16(F16Stevens.MASS, F16Stevens.INERTIA, 0.35)
atmosphere = F16StevensAtmosphere(get_height(dss_base))
gravity = LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M)

dss = SixDOFBodyEuler(dss_base)

# Check getter methods for DSState
getter_methods = [
    get_earth_position,
    get_height,
    get_euler_angles,
    get_tas,
    get_α,
    get_β,
    get_tasαβ,
    get_body_velocity,
    get_horizon_velocity,
    get_ang_vel_body,
    get_euler_angles_rates,
    get_engine_power,
    get_ang_vel_body,
]

@testset "$(gm)" for gm in getter_methods
    @test isapprox(gm(dss), gm(dss_base))
end

# Check getter methods for DSStateDot
args = [
    0,
    F16Stevens.MASS,
    F16Stevens.INERTIA,
    [5000, 6000, 7000],
    [2000, 3000, 4000],
    [F16Stevens.HX, F16Stevens.HX*0.1, F16Stevens.HX*0.5],
    0.1
]

dssd_base = state_eqs(dss_base, args...)
dssd = state_eqs(dss, args...)

getter_methods = [
    get_earth_position,
    get_height,
    get_euler_angles,
    get_tas,
    get_α,
    get_β,
    get_tasαβ,
    get_body_velocity,
    get_horizon_velocity,
    get_ang_vel_body,
    get_euler_angles_rates,
    get_engine_power,
    get_ang_vel_body,
    get_tas_dot,
    get_α_dot,
    get_β_dot,
    get_tasαβ_dot,
    get_accel_body,
    get_ang_accel_body,
    get_engine_power_dot,
]

@testset "$(gm)" for gm in getter_methods
    @test isapprox(gm(dssd), gm(dssd_base))
end

# Check trim results
@testset "Trim γ=$(cond[1]) ψ_dot=$(cond[2])" for cond in (
    (0.0, 0.0), (0.0, 0.3), (0.2, 0.0), (0.2, 0.3)
    )

    γ = cond[1]
    ψ_dot = cond[2]

    args = [
        [0.8349601, -1.481766, 0.09553108, -0.4118124],
        ac,
        atmosphere,
        gravity,
        γ,
        ψ_dot,
    ]

    dssd_base, controls_trim_base, outputs_trim_base, cost_base = trim(
        dss_base, args...
    )

    dssd, controls_trim, outputs_trim, cost = trim(
        dss_base, args...
    )

    @test isapprox(get_value.(controls_trim, 0.0), get_value.(controls_trim_base, 0.0))
    @test isapprox(get_α(dssd), get_α(dssd_base))
    @test isapprox(get_β(dssd), get_β(dssd_base))
end

# Check simulation
@testset "Simulation γ=$(cond[1]) ψ_dot=$(cond[2])" for cond in (
    (0.0, 0.0), (0.0, 0.3), (0.2, 0.0), (0.2, 0.3)
    )

    γ = cond[1]
    ψ_dot = cond[2]

    args_trim = [
        [0.8349601, -1.481766, 0.09553108, -0.4118124],
        ac,
        atmosphere,
        gravity,
        γ,
        ψ_dot,
    ]

    dssd_base_trim, controls_trim_base, outputs_trim_base, cost_base = trim(
        dss_base, args_trim...
    )

    t0 = 0.0
    t1 = 15.0
    dt = 0.01

    args_sim = [
        controls_trim_base,
        ac,
        F16StevensAtmosphere,
        gravity,
    ]

    kwargs_sim = Dict(:solver=>RK4(), :solve_args=>Dict(:reltol=>1e-10, :saveat=>dt))

    r_base = simulate(
        t0,
        t1,
        get_ds_state(dssd_base_trim),
        args_sim...;
        kwargs_sim...
    )

    r = simulate(
        t0,
        t1,
        SixDOFBodyEuler(get_ds_state(dssd_base_trim)),
        args_sim...;
        kwargs_sim...
    )

    # Check final state
    names = [Symbol("x$ii") for ii in 1:13]
    final_dss_base = SixDOFAeroEuler(Array(r_base[end, names]))
    final_dss = SixDOFBodyEuler(Array(r[end, names]))

    @test isapprox(get_tas(final_dss), get_tas(final_dss_base), atol=1e-5)
    @test isapprox(get_α(final_dss), get_α(final_dss_base), atol=1e-5)
    @test isapprox(get_β(final_dss), get_β(final_dss_base), atol=1e-5)
    @test isapprox(get_body_velocity(final_dss), get_body_velocity(final_dss_base), atol=1e-4)
    @test isapprox(get_ang_vel_body(final_dss), get_ang_vel_body(final_dss_base), atol=1e-5)
    @test isapprox(get_euler_angles(final_dss), get_euler_angles(final_dss_base), atol=1e-5)
    @test isapprox(get_earth_position(final_dss), get_earth_position(final_dss_base), atol=1e-3)
end
