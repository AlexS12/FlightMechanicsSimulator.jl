using Test

using DataFrames
using FlightMechanicsUtils
using OrdinaryDiffEq

using FlightMechanicsSimulator


# ---------- PROPAGATE COORDINATED TURN ----------
# Check that this is a trimmed condition
# Stevens, B. L., Lewis, F. L., & Johnson, E. N. (2015). Aircraft control
# and simulation: dynamics, controls design, and autonomous systems. John Wiley
# & Sons.
# Example 3.6-2 (page 191)

# Trimmed conditions: page 192
# C     X(1)  -> vt (m/s)
# C     X(2)  -> alpha (rad)
# C     X(3)  -> beta (rad)
# C     X(4)  -> phi (rad)
# C     X(5)  -> theta (rad)
# C     X(6)  -> psi (rad)
# C     X(7)  -> P (rad/s)
# C     X(8)  -> Q (rad/s)
# C     X(9)  -> R (rad/s)
# C     X(10) -> North (m)
# C     X(11) -> East (m)
# C     X(12) -> Altitude (m)
# C     X(13) -> Pow
x_stev = [
    502.0 * FT2M,
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
]
controls_stev = [0.8349601, -1.481766, 0.09553108, -0.4118124]
xcg = 0.35

# RETRIM to refine flying condition
dssd, controls_trim, outputs_trim, cost = trim(
    SixDOFAeroEuler(x_stev),
    controls_stev,
    F16(F16Stevens.MASS, F16Stevens.INERTIA, xcg),
    F16StevensAtmosphere(x_stev[12]),
    LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M),
    0.0,
    0.3,
)

x_trim = get_x(dssd)

dt = 0.01  # s
t0 = 0.0  # s
t1 = 180.0  # s

results = simulate(
    t0,
    t1,
    SixDOFAeroEuler(x_trim),
    controls_trim,
    F16(F16Stevens.MASS, F16Stevens.INERTIA, xcg),
    F16StevensAtmosphere,
    LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M);
    solver=RK4(),
    solve_args=Dict(:reltol=>1e-6),
    )

sol, log_out = FlightMechanicsSimulator._simulate(
    t0,
    t1,
    SixDOFAeroEuler(x_trim),
    controls_trim,
    F16(F16Stevens.MASS, F16Stevens.INERTIA, xcg),
    F16StevensAtmosphere,
    LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M);
    solver=RK4(),
    solve_args=Dict(:reltol=>1e-6),
    )

@testset "SimulationLogs" begin
    for (ii, name) in enumerate(get_x_names(dssd))
        ode_sol = [sol.u[jj][ii] for jj in 1:size(sol.t)[1]]
        @test isapprox(
            ode_sol - results[!, name],
            zeros(size(sol.t))
        )
    end
end
