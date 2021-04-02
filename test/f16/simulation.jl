using Test

using CSV
using DataFrames
using FlightMechanicsUtils
using OrdinaryDiffEq
using StaticArrays

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
x_dot, outputs = f(
    time,
    SixDOFAeroEuler(x_stev),
    controls_stev,
    F16(F16Stevens.MASS, F16Stevens.INERTIA, xcg),
    F16StevensAtmosphere(x_stev[12]),
    LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M),
)

# Linear acceleration
@test isapprox(x_dot[1:3], zeros(3), atol = 5e-4)
# Angular acceleration
@test isapprox(x_dot[7:9], zeros(3), atol = 1e-5)
# Cost function
cost =
    x_dot[1]^2 +
    100 * (x_dot[2]^2 + x_dot[3]^2) +
    10 * (x_dot[7]^2 + x_dot[8]^2 + x_dot[9]^2)


@test isapprox(cost, 0.0, atol = 1e-5)

# Desired turn rate
@test isapprox(x_dot[6], 0.3, atol = 1e-6)
# No vertical changes
@test isapprox(x_dot[12], 0.0, atol = 1e-4)

# RETRIM to refine flying condition
x_trim, controls_trim, x_dot_trim, outputs_trim, cost = trim(
    SixDOFAeroEuler(x_stev),
    controls_stev,
    F16(F16Stevens.MASS, F16Stevens.INERTIA, xcg),
    F16StevensAtmosphere(x_stev[12]),
    LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M),
    0.0,
    0.3,
)

# Linear acceleration
@test isapprox(x_dot_trim[1:3], zeros(3), atol=1e-12)
# Angular acceleration
@test isapprox(x_dot_trim[7:9], zeros(3), atol=1e-12)
# Desired turn rate
@test isapprox(x_dot_trim[6], 0.3, atol = 1e-6)
# No vertical changes
@test isapprox(x_dot_trim[12], 0.0, atol = 1e-4)

dt = 0.01  # s
t0 = 0.0  # s
t1 = 180.0  # s

x = x_trim
# Transform to Input for simulate
controls = ConstantInput.(controls_trim)

results = simulate(
    t0,
    t1,
    SixDOFAeroEuler(x),
    controls,
    F16(F16Stevens.MASS, F16Stevens.INERTIA, xcg),
    F16StevensAtmosphere,
    LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M);
    solver=RK4(),
    solve_args=Dict(:reltol=>1e-10, :saveat=>dt),
    )

# Check X, Y against Stevens
# Stevens, B. L., Lewis, F. L., & Johnson, E. N. (2015). Aircraft control
# and simulation: dynamics, controls design, and autonomous systems. John Wiley
# & Sons.
# Example 3.6-6 (page 198)
xy_trajectory_data = [
    # time (s)  X (ft)  Y(ft)
    0.00e0      0.00e0    0.00e0
    1.00e1      2.36e2    3.33e3
    2.00e1     -4.68e2    6.65e1
    3.00e1      6.90e2    3.20e3
    4.00e1     -8.97e2    2.61e2
    5.00e1      1.09e3    2.94e3
    6.00e1     -1.26e3    5.68e2
    7.00e1      1.40e3    2.59e3
    8.00e1     -1.51e3    9.62e2
    9.00e1      1.60e3    2.16e3
    1.00e2     -1.66e3    1.41e3
    1.10e2      1.67e3    1.70e3
    1.20e2     -1.66e3    1.89e3
    1.30e2      1.61e3    1.22e3
    1.40e2     -1.53e3    2.35e3
    1.50e2      1.42e3    7.87e2
    1.60e2     -1.28e3    2.76e3
    1.70e2      1.11e3    4.22e2
    1.80e2     -9.21e2    3.07e3
 ]

 for case in eachrow(xy_trajectory_data)
    r = results[findall(in(case[1]), results.time), :]
    @test isapprox(case[2], r[1, :x] * M2FT, atol=20)
    @test isapprox(case[3], r[1, :y] * M2FT, atol=20)
 end

# Select last time step state
 x = results[end, :]

# Check that TAS, α, β, θ, ϕ, p, q, r, alt, pow remain constant
@test isapprox(x[:tas], x_trim[1])  # TAS
@test isapprox(x[:α], x_trim[2])  # α
@test isapprox(x[:β], x_trim[3])  # β
@test isapprox(x[:ϕ], x_trim[4])  # ϕ
@test isapprox(x[:θ], x_trim[5])  # θ
@test isapprox(x[:p], x_trim[7])  # p
@test isapprox(x[:q], x_trim[8])  # q
@test isapprox(x[:r], x_trim[9])  # r
@test isapprox(x[:z], x_trim[12], atol=1e-12)  # alt
@test isapprox(x[:pow], x_trim[13])  # pow
