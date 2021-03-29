using Test
using FlightMechanicsSimulator
using FlightMechanicsUtils


#  Stevens, B. L., Lewis, F. L., & Johnson, E. N. (2015). Aircraft control
#  and simulation: dynamics, controls design, and autonomous systems. John Wiley
#  & Sons. (page 193 table 3.6-2)
trim_test_data = [
#   TAS   thtl    AOA     DE       thtl_tol    AOA_tol    DE_tol
#   ft/s  unit    deg     deg      unit        deg        deg
    130   0.816   45.6    20.1     0.0005      0.05       0.15
    140   0.736   40.3   -1.36     0.001       0.05       0.05
    150   0.619   34.6    0.173    0.0005      0.05       0.05
    170   0.464   27.2    0.621    0.001       0.05       0.05
    640   0.23    0.742  -0.871    0.0005      0.015      0.0005
    800   0.378   -0.045 -0.943    0.0005      0.001      0.001
    200   0.287   19.7    0.723    0.0005      0.05       0.05
    260   0.148   11.6   -0.09     0.0005      0.05       0.05
    300   0.122   8.49   -0.591    0.0005      0.01       0.005
    350   0.107   5.87   -0.539    0.001       0.005      0.005
    400   0.108   4.16   -0.591    0.0005      0.005      0.005
    440   0.113   3.19   -0.671    0.0005      0.005      0.005
    500   0.137   2.14   -0.756    0.001       0.01       0.005
    540   0.16    1.63   -0.798    0.0005      0.005      0.005
    600   0.2     1.04   -0.846    0.0005      0.01       0.005
    700   0.282   0.382  -0.9      0.0005      0.001      0.0005
]

xcg = 0.35

for case in eachrow(trim_test_data)
    local x, controls, x_trim, controls_trim, x_dot_trim, outputs_trim, cost
    x = [
        case[1]*FT2M,  #-> vt (m/s)
        deg2rad(10.),  # -> alpha (rad)
        0.0,  # -> beta (rad)
        0.0,  # -> phi (rad)
        deg2rad(10.),  # -> theta (rad)
        0.0,  #  -> psi (rad)
        0.0,  # -> P (rad/s)
        0.0,  # -> Q (rad/s)
        0.0,  # -> R (rad/s)
        0.0,  # -> North (m)
        0.0,  # -> East (m)
        0.0,  # -> Altitude (m)
        50.0,  # -> Pow
    ]

    controls = [
        0.5,  # thtl
        0.0,  # elev
        0.0,  # ail
        0.0,  # rudder
    ]

    # TRIM
    x_trim, controls_trim, x_dot_trim, outputs_trim, cost = F16.trimmer(
        f,
        x,
        controls,
        F16.F16Stevens(F16.MASS, F16.INERTIA, xcg),
        F16StevensAtmosphere,
        LHDownGravity(FlightMechanicsSimulator.F16.GD*FT2M),
        0.0,
        0.0,
    )

    @test isapprox(cost, zeros(6), atol=1e-12)
    @test isapprox(controls_trim[1], case[2], atol=case[5])  # THTL
    @test isapprox(rad2deg(x_trim[2]), case[3], atol=case[6])  # AOA
    @test isapprox(controls_trim[2], case[4], atol=case[7])  # DE
end

#  Stevens, B. L., Lewis, F. L., & Johnson, E. N. (2015). Aircraft control
#  and simulation: dynamics, controls design, and autonomous systems. John Wiley
#  & Sons. (page 195 table 3.6-3)

# NOMINAL (first column)
xcg = 0.35

x = [
    502*FT2M,  #-> vt (m/s)
    deg2rad(10.),  # -> alpha (rad)
    0.0,  # -> beta (rad)
    0.0,  # -> phi (rad)
    deg2rad(10.),  # -> theta (rad)
    0.0,  #  -> psi (rad)
    0.0,  # -> P (rad/s)
    0.0,  # -> Q (rad/s)
    0.0,  # -> R (rad/s)
    0.0,  # -> North (m)
    0.0,  # -> East (m)
    0.0,  # -> Altitude (m)
    50.0,  # -> Pow
]

controls = [
    0.5,  # thtl
    0.0,  # elev
    0.0,  # ail
    0.0,  # rudder
]

x_trim, controls_trim, x_dot_trim, outputs_trim, cost = F16.trimmer(
        f,
        x,
        controls,
        F16.F16Stevens(F16.MASS, F16.INERTIA, xcg),
        F16StevensAtmosphere,
        LHDownGravity(FlightMechanicsSimulator.F16.GD*FT2M),
        0.0,
        0.0,
    )

@test isapprox(cost, zeros(6), atol=1e-12)
@test isapprox(x_trim[2], 0.03691, atol=0.00005)  # AOA
@test isapprox(x_trim[3], -4e-9, atol=1e-8)  # AOS
@test isapprox(x_trim[4], 0)  # PHI
@test isapprox(x_trim[5], 0.03691, atol=0.00005)  # THETA
@test isapprox(x_trim[7], 0)  # P
@test isapprox(x_trim[8], 0)  # Q
@test isapprox(x_trim[9], 0)  # R
@test isapprox(controls_trim[1], 0.1385, atol=0.0001)  # THTL
@test isapprox(controls_trim[2], -0.7588, atol=0.0002)  # DE
@test isapprox(controls_trim[3], -1.2e-7, atol=1e-6)  # DA
@test isapprox(controls_trim[4], 6.2e-7, atol=1e-6)  # DR

# XCG = 0.3 (second column)
xcg = 0.3

x = [
    502*FT2M,  #-> vt (m/s)
    deg2rad(10.),  # -> alpha (rad)
    0.0,  # -> beta (rad)
    0.0,  # -> phi (rad)
    deg2rad(10.),  # -> theta (rad)
    0.0,  #  -> psi (rad)
    0.0,  # -> P (rad/s)
    0.0,  # -> Q (rad/s)
    0.0,  # -> R (rad/s)
    0.0,  # -> North (m)
    0.0,  # -> East (m)
    0.0,  # -> Altitude (m)
    50.0,  # -> Pow
]

controls = [
    0.5,  # thtl
    0.0,  # elev
    0.0,  # ail
    0.0,  # rudder
]

x_trim, controls_trim, x_dot_trim, outputs_trim, cost = F16.trimmer(
        f,
        x,
        controls,
        F16.F16Stevens(F16.MASS, F16.INERTIA, xcg),
        F16StevensAtmosphere,
        LHDownGravity(FlightMechanicsSimulator.F16.GD*FT2M),
        0.0,
        0.0,
    )

@test isapprox(cost, zeros(6), atol=1e-12)
@test isapprox(x_trim[2], 0.03936, atol=0.00005)  # AOA
@test isapprox(x_trim[3], 4.1e-9, atol=1e-8)  # AOS
@test isapprox(x_trim[4], 0)  # PHI
@test isapprox(x_trim[5], 0.03936, atol=0.00005)  # THETA
@test isapprox(x_trim[7], 0)  # P
@test isapprox(x_trim[8], 0)  # Q
@test isapprox(x_trim[9], 0)  # R
@test isapprox(controls_trim[1], 0.1485, atol=0.00005)  # THTL
@test isapprox(controls_trim[2], -1.931, atol=0.0001)  # DE
@test isapprox(controls_trim[3], -7e-8, atol=1e-6)  # DA
@test isapprox(controls_trim[4], 8.3e-7, atol=1e-6)  # DR

# XCG = 0.38 (third column)
xcg = 0.38

x = [
    502*FT2M,  #-> vt (m/s)
    deg2rad(10.),  # -> alpha (rad)
    0.0,  # -> beta (rad)
    0.0,  # -> phi (rad)
    deg2rad(10.),  # -> theta (rad)
    0.0,  #  -> psi (rad)
    0.0,  # -> P (rad/s)
    0.0,  # -> Q (rad/s)
    0.0,  # -> R (rad/s)
    0.0,  # -> North (m)
    0.0,  # -> East (m)
    0.0,  # -> Altitude (m)
    50.0,  # -> Pow
]

controls = [
    0.5,  # thtl
    0.0,  # elev
    0.0,  # ail
    0.0,  # rudder
]

x_trim, controls_trim, x_dot_trim, outputs_trim, cost = F16.trimmer(
        f,
        x,
        controls,
        F16.F16Stevens(F16.MASS, F16.INERTIA, xcg),
        F16StevensAtmosphere,
        LHDownGravity(FlightMechanicsSimulator.F16.GD*FT2M),
        0.0,
        0.0,
    )

@test isapprox(cost, zeros(6), atol=1e-12)
@test isapprox(x_trim[2], 0.03544, atol=0.00005)  # AOA
@test isapprox(x_trim[3], 3.1e-8, atol=1e-7)  # AOS
@test isapprox(x_trim[4], 0)  # PHI
@test isapprox(x_trim[5], 0.03544, atol=0.00005)  # THETA
@test isapprox(x_trim[7], 0)  # P
@test isapprox(x_trim[8], 0)  # Q
@test isapprox(x_trim[9], 0)  # R
@test isapprox(controls_trim[1], 0.1325, atol=0.0001)  # THTL
@test isapprox(controls_trim[2], -0.05590, atol=0.0005)  # DE
@test isapprox(controls_trim[3], -5.1e-8, atol=1e-6)  # DA
@test isapprox(controls_trim[4], 4.3e-6, atol=1e-5)  # DR


# Coordinated turn (fourth column)
xcg = 0.3

x = [
    502*FT2M,  #-> vt (m/s)
    deg2rad(10.),  # -> alpha (rad)
    0.0,  # -> beta (rad)
    0.0,  # -> phi (rad)
    deg2rad(10.),  # -> theta (rad)
    0.0,  #  -> psi (rad)
    0.0,  # -> P (rad/s)
    0.0,  # -> Q (rad/s)
    0.0,  # -> R (rad/s)
    0.0,  # -> North (m)
    0.0,  # -> East (m)
    0.0,  # -> Altitude (m)
    50.0,  # -> Pow
]

controls = [
    0.5,  # thtl
    0.0,  # elev
    0.0,  # ail
    0.0,  # rudder
]

x_trim, controls_trim, x_dot_trim, outputs_trim, cost = F16.trimmer(
        f,
        x,
        controls,
        F16.F16Stevens(F16.MASS, F16.INERTIA, xcg),
        F16StevensAtmosphere,
        LHDownGravity(FlightMechanicsSimulator.F16.GD*FT2M),
        0.0,
        0.3,  # rad/s
    )

@test isapprox(cost, zeros(6), atol=1e-12)
@test isapprox(x_trim[2], 0.2485, atol=0.0005)  # AOA
@test isapprox(x_trim[3], 4.8e-4, atol=0.00005)  # AOS
@test isapprox(x_trim[4], 1.367, atol=0.0005)  # PHI
@test isapprox(x_trim[5], 0.05185, atol=0.00005)  # THETA
@test isapprox(x_trim[7], -0.01555, atol=0.00001)  # P
@test isapprox(x_trim[8], 0.2934, atol=0.00005)  # Q
@test isapprox(x_trim[9], 0.06071, atol=0.000005)  # R
@test isapprox(controls_trim[1], 0.8499, atol=0.0005)  # THTL
@test isapprox(controls_trim[2], -6.256, atol=0.001)  # DE
@test isapprox(controls_trim[3], 0.09891, atol=0.00005)  # DA
@test isapprox(controls_trim[4], -0.4218, atol=0.0005)  # DR
