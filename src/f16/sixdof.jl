function f(time, X, XCG, controls)

    # C     X(1)  -> vt (ft/s)
    # C     X(2)  -> alpha (rad)
    # C     X(3)  -> beta (rad)
    # C     X(4)  -> phi (rad)
    # C     X(5)  -> theta (rad)
    # C     X(6)  -> psi (rad)
    # C     X(7)  -> P (rad/s)
    # C     X(8)  -> Q (rad/s)
    # C     X(9)  -> R (rad/s)
    # C     X(10) -> North (ft)
    # C     X(11) -> East (ft)
    # C     X(12) -> Altitude (ft)
    # C     X(13) -> Pow

    outputs = Array{Float64}(undef, 7)

    # Assign state & control variables
    VT = X[1]
    ALPHA = X[2] * RAD2DEG
    BETA = X[3] * RAD2DEG
    PHI = X[4]
    THETA = X[5]
    PSI = X[6]
    P = X[7]
    Q = X[8]
    R = X[9]
    ALT = X[12]
    POW = X[13]

    # Air data computer
    AMACH, QBAR = adc(VT, ALT)
    # Engine model
    THTL = controls[1]
    CPOW = tgear(THTL)
    xd_13 = pdot(POW, CPOW)

    # Calculate forces and moments
    T, TY, TZ, MTX, MTY, MTZ = calculate_prop_forces_moments(X, controls)
    CXT, CYT, CZT, CLT, CMT, CNT = calculate_aero_forces_moments(X, controls, XCG)

    STH = sin(THETA)
    CTH = cos(THETA)
    SPH = sin(PHI)
    CPH = cos(PHI)

    QS = QBAR * S

    # Total forces & moments
    Fx = -MASS * GD * STH + (QS * CXT + T)
    Fy = MASS * GD * CTH * SPH + QS * CYT
    Fz = MASS * GD * CTH * CPH + QS * CZT

    L = QS * B * CLT
    M = QS * CBAR * CMT
    N = QS * B * CNT

    inertia = [
        AXX 0.0 AXZ;
        0.0 AYY 0.0;
        AXZ 0.0 AZZ
        ]

    forces = [Fx, Fy, Fz]
    moments = [L, M, N]
    h = [HX, 0, 0]

    x_dot = sixdof_aero_earth_euler_fixed_mass(time, X, MASS, inertia, forces, moments, h)

    x_dot = [x_dot..., xd_13]
    # Outputs
    RMQS = QS / MASS

    AX = (QS * CXT + T) / GD  # <<-- ASM: Definition missing
    AY = RMQS * CYT
    AZ = RMQS * CZT

    AN = -AZ / GD
    ALAT = AY / GD

    outputs[1] = AN
    outputs[2] = ALAT
    outputs[3] = AX
    outputs[4] = QBAR
    outputs[5] = AMACH
    outputs[6] = Q
    outputs[7] = ALPHA

    return x_dot, outputs

end
