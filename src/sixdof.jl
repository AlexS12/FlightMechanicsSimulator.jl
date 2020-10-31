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
    ALPHA = X[2] * RTOD
    BETA = X[3] * RTOD
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
    fx = -MASS * GD * STH + (QS * CXT + T)
    fy = MASS * GD * CTH * SPH + QS * CYT
    fz = MASS * GD * CTH * CPH + QS * CZT

    mx = QS * B * CLT
    my = QS * CBAR * CMT
    mz = QS * B * CNT

    inertia = [
        AXX 0.0 AXZ;
        0.0 AYY 0.0;
        AXZ 0.0 AZZ
        ]

    forces = [fx, fy, fz]
    moments = [mx, my, mz]
    h = [HX, 0, 0]

    xd = sixdof_aero_earth_euler_fixed_mass(time, X, MASS, inertia, forces, moments, h)

    xd = [xd..., xd_13]
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

    return xd, outputs

end


function sixdof_aero_earth_euler_fixed_mass(time, X, MASS, inertia, forces, moments, h)
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

    xd = Array{Float64}(undef, 12)

    # Assign state & control variables
    VT = X[1]
    ALPHA = X[2] * RTOD
    BETA = X[3] * RTOD
    PHI = X[4]
    THETA = X[5]
    PSI = X[6]
    P = X[7]
    Q = X[8]
    R = X[9]
    ALT = X[12]

    # Unpack forces
    fx, fy, fz = forces
    # Unpack moments
    mx, my, mz = moments
    # Unpack angular momentum contributions
    HX, HY, HZ = h
    # Unpack inertia
    AXX, AYY, AZZ = inertia[1, 1], inertia[2, 2], inertia[3, 3]
    AXY, AXZ = inertia[1, 2], inertia[1, 3]

    # Get ready for state equations
    CBTA = cos(X[3])
    U = VT * cos(X[2]) * CBTA
    V = VT * sin(X[3])
    W = VT * sin(X[2]) * CBTA

    STH = sin(THETA)
    CTH = cos(THETA)
    SPH = sin(PHI)
    CPH = cos(PHI)
    SPSI = sin(PSI)
    CPSI = cos(PSI)

    # Force equations
    UDOT = R * V - Q * W + fx / MASS
    VDOT = P * W - R * U + fy / MASS
    WDOT = Q * U - P * V + fz / MASS
    DUM = (U * U + W * W)

    xd[1] = (U * UDOT + V * VDOT + W * WDOT) / VT
    xd[2] = (U * WDOT - W * UDOT) / DUM
    xd[3] = (VT * VDOT - V * xd[1]) * CBTA / DUM

    # Kinematics
    xd[4] = P + (STH / CTH) * (Q * SPH + R * CPH)
    xd[5] = Q * CPH - R * SPH
    xd[6] = (Q * SPH + R * CPH) / CTH

    # Moments
    PQ = P * Q
    QR = Q * R
    QHX = Q * HX

    xd[7] = (XPQ * PQ - XQR * QR + AZZ * mx + AXZ * (mz + QHX)) / GAM
    xd[8] = (YPR * P * R - AXZ * (P^2 - R^2) + my - R * HX) / AYY
    xd[9] = (ZPQ * PQ - XPQ * QR + AXZ * mx + AXX * (mz + QHX)) / GAM

    # Navigation
    T1 = SPH * CPSI
    T2 = CPH * STH
    T3 = SPH * SPSI
    S1 = CTH * CPSI
    S2 = CTH * SPSI
    S3 = T1 * STH - CPH * SPSI
    S4 = T3 * STH + CPH * CPSI
    S5 = SPH * CTH
    S6 = T2 * CPSI + T3
    S7 = T2 * SPSI - T1
    S8 = CPH * CTH

    xd[10] = U * S1 + V * S3 + W * S6  # North speed
    xd[11] = U * S2 + V * S4 + W * S7  # East speed
    xd[12] = U * STH - V * S5 - W * S8  # Vertical speed

    return xd

end
