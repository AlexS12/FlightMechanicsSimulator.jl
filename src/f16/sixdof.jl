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


# TODO: doc
function sixdof_aero_earth_euler_fixed_mass(time, x, mass, inertia, forces, moments, h)
    # C     x(1)  -> tas (ft/s)
    # C     x(2)  -> α (rad)
    # C     x(3)  -> β (rad)
    # C     x(4)  -> ϕ (rad)
    # C     x(5)  -> θ (rad)
    # C     x(6)  -> ψ (rad)
    # C     x(7)  -> p (rad/s)
    # C     x(8)  -> q (rad/s)
    # C     x(9)  -> r (rad/s)
    # C     x(10) -> North (ft)
    # C     x(11) -> East (ft)
    # C     x(12) -> Altitude (ft)

    x_dot = Array{Float64}(undef, 12)

    # Assign state
    tas = x[1]
    α = x[2]
    β = x[3]
    ϕ = x[4]
    θ = x[5]
    ψ = x[6]
    p = x[7]
    q = x[8]
    r = x[9]

    # Unpack forces
    Fx, Fy, Fz = forces
    # Unpack moments
    L, M, N = moments
    # Unpack angular momentum contributions
    hx, hy, hz = h
    # Unpack inertia
    Ixx, Iyy, Izz = inertia[1, 1], inertia[2, 2], inertia[3, 3]
    Ixz = inertia[1, 3]

    # Get ready for state equations
    # TODO: use wind to body
    cβ = cos(β)
    u = tas * cos(α) * cβ
    v = tas * sin(β)
    w = tas * sin(α) * cβ

    sψ, cψ = sin(ψ), cos(ψ)
    sθ, cθ = sin(θ), cos(θ)
    sϕ, cϕ = sin(ϕ), cos(ϕ)

    # Force equations
    udot = r * v - q * w + Fx / mass
    vdot = p * w - r * u + Fy / mass
    wdot = q * u - p * v + Fz / mass
    dum = (u * u + w * w)

    # TODO: use uvw_dot_to_tasαβ_dot
    x_dot[1] = (u * udot + v * vdot + w * wdot) / tas
    x_dot[2] = (u * wdot - w * udot) / dum
    x_dot[3] = (tas * vdot - v * x_dot[1]) * cβ / dum

    # Kinematics
    x_dot[4] = p + (sθ / cθ) * (q * sϕ + r * cϕ)
    x_dot[5] = q * cϕ - r * sϕ
    x_dot[6] = (q * sϕ + r * cϕ) / cθ

    # Moments
    pq = p * q
    qr = q * r

    rhy_qhz = (r * hy - q * hz)
    qhx_phy = (q * hx - p * hy)

    # If inertia is constant this terms are constant too.
    # TODO: think about passing them as arguments to improve speed.
    IxzS = Ixz^2
    xpq = Ixz * (Ixx - Iyy + Izz)
    gam = Ixx * Izz - IxzS
    xqr = Izz * (Izz - Iyy) + IxzS
    zpq = (Ixx - Iyy) * Ixx + IxzS
    ypr = Izz - Ixx

    x_dot[7] = (xpq * pq - xqr * qr + Izz * (L + rhy_qhz) + Ixz * (N + qhx_phy)) / gam
    x_dot[8] = (ypr * p * r - Ixz * (p^2 - r^2) + M - r * hx + p * hz) / Iyy
    x_dot[9] = (zpq * pq - xpq * qr + Ixz * (L + rhy_qhz) + Ixx * (N + qhx_phy)) / gam

    # Navigation
    t1 = sϕ * cψ
    t2 = cϕ * sθ
    t3 = sϕ * sψ
    s1 = cθ * cψ
    s2 = cθ * sψ
    s3 = t1 * sθ - cϕ * sψ
    s4 = t3 * sθ + cϕ * cψ
    s5 = sϕ * cθ
    s6 = t2 * cψ + t3
    s7 = t2 * sψ - t1
    s8 = cϕ * cθ

    x_dot[10] = u * s1 + v * s3 + w * s6  # North speed
    x_dot[11] = u * s2 + v * s4 + w * s7  # East speed
    x_dot[12] = u * sθ - v * s5 - w * s8  # Vertical speed

    return x_dot

end
