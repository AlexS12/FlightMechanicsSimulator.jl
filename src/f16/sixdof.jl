function f(time, x, mass, xcg, controls)

    # C     x(1)  -> vt (ft/s)
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
    # C     x(13) -> pow

    outputs = Array{Float64}(undef, 7)

    # Assign state
    vt = x[1]
    α = x[2] * RAD2DEG
    β = x[3] * RAD2DEG
    ϕ = x[4]
    θ = x[5]
    ψ = x[6]
    p = x[7]
    q = x[8]
    r = x[9]
    height = x[12]
    pow = x[13]

    # Air data computer
    # TODO: call atmosphere and pass outputs to adc instead of calling atmosphere inside adc
    amach, qbar = adc(vt, height)
    # Engine model
    thtl = controls[1]
    cpow = tgear(thtl)
    xd_13 = pdot(pow, cpow)

    # Calculate forces and moments
    Tx, Ty, Tz, LT, MT, NT = calculate_prop_forces_moments(x, controls)
    CXA, CYA, CZA, LA, MA, NA = calculate_aero_forces_moments(x, controls, xcg)
    Fgx, Fgy, Fgz = calculate_gravity_forces(GD, mass, θ, ϕ)

    qbarS = qbar * S

    # Total forces & moments
    Fx = Fgx + qbarS * CXA + Tx
    Fy = Fgy + qbarS * CYA
    Fz = Fgz + qbarS * CZA

    L = qbarS * B * LA
    M = qbarS * CBAR * MA
    N = qbarS * B * NA

    inertia = [
        AXX 0.0 AXZ;
        0.0 AYY 0.0;
        AXZ 0.0 AZZ
        ]

    forces = [Fx, Fy, Fz]
    moments = [L, M, N]
    h = [HX, 0, 0]

    x_dot = sixdof_aero_earth_euler_fixed_mass(time, x, mass, inertia, forces, moments, h)

    x_dot = [x_dot..., xd_13]
    # Outputs
    rmqs = qbarS / mass

    ax = (qbarS * CXA + Tx) / GD  # <<-- ASM: Definition missing
    ay = rmqs * CYA
    az = rmqs * CZA

    an = -az / GD
    alat = ay / GD

    outputs[1] = an
    outputs[2] = alat
    outputs[3] = ax
    outputs[4] = qbar
    outputs[5] = amach
    outputs[6] = q
    outputs[7] = α

    return x_dot, outputs

end
