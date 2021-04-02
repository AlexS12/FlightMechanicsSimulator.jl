struct SixDOFAeroEuler{T}<:DSState{T}
    x::SVector{13, T}
end

SixDOFAeroEuler(x::AbstractVector) = SixDOFAeroEuler(SVector{13, eltype(x)}(x))

get_x_names(dss::SixDOFAeroEuler) = [:tas, :α, :β, :ϕ, :θ, :ψ, :p, :q, :r, :x, :y, :z, :pow]


function state_eqs(dss::SixDOFAeroEuler, time, mass, inertia, forces, moments, h, pow_dot)

    xdot = sixdof_aero_earth_euler_fixed_mass(
         time,
         get_x(dss),
         mass,
         inertia,
         forces,
         moments,
         h
    )

    xdot = [xdot..., pow_dot]

    return DSStateDot(dss, xdot)
end


# TODO: doc
function sixdof_aero_earth_euler_fixed_mass!(
    time, x, mass, inertia, forces, moments, h, x_dot
    )
    # C     x(1)  -> tas (m/s)
    # C     x(2)  -> α (rad)
    # C     x(3)  -> β (rad)
    # C     x(4)  -> ϕ (rad)
    # C     x(5)  -> θ (rad)
    # C     x(6)  -> ψ (rad)
    # C     x(7)  -> p (rad/s)
    # C     x(8)  -> q (rad/s)
    # C     x(9)  -> r (rad/s)
    # C     x(10) -> North (m)
    # C     x(11) -> East (m)
    # C     x(12) -> Altitude (m)

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
    u, v, w = wind2body(tas, 0, 0, α, β)

    sψ, cψ = sin(ψ), cos(ψ)
    sθ, cθ = sin(θ), cos(θ)
    sϕ, cϕ = sin(ϕ), cos(ϕ)

    # Force equations
    udot = r * v - q * w + Fx / mass
    vdot = p * w - r * u + Fy / mass
    wdot = q * u - p * v + Fz / mass

    x_dot[1], x_dot[2], x_dot[3] = uvw_dot_to_tasαβ_dot(u, v, w, udot, vdot, wdot)

    # Kinematics
    x_dot[6], x_dot[5], x_dot[4] = pqr_2_ψθϕ_dot(p, q, r, θ, ϕ)

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

end


# TODO: doc
function sixdof_aero_earth_euler_fixed_mass(time, x, mass, inertia, forces, moments, h)
    x_dot = Array{Float64}(undef, 12)
    sixdof_aero_earth_euler_fixed_mass!(time, x, mass, inertia, forces, moments, h, x_dot)
    return x_dot
end
