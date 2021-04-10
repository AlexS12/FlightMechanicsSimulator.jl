# TODO: DOC
struct SixDOFBodyEuler{T}<:DSState{T}
    x::SVector{13, T}
end


SixDOFBodyEuler(x::AbstractVector) = SixDOFBodyEuler(SVector{13, eltype(x)}(x))

SixDOFBodyEuler(dss::DSState) = SixDOFBodyEuler(
    [
        get_body_velocity(dss)...,
        get_euler_angles(dss)[3:-1:1]...,
        get_ang_vel_body(dss)...,
        get_earth_position(dss)...,
        get_engine_power(dss),
    ]
)


get_x_names(dss::SixDOFBodyEuler) = [:u, :v, :w, :ϕ, :θ, :ψ, :p, :q, :r, :x, :y, :z, :pow]

# Mandatory getters for DSState
get_earth_position(dss::SixDOFBodyEuler) = get_x(dss)[10:12]
get_height(dss::SixDOFBodyEuler) = -get_x(dss)[12]

get_euler_angles(dss::SixDOFBodyEuler) = get_x(dss)[6:-1:4]

get_tas(dss::SixDOFBodyEuler) = get_tasαβ(dss)[1]
get_α(dss::SixDOFBodyEuler) = get_tasαβ(dss)[2]
get_β(dss::SixDOFBodyEuler) = get_tasαβ(dss)[3]
get_tasαβ(dss::SixDOFBodyEuler) = uvw_to_tasαβ(get_body_velocity(dss)...)
get_body_velocity(dss::SixDOFBodyEuler) = get_x(dss)[1:3]
get_horizon_velocity(dss::SixDOFBodyEuler) = body2horizon(
    get_body_velocity(dss)..., get_euler_angles(dss)...
)

get_ang_vel_body(dss::SixDOFBodyEuler) = get_x(dss)[7:9]
get_euler_angles_rates(dss::SixDOFBodyEuler) = pqr_2_ψθϕ_dot(
    get_ang_vel_body(dss)..., get_euler_angles(dss)[2:3]...
)

get_engine_power(dss::SixDOFBodyEuler) = get_x(dss)[13]

# Mandatory getters for DSStateDot
get_tas_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFBodyEuler, N, T} = get_tasαβ_dot(dssd)[1]
get_α_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFBodyEuler, N, T} = get_tasαβ_dot(dssd)[2]
get_β_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFBodyEuler, N, T} = get_tasαβ_dot(dssd)[3]
get_tasαβ_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFBodyEuler, N, T} = uvw_dot_to_tasαβ_dot(
    get_x(dssd)[1:3]..., get_xdot(dssd)[1:3]...
)
get_accel_body(dssd::DSStateDot{S, N, T}) where {S<:SixDOFBodyEuler, N, T} = get_xdot(dssd)[1:3]
get_ang_accel_body(dssd::DSStateDot{S, N, T}) where {S<:SixDOFBodyEuler, N, T} = get_xdot(dssd)[7:9]
get_engine_power_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFBodyEuler, N, T} = get_xdot(dssd)[13]


function state_eqs(dss::SixDOFBodyEuler, time, mass, inertia, forces, moments, h, pow_dot)

    xdot = sixdof_body_earth_euler_fixed_mass(
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


function sixdof_body_earth_euler_fixed_mass(time, x, mass, inertia, forces, moments, h)

    m = mass
    u, v, w, ϕ, θ, ψ, p, q, r, xe, ye, ze = x

    Fx, Fy, Fz = forces
    L, M, N = moments

    sψ, cψ = sin(ψ), cos(ψ)
    sθ, cθ = sin(θ), cos(θ)
    sϕ, cϕ = sin(ϕ), cos(ϕ)

    # Linear momentum equations
    u_dot = Fx / m + r * v - q * w
    v_dot = Fy / m - r * u + p * w
    w_dot = Fz / m + q * u - p * v

    # Angular momentum equations
    Ix = inertia[1, 1]
    Iy = inertia[2, 2]
    Iz = inertia[3, 3]
    Jxz = -inertia[1, 3]

    Jxz2 = Jxz*Jxz
    Γ = (Ix*Iz - Jxz2)
    temp = (Ix + Iz - Iy)

    # Engine angular momentum contribution
    hx, hy, hz = h

    rhy_qhz = (r*hy - q*hz)
    qhx_phy = (q*hx - p*hy)

    pe_dot = Iz * rhy_qhz + Jxz * qhx_phy
    qe_dot = -r*hx + p*hz
    re_dot = Jxz * rhy_qhz + Ix * qhx_phy

    # Angular momentum equations
    p_dot = L*Iz + N*Jxz - q*r*(Iz*Iz - Iz*Iy + Jxz2) + p*q * Jxz * temp + pe_dot
    p_dot /= Γ
    q_dot = M + (Iz - Ix) * p*r - Jxz * (p*p - r*r) + qe_dot
    q_dot /= Iy
    r_dot = L*Jxz + N*Ix + p*q * (Ix*Ix - Ix*Iy + Jxz2) - q*r * Jxz * temp + re_dot
    r_dot /= Γ

    # Angular Kinematic equations
    ψ_dot = (q * sϕ + r * cϕ) / cθ
    θ_dot = q * cϕ - r * sϕ
    # ϕ_dot = p + (q * sϕ + r * cϕ) * tan(θ)
    ϕ_dot = p + ψ_dot * sθ

    # Linear kinematic equations
    xe_dot =  cθ*cψ * u + (sϕ*sθ*cψ - cϕ*sψ) * v + (cϕ*sθ*cψ + sϕ*sψ) * w
    ye_dot =  cθ*sψ * u + (sϕ*sθ*sψ + cϕ*cψ) * v + (cϕ*sθ*sψ - sϕ*cψ) * w
    ze_dot = -sθ    * u +  sϕ*cθ             * v +  cϕ*cθ             * w

    x_dot = [u_dot, v_dot, w_dot, ϕ_dot, θ_dot, ψ_dot, p_dot, q_dot, r_dot, xe_dot, ye_dot, ze_dot]

    return x_dot
end

# TODO: mutating version
