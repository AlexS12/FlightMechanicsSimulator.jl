"""
    simulate(
        tini, tfin, dss, controls, aircraft, atmosphere, gravity;
        solver=TSit5(), solve_args=Dict()
    )

"""
function simulate(tini, tfin, dss, controls, aircraft, atmosphere, gravity;
    solver=TSit5(), solve_args=Dict()
    )

    tspan = (tini, tfin)
    p = [dss, controls, aircraft, atmosphere, gravity]

    x0 = get_x(dss)
    prob = ODEProblem{false}(f, x0, tspan, p)
    sol = solve(prob, solver; solve_args...)

    df = DataFrame(sol')
    rename!(df, get_x_names(dss))
    df[!, :time] = sol.t

    return df
end


function f(x, p, t)
    dss = p[1]
    controls = p[2]
    aircraft = p[3]
    atmosphere = p[4]
    gravity = p[5]

    controls_arr = get_value.(controls, t)
    # TODO: receive as argument in p
    dynamic_system_state = typeof(dss)(
        SVector{length(x)}(x)
    )

    atmosphere = atmosphere(x[12])

    x_dot, outputs = f(time, dynamic_system_state, controls_arr, aircraft, atmosphere, gravity)

    return x_dot
end


function f(
    time,
    dynamic_system::DynamicSystemState,
    controls,
    aircraft::Aircraft,
    atmosphere::Atmosphere,
    gravity::Gravity
    )

    # C     x(1)  -> vt (m/s)
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
    # C     x(13) -> pow
    ac = aircraft

    mass = get_mass(ac)
    xcg = get_xcg_mac(ac)
    inertia = get_inertia_tensor(ac)

    S = get_surface(ac)
    c = get_chord(ac)
    b = get_wing_span(ac)

    # Assign state
    x = get_x(dynamic_system)

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

    T = get_temperature(atmosphere)
    ρ = get_density(atmosphere)
    a = get_sound_velocity(atmosphere)
    p = get_pressure(atmosphere)

    # Air data computer
    # TODO: this will be inside AeroState
    amach = vt / a  # mach number
    qbar = 0.5 * ρ * vt^2  # dynamic pressure

    # Calculate forces and moments
    # Propulsion
    Tx, Ty, Tz, LT, MT, NT = calculate_prop_forces_moments(ac, x, amach, controls)
    h = calculate_prop_gyro_effects(ac)

    # Engine dynamic model
    thtl = controls[1]
    pdot = calculate_pdot(ac, thtl, pow)

    # Aerodynamics
    Fax, Fay, Faz, La, Ma, Na = calculate_aero_forces_moments(ac, x, controls, xcg, qbar, S, b, c)
    # Gravity
    Fgx, Fgy, Fgz = get_gravity_body(gravity, θ, ϕ) .* mass

    # Total forces & moments
    Fx = Fgx + Fax + Tx
    Fy = Fgy + Fay + Ty
    Fz = Fgz + Faz + Tz

    L = La
    M = Ma
    N = Na

    forces = [Fx, Fy, Fz]
    moments = [L, M, N]

    dynamic_system_state_dot = state_eqs(
        dynamic_system, time, mass, inertia, forces, moments, h, pdot
    )

    x_dot = get_xdot(dynamic_system_state_dot)

    # Outputs
    gravity_down = get_gravity_accel(gravity)
    outputs = calculate_outputs(x, amach, qbar, S, mass, gravity_down, [Fax, Fay, Faz], [Tx, Ty, Tz])

    return x_dot, outputs

end


function calculate_outputs(x, amach, qbar, S, mass, g, Fa, Fp)

    outputs = Array{Float64}(undef, 7)

    # vt = x[1]
    α = x[2] * RAD2DEG
    # β = x[3] * RAD2DEG
    # ϕ = x[4]
    # θ = x[5]
    # ψ = x[6]
    # p = x[7]
    q = x[8]
    # r = x[9]
    # height = x[12]
    # pow = x[13]

    Fax, Fay, Faz = Fa
    Tx, Ty, Tz = Fp

    qbarS = qbar * S

    rmqs = qbarS / mass

    ax = (Fax + Tx) / mass
    ay = (Fay + Ty) / mass
    az = (Faz + Tz) / mass

    along = ax / g
    an = -az / g
    alat = ay / g

    outputs[1] = an
    outputs[2] = alat
    outputs[3] = along
    outputs[4] = qbar
    outputs[5] = amach
    outputs[6] = q
    outputs[7] = α

    return outputs
end
