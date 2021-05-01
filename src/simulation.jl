"""
    simulate(
        tini, tfin, dss, controls, aircraft, atmosphere, gravity;
        solver=Tsit5(), solve_args=Dict()
    )

"""
function simulate(tini, tfin, dss, controls, aircraft, atmosphere, gravity;
    solver=Tsit5(), solve_args=Dict()
    )

    sol, out = _simulate(tini, tfin, dss, controls, aircraft, atmosphere, gravity;
        solver, solve_args
    )

    df = DataFrame(sol, out)

    return df
end


function _simulate(tini, tfin, dss, controls, aircraft, atmosphere, gravity;
    solver=Tsit5(), solve_args=Dict()
    )

    tspan = (tini, tfin)
    p = [dss, controls, aircraft, atmosphere, gravity]

    x0 = get_x(dss)
    prob = ODEProblem{false}(f, x0, tspan, p)
    sol = solve(prob, solver; solve_args...)
    out = get_log(sol)

    return sol, out
end


function f(x, p, t)
    dss = p[1]
    controls = p[2]
    aircraft = p[3]
    atmosphere = p[4]
    gravity = p[5]

    controls_arr = get_value.(controls, t)

    dss = typeof(dss)(x)

    atmosphere = atmosphere(get_height(dss))

    dssd, outputs = f(time, dss, controls_arr, aircraft, atmosphere, gravity)

    return get_xdot(dssd)
end


function f(
    time,
    dss::DSState,
    controls,
    aircraft::Aircraft,
    atmosphere::Atmosphere,
    gravity::Gravity
    )

    ac = aircraft

    mass = get_mass(ac)
    xcg = get_xcg_mac(ac)
    inertia = get_inertia_tensor(ac)

    S = get_surface(ac)
    c = get_chord(ac)
    b = get_wing_span(ac)

    vt = get_tas(dss)
    α = get_α(dss) * RAD2DEG
    β = get_β(dss) * RAD2DEG
    ψ, θ, ϕ = get_euler_angles(dss)
    p, q, r = get_ang_vel_body(dss)
    height = get_height(dss)
    pow = get_engine_power(dss)

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
    Tx, Ty, Tz, LT, MT, NT = calculate_prop_forces_moments(ac, dss, amach, controls)
    h = calculate_prop_gyro_effects(ac)

    # Engine dynamic model
    thtl = controls[1]
    pdot = calculate_pdot(ac, thtl, pow)

    # Aerodynamics
    Fax, Fay, Faz, La, Ma, Na = calculate_aero_forces_moments(ac, dss, controls, xcg, qbar, S, b, c)
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

    dssd = state_eqs(
        dss, time, mass, inertia, forces, moments, h, pdot
    )

    # Outputs
    gravity_down = get_gravity_accel(gravity)
    outputs = calculate_outputs(dss, amach, qbar, S, mass, gravity_down, [Fax, Fay, Faz], [Tx, Ty, Tz])

    # LOGGING
    @log atmosphere
    @log dssd
    @log FT = [Tx, Ty, Tz]
    @log MT = [LT, MT, NT]
    @log HT = h
    @log Fa = [Fax, Fay, Faz]
    @log Ma = [La, Ma, Na]
    @log Fg = [Fgx, Fgy, Fgz]
    @log F = forces
    @log M = moments
    @log controls
    @log amach
    @log qbar

    return dssd, outputs

end


function calculate_outputs(dss::DSState, amach, qbar, S, mass, g, Fa, Fp)

    outputs = Array{Float64}(undef, 7)

    α = get_α(dss) * RAD2DEG
    p, q, r = get_ang_vel_body(dss)


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


function DataFrame(sol::ODESolution, out::SimulationLog)

    # Get variables from SimulationLogs
    out = get_log(sol)

    rv = get_tasαβ.(out.dssd)
    tas, α, β = [getindex.(rv, ii) for ii in 1:3]
    rv = get_body_velocity.(out.dssd)
    u, v, w = [getindex.(rv, ii) for ii in 1:3]
    rv = get_euler_angles.(out.dssd)
    ψ, θ, ϕ = [getindex.(rv, ii) for ii in 1:3]
    rv = get_ang_vel_body.(out.dssd)
    p, q, r = [getindex.(rv, ii) for ii in 1:3]
    rv = get_earth_position.(out.dssd)
    x, y, z = [getindex.(rv, ii) for ii in 1:3]
    h = get_height.(out.dssd)
    pow = get_engine_power.(out.dssd)
    pow_dot = get_engine_power_dot.(out.dssd)
    rv = get_euler_angles_rates.(out.dssd)
    ψ_dot, θ_dot, ϕ_dot = [getindex.(rv, ii) for ii in 1:3]
    rv = get_accel_body.(out.dssd)
    u_dot, v_dot, w_dot = [getindex.(rv, ii) for ii in 1:3]
    rv = get_ang_accel_body.(out.dssd)
    p_dot, q_dot, r_dot = [getindex.(rv, ii) for ii in 1:3]
    rv = get_tasαβ_dot.(out.dssd)
    tas_dot, α_dot, β_dot = [getindex.(rv, ii) for ii in 1:3]

    T_atm = get_temperature.(out.atmosphere)
    ρ_atm = get_density.(out.atmosphere)
    p_atm = get_pressure.(out.atmosphere)
    a_atm = get_sound_velocity.(out.atmosphere)

    mach = out.amach
    qbar = out.qbar

    Fxa, Fya, Fza = [getindex.(out.Fa, ii) for ii in 1:3]
    Fxt, Fyt, Fzt = [getindex.(out.FT, ii) for ii in 1:3]
    Fxg, Fyg, Fzg = [getindex.(out.Fg, ii) for ii in 1:3]
    Fx, Fy, Fz = [getindex.(out.F, ii) for ii in 1:3]
    Ma, La, Na = [getindex.(out.Ma, ii) for ii in 1:3]
    MT, LT, NT = [getindex.(out.MT, ii) for ii in 1:3]
    HTx, HTy, HTz = [getindex.(out.HT, ii) for ii in 1:3]
    M, L, N = [getindex.(out.M, ii) for ii in 1:3]

    thtl, elev, ail, rud = [getindex.(out.controls, ii) for ii in 1:4]

    time = sol.t

    df_log = DataFrame(;
        time,
        tas, α, β,
        u, v, w,
        ψ, θ, ϕ,
        p, q, r,
        ψ_dot, θ_dot, ϕ_dot,
        x, y, z,
        h,
        pow, pow_dot,
        u_dot, v_dot, w_dot,
        p_dot, q_dot, r_dot,
        tas_dot, α_dot, β_dot,
        T_atm, p_atm, ρ_atm, a_atm,
        mach, qbar,
        thtl, elev, ail, rud,
        Fxa, Fya, Fza,
        Fxt, Fyt, Fzt,
        Fxg, Fyg, Fzg,
        Fx, Fy, Fz,
        Ma, La, Na,
        MT, LT, NT,
        HTx, HTy, HTz,
        M, L, N,
    )

    # DataFrame from solution
    df_sol = DataFrame(sol', :auto)

    df = hcat(df_log, df_sol)

    return df
end
