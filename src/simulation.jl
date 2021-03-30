"""
    simulate(tini, tfin, dt, x0, mass, xcg, controls)

Propagate a simulation from tini to tfin with dt time step.
- tini: initial time (s)
- tfin: final simulation time (s)
- dt: time step (s)
- x0: initial state. Array{13, Number} according to `F16Stevens.f`
- controls: inputs. Array{4, Input} according to `F16Stevens.f`
- aircraft: aircraft instance.
"""
function simulate(tini, tfin, dt, x0, controls, aircraft, atmosphere, gravity;
    solver=TSit5(), solve_args=Dict()
    )

    tspan = (tini, tfin)
    p = [controls, aircraft, atmosphere, gravity]

    prob = ODEProblem{false}(f, x0, tspan, p)
    sol = solve(prob, solver; solve_args...)

    results = hcat([[sol.t[ii]; sol.u[ii]] for ii in 1:length(sol.t)]...)'

    return results
end


function f(x, p, t)
    controls = p[1]
    aircraft = p[2]
    atmosphere = p[3]
    gravity=p[4]

    controls_arr = get_value.(controls, t)

    x_dot, outputs = f(time, x, controls_arr, aircraft, atmosphere, gravity)

    return x_dot
end


function f(time, x, controls, aircraft, atmosphere, gravity)

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

    atmosphere = atmosphere(height)
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

    x_dot = sixdof_aero_earth_euler_fixed_mass(time, x, mass, inertia, forces, moments, h)

    # Engine dynamic model
    thtl = controls[1]
    pdot = calculate_pdot(ac, thtl, pow)

    x_dot = [x_dot..., pdot]

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
