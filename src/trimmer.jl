
function trim(
    dss_guess, controls_guess, aircraft, atmosphere, gravity, γ=0.0, ψ_dot=0.0;
    show_trace=false,
    ftol=1e-16,
    iterations=5000
)
    # THTL = controls[1]
    # EL = controls[2]
    # AIL = controls[3]
    # RDR = controls[4]

    # TRIMMING SOLUTION
    sol_gues = [
        get_α(dss_guess),  # alpha (rad)
        get_β(dss_guess),  # beta (rad)
        controls_guess..., # thtl  (0-1), el (deg), ail (deg), rdr (deg)
    ]

    # CONSTS
    consts = [
        get_tas(dss_guess),  # TAS (m/s)
        get_euler_angles(dss_guess)[1],  # psi (rad)
        get_earth_position(dss_guess)...,  # north, east, down (m)
        ψ_dot,  # ψ_dot (rad/s)
        γ,  # γ (rad)
        aircraft,
        atmosphere,
        gravity,
    ]

    ds_type = typeof(dss_guess)

    f_opt(sol) = trim_cost_function(ds_type, sol, consts; full_output=false)

    result = nlsolve(
        f_opt, sol_gues;
        ftol=ftol, show_trace=show_trace, iterations=iterations
    )

     if show_trace
            println(result)
     end

    # TODO: return a DSStateDot instead of x
    sol = result.zero
    x, controls, xd, outputs, cost = trim_cost_function(
        ds_type, sol, consts; full_output=true
    )
    return x, controls, xd, outputs, cost
end


function trim_cost_function(::Type{T}, sol, consts; full_output=false) where {T<:SixDOFAeroEuler}

    tas = consts[1]
    ψ = consts[2]
    x = consts[3]
    y = consts[4]
    alt = consts[5]
    ψ_dot = consts[6]
    γ = consts[7]
    aircraft = consts[8]
    atmosphere = consts[9]
    gravity = consts[10]

    α = sol[1]
    β = sol[2]
    thtl = sol[3]
    controls = sol[3:6]
    gd = get_gravity_accel(gravity)

    ϕ, θ, p, q, r = apply_trimmer_constrains(tas, α, β, γ, ψ_dot, gd)
    # Construct state vector
    # TODO: implement a method to construct x for each ds?
    x = [tas, α, β, ϕ, θ, ψ, p, q, r, x, y, alt, tgear(aircraft, thtl)]

    dss = SixDOFAeroEuler(x)
    x_dot, outputs = f(time, dss, controls, aircraft, atmosphere, gravity)
    dssd = DSStateDot(dss, x_dot)

    cost = [x_dot[1:3]..., x_dot[7:9]...]

    if full_output
        return x, controls, x_dot, outputs, cost
    else
        return cost
    end
end


function apply_trimmer_constrains(tas, α, β, γ, ψ_dot, gd)
    # Coordinated turn bank --> phi
    ϕ = coordinated_turn_bank(ψ_dot, α, β, tas, γ, gd)
    # Climb -> theta
    θ = rate_of_climb_constrain_no_wind(γ, α, β, ϕ)
    # Angular kinemtic -> p, q, r
    p, q, r = ψθϕ_dot_2_pqr(ψ_dot, 0, 0, θ, ϕ)

    return ϕ, θ, p, q, r
end
