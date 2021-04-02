
function trim(
    x_guess, controls_guess, aircraft, atmosphere, gravity, γ=0.0, ψ_dot=0.0;
    show_trace=false,
    ftol=1e-16,
    iterations=5000
)

    #  STATE VECTOR
    # C     X(1)  -> vt (m/s)
    # C     X(2)  -> alpha (rad)
    # C     X(3)  -> beta (rad)
    # C     X(4)  -> phi (rad)
    # C     X(5)  -> theta (rad)
    # C     X(6)  -> psi (rad)
    # C     X(7)  -> P (rad/s)
    # C     X(8)  -> Q (rad/s)
    # C     X(9)  -> R (rad/s)
    # C     X(10) -> North (m)
    # C     X(11) -> East (m)
    # C     X(12) -> Altitude (m)
    # C     X(13) -> POW (0-100)

    # THTL = controls[1]
    # EL = controls[2]
    # AIL = controls[3]
    # RDR = controls[4]

    # TRIMMING SOLUTION
    sol_gues = [
        x_guess[2],  # alpha (rad)
        x_guess[3],  # beta (rad)
        controls_guess..., # thtl  (0-1), el (deg), ail (deg), rdr (deg)
    ]

    # CONSTS
    consts = [
        x_guess[1],  # TAS (m/s)
        x_guess[6],  # psi (rad)
        x_guess[10],  # north (m)
        x_guess[11],  # east (m)
        x_guess[12],  # alt (m)
        ψ_dot,  # ψ_dot (rad/s)
        γ,  # γ (rad)
        aircraft,
        atmosphere,
        gravity,
    ]

    f_opt(sol) = trim_cost_function(sol, consts; full_output=false)

    result = nlsolve(
        f_opt, sol_gues;
        ftol=ftol, show_trace=show_trace, iterations=iterations
    )

     if show_trace
            println(result)
     end

    sol = result.zero
    x, controls, xd, outputs, cost = trim_cost_function(sol, consts; full_output=true)
    return x, controls, xd, outputs, cost
end


function trim_cost_function(sol, consts; full_output=false)

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
    x = [tas, α, β, ϕ, θ, ψ, p, q, r, x, y, alt, tgear(aircraft, thtl)]
    # TODO: allow trimmer to use other dynamic systems for trimming the a/c
    dss = SixDOFAeroEuler(x)

    x_dot, outputs = f(time, dss, controls, aircraft, atmosphere, gravity)

    dss_dot = DSStateDot(dss, x_dot)

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
