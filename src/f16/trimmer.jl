# using Optim
using NLsolve


function trimmer(fun, x_guess, controls_guess, γ=0.0, ψ_dot=0.0, mass=MASS, xcg=0.35; show_trace=false, ftol=1e-16, iterations=5000)

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
    # alpha (rad)
    # beta (rad)
    # thtl  (0-1)
    # el  (deg)
    # ail (deg)
    # rdr (deg)
    sol_gues = [x_guess[2], x_guess[3], controls_guess...]

    # CONSTS
    # MASS
    # XCG
    # TAS (m/s)
    # psi (rad)
    # north (m)
    # east (m)
    # alt (m)
    # ψ_dot (rad/s)
    # γ (rad)
    consts = [mass, xcg, x_guess[1], x_guess[6], x_guess[10], x_guess[11], x_guess[12], ψ_dot, γ]

    f_opt(sol) = trim_cost_function(sol, consts, fun; full_output=false)

    result = nlsolve(f_opt, sol_gues; ftol=ftol, show_trace=show_trace, iterations=iterations)

     if show_trace
            println(result)
     end

    sol = result.zero
    x, controls, xd, outputs, cost = trim_cost_function(sol, consts, fun; full_output=true)
    return x, controls, xd, outputs, cost
end


function trim_cost_function(sol, consts, fun; full_output=false)

    mass = consts[1]
    xcg = consts[2]
    tas = consts[3]
    ψ = consts[4]
    x = consts[5]
    y = consts[6]
    alt = consts[7]
    ψ_dot = consts[8]
    γ = consts[9]

    α = sol[1]
    β = sol[2]
    thtl = sol[3]
    controls = sol[3:6]

    x = calculate_state_with_constrains(tas, α, β, γ, ψ_dot, x, y, alt, ψ, thtl)

    x_dot, outputs = fun(time, x, mass, xcg, controls)

    cost = [x_dot[1:3]..., x_dot[7:9]...]

    if full_output
        return x, controls, x_dot, outputs, cost
    else
        return cost
    end
end


function calculate_state_with_constrains(tas, α, β, γ, ψ_dot, x, y, alt, ψ, thtl)
    # Coordinated turn bank --> phi
    # TODO: should use gD and not GD*FT2M. But tests against Stevens would fail
    ϕ = coordinated_turn_bank(ψ_dot, α, β, tas, γ, GD * FT2M)
    # Climb -> theta
    θ = rate_of_climb_constrain_no_wind(γ, α, β, ϕ)
    # Angular kinemtic -> p, q, r
    p, q, r = ψθϕ_dot_2_pqr(ψ_dot, 0, 0, θ, ϕ)

    x = [tas, α, β, ϕ, θ, ψ, p, q, r, x, y, alt, tgear(thtl)]
end
