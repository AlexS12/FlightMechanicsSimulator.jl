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
    psi = consts[4]
    x = consts[5]
    y = consts[6]
    alt = consts[7]
    turn_rate = consts[8]
    gamma = consts[9]

    alpha = sol[1]
    beta = sol[2]
    thtl = sol[3]
    controls = sol[3:6]

    x = calculate_state_with_constrains(tas, alpha, beta, gamma, turn_rate, x, y, alt, psi, thtl)

    x_dot, outputs = fun(time, x, mass, xcg, controls)

    cost = [x_dot[1:3]..., x_dot[7:9]...]

    if full_output
        return x, controls, x_dot, outputs, cost
    else
        return cost
    end
end


function calculate_state_with_constrains(tas, alpha, beta, gamma, turn_rate, x, y, alt, psi, thtl)
    # Coordinated turn bank --> phi
    # TODO: should use gD and not GD*FT2M. But tests against Stevens would fail
    G = turn_rate * tas / (GD * FT2M)

    if abs(gamma) < 1e-8
        phi = G * cos(beta) / (cos(alpha) - G * sin(alpha) * sin(beta))
        phi = atan(phi)
    else
        a = 1 - G * tan(alpha) * sin(beta)
        b = sin(gamma) / cos(beta)
        c = 1 + G^2 * cos(beta)^2

        sq = sqrt(c * (1 - b^2) + G^2 * sin(beta)^2)

        num = (a - b^2) + b * tan(alpha) * sq
        den = a ^ 2 - b^2 * (1 + c * tan(alpha)^2)

        phi = atan(G * cos(beta) / cos(alpha) * num / den)
    end

    # Climb -> theta
    a = cos(alpha) * cos(beta)
    b = sin(phi) * sin(beta) + cos(phi) * sin(alpha) * cos(beta)
    sq = sqrt(a^2 - sin(gamma)^2 + b^2)
    theta = (a * b + sin(gamma) * sq) / (a^2 - sin(gamma)^2)
    theta = atan(theta)

    # Angular kinemtic -> p, q, r
    p = - turn_rate * sin(theta)
    q = turn_rate * sin(phi) * cos(theta)
    r = turn_rate * cos(theta) * cos(phi)

    x = [
        tas,
        alpha,
        beta,
        phi,
        theta,
        psi,
        p,
        q,
        r,
        x,
        y,
        alt,
        tgear(thtl)
    ]
end
