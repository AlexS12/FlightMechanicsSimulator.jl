struct TrimConditions
    # TODO: improve typing
    tas::Number  # m/s
    ψ::Number  # rad
    position::AbstractArray  # [x, y, z] earth (m)
    ψ_dot::Number  # rad/s
    γ::Number  # rad
    aircraft::Aircraft
    atmosphere::Atmosphere
    gravity::Gravity
end


struct TrimSolution
    # TODO: improve typing
    α::Number
    β::Number
    controls::AbstractArray
    θ::Number
    ϕ::Number
    p::Number
    q::Number
    r::Number
end



function trim(
    dss_guess, controls_guess, aircraft, atmosphere, gravity, γ=0.0, ψ_dot=0.0;
    show_trace=false,
    ftol=1e-16,
    iterations=5000
)

    # TRIMMING SOLUTION
    sol_gues = [
        get_α(dss_guess),  # alpha (rad)
        get_β(dss_guess),  # beta (rad)
        controls_guess..., # thtl  (0-1), el (deg), ail (deg), rdr (deg)
    ]

    # CONSTS
    trim_conditions = TrimConditions(
        get_tas(dss_guess),  # TAS (m/s)
        get_euler_angles(dss_guess)[1],  # psi (rad)
        get_earth_position(dss_guess),  # north, east, down (m)
        ψ_dot,  # ψ_dot (rad/s)
        γ,  # γ (rad)
        aircraft,
        atmosphere,
        gravity,
    )

    ds_type = typeof(dss_guess)

    f_opt(sol) = trim_cost_function(ds_type, sol, trim_conditions; full_output=false)

    result = nlsolve(
        f_opt, sol_gues;
        ftol=ftol, show_trace=show_trace, iterations=iterations
    )

     if show_trace
            println(result)
     end

    sol = result.zero
    dssd, controls, outputs, cost = trim_cost_function(
        ds_type, sol, trim_conditions; full_output=true
    )
    return dssd, controls, outputs, cost
end


function trim_cost_function(ds::Type{T}, sol, trim_conditions; full_output=false) where {T<:DSState}

    tas = trim_conditions.tas
    ψ = trim_conditions.ψ
    x, y, z = trim_conditions.position
    ψ_dot = trim_conditions.ψ_dot
    γ = trim_conditions.γ
    aircraft = trim_conditions.aircraft
    atmosphere = trim_conditions.atmosphere
    gravity = trim_conditions.gravity

    gd = get_gravity_accel(gravity)

    ϕ, θ, p, q, r = apply_trimmer_constrains(
        tas,  # m/s
        sol[1],  # α (rad)
        sol[2],  # β (rad)
        γ,  # rad
        ψ_dot,  # rad/s
        gd,  # m/s²
    )

    trim_solution = TrimSolution(
        sol[1],  # α (rad)
        sol[2],  # β (rad)
        sol[3:6],  # controls  -> thtl  (0-1), el (deg), ail (deg), rdr (deg)
        θ,  # rad
        ϕ,  # rad
        p,  # rad/s
        q,  # rad/s
        r,  # rad/s
    )
    # Construct state vector
    dss = ds(trim_conditions, trim_solution, aircraft)
    dssd, outputs = f(time, dss, trim_solution.controls, aircraft, atmosphere, gravity)

    x_dot = get_xdot(dssd)
    cost = [x_dot[1:3]..., x_dot[7:9]...]

    if full_output
        return dssd, trim_solution.controls, outputs, cost
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


function SixDOFAeroEuler{T}(tc::TrimConditions, ts::TrimSolution, ac::Aircraft) where {T}
    x = [
        tc.tas,
        ts.α,
        ts.β,
        ts.ϕ,
        ts.θ,
        tc.ψ,
        ts.p,
        ts.q,
        ts.r,
        tc.position...,
        tgear(ac, ts.controls[1])
    ]
    return SixDOFAeroEuler(x)
end


function SixDOFBodyEuler{T}(tc::TrimConditions, ts::TrimSolution, ac::Aircraft) where {T}
    x = [
        wind2body(tc.tas, 0, 0, ts.α, ts.β)...,
        ts.ϕ,
        ts.θ,
        tc.ψ,
        ts.p,
        ts.q,
        ts.r,
        tc.position...,
        tgear(ac, ts.controls[1])
    ]
    return SixDOFBodyEuler(x)
end
