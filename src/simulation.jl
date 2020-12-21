function simulate(tini, tfin, dt, x0, mass, xcg, controls)

    t = tini
    x = x0

    results = []

    while t < tfin + dt / 2.0
        # Store results from previous step
        push!(results, vcat([t], x))
        # Get control values for current timestep
        controls_arr = get_value.(controls, t)
        # Propagate
        x = F16.rk4(F16.f, dt, x, t, mass, xcg, controls_arr)
        # Prepare next time step
        t += dt
    end
    # Append last result
    push!(results, vcat([t], x))
    # Concat results
    results = hcat(results...)'

    return results
end
