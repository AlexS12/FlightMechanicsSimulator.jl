function simulate(tini, tfin, dt, x0, mass, xcg, controls)

    t = tini
    x = x0

    results = []

    while t < tfin + dt / 2.0
        # Store results from previous step
        push!(results, vcat([t], x))
        # Propagate
        x = F16.rk4(F16.f, dt, x, t, mass, xcg, controls)
        # Prepare next time step
        t += dt
    end
    # Append last result
    push!(results, vcat([t], x))
    # Concat results
    results = hcat(results...)'

    return results
end
