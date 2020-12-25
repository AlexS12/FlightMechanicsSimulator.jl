"""
    simulate(tini, tfin, dt, x0, mass, xcg, controls)

Propagate a simulation from tini to tfin with dt time step.
- tini: initial time (s)
- tfin: final simulation time (s)
- dt: time step (s)
- x0: initial state. Array{13, Number} according to `F16.f`
- mass: aircraft total mass (lb). Constant for the simulation.
- xcg: aircraft CG position MAC [0-1]. Constant for the simulation.
- controls: inputs. Array{4, Input} according to `F16.f`
"""
function simulate(tini, tfin, dt, x0, mass, xcg, controls)

    t = tini
    x = x0

    results = []
    # Append initial condition
    push!(results, vcat([t], x))

    while t < tfin + dt / 2.0
        # Simulate next time step
        x = simulate_timestep(F16.f, dt, x, t, mass, xcg, controls)
        # Store results from previous step
        push!(results, vcat([t], x))
        # Prepare next time step
        t += dt
    end

    # Concat results
    results = hcat(results...)'

    return results
end


function simulate_timestep(f, dt, x, t, mass, xcg, controls)
    # Get control values for current timestep
    controls_arr = get_value.(controls, t)
    # Propagate
    x = F16.rk4(f, dt, x, t, mass, xcg, controls_arr)
end
