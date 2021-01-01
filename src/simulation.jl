using OrdinaryDiffEq


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

    tspan = (tini, tfin)
    p = [mass, xcg, controls]

    prob = ODEProblem{false}(F16.f, x0, tspan, p)
    sol = solve(prob, RK4(), reltol=1e-10, saveat=dt)

    results = hcat([[sol.t[ii]; sol.u[ii]] for ii in 1:length(sol.t)]...)'

    return results
end
