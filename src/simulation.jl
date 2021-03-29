using OrdinaryDiffEq


"""
    simulate(tini, tfin, dt, x0, mass, xcg, controls)

Propagate a simulation from tini to tfin with dt time step.
- tini: initial time (s)
- tfin: final simulation time (s)
- dt: time step (s)
- x0: initial state. Array{13, Number} according to `F16.f`
- controls: inputs. Array{4, Input} according to `F16.f`
- aircraft: aircraft instance.
"""
function simulate(tini, tfin, dt, x0, controls, aircraft, atmosphere, gravity;
    solver=TSit5(), solve_args=Dict()
    )

    tspan = (tini, tfin)
    p = [controls, aircraft, atmosphere, gravity]

    prob = ODEProblem{false}(F16.f, x0, tspan, p)
    sol = solve(prob, solver; solve_args...)

    results = hcat([[sol.t[ii]; sol.u[ii]] for ii in 1:length(sol.t)]...)'

    return results
end
