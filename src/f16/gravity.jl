function calculate_gravity_forces(g, mass, θ, ϕ)
    fxb, fyb, fzb = horizon2body(0, 0, mass * g, 0, θ, ϕ)
    return fxb, fyb, fzb
end
