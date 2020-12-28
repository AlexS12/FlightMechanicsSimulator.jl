function calculate_gravity_forces(g, mass, θ, ϕ)
    fxb = -mass * g * sin(θ)
    fyb = mass * g * cos(θ) * sin(ϕ)
    fzb = mass * g * cos(θ) * cos(ϕ)
    return fxb, fyb, fzb
end