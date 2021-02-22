# TODO: missing layers (do not forget tests)
# TODO: documentation
function atmosphere_isa(alt)

    α = -0.0065  # K/m
    T0 = 288.15  # K
    p0 = 101325.0  # Pa

    T = T0 + α * alt
    p = p0 * (T0 / T) ^ (gD / (R_AIR * α))

    rho = p / (R_AIR * T)
    a   = sqrt(γ_AIR * R_AIR * T)

    return T, p, rho, a
end
