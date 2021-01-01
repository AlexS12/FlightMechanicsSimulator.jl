function atmosphere(alt)
    # alt (ft)
    tfac = 1.0 - 0.703e-5 * alt
    # temperature (Rankine)
    alt >= 35000.0 ? T = 390.0 : T = 519.0 * tfac
    # density (slug / ft^3)
    ρ = R0 * tfac^4.14
    # sound velocity (ft/s)
    a = sqrt(1.4 * 1716.3 * T)
    # pressure (psf)
    p = 1715.0 * ρ * T
    return T, ρ, a, p
end


function adc(vt, T, ρ, a, p)

    amach = vt / a  # mach number
    qbar = 0.5 * ρ * vt^2  # dynamic pressure

    return amach, qbar

end
