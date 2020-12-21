function adc(vt, alt)

    tfac = 1.0 - 0.703e-5 * alt

    alt >= 35000.0 ? t = 390.0 : t = 519.0 * tfac  # temperature

    rho = R0 * tfac^4.14  # density
    amach = vt / sqrt(1.4 * 1716.3 * t)  # mach number
    qbar = 0.5 * rho * vt^2  # dynamic pressure
    # ps = 1715.0 * rho * t  # static pressure

    return amach, qbar

end
