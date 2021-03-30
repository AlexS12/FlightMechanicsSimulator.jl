function atmosphere(alt)

    R0 = 2.377e-3  # Sea level density  (slug/ft³)

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


function atmosphere_f16(alt)
    T, ρ, a, p = atmosphere(alt * M2FT)
    return T * RANK2KEL, ρ * SLUGFT32KGM3, a * FT2M, p * PSF2PA
end
