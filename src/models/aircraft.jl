abstract type Aircraft end


# Mass properties
"""
    get_mass(ac::Aircraft)

Get aircraft mass (kg).
"""
function get_mass(ac::Aircraft) end


"""
    get_inertia_tensor(ac::Aircraft)

Get aircraft inertia tensor (kg·m²).
"""
function get_inertia_tensor(ac::Aircraft) end


"""
    get_xcg_mac(ac::Aircraft)

Get aircraft center of gravity position (units of mean aerodynamic chord).
"""
function get_xcg_mac(ac::Aircraft) end

# Geometric properties
"""
    get_surface(ac:Aircraft)

Get aircraft wing surface, ``S`` (m²).
"""
function get_surface(ac::Aircraft) end


"""
    get_chord(ac::Aircraft)

Get aircraft mean aerodynamic chord, ``c``, (m).
"""
function get_chord(ac::Aircraft) end


"""
    get_wing_span(ac::Aircraft)

Get aircraft wing span, ``b``, (m).
"""
function get_wing_span(ac::Aircraft) end


"""
    get_aspect_ratio(ac::Aircraft)

Get aspect ratio, ``b^2/S``.
"""
function get_aspect_ratio(ac::Aircraft)
    b = get_wing_span(ac)
    S = get_surface(ac)
    return b^2 / S
end

# Forces and moments
# TODO: document once headers are standarized. Right now they are too specific to F16Stevens
function calculate_prop_forces_moments end
function calculate_prop_gyro_effects end
function calculate_pdot end

function calculate_aero_forces_moments end

# Flight Control system
function tgear end
