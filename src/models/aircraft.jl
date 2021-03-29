abstract type Aircraft end

# TODO: DOC
# Mass properties
function get_mass(ac::Aircraft) end
function get_inertia_tensor(ac::Aircraft) end


"""
    get_xcg_mac(ac::Aircraft)

Get aircraft center of gravity position (units of mean aerodynamic chord).
"""
function get_xcg_mac(ac::Aircraft) end

# Geometric properties
function get_surface(ac::Aircraft) end
function get_chord(ac::Aircraft) end
function get_wing_span(ac::Aircraft) end

# Forces and moments
function calculate_prop_forces_moments end
function calculate_prop_gyro_effects end
function calculate_pdot end

function calculate_aero_forces_moments end

# Flight Control system
function tgear end
