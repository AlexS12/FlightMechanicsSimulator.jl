abstract type Aircraft end

# TODO: DOC
# Mass properties
function get_mass(ac::Aircraft) end
function get_inertia_tensor(ac::Aircraft) end
function get_cg_mac(ac::Aircraft) end

# Geometric properties
function get_surface(ac::Aircraft) end
function get_chord(ac::Aircraft) end
function get_wing_span(ac::Aircraft) end
