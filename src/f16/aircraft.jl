# F16 constants
const MASS = 20500.0 * LB2KG # kgf

const AXX = 9496.0 * SLUG2KG * FT2M^2 # slug·ft² -> kg·m²
const AYY = 55814.0 * SLUG2KG * FT2M^2  # slug·ft² -> kg·m²
const AZZ = 63100.0 * SLUG2KG * FT2M^2 # slug·ft² -> kg·m²
const AXZ = 982.0 * SLUG2KG * FT2M^2 # slug·ft² -> kg·m²

const INERTIA = @SMatrix [
    AXX 0.0 AXZ;
    0.0 AYY 0.0;
    AXZ 0.0 AZZ
]  # Kg·m²

const XCGR = 0.35  # units MAC

const S = 300 * FT2M^2  # m²
const CBAR = 11.32 * FT2M # m
const B = 30 * FT2M  # m

const HX = 160.0 * SLUG2KG * FT2M^2 # slug·ft² -> kg·m²

const DE_MAX = 25.0  # deg
const DA_MAX = 20.0  # deg  #XXX: In Stevens' book says 21.5 deg (Appendix A Section A.4)
const DR_MAX = 30.0  # deg

# Get the GD value used in Stevens to pass tests
# TODO: Note that this is different from gD * M2FT!!
const GD = floor(gD * M2FT, digits=2)


# TODO: doc
struct F16Stevens{T}<:Aircraft
    mass::T
    inertia::SMatrix{3, 3, T}
    cg_mac::T
end


F16Stevens() = F16Stevens(MASS, INERTIA, XCGR)


get_mass(ac::F16Stevens) = ac.mass
get_inertia_tensor(ac::F16Stevens) = ac.inertia
get_xcg_mac(ac::F16Stevens) = ac.cg_mac

get_surface(ac::F16Stevens) = S
get_chord(ac::F16Stevens) = CBAR
get_wing_span(ac::F16Stevens) = B
