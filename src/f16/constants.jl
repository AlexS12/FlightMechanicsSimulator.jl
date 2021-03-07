# F16 constants
const AXX = 9496.0 * SLUG2KG * FT2M^2 # slug·ft² -> kg·m²
const AYY = 55814.0 * SLUG2KG * FT2M^2  # slug·ft² -> kg·m²
const AZZ = 63100.0 * SLUG2KG * FT2M^2 # slug·ft² -> kg·m²
const AXZ = 982.0 * SLUG2KG * FT2M^2 # slug·ft² -> kg·m²

# Get the GD value used in Stevens to pass tests
# TODO: Note that this is different from gD * M2FT!!
const GD = floor(gD * M2FT, digits=2)

const MASS = 20500.0 * LB2KG # kgf

const S = 300 * FT2M^2  # m²
const B = 30 * FT2M  # m
const CBAR = 11.32 * FT2M # m
const XCGR = 0.35  # units MAC
const HX = 160.0 * SLUG2KG * FT2M^2 # slug·ft² -> kg·m²

const DE_MAX = 25.0  # deg
const DA_MAX = 20.0  # deg  #XXX: In Stevens' book says 21.5 deg (Appendix A Section A.4)
const DR_MAX = 30.0  # deg
