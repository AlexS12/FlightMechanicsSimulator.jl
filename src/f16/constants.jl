# F16 constants
const AXX = 9496.0  # slug·ft²
const AYY = 55814.0  # slug·ft²
const AZZ = 63100.0  # slug·ft²
const AXZ = 982.0  # slug·ft²

# Get the GD value used in Stevens to pass tests
const GD = floor(gD * M2FT, digits=2)

const WEIGHT = 20500.0  # lbf
const MASS = WEIGHT / GD  # lb

const S = 300  # ft^2
const B = 30  # ft
const CBAR = 11.32  # ft
const XCGR = 0.35  # units MAC
const HX = 160.0  # slug·ft²

const DE_MAX = 25.0  # deg
const DA_MAX = 20.0  # deg  #XXX: In Stevens' book says 21.5 deg (Appendix A Section A.4)
const DR_MAX = 30.0  # deg
