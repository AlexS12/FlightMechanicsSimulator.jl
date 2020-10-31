const AXX = 9496.0  # slug·ft²
const AYY = 55814.0  # slug·ft²
const AZZ = 63100.0  # slug·ft²
const AXZ = 982.0  # slug·ft²

const AXZS = AXZ^2
const XPQ = AXZ * (AXX - AYY + AZZ)
const GAM = AXX * AZZ - AXZ^2

const XQR = AZZ * (AZZ - AYY) + AXZS
const ZPQ = (AXX - AYY) * AXX + AXZS

const YPR= AZZ - AXX

const WEIGHT= 20500.0  #
const GD= 32.17  # ft/s^2
const MASS= WEIGHT / GD  # lb

const S = 300  # ft^2
const B = 30  # ft
const CBAR = 11.32  # ft
const XCGR = 0.35  # units MAC
const HX = 160.0  # slug·ft²

const DE_MAX = 25.0  # deg
const DA_MAX = 20.0  # deg  #XXX: In Stevens' book says 21.5 deg (Appendix A Section A.4)
const DR_MAX = 30.0  # deg

const RTOD = 57.29578
const DEG2RAD = 1 / RTOD

const R0 = 2.377e-3  # Sea level density  (slug/ft²)
