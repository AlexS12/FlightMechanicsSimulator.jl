# Constants

# Down component of gravity acceleration at Earth surface at 45º geodetic latitude.
# Stevens, B. L., Lewis, F. L., & Johnson, E. N. (2015). Aircraft control
#  and simulation: dynamics, controls design, and autonomous systems. John Wiley
#  & Sons. Equation (page 33)
const gD = 9.80665  # m/s²

const R0 = 2.377e-3  # Sea level density  (slug/ft³)

# Unit conversions
const RAD2DEG = 57.29578  # Radians (rad) to degrees (deg)
const DEG2RAD = 1 / RAD2DEG

const FT2M = 0.3048  # Feet (ft) to meter (m)
const M2FT = 1 / FT2M
