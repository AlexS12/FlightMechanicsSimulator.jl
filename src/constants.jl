# Constants

# Down component of gravity acceleration at Earth surface at 45º geodetic latitude.
# Stevens, B. L., Lewis, F. L., & Johnson, E. N. (2015). Aircraft control
#  and simulation: dynamics, controls design, and autonomous systems. John Wiley
#  & Sons. Equation (page 33)
const gD = 9.80665  # m/s²

# Air Constants
const γ_AIR = 1.4  # Adiabatic index or ratio of specific heats (dry air at 20º C)
const R_AIR = 287.05287  # Specific gas constant for dry air (J/(Kg·K))

# Unit conversions
const RAD2DEG = 57.29578  # Radians (rad) to degrees (deg)
const DEG2RAD = 1 / RAD2DEG

const FT2M = 0.3048  # Feet (ft) to meter (m)
const M2FT = 1 / FT2M

const KEL2RANK = 1.8  # Kelvin to Rankine
const RANK2KEL = 1 / KEL2RANK

const KG2LB = 2.20462262185  # Kilogram (Kg) to pound (lb)
const LB2KG = 1 / KG2LB

const SLUG2KG = 14.5939029372  # Slug to Kilogram (kg)
const KG2SLUG = 1 / SLUG2KG

const PSF2PA = LB2KG * gD / FT2M^2  # Pounds per square foot (PSF) to Pascal (Pa)
const PA2PSF = 1 / PSF2PA

const SLUGFT32KGM3 = SLUG2KG / FT2M^3  # slug/ft³ to Kg/m³
const KGM32SLUGFT3 = 1 / SLUGFT32KGM3
