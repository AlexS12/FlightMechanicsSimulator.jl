module F16

using StaticArrays
using FlightMechanicsUtils
using ..FlightMechanicsSimulator

import ..FlightMechanicsSimulator: get_mass, get_inertia_tensor, get_cg_mac
import ..FlightMechanicsSimulator: get_chord, get_surface, get_wing_span
import ..FlightMechanicsSimulator: calculate_prop_forces_moments
import ..FlightMechanicsSimulator: calculate_prop_gyro_effects
import ..FlightMechanicsSimulator: calculate_pdot
import ..FlightMechanicsSimulator: calculate_aero_forces_moments
import ..FlightMechanicsSimulator: tgear

include("aircraft.jl")
include("adc.jl")
include("engine.jl")
include("aero.jl")
include("trimmer.jl")

end
