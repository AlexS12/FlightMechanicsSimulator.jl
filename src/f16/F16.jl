module F16

using StaticArrays
using FlightMechanicsUtils
using ..FlightMechanicsSimulator

import ..FlightMechanicsSimulator: get_mass, get_inertia_tensor, get_cg_mac
import ..FlightMechanicsSimulator: get_chord, get_surface, get_wing_span

include("aircraft.jl")
include("adc.jl")
include("engine.jl")
include("aero.jl")
include("sixdof.jl")
include("trimmer.jl")

end
