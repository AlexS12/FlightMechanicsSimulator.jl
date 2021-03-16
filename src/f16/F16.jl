module F16

using ..FlightMechanicsSimulator
using FlightMechanicsUtils

include("constants.jl")
include("adc.jl")
include("engine.jl")
include("aero.jl")
include("gravity.jl")
include("sixdof.jl")
include("trimmer.jl")

end
