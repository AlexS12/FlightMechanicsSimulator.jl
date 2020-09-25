module FlightMechanicsSimulator

module Fortran
    include("fortran_wrapper.jl")
end

include("constants.jl")
include("adc.jl")
include("engine.jl")
include("aero.jl")
include("sixdof.jl")
include("rk4.jl")
include("trimmer.jl")

end
