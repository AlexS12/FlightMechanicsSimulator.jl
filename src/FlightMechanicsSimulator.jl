module FlightMechanicsSimulator

export GD, R0, RAD2DEG, DEG2RAD
include("constants.jl")

export sixdof_aero_earth_euler_fixed_mass
include("dynamic_systems.jl")

include("f16/F16.jl")  # F16 submodule

end
