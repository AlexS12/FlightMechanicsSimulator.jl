module FlightMechanicsSimulator

export gD, R0, RAD2DEG, DEG2RAD, M2FT, FT2M
include("constants.jl")

export sixdof_aero_earth_euler_fixed_mass
include("dynamic_systems.jl")

export F16
include("f16/F16.jl")  # F16 submodule

export simulate
include("simulation.jl")

export ConstantInput, StepInput, DoubletInput, RampInput, SinusoidalInput, get_value
include("models/inputs.jl")

end
