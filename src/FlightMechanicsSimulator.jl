module FlightMechanicsSimulator

using FlightMechanicsUtils


export sixdof_aero_earth_euler_fixed_mass
include("dynamic_systems.jl")

export F16
include("f16/F16.jl")  # F16 submodule

export simulate
include("simulation.jl")

export ConstantInput, StepInput, DoubletInput, RampInput, SinusoidalInput, get_value
include("models/inputs.jl")

export ISA1976, F16StevensAtmosphere
export get_density, get_height, get_pressure, get_sound_velocity, get_temperature
include("models/atmosphere.jl")

end
