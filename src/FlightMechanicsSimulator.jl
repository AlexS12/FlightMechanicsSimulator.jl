module FlightMechanicsSimulator

using NLsolve
using StaticArrays

using FlightMechanicsUtils


export sixdof_aero_earth_euler_fixed_mass
include("dynamic_systems.jl")

export trim
include("trimmer.jl")

export simulate
export f
include("simulation.jl")

export ConstantInput, StepInput, DoubletInput, RampInput, SinusoidalInput, get_value
include("models/inputs.jl")

export ISA1976, F16StevensAtmosphere
export get_density, get_height, get_pressure, get_sound_velocity, get_temperature
include("models/atmosphere.jl")

export LHDownGravity
export get_gravity_accel, get_gravity_body, get_gravity_horizon
include("models/gravity.jl")

export Aircraft
export get_mass, get_inertia_tensor, get_cg_mac
export get_chord, get_surface, get_wing_span
export calculate_prop_forces_moments, calculate_prop_gyro_effects, calculate_pdot
export calculate_aero_forces_moments
export tgear
include("models/aircraft.jl")

export F16
include("f16/F16.jl")  # F16 submodule

end
