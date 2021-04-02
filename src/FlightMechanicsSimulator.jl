module FlightMechanicsSimulator

using DataFrames
using NLsolve
using OrdinaryDiffEq
using StaticArrays

using FlightMechanicsUtils

import Base: eltype


export trim
include("trimmer.jl")

export ConstantInput, StepInput, DoubletInput, RampInput, SinusoidalInput, get_value
include("models/inputs.jl")

export ISA1976, F16StevensAtmosphere
export get_density, get_height, get_pressure, get_sound_velocity, get_temperature
include("models/atmosphere.jl")

export LHDownGravity
export get_gravity_accel, get_gravity_body, get_gravity_horizon
include("models/gravity.jl")

export Aircraft
export get_mass, get_inertia_tensor, get_xcg_mac
export get_chord, get_surface, get_wing_span, get_aspect_ratio
export calculate_prop_forces_moments, calculate_prop_gyro_effects, calculate_pdot
export calculate_aero_forces_moments
export tgear
include("models/aircraft.jl")

export DSStateDot
export get_x, get_n_states, state_eqs
export get_xdot, get_ds_state
include("models/dynamic_system.jl")

export SixDOFAeroEuler
export sixdof_aero_earth_euler_fixed_mass
include("dynamics/sixdof_aero_euler.jl")

export simulate
export f
include("simulation.jl")

export F16
export F16Stevens
include("F16Stevens/F16Stevens.jl")  # F16 submodule
using .F16Stevens

end
