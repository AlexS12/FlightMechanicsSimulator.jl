module FlightMechanicsSimulator

using DataFrames
using NLsolve
using OrdinaryDiffEq
using StaticArrays

using FlightMechanicsUtils

import Base: eltype


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

export DSState, DSStateDot
export get_x, get_n_states, get_x_names, state_eqs
export get_xdot, get_ds_state
export get_earth_position, get_height
export get_euler_angles
export get_tas, get_α, get_β, get_tasαβ
export get_body_velocity, get_horizon_velocity
export get_ang_vel_body, get_euler_angles_rates
export get_engine_power
export get_tas_dot, get_α_dot, get_β_dot, get_tasαβ_dot
export get_accel_body
export get_ang_accel_body
export get_engine_power_dot
include("models/dynamic_system.jl")

export SixDOFAeroEuler
export sixdof_aero_earth_euler_fixed_mass
include("dynamics/sixdof_aero_euler.jl")

export SixDOFBodyEuler
export sixdof_body_earth_euler_fixed_mass
include("dynamics/sixdof_body_euler.jl")

export trim
include("trimmer.jl")

export simulate
export f
include("simulation.jl")

export F16
export F16Stevens
include("F16Stevens/F16Stevens.jl")  # F16 submodule
using .F16Stevens

end
