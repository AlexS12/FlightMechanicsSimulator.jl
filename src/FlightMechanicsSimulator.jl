module FlightMechanicsSimulator

export gD, R_AIR, γ_AIR,
  RAD2DEG, DEG2RAD,
  M2FT, FT2M,
  KEL2RANK, RANK2KEL,
  KG2LB, LB2KG, SLUG2KG, KG2SLUG,
  PSF2PA, PA2PSF,
  SLUGFT32KGM3, KGM32SLUGFT3
include("constants.jl")

export atmosphere_isa
include("atmosphere.jl")

export sixdof_aero_earth_euler_fixed_mass
include("dynamic_systems.jl")

export F16
include("f16/F16.jl")  # F16 submodule

export simulate
include("simulation.jl")

export ConstantInput, StepInput, DoubletInput, RampInput, SinusoidalInput, get_value
include("models/inputs.jl")

end
