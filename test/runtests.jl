using SafeTestsets

@safetestset "F16 ADC" begin include("adc.jl") end
@safetestset "F16 Aerodynamics" begin include("aero.jl") end
@safetestset "F16 Engine" begin include("engine.jl") end
@safetestset "SixDOF" begin include("sixdof.jl") end
@safetestset "Simulation" begin include("simulation.jl") end
@safetestset "Trimmer" begin include("trimmer.jl") end
@safetestset "Inputs" begin include("models/inputs.jl") end
@safetestset "Atmosphere" begin include("models/atmosphere.jl") end
@safetestset "Gravity" begin include("models/gravity.jl") end
