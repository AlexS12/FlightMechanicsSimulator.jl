using SafeTestsets

@safetestset "Atmosphere" begin include("models/atmosphere.jl") end
@safetestset "Gravity" begin include("models/gravity.jl") end
@safetestset "Inputs" begin include("models/inputs.jl") end
@safetestset "F16 ADC" begin include("adc.jl") end
@safetestset "F16 Aerodynamics" begin include("aero.jl") end
@safetestset "F16 Engine" begin include("engine.jl") end
@safetestset "F16 SixDOF" begin include("sixdof.jl") end
@safetestset "F16 Trimmer" begin include("trimmer.jl") end
@safetestset "F16 Simulation" begin include("simulation.jl") end
