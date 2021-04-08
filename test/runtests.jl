using SafeTestsets

@safetestset "Atmosphere" begin include("models/atmosphere.jl") end
@safetestset "Gravity" begin include("models/gravity.jl") end
@safetestset "Inputs" begin include("models/inputs.jl") end
@safetestset "Dynamic Systems" begin include("models/dynamic_system.jl") end
@safetestset "F16 Aircraft" begin include("f16/aircraft.jl") end
@safetestset "F16 ADC" begin include("f16/adc.jl") end
@safetestset "F16 Aerodynamics" begin include("f16/aero.jl") end
@safetestset "F16 Engine" begin include("f16/engine.jl") end
@safetestset "F16 SixDOF" begin include("f16/sixdof.jl") end
@safetestset "F16 Trimmer" begin include("f16/trimmer.jl") end
@safetestset "F16 Simulation" begin include("f16/simulation.jl") end
