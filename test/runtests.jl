using SafeTestsets

@safetestset "ADC" begin include("adc.jl") end
@safetestset "Aerodynamics" begin include("aero.jl") end
@safetestset "Engine" begin include("engine.jl") end
@safetestset "SixDOF" begin include("sixdof.jl") end
@safetestset "Propagator" begin include("rk4.jl") end
@safetestset "Trimmer" begin include("trimmer.jl") end
@safetestset "Inputs" begin include("models/inputs.jl") end
@safetestset "Atmosphere" begin include("models/atmosphere.jl") end
@safetestset "Gravity" begin include("models/gravity.jl") end
