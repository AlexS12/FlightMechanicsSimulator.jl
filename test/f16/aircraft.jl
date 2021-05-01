using Test

using StaticArrays

using FlightMechanicsSimulator



@testset "F16 Aircraft" begin
    ac1 = F16()

    @test isa(ac1, FlightMechanicsSimulator.Aircraft)

    @test isapprox(get_mass(ac1), F16Stevens.MASS)
    @test isapprox(get_inertia_tensor(ac1), F16Stevens.INERTIA)
    @test isapprox(get_xcg_mac(ac1), F16Stevens.XCGR)
    @test isapprox(get_wing_span(ac1), F16Stevens.B)
    @test isapprox(get_surface(ac1), F16Stevens.S)
    @test isapprox(get_chord(ac1), F16Stevens.CBAR)
    @test isapprox(get_aspect_ratio(ac1), F16Stevens.B^2 / F16Stevens.S)

    inertia = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    ac2 = F16(12345.0, inertia, 0.5)

    @test isapprox(get_mass(ac2), 12345.0)
    @test isapprox(get_inertia_tensor(ac2), inertia)
    @test isapprox(get_xcg_mac(ac2), 0.5)

end
