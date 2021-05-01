using Test

using FlightMechanicsUtils

using FlightMechanicsSimulator


@testset "LHDownGravity" begin
    # Check constructor with default value
    grav_1 = LHDownGravity()
    # Check constructor with custom value
    grav_2 = LHDownGravity(F16Stevens.GD * FT2M)

    # Check return type
    @test isa(grav_1, LHDownGravity)
    # Check inheritance
    @test isa(grav_1, FlightMechanicsSimulator.Gravity)

    # Check getters
    @test isapprox(get_gravity_accel(grav_1), 9.80665)
    @test isapprox(get_gravity_accel(grav_2), 9.805416)

    @test isapprox(get_gravity_horizon(grav_1), [0, 0, 9.80665])
    @test isapprox(get_gravity_horizon(grav_2), [0, 0, 9.805416])
end
