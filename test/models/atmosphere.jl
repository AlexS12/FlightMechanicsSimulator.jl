using Test

using FlightMechanicsSimulator


@testset "ISA1976" begin
    # Check constructor
    atmosphere = ISA1976(0.0)

    # Check return type is the expected one.
    @test isa(atmosphere, ISA1976)
    # Check inheritance.
    @test isa(atmosphere, FlightMechanicsSimulator.Atmosphere)

    # Check getters
    @test isapprox(get_temperature(atmosphere), 288.15)
    @test isapprox(get_pressure(atmosphere), 101325)
    @test isapprox(get_density(atmosphere), 1.225)
    @test isapprox(get_sound_velocity(atmosphere), 340.293988026089)
    @test isapprox(get_height(atmosphere), 0.0)

    # Check against other constructor
    atmosphere2 = ISA1976(
        atmosphere.temperature,
        atmosphere.pressure,
        atmosphere.density,
        atmosphere.sound_velocity,
        atmosphere.height,
    )
    @test atmosphere2 === atmosphere
end


@testset "F16StevensAtmosphere" begin
    # Check constructor
    atmosphere = F16StevensAtmosphere(0.0)

    # Check return type is the expected one.
    @test isa(atmosphere, F16StevensAtmosphere)
    # Check inheritance.
    @test isa(atmosphere, FlightMechanicsSimulator.Atmosphere)

    # Check getters
    @test isapprox(get_temperature(atmosphere), 288, atol=1)
    @test isapprox(get_pressure(atmosphere), 101325, atol=100)
    @test isapprox(get_density(atmosphere), 1.225, atol=0.001)
    @test isapprox(get_sound_velocity(atmosphere), 340, atol=0.5)
    @test isapprox(get_height(atmosphere), 0.0)

    # Check against other constructor
    atmosphere2 = F16StevensAtmosphere(
        atmosphere.temperature,
        atmosphere.pressure,
        atmosphere.density,
        atmosphere.sound_velocity,
        atmosphere.height,
    )
    @test atmosphere2 === atmosphere
end
