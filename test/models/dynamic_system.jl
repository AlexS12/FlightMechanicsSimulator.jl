using Test

using FlightMechanicsUtils

using FlightMechanicsSimulator


x_sixdofaeroeuler = SixDOFAeroEuler([
    502.0 * FT2M,
    0.2392628,
    0.0005061803,
    1.366289,
    0.05000808,
    0.2340769,
    -0.01499617,
    0.2933811,
    0.06084932,
    0.0 * FT2M,
    0.0 * FT2M,
    0.0 * FT2M,
    64.12363,
])

x_sixdofbodyeuler = SixDOFBodyEuler([
    wind2body(502.0 * FT2M, 0, 0,  0.2392628, 0.0005061803)...,
    1.366289,
    0.05000808,
    0.2340769,
    -0.01499617,
    0.2933811,
    0.06084932,
    0.0 * FT2M,
    0.0 * FT2M,
    0.0 * FT2M,
    64.12363,
])

@testset "Interface $(typeof(dss))" for dss in [x_sixdofaeroeuler, x_sixdofbodyeuler]

    methods_scalar_output = [
        get_height,
        get_tas,
        get_α,
        get_β,
        get_engine_power,
    ]

    @testset "$im" for im in methods_scalar_output
        @test hasmethod(im, (typeof(dss),))
        @test isa(im(dss), Number)
    end

    methods_vec3_output = [
        get_earth_position,
        get_tasαβ,
        get_euler_angles,
        get_body_velocity,
        get_horizon_velocity,
        get_ang_vel_body,
        get_euler_angles_rates,
    ]

    @testset "$im" for im in methods_vec3_output
        @test hasmethod(im, (typeof(dss),))
        rv = im(dss)
        @test isa(rv, AbstractArray)
        @test length(rv) == 3
        @test eltype(rv) <: Number
    end

    @testset "$im" for im in [get_n_states]
        @test hasmethod(im, (typeof(dss),))
        rv = im(dss)
        @test isa(rv, Integer)
        @test length(dss.x) == rv
    end

    @testset "$im" for im in [get_x]
        @test hasmethod(im, (typeof(dss),))
        rv = im(dss)
        @test isa(rv, AbstractArray)
        @test length(rv) == get_n_states(dss)
        @test eltype(rv) <: Number
    end

    @testset "$im" for im in [get_x_names]
        @test hasmethod(im, (typeof(dss),))
        rv = im(dss)
        @test isa(rv, AbstractArray)
        @test length(rv) == get_n_states(dss)
        @test eltype(rv) == Symbol
    end

    @testset "state_eqs" begin
        @test hasmethod(
            state_eqs,
            (
                typeof(dss),
                Number,
                Number,
                AbstractMatrix,
                AbstractVector,
                AbstractVector,
                AbstractVector,
                Number,
            ),
        )

        time = 0.0
        mass = 10500.0
        inertia = [555.5 0.0 1234.5; 0.0 5674.0 0.0; 987.2 0.0  7895.3]
        forces = [1.0, 2.0, 3.0]
        moments = [4.0, 5.0, 6.0]
        hx = [12.0, 20.0, 30.0]
        pdot = 12.5

        dssd = state_eqs(dss, time, mass, inertia, forces, moments, hx, pdot)

        @test isa(dssd, DSStateDot)
    end

    @testset "DSStateDot{$(typeof(dss)), $(get_n_states(dss)), $(eltype(dss))}" begin
        # Test constructor
        dssd = DSStateDot(dss, get_x(dss))

        @testset "$im" for im in [get_ds_state]
            @test hasmethod(im, (typeof(dssd),))
            rv = im(dssd)
            @test isa(rv, typeof(dss))
        end

        @testset "$im" for im in [get_xdot]
            @test hasmethod(im, (typeof(dssd),))
            rv = im(dssd)
            @test isa(rv, AbstractArray)
            @test length(rv) == get_n_states(dss)
            @test eltype(rv) <: Number
        end

        @testset "$im" for im in [get_n_states]
            @test hasmethod(im, (typeof(dss),))
            rv = im(dssd)
            @test isa(rv, Integer)
            @test length(dss.x) == rv
        end

        @testset "$im" for im in [get_x]
            @test hasmethod(im, (typeof(dss),))
            rv = im(dssd)
            @test isa(rv, AbstractArray)
            @test length(rv) == get_n_states(dss)
            @test eltype(rv) <: Number
        end

        @testset "$im" for im in [get_x_names]
            @test hasmethod(im, (typeof(dss),))
            rv = im(dssd)
            @test isa(rv, AbstractArray)
            @test length(rv) == get_n_states(dss)
            @test eltype(rv) == Symbol
        end

        methods_scalar_output = [
            get_height,
            get_tas,
            get_α,
            get_β,
            get_engine_power,
            get_tas_dot,
            get_α_dot,
            get_β_dot,
            get_engine_power_dot,
        ]

        @testset "$im" for im in methods_scalar_output
            @test hasmethod(im, (typeof(dssd),))
            @test isa(im(dssd), Number)
        end

        methods_vec3_output = [
            get_earth_position,
            get_tasαβ,
            get_euler_angles,
            get_body_velocity,
            get_horizon_velocity,
            get_ang_vel_body,
            get_euler_angles_rates,
            get_tasαβ_dot,
            get_accel_body,
            get_ang_vel_body,
        ]

        @testset "$im" for im in methods_vec3_output
            @test hasmethod(im, (typeof(dssd),))
            rv = im(dssd)
            @test isa(rv, AbstractArray)
            @test length(rv) == 3
            @test eltype(rv) <: Number
        end
    end
end
