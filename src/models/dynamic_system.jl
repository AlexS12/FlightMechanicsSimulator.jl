abstract type DynamicSystemState end


get_x(dss::DynamicSystemState) = dss.x
get_n_states(dss::DynamicSystemState) = length(dss.x)
function state_eqs(dss::DynamicSystemState,time, mass, inertia, forces, moments, h) end


abstract type DynamicSystemStateDot end


get_xdot(dssd::DynamicSystemStateDot) = dssd.xdot
get_dynamic_system_state(dssd::DynamicSystemStateDot) = dssd.dss


struct SixDOFAeroEuler{T}<:DynamicSystemState
    x::SVector{12, T}
end

SixDOFAeroEuler(x::AbstractVector) = SixDOFAeroEuler(SVector{12}(x))


get_n_states(dss::SixDOFAeroEuler) = 12


function state_eqs(dss::SixDOFAeroEuler, time, mass, inertia, forces, moments, h)
     xdot = sixdof_aero_earth_euler_fixed_mass(
         time,
         get_x(dss),
         mass,
         inertia,
         forces,
         moments,
         h
    )
    return SixDOFAeroEulerDot(dss, SVector{12}(xdot))
end


struct SixDOFAeroEulerDot{S<:DynamicSystemState, T}<:DynamicSystemStateDot
    dss::S
    # TODO: SVector could probably be inefered from T
    xdot::SVector{12, T}
end

SixDOFAeroEulerDot(dss, x::AbstractVector) = SixDOFAeroEulerDot(dss, SVector{12}(x))
