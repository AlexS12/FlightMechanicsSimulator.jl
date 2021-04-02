abstract type DSState end


get_x(dss::DSState) = dss.x
get_n_states(dss::DSState) = length(dss.x)
function state_eqs(dss::DSState, time, mass, inertia, forces, moments, h) end


struct DSStateDot{S, N, T}
    dss::S
    xdot::SVector{N, T}
    DSStateDot(dss::S, xdot::V) where {S<:DSState, V<:SVector} =
        new{typeof(dss), get_n_states(dss), eltype(dss)}(
            dss,
            SVector{get_n_states(dss), eltype(dss)}(xdot)
        )
end

function DSStateDot(dss::DSState, xdot::AbstractArray)
    return DSStateDot(dss, SVector{get_n_states(dss), eltype(dss)}(xdot))
end


get_xdot(dssd::DSStateDot) = dssd.xdot
get_ds_state(dssd::DSStateDot) = dssd.dss

get_x(dssd::DSStateDot) = get_x(dssd.dss)
get_n_states(dssd::DSStateDot) = get_n_states(dssd.dss)


struct SixDOFAeroEuler{T}<:DSState
    x::SVector{13, T}
end

SixDOFAeroEuler(x::AbstractVector) = SixDOFAeroEuler(SVector{13, eltype(x)}(x))

eltype(::Type{SixDOFAeroEuler{T}}) where {T} = T

get_x_names(dss::SixDOFAeroEuler) = [:tas, :α, :β, :ϕ, :θ, :ψ, :p, :q, :r, :x, :y, :z, :pow]


function state_eqs(dss::SixDOFAeroEuler, time, mass, inertia, forces, moments, h, pow_dot)

    xdot = sixdof_aero_earth_euler_fixed_mass(
         time,
         get_x(dss),
         mass,
         inertia,
         forces,
         moments,
         h
    )

    xdot = [xdot..., pow_dot]

    return DSStateDot(dss, xdot)
end
