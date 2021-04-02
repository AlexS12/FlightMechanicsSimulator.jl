abstract type DSState{T} end

eltype(::Type{<:DSState{T}}) where {T} = T

get_x(dss::DSState) = dss.x
get_n_states(dss::DSState) = length(dss.x)

function get_x_names(dss::DSState) end
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
