abstract type DSState{T} end

eltype(::Type{<:DSState{T}}) where {T} = T

get_x(dss::DSState) = dss.x
get_n_states(dss::DSState) = length(dss.x)

function get_x_names(dss::DSState) end
function state_eqs(dss::DSState, time, mass, inertia, forces, moments, h) end

function get_earth_position(dss::DSState) end
function get_height(dss::DSState) end

function get_euler_angles(dss::DSState) end

function get_tas(dss::DSState) end
function get_α(dss::DSState) end
function get_β(dss::DSState) end
function get_tasαβ(dss::DSState) end
function get_body_velocity(dss::DSState) end
function get_horizon_velocity(dss::DSState) end

function get_ang_vel_body(dss::DSState) end
function get_euler_angles_rates(dss::DSState) end

function get_engine_power(dss::DSState) end


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
get_x_names(dssd::DSStateDot) = get_x_names(dssd.dss)
get_n_states(dssd::DSStateDot) = get_n_states(dssd.dss)

get_earth_position(dssd::DSStateDot) = get_earth_position(dssd.dss)
get_height(dssd::DSStateDot) = get_height(dssd.dss)
get_euler_angles(dssd::DSStateDot) = get_euler_angles(dssd.dss)
get_tas(dssd::DSStateDot) = get_tas(dssd.dss)
get_α(dssd::DSStateDot) = get_α(dssd.dss)
get_β(dssd::DSStateDot) = get_β(dssd.dss)
get_tasαβ(dssd::DSStateDot) = get_tasαβ(dssd.dss)
get_body_velocity(dssd::DSStateDot) = get_body_velocity(dssd.dss)
get_horizon_velocity(dssd::DSStateDot) = get_horizon_velocity(dssd.dss)
get_ang_vel_body(dssd::DSStateDot) = get_ang_vel_body(dssd.dss)
get_euler_angles_rates(dssd::DSStateDot) = get_euler_angles_rates(dssd.dss)
get_engine_power(dssd::DSStateDot) = get_engine_power(dssd.dss)

function get_tas_dot(dssd::DSStateDot) end
function get_α_dot(dssd::DSStateDot) end
function get_β_dot(dssd::DSStateDot) end
function get_tasαβ_dot(dssd::DSStateDot) end
function get_accel_body(dssd::DSStateDot) end

function get_ang_accel_body(dssd::DSStateDot) end

function get_engine_power_dot(dssd::DSStateDot) end
