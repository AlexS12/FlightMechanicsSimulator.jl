# --------------------------------  Dynamic System State --------------------------------
"""
    DSState{T}
    DSState(x::AbstractVector)

Dynamic system state of type `T` (type of the elements of the state vector, ie. `Float64`).
"""
abstract type DSState{T} end


eltype(::Type{<:DSState{T}}) where {T} = T

"""
    get_x(dss::DSState)
    get_x(dssd::DSStateDot)

Return the state vector `x` of the dynamic system.
"""
get_x(dss::DSState) = dss.x

"""
    get_n_states(dss::DSState)
    get_n_states(dssd::DSStateDot)

Return the length of the state vector `x` of the dynamic system.
"""
get_n_states(dss::DSState) = length(dss.x)


# ----------------------- Methods to be implemented by each subtype -----------------------
"""
    get_x_names(dss::DSState)
    get_x_names(dssd::DSStateDot)

Return an array of symbols with the name for each variable in the state vector.
"""
function get_x_names end

"""
    state_eqs(dss::DSState, time, mass, inertia, forces, moments, h)

Calculate the derivative of the state vector using the state equations.

# Arguments

- `dss::DSState`: dynamic system state.
- `time::Number`: time (s).
- `mass::Number`: mass (kg).
- `inertia::AbstractMatrix`: 3x3 inertia tensor (kg·m^2).
- `forces::AbstractVector`: body axis forces [Fx, Fy, Fz] (N).
- `moments::AbstractVector`: body axis moments (L, M, N) (N·m).
- `h::AbstractVector`:  Additional angular momentum contributions such as those coming from
 spinning rotors (kg·m²/s).

 # Returns
 - `dssd::DSStateDot`: dynamic system state dot.
"""
function state_eqs end

"""
    get_earth_positon(dss::DSState)
    get_earth_positon(dssd::DSStateDot)

Return Earth position [x, y, z] (m).
"""
function get_earth_position end

"""
    get_height(dss::DSState)
    get_height(dssd::DSStateDot)

Return height (m).
"""
function get_height end

"""
    get_euler_angles(dss::DSState)
    get_euler_angles(dssd::DSStateDot)

Return Euler angles (ψ, θ, ϕ) (rad).
"""
function get_euler_angles end

"""
    get_tas(Dss::DState)
    get_tas(dssd::DSStateDot)

Return true air speed (m/s).
"""
function get_tas end

"""
    get_α(dss::DSState)
    get_α(dssd::DSStateDot)

Get angle of attack (rad).
"""
function get_α end

"""
    get_β(dss::DSState)
    get_β(dssd::DSStateDot)

Get angle of side slip (rad).
"""
function get_β end

"""
    get_tasαβ(dss:DSState)
    get_tasαβ(dssd::DSStateDot)

Get [tas (m/s), α (rad), β (rad)].
"""
function get_tasαβ end

"""
    get_body_velocity(dss::DState)
    get_body_velocity(dssd::DSStateDot)

Get velocity in body axis [u, v, w] (m/s).
"""
function get_body_velocity end

"""
    get_horizon_velocity(dss::DSState)
    get_horizon_velocity(dssd::DSStateDot)

Get velocity in horizon axis [Vn, Ve, Vd] (m/s).
"""
function get_horizon_velocity end

"""
    get_ang_vel_body(dss::DSState)
    get_ang_vel_body(dssd::DSStateDot)

Get angular velocity in body axis [p, q, r] (rad/s).
"""
function get_ang_vel_body end

"""
    get_euler_angles_rates(dss::DSState)
    get_euler_angles_rates(dssd::DSStateDot)

Get euler angles rates [ψ_dot, θ_dot, ϕ_dot] (rad/s).
"""
function get_euler_angles_rates end

"""
    get_engine_power(dss:DSState)
    get_engine_power(dssd::DSStateDot)

Get engine power (%).
"""
function get_engine_power end


# ------------------------------- Dynamic System State  Dot ------------------------------
"""
    DSStateDot{S, N, T}
    DSStateDot(dss::DSState, xdot::AbstractArray)

Time derivative, `x_dot`, of the state vector `x`.

- `S`: type of the associated `DSState`.
- `N`: number of elements of the `x_dot` vector.
- `T`: type of the elements of `x_dot` vector.
"""
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


"""
    get_xdot(dssd::DSStateDot)

Get the time derivative of state vector, `x_dot`.
"""
get_xdot(dssd::DSStateDot) = dssd.xdot

"""
    get_ds_state(dssd::DSStateDot)

Get the `DSState`.
"""
get_ds_state(dssd::DSStateDot) = dssd.dss

# --------------------- Methoos reexported from coposition with DSState --------------------
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


# ----------------------- Methods to be implemented by each subtype -----------------------
"""
    get_tas_dot(dssd::DSStateDot)

Get derivative of true air speed (m/s²).
"""
function get_tas_dot end

"""
    get_α_dot(dssd::DSStateDot)

Get angle of attack derivative (rad/s).
"""
function get_α_dot end

"""
    get_β_dot(dssd::DSStateDot)

Get derivative of angles of sideslip (rad/s).
"""
function get_β_dot end

"""
    get_tasαβ_dot(dssd::DSStateDot)

Get [tas_dot (m/s²), α_dot (rad/s), β_dot (rad/s)].
"""
function get_tasαβ_dot end

"""
    get_accel_body(dssd::DSStateDot)

Get accelerations in body axis [u_dot, v_dot, w_dot] (m/s²).
"""
function get_accel_body end

"""
    get_ang_accel_body(dssd::DSStateDot)

Get angular acceleration in body axis [p_dot, q_dot, r_dot] (rad/s²).
"""
function get_ang_accel_body end

"""
    get_engine_power_dot(dssd::DSStateDot)

Get derivative of engine power (s⁻¹).
"""
function get_engine_power_dot end


"""
    get_trimmer_cost(dssd::DSStateDot)

Get the cost function evaluation for the given `DSStateDot`.
"""
function get_trimmer_cost end
