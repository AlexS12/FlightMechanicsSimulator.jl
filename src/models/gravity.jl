# ---------------------------------------- Gravity ----------------------------------------
# Define abstract type and interface.
abstract type Gravity end


# TODO: make an interface which receives position. Right now, for flat Earth and constant
# gravity does not make sense. However, will be usefull in the future.
"""
    get_gravity_accel(g::Gravity)

Get gravity acceleration value (m/s²).
"""
function get_gravity_accel(g::Gravity) end


"""
    get_gravity_horizon(g::Gravity)

Get gravity vector in local horizon coordinates (North, East, Down) (m/s²).
"""
function get_gravity_horizon(g::Gravity) end


"""
    get_gravity_horizon(g::Gravity, θ, ϕ)

Get gravity vector in body axis (m/s²).
"""
function get_gravity_body(g::Gravity, θ, ϕ)
    ax, ay, az = horizon2body(get_gravity_horizon(g)..., 0.0, θ, ϕ)
    return @SVector [ax, ay, az]
end


# ------------------------------ Local Horizon Down Gravity -------------------------------
"""
    LHDownGravity{T<:Number} <: Gravity

Local Horizon Down gravity.

# Fields
- `gD::T`: down gravity acceleration value (m/s²).
"""
struct LHDownGravity{T<:Number} <: Gravity
    """Gravity acceleration (m/s²)"""
    gD::T
end


"""
    LHDownGravity([gD])

If no value of acceleration is given, default acceleration is: `gD=$(gD) m/s²`.
"""
LHDownGravity(gD=gD) = LHDownGravity(gD)


get_gravity_accel(g::LHDownGravity) = g.gD
get_gravity_horizon(g::LHDownGravity) = @SVector [0, 0, g.gD]
