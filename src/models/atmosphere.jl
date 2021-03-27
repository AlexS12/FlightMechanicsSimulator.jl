abstract type Atmosphere end


"""
    get_temperature(a::Atmosphere)

Get atmosphere temperature (K).
"""
get_temperature(a::Atmosphere) = a.temperature


"""
    get_pressure(a::Atmosphere)

Get atmosphere pressure (Pa).
"""
get_pressure(a::Atmosphere) = a.pressure


"""
    get_density(a::Atmosphere)

Get atmosphere density (kg/m³).
"""
get_density(a::Atmosphere) = a.density


"""
    get_sound_velocity(a::Atmospere)

Get sound velocity (m/s).
"""
get_sound_velocity(a::Atmosphere) = a.sound_velocity


"""
    get_height(a::Atmospere)

Get geopotential height (m).
"""
get_height(a::Atmosphere) = a.height


"""
    F16StevensAtmosphere{T<:Number}<:Atmosphere

ISA atmosphere model as coded in [1].

# Fields
- `temperature::T`: temperature (K).
- `pressure::T`: pressure (Pa).
- `density::T`: density (kg/m³).
- `sound_velocity::T`: sound velocity (m/s).
- `height::T`: height (m).

# Notes
- Differs from [`F16StevensAtmosphere`](@ref). Check `FlightMechanicsUtils.jl` for more
  information.

# References
1. Stevens, B. L., Lewis, F. L., (1992). Aircraft control and simulation: dynamics, controls
   design, and autonomous systems. John Wiley & Sons. Section A.6 (page 715)
"""
struct F16StevensAtmosphere{T<:Number}<:Atmosphere
    """Temperature (K)"""
    temperature::T
    """Pressure (Pa)"""
    pressure::T
    """Density (kg/m³)"""
    density::T
    """Sound velocity (m/s)"""
    sound_velocity::T
    """Geopotential Height (m)"""
    height::T
end


"""
    F16StevensAtmosphere(h)

Calculate atmosphere for given height, `h` (m).

# Returns
- `F16StevensAtmosphere`
"""
function F16StevensAtmosphere(h)
    T, ρ, a, p = F16.atmosphere_f16(h)
    return F16StevensAtmosphere(T, p, ρ, a, h)
end

"""
    ISA1976{T<:Number}<:Atmosphere

ISA 1976 atmosphere implementation.

# Fields
- `temperature::T`: temperature (K).
- `pressure::T`: pressure (Pa).
- `density::T`: density (kg/m³).
- `sound_velocity::T`: sound velocity (m/s).
- `height::T`: height (m).
"""
struct ISA1976{T<:Number}<:Atmosphere
    """Temperature (K)"""
    temperature::T
    """Pressure (Pa)"""
    pressure::T
    """Density (kg/m³)"""
    density::T
    """Sound velocity (m/s)"""
    sound_velocity::T
    """Geopotential Height (m)"""
    height::T
end


"""
    ISA1976(h)

Calculate atmosphere for given height, `h` (m).

# Returns
- `ISA1976`
"""
ISA1976(h) = ISA1976(atmosphere_isa(h)..., h)
