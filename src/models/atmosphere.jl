abstract type Atmosphere end


get_temperature(a::Atmosphere) = a.temperature
get_pressure(a::Atmosphere) = a.pressure
get_density(a::Atmosphere) = a.density
get_sound_velocity(a::Atmosphere) = a.sound_velocity


function calculate_atmosphere end


struct F16StevensAtmosphere{T<:Number}<:Atmosphere
    """Temperature (K)"""
    temperature::T
    """Pressure (Pa)"""
    pressure::T
    """Density (kg/m³)"""
    density::T
    """Sound velocity (m/s)"""
    sound_velocity::T
end


function F16StevensAtmosphere(h)
    T, ρ, a, p = F16.atmosphere_f16(h)
    return F16StevensAtmosphere(T, p, ρ, a)
end

struct ISA1976{T<:Number}<:Atmosphere
    """Temperature (K)"""
    temperature::T
    """Pressure (Pa)"""
    pressure::T
    """Density (kg/m³)"""
    density::T
    """Sound velocity (m/s)"""
    sound_velocity::T
end


ISA1976(h) = ISA1976(atmosphere_isa(h)...)
