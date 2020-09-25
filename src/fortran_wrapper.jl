using Libdl


const DLL = "../lib/f16_fortran/bin/f16"


function adc(vt::Float64, alt::Float64)

    amach = Array{Float64}(undef, 1)
    qbar = Array{Float64}(undef, 1)

    ccall(
        (:adc_, DLL),
        Nothing,
        (Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}),
        vt,
        alt,
        amach,
        qbar,
    )

    return amach[1], qbar[1]
end


function damp(α::Float64)
    D = Array{Float64}(undef, 9)
    ccall((:damp_, DLL), Nothing, (Ref{Float64}, Ptr{Float64}), α, D)
    return D
end

CX(α::Float64, de::Float64) =
    ccall((:cx_, DLL), Float64, (Ref{Float64}, Ref{Float64}), α, de)

CY(β::Float64, da::Float64, dr::Float64) =
    ccall((:cy_, DLL), Float64, (Ref{Float64}, Ref{Float64}, Ref{Float64}), β, da, dr)

CZ(α::Float64, β::Float64, de::Float64) =
    ccall((:cz_, DLL), Float64, (Ref{Float64}, Ref{Float64}, Ref{Float64}), α, β, de)

CM(α::Float64, de::Float64) =
    ccall((:cm_, DLL), Float64, (Ref{Float64}, Ref{Float64}), α, de)

CL(α::Float64, β::Float64) =
    ccall((:cl_, DLL), Float64, (Ref{Float64}, Ref{Float64}), α, β)

CN(α::Float64, β::Float64) =
    ccall((:cn_, DLL), Float64, (Ref{Float64}, Ref{Float64}), α, β)

DLDA(α::Float64, β::Float64) =
    ccall((:dlda_, DLL), Float64, (Ref{Float64}, Ref{Float64}), α, β)

DLDR(α::Float64, β::Float64) =
    ccall((:dldr_, DLL), Float64, (Ref{Float64}, Ref{Float64}), α, β)

DNDA(α::Float64, β::Float64) =
    ccall((:dnda_, DLL), Float64, (Ref{Float64}, Ref{Float64}), α, β)

DNDR(α::Float64, β::Float64) =
    ccall((:dndr_, DLL), Float64, (Ref{Float64}, Ref{Float64}), α, β)

tgear(thtl::Float64) = ccall((:tgear_, DLL), Float64, (Ref{Float64},), thtl)

rtau(dp::Float64) = ccall((:rtau_, DLL), Float64, (Ref{Float64},), dp)

pdot(p3::Float64, p1::Float64) =
    ccall((:pdot_, DLL), Float64, (Ref{Float64}, Ref{Float64}), p3, p1)

thrust(pow::Float64, alt::Float64, rmach::Float64) = ccall(
    (:thrust_, DLL),
    Float64,
    (Ref{Float64}, Ref{Float64}, Ref{Float64}),
    pow,
    alt,
    rmach,
)


function rk4(
    fun,
    dt::Float64,
    x::Array{Float64,1},
    time::Float64,
    xcg::Float64,
    controls::Array{Float64,1},
)

    x_ = copy(x)
    nx = length(x_)

    ccall(
        (:rk4_, DLL),
        Nothing,
        (
            Ptr{Nothing},
            Ref{Float64},
            Ptr{Float64},
            Ref{Int64},
            Ref{Float64},
            Ref{Float64},
            Ptr{Float64},
        ),
        fun,
        dt,
        x_,
        nx,
        time,
        xcg,
        controls,
    )

    return x_
end


function f(
    time::Float64,
    x::Array{Float64,1},
    xcg::Float64,
    controls::Array{Float64,1},
)

    xd = Array{Float64}(undef, 13)
    outputs = Array{Float64}(undef, 7)

    ccall(
        (:f_, DLL),
        Nothing,
        (
            Ref{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ref{Float64},
            Ptr{Float64},
            Ptr{Float64},
        ),
        time,
        x,
        xd,
        xcg,
        controls,
        outputs,
    )

    return xd, outputs
end
