function rk4(f, DT, XX, TIME, XCG, CONTROLS)

    NX = length(XX)
    XA = Array{Float64}(undef, NX)
    X = Array{Float64}(undef, NX)

    # 1
    XD, outputs = f(TIME, XX, XCG, CONTROLS)

    for M = 1:NX
        XA[M] = XD[M] * DT
        X[M] = XX[M] + 0.5 * XA[M]
    end

    # 2
    XD, outputs = f(TIME, XX, XCG, CONTROLS)
    for M = 1:NX
        Q = XD[M] * DT
        X[M] = XX[M] + 0.5 * Q
        XA[M] = XA[M] + Q + Q
    end

    # 3
    XD, outputs = f(TIME, XX, XCG, CONTROLS)
    for M = 1:NX
        Q = XD[M] * DT
        X[M] = XX[M] + Q
        XA[M] = XA[M] + Q + Q
    end

    # 4
    XD, outputs = f(TIME, XX, XCG, CONTROLS)
    for M = 1:NX
        X[M] = XX[M] + (XA[M] + XD[M] * DT) / 6.0
    end

    return X
end
