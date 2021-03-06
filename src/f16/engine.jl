function tgear(thtl)  # Power command v. thtl. relationship
    if thtl <= 0.77
        cpow = 64.94 * thtl
    else
        cpow = 217.38 * thtl - 117.38
    end
    return cpow
end


function pdot(p3, p1) # PDOT= rate of change of power
    if p1 >= 50.0  # P3= actual power, P1= power command
        if p3 >= 50.0
            t = 5.0
            p2 = p1
        else
            p2 = 60.0
            t = rtau(p2 - p3)
        end
    else
        if p3 >= 50.0
            t = 5.0
            p2 = 40.0
        else
            p2 = p1
            t = rtau(p2 - p3)
        end
    end
    pd = t * (p2 - p3)
    return pd
end


function rtau(dp)
    if dp <= 25.0
        rt = 1.0  # Reciprocal time constant
    elseif dp >= 50.0
        rt = 0.1
    else
        rt = 1.9 - 0.036 * dp
    end
    return rt
end


idle_data =
    [
        1060.0 670.0 880.0 1140.0 1500.0 1860.0
        635.0 425.0 690.0 1010.0 1330.0 1700.0
        60.0 25.0 345.0 755.0 1130.0 1525.0
        -1020.0 -710.0 -300.0 350.0 910.0 1360.0
        -2700.0 -1900.0 -1300.0 -247.0 600.0 1100.0
        -3600.0 -1400.0 -595.0 -342.0 -200.0 700.0
    ]'

mil_data =
    [
        12680.0 9150.0 6200.0 3950.0 2450.0 1400.0
        12680.0 9150.0 6313.0 4040.0 2470.0 1400.0
        12610.0 9312.0 6610.0 4290.0 2600.0 1560.0
        12640.0 9839.0 7090.0 4660.0 2840.0 1660.0
        12390.0 10176.0 7750.0 5320.0 3250.0 1930.0
        11680.0 9848.0 8050.0 6100.0 3800.0 2310.0
    ]'

max_data =
    [
        20000.0 15000.0 10800.0 7000.0 4000.0 2500.0
        21420.0 15700.0 11225.0 7323.0 4435.0 2600.0
        22700.0 16860.0 12250.0 8154.0 5000.0 2835.0
        24240.0 18910.0 13760.0 9285.0 5700.0 3215.0
        26070.0 21075.0 15975.0 11115.0 6860.0 3950.0
        28886.0 23319.0 18300.0 13484.0 8642.0 5057.0
    ]'


function thrust(pow, alt, rmach)  # Engine thrust model

    if alt <= 0.
        alt = 0.0
    end

    h = 0.0001 * alt
    i = floor(Int, h)
    if i >= 5
        i = 4
    end

    dh = h - float(i)

    rm = 5.0 * rmach
    m = floor(Int, rm)
    if m >= 5
        m = 4
    end

    dm = rm - float(m)
    cdh = 1.0 - float(dh)

    i = i + 1
    m = m + 1

    s = mil_data[i, m] * cdh + mil_data[i+1, m] * dh
    t = mil_data[i, m+1] * cdh + mil_data[i+1, m+1] * dh
    tmil = s + (t - s) * dm

    if pow < 50.0
        s = idle_data[i, m] * cdh + idle_data[i+1, m] * dh
        t = idle_data[i, m+1] * cdh + idle_data[i+1, m+1] * dh
        tidl = s + (t - s) * dm
        thrst = tidl + (tmil - tidl) * pow * 0.02
    else
        s = max_data[i, m] * cdh + max_data[i+1, m] * dh
        t = max_data[i, m+1] * cdh + max_data[i+1, m+1] * dh
        tmax = s + (t - s) * dm
        thrst = tmil + (tmax - tmil) * (pow - 50.0) * 0.02
    end

    return thrst
end


function calculate_prop_forces_moments(x, mach, controls)

    # Assign state & control variables
    vt = x[1] * M2FT
    alt = x[12] * M2FT
    pow = x[13]

    thtl = controls[1]

    thrust_ = thrust(pow, alt, mach)

    return [thrust_, 0.0, 0.0, 0.0, 0.0, 0.0]
end

calculate_prop_gyro_effects() = [HX, 0, 0]
