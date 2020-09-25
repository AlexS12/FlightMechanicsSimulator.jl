using Test
using FlightMechanicsSimulator


#  Stevens, B. L., Lewis, F. L., & Johnson, E. N. (2015). Aircraft control
#  and simulation: dynamics, controls design, and autonomous systems. John Wiley
#  & Sons. (page 185 table 3.5-2)

# xcg = 0.4 * CMA

# INPUTS
# U(1) = THTL =   0.9  [-]
# U(2) = ELEV =  20    [DEG]
# U(3) = DAIL = -15    [DEG]
# U(4) = RDR  = -20    [DEG]

# STATE
# IDX     name      X        UNITS        XDOT           XDOT MORELLI
# 1       VT       500       [ft/s]      -75.23724       -77.57521
# 2       ALPHA      0.5     [RAD]       -0.8813491      -0.88123
# 3       BETA      -0.2     [RAD]       -0.4759990      -0.45276
# 4       PHI       -1       [RAD]        2.505734        0.70000
# 5       THETA      1       [RAD]        0.3250820       0.32508
# 6       PSI       -1       [RAD]        2.145926        2.14593
# 7       P          0.7     [RAD/S]     12.62679        12.91108
# 8       Q         -0.8     [RAD/S]      0.9649671      0.97006
# 9       R          0.9     [RAD/S]      0.5809759      -0.55450
# 10      X       1000       [FT]       342.4439         342.44390
# 11      Y        900       [FT]      -266.7707         -266.77068
# 12      ALT    10000       [FT]       248.1241          248.12412
# 13      POW       90       [%]        -58.6899         -58.6900

t = 0.0  # s

x = [500.0, 0.5, -0.2, -1.0, 1.0, -1.0, 0.7, -0.8, 0.9, 1000.0, 900.0, 10000.0, 90.0]
xcg = 0.40
controls = [0.9, 20.0, -15.0, -20.0]

xd_exp = [
    -75.23724,
    -0.8813491,
    -0.4759990,
    2.505734,
    0.3250820,
    2.145926,
    12.62679,
    0.9649671,
    0.5809759,
    342.4439,
    -266.7707,
    248.1241,
    -58.6899,
]

xd_1, outputs_1 = FlightMechanicsSimulator.Fortran.f(t, x, xcg, controls)
xd_2, outputs_2 = FlightMechanicsSimulator.f(t, x, xcg, controls)

# Check Fortran against Julia
@test isapprox(xd_1, xd_2)
@test isapprox(outputs_1, outputs_2)
# Check against Stevens
@test isapprox(xd_1, xd_exp, atol=0.05)

# -----------------------------------
t = 0.0;

vt_test = [50.0, 75.0]  # ft/s
α_test = deg2rad.([1.0, 5.0])
β_test = deg2rad.([-5.0, 5.0])
ϕ_test = deg2rad.([-30.0, 25.0])
θ_test = deg2rad.([-15.0, 25.0])
ψ_test = deg2rad.([45.0, 175.0])
p_test = deg2rad.([-15.0, 30.0])
q_test = deg2rad.([-5.0, 10.0])
r_test = deg2rad.([-20.0, 30.0])
norh_ft = 0.0  # ft
east_ft = 0.0  # ft
alt_test = [5000.0, 45000.0]  # ft
pow_test = [10., 80.]  # %

de_test = [-25.0, 20.0]  # deg
da_test = [-15.0, 10.0]  # deg
dr_test = [-25.0, 20.0]  # deg
thtl_test = [0.2, 0.6]

xcg_test = [0.35, 0.25]

for vt_fts in vt_test
    for α_rad in α_test
        for β_rad in β_test
            for ϕ_rad in ϕ_test
                for θ_rad in ϕ_test
                    for ψ_rad in ψ_test
                        for p_rads in p_test
                            for q_rads in q_test
                                for r_rads in r_test
                                    for alt_ft in alt_test
                                        for pow in pow_test
                                            for thtl in thtl_test
                                                for el in de_test
                                                    for ail in da_test
                                                        for rdr in dr_test
                                                            for pow in thtl_test
                                                                for xcg in xcg_test
                                                                    local x
                                                                    local controls
                                                                    x = [
                                                                        vt_fts,
                                                                        α_rad,
                                                                        β_rad,
                                                                        ϕ_rad,
                                                                        θ_rad,
                                                                        ψ_rad,
                                                                        p_rads,
                                                                        q_rads,
                                                                        r_rads,
                                                                        norh_ft,
                                                                        east_ft,
                                                                        alt_ft,
                                                                        pow,
                                                                    ]

                                                                    controls =
                                                                        [thtl, el, ail, rdr]

                                                                    xd1, outputs1 =
                                                                    FlightMechanicsSimulator.Fortran.f(
                                                                            t,
                                                                            x,
                                                                            xcg,
                                                                            controls,
                                                                        )

                                                                    xd2, outputs2 =
                                                                        FlightMechanicsSimulator.f(
                                                                            t,
                                                                            x,
                                                                            xcg,
                                                                            controls,
                                                                        )

                                                                    @test isapprox(
                                                                        xd1,
                                                                        xd2,
                                                                        nans = true,
                                                                        atol = 1e-10,
                                                                    )
                                                                    @test isapprox(
                                                                        outputs1,
                                                                        outputs2,
                                                                        nans = true,
                                                                        atol = 1e-10,
                                                                    )
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
