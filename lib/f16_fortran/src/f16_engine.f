C
C  This function computes the engine
C  power level command, POW, for an input throttle setting, THTL,
C  for the F-16 engine model.
C
C  THTL = THROTTLE SETTING.  ( 0 <= THTL <= 1.0 )
C
C
      FUNCTION TGEAR(THTL) ! Power command v. thtl. relationship
      IF(THTL.LE.0.77) THEN
      TGEAR = 64.94*THTL
      ELSE
      TGEAR = 217.38*THTL-117.38
      END IF
      RETURN
      END

C
C  This function computes the rate of change in engine power level
C  using a first order lag as a function of actual power, POW, 
C  and commanded power, CPOW, for the F-16 engine model.
C
C  POW =  ENGINE POWER LEVEL, PERCENT.  ( 0 <= POW <= 100. )
C  CPOW = COMMANDED ENGINE POWER LEVEL, PERCENT.  ( 0 <= CPOW <= 100. )
C
C
      FUNCTION PDOT(P3,P1) ! PDOT= rate of change of power
      IF (P1.GE.50.0) THEN ! P3= actual power, P1= power command
      IF (P3.GE.50.0) THEN
      T=5.0
      P2=P1
      ELSE
      P2=60.0
      T=RTAU(P2-P3)
      END IF
      ELSE
      IF (P3.GE.50.0) THEN
      T=5.0
      P2=40.0
      ELSE
      P2=P1
      T=RTAU(P2-P3)
      END IF
      END IF
      PDOT=T*(P2-P3)
      RETURN
      END


C
C  This function computes the thrust lag reciprocal time constant
C  for the F-16 engine model.
C
C  DP = CHANGE IN POWER LEVEL, PERCENT  ( 0 <= DP <= 100. )
C
C
      FUNCTION RTAU(DP) 
      ! used by function PDOT
      IF (DP.LE.25.0) THEN
      RTAU=1.0
      ! reciprocal time constant
      ELSE IF (DP.GE.50.0)THEN
      RTAU=0.1
      ELSE
      RTAU=1.9-.036*DP
      END IF
      RETURN
      END

C
C  This function computes the thrust 
C  for the F-16 model.
C
C  POW = ENGINE POWER LEVEL, PERCENT.  ( 0 <= POW <= 100. )
C  ALT = ALTITUDE, FT.        ( 0 <= ALT <= 50000. )
C  RMACH = MACH NUMBER.       ( 0 <= RMACH <= 1.0 )
C
C
      FUNCTION THRUST(POW,ALT,RMACH)

      REAL A(0:5,0:5), B(0:5,0:5), C(0:5,0:5)
C
C  IDLE POWER DATA.
C
      DATA A/  1060.,  670.,  880., 1140., 1500., 1860.,
     $          635.,  425.,  690., 1010., 1330., 1700.,
     $           60.,   25.,  345.,  755., 1130., 1525.,
     $        -1020., -710., -300.,  350.,  910., 1360.,
     $        -2700.,-1900.,-1300., -247.,  600., 1100.,
     $        -3600.,-1400., -595., -342., -200.,  700./
C
C  MIL POWER DATA.
C
      DATA B/ 12680., 9150., 6200., 3950., 2450., 1400.,
     $        12680., 9150., 6313., 4040., 2470., 1400.,
     $        12610., 9312., 6610., 4290., 2600., 1560.,
     $        12640., 9839., 7090., 4660., 2840., 1660.,
     $        12390.,10176., 7750., 5320., 3250., 1930.,
     $        11680., 9848., 8050., 6100., 3800., 2310./
C
C  MAX POWER DATA.
C
      DATA C/ 20000.,15000.,10800., 7000., 4000., 2500.,
     $        21420.,15700.,11225., 7323., 4435., 2600.,
     $        22700.,16860.,12250., 8154., 5000., 2835.,
     $        24240.,18910.,13760., 9285., 5700., 3215.,
     $        26070.,21075.,15975.,11115., 6860., 3950.,
     $        28886.,23319.,18300.,13484., 8642., 5057./
C
C
C  ROW INDEX FOR ALTITUDE.
C
      H=0.0001*ALT
      I=INT(H)
      IF (I.GE.5) I=4
      DH=H-FLOAT(I)
C
C  COLUMN INDEX FOR MACH NUMBER.
C
      RM=5.*RMACH
      M=INT(RM)
      IF (M.GE.5) M=4
      DM=RM-FLOAT(M)
      CDH=1.0-DH
C
C  COMPUTE MIL THRUST.
C
C  ALTITUDE INTERPOLATION.
      S= B(I,M)*CDH + B(I+1,M)*DH
      T= B(I,M+1)*CDH + B(I+1,M+1)*DH
C  MACH NUMBER INTERPOLATION.
      TMIL= S + (T-S)*DM
C
C  INTERPOLATE WITH IDLE OR MAX THRUST, DEPENDING ON POWER LEVEL COMMAND.
C
      IF (POW.LT.50.0) THEN
C
C  COMPUTE IDLE THRUST.
C
C  ALTITUDE INTERPOLATION.
        S= A(I,M)*CDH + A(I+1,M)*DH
        T= A(I,M+1)*CDH + A(I+1,M+1)*DH
C  MACH NUMBER INTERPOLATION.
        TIDL= S + (T-S)*DM
        THRUST= TIDL + (TMIL-TIDL)*POW/50.0
      ELSE
C
C  COMPUTE MAX THRUST.
C
C  ALTITUDE INTERPOLATION.
        S= C(I,M)*CDH + C(I+1,M)*DH
        T= C(I,M+1)*CDH + C(I+1,M+1)*DH
C  MACH NUMBER INTERPOLATION.
        TMAX= S + (T-S)*DM
        THRUST= TMIL + (TMAX-TMIL)*(POW-50.0)*0.02
      END IF
C
      RETURN
      END