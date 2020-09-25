      SUBROUTINE F(TIME,X,XD, XCG, CONTROLS, OUTPUT)
      REAL X(13), XD(13), D(9), MASS
      ! COMMON/PARAM/XCG
      ! COMMON/CONTROLS/THTL,EL,AIL,RDR
      ! COMMON/OUTPUT/AN,ALAT,AX,QBAR,AMACH,Q,ALPHA
      REAL  XCG
      REAL  CONTROLS(4)
      REAL  OUTPUT(7)

      PARAMETER (AXX=9496.0, AYY= 55814.0, AZZ=63100.0, AXZ= 982.0)
      PARAMETER (AXZS=AXZ**2, XPQ=AXZ*(AXX-AYY+AZZ),GAM=AXX*AZZ-AXZ**2)
      PARAMETER (XQR= AZZ*(AZZ-AYY)+AXZS, ZPQ=(AXX-AYY)*AXX+AXZS)
      PARAMETER ( YPR= AZZ - AXX )
      PARAMETER (WEIGHT= 20500.0, GD= 32.17, MASS= weight/gd)
      DATA S,B,CBAR,XCGR,HX/300,30,11.32,0.35,160.0/
      DATA RTOD / 57.29578/

      !! UNPACK CONTROLS
      THTL = CONTROLS(1)
      EL = CONTROLS(2)
      AIL = CONTROLS(3)
      RDR = CONTROLS(4)
C
C Assign state & control variables
C
C     X(1)  -> vt (ft/s)
C     X(2)  -> alpha (rad)
C     X(3)  -> beta (rad)
C     X(4)  -> phi (rad)
C     X(5)  -> theta (rad)
C     X(6)  -> psi (rad)
C     X(7)  -> P (rad/s)
C     X(8)  -> Q (rad/s)
C     X(9)  -> R (rad/s)
C     X(10) -> North (ft)
C     X(11) -> East (ft)
C     X(12) -> Altitude (ft)
C     X(13) -> Pow

      VT= X(1); ALPHA= X(2)*RTOD; BETA= X(3)*RTOD
      PHI=X(4); THETA= X(5); PSI= X(6)
      P= X(7); Q= X(8); R= X(9); ALT= X(12); POW= X(13)
C
C Air data computer and engine model
C
      CALL ADC(VT,ALT,AMACH,QBAR); CPOW= TGEAR(THTL)
      XD(13) = PDOT(POW,CPOW);
      T= THRUST(POW,ALT,AMACH)
C
C Look-up tables and component buildup
C
      CXT = CX (ALPHA,EL)
      CYT = CY (BETA,AIL,RDR)
      CZT = CZ (ALPHA,BETA,EL)
      DAIL= AIL/20.0; DRDR= RDR/30.0
      CLT = CL(ALPHA,BETA) + DLDA(ALPHA,BETA)*DAIL +
     1      DLDR(ALPHA,BETA)*DRDR
      CMT = CM(ALPHA,EL)
      CNT = CN(ALPHA,BETA) + DNDA(ALPHA,BETA)*DAIL +
     1      DNDR(ALPHA,BETA)*DRDR
C
C Add damping derivatives :
C
      CBTA = COS(X(3)); U=VT*COS(X(2))*CBTA
      V= VT * SIN(X(3)); W=VT*SIN(X(2))*CBTA
      TVT= 0.5/VT; B2V= B*TVT; CQ= CBAR*Q*TVT
      CALL DAMP(ALPHA,D)
      CXT= CXT + CQ * D(1)
      CYT= CYT + B2V * ( D(2)*R + D(3)*P )
      CZT= CZT + CQ * D(4)
      CLT= CLT + B2V * ( D(5)*R + D(6)*P )
      CMT= CMT + CQ * D(7) + CZT * (XCGR-XCG)
      CNT= CNT + B2V*(D(8)*R + D(9)*P) - CYT*(XCGR-XCG) * CBAR/B
C
C Get ready for state equations
C
      STH = SIN(THETA); CTH= COS(THETA); SPH= SIN(PHI)
      CPH = COS(PHI) ; SPSI= SIN(PSI); CPSI= COS(PSI)
      QS = QBAR * S ; QSB= QS * B; RMQS= QS/MASS
      GCTH = GD * CTH ; QSPH= Q * SPH
      AY = RMQS*CYT ; AZ= RMQS * CZT
C
C Force equations
C
      UDOT = R*V - Q*W - GD*STH + (QS * CXT + T)/MASS
      VDOT = P*W - R*U + GCTH * SPH + AY
      WDOT = Q*U - P*V + GCTH * CPH + AZ
      DUM = (U*U + W*W)
      xd(1) = (U*UDOT + V*VDOT + W*WDOT)/VT
      xd(2) = (U*WDOT - W*UDOT) / DUM
      xd(3) = (VT*VDOT- V*XD(1)) * CBTA / DUM
C
C Kinematics
C
      xd(4) = P + (STH/CTH)*(QSPH + R*CPH)
      xd(5) = Q*CPH - R*SPH
      xd(6) = (QSPH + R*CPH)/CTH
C
C Moments
C
      ROLL = QSB*CLT
      PITCH = QS *CBAR*CMT
      YAW = QSB*CNT
      PQ = p*Q
      QR = Q*R
      QHX = Q*HX

      xd(7) = ( XPQ*PQ - XQR*QR + AZZ*ROLL + AXZ*(YAW + QHX) )/GAM
      xd(8) = ( YPR*P*R - AXZ*(P**2 - R**2) + PITCH - R*HX )/AYY
      xd(9) = ( ZPQ*PQ - XPQ*QR + AXZ*ROLL + AXX*(YAW + QHX) )/GAM
C
C Navigation
C
      T1= SPH * CPSI; T2= CPH * STH; T3= SPH * SPSI
      S1= CTH * CPSI; S2= CTH * SPSI; S3= T1 * STH - CPH * SPSI
      S4= T3 * STH + CPH * CPSI; S5= SPH * CTH; S6= T2*CPSI + T3
      S7= T2 * SPSI - T1; S8= CPH * CTH
C
      xd(10) = U * S1 + V * S3 + W * S6 ! North speed
      xd(11) = U * S2 + V * S4 + W * S7 ! East speed
      xd(12) = U * STH -V * S5 - W * S8 ! Vertical speed
C
C Outputs
C
      AN = -AZ/GD; ALAT= AY/GD
      AX=(QS * CXT + T)/GD  ! <<-- ASM: Definition missing

      !! PACK OUTPUTS
      OUTPUT(1) = AN
      OUTPUT(2) = ALAT
      OUTPUT(3) = AX
      OUTPUT(4) = QBAR
      OUTPUT(5) = AMACH
      OUTPUT(6) = Q
      OUTPUT(7) = ALPHA
      RETURN
      END
