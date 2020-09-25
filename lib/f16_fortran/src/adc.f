      SUBROUTINE ADC(VT,ALT,AMACH,QBAR)
      DATA R0/2.377E-3/
      TFAC = 1.0 - 0.703E-5 * ALT
      T = 519.0 * TFAC
      IF (ALT .GE. 35000.0) T= 390.0
      RHO = R0 * (TFAC**4.14)
      AMACH= VT/SQRT(1.4*1716.3*T)
      QBAR = 0.5*RHO*VT*VT
C     PS = 1715.0 * RHO * T
      RETURN
      END