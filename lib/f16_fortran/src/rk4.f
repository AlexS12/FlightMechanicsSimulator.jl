      SUBROUTINE RK4(F,DT,XX,NX,TIME,XCG,CONTROLS)

      EXTERNAL  F
      INTEGER(8) ::  NN, NX, M
      PARAMETER (NN=13)
      REAL(8) :: XX(NN),XD(NN),X(NN),XA(NN)
      REAL(8) :: TIME
      REAL(8) :: XCG
      REAL(8) :: CONTROLS(4)
      REAL(8) :: OUTPUT(7)


      CALL F(TIME, XX, XD, XCG, CONTROLS, OUTPUT)

      DO M=1,NX
        XA (M)=XD(M)*DT
        X(M)=XX(M)+0.5d0*XA(M)
      enddo

      CALL F(TIME, XX, XD, XCG, CONTROLS, OUTPUT)

      DO M=1,NX
        Q=XD(M)*DT
        X(M)=XX(M)+0.5*Q
        XA(M)=XA(M)+Q+Q
      enddo

      CALL F(TIME, XX, XD, XCG, CONTROLS, OUTPUT)

      DO M=1,NX
        Q=XD(M)*DT
        X(M)=XX(M)+Q
        XA(M)=XA(M)+Q+Q
      enddo

      CALL F(TIME, XX, XD, XCG, CONTROLS, OUTPUT)

      DO M=1,NX
        XX(M)=XX(M)+(XA(M)+XD(M)*DT)/6.d0
      enddo

      RETURN
      END