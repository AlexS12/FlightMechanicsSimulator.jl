      SUBROUTINE DAMP(ALPHA, D)
C
C  This subroutine computes the damping coefficients
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C
C  D = [CXq, CYr, CYp, CZq, Clr, Clp, Cmq, Cnr, Cnp]
C
      REAL A(-2:9,9),D(9)
C
      DATA A/
     & -.267, -.110, .308, 1.34, 2.08, 2.91, 2.76,
     & 2.05, 1.50, 1.49, 1.83, 1.21,
     & .882, .852, .876, .958, .962, .974, .819,
     & .483, .590, 1.21, -.493, -1.04,
     & -.108, -.108, -.188, .110, .258, .226, .344,
C    & .362, .611, .529, .298, -2.27,   <<--
     & .362, .611, .529, .298, -0.227, !<<--  Morelli
     & -8.80, -25.8, -28.9, -31.4, -31.2, -30.7, -27.7,
     & -28.2, -29.0, -29.8, -38.3, -35.3,
     & -.126, -.026, .063, .113, .208, .230, .319,
     & .437, .680, .100, .447, -.330,
     & -.360, -.359, -.443, -.420, -.383, -.375, -.329,
     & -.294, -.230, -.210, -.120, -.100,
C    & -7.21, -.540, -5.23, -5.26, -6.11, -6.64, -5.69,  <<--
     & -7.21, -5.40, -5.23, -5.26, -6.11, -6.64, -5.69, !<<--  Morelli
     & -6.00, -6.20, -6.40, -6.60, -6.00,
     & -.380, -.363, -.378, -.386, -.370, -.453, -.550,
     & -.582, -.595, -.637, -1.02, -.840,
     & .061, .052, .052, -.012, -.013, -.024, .050,
     & .150, .130, .158, .240, .150/
C
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      DO 10  I=1,9
      D(I) = A(K,I) + ABS(DA)*(A(L,I)-A(K,I))
 10   CONTINUE
C
      RETURN
      END
C

C
      FUNCTION CX(ALPHA,EL)
C
C  This function computes the X body axis aerodynamic force coefficient
C  for the F-16 aerodynamic model.
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  EL = ELEVATOR DEFLECTION, DEG  ( -25 <= EL <= 25 )
C
C
      REAL A(-2:9,-2:2)
C
      DATA A/  -.099,-.081,-.081,-.063,-.025,.044,
     $.097,.113,.145,.167,.174,.166,
     $-.048,-.038,-.040,-.021,.016,.083,
     $.127,.137,.162,.177,.179,.167,
     $-.022,-.020,-.021,-.004,.032,.094,
     $.128,.130,.154,.161,.155,.138,
     $-.040,-.038,-.039,-.025,.006,.062,
     $.087,.085,.100,.110,.104,.091,
     $-.083,-.073,-.076,-.072,-.046,.012,
     $.024,.025,.043,.053,.047,.040 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = EL/12.0
      M = INT(S)
      IF(M.LE.-2) M=-1
      IF(M.GE.2) M=1
      DE = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DE))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      CX = V + (W-V)*ABS(DE)
C
      RETURN
      END
C
C
      FUNCTION CY(BETA,AIL,RDR)
C
C  This function computes the Y body axis aerodynamic force coefficient
C  for the F-16 aerodynamic model. 
C
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C  AIL = AILERON DEFLECTION, DEG  ( -21.5 <= AIL <= 21.5 )
C  RDR = RUDDER DEFLECTION, DEG  ( -30 <= RDR <= 30 )
C
C
      CY = -.02*BETA + .021*(AIL/20.0) + .086*(RDR/30.0)
C
      RETURN
      END
C
C
      FUNCTION CZ(ALPHA,BETA,EL)
C
C  This function computes the Z body axis aerodynamic force coefficient
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C  EL = ELEVATOR DEFLECTION, DEG  ( -25 <= EL <= 25 )
C
C
      REAL A(-2:9)
C
      DATA A /  .770,.241,-.100,-.416,-.731,-1.053,
     &         -1.366,-1.646,-1.917,-2.120,-2.248,-2.229 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = A(K) + ABS(DA)*(A(L)-A(K))
      CZ = S*(1-(BETA/57.3)**2) - .19*(EL/25.0)
C
      RETURN
      END
C
C
      FUNCTION CM(ALPHA,EL)
C
C  This function computes the Y body axis aerodynamic moment coefficient
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  EL = ELEVATOR DEFLECTION, DEG  ( -25 <= EL <= 25 )
C
C
      REAL A(-2:9,-2:2)
C
      DATA A/  .205,.168,.186,.196,.213,.251,
     $.245,.238,.252,.231,.198,.192,
     $.081,.077,.107,.110,.110,.141,
     $.127,.119,.133,.108,.081,.093,
     $-.046,-.020,-.009,-.005,-.006,.010,
     $.006,-.001,.014,.000,-.013,.032,
     $-.174,-.145,-.121,-.127,-.129,-.102,
     $-.097,-.113,-.087,-.084,-.069,-.006,
     $-.259,-.202,-.184,-.193,-.199,-.150,
     $-.160,-.167,-.104,-.076,-.041,-.005/
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = EL/12.0
      M = INT(S)
      IF(M.LE.-2) M=-1
      IF(M.GE.2) M=1
      DE = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DE))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      CM = V + (W-V)*ABS(DE)
C
      RETURN
      END
C
C
      FUNCTION CL(ALPHA,BETA)
C
C  This function computes the X body axis aerodynamic moment coefficient
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      REAL A(-2:9,0:6)
C
C      DATA A/  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
C     $-.001,-.004,-.008,-.012,-.016,-.022,
C     $-.022,-.021,-.015,-.008,-.013,-.015,
C     $-.003,-.009,-.017,-.024,-.030,-.041,
C     $-.045,-.040,-.016,-.002,-.010,-.019,
C     $-.001,-.010,-.020,-.030,-.039,-.054,
C     $-.057,-.054,-.023,-.006,-.014,-.027,
C     $.000,-.010,-.022,-.034,-.047,-.060,
C     $-.069,-.067,-.033,-.036,-.035,-.035,
C     $.007,-.010,-.023,-.034,-.049,-.063,
C     $-.081,-.079,-.060,-.058,-.062,-.059,
C     $.009,-.011,-.023,-.037,-.050,-.068,
C     $-.089,-.088,-.091,-.076,-.077,-.076 /
C
      DATA A/12*0,
     & -.001, -.004, -.008, -.012, -.016, -.019, -.020,
     & -.020, -.015, -.008, -.013, -.015,
     & -.003, -.009, -.017, -.024, -.030, -.034, -.040,
     & -.037, -.016, -.002, -.010, -.019,
     & -.001, -.010, -.020, -.030, -.039, -.044, -.050,
     & -.049, -.023, -.006, -.014, -.027,
     & .000, -.010, -.022, -.034, -.047, -.046, -.059,
     & -.061, -.033, -.036, -.035, -.035,
     & .007, -.010, -.023, -.034, -.049, -.046, -.068,
     & -.071, -.060, -.058, -.062, -.059,
     & .009, -.011, -.023, -.037, -.050, -.047, -.074,
     & -.079, -.091, -.076, -.077, -.076/

      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = .2*ABS(BETA)
      M = INT(S)
      IF(M.EQ.0) M=1
      IF(M.GE.6) M=5
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      CL = (V + (W-V)*ABS(DB))*SIGN(1.0,BETA)
C
      RETURN
      END
C
C
      FUNCTION CN(ALPHA,BETA)
C
C  This function computes the Z body axis aerodynamic moment coefficient
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      REAL A(-2:9,0:6)
C
      DATA A/  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     $.018,.019,.018,.019,.019,.018,
     $.013,.007,.004,-.014,-.017,-.033,
     $.038,.042,.042,.042,.043,.039,
     $.030,.017,.004,-.035,-.047,-.057,
     $.056,.057,.059,.058,.058,.053,
     $.032,.012,.002,-.046,-.071,-.073,
     $.064,.077,.076,.074,.073,.057,
     $.029,.007,.012,-.034,-.065,-.041,
     $.074,.086,.093,.089,.080,.062,
     $.049,.022,.028,-.012,-.002,-.013,
     $.079,.090,.106,.106,.096,.080,
     $.068,.030,.064,.015,.011,-.001 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = .2*ABS(BETA)
      M = INT(S)
      IF(M.EQ.0) M=1
      IF(M.GE.6) M=5
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      CN = (V + (W-V)*ABS(DB))*SIGN(1.0,BETA)
C
      RETURN
      END
C
C
      FUNCTION DLDA(ALPHA,BETA)
C
C  This function computes the rolling moment due to aileron deflection
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      REAL A(-2:9,-3:3)
C
      DATA A/  -.041,-.052,-.053,-.056,-.050,-.056,
     $-.082,-.059,-.042,-.038,-.027,-.017,
     $-.041,-.053,-.053,-.053,-.050,-.051,
     $-.066,-.043,-.038,-.027,-.023,-.016,
     $-.042,-.053,-.052,-.051,-.049,-.049,
     $-.043,-.035,-.026,-.016,-.018,-.014,
     $-.040,-.052,-.051,-.052,-.048,-.048,
     $-.042,-.037,-.031,-.026,-.017,-.012,
     $-.043,-.049,-.048,-.049,-.043,-.042,
     $-.042,-.036,-.025,-.021,-.016,-.011,
     $-.044,-.048,-.048,-.047,-.042,-.041,
     $-.020,-.028,-.013,-.014,-.011,-.010,
     $-.043,-.049,-.047,-.045,-.042,-.037,
     $-.003,-.013,-.010,-.003,-.007,-.008 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = 0.1*BETA
      M = INT(S)
      IF(M.EQ.-3) M=-2
      IF(M.GE.3) M=2
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      DLDA = V + (W-V)*ABS(DB)
C
      RETURN
      END
C
C
      FUNCTION DLDR(ALPHA,BETA)
C
C  This function computes the rolling moment due to rudder deflection
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      REAL A(-2:9,-3:3)
C
!       DATA A/  .005,.017,.014,.010,-.005,.009,
!      $.019,.005,.000,-.005,-.011,.008,
!      $.007,.016,.014,.014,.013,.009,
!      $.012,.005,.000,.004,.009,.007,
!      $.013,.013,.011,.012,.011,.009,
!      $.008,.005,.000,.005,.003,.005,
!      $.018,.015,.015,.014,.014,.014,
!      $.014,.015,.013,.011,.006,.001,
!      $.015,.014,.013,.013,.012,.011,
!      $.011,.010,.008,.008,.007,.003,
!      $.021,.011,.010,.011,.010,.009,
!      $.008,.010,.006,.005,.000,.001,
!      $.023,.010,.011,.011,.011,.010,
!      $.008,.010,.006,.014,.020,.000 /
C

      DATA A/ .005, .017, .014, .010, -.005, .009, .019,
     & .005, -.000, -.005, -.011, .008,
     & .007, .016, .014, .014, .013, .009, .012,
     & .005, .000, .004, .009, .007,
     & .013, .013, .011, .012, .011, .009, .008,
     & .005, -.002, .005, .003, .005,
     & .018, .015, .015, .014, .014, .014, .014,
     & .015, .013, .011, .006, .001,
     & .015, .014, .013, .013, .012, .011, .011,
     & .010, .008, .008, .007, .003,
     & .021, .011, .010, .011, .010, .009, .008,
     & .010, .006, .005, .000, .001,
     & .023, .010, .011, .011, .011, .010, .008,
     & .010, .006, .014, .020, .000/
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = 0.1*BETA
      M = INT(S)
      IF(M.EQ.-3) M=-2
      IF(M.GE.3) M=2
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      DLDR = V + (W-V)*ABS(DB)
C
      RETURN
      END
C
C
      FUNCTION DNDA(ALPHA,BETA)
C
C  This function computes the yawing moment due to aileron deflection
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      REAL A(-2:9,-3:3)
C
!       DATA A/  .001,-.027,-.017,-.013,-.012,-.016,
!      $.001,.017,.011,.017,.008,.016,
!      $.002,-.014,-.016,-.016,-.014,-.019,
!      $-.021,.002,.012,.016,.015,.011,
!      $-.006,-.008,-.006,-.006,-.005,-.008,
!      $-.005,.007,.004,.007,.006,.006,
!      $-.011,-.011,-.010,-.009,-.008,-.006,
!      $.000,.004,.007,.010,.004,.010,
!      $-.015,-.015,-.014,-.012,-.011,-.008,
!      $-.002,.002,.006,.012,.011,.011,
!      $-.024,-.010,-.004,-.002,-.001,.003,
!      $.014,.006,-.001,.004,.004,.006,
!      $-.022,.002,-.003,-.005,-.003,-.001,
!      $-.009,-.009,-.001,.003,-.002,.001 /
C
      DATA A/  .001,-.027,-.017,-.013,-.012,-.016,
     $.001,.017,.011,.017,.008,.016,
     $.002,-.014,-.016,-.016,-.014,-.019,
     $-.021,.002,.012,.015,.015,.011,
     $-.006,-.008,-.006,-.006,-.005,-.008,
     $-.005,.007,.004,.007,.006,.006,
     $-.011,-.011,-.010,-.009,-.008,-.006,
     $.000,.004,.007,.010,.004,.010,
     $-.015,-.015,-.014,-.012,-.011,-.008,
     $-.002,.002,.006,.012,.011,.011,
     $-.024,-.010,-.004,-.002,-.001,.003,
     $.014,.006,-.001,.004,.004,.006,
     $-.022,.002,-.003,-.005,-.003,-.001,
     $-.009,-.009,-.001,.003,-.002,.001 /

      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = 0.1*BETA
      M = INT(S)
      IF(M.EQ.-3) M=-2
      IF(M.GE.3) M=2
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      DNDA = V + (W-V)*ABS(DB)
C
      RETURN
      END
C
C
      FUNCTION DNDR(ALPHA,BETA)
C
C  This function computes the yawing moment due to rudder deflection
C  for the F-16 aerodynamic model.  
C
C  ALPHA = ANGLE OF ATTACK, DEG  ( -10 <= ALPHA <= 45 )
C  BETA = SIDESLIP ANGLE, DEG  ( -30 <= BETA <= 30 )
C
C
      REAL A(-2:9,-3:3)
C
      DATA A/  -.018,-.052,-.052,-.052,-.054,-.049,
     $-.059,-.051,-.030,-.037,-.026,-.013,
     $-.028,-.051,-.043,-.046,-.045,-.049,
     $-.057,-.052,-.030,-.033,-.030,-.008,
     $-.037,-.041,-.038,-.040,-.040,-.038,
     $-.037,-.030,-.027,-.024,-.019,-.013,
     $-.048,-.045,-.045,-.045,-.044,-.045,
     $-.047,-.048,-.049,-.045,-.033,-.016,
     $-.043,-.044,-.041,-.041,-.040,-.038,
     $-.034,-.035,-.035,-.029,-.022,-.009,
     $-.052,-.034,-.036,-.036,-.035,-.028,
     $-.024,-.023,-.020,-.016,-.010,-.014,
     $-.062,-.034,-.027,-.028,-.027,-.027,
     $-.023,-.023,-.019,-.009,-.025,-.010 /
C
      S = 0.2*ALPHA
      K = INT(S)
      IF(K.LE.-2) K=-1
      IF(K.GE.9) K=8
      DA = S - FLOAT(K)
      L = K + INT(SIGN(1.1,DA))
      S = 0.1*BETA
      M = INT(S)
      IF(M.EQ.-3) M=-2
      IF(M.GE.3) M=2
      DB = S - FLOAT(M)
      N = M + INT(SIGN(1.1,DB))
      V = A(K,M) + ABS(DA)*(A(L,M)-A(K,M))
      W = A(K,N) + ABS(DA)*(A(L,N)-A(K,N))
      DNDR = V + (W-V)*ABS(DB)
C
      RETURN
      END
C