!********************************************************************
        double precision   FUNCTION erf_spline(x)
!********************************************************************
!
!     Written by Jesper Kielberg Pedersen, Mar. 2003
!
!     Purpose : Calculate value of error-function in point X
!               using splines
!
!               Based on f90-code from :
!               Naval Surface Warfare Center Mathematical Library
!               (http://www.math.iastate.edu/burkardt/f_src/nswc/nswc.html)
!*****************************************************************************
      implicit none
      double precision x, C, AX,HALF,T, ONE,TOP, BOT, FN_VAL, ZERO,X2
      PARAMETER (C = .564189583547756D0, ONE = 1.0D0, HALF = 0.5D0, ZERO = 0.0D0)
      DOUBLE PRECISION A(5), B(3), P(8), Q(8), R(5), S(4)

      SAVE A,B,P,Q,R,S
      DATA A / .771058495001320D-04, -.133733772997339D-02, &
     &         .323076579225834D-01,  .479137145607681D-01, &
     &         .128379167095513D+00 /
      DATA B / .301048631703895D-02,  .538971687740286D-01, &
     &         .375795757275549D+00 /
      DATA P / -1.36864857382717D-07, 5.64195517478974D-01, &
     &          7.21175825088309D+00, 4.31622272220567D+01, &
     &          1.52989285046940D+02, 3.39320816734344D+02, &
     &          4.51918953711873D+02, 3.00459261020162D+02 /
      DATA Q /  1.00000000000000D+00, 1.27827273196294D+01, &
     &          7.70001529352295D+01, 2.77585444743988D+02, &
     &          6.38980264465631D+02, 9.31354094850610D+02, &
     &          7.90950925327898D+02, 3.00459260956983D+02 /
      DATA R /  2.10144126479064D+00, 2.62370141675169D+01, &
     &          2.13688200555087D+01, 4.65807828718470D+00, &
     &          2.82094791773523D-01 /
      DATA S /  9.41537750555460D+01, 1.87114811799590D+02, &
     &          9.90191814623914D+01, 1.80124575948747D+01 /

      AX = ABS(X)

      IF (AX .LE. HALF) THEN
        T = X*X
        TOP = ((((A(1)*T + A(2))*T + A(3))*T + A(4))*T + A(5)) + ONE
        BOT = ((B(1)*T + B(2))*T + B(3))*T + ONE
        FN_VAL = X*(TOP/BOT)
      ELSE IF (AX .LE. 4.0D0) THEN
        TOP = ((((((P(1)*AX + P(2))*AX + P(3))*AX + P(4))*AX + P(5))*AX &
     &        + P(6))*AX + P(7))*AX + P(8)
        BOT = ((((((Q(1)*AX + Q(2))*AX + Q(3))*AX + Q(4))*AX + Q(5))*AX &
     &        + Q(6))*AX + Q(7))*AX + Q(8)
        FN_VAL = HALF + (HALF - DEXP(-X*X)*TOP/BOT)
        IF (X .LT. ZERO) FN_VAL = -FN_VAL
      ELSE IF (AX .LT. 5.8D0) THEN
        X2 = X*X
        T = ONE / X2
        TOP = (((R(1)*T + R(2))*T + R(3))*T + R(4))*T + R(5)
        BOT = (((S(1)*T + S(2))*T + S(3))*T + S(4))*T + ONE
        FN_VAL = (C - TOP/(X2*BOT)) / AX
        FN_VAL = HALF + (HALF - DEXP(-X2)*FN_VAL)
        IF (X .LT. ZERO) FN_VAL = -FN_VAL
      ELSE
        FN_VAL = SIGN(ONE, X)
      END IF
!
      erf_spline = FN_VAL
      RETURN
      END
!*****************************************************************************
