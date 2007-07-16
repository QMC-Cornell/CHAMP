        PROGRAM MITTIKB

C       ============================================================
C       Purpose: This program computes the integral of [I0(t)-1]/t
C                with respect to t from 0 to x and K0(t)/t with
C                respect to t from x to ì using subroutine ITTIKB
C       Input :  x   --- Upper limit of the integral
C       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
C                TTK --- Integration of K0(t)/t from x to ì
C       Example:
C                   x     [1-I0(t)]/tdt      K0(t)/tdt
C                ---------------------------------------
C                  5.0     .710478D+01     .586361D-03
C                 10.0     .340811D+03     .156293D-05
C                 15.0     .254373D+05     .598363D-08
C                 20.0     .236735D+07     .267906D-10
C                 25.0     .246534D+09     .131007D-12
C       ============================================================

        DOUBLE PRECISION X,TTI,TTK
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x     [1-I0(t)]/tdt      K0(t)/tdt'
        WRITE(*,*)'---------------------------------------'
        CALL ITTIKB(X,TTI,TTK)
        WRITE(*,10)X,TTI,TTK
10      FORMAT(1X,F5.1,2D16.6)
        END


        SUBROUTINE ITTIKB(X,TTI,TTK)

C       =========================================================
C       Purpose: Integrate [I0(t)-1]/t with respect to t from 0
C                to x, and K0(t)/t with respect to t from x to ì
C       Input :  x   --- Variable in the limits  ( x ò 0 )
C       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
C                TTK --- Integration of K0(t)/t from x to ì
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IF (X.EQ.0.0D0) THEN
           TTI=0.0D0
        ELSE IF (X.LE.5.0D0) THEN
           X1=X/5.0D0
           T=X1*X1
           TTI=(((((((.1263D-3*T+.96442D-3)*T+.968217D-2)*T
     &         +.06615507D0)*T+.33116853D0)*T+1.13027241D0)
     &         *T+2.44140746D0)*T+3.12499991D0)*T
        ELSE
           T=5.0D0/X
           TTI=(((((((((2.1945464D0*T-3.5195009D0)*T
     &         -11.9094395D0)*T+40.394734D0)*T-48.0524115D0)
     &         *T+28.1221478D0)*T-8.6556013D0)*T+1.4780044D0)
     &         *T-.0493843D0)*T+.1332055D0)*T+.3989314D0
           TTI=TTI*DEXP(X)/(DSQRT(X)*X)
        ENDIF
        IF (X.EQ.0.0D0) THEN
           TTK=1.0D+300
        ELSE IF (X.LE.2.0D0) THEN
           T1=X/2.0D0
           T=T1*T1
           TTK=(((((.77D-6*T+.1544D-4)*T+.48077D-3)*T
     &         +.925821D-2)*T+.10937537D0)*T+.74999993D0)*T
           E0=EL+DLOG(X/2.0D0)
           TTK=PI*PI/24.0D0+E0*(.5D0*E0+TTI)-TTK
        ELSE IF (X.LE.4.0D0) THEN
           T=2.0D0/X
           TTK=(((.06084D0*T-.280367D0)*T+.590944D0)*T
     &         -.850013D0)*T+1.234684D0
           TTK=TTK*DEXP(-X)/(DSQRT(X)*X)
        ELSE
           T=4.0D0/X
           TTK=(((((.02724D0*T-.1110396D0)*T+.2060126D0)*T
     &         -.2621446D0)*T+.3219184D0)*T-.5091339D0)*T
     &         +1.2533141D0
           TTK=TTK*DEXP(-X)/(DSQRT(X)*X)
        ENDIF
        RETURN
        END
