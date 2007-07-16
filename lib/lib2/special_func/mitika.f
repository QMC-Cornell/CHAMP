        PROGRAM MITIKA

C       ============================================================
C       Purpose: This program evaluates the integral of modified
C                Bessel functions I0(t) and K0(t) with respect to t
C                from 0 to x using subroutine ITIKA
C       Input :  x  --- Upper limit of the integral  ( x ò 0 )
C       Output:  TI --- Integration of I0(t) from 0 to x
C                TK --- Integration of K0(t) from 0 to x
C       Example:
C                    x         I0(t)dt         K0(t)dt
C                 --------------------------------------
C                   5.0    .31848668D+02     1.56738739
C                  10.0    .29930445D+04     1.57077931
C                  15.0    .35262048D+06     1.57079623
C                  20.0    .44758593D+08     1.57079633
C                  25.0    .58991731D+10     1.57079633
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x         I0(t)dt         K0(t)dt'
        WRITE(*,*)' --------------------------------------'
        CALL ITIKA(X,TI,TK)
        WRITE(*,10)X,TI,TK
10      FORMAT(1X,F5.1,D17.8,F15.8)
        END


        SUBROUTINE ITIKA(X,TI,TK)

C       =======================================================
C       Purpose: Integrate modified Bessel functions I0(t) and
C                K0(t) with respect to t from 0 to x
C       Input :  x  --- Upper limit of the integral  ( x ò 0 )
C       Output:  TI --- Integration of I0(t) from 0 to x
C                TK --- Integration of K0(t) from 0 to x
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        DATA A/.625D0,1.0078125D0,
     &       2.5927734375D0,9.1868591308594D0,
     &       4.1567974090576D+1,2.2919635891914D+2,
     &       1.491504060477D+3,1.1192354495579D+4,
     &       9.515939374212D+4,9.0412425769041D+5/
        IF (X.EQ.0.0D0) THEN
           TI=0.0D0
           TK=0.0D0
           RETURN
        ELSE IF (X.LT.20.0D0) THEN
           X2=X*X
           TI=1.0D0
           R=1.0D0
           DO 10 K=1,50
              R=.25D0*R*(2*K-1.0D0)/(2*K+1.0D0)/(K*K)*X2
              TI=TI+R
              IF (DABS(R/TI).LT.1.0D-12) GO TO 15
10         CONTINUE
15         TI=TI*X
        ELSE
           TI=1.0D0
           R=1.0D0
           DO 20 K=1,10
              R=R/X
20            TI=TI+A(K)*R
           RC1=1.0D0/DSQRT(2.0D0*PI*X)
           TI=RC1*DEXP(X)*TI
        ENDIF
        IF (X.LT.12.0D0) THEN
           E0=EL+DLOG(X/2.0D0)
           B1=1.0D0-E0
           B2=0.0D0
           RS=0.0D0
           R=1.0D0
           DO 25 K=1,50
              R=.25D0*R*(2*K-1.0D0)/(2*K+1.0D0)/(K*K)*X2
              B1=B1+R*(1.0D0/(2*K+1)-E0)
              RS=RS+1.0D0/K
              B2=B2+R*RS
              TK=B1+B2
              IF (DABS((TK-TW)/TK).LT.1.0D-12) GO TO 30
25            TW=TK
30         TK=TK*X
        ELSE
           TK=1.0D0
           R=1.0D0
           DO 35 K=1,10
              R=-R/X
35            TK=TK+A(K)*R
           RC2=DSQRT(PI/(2.0D0*X))
           TK=PI/2.0D0-RC2*TK*DEXP(-X)
        ENDIF
        RETURN
        END
