        PROGRAM MITSL0

C       ===========================================================
C       Purpose: This program evaluates the integral of modified
C                Struve function L0(t) with respect to t from 0
C                to x using subroutine ITSL0
C       Input :  x   --- Upper limit  ( x ò 0 )
C       Output:  TL0 --- Integration of L0(t) from 0 to x
C       Example:
C                      x        L0(t)dt
C                   -----------------------
C                     0.0    .0000000D+00
C                     5.0    .3003079D+02
C                    10.0    .2990773D+04
C                    15.0    .3526179D+06
C                    20.0    .4475860D+08
C                    30.0    .7955389D+12
C                    40.0    .1508972D+17
C                    50.0    .2962966D+21
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x        L0(t)dt'
        WRITE(*,*)'-----------------------'
        CALL ITSL0(X,TL0)
        WRITE(*,10)X,TL0
10      FORMAT(1X,F5.1,D16.7)
        END


        SUBROUTINE ITSL0(X,TL0)

C       ===========================================================
C       Purpose: Evaluate the integral of modified Struve function
C                L0(t) with respect to t from 0 to x
C       Input :  x   --- Upper limit  ( x ò 0 )
C       Output:  TL0 --- Integration of L0(t) from 0 to x
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(18)
        PI=3.141592653589793D0
        R=1.0D0
        IF (X.LE.20.0) THEN
           S=0.5D0
           DO 10 K=1,100
              RD=1.0D0
              IF (K.EQ.1) RD=0.5D0
              R=R*RD*K/(K+1.0D0)*(X/(2.0D0*K+1.0D0))**2
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 15
10         CONTINUE
15         TL0=2.0D0/PI*X*X*S
        ELSE
           S=1.0D0
           DO 20 K=1,10
              R=R*K/(K+1.0D0)*((2.0D0*K+1.0D0)/X)**2
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 25
20            CONTINUE
25         EL=.57721566490153D0
           S0=-S/(PI*X*X)+2.0D0/PI*(DLOG(2.0D0*X)+EL)
           A0=1.0D0
           A1=5.0D0/8.0D0
           A(1)=A1
           DO 30 K=1,10
              AF=((1.5D0*(K+.50D0)*(K+5.0D0/6.0D0)*A1-.5D0*
     &            (K+.5D0)**2*(K-.5D0)*A0))/(K+1.0D0)
              A(K+1)=AF
              A0=A1
30            A1=AF
           TI=1.0D0
           R=1.0D0
           DO 35 K=1,11
              R=R/X
35            TI=TI+A(K)*R
           TL0=TI/DSQRT(2*PI*X)*DEXP(X)+S0
        ENDIF
        RETURN
        END
