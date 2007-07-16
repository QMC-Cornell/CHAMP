        PROGRAM ME1XB

C       =========================================================
C       Purpose: This program computes the exponential integral
C                E1(x) using subroutine E1XB
C       Input :  x  --- Argument of E1(x)  ( x > 0 )
C       Output:  E1 --- E1(x)
C       Example:
C                  x          E1(x)
C                -------------------------
C                 0.0     .1000000000+301
C                 1.0     .2193839344E+00
C                 2.0     .4890051071E-01
C                 3.0     .1304838109E-01
C                 4.0     .3779352410E-02
C                 5.0     .1148295591E-02
C       =========================================================

        DOUBLE PRECISION E1,X
        WRITE(*,*)'Please enter x '
        READ(*,*) X
        WRITE(*,*)'   x          E1(x)'
        WRITE(*,*)' -------------------------'
        CALL E1XB(X,E1)
        WRITE(*,10)X,E1
10      FORMAT(1X,F5.1,E20.10)
        END


        SUBROUTINE E1XB(X,E1)

C       ============================================
C       Purpose: Compute exponential integral E1(x)
C       Input :  x  --- Argument of E1(x)
C       Output:  E1 --- E1(x)  ( x > 0 )
C       ============================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0) THEN
           E1=1.0D+300
        ELSE IF (X.LE.1.0) THEN
           E1=1.0D0
           R=1.0D0
           DO 10 K=1,25
              R=-R*K*X/(K+1.0D0)**2
              E1=E1+R
              IF (DABS(R).LE.DABS(E1)*1.0D-15) GO TO 15
10         CONTINUE
15         GA=0.5772156649015328D0
           E1=-GA-DLOG(X)+X*E1
        ELSE
           M=20+INT(80.0/X)
           T0=0.0D0
           DO 20 K=M,1,-1
              T0=K/(1.0D0+K/(X+T0))
20         CONTINUE
           T=1.0D0/(X+T0)
           E1=DEXP(-X)*T
        ENDIF
        RETURN
        END
