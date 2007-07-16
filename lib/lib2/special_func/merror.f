        PROGRAM MERROR

C       ===================================================
C       Purpose: This program computes the error function
C                erf(x) using subroutine ERROR
C       Input:   x   --- Argument of erf(x)
C       Output:  ERR --- erf(x)
C       Example:
C                  x         erf(x)
C                ---------------------
C                 1.0       .84270079
C                 2.0       .99532227
C                 3.0       .99997791
C                 4.0       .99999998
C                 5.0      1.00000000
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x         erf(x)'
        WRITE(*,*)'---------------------'
        CALL ERROR(X,ER)
        WRITE(*,10)X,ER
10      FORMAT(1X,F5.2,F15.8)
        END


        SUBROUTINE ERROR(X,ERR)

C       =========================================
C       Purpose: Compute error function erf(x)
C       Input:   x   --- Argument of erf(x)
C       Output:  ERR --- erf(x)
C       =========================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        X2=X*X
        IF (DABS(X).LT.3.5D0) THEN
           ER=1.0D0
           R=1.0D0
           DO 10 K=1,50
              R=R*X2/(K+0.5D0)
              ER=ER+R
              IF (DABS(R).LE.DABS(ER)*EPS) GO TO 15
10         CONTINUE
15         C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
           ERR=C0*ER
        ELSE
           ER=1.0D0
           R=1.0D0
           DO 20 K=1,12
              R=-R*(K-0.5D0)/X2
20            ER=ER+R
           C0=DEXP(-X2)/(DABS(X)*DSQRT(PI))
           ERR=1.0D0-C0*ER
           IF (X.LT.0.0) ERR=-ERR
        ENDIF
        RETURN
        END
