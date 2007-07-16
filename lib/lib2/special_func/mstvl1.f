        PROGRAM MSTVL1

C       =====================================================
C       Purpose: This program computes the modified Struve
C                function L1(x) using subroutine STVL1
C       Input :  x   --- Argument of L1(x) ( x ò 0 )
C       Output:  SL1 --- L1(x)
C       Example:
C                     x        L1(x)
C                 -----------------------
C                   0.0   .00000000D+00
C                   5.0   .23728216D+02
C                  10.0   .26703583D+04
C                  15.0   .32812429D+06
C                  20.0   .42454973D+08
C                  30.0   .76853204D+12
C                  40.0   .14707396D+17
C                  50.0   .29030786D+21
C       =====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x        L1(x)'
        WRITE(*,*)'-----------------------'
        CALL STVL1(X,SL1)
        WRITE(*,10)X,SL1
10      FORMAT(1X,F5.1,D16.8)
        END


        SUBROUTINE STVL1(X,SL1)

C       ================================================
C       Purpose: Compute modified Struve function L1(x)
C       Input :  x   --- Argument of L1(x) ( x ò 0 )
C       Output:  SL1 --- L1(x)
C       ================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        R=1.0D0
        IF (X.LE.20.0D0) THEN
           S=0.0D0
           DO 10 K=1,60
              R=R*X*X/(4.0D0*K*K-1.0D0)
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 15
10         CONTINUE
15         SL1=2.0D0/PI*S
        ELSE
           S=1.0D0
           KM=INT(.50*X)
           IF (X.GT.50) KM=25
           DO 20 K=1,KM
              R=R*(2.0D0*K+3.0D0)*(2.0D0*K+1.0D0)/(X*X)
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 25
20            CONTINUE
25         SL1=2.0D0/PI*(-1.0D0+1.0D0/(X*X)+3.0D0*S/X**4)
           A1=DEXP(X)/DSQRT(2.0D0*PI*X)
           R=1.0D0
           BI1=1.0D0
           DO 30 K=1,16
              R=-0.125D0*R*(4.0D0-(2.0D0*K-1.0D0)**2)/(K*X)
              BI1=BI1+R
              IF (DABS(R/BI1).LT.1.0D-12) GO TO 35
30         CONTINUE
35         SL1=SL1+A1*BI1
        ENDIF
        RETURN
        END
