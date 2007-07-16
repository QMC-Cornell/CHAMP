        PROGRAM MSTVL0

C       =====================================================
C       Purpose: This program computes modified Struve
C                function L0(x) using subroutine STVL0
C       Input :  x   --- Argument of L0(x) ( x ò 0 )
C       Output:  SL0 --- L0(x)
C       Example:
C                   x        L0(x)
C               ------------------------
C                  0.0   .00000000D+00
C                  5.0   .27105917D+02
C                 10.0   .28156522D+04
C                 15.0   .33964933D+06
C                 20.0   .43558283D+08
C                 30.0   .78167230D+12
C                 40.0   .14894775D+17
C                 50.0   .29325538D+21
C       =====================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x        L0(x)'
        WRITE(*,*)'-----------------------'
        CALL STVL0(X,SL0)
        WRITE(*,10)X,SL0
10      FORMAT(1X,F5.1,D16.8)
        END


        SUBROUTINE STVL0(X,SL0)

C       ================================================
C       Purpose: Compute modified Struve function L0(x)
C       Input :  x   --- Argument of L0(x) ( x ò 0 )
C       Output:  SL0 --- L0(x)
C       ================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        S=1.0D0
        R=1.0D0
        IF (X.LE.20.0D0) THEN
           A0=2.0D0*X/PI
           DO 10 K=1,60
              R=R*(X/(2.0D0*K+1.0D0))**2
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 15
10         CONTINUE
15         SL0=A0*S
        ELSE
           KM=INT(.5*(X+1.0))
           IF (X.GE.50.0) KM=25
           DO 20 K=1,KM
              R=R*((2.0D0*K-1.0D0)/X)**2
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 25
20         CONTINUE
25         A1=DEXP(X)/DSQRT(2.0D0*PI*X)
           R=1.0D0
           BI0=1.0D0
           DO 30 K=1,16
              R=0.125D0*R*(2.0D0*K-1.0D0)**2/(K*X)
              BI0=BI0+R
              IF (DABS(R/BI0).LT.1.0D-12) GO TO 35
30         CONTINUE
35         BI0=A1*BI0
           SL0=-2.0D0/(PI*X)*S+BI0
        ENDIF
        RETURN
        END
