        PROGRAM MSTVH0

C       ====================================================
C       Purpose: This program computes Struve function
C                H0(x) using subroutine STVH0
C       Input :  x   --- Argument of H0(x) ( x ò 0 )
C       Output:  SH0 --- H0(x)
C       Example:
C                   x          H0(x)
C                ----------------------
C                  0.0       .00000000
C                  5.0      -.18521682
C                 10.0       .11874368
C                 15.0       .24772383
C                 20.0       .09439370
C                 25.0      -.10182519
C       ====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x          H0(x)'
        WRITE(*,*)'----------------------'
        CALL STVH0(X,SH0)
        WRITE(*,10)X,SH0
10      FORMAT(1X,F5.1,E16.8)
        END


        SUBROUTINE STVH0(X,SH0)

C       =============================================
C       Purpose: Compute Struve function H0(x)
C       Input :  x   --- Argument of H0(x) ( x ò 0 )
C       Output:  SH0 --- H0(x)
C       =============================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        S=1.0D0
        R=1.0D0
        IF (X.LE.20.0D0) THEN
           A0=2.0*X/PI
           DO 10 K=1,60
              R=-R*X/(2.0D0*K+1.0D0)*X/(2.0D0*K+1.0D0)
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         SH0=A0*S
        ELSE
           KM=INT(.5*(X+1.0))
           IF (X.GE.50.0) KM=25
           DO 20 K=1,KM
              R=-R*((2.0D0*K-1.0D0)/X)**2
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 25
20         CONTINUE
25         T=4.0D0/X
           T2=T*T
           P0=((((-.37043D-5*T2+.173565D-4)*T2-.487613D-4)
     &        *T2+.17343D-3)*T2-.1753062D-2)*T2+.3989422793D0
           Q0=T*(((((.32312D-5*T2-.142078D-4)*T2+.342468D-4)*
     &        T2-.869791D-4)*T2+.4564324D-3)*T2-.0124669441D0)
           TA0=X-.25D0*PI
           BY0=2.0D0/DSQRT(X)*(P0*DSIN(TA0)+Q0*DCOS(TA0))
           SH0=2.0D0/(PI*X)*S+BY0
        ENDIF
        RETURN
        END
