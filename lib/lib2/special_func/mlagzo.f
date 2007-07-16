        PROGRAM MLAGZO

C       ===========================================================
C       Purpose : This program computes the zeros of Laguerre
C                 polynomial Ln(x) in the interval [0,ì] and the
C                 corresponding weighting coefficients for Gauss-
C                 Laguerre integration using subroutine LAGZO
C       Input :   n    --- Order of the Laguerre polynomial
C                 X(n) --- Zeros of the Laguerre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(100),W(100)
        WRITE(*,*)'Please enter the order of Ln(x), n '
        READ(*,*)N
        WRITE(*,20)N
        CALL LAGZO(N,X,W)
        WRITE(*,*)'  Nodes and weights for Gauss-Lagurre integration'
        WRITE(*,*)
        WRITE(*,*)'  i             xi                      Wi'
        WRITE(*,*)' -----------------------------------------',
     &            '------------'
        DO 10 J=1,N
10         WRITE(*,30)J,X(J),W(J)
20      FORMAT(1X,'n =',I3)
30      FORMAT(1X,I3,3X,D22.13,3X,D22.13)
        END


        SUBROUTINE LAGZO(N,X,W)

C       =========================================================
C       Purpose : Compute the zeros of Laguerre polynomial Ln(x)
C                 in the interval [0,ì], and the corresponding
C                 weighting coefficients for Gauss-Laguerre
C                 integration
C       Input :   n    --- Order of the Laguerre polynomial
C                 X(n) --- Zeros of the Laguerre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),W(N)
        HN=1.0D0/N
        DO 35 NR=1,N
           IF (NR.EQ.1) Z=HN
           IF (NR.GT.1) Z=X(NR-1)+HN*NR**1.27
           IT=0
10         IT=IT+1
           Z0=Z
           P=1.0D0
           DO 15 I=1,NR-1
15            P=P*(Z-X(I))
           F0=1.0D0
           F1=1.0D0-Z
           DO 20 K=2,N
              PF=((2.0D0*K-1.0D0-Z)*F1-(K-1.0D0)*F0)/K
              PD=K/Z*(PF-F1)
              F0=F1
20            F1=PF
           FD=PF/P
           Q=0.0D0
           DO 30 I=1,NR-1
              WP=1.0D0
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 WP=WP*(Z-X(J))
25            CONTINUE
              Q=Q+WP
30         CONTINUE
           GD=(PD-Q*FD)/P
           Z=Z-FD/GD
           IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
           X(NR)=Z
           W(NR)=1.0D0/(Z*PD*PD)
35      CONTINUE
        RETURN
        END
