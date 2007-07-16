        PROGRAM MLEGZO

C       ============================================================
C       Purpose : This program computes the zeros of Legendre
C                 polynomial Pn(x) in the interval [-1,1] and the
C                 corresponding weighting coefficients for Gauss-
C                 Legendre integration using subroutine LEGZO
C       Input :   n    --- Order of the Legendre polynomial
C       Output:   X(n) --- Zeros of the Legendre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(120),W(120)
        WRITE(*,*)'Please enter the order of Pn(x), n '
        READ(*,*)N
        WRITE(*,15)N
        CALL LEGZO(N,X,W)
        WRITE(*,*)'  Nodes and weights for Gauss-Legendre integration'
        WRITE(*,*)
        WRITE(*,*)'  i              xi                   Wi'
        WRITE(*,*)' ------------------------------------------------'
        DO 10 I=1,N
10         WRITE(*,20)I,X(I),W(I)
15      FORMAT(1X,'n =',I3)
20      FORMAT(1X,I3,1X,F22.13,D22.13)
        END


        SUBROUTINE LEGZO(N,X,W)

C       =========================================================
C       Purpose : Compute the zeros of Legendre polynomial Pn(x)
C                 in the interval [-1,1], and the corresponding
C                 weighting coefficients for Gauss-Legendre
C                 integration
C       Input :   n    --- Order of the Legendre polynomial
C       Output:   X(n) --- Zeros of the Legendre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),W(N)
        N0=(N+1)/2
        DO 45 NR=1,N0
           Z=DCOS(3.1415926D0*(NR-0.25D0)/N)
10         Z0=Z
           P=1.0D0
           DO 15 I=1,NR-1
15            P=P*(Z-X(I))
           F0=1.0D0
           IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0.0D0
           F1=Z
           DO 20 K=2,N
              PF=(2.0D0-1.0D0/K)*Z*F1-(1.0D0-1.0D0/K)*F0
              PD=K*(F1-Z*PF)/(1.0D0-Z*Z)
              F0=F1
20            F1=PF
           IF (Z.EQ.0.0) GO TO 40
           FD=PF/P
           Q=0.0D0
           DO 35 I=1,NR-1
              WP=1.0D0
              DO 30 J=1,NR-1
                 IF (J.NE.I) WP=WP*(Z-X(J))
30            CONTINUE
35            Q=Q+WP
           GD=(PD-Q*FD)/P
           Z=Z-FD/GD
           IF (DABS(Z-Z0).GT.DABS(Z)*1.0D-15) GO TO 10
40         X(NR)=Z
           X(N+1-NR)=-Z
           W(NR)=2.0D0/((1.0D0-Z*Z)*PD*PD)
45         W(N+1-NR)=W(NR)
        RETURN
        END
