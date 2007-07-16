        PROGRAM MHERZO

C       ===========================================================
C       Purpose : This program computes the zeros of Hermite
C                 polynomial Ln(x) in the interval [-ì,ì] and the
C                 corresponding weighting coefficients for Gauss-
C                 Hermite integration using subroutine HERZO
C       Input :   n    --- Order of the Hermite polynomial
C                 X(n) --- Zeros of the Hermite polynomial
C                 W(n) --- Corresponding weighting coefficients
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(100),W(100)
        WRITE(*,*)'Please enter the order of Hn(x), n '
        READ(*,*)N
        WRITE(*,20)N
        CALL HERZO(N,X,W)
        WRITE(*,*)'  Nodes and weights for Gauss-Hermite integration'
        WRITE(*,*)
        WRITE(*,*)'  i             xi                      Wi'
        WRITE(*,*)' -----------------------------------------',
     &            '------------'
        DO 10 J=1,N
10         WRITE(*,30)J,X(J),W(J)
20      FORMAT(1X,'n =',I3)
30      FORMAT(1X,I3,3X,D22.13,3X,D22.13)
        END


        SUBROUTINE HERZO(N,X,W)

C       ========================================================
C       Purpose : Compute the zeros of Hermite polynomial Ln(x)
C                 in the interval [-ì,ì], and the corresponding
C                 weighting coefficients for Gauss-Hermite
C                 integration
C       Input :   n    --- Order of the Hermite polynomial
C                 X(n) --- Zeros of the Hermite polynomial
C                 W(n) --- Corresponding weighting coefficients
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),W(N)
        HN=1.0D0/N
        ZL=-1.1611D0+1.46D0*N**0.5
        DO 40 NR=1,N/2
           IF (NR.EQ.1) Z=ZL
           IF (NR.NE.1) Z=Z-HN*(N/2+1-NR)
           IT=0
10         IT=IT+1
           Z0=Z
           F0=1.0D0
           F1=2.0D0*Z
           DO 15 K=2,N
              HF=2.0D0*Z*F1-2.0D0*(K-1.0D0)*F0
              HD=2.0D0*K*F1
              F0=F1
15            F1=HF
           P=1.0D0
           DO 20 I=1,NR-1
20            P=P*(Z-X(I))
           FD=HF/P
           Q=0.0D0
           DO 30 I=1,NR-1
              WP=1.0D0
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 WP=WP*(Z-X(J))
25            CONTINUE
30            Q=Q+WP
           GD=(HD-Q*FD)/P
           Z=Z-FD/GD
           IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
           X(NR)=Z
           X(N+1-NR)=-Z
           R=1.0D0
           DO 35 K=1,N
35            R=2.0D0*R*K
           W(NR)=3.544907701811D0*R/(HD*HD)
40         W(N+1-NR)=W(NR)
        IF (N.NE.2*INT(N/2)) THEN
           R1=1.0D0
           R2=1.0D0
           DO 45 J=1,N
              R1=2.0D0*R1*J
              IF (J.GE.(N+1)/2) R2=R2*J
45         CONTINUE
           W(N/2+1)=0.88622692545276D0*R1/(R2*R2)
           X(N/2+1)=0.0D0
        ENDIF
        RETURN
        END
