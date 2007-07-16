        PROGRAM MLPNI

C       ========================================================
C       Purpose: This program computes the Legendre polynomials
C                Pn(x), Pn'(x) and the integral of Pn(t) from 0
C                to x using subroutine LPNI
C       Input :  x --- Argument of Pn(x)
C                n --- Degree of Pn(x) ( n = 0,1,... )
C       Output:  PN(n) --- Pn(x)
C                PD(n) --- Pn'(x)
C                PL(n) --- Integral of Pn(t) from 0 to x
C       Example: x = 0.50
C                n       Pn(x)         Pn'(x)        Pn(t)dt
C               ---------------------------------------------
C                0    1.00000000     .00000000     .50000000
C                1     .50000000    1.00000000     .12500000
C                2    -.12500000    1.50000000    -.18750000
C                3    -.43750000     .37500000    -.14843750
C                4    -.28906250   -1.56250000     .05859375
C                5     .08984375   -2.22656250     .11816406
C       ========================================================

        DOUBLE PRECISION PN,PD,PL,X
        DIMENSION PN(0:100),PD(0:100),PL(0:100)
        WRITE(*,*)'  Please enter Nmax and x'
        READ(*,*)N,X
        WRITE(*,30)X
        WRITE(*,*)
        WRITE(*,*)'  n        Pn(x)          Pn''(x)         Pn(t)dt'
        WRITE(*,*)' ---------------------------------------------------'
        CALL LPNI(N,X,PN,PD,PL)
        DO 10 K=0,N
10         WRITE(*,20)K,PN(K),PD(K),PL(K)
20      FORMAT(1X,I3,3E16.8)
30      FORMAT(3X,'x =',F5.2)
        END


        SUBROUTINE LPNI(N,X,PN,PD,PL)

C       =====================================================
C       Purpose: Compute Legendre polynomials Pn(x), Pn'(x)
C                and the integral of Pn(t) from 0 to x
C       Input :  x --- Argument of Pn(x)
C                n --- Degree of Pn(x) ( n = 0,1,... )
C       Output:  PN(n) --- Pn(x)
C                PD(n) --- Pn'(x)
C                PL(n) --- Integral of Pn(t) from 0 to x
C       =====================================================

        IMPLICIT DOUBLE PRECISION (P,R,X)
        DIMENSION PN(0:N),PD(0:N),PL(0:N)
        PN(0)=1.0D0
        PN(1)=X
        PD(0)=0.0D0
        PD(1)=1.0D0
        PL(0)=X
        PL(1)=0.5D0*X*X
        P0=1.0D0
        P1=X
        DO 15 K=2,N
           PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0
           PN(K)=PF
           IF (DABS(X).EQ.1.0D0) THEN
              PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
           ENDIF
           PL(K)=(X*PN(K)-PN(K-1))/(K+1.0D0)
           P0=P1
           P1=PF
           IF (K.EQ.2*INT(K/2)) GO TO 15
           R=1.0D0/(K+1.0D0)
           N1=(K-1)/2
           DO 10 J=1,N1
10            R=(0.5D0/J-1.0D0)*R
           PL(K)=PL(K)+R
15      CONTINUE
        RETURN
        END
