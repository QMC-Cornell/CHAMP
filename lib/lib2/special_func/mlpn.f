        PROGRAM MLPN

C       ========================================================
C       Purpose: This program computes the Legendre polynomials
C                Pn(x) and their derivatives Pn'(x) using
C                subroutine LPN
C       Input :  x --- Argument of Pn(x)
C                n --- Degree of Pn(x) ( n = 0,1,...)
C       Output:  PN(n) --- Pn(x)
C                PD(n) --- Pn'(x)
C       Example:    x = 0.5
C                  n          Pn(x)            Pn'(x)
C                ---------------------------------------
C                  0       1.00000000        .00000000
C                  1        .50000000       1.00000000
C                  2       -.12500000       1.50000000
C                  3       -.43750000        .37500000
C                  4       -.28906250      -1.56250000
C                  5        .08984375      -2.22656250
C       ========================================================

        DOUBLE PRECISION PN,PD,X
        DIMENSION PN(0:100),PD(0:100)
        WRITE(*,*)'  Please enter Nmax and x '
        READ(*,*)N,X
        WRITE(*,30)X
        WRITE(*,*)
        CALL LPN(N,X,PN,PD)
        WRITE(*,*)'  n         Pn(x)           Pn''(x)'
        WRITE(*,*)'---------------------------------------'
        DO 10 K=0,N
10         WRITE(*,20)K,PN(K),PD(K)
20      FORMAT(1X,I3,2E17.8)
30      FORMAT(3X,'x =',F5.1)
        END


        SUBROUTINE LPN(N,X,PN,PD)

C       ===============================================
C       Purpose: Compute Legendre polynomials Pn(x)
C                and their derivatives Pn'(x)
C       Input :  x --- Argument of Pn(x)
C                n --- Degree of Pn(x) ( n = 0,1,...)
C       Output:  PN(n) --- Pn(x)
C                PD(n) --- Pn'(x)
C       ===============================================

        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PN(0:N),PD(0:N)
        PN(0)=1.0D0
        PN(1)=X
        PD(0)=0.0D0
        PD(1)=1.0D0
        P0=1.0D0
        P1=X
        DO 10 K=2,N
           PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0
           PN(K)=PF
           IF (DABS(X).EQ.1.0D0) THEN
              PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
           ENDIF
           P0=P1
10         P1=PF
        RETURN
        END
