        PROGRAM MBERNOB

C       ===========================================================
C       Purpose: This program computes Bernoulli number Bn using
C                subroutine BERNOB
C       Example: Compute Bernouli number Bn for n = 0,1,...,10
C                Computed results:

C                   n            Bn
C                 --------------------------
C                   0     .100000000000D+01
C                   1    -.500000000000D+00
C                   2     .166666666667D+00
C                   4    -.333333333333D-01
C                   6     .238095238095D-01
C                   8    -.333333333333D-01
C                  10     .757575757576D-01
C       ===========================================================

        DOUBLE PRECISION B
        DIMENSION B(0:200)
        WRITE(*,*)'  Please enter Nmax'
        READ(*,*)N
        CALL BERNOB(N,B)
        WRITE(*,*)'   n            Bn'
        WRITE(*,*)' --------------------------'
        WRITE(*,20)0,B(0)
        WRITE(*,20)1,B(1)
        DO 10 K=2,N,2
10         WRITE(*,20)K,B(K)
20      FORMAT(2X,I3,D22.12)
        END


        SUBROUTINE BERNOB(N,BN)

C       ======================================
C       Purpose: Compute Bernoulli number Bn
C       Input :  n --- Serial number
C       Output:  BN(n) --- Bn
C       ======================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BN(0:N)
        TPI=6.283185307179586D0
        BN(0)=1.0D0
        BN(1)=-0.5D0
        BN(2)=1.0D0/6.0D0
        R1=(2.0D0/TPI)**2
        DO 20 M=4,N,2
           R1=-R1*(M-1)*M/(TPI*TPI)
           R2=1.0D0
           DO 10 K=2,10000
              S=(1.0D0/K)**M
              R2=R2+S
              IF (S.LT.1.0D-15) GOTO 20
10         CONTINUE
20         BN(M)=R1*R2
        RETURN
        END
