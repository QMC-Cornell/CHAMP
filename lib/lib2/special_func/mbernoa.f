        PROGRAM MBERNOA

C       ===========================================================
C       Purpose: This program computes Bernoulli number Bn using
C                subroutine BERNOA
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
        CALL BERNOA(N,B)
        WRITE(*,*)'   n            Bn'
        WRITE(*,*)' --------------------------'
        WRITE(*,20)0,B(0)
        WRITE(*,20)1,B(1)
        DO 10 K=2,N,2
10         WRITE(*,20)K,B(K)
20      FORMAT(2X,I3,D22.12)
        END


        SUBROUTINE BERNOA(N,BN)

C       ======================================
C       Purpose: Compute Bernoulli number Bn
C       Input :  n --- Serial number
C       Output:  BN(n) --- Bn
C       ======================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BN(0:N)
        BN(0)=1.0D0
        BN(1)=-0.5D0
        DO 30 M=2,N
           S=-(1.0D0/(M+1.0D0)-0.5D0)
           DO 20 K=2,M-1
              R=1.0D0
              DO 10 J=2,K
10               R=R*(J+M-K)/J
20            S=S-R*BN(K)
30         BN(M)=S
        DO 40 M=3,N,2
40         BN(M)=0.0D0
        RETURN
        END
