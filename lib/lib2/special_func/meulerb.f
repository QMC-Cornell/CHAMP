        PROGRAM MEULERB

C       ==========================================================
C       Purpose: This program computes Euler number En using
C                subroutine EULERB
C       Example: Compute Euler number En for n = 0,2,...,10
C                Computed results:

C                   n            En
C                 --------------------------
C                   0     .100000000000D+01
C                   2    -.100000000000D+01
C                   4     .500000000000D+01
C                   6    -.610000000000D+02
C                   8     .138500000000D+04
C                  10    -.505210000000D+05
C       ==========================================================

        DOUBLE PRECISION E
        DIMENSION E(0:200)
        WRITE(*,*)'  Please enter Nmax '
        READ(*,*)N
        CALL EULERB(N,E)
        WRITE(*,*)'   n            En'
        WRITE(*,*)' --------------------------'
        DO 10 K=0,N,2
10         WRITE(*,20)K,E(K)
20      FORMAT(2X,I3,D22.12)
        END


        SUBROUTINE EULERB(N,EN)

C       ======================================
C       Purpose: Compute Euler number En
C       Input :  n --- Serial number
C       Output:  EN(n) --- En
C       ======================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:N)
        HPI=2.0D0/3.141592653589793D0
        EN(0)=1.0D0
        EN(2)=-1.0D0
        R1=-4.0D0*HPI**3
        DO 20 M=4,N,2
           R1=-R1*(M-1)*M*HPI*HPI
           R2=1.0D0
           ISGN=1.0D0
           DO 10 K=3,1000,2
              ISGN=-ISGN
              S=(1.0D0/K)**(M+1)
              R2=R2+ISGN*S
              IF (S.LT.1.0D-15) GOTO 20
10         CONTINUE
20         EN(M)=R1*R2
        RETURN
        END
