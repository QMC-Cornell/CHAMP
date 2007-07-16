        PROGRAM MEULERA

C       ==========================================================
C       Purpose: This program computes Euler number En using
C                subroutine EULERA
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
        CALL EULERA(N,E)
        WRITE(*,*)'   n            En'
        WRITE(*,*)' --------------------------'
        DO 10 K=0,N,2
10         WRITE(*,20)K,E(K)
20      FORMAT(2X,I3,D22.12)
        END


        SUBROUTINE EULERA(N,EN)

C       ======================================
C       Purpose: Compute Euler number En
C       Input :  n --- Serial number
C       Output:  EN(n) --- En
C       ======================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:N)
        EN(0)=1.0D0
        DO 30 M=1,N/2
           S=1.0D0
           DO 20 K=1,M-1
              R=1.0D0
              DO 10 J=1,2*K
10               R=R*(2.0D0*M-2.0D0*K+J)/J
20            S=S+R*EN(2*K)
30         EN(2*M)=-S
        RETURN
        END
