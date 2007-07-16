        PROGRAM MCJK

C       ============================================================
C       Purpose: This program computes the expansion coefficients
C                for the asymptotic expansion of Bessel functions
C                with large orders using subroutine CJK
C       Input :  Km   --- Maximum k
C       Output:  A(L) --- Cj(k) where j and k are related to L by
C                         L=j+1+[k*(k+1)]/2; j,k=0,1,2,...,Km
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(231)
        WRITE(*,*)'Please enter Km ( ó 20 )'
        READ(*,*)KM
        LM=KM+1+(KM*(KM+1))/2
        CALL CJK(KM,A)
        DO 10 K=1,LM
10         WRITE(*,15)K,A(K)
15      FORMAT(1X,I3,D25.14)
        END


        SUBROUTINE CJK(KM,A)

C       ========================================================
C       Purpose: Compute the expansion coefficients for the
C                asymptotic expansion of Bessel functions
C                with large orders
C       Input :  Km   --- Maximum k
C       Output:  A(L) --- Cj(k) where j and k are related to L
C                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(*)
        A(1)=1.0D0
        F0=1.0D0
        G0=1.0D0
        DO 10 K=0,KM-1
           L1=(K+1)*(K+2)/2+1
           L2=(K+1)*(K+2)/2+K+2
           F=(0.5D0*K+0.125D0/(K+1))*F0
           G=-(1.5D0*K+0.625D0/(3.0*(K+1.0D0)))*G0
           A(L1)=F
           A(L2)=G
           F0=F
10         G0=G
        DO 15 K=1,KM-1
           DO 15 J=1,K
              L3=K*(K+1)/2+J+1
              L4=(K+1)*(K+2)/2+J+1
              A(L4)=(J+0.5D0*K+0.125D0/(2.0*J+K+1.0))*A(L3)
     &             -(J+0.5D0*K-1.0+0.625D0/(2.0*J+K+1.0))*A(L3-1)
15         CONTINUE
        RETURN
        END
