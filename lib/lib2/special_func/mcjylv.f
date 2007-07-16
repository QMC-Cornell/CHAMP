        PROGRAM MCJYLV

C       ========================================================
C       Purpose: This program computes Bessel functions Jv(z)
C                and Yv(z) and their derivatives with a large
C                order and complex argument using subroutine
C                CJYLV
C       Input:   v --- Order of Jv(z) and Yv(z)
C                z --- Complex argument
C       Output:  CBJV --- Jv(z)
C                CDJV --- Jv'(z)
C                CBYV --- Yv(z)
C                CDYV --- Yv'(z)
C       Examples:
C                v = 100.00,    z = 4.00 + 2.00 i

C                Jv(z) = -.6444792518-123 + .6619157435-123 i
C                Jv'(z)= -.6251103777-122 + .1967638668-121 i
C                Yv(z) =  .2403065353+121 + .2472039414+121 i
C                Yv'(z)= -.7275814786+122 - .2533588851+122 i

C                v =100.5,     z = 4.00 + 2.00 i

C                Jv(z) = -.1161315754-123 + .7390127781-124 i
C                Jv'(z)= -.1588519437-122 + .2652227059-122 i
C                Yv(z) =  .1941381412+122 + .1237578195+122 i
C                Yv'(z)= -.5143285247+123 - .5320026773+122 i
C       ========================================================

        IMPLICIT DOUBLE PRECISION (V,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        WRITE(*,*)'Please enter v,x and y ( z = x+iy )'
        READ(*,*)V,X,Y
        WRITE(*,10)V,X,Y
        Z=CMPLX(X,Y)
        CALL CJYLV(V,Z,CBJV,CDJV,CBYV,CDYV)
        WRITE(*,*)
        WRITE(*,20)CBJV
        WRITE(*,30)CDJV
        WRITE(*,*)
        WRITE(*,40)CBYV
        WRITE(*,50)CDYV
10      FORMAT(8X,'v = ',F6.2,',    ','z =',F7.2,' + i ',F7.2)
20      FORMAT(8X,'Jv(z) =',D17.10,' + i',D17.10)
30      FORMAT(8X,'Jv''(z)=',D17.10,' + i',D17.10)
40      FORMAT(8X,'Yv(z) =',D17.10,' + i',D17.10)
50      FORMAT(8X,'Yv''(z)=',D17.10,' + i',D17.10)
        END


        SUBROUTINE CJYLV(V,Z,CBJV,CDJV,CBYV,CDYV)

C       ===================================================
C       Purpose: Compute Bessel functions Jv(z) and Yv(z)
C                and their derivatives with a complex
C                argument and a large order
C       Input:   v --- Order of Jv(z) and Yv(z)
C                z --- Complex argument
C       Output:  CBJV --- Jv(z)
C                CDJV --- Jv'(z)
C                CBYV --- Yv(z)
C                CDYV --- Yv'(z)
C       Routine called:
C                CJK to compute the expansion coefficients
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CF(12),A(91)
        KM=12
        CALL CJK(KM,A)
        PI=3.141592653589793D0
        DO 30 L=1,0,-1
           V0=V-L
           CWS=CDSQRT(1.0D0-(Z/V0)*(Z/V0))
           CETA=CWS+CDLOG(Z/V0/(1.0D0+CWS))
           CT=1.0D0/CWS
           CT2=CT*CT
           DO 15 K=1,KM
              L0=K*(K+1)/2+1
              LF=L0+K
              CF(K)=A(LF)
              DO 10 I=LF-1,L0,-1
10               CF(K)=CF(K)*CT2+A(I)
15            CF(K)=CF(K)*CT**K
           VR=1.0D0/V0
           CSJ=(1.0D0,0.0D0)
           DO 20 K=1,KM
20            CSJ=CSJ+CF(K)*VR**K
           CBJV=CDSQRT(CT/(2.0D0*PI*V0))*CDEXP(V0*CETA)*CSJ
           IF (L.EQ.1) CFJ=CBJV
           CSY=(1.0D0,0.0D0)
           DO 25 K=1,KM
25            CSY=CSY+(-1)**K*CF(K)*VR**K
           CBYV=-CDSQRT(2.0D0*CT/(PI*V0))*CDEXP(-V0*CETA)*CSY
           IF (L.EQ.1) CFY=CBYV
30      CONTINUE
        CDJV=-V/Z*CBJV+CFJ
        CDYV=-V/Z*CBYV+CFY
        RETURN
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

