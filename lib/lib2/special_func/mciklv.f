        PROGRAM MCIKLV

C       =========================================================
C       Purpose: This program computes modified Bessel functions
C                Iv(z) and Kv(z) and their derivatives for a
C                large order and a complex argument using
C                subroutine CIKLV
C       Input:   v --- Order of Iv(z) and Kv(z)
C                z --- Complex argument
C       Output:  CBIV --- Iv(z)
C                CDIV --- Iv'(z)
C                CBKV --- Kv(z)
C                CDKV --- Kv'(z)
C       Examples:
C                v =100.00,    z =   4.00 + i   2.00

C       Iv(z) = -.7373606617-123 + .6461109082-123 i
C       Iv'(z)= -.8307094243-122 + .2030132500-121 i
C       Kv(z) = -.3836166007+121 - .3356017795+121 i
C       Kv'(z)=  .1103271276+123 + .2886519240+122 i

C                v =100.50,    z =   4.00 + i   2.00
C       Iv(z) = -.1289940051-123 + .6845756182-124 i
C       Iv'(z)= -.1907996261-122 + .2672465997-122 i
C       Kv(z) = -.3008779281+122 - .1593719779+122 i
C       Kv'(z)=  .7653781978+123 + .1857772148+122 i
C       =========================================================

        IMPLICIT DOUBLE PRECISION (V,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        WRITE(*,*)'Please enter v,x,y ( z = x+iy )'
        READ(*,*)V,X,Y
        WRITE(*,10)V,X,Y
        Z=CMPLX(X,Y)
        CALL CIKLV(V,Z,CBIV,CDIV,CBKV,CDKV)
        WRITE(*,*)
        WRITE(*,20)CBIV
        WRITE(*,30)CDIV
        WRITE(*,*)
        WRITE(*,40)CBKV
        WRITE(*,50)CDKV
10      FORMAT(8X,'v =',F6.2,',    ','z =',F7.2,' + i',F7.2)
20      FORMAT(8X,'Iv(z) =',D17.10,' + i ',D17.10)
30      FORMAT(8X,'Iv''(z)=',D17.10,' + i ',D17.10)
40      FORMAT(8X,'Kv(z) =',D17.10,' + i ',D17.10)
50      FORMAT(8X,'Kv''(z)=',D17.10,' + i ',D17.10)
        END



        SUBROUTINE CIKLV(V,Z,CBIV,CDIV,CBKV,CDKV)

C       =====================================================
C       Purpose: Compute modified Bessel functions Iv(z) and
C                Kv(z) and their derivatives with a complex
C                argument and a large order
C       Input:   v --- Order of Iv(z) and Kv(z)
C                z --- Complex argument
C       Output:  CBIV --- Iv(z)
C                CDIV --- Iv'(z)
C                CBKV --- Kv(z)
C                CDKV --- Kv'(z)
C       Routine called:
C                CJK to compute the expansion coefficients
C       ====================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CF(12),A(91)
        PI=3.141592653589793D0
        KM=12
        CALL CJK(KM,A)
        DO 30 L=1,0,-1
           V0=V-L
           CWS=CDSQRT(1.0D0+(Z/V0)*(Z/V0))
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
           CSI=(1.0D0,0.0D0)
           DO 20 K=1,KM
20            CSI=CSI+CF(K)*VR**K
           CBIV=CDSQRT(CT/(2.0D0*PI*V0))*CDEXP(V0*CETA)*CSI
           IF (L.EQ.1) CFI=CBIV
           CSK=(1.0D0,0.0D0)
           DO 25 K=1,KM
25            CSK=CSK+(-1)**K*CF(K)*VR**K
           CBKV=CDSQRT(PI*CT/(2.0D0*V0))*CDEXP(-V0*CETA)*CSK
           IF (L.EQ.1) CFK=CBKV
30      CONTINUE
        CDIV=CFI-V/Z*CBIV
        CDKV=-CFK-V/Z*CBKV
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

