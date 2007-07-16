        PROGRAM MCIKVB

C       =============================================================
C       Purpose: This program computes the modified Bessel functions
C                Iv(z), Kv(z) and their derivatives for an arbitrary
C                order and complex argument using subroutine CIKVB
C       Input :  z --- Complex argument z
C                v --- Real order of Iv(z) and Kv(z)
C                      ( v = n+v0,  0 ף n ף 250, 0 ף v0 < 1 )
C       Output:  CBI(n) --- In+v0(z)
C                CDI(n) --- In+v0'(z)
C                CBK(n) --- Kn+v0(z)
C                CDK(n) --- Kn+v0'(z)
C       Example: v= n+v0,   v0 = .25,   z = 4.0+ i 2.0

C    n     Re[Iv(z)]      Im[Iv(z)]     Re[Iv'(z)]     Im[Iv'(z)]
C  -----------------------------------------------------------------
C    0  -.19336550D+01  .10328998D+02 -.23119621D+01  .91612230D+01
C    1  -.24735044D+01  .85964317D+01 -.23898329D+01  .78707023D+01
C    2  -.28460107D+01  .54124063D+01 -.24105909D+01  .55204965D+01
C    3  -.23476775D+01  .24445612D+01 -.21145027D+01  .30604463D+01
C    4  -.13829947D+01  .70848630D+00 -.14732387D+01  .12545751D+01
C    5  -.59879982D+00  .64588999D-01 -.78816416D+00  .32629794D+00

C    n     Re[Kv(z)]      Im[Kv(z)]     Re[Kv'(z)]     Im[Kv'(z)]
C   ----------------------------------------------------------------
C    0  -.64820386D-02 -.84715754D-02  .75118612D-02  .89920077D-02
C    1  -.80477525D-02 -.92535355D-02  .96506687D-02  .97789903D-02
C    2  -.12819299D-01 -.11086405D-01  .16310878D-01  .11358076D-01
C    3  -.24574004D-01 -.13462616D-01  .33167751D-01  .11850554D-01
C    4  -.53516204D-01 -.12614703D-01  .75424026D-01  .14407268D-02
C    5  -.12627405D+00  .10581162D-01  .18054884D+00 -.64789392D-01
C       =============================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)
        WRITE(*,*)'  Please enter v, x and y ( z=x+iy )'
        READ(*,*)V,X,Y
        Z=CMPLX(X,Y)
        N=INT(V)
        V0=V-N
        WRITE(*,25)V0,X,Y
        IF (N.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL CIKVB(V,Z,VM,CBI,CDI,CBK,CDK)
        NM=INT(VM)
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Iv(z)]       Im[Iv(z)] ',
     &            '     Re[Iv''(z)]      Im[Iv''(z)] '
        WRITE(*,*)' ---------------------------------------',
     &            '------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20) K,CBI(K),CDI(K)
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Kv(z)]       Im[Kv(z)] ',
     &            '     Re[Kv''(z)]      Im[Kv''(z)] '
        WRITE(*,*)' ---------------------------------------',
     &            '------------------------------'
        DO 15 K=0,NM,NS
15         WRITE(*,20)K,CBK(K),CDK(K)
20      FORMAT(1X,I4,1X,4D16.8)
25      FORMAT(8X,'v= n+v0',',   ','v0 =',F6.2,',   ','z =',F5.1,
     &        '+ i',F5.1)
        END


        SUBROUTINE CIKVB(V,Z,VM,CBI,CDI,CBK,CDK)

C       ===========================================================
C       Purpose: Compute the modified Bessel functions Iv(z), Kv(z)
C                and their derivatives for an arbitrary order and
C                complex argument
C       Input :  z --- Complex argument z
C                v --- Real order of Iv(z) and Kv(z)
C                      ( v =n+v0, n = 0,1,2,..., 0 ף v0 < 1 )
C       Output:  CBI(n) --- In+v0(z)
C                CDI(n) --- In+v0'(z)
C                CBK(n) --- Kn+v0(z)
C                CDK(n) --- Kn+v0'(z)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA for computing the gamma function
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:*),CDI(0:*),CBK(0:*),CDK(0:*)
        Z1=Z
        Z2=Z*Z
        A0=CDABS(Z)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        N=INT(V)
        V0=V-N
        PIV=PI*V0
        VT=4.0D0*V0*V0
        IF (N.EQ.0) N=1
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBI(K)=0.0D0
              CDI(K)=0.0D0
              CBK(K)=-1.0D+300
10            CDK(K)=1.0D+300
           IF (V0.EQ.0.0) THEN
              CBI(0)=(1.0D0,0.0D0)
              CDI(1)=(0.5D0,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        K0=14
        IF (A0.GE.35.0) K0=10
        IF (A0.GE.50.0) K0=8
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LT.18.0) THEN
           IF (V0.EQ.0.0) THEN
              CA1=(1.0D0,0.0D0)
           ELSE
              V0P=1.0D0+V0
              CALL GAMMA(V0P,GAP)
              CA1=(0.5D0*Z1)**V0/GAP
           ENDIF
           CI0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 15 K=1,50
              CR=0.25D0*CR*Z2/(K*(K+V0))
              CI0=CI0+CR
              IF (CDABS(CR/CI0).LT.1.0D-15) GO TO 20
15         CONTINUE
20         CBI0=CI0*CA1
        ELSE
           CA=CDEXP(Z1)/CDSQRT(2.0D0*PI*Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 25 K=1,K0
              CR=-0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
25            CS=CS+CR
           CBI0=CA*CS
        ENDIF
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 30 K=M,0,-1
           CF=2.0D0*(V0+K+1.0D0)/Z1*CF1+CF2
           IF (K.LE.N) CBI(K)=CF
           CF2=CF1
30         CF1=CF
        CS=CBI0/CF
        DO 35 K=0,N
35         CBI(K)=CS*CBI(K)
        IF (A0.LE.9.0) THEN
           IF (V0.EQ.0.0) THEN
              CT=-CDLOG(0.5D0*Z1)-0.5772156649015329D0
              CS=(0.0D0,0.0D0)
              W0=0.0D0
              CR=(1.0D0,0.0D0)
              DO 40 K=1,50
                 W0=W0+1.0D0/K
                 CR=0.25D0*CR/(K*K)*Z2
                 CP=CR*(W0+CT)
                 CS=CS+CP
                 IF (K.GE.10.AND.CDABS(CP/CS).LT.1.0D-15) GO TO 45
40            CONTINUE
45            CBK0=CT+CS
           ELSE
              V0N=1.0D0-V0
              CALL GAMMA(V0N,GAN)
              CA2=1.0D0/(GAN*(0.5D0*Z1)**V0)
              CA1=(0.5D0*Z1)**V0/GAP
              CSU=CA2-CA1
              CR1=(1.0D0,0.0D0)
              CR2=(1.0D0,0.0D0)
              DO 50 K=1,50
                 CR1=0.25D0*CR1*Z2/(K*(K-V0))
                 CR2=0.25D0*CR2*Z2/(K*(K+V0))
                 CP=CA2*CR1-CA1*CR2
                 CSU=CSU+CP
                 IF (K.GE.10.AND.CDABS(CP/CSU).LT.1.0D-15) GO TO 55
50            CONTINUE
55            CBK0=0.5D0*PI*CSU/DSIN(PIV)
           ENDIF
        ELSE
           CB=CDEXP(-Z1)*CDSQRT(0.5D0*PI/Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 60 K=1,K0
              CR=0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
60            CS=CS+CR
           CBK0=CB*CS
        ENDIF
        CBK(0)=CBK0
        IF (REAL(Z).LT.0.0) THEN
           DO 65 K=0,N
              CVK=CDEXP((K+V0)*PI*CI)
              IF (DIMAG(Z).LT.0.0D0) THEN
                 CBK(K)=CVK*CBK(K)+PI*CI*CBI(K)
                 CBI(K)=CBI(K)/CVK
              ELSE IF (DIMAG(Z).GT.0.0) THEN
                 CBK(K)=CBK(K)/CVK-PI*CI*CBI(K)
                 CBI(K)=CVK*CBI(K)
              ENDIF
65         CONTINUE
        ENDIF
        DO 70 K=1,N
           CKK=(1.0D0/Z-CBI(K)*CBK(K-1))/CBI(K-1)
           CBK(K)=CKK
70      CONTINUE
        CDI(0)=V0/Z*CBI(0)+CBI(1)
        CDK(0)=V0/Z*CBK(0)-CBK(1)
        DO 80 K=1,N
           CDI(K)=-(K+V0)/Z*CBI(K)+CBI(K-1)
80         CDK(K)=-(K+V0)/Z*CBK(K)-CBK(K-1)
        VM=N+V0
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)

C       ==================================================
C       Purpose: Compute gamma function ג(x)
C       Input :  x  --- Argument of ג(x)
C                       ( x is not equal to 0,-1,-2,תתת)
C       Output:  GA --- ג(x)
C       ==================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)

C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that the magnitude of
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)

C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
