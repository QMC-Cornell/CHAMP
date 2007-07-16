        PROGRAM MCJYVB

C       =============================================================
C       Purpose: This program computes Bessel functions Jv(z), Yv(z),
C                and their derivatives for a complex argument using
C                subroutine CJYVB
C       Input :  z --- Complex argument
C                v --- Order of Jv(z) and Yv(z)
C                      ( v = n+v0, 0 ף n ף 250, 0 ף v0 < 1 )
C       Output:  CBJ(n) --- Jn+v0(z)
C                CDJ(n) --- Jn+v0'(z)
C                CBY(n) --- Yn+v0(z)
C                CDY(n) --- Yn+v0'(z)
C       Example:
C                v = n +v0,  v0 = 1/3,   z = 4.0 + i 2.0

C     n     Re[Jv(z)]       Im[Jv(z)]      Re[Jv'(z)]      Im[Jv'(z)]
C    -------------------------------------------------------------------
C     0  -.13829878D+01  -.30855145D+00  -.18503756D+00   .13103689D+01
C     1   .82553327D-01  -.12848394D+01  -.12336901D+01   .45079506D-01
C     2   .10843924D+01  -.39871046D+00  -.33046401D+00  -.84574964D+00
C     3   .74348135D+00   .40665987D+00   .45318486D+00  -.42198992D+00
C     4   .17802266D+00   .44526939D+00   .39624497D+00   .97902890D-01
C     5  -.49008598D-01   .21085409D+00   .11784299D+00   .19422044D+00

C     n     Re[Yv(z)]      Im[Yv(z)]       Re[Yv'(z)]      Im[Yv'(z)]
C    -------------------------------------------------------------------
C     0   .34099851D+00  -.13440666D+01  -.13544477D+01  -.15470699D+00
C     1   .13323787D+01   .53735934D-01  -.21467271D-01  -.11807457D+01
C     2   .38393305D+00   .10174248D+01   .91581083D+00  -.33147794D+00
C     3  -.49924295D+00   .71669181D+00   .47786442D+00   .37321597D+00
C     4  -.57179578D+00   .27099289D+00  -.12111686D+00   .23405313D+00
C     5  -.25700924D+00   .24858555D+00  -.43023156D+00  -.13123662D+00
C       =============================================================

        IMPLICIT DOUBLE PRECISION (V,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:250),CDJ(0:250),CBY(0:250),CDY(0:250)
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
        CALL CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)
        NM=INT(VM)
        WRITE(*,*)
        WRITE(*,*)'  n       Re[Jv(z)]       Im[Jv(z)]',
     &            '       Re[Jv''(z)]      Im[Jv''(z)]'
        WRITE(*,*)' ----------------------------------',
     &            '-----------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20) K,CBJ(K),CDJ(K)
        WRITE(*,*)
        WRITE(*,*)'  n       Re[Yv(z)]       Im[Yv(z)]',
     &            '       Re[Yv''(z)]      Im[Yv''(z)]'
        WRITE(*,*)' ----------------------------------',
     &            '-----------------------------------'
        DO 15 K=0,NM,NS
15         WRITE(*,20) K,CBY(K),CDY(K)
20      FORMAT(1X,I3,2X,4D16.8)
25      FORMAT(8X,'v = n+v0',',  v0 =',F5.2,',  z =',F7.2,' +',F7.2,'i')
        END


        SUBROUTINE CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)

C       ===========================================================
C       Purpose: Compute Bessel functions Jv(z), Yv(z) and their
C                derivatives for a complex argument
C       Input :  z --- Complex argument
C                v --- Order of Jv(z) and Yv(z)
C                      ( v = n+v0, n = 0,1,2,..., 0 ף v0 < 1 )
C       Output:  CBJ(n) --- Jn+v0(z)
C                CDJ(n) --- Jn+v0'(z)
C                CBY(n) --- Yn+v0(z)
C                CDY(n) --- Yn+v0'(z)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA for computing the gamma function
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A,B,G,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:*),CDJ(0:*),CBY(0:*),CDY(0:*)
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z1=Z
        Z2=Z*Z
        N=INT(V)
        V0=V-N
        PV0=PI*V0
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           IF (V0.EQ.0.0) THEN
              CBJ(0)=(1.0D0,0.0D0)
              CDJ(1)=(0.5D0,0.0D0)
           ELSE
              CDJ(0)=(1.0D+300,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        IF (REAL(Z).LT.0.0D0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CJV0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 15 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+V0))
              CJV0=CJV0+CR
              IF (CDABS(CR).LT.CDABS(CJV0)*1.0D-15) GO TO 20
15         CONTINUE
20         VG=1.0D0+V0
           CALL GAMMA(VG,GA)
           CA=(0.5D0*Z1)**V0/GA
           CJV0=CJV0*CA
        ELSE
           K0=11
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           VV=4.0D0*V0*V0
           CPZ=(1.0D0,0.0D0)
           CRP=(1.0D0,0.0D0)
           DO 25 K=1,K0
              CRP=-0.78125D-2*CRP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &            (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*Z2)
25            CPZ=CPZ+CRP
           CQZ=(1.0D0,0.0D0)
           CRQ=(1.0D0,0.0D0)
           DO 30 K=1,K0
              CRQ=-0.78125D-2*CRQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &            (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*Z2)
30            CQZ=CQZ+CRQ
           CQZ=0.125D0*(VV-1.0)*CQZ/Z1
           ZK=Z1-(0.5D0*V0+0.25D0)*PI
           CA0=CDSQRT(RP2/Z1)
           CCK=CDCOS(ZK)
           CSK=CDSIN(ZK)
           CJV0=CA0*(CPZ*CCK-CQZ*CSK)
           CYV0=CA0*(CPZ*CSK+CQZ*CCK)
        ENDIF
        IF (A0.LE.12.0) THEN
           IF (V0.NE.0.0) THEN
              CJVN=(1.0D0,0.0D0)
              CR=(1.0D0,0.0D0)
              DO 35 K=1,40
                 CR=-0.25D0*CR*Z2/(K*(K-V0))
                 CJVN=CJVN+CR
                 IF (CDABS(CR).LT.CDABS(CJVN)*1.0D-15) GO TO 40
35            CONTINUE
40            VG=1.0D0-V0
              CALL GAMMA(VG,GB)
              CB=(2.0D0/Z1)**V0/GB
              CJU0=CJVN*CB
              CYV0=(CJV0*DCOS(PV0)-CJU0)/DSIN(PV0)
           ELSE
              CEC=CDLOG(Z1/2.0D0)+.5772156649015329D0
              CS0=(0.0D0,0.0D0)
              W0=0.0D0
              CR0=(1.0D0,0.0D0)
              DO 45 K=1,30
                 W0=W0+1.0D0/K
                 CR0=-0.25D0*CR0/(K*K)*Z2
45               CS0=CS0+CR0*W0
              CYV0=RP2*(CEC*CJV0-CS0)
           ENDIF
        ENDIF
        IF (N.EQ.0) N=1
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 50 K=M,0,-1
           CF=2.0D0*(V0+K+1.0D0)/Z1*CF1-CF2
           IF (K.LE.N) CBJ(K)=CF
           CF2=CF1
50         CF1=CF
        CS=CJV0/CF
        DO 55 K=0,N
55         CBJ(K)=CS*CBJ(K)
        IF (REAL(Z).LT.0.0D0) THEN
           CFAC0=CDEXP(PV0*CI)
           IF (DIMAG(Z).LT.0.0D0) THEN
              CYV0=CFAC0*CYV0-2.0D0*CI*DCOS(PV0)*CJV0
           ELSE IF (DIMAG(Z).GT.0.0D0) THEN
              CYV0=CYV0/CFAC0+2.0D0*CI*DCOS(PV0)*CJV0
           ENDIF
           DO 60 K=0,N
              IF (DIMAG(Z).LT.0.0D0) THEN
                 CBJ(K)=CDEXP(-PI*(K+V0)*CI)*CBJ(K)
              ELSE IF (DIMAG(Z).GT.0.0D0) THEN
                 CBJ(K)=CDEXP(PI*(K+V0)*CI)*CBJ(K)
              ENDIF
60         CONTINUE
           Z1=Z1
        ENDIF
        CBY(0)=CYV0
        DO 65 K=1,N
           CYY=(CBJ(K)*CBY(K-1)-2.0D0/(PI*Z))/CBJ(K-1)
           CBY(K)=CYY
65      CONTINUE
        CDJ(0)=V0/Z*CBJ(0)-CBJ(1)
        DO 70 K=1,N
70         CDJ(K)=-(K+V0)/Z*CBJ(K)+CBJ(K-1)
        CDY(0)=V0/Z*CBY(0)-CBY(1)
        DO 75 K=1,N
75         CDY(K)=CBY(K-1)-(K+V0)/Z*CBY(K)
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
