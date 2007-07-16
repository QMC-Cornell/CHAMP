        PROGRAM MCJYVA

C       ===============================================================
C       Purpose: This program computes Bessel functions Jv(z), Yv(z),
C                and their derivatives for a complex argument using
C                subroutine CJYVA
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
C    ------------------------------------------------------------------
C     0  -.13829878D+01  -.30855145D+00  -.18503756D+00   .13103689D+01
C     1   .82553327D-01  -.12848394D+01  -.12336901D+01   .45079506D-01
C     2   .10843924D+01  -.39871046D+00  -.33046401D+00  -.84574964D+00
C     3   .74348135D+00   .40665987D+00   .45318486D+00  -.42198992D+00
C     4   .17802266D+00   .44526939D+00   .39624497D+00   .97902890D-01
C     5  -.49008598D-01   .21085409D+00   .11784299D+00   .19422044D+00

C     n     Re[Yv(z)]      Im[Yv(z)]       Re[Yv'(z)]      Im[Yv'(z)]
C    ------------------------------------------------------------------
C     0   .34099851D+00  -.13440666D+01  -.13544477D+01  -.15470699D+00
C     1   .13323787D+01   .53735934D-01  -.21467271D-01  -.11807457D+01
C     2   .38393305D+00   .10174248D+01   .91581083D+00  -.33147794D+00
C     3  -.49924295D+00   .71669181D+00   .47786442D+00   .37321597D+00
C     4  -.57179578D+00   .27099289D+00  -.12111686D+00   .23405313D+00
C     5  -.25700924D+00   .24858555D+00  -.43023156D+00  -.13123662D+00
C       ===============================================================

        IMPLICIT DOUBLE PRECISION (V,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        COMMON CBJ(0:251),CDJ(0:251),CBY(0:251),CDY(0:251)
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
        CALL CJYVA(V,Z,VM,CBJ,CDJ,CBY,CDY)
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


        SUBROUTINE CJYVA(V,Z,VM,CBJ,CDJ,CBY,CDY)

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
        PV1=PI*(1.0D0+V0)
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
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           DO 25 L=0,1
              VL=V0+L
              CJVL=(1.0D0,0.0D0)
              CR=(1.0D0,0.0D0)
              DO 15 K=1,40
                 CR=-0.25D0*CR*Z2/(K*(K+VL))
                 CJVL=CJVL+CR
                 IF (CDABS(CR).LT.CDABS(CJVL)*1.0D-15) GO TO 20
15            CONTINUE
20            VG=1.0D0+VL
              CALL GAMMA(VG,GA)
              CA=(0.5D0*Z1)**VL/GA
              IF (L.EQ.0) CJV0=CJVL*CA
              IF (L.EQ.1) CJV1=CJVL*CA
25         CONTINUE
        ELSE
           K0=11
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           DO 40 J=0,1
              VV=4.0D0*(J+V0)*(J+V0)
              CPZ=(1.0D0,0.0D0)
              CRP=(1.0D0,0.0D0)
              DO 30 K=1,K0
                 CRP=-0.78125D-2*CRP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &               (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*Z2)
30               CPZ=CPZ+CRP
              CQZ=(1.0D0,0.0D0)
              CRQ=(1.0D0,0.0D0)
              DO 35 K=1,K0
                 CRQ=-0.78125D-2*CRQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &               (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*Z2)
35               CQZ=CQZ+CRQ
              CQZ=0.125D0*(VV-1.0)*CQZ/Z1
              ZK=Z1-(0.5D0*(J+V0)+0.25D0)*PI
              CA0=CDSQRT(RP2/Z1)
              CCK=CDCOS(ZK)
              CSK=CDSIN(ZK)
              IF (J.EQ.0) THEN
                 CJV0=CA0*(CPZ*CCK-CQZ*CSK)
                 CYV0=CA0*(CPZ*CSK+CQZ*CCK)
              ELSE IF (J.EQ.1) THEN
                 CJV1=CA0*(CPZ*CCK-CQZ*CSK)
                 CYV1=CA0*(CPZ*CSK+CQZ*CCK)
              ENDIF
40         CONTINUE
        ENDIF
        IF (A0.LE.12.0) THEN
           IF (V0.NE.0.0) THEN
              DO 55 L=0,1
                 VL=V0+L
                 CJVL=(1.0D0,0.0D0)
                 CR=(1.0D0,0.0D0)
                 DO 45 K=1,40
                    CR=-0.25D0*CR*Z2/(K*(K-VL))
                    CJVL=CJVL+CR
                    IF (CDABS(CR).LT.CDABS(CJVL)*1.0D-15) GO TO 50
45               CONTINUE
50               VG=1.0D0-VL
                 CALL GAMMA(VG,GB)
                 CB=(2.0D0/Z1)**VL/GB
                 IF (L.EQ.0) CJU0=CJVL*CB
                 IF (L.EQ.1) CJU1=CJVL*CB
55            CONTINUE
              CYV0=(CJV0*DCOS(PV0)-CJU0)/DSIN(PV0)
              CYV1=(CJV1*DCOS(PV1)-CJU1)/DSIN(PV1)
           ELSE
              CEC=CDLOG(Z1/2.0D0)+.5772156649015329D0
              CS0=(0.0D0,0.0D0)
              W0=0.0D0
              CR0=(1.0D0,0.0D0)
              DO 60 K=1,30
                 W0=W0+1.0D0/K
                 CR0=-0.25D0*CR0/(K*K)*Z2
60               CS0=CS0+CR0*W0
              CYV0=RP2*(CEC*CJV0-CS0)
              CS1=(1.0D0,0.0D0)
              W1=0.0D0
              CR1=(1.0D0,0.0D0)
              DO 65 K=1,30
                 W1=W1+1.0D0/K
                 CR1=-0.25D0*CR1/(K*(K+1))*Z2
65               CS1=CS1+CR1*(2.0D0*W1+1.0D0/(K+1.0D0))
              CYV1=RP2*(CEC*CJV1-1.0D0/Z1-0.25D0*Z1*CS1)
           ENDIF
        ENDIF
        IF (REAL(Z).LT.0.0D0) THEN
           CFAC0=CDEXP(PV0*CI)
           CFAC1=CDEXP(PV1*CI)
           IF (DIMAG(Z).LT.0.0D0) THEN
              CYV0=CFAC0*CYV0-2.0D0*CI*DCOS(PV0)*CJV0
              CYV1=CFAC1*CYV1-2.0D0*CI*DCOS(PV1)*CJV1
              CJV0=CJV0/CFAC0
              CJV1=CJV1/CFAC1
           ELSE IF (DIMAG(Z).GT.0.0D0) THEN
              CYV0=CYV0/CFAC0+2.0D0*CI*DCOS(PV0)*CJV0
              CYV1=CYV1/CFAC1+2.0D0*CI*DCOS(PV1)*CJV1
              CJV0=CFAC0*CJV0
              CJV1=CFAC1*CJV1
           ENDIF
        ENDIF
        CBJ(0)=CJV0
        CBJ(1)=CJV1
        IF (N.GE.2.AND.N.LE.INT(0.25*A0)) THEN
           CF0=CJV0
           CF1=CJV1
           DO 70 K=2,N
              CF=2.0D0*(K+V0-1.0D0)/Z*CF1-CF0
              CBJ(K)=CF
              CF0=CF1
70            CF1=CF
        ELSE IF (N.GE.2) THEN
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              N=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 75 K=M,0,-1
              CF=2.0D0*(V0+K+1.0D0)/Z*CF1-CF2
              IF (K.LE.N) CBJ(K)=CF
              CF2=CF1
75            CF1=CF
           IF (CDABS(CJV0).GT.CDABS(CJV1)) CS=CJV0/CF
           IF (CDABS(CJV0).LE.CDABS(CJV1)) CS=CJV1/CF2
           DO 80 K=0,N
80            CBJ(K)=CS*CBJ(K)
        ENDIF
        CDJ(0)=V0/Z*CBJ(0)-CBJ(1)
        DO 85 K=1,N
85         CDJ(K)=-(K+V0)/Z*CBJ(K)+CBJ(K-1)
        CBY(0)=CYV0
        CBY(1)=CYV1
        YA0=CDABS(CYV0)
        LB=0
        CG0=CYV0
        CG1=CYV1
        DO 90 K=2,N
           CYK=2.0D0*(V0+K-1.0D0)/Z*CG1-CG0
           IF (CDABS(CYK).GT.1.0D+290) GO TO 90
           YAK=CDABS(CYK)
           YA1=CDABS(CG0)
           IF (YAK.LT.YA0.AND.YAK.LT.YA1) LB=K
           CBY(K)=CYK
           CG0=CG1
           CG1=CYK
90      CONTINUE
        IF (LB.LE.4.OR.DIMAG(Z).EQ.0.0D0) GO TO 125
95      IF (LB.EQ.LB0) GO TO 125
        CH2=(1.0D0,0.0D0)
        CH1=(0.0D0,0.0D0)
        LB0=LB
        DO 100 K=LB,1,-1
           CH0=2.0D0*(K+V0)/Z*CH1-CH2
           CH2=CH1
100        CH1=CH0
        CP12=CH0
        CP22=CH2
        CH2=(0.0D0,0.0D0)
        CH1=(1.0D0,0.0D0)
        DO 105 K=LB,1,-1
           CH0=2.0D0*(K+V0)/Z*CH1-CH2
           CH2=CH1
105        CH1=CH0
        CP11=CH0
        CP21=CH2
        IF (LB.EQ.N) CBJ(LB+1)=2.0D0*(LB+V0)/Z*CBJ(LB)-CBJ(LB-1)
        IF (CDABS(CBJ(0)).GT.CDABS(CBJ(1))) THEN
           CBY(LB+1)=(CBJ(LB+1)*CYV0-2.0D0*CP11/(PI*Z))/CBJ(0)
           CBY(LB)=(CBJ(LB)*CYV0+2.0D0*CP12/(PI*Z))/CBJ(0)
        ELSE
           CBY(LB+1)=(CBJ(LB+1)*CYV1-2.0D0*CP21/(PI*Z))/CBJ(1)
           CBY(LB)=(CBJ(LB)*CYV1+2.0D0*CP22/(PI*Z))/CBJ(1)
        ENDIF
        CYL2=CBY(LB+1)
        CYL1=CBY(LB)
        DO 110 K=LB-1,0,-1
           CYLK=2.0D0*(K+V0+1.0D0)/Z*CYL1-CYL2
           CBY(K)=CYLK
           CYL2=CYL1
110        CYL1=CYLK
        CYL1=CBY(LB)
        CYL2=CBY(LB+1)
        DO 115 K=LB+1,N-1
           CYLK=2.0D0*(K+V0)/Z*CYL2-CYL1
           CBY(K+1)=CYLK
           CYL1=CYL2
115        CYL2=CYLK
        DO 120 K=2,N
           WA=CDABS(CBY(K))
           IF (WA.LT.CDABS(CBY(K-1))) LB=K
120     CONTINUE
        GO TO 95
125     CDY(0)=V0/Z*CBY(0)-CBY(1)
        DO 130 K=1,N
130        CDY(K)=CBY(K-1)-(K+V0)/Z*CBY(K)
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
