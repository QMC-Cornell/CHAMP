        PROGRAM MCJYNA

C       ================================================================
C       Purpose: This program computes Bessel functions Jn(z), Yn(z)
C                and their derivatives for a complex argument using
C                subroutine CJYNA
C       Input :  z --- Complex argument of Jn(z) and Yn(z)
C                n --- Order of Jn(z) and Yn(z)
C                      ( n = 0,1,תתת, n ף 250 )
C       Output:  CBJ(n) --- Jn(z)
C                CDJ(n) --- Jn'(z)
C                CBY(n) --- Yn(z)
C                CDY(n) --- Yn'(z)
C       Eaxmple: z = 4.0 + i 2.0

C     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
C    -------------------------------------------------------------------
C     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
C     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
C     2   .93050039D+00  -.77959350D+00  -.72363400D+00  -.72836666D+00
C     3   .93991546D+00   .23042918D+00   .29742236D+00  -.63587637D+00
C     4   .33565567D+00   .49215925D+00   .47452722D+00  -.29035945D-01
C     5  -.91389835D-02   .28850107D+00   .20054412D+00   .19908868D+00

C     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
C   --------------------------------------------------------------------
C     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
C     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
C     2   .79074211D+00   .86842120D+00   .78932897D+00  -.70142425D+00
C     3  -.29934789D+00   .89064431D+00   .70315755D+00   .24423024D+00
C     4  -.61557299D+00   .37996071D+00   .41126221D-01   .34044655D+00
C     5  -.38160033D+00   .20975121D+00  -.33884827D+00  -.20590670D-01

C                z = 20.0 + i 10.0 ,      Nmax = 5

C     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
C   --------------------------------------------------------------------
C     0   .15460268D+04  -.10391216D+04  -.10601232D+04  -.15098284D+04
C     1   .10601232D+04   .15098284D+04   .14734253D+04  -.10783122D+04
C     2  -.14008238D+04   .11175029D+04   .11274890D+04   .13643952D+04
C     3  -.11948548D+04  -.12189620D+04  -.11843035D+04   .11920871D+04
C     4   .96778325D+03  -.12666712D+04  -.12483664D+04  -.93887194D+03
C     5   .13018781D+04   .65878188D+03   .64152944D+03  -.12682398D+04

C     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
C   --------------------------------------------------------------------
C     0   .10391216D+04   .15460268D+04   .15098284D+04  -.10601232D+04
C     1  -.15098284D+04   .10601232D+04   .10783122D+04   .14734253D+04
C     2  -.11175029D+04  -.14008238D+04  -.13643952D+04   .11274890D+04
C     3   .12189620D+04  -.11948548D+04  -.11920871D+04  -.11843035D+04
C     4   .12666712D+04   .96778324D+03   .93887194D+03  -.12483664D+04
C     5  -.65878189D+03   .13018781D+04   .12682398D+04   .64152944D+03
C       ================================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        COMMON CBJ(0:251),CDJ(0:251),CBY(0:251),CDY(0:251)
        WRITE(*,*)'  Please enter n, x,y (z=x+iy) '
        READ(*,*)N,X,Y
        Z=CMPLX(X,Y)
        WRITE(*,40)X,Y,N
        IF (N.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL CJYNA(N,Z,NM,CBJ,CDJ,CBY,CDY)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Jn(z)]       Im[Jn(z)]',
     &            '       Re[Jn''(z)]      Im[Jn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,30)K,CBJ(K),CDJ(K)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Yn(z)]       Im[Yn(z)]',
     &            '       Re[Yn''(z)]      Im[Yn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        DO 20 K=0,NM,NS
20         WRITE(*,30)K,CBY(K),CDY(K)
30      FORMAT(1X,I4,4D16.8)
40      FORMAT(3X,3Hz =,F5.1,' + i ',F5.1,' ,',6X,6HNmax =,I3)
        END


        SUBROUTINE CJYNA(N,Z,NM,CBJ,CDJ,CBY,CDY)

C       =======================================================
C       Purpose: Compute Bessel functions Jn(z), Yn(z) and
C                their derivatives for a complex argument
C       Input :  z --- Complex argument of Jn(z) and Yn(z)
C                n --- Order of Jn(z) and Yn(z)
C       Output:  CBJ(n) --- Jn(z)
C                CDJ(n) --- Jn'(z)
C                CBY(n) --- Yn(z)
C                CDY(n) --- Yn'(z)
C                NM --- Highest order computed
C       Rouitines called:
C            (1) CJY01 to calculate J0(z), J1(z), Y0(z), Y1(z)
C            (2) MSTA1 and MSTA2 to calculate the starting
C                point for backward recurrence
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:N),CDJ(0:N),CBY(0:N),CDY(0:N)
        PI=3.141592653589793D0
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 5 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
5             CDY(K)=(1.0D+300,0.0D0)
           CBJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        CALL CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
        CBJ(0)=CBJ0
        CBJ(1)=CBJ1
        CBY(0)=CBY0
        CBY(1)=CBY1
        CDJ(0)=CDJ0
        CDJ(1)=CDJ1
        CDY(0)=CDY0
        CDY(1)=CDY1
        IF (N.LE.1) RETURN
        IF (N.LT.INT(0.25*A0)) THEN
           CJ0=CBJ0
           CJ1=CBJ1
           DO 70 K=2,N
              CJK=2.0D0*(K-1.0D0)/Z*CJ1-CJ0
              CBJ(K)=CJK
              CJ0=CJ1
70            CJ1=CJK
        ELSE
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 75 K=M,0,-1
              CF=2.0D0*(K+1.0D0)/Z*CF1-CF2
              IF (K.LE.NM) CBJ(K)=CF
              CF2=CF1
75            CF1=CF
           IF (CDABS(CBJ0).GT.CDABS(CBJ1)) THEN
              CS=CBJ0/CF
           ELSE
              CS=CBJ1/CF2
           ENDIF
           DO 80 K=0,NM
80            CBJ(K)=CS*CBJ(K)
        ENDIF
        DO 85 K=2,NM
85         CDJ(K)=CBJ(K-1)-K/Z*CBJ(K)
        YA0=CDABS(CBY0)
        LB=0
        CG0=CBY0
        CG1=CBY1
        DO 90 K=2,NM
           CYK=2.0D0*(K-1.0D0)/Z*CG1-CG0
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
           CH0=2.0D0*K/Z*CH1-CH2
           CH2=CH1
100        CH1=CH0
        CP12=CH0
        CP22=CH2
        CH2=(0.0D0,0.0D0)
        CH1=(1.0D0,0.0D0)
        DO 105 K=LB,1,-1
           CH0=2.0D0*K/Z*CH1-CH2
           CH2=CH1
105        CH1=CH0
        CP11=CH0
        CP21=CH2
        IF (LB.EQ.NM) CBJ(LB+1)=2.0D0*LB/Z*CBJ(LB)-CBJ(LB-1)
        IF (CDABS(CBJ(0)).GT.CDABS(CBJ(1))) THEN
           CBY(LB+1)=(CBJ(LB+1)*CBY0-2.0D0*CP11/(PI*Z))/CBJ(0)
           CBY(LB)=(CBJ(LB)*CBY0+2.0D0*CP12/(PI*Z))/CBJ(0)
        ELSE
           CBY(LB+1)=(CBJ(LB+1)*CBY1-2.0D0*CP21/(PI*Z))/CBJ(1)
           CBY(LB)=(CBJ(LB)*CBY1+2.0D0*CP22/(PI*Z))/CBJ(1)
        ENDIF
        CYL2=CBY(LB+1)
        CYL1=CBY(LB)
        DO 110 K=LB-1,0,-1
           CYLK=2.0D0*(K+1.0D0)/Z*CYL1-CYL2
           CBY(K)=CYLK
           CYL2=CYL1
110        CYL1=CYLK
        CYL1=CBY(LB)
        CYL2=CBY(LB+1)
        DO 115 K=LB+1,NM-1
           CYLK=2.0D0*K/Z*CYL2-CYL1
           CBY(K+1)=CYLK
           CYL1=CYL2
115        CYL2=CYLK
        DO 120 K=2,NM
           WA=CDABS(CBY(K))
           IF (WA.LT.CDABS(CBY(K-1))) LB=K
120     CONTINUE
        GO TO 95
125     CONTINUE
        DO 130 K=2,NM
130        CDY(K)=CBY(K-1)-K/Z*CBY(K)
        RETURN
        END


        SUBROUTINE CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)

C       ===========================================================
C       Purpose: Compute complex Bessel functions J0(z), J1(z)
C                Y0(z), Y1(z), and their derivatives
C       Input :  z --- Complex argument
C       Output:  CBJ0 --- J0(z)
C                CDJ0 --- J0'(z)
C                CBJ1 --- J1(z)
C                CDJ1 --- J1'(z)
C                CBY0 --- Y0(z)
C                CDY0 --- Y0'(z)
C                CBY1 --- Y1(z)
C                CDY1 --- Y1'(z)
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION A(12),B(12),A1(12),B1(12)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        RP2=2.0D0/PI
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBJ0=(1.0D0,0.0D0)
           CBJ1=(0.0D0,0.0D0)
           CDJ0=(0.0D0,0.0D0)
           CDJ1=(0.5D0,0.0D0)
           CBY0=-(1.0D300,0.0D0)
           CBY1=-(1.0D300,0.0D0)
           CDY0=(1.0D300,0.0D0)
           CDY1=(1.0D300,0.0D0)
           RETURN
        ENDIF
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CBJ0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,40
              CR=-0.25D0*CR*Z2/(K*K)
              CBJ0=CBJ0+CR
              IF (CDABS(CR/CBJ0).LT.1.0D-15) GO TO 15
10         CONTINUE
15         CBJ1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+1.0D0))
              CBJ1=CBJ1+CR
              IF (CDABS(CR/CBJ1).LT.1.0D-15) GO TO 25
20         CONTINUE
25         CBJ1=0.5D0*Z1*CBJ1
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(0.0D0,0.0D0)
           DO 30 K=1,40
              W0=W0+1.0D0/K
              CR=-0.25D0*CR/(K*K)*Z2
              CP=CR*W0
              CS=CS+CP
              IF (CDABS(CP/CS).LT.1.0D-15) GO TO 35
30         CONTINUE
35         CBY0=RP2*(CDLOG(Z1/2.0D0)+EL)*CBJ0-RP2*CS
           W1=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(1.0D0,0.0D0)
           DO 40 K=1,40
              W1=W1+1.0D0/K
              CR=-0.25D0*CR/(K*(K+1))*Z2
              CP=CR*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS=CS+CP
              IF (CDABS(CP/CS).LT.1.0D-15) GO TO 45
40         CONTINUE
45         CBY1=RP2*((CDLOG(Z1/2.0D0)+EL)*CBJ1-1.0D0/Z1-.25D0*Z1*CS)
        ELSE
           DATA A/-.703125D-01,.112152099609375D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01,
     &            -.1100171402692467D+03,.3038090510922384D+04,
     &            -.1188384262567832D+06,.6252951493434797D+07,
     &            -.4259392165047669D+09,.3646840080706556D+11,
     &            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .732421875D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02,
     &             .5513358961220206D+03,-.1825775547429318D+05,
     &             .8328593040162893D+06,-.5006958953198893D+08,
     &             .3836255180230433D+10,-.3649010818849833D+12,
     &             .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875D+00,-.144195556640625D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01,
     &             .1215978918765359D+03,-.3302272294480852D+04,
     &             .1276412726461746D+06,-.6656367718817688D+07,
     &             .4502786003050393D+09,-.3833857520742790D+11,
     &             .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02,
     &             -.6038440767050702D+03,.1971837591223663D+05,
     &             -.8902978767070678D+06,.5310411010968522D+08,
     &             -.4043620325107754D+10,.3827011346598605D+12,
     &             -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           CT1=Z1-0.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 50 K=1,K0
50            CP0=CP0+A(K)*Z1**(-2*K)
           CQ0=-0.125D0/Z1
           DO 55 K=1,K0
55            CQ0=CQ0+B(K)*Z1**(-2*K-1)
           CU=CDSQRT(RP2/Z1)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CT2=Z1-0.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 60 K=1,K0
60            CP1=CP1+A1(K)*Z1**(-2*K)
           CQ1=0.375D0/Z1
           DO 65 K=1,K0
65            CQ1=CQ1+B1(K)*Z1**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
        ENDIF
        IF (REAL(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBY0=CBY0-2.0D0*CI*CBJ0
           IF (DIMAG(Z).GT.0.0) CBY0=CBY0+2.0D0*CI*CBJ0
           IF (DIMAG(Z).LT.0.0) CBY1=-(CBY1-2.0D0*CI*CBJ1)
           IF (DIMAG(Z).GT.0.0) CBY1=-(CBY1+2.0D0*CI*CBJ1)
           CBJ1=-CBJ1
        ENDIF
        CDJ0=-CBJ1
        CDJ1=CBJ0-1.0D0/Z*CBJ1
        CDY0=-CBY1
        CDY1=CBY0-1.0D0/Z*CBY1
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
