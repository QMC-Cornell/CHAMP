        PROGRAM MCYZO

C       ===========================================================
C       Purpose : This program evaluates the complex zeros of
C                 Y0(z), Y0'(z), Y1(z) and Y1'(z), and their
C                 associated values at the zeros using the
C                 modified Newton's iteration method
C       Input:    NT --- Total number of roots/zeros
C                 KF --- Function choice code
C                        KF=0 for  Y0(z) & Y1(z0)
C                        KF=1 for  Y1(z) & Y0(z1)
C                        KF=2 for  Y1'(z) & Y1(z1')
C                 KC --- Choice code
C                        KC=0 for complex roots
C                        KC=1 for real roots
C       Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z)
C                 ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z)
C                           at the L-th zero
C       Examples: NT = 5

C   No.      z0, Zeros of Y0(z)                Y1(z0)
C  -----------------------------------------------------------------
C    1   -2.403016632 + i .5398823130   .1007476893 - i .8819677101
C    2   -5.519876702 + i .5471800106  -.0292464182 + i .5871695027
C    3   -8.653672403 + i .5484120673   .0149080637 - i .4694587524
C    4  -11.791512030 + i .5488191184  -.0093736817 + i .4023045429
C    5  -14.930906564 + i .5490008289   .0065788031 - i .3575673214

C   No.      z1, Zeros of Y1(z)                 Y0(z1)
C  -----------------------------------------------------------------
C    1    -.502743273 + i .7862437145  -.4595276847 + i1.3171019361
C    2   -3.833535193 + i .5623565382   .0483019087 - i .6925128842
C    3   -7.015903683 + i .5533930459  -.0201269494 + i .5186425332
C    4  -10.173573834 + i .5512733877   .0116140017 - i .4320329636
C    5  -13.323739307 + i .5504585830  -.0077719300 + i .3779698048

C   No.      z1', Zeros of Y1'(z)                Y1(z1')
C   ----------------------------------------------------------------
C    1     .576785129 + i .9039847922  -.7634970879 + i .5892448647
C    2   -1.940477342 - i .7211859189   .1620640057 + i .9520278864
C    3   -5.333478617 - i .5672196368  -.0317940081 - i .5968536736
C    4   -8.536768577 - i .5560607040   .0154177166 + i .4726011652
C    5  -11.706175219 - i .5528590607  -.0095443768 - i .4037533396
C       ============================================================

        IMPLICIT COMPLEX *16 (Z)
        DIMENSION ZO(50),ZV(50)
        WRITE(*,*)'Please Enter NT, KF and KC'
        WRITE(*,*)'  NT --- Total number of the roots'
        WRITE(*,*)'  KF  --- Function choice code'
        WRITE(*,*)'          KF=0 for Y0(z) & Y1(z0)'
        WRITE(*,*)'          KF=1 for Y1(z) & Y0(z1)'
        WRITE(*,*)'          KF=2 for Y1''(z) & Y1(z1'')'
        WRITE(*,*)'  KC  --- Choice code'
        WRITE(*,*)'          KC=0 for complex roots'
        WRITE(*,*)'          KC=1 for real roots'
        READ(*,*)NT,KF,KC
        WRITE(*,20)NT,KF,KC
        WRITE(*,*)
        WRITE(*,15)
        CALL CYZO(NT,KF,KC,ZO,ZV)
        WRITE(*,*)
        IF (KF.EQ.0) THEN
           WRITE(*,*)' No.          z0, Zeros of Y0(z)',
     &               '                 Y1(z0)'
        ELSE IF (KF.EQ.1) THEN
           WRITE(*,*)' No.          z1, Zeros of Y1(z)',
     &               '                 Y0(z1)'
        ELSE IF (KF.EQ.2) THEN
           WRITE(*,*)' No.        z1'', Zeros of Y1''(z)',
     &               '                Y1(z1'')'
        ENDIF
        WRITE(*,*)'--------------------------------------',
     &            '----------------------------'
        DO 10 I=1,NT
10         WRITE(*,25)I,ZO(I),ZV(I)
15      FORMAT(20X,'*****    Please Wait !    *****')
20      FORMAT(2X,'NT=',I3,',  ','KF=',I3,',  ','KC=',I3)
25      FORMAT(1X,I3,2X,F15.9,3F15.10)
        END


        SUBROUTINE CYZO(NT,KF,KC,ZO,ZV)

C       ===========================================================
C       Purpose : Compute the complex zeros of Y0(z), Y1(z) and
C                 Y1'(z), and their associated values at the zeros
C                 using the modified Newton's iteration method
C       Input:    NT --- Total number of zeros/roots
C                 KF --- Function choice code
C                        KF=0 for  Y0(z) & Y1(z0)
C                        KF=1 for  Y1(z) & Y0(z1)
C                        KF=2 for  Y1'(z) & Y1(z1')
C                 KC --- Choice code
C                        KC=0 for complex roots
C                        KC=1 for real roots
C       Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z)
C                 ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z)
C                           at the L-th zero
C       Routine called: CY01 for computing Y0(z) and Y1(z), and
C                       their derivatives
C       ===========================================================
        IMPLICIT DOUBLE PRECISION (H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION ZO(NT),ZV(NT)
        IF (KC.EQ.0) THEN
           X=-2.4D0
           Y=0.54D0
           H=3.14
        ELSE IF (KC.EQ.1) THEN
           X=0.89
           Y=0.0
           H=-3.14
        ENDIF
        IF (KF.EQ.1) X=-0.503
        IF (KF.EQ.2) X=0.577
        ZERO=CMPLX(X,Y)
        Z=ZERO
        DO 35 NR=1,NT
10         IF (NR.NE.1) Z=ZO(NR-1)-H
           IT=0
15         IT=IT+1
           CALL CY01(KF,Z,ZF,ZD)
           ZP=(1.0D0,0.0D0)
           DO 20 I=1,NR-1
20            ZP=ZP*(Z-ZO(I))
           ZFD=ZF/ZP
           ZQ=(0.0D0,0.0D0)
           DO 30 I=1,NR-1
              ZW=(1.0D0,0.0D0)
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 ZW=ZW*(Z-ZO(J))
25            CONTINUE
              ZQ=ZQ+ZW
30         CONTINUE
           ZGD=(ZD-ZQ*ZFD)/ZP
           Z=Z-ZFD/ZGD
           W0=W
           W=CDABS(Z)
           IF (IT.LE.50.AND.DABS((W-W0)/W).GT.1.0D-12) GO TO 15
           ZO(NR)=Z
35      CONTINUE
        DO 40 I=1,NT
           Z=ZO(I)
           IF (KF.EQ.0.OR.KF.EQ.2) THEN
              CALL CY01(1,Z,ZF,ZD)
              ZV(I)=ZF
           ELSE IF (KF.EQ.1) THEN
              CALL CY01(0,Z,ZF,ZD)
              ZV(I)=ZF
           ENDIF
40      CONTINUE
        RETURN
        END


        SUBROUTINE CY01(KF,Z,ZF,ZD)

C       ===========================================================
C       Purpose: Compute complex Bessel functions Y0(z), Y1(z)
C                and their derivatives
C       Input :  z  --- Complex argument of Yn(z) ( n=0,1 )
C                KF --- Function choice code
C                    KF=0 for ZF=Y0(z) and ZD=Y0'(z)
C                    KF=1 for ZF=Y1(z) and ZD=Y1'(z)
C                    KF=2 for ZF=Y1'(z) and ZD=Y1''(z)
C       Output:  ZF --- Y0(z) or Y1(z) or Y1'(z)
C                ZD --- Y0'(z) or Y1'(z) or Y1''(z)
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
           CBY0=-(1.0D300,0.0D0)
           CBY1=-(1.0D300,0.0D0)
           CDY0=(1.0D300,0.0D0)
           CDY1=(1.0D300,0.0D0)
           GO TO 70
        ENDIF
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CBJ0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,40
              CR=-0.25D0*CR*Z2/(K*K)
              CBJ0=CBJ0+CR
              IF (CDABS(CR).LT.CDABS(CBJ0)*1.0D-15) GO TO 15
10         CONTINUE
15         CBJ1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+1.0D0))
              CBJ1=CBJ1+CR
              IF (CDABS(CR).LT.CDABS(CBJ1)*1.0D-15) GO TO 25
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
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 35
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
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 45
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
           CT1=Z1-.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 50 K=1,K0
50            CP0=CP0+A(K)*Z1**(-2*K)
           CQ0=-0.125D0/Z1
           DO 55 K=1,K0
55            CQ0=CQ0+B(K)*Z1**(-2*K-1)
           CU=CDSQRT(RP2/Z1)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CT2=Z1-.75D0*PI
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
        CDY0=-CBY1
        CDY1=CBY0-1.0D0/Z*CBY1
70      IF (KF.EQ.0) THEN
           ZF=CBY0
           ZD=CDY0
        ELSE IF (KF.EQ.1) THEN
           ZF=CBY1
           ZD=CDY1
        ELSE IF (KF.EQ.2) THEN
           ZF=CDY1
           ZD=-CDY1/Z-(1.0D0-1.0D0/(Z*Z))*CBY1
        ENDIF
        RETURN
        END

