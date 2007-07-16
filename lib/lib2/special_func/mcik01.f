        PROGRAM MCIK01

C       =============================================================
C       Purpose: This program computes the modified Bessel functions
C                I0(z), I1(z), K0(z), K1(z), and their derivatives
C                for a complex argument using subroutine CIK01
C       Input :  z --- Complex argument
C       Output:  CBI0 --- I0(z)
C                CDI0 --- I0'(z)
C                CBI1 --- I1(z)
C                CDI1 --- I1'(z)
C                CBK0 --- K0(z)
C                CDK0 --- K0'(z)
C                CBK1 --- K1(z)
C                CDK1 --- K1'(z)
C       Example: z = 20.0 + i 10.0

C     n      Re[In(z)]      Im[In(z)]      Re[In'(z)]     Im[In'(z)]
C    -----------------------------------------------------------------
C     0   -.38773811D+08 -.13750292D+08 -.37852037D+08 -.13869150D+08
C     1   -.37852037D+08 -.13869150D+08 -.36982347D+08 -.13952566D+08

C     n      Re[Kn(z)]      Im[Kn(z)]      Re[Kn'(z)]     Im[Kn'(z)]
C    -----------------------------------------------------------------
C     0   -.37692389D-09  .39171613D-09  .38056380D-09 -.40319029D-09
C     1   -.38056380D-09  .40319029D-09  .38408264D-09 -.41545502D-09
C       =============================================================

        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        WRITE(*,*)'  Please enter x and y (z=x+iy) '
        READ(*,*)X,Y
        Z=CMPLX(X,Y)
        WRITE(*,30)X,Y
        CALL CIK01(Z,CBI0,CDI0,CBI1,CDI1,CBK0,CDK0,CBK1,CDK1)
        WRITE(*,*)
        WRITE(*,*)'  n      Re[In(z)]      Im[In(z)]',
     &            '      Re[In''(z)]     Im[In''(z)]'
        WRITE(*,*)' -------------------------------',
     &            '----------------------------------'
        WRITE(*,10)CBI0,CDI0
        WRITE(*,20)CBI1,CDI1
        WRITE(*,*)
        WRITE(*,*)'  n      Re[Kn(z)]      Im[Kn(z)]',
     &            '      Re[Kn''(z)]     Im[Kn''(z)]'
        WRITE(*,*)' -------------------------------',
     &            '----------------------------------'
        WRITE(*,10)CBK0,CDK0
        WRITE(*,20)CBK1,CDK1
10      FORMAT(3X,'0',2X,4D15.7)
20      FORMAT(3X,'1',2X,4D15.7)
30      FORMAT(3X,3Hz =,F7.2,' + i',F7.2)
        END


        SUBROUTINE CIK01(Z,CBI0,CDI0,CBI1,CDI1,CBK0,CDK0,CBK1,CDK1)

C       ==========================================================
C       Purpose: Compute modified Bessel functions I0(z), I1(z),
C                K0(z), K1(z), and their derivatives for a
C                complex argument
C       Input :  z --- Complex argument
C       Output:  CBI0 --- I0(z)
C                CDI0 --- I0'(z)
C                CBI1 --- I1(z)
C                CDI1 --- I1'(z)
C                CBK0 --- K0(z)
C                CDK0 --- K0'(z)
C                CBK1 --- K1(z)
C                CDK1 --- K1'(z)
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION A(12),B(12),A1(10)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBI0=(1.0D0,0.0D0)
           CBI1=(0.0D0,0.0D0)
           CDI0=(0.0D0,0.0D0)
           CDI1=(0.5D0,0.0D0)
           CBK0=(1.0D+300,0.0D0)
           CBK1=(1.0D+300,0.0D0)
           CDK0=-(1.0D+300,0.0D0)
           CDK1=-(1.0D+300,0.0D0)
           RETURN
        ENDIF
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LE.18.0) THEN
           CBI0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,50
              CR=0.25D0*CR*Z2/(K*K)
              CBI0=CBI0+CR
              IF (CDABS(CR/CBI0).LT.1.0D-15) GO TO 15
10         CONTINUE
15         CBI1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,50
              CR=0.25D0*CR*Z2/(K*(K+1))
              CBI1=CBI1+CR
              IF (CDABS(CR/CBI1).LT.1.0D-15) GO TO 25
20         CONTINUE
25         CBI1=0.5D0*Z1*CBI1
        ELSE
           DATA A/0.125D0,7.03125D-2,
     &            7.32421875D-2,1.1215209960938D-1,
     &            2.2710800170898D-1,5.7250142097473D-1,
     &            1.7277275025845D0,6.0740420012735D0,
     &            2.4380529699556D01,1.1001714026925D02,
     &            5.5133589612202D02,3.0380905109224D03/
           DATA B/-0.375D0,-1.171875D-1,
     &            -1.025390625D-1,-1.4419555664063D-1,
     &            -2.7757644653320D-1,-6.7659258842468D-1,
     &            -1.9935317337513D0,-6.8839142681099D0,
     &            -2.7248827311269D01,-1.2159789187654D02,
     &            -6.0384407670507D02,-3.3022722944809D03/
           K0=12
           IF (A0.GE.35.0) K0=9
           IF (A0.GE.50.0) K0=7
           CA=CDEXP(Z1)/CDSQRT(2.0D0*PI*Z1)
           CBI0=(1.0D0,0.0D0)
           ZR=1.0D0/Z1
           DO 30 K=1,K0
30            CBI0=CBI0+A(K)*ZR**K
           CBI0=CA*CBI0
           CBI1=(1.0D0,0.0D0)
           DO 35 K=1,K0
35            CBI1=CBI1+B(K)*ZR**K
           CBI1=CA*CBI1
        ENDIF
        IF (A0.LE.9.0) THEN
           CS=(0.0D0,0.0D0)
           CT=-CDLOG(0.5D0*Z1)-0.5772156649015329D0
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           DO 40 K=1,50
              W0=W0+1.0D0/K
              CR=0.25D0*CR/(K*K)*Z2
              CS=CS+CR*(W0+CT)
              IF (CDABS((CS-CW)/CS).LT.1.0D-15) GO TO 45
40            CW=CS
45         CBK0=CT+CS
        ELSE
           DATA A1/0.125D0,0.2109375D0,
     &             1.0986328125D0,1.1775970458984D01,
     &             2.1461706161499D02,5.9511522710323D03,
     &             2.3347645606175D05,1.2312234987631D07,
     &             8.401390346421D08,7.2031420482627D10/
           CB=0.5D0/Z1
           ZR2=1.0D0/Z2
           CBK0=(1.0D0,0.0D0)
           DO 50 K=1,10
50            CBK0=CBK0+A1(K)*ZR2**K
           CBK0=CB*CBK0/CBI0
        ENDIF
        CBK1=(1.0D0/Z1-CBI1*CBK0)/CBI0
        IF (REAL(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBK0=CBK0+CI*PI*CBI0
           IF (DIMAG(Z).GT.0.0) CBK0=CBK0-CI*PI*CBI0
           IF (DIMAG(Z).LT.0.0) CBK1=-CBK1+CI*PI*CBI1
           IF (DIMAG(Z).GT.0.0) CBK1=-CBK1-CI*PI*CBI1
           CBI1=-CBI1
        ENDIF
        CDI0=CBI1
        CDI1=CBI0-1.0D0/Z*CBI1
        CDK0=-CBK1
        CDK1=-CBK0-1.0D0/Z*CBK1
        RETURN
        END
