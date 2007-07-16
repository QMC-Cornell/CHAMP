        PROGRAM MFCSZO

C       ===========================================================
C       Purpose : This program computes the complex zeros of the
C                 Fresnel integral C(z) or S(z) using subroutine
C                 FCSZO
C       Input :   KF  --- Function code
C                         KF=1 for C(z) or KF=2 for S(z)
C                 NT  --- Total number of zeros
C       Output:   ZO(L) --- L-th zero of C(z) or S(z)
C       Example:  NT=10

C       n     Complex zeros of C(z)        Complex zeros of S(z)
C     ------------------------------------------------------------
C       1    1.7436675 + i .30573506      2.0092570 + i .28854790
C       2    2.6514596 + i .25290396      2.8334772 + i .24428524
C       3    3.3203593 + i .22395346      3.4675331 + i .21849268
C       4    3.8757345 + i .20474747      4.0025782 + i .20085103
C       5    4.3610635 + i .19066973      4.4741893 + i .18768859
C       6    4.7976077 + i .17970801      4.9006784 + i .17732036
C       7    5.1976532 + i .17081930      5.2929467 + i .16884418
C       8    5.5690602 + i .16339854      5.6581068 + i .16172492
C       9    5.9172173 + i .15706585      6.0011034 + i .15562108
C      10    6.2460098 + i .15156826      6.3255396 + i .15030246
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION ZO(100)
        WRITE(*,*)'Please Enter KF and NT '
        READ(*,*)KF,NT
        WRITE(*,20)KF,NT
        WRITE(*,*)' *****     Please Wait !     *****'
        CALL FCSZO(KF,NT,ZO)
        WRITE(*,*)
        IF (KF.EQ.1) WRITE(*,*)'  n        Complex zeros of C(z)'
        IF (KF.EQ.2) WRITE(*,*)'  n        Complex zeros of S(z)'
        WRITE(*,*)'-----------------------------------'
        DO 10 I=1,NT
10         WRITE(*,30) I,ZO(I)
20      FORMAT(2X,'KF=',I2,',     ','NT=',I3)
30      FORMAT(1X,I3,F13.8,2X,2H+i,F13.8)
        END


        SUBROUTINE FCSZO(KF,NT,ZO)

C       ===============================================================
C       Purpose: Compute the complex zeros of Fresnel integral C(z)
C                or S(z) using modified Newton's iteration method
C       Input :  KF  --- Function code
C                        KF=1 for C(z) or KF=2 for S(z)
C                NT  --- Total number of zeros
C       Output:  ZO(L) --- L-th zero of C(z) or S(z)
C       Routines called:
C            (1) CFC for computing Fresnel integral C(z)
C            (2) CFS for computing Fresnel integral S(z)
C       ==============================================================

        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION ZO(NT)
        PI=3.141592653589793D0
        DO 35 NR=1,NT
           IF (KF.EQ.1) PSQ=DSQRT(4.0D0*NR-1.0D0)
           IF (KF.EQ.2) PSQ=2.0D0*NR**(0.5)
           PX=PSQ-DLOG(PI*PSQ)/(PI*PI*PSQ**3.0)
           PY=DLOG(PI*PSQ)/(PI*PSQ)
           Z=CMPLX(PX,PY)
           IF (KF.EQ.2) THEN
              IF (NR.EQ.2) Z=(2.8334,0.2443)
              IF (NR.EQ.3) Z=(3.4674,0.2185)
              IF (NR.EQ.4) Z=(4.0025,0.2008)
           ENDIF
           IT=0
15         IT=IT+1
           IF (KF.EQ.1) CALL CFC(Z,ZF,ZD)
           IF (KF.EQ.2) CALL CFS(Z,ZF,ZD)
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
30            ZQ=ZQ+ZW
           ZGD=(ZD-ZQ*ZFD)/ZP
           Z=Z-ZFD/ZGD
           W0=W
           W=CDABS(Z)
           IF (IT.LE.50.AND.DABS((W-W0)/W).GT.1.0D-12) GO TO 15
35         ZO(NR)=Z
        RETURN
        END


        SUBROUTINE CFC(Z,ZF,ZD)

C       =========================================================
C       Purpose: Compute complex Fresnel integral C(z) and C'(z)
C       Input :  z --- Argument of C(z)
C       Output:  ZF --- C(z)
C                ZD --- C'(z)
C       =========================================================

        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,S,Z)
        EPS=1.0D-14
        PI=3.141592653589793D0
        W0=CDABS(Z)
        ZP=0.5D0*PI*Z*Z
        ZP2=ZP*ZP
        Z0=(0.0D0,0.0D0)
        IF (Z.EQ.Z0) THEN
           C=Z0
        ELSE IF (W0.LE.2.5) THEN
           CR=Z
           C=CR
           DO 10 K=1,80
              CR=-.5D0*CR*(4.0D0*K-3.0D0)/K/(2.0D0*K-1.0D0)
     &          /(4.0D0*K+1.0D0)*ZP2
              C=C+CR
              WA=CDABS(C)
              IF (DABS((WA-WA0)/WA).LT.EPS.AND.K.GT.10) GO TO 30
10            WA0=WA
        ELSE IF (W0.GT.2.5.AND.W0.LT.4.5) THEN
           M=85
           C=Z0
           CF1=Z0
           CF0=(1.0D-100,0.0D0)
           DO 15 K=M,0,-1
              CF=(2.0D0*K+3.0D0)*CF0/ZP-CF1
              IF (K.EQ.INT(K/2)*2) C=C+CF
              CF1=CF0
15            CF0=CF
           C=CDSQRT(2.0D0/(PI*ZP))*CDSIN(ZP)/CF*C
        ELSE
           CR=(1.0D0,0.0D0)
           CF=(1.0D0,0.0D0)
           DO 20 K=1,20
              CR=-.25D0*CR*(4.0D0*K-1.0D0)*(4.0D0*K-3.0D0)/ZP2
20            CF=CF+CR
           CR=1.0D0/(PI*Z*Z)
           CG=CR
           DO 25 K=1,12
              CR=-.25D0*CR*(4.0D0*K+1.0D0)*(4.0D0*K-1.0D0)/ZP2
25            CG=CG+CR
           C=.5D0+(CF*CDSIN(ZP)-CG*CDCOS(ZP))/(PI*Z)
        ENDIF
30      ZF=C
        ZD=CDCOS(0.5*PI*Z*Z)
        RETURN
        END


        SUBROUTINE CFS(Z,ZF,ZD)

C       =========================================================
C       Purpose: Compute complex Fresnel Integral S(z) and S'(z)
C       Input :  z  --- Argument of S(z)
C       Output:  ZF --- S(z)
C                ZD --- S'(z)
C       =========================================================

        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,S,Z)
        EPS=1.0D-14
        PI=3.141592653589793D0
        W0=CDABS(Z)
        ZP=0.5D0*PI*Z*Z
        ZP2=ZP*ZP
        Z0=(0.0D0,0.0D0)
        IF (Z.EQ.Z0) THEN
           S=Z0
        ELSE IF (W0.LE.2.5) THEN
           S=Z*ZP/3.0D0
           CR=S
           DO 10 K=1,80
              CR=-.5D0*CR*(4.0D0*K-1.0D0)/K/(2.0D0*K+1.0D0)
     &          /(4.0D0*K+3.0D0)*ZP2
              S=S+CR
              WB=CDABS(S)
              IF (DABS(WB-WB0).LT.EPS.AND.K.GT.10) GO TO 30
10            WB0=WB
        ELSE IF (W0.GT.2.5.AND.W0.LT.4.5) THEN
           M=85
           S=Z0
           CF1=Z0
           CF0=(1.0D-100,0.0D0)
           DO 15 K=M,0,-1
              CF=(2.0D0*K+3.0D0)*CF0/ZP-CF1
              IF (K.NE.INT(K/2)*2) S=S+CF
              CF1=CF0
15            CF0=CF
           S=CDSQRT(2.0D0/(PI*ZP))*CDSIN(ZP)/CF*S
        ELSE
           CR=(1.0D0,0.0D0)
           CF=(1.0D0,0.0D0)
           DO 20 K=1,20
              CR=-.25D0*CR*(4.0D0*K-1.0D0)*(4.0D0*K-3.0D0)/ZP2
20            CF=CF+CR
           CR=1.0D0/(PI*Z*Z)
           CG=CR
           DO 25 K=1,12
              CR=-.25D0*CR*(4.0D0*K+1.0D0)*(4.0D0*K-1.0D0)/ZP2
25            CG=CG+CR
           S=.5D0-(CF*CDCOS(ZP)+CG*CDSIN(ZP))/(PI*Z)
        ENDIF
30      ZF=S
        ZD=CDSIN(0.5*PI*Z*Z)
        RETURN
        END
