        PROGRAM MAIRYB

C       ============================================================
C       Purpose: This program computes Airy functions and their
C                derivatives using subroutine AIRYB
C       Input:   x  --- Argument of Airy function
C       Output:  AI --- Ai(x)
C                BI --- Bi(x)
C                AD --- Ai'(x)
C                BD --- Bi'(x)
C       Example:

C   x       Ai(x)          Bi(x)          Ai'(x)         Bi'(x)
C  ----------------------------------------------------------------
C   0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00
C  10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10
C  20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26
C  30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48

C   x       Ai(-x)         Bi(-x)         Ai'(-x)        Bi'(-x)
C  ----------------------------------------------------------------
C   0       .35502805      .61492663     -.25881940      .44828836
C  10       .04024124     -.31467983      .99626504      .11941411
C  20      -.17640613     -.20013931      .89286286     -.79142903
C  30      -.08796819     -.22444694     1.22862060     -.48369473
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        CALL AIRYB(X,AI,BI,AD,BD)
        WRITE(*,30)
        WRITE(*,40)
        WRITE(*,10)X,AI,BI,AD,BD
        WRITE(*,*)
        CALL AIRYB(-X,AI,BI,AD,BD)
        WRITE(*,50)
        WRITE(*,40)
        WRITE(*,20)X,AI,BI,AD,BD
10      FORMAT(1X,F5.1,4D16.8)
20      FORMAT(1X,F5.1,4D16.8)
30      FORMAT(4X,'x',8X,'Ai(x)',11X,'Bi(x)',11X,'Ai''(x)',
     &         10X,'Bi''(x)')
40      FORMAT(2X,'----------------------------------',
     &        '-----------------------------------')
50      FORMAT(4X,'x',8X,'Ai(-x)',10X,'Bi(-x)',10X,
     &        'Ai''(-x)',9X,'Bi''(-x)')
        END


        SUBROUTINE AIRYB(X,AI,BI,AD,BD)

C       =======================================================
C       Purpose: Compute Airy functions and their derivatives
C       Input:   x  --- Argument of Airy function
C       Output:  AI --- Ai(x)
C                BI --- Bi(x)
C                AD --- Ai'(x)
C                BD --- Bi'(x)
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION CK(41),DK(41)
        EPS=1.0D-15
        PI=3.141592653589793D0
        C1=0.355028053887817D0
        C2=0.258819403792807D0
        SR3=1.732050807568877D0
        XA=DABS(X)
        XQ=DSQRT(XA)
        IF (X.GT.0.0D0) XM=5.0
        IF (X.LE.0.0D0) XM=8.0
        IF (X.EQ.0.0D0) THEN
           AI=C1
           BI=SR3*C1
           AD=-C2
           BD=SR3*C2
           RETURN
        ENDIF
        IF (XA.LE.XM) THEN
           FX=1.0D0
           R=1.0D0
           DO 10 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K-1.0D0)*X
              FX=FX+R
              IF (DABS(R).LT.DABS(FX)*EPS) GO TO 15
10         CONTINUE
15         GX=X
           R=X
           DO 20 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K+1.0D0)*X
              GX=GX+R
              IF (DABS(R).LT.DABS(GX)*EPS) GO TO 25
20         CONTINUE
25         AI=C1*FX-C2*GX
           BI=SR3*(C1*FX+C2*GX)
           DF=0.5D0*X*X
           R=DF
           DO 30 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K+2.0D0)*X
              DF=DF+R
              IF (DABS(R).LT.DABS(DF)*EPS) GO TO 35
30         CONTINUE
35         DG=1.0D0
           R=1.0D0
           DO 40 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K-2.0D0)*X
              DG=DG+R
              IF (DABS(R).LT.DABS(DG)*EPS) GO TO 45
40         CONTINUE
45         AD=C1*DF-C2*DG
           BD=SR3*(C1*DF+C2*DG)
        ELSE
           XE=XA*XQ/1.5D0
           XR1=1.0D0/XE
           XAR=1.0D0/XQ
           XF=DSQRT(XAR)
           RP=0.5641895835477563D0
           R=1.0D0
           DO 50 K=1,40
              R=R*(6.0D0*K-1.0D0)/216.0D0*(6.0D0*K-3.0D0)
     &          /K*(6.0D0*K-5.0D0)/(2.0D0*K-1.0D0)
              CK(K)=R
50            DK(K)=-(6.0D0*K+1.0D0)/(6.0D0*K-1.0D0)*CK(K)
           KM=INT(24.5-XA)
           IF (XA.LT.6.0) KM=14
           IF (XA.GT.15.0) KM=10
           IF (X.GT.0.0D0) THEN
              SAI=1.0D0
              SAD=1.0D0
              R=1.0D0
              DO 55 K=1,KM
                 R=-R*XR1
                 SAI=SAI+CK(K)*R
55               SAD=SAD+DK(K)*R
              SBI=1.0D0
              SBD=1.0D0
              R=1.0D0
              DO 60 K=1,KM
                 R=R*XR1
                 SBI=SBI+CK(K)*R
60               SBD=SBD+DK(K)*R
              XP1=DEXP(-XE)
              AI=0.5D0*RP*XF*XP1*SAI
              BI=RP*XF/XP1*SBI
              AD=-.5D0*RP/XF*XP1*SAD
              BD=RP/XF/XP1*SBD
           ELSE
              XCS=DCOS(XE+PI/4.0D0)
              XSS=DSIN(XE+PI/4.0D0)
              SSA=1.0D0
              SDA=1.0D0
              R=1.0D0
              XR2=1.0D0/(XE*XE)
              DO 65 K=1,KM
                 R=-R*XR2
                 SSA=SSA+CK(2*K)*R
65               SDA=SDA+DK(2*K)*R
              SSB=CK(1)*XR1
              SDB=DK(1)*XR1
              R=XR1
              DO 70 K=1,KM
                 R=-R*XR2
                 SSB=SSB+CK(2*K+1)*R
70               SDB=SDB+DK(2*K+1)*R
              AI=RP*XF*(XSS*SSA-XCS*SSB)
              BI=RP*XF*(XCS*SSA+XSS*SSB)
              AD=-RP/XF*(XCS*SDA+XSS*SDB)
              BD=RP/XF*(XSS*SDA-XCS*SDB)
           ENDIF
        ENDIF
        RETURN
        END
