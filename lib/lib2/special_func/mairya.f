        PROGRAM MAIRYA

C       ============================================================
C       Purpose: This program computes Airy functions and their
C                derivatives using subroutine AIRYA
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
        CALL AIRYA(X,AI,BI,AD,BD)
        WRITE(*,30)
        WRITE(*,40)
        WRITE(*,10)X,AI,BI,AD,BD
        WRITE(*,*)
        CALL AIRYA(-X,AI,BI,AD,BD)
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


        SUBROUTINE AIRYA(X,AI,BI,AD,BD)

C       ======================================================
C       Purpose: Compute Airy functions and their derivatives
C       Input:   x  --- Argument of Airy function
C       Output:  AI --- Ai(x)
C                BI --- Bi(x)
C                AD --- Ai'(x)
C                BD --- Bi'(x)
C       Routine called:
C                AJYIK for computing Jv(x), Yv(x), Iv(x) and
C                Kv(x) with v=1/3 and 2/3
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PIR=0.318309886183891D0
        C1=0.355028053887817D0
        C2=0.258819403792807D0
        SR3=1.732050807568877D0
        Z=XA**1.5/1.5D0
        XQ=DSQRT(XA)
        CALL AJYIK(Z,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2)
        IF (X.EQ.0.0D0) THEN
           AI=C1
           BI=SR3*C1
           AD=-C2
           BD=SR3*C2
        ELSE IF (X.GT.0.0D0) THEN
           AI=PIR*XQ/SR3*VK1
           BI=XQ*(PIR*VK1+2.0D0/SR3*VI1)
           AD=-XA/SR3*PIR*VK2
           BD=XA*(PIR*VK2+2.0D0/SR3*VI2)
        ELSE
           AI=0.5D0*XQ*(VJ1-VY1/SR3)
           BI=-0.5D0*XQ*(VJ1/SR3+VY1)
           AD=0.5D0*XA*(VJ2+VY2/SR3)
           BD=0.5D0*XA*(VJ2/SR3-VY2)
        ENDIF
        RETURN
        END


        SUBROUTINE AJYIK(X,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2)

C       =======================================================
C       Purpose: Compute Bessel functions Jv(x) and Yv(x),
C                and modified Bessel functions Iv(x) and
C                Kv(x), and their derivatives with v=1/3,2/3
C       Input :  x --- Argument of Jv(x),Yv(x),Iv(x) and
C                      Kv(x) ( x ò 0 )
C       Output:  VJ1 --- J1/3(x)
C                VJ2 --- J2/3(x)
C                VY1 --- Y1/3(x)
C                VY2 --- Y2/3(x)
C                VI1 --- I1/3(x)
C                VI2 --- I2/3(x)
C                VK1 --- K1/3(x)
C                VK2 --- K2/3(x)
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0D0) THEN
           VJ1=0.0D0
           VJ2=0.0D0
           VY1=-1.0D+300
           VY2=1.0D+300
           VI1=0.0D0
           VI2=0.0D0
           VK1=-1.0D+300
           VK2=-1.0D+300
           RETURN
        ENDIF
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        GP1=.892979511569249D0
        GP2=.902745292950934D0
        GN1=1.3541179394264D0
        GN2=2.678938534707747D0
        VV0=0.444444444444444D0
        UU0=1.1547005383793D0
        X2=X*X
        K0=12
        IF (X.GE.35.0) K0=10
        IF (X.GE.50.0) K0=8
        IF (X.LE.12.0) THEN
           DO 25 L=1,2
              VL=L/3.0D0
              VJL=1.0D0
              R=1.0D0
              DO 15 K=1,40
                 R=-0.25D0*R*X2/(K*(K+VL))
                 VJL=VJL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 20
15            CONTINUE
20            A0=(0.5D0*X)**VL
              IF (L.EQ.1) VJ1=A0/GP1*VJL
              IF (L.EQ.2) VJ2=A0/GP2*VJL
25         CONTINUE
        ELSE
           DO 40 L=1,2
              VV=VV0*L*L
              PX=1.0D0
              RP=1.0D0
              DO 30 K=1,K0
                 RP=-0.78125D-2*RP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &              (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*X2)
30               PX=PX+RP
              QX=1.0D0
              RQ=1.0D0
              DO 35 K=1,K0
                 RQ=-0.78125D-2*RQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &              (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*X2)
35               QX=QX+RQ
              QX=0.125D0*(VV-1.0)*QX/X
              XK=X-(0.5D0*L/3.0D0+0.25D0)*PI
              A0=DSQRT(RP2/X)
              CK=DCOS(XK)
              SK=DSIN(XK)
              IF (L.EQ.1) THEN
                 VJ1=A0*(PX*CK-QX*SK)
                 VY1=A0*(PX*SK+QX*CK)
              ELSE IF (L.EQ.2) THEN
                 VJ2=A0*(PX*CK-QX*SK)
                 VY2=A0*(PX*SK+QX*CK)
              ENDIF
40         CONTINUE
        ENDIF
        IF (X.LE.12.0D0) THEN
           DO 55 L=1,2
              VL=L/3.0D0
              VJL=1.0D0
              R=1.0D0
              DO 45 K=1,40
                 R=-0.25D0*R*X2/(K*(K-VL))
                 VJL=VJL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 50
45            CONTINUE
50            B0=(2.0D0/X)**VL
              IF (L.EQ.1) UJ1=B0*VJL/GN1
              IF (L.EQ.2) UJ2=B0*VJL/GN2
55         CONTINUE
           PV1=PI/3.0D0
           PV2=PI/1.5D0
           VY1=UU0*(VJ1*DCOS(PV1)-UJ1)
           VY2=UU0*(VJ2*DCOS(PV2)-UJ2)
        ENDIF
        IF (X.LE.18.0) THEN
           DO 70 L=1,2
              VL=L/3.0D0
              VIL=1.0D0
              R=1.0D0
              DO 60 K=1,40
                 R=0.25D0*R*X2/(K*(K+VL))
                 VIL=VIL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 65
60            CONTINUE
65            A0=(0.5D0*X)**VL
              IF (L.EQ.1) VI1=A0/GP1*VIL
              IF (L.EQ.2) VI2=A0/GP2*VIL
70         CONTINUE
        ELSE
           C0=DEXP(X)/DSQRT(2.0D0*PI*X)
           DO 80 L=1,2
              VV=VV0*L*L
              VSL=1.0D0
              R=1.0D0
              DO 75 K=1,K0
                 R=-0.125D0*R*(VV-(2.0D0*K-1.0D0)**2.0)/(K*X)
75               VSL=VSL+R
              IF (L.EQ.1) VI1=C0*VSL
              IF (L.EQ.2) VI2=C0*VSL
80         CONTINUE
        ENDIF
        IF (X.LE.9.0D0) THEN
           DO 95 L=1,2
              VL=L/3.0D0
               IF (L.EQ.1) GN=GN1
               IF (L.EQ.2) GN=GN2
               A0=(2.0D0/X)**VL/GN
               SUM=1.0D0
               R=1.0D0
               DO 85 K=1,60
                  R=0.25D0*R*X2/(K*(K-VL))
                  SUM=SUM+R
                  IF (DABS(R).LT.1.0D-15) GO TO 90
85             CONTINUE
90            IF (L.EQ.1) VK1=0.5D0*UU0*PI*(SUM*A0-VI1)
              IF (L.EQ.2) VK2=0.5D0*UU0*PI*(SUM*A0-VI2)
95         CONTINUE
        ELSE
           C0=DEXP(-X)*DSQRT(0.5D0*PI/X)
           DO 105 L=1,2
              VV=VV0*L*L
              SUM=1.0D0
              R=1.0D0
              DO 100 K=1,K0
                 R=0.125D0*R*(VV-(2.0*K-1.0)**2.0)/(K*X)
100              SUM=SUM+R
              IF (L.EQ.1) VK1=C0*SUM
              IF (L.EQ.2) VK2=C0*SUM
105        CONTINUE
        ENDIF
        RETURN
        END


