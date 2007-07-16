        PROGRAM MPBWA

C       ============================================================
C       Purpose: This program computes the parabolic cylinder
C                functions W(a,ñx) and their derivatives using
C                subroutine PBWA
C       Input  : a --- Parameter  ( 0 ó |a| ó 5 )
C                x --- Argument of W(a,ñx)  ( 0 ó |x| ó 5 )
C       Output : W1F --- W(a,x)
C                W1D --- W'(a,x)
C                W2F --- W(a,-x)
C                W2D --- W'(a,-x)
C       Example: x = 5.0
C                 a      W(a,x)     W'(a,x)    W(a,-x)   W'(a,-x)
C              ----------------------------------------------------
C                0.5   .1871153    .1915744  -.8556585   4.4682493
C                1.5  -.0215853    .0899870 -8.8586002  -9.3971967
C                0.0   .3009549   -.7148233   .6599634   1.7552224
C               -0.5  -.1934088  -1.3474400   .6448148   -.6781011
C               -1.5  -.5266539    .8219516  -.2822774  -1.4582283
C               -5.0   .0893618  -1.8118641   .5386084    .2698553
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter a and x '
        READ(*,*)A,X
        WRITE(*,10)A,X
        WRITE(*,*)
        WRITE(*,*)'   a       W(a,x)          W''(a,x)',
     &            '         W(a,-x)         W''(a,-x)'
        WRITE(*,*)' -----------------------------------------',
     &            '----------------------------'
        CALL PBWA(A,X,W1F,W1D,W2F,W2D)
        WRITE(*,20)A,W1F,W1D,W2F,W2D
10      FORMAT(1X,'a=',F5.1,3X,'x=',F5.1)
20      FORMAT(1X,F5.1,4D16.8)
        END


        SUBROUTINE PBWA(A,X,W1F,W1D,W2F,W2D)

C       ======================================================
C       Purpose: Compute parabolic cylinder functions W(a,ñx)
C                and their derivatives
C       Input  : a --- Parameter  ( 0 ó |a| ó 5 )
C                x --- Argument of W(a,ñx)  ( 0 ó |x| ó 5 )
C       Output : W1F --- W(a,x)
C                W1D --- W'(a,x)
C                W2F --- W(a,-x)
C                W2D --- W'(a,-x)
C       Routine called:
C               CGAMA for computing complex gamma function
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION H(100),D(100)
        EPS=1.0D-15
        P0=0.59460355750136D0
        IF (A.EQ.0.0D0) THEN
           G1=3.625609908222D0
           G2=1.225416702465D0
        ELSE
           X1=0.25D0
           Y1=0.5D0*A
           CALL CGAMA(X1,Y1,1,UGR,UGI)
           G1=DSQRT(UGR*UGR+UGI*UGI)
           X2=0.75D0
           CALL CGAMA(X2,Y1,1,VGR,VGI)
           G2=DSQRT(VGR*VGR+VGI*VGI)
        ENDIF
        F1=DSQRT(G1/G2)
        F2=DSQRT(2.0D0*G2/G1)
        H0=1.0D0
        H1=A
        H(1)=A
        DO 10 L1=4,200,2
           M=L1/2
           HL=A*H1-0.25D0*(L1-2.0D0)*(L1-3.0D0)*H0
           H(M)=HL
           H0=H1
10         H1=HL
        Y1F=1.0D0
        R=1.0D0
        DO 15 K=1,100
           R=0.5D0*R*X*X/(K*(2.0D0*K-1.0D0))
           R1=H(K)*R
           Y1F=Y1F+R1
           IF (DABS(R1/Y1F).LE.EPS.AND.K.GT.30) GO TO 20
15      CONTINUE
20      Y1D=A
        R=1.0D0
        DO 25 K=1,100
           R=0.5D0*R*X*X/(K*(2.0D0*K+1.0D0))
           R1=H(K+1)*R
           Y1D=Y1D+R1
           IF (DABS(R1/Y1D).LE.EPS.AND.K.GT.30) GO TO 30
25      CONTINUE
30      Y1D=X*Y1D
        D1=1.0D0
        D2=A
        D(1)=1.0D0
        D(2)=A
        DO 40 L2=5,160,2
           M=(L2+1)/2
           DL=A*D2-0.25D0*(L2-2.0D0)*(L2-3.0D0)*D1
           D(M)=DL
           D1=D2
40         D2=DL
        Y2F=1.0D0
        R=1.0D0
        DO 45 K=1,100
           R=0.5D0*R*X*X/(K*(2.0D0*K+1.0D0))
           R1=D(K+1)*R
           Y2F=Y2F+R1
           IF (DABS(R1/Y2F).LE.EPS.AND.K.GT.30) GO TO 50
45      CONTINUE
50      Y2F=X*Y2F
        Y2D=1.0D0
        R=1.0D0
        DO 55 K=1,100
           R=0.5D0*R*X*X/(K*(2.0D0*K-1.0D0))
           R1=D(K+1)*R
           Y2D=Y2D+R1
           IF (DABS(R1/Y2D).LE.EPS.AND.K.GT.30) GO TO 60
55      CONTINUE
60      W1F=P0*(F1*Y1F-F2*Y2F)
        W2F=P0*(F1*Y1F+F2*Y2F)
        W1D=P0*(F1*Y1D-F2*Y2D)
        W2D=P0*(F1*Y1D+F2*Y2D)
        RETURN
        END


        SUBROUTINE CGAMA(X,Y,KF,GR,GI)

C       =========================================================
C       Purpose: Compute complex gamma function â(z) or Ln[â(z)]
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                KF --- Function code
C                       KF=0 for Ln[â(z)]
C                       KF=1 for â(z)
C       Output:  GR --- Real part of Ln[â(z)] or â(z)
C                GI --- Imaginary part of Ln[â(z)] or â(z)
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        PI=3.141592653589793D0
        DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/
        IF (Y.EQ.0.0D0.AND.X.EQ.INT(X).AND.X.LE.0.0D0) THEN
           GR=1.0D+300
           GI=0.0D0
           RETURN
        ELSE IF (X.LT.0.0D0) THEN
           X1=X
           Y1=Y
           X=-X
           Y=-Y
        ENDIF
        X0=X
        IF (X.LE.7.0) THEN
           NA=INT(7-X)
           X0=X+NA
        ENDIF
        Z1=DSQRT(X0*X0+Y*Y)
        TH=DATAN(Y/X0)
        GR=(X0-.5D0)*DLOG(Z1)-TH*Y-X0+0.5D0*DLOG(2.0D0*PI)
        GI=TH*(X0-0.5D0)+Y*DLOG(Z1)-Y
        DO 10 K=1,10
           T=Z1**(1-2*K)
           GR=GR+A(K)*T*DCOS((2.0D0*K-1.0D0)*TH)
10         GI=GI-A(K)*T*DSIN((2.0D0*K-1.0D0)*TH)
        IF (X.LE.7.0) THEN
           GR1=0.0D0
           GI1=0.0D0
           DO 15 J=0,NA-1
              GR1=GR1+.5D0*DLOG((X+J)**2+Y*Y)
15            GI1=GI1+DATAN(Y/(X+J))
           GR=GR-GR1
           GI=GI-GI1
        ENDIF
        IF (X1.LT.0.0D0) THEN
           Z1=DSQRT(X*X+Y*Y)
           TH1=DATAN(Y/X)
           SR=-DSIN(PI*X)*DCOSH(PI*Y)
           SI=-DCOS(PI*X)*DSINH(PI*Y)
           Z2=DSQRT(SR*SR+SI*SI)
           TH2=DATAN(SI/SR)
           IF (SR.LT.0.0D0) TH2=PI+TH2
           GR=DLOG(PI/(Z1*Z2))-GR
           GI=-TH1-TH2-GI
           X=X1
           Y=Y1
        ENDIF
        IF (KF.EQ.1) THEN
           G0=DEXP(GR)
           GR=G0*DCOS(GI)
           GI=G0*DSIN(GI)
        ENDIF
        RETURN
        END
