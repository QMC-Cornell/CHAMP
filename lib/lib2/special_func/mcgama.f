        PROGRAM MCGAMA

C       ==========================================================
C       Purpose: This program computes the gamma function â(z)
C                or ln[â(z)] for a complex argument using
C                subroutine CGAMA
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                KF --- Function code
C                       KF=0 for ln[â(z)]
C                       KF=1 for â(z)
C       Output:  GR --- Real part of ln[â(z)] or â(z)
C                GI --- Imaginary part of ln[â(z)] or â(z)
C       Examples:

C         x         y           Re[â(z)]           Im[â(z)]
C       --------------------------------------------------------
C        2.50      5.00     .2267360319D-01    -.1172284404D-01
C        5.00     10.00     .1327696517D-01     .3639011746D-02
C        2.50     -5.00     .2267360319D-01     .1172284404D-01
C        5.00    -10.00     .1327696517D-01    -.3639011746D-02

C         x         y          Re[lnâ(z)]         Im[lnâ(z)]
C      ---------------------------------------------------------
C        2.50      5.00    -.3668103262D+01     .5806009801D+01
C        5.00     10.00    -.4285507444D+01     .1911707090D+02
C        2.50     -5.00    -.3668103262D+01    -.5806009801D+01
C        5.00    -10.00    -.4285507444D+01    -.1911707090D+02
C       ==========================================================

        DOUBLE PRECISION X,Y,GR,GI
        WRITE(*,*)'  Please enter KF, x and y'
        READ(*,*)KF,X,Y
        WRITE(*,*)
        IF (KF.EQ.1) THEN
            WRITE(*,*)'       x         y           Re[â(z)]',
     &                '           Im[â(z)]'
        ELSE
            WRITE(*,*)'       x         y          Re[lnâ(z)]',
     &                '         Im[lnâ(z)]'
        ENDIF
        WRITE(*,*)'    ------------------------------------',
     &            '---------------------'
        CALL CGAMA(X,Y,KF,GR,GI)
        WRITE(*,10)X,Y,GR,GI
10      FORMAT(1X,2F10.2,2D20.10)
        END


        SUBROUTINE CGAMA(X,Y,KF,GR,GI)

C       =========================================================
C       Purpose: Compute the gamma function â(z) or ln[â(z)]
C                for a complex argument
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                KF --- Function code
C                       KF=0 for ln[â(z)]
C                       KF=1 for â(z)
C       Output:  GR --- Real part of ln[â(z)] or â(z)
C                GI --- Imaginary part of ln[â(z)] or â(z)
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
