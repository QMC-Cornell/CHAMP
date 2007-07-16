        PROGRAM MFFK

C       ==============================================================
C       Purpose: This program computes the modified Fresnel integrals
C                FÒ(x) and KÒ(x) using subroutine FFK
C       Input :  x   --- Argument of FÒ(x) and KÒ(x)
C                KS  --- Sign code
C                        KS=0 for calculating F+(x) and K+(x)
C                        KS=1 for calculating F_(x) and K_(x)
C       Output:  FR  --- Re[FÒ(x)]
C                FI  --- Im[FÒ(x)]
C                FM  --- |FÒ(x)|
C                FA  --- Arg[FÒ(x)]  (Degs.)
C                GR  --- Re[KÒ(x)]
C                GI  --- Im[KÒ(x)]
C                GM  --- |KÒ(x)|
C                GA  --- Arg[KÒ(x)]  (Degs.)
C       Example:

C         x     Re[FÒ(x)]   ÒIm[FÒ(x)]    Mod[FÒ(x)]  ÒArg[FÒ(x)]
C       ----------------------------------------------------------
C        0.0    .62665707    .62665707    .88622693    45.000000
C        2.0    .16519561   -.17811942    .24293233   -47.155835
C        4.0    .03219674   -.12047678    .12470479   -75.037684
C        6.0    .08245304   -.01180212    .08329342    -8.145843
C        8.0   -.05729996    .02493542    .06249048   156.482601
C       10.0    .02553188    .04298617    .04999688    59.291561

C         x     Re[KÒ(x)]   ÒIm[KÒ(x)]    Mod[KÒ(x)]  ÒArg[KÒ(x)]
C       ----------------------------------------------------------
C        0.0    .50000000    .00000000    .50000000     0.000000
C        2.0    .10702394    .08562295    .13705989    38.661047
C        4.0    .05126306    .04818949    .07035714    43.229843
C        6.0    .03368650    .03276566    .04699328    44.206095
C        8.0    .02512396    .02473472    .03525648    44.552712
C       10.0    .02004532    .01984592    .02820772    44.713609
C       ===============================================================

        IMPLICIT DOUBLE PRECISION (F,G,X)
        WRITE(*,*)'Please enter x'
        READ(*,*) X
        WRITE(*,*)
        WRITE(*,*)'   x      Re[FÒ(x)]    ÒIm[FÒ(x)]     ',
     &             'Mod[FÒ(x)]   ÒArg[FÒ(x)]'
        WRITE(*,*)' ---------------------------------------',
     &            '-----------------------'
        CALL FFK(0,X,FR,FI,FM,FA,GR,GI,GM,GA)
        WRITE(*,10)X,FR,FI,FM,FA
        WRITE(*,*)
        WRITE(*,*)'   x      Re[KÒ(x)]    ÒIm[KÒ(x)]     ',
     &             'Mod[KÒ(x)]   ÒArg[KÒ(x)]'
        WRITE(*,*)' ---------------------------------------',
     &            '-----------------------'
        WRITE(*,10)X,GR,GI,GM,GA
10      FORMAT(1X,F5.1,3F14.8,F14.6)
        END


        SUBROUTINE FFK(KS,X,FR,FI,FM,FA,GR,GI,GM,GA)

C       =======================================================
C       Purpose: Compute modified Fresnel integrals FÒ(x)
C                and KÒ(x)
C       Input :  x   --- Argument of FÒ(x) and KÒ(x)
C                KS  --- Sign code
C                        KS=0 for calculating F+(x) and K+(x)
C                        KS=1 for calculating F_(x) and K_(x)
C       Output:  FR  --- Re[FÒ(x)]
C                FI  --- Im[FÒ(x)]
C                FM  --- |FÒ(x)|
C                FA  --- Arg[FÒ(x)]  (Degs.)
C                GR  --- Re[KÒ(x)]
C                GI  --- Im[KÒ(x)]
C                GM  --- |KÒ(x)|
C                GA  --- Arg[KÒ(x)]  (Degs.)
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        SRD= 57.29577951308233D0
        EPS=1.0D-15
        PI=3.141592653589793D0
        PP2=1.2533141373155D0
        P2P=.7978845608028654D0
        XA=DABS(X)
        X2=X*X
        X4=X2*X2
        IF (X.EQ.0.0D0) THEN
           FR=.5D0*DSQRT(0.5D0*PI)
           FI=(-1)**KS*FR
           FM=DSQRT(0.25D0*PI)
           FA=(-1)**KS*45.0D0
           GR=.5D0
           GI=0.0D0
           GM=.5D0
           GA=0.0D0
        ELSE
           IF (XA.LE.2.5D0) THEN
              XR=P2P*XA
              C1=XR
              DO 10 K=1,50
                 XR=-.5D0*XR*(4.0D0*K-3.0D0)/K/(2.0D0*K-1.0D0)
     &              /(4.0D0*K+1.0D0)*X4
                 C1=C1+XR
                 IF (DABS(XR/C1).LT.EPS) GO TO 15
10            CONTINUE
15            S1=P2P*XA*XA*XA/3.0D0
              XR=S1
              DO 20 K=1,50
                 XR=-.5D0*XR*(4.0D0*K-1.0D0)/K/(2.0D0*K+1.0D0)
     &              /(4.0D0*K+3.0D0)*X4
                 S1=S1+XR
                 IF (DABS(XR/S1).LT.EPS) GO TO 40
20            CONTINUE
           ELSE IF (XA.LT.5.5D0) THEN
              M=INT(42+1.75*X2)
              XSU=0.0D0
              XC=0.0D0
              XS=0.0D0
              XF1=0.0D0
              XF0=1D-100
              DO 25 K=M,0,-1
                 XF=(2.0D0*K+3.0D0)*XF0/X2-XF1
                 IF (K.EQ.2*INT(K/2))  THEN
                    XC=XC+XF
                 ELSE
                    XS=XS+XF
                 ENDIF
                 XSU=XSU+(2.0D0*K+1.0D0)*XF*XF
                 XF1=XF0
25               XF0=XF
              XQ=DSQRT(XSU)
              XW=P2P*XA/XQ
              C1=XC*XW
              S1=XS*XW
           ELSE
              XR=1.0D0
              XF=1.0D0
              DO 30 K=1,12
                 XR=-.25D0*XR*(4.0D0*K-1.0D0)*(4.0D0*K-3.0D0)/X4
30               XF=XF+XR
              XR=1.0D0/(2.0D0*XA*XA)
              XG=XR
              DO 35 K=1,12
                 XR=-.25D0*XR*(4.0D0*K+1.0D0)*(4.0D0*K-1.0D0)/X4
35               XG=XG+XR
              C1=.5D0+(XF*DSIN(X2)-XG*DCOS(X2))/DSQRT(2.0D0*PI)/XA
              S1=.5D0-(XF*DCOS(X2)+XG*DSIN(X2))/DSQRT(2.0D0*PI)/XA
           ENDIF
40         FR=PP2*(.5D0-C1)
           FI0=PP2*(.5D0-S1)
           FI=(-1)**KS*FI0
           FM=DSQRT(FR*FR+FI*FI)
           IF (FR.GE.0.0) THEN
              FA=SRD*DATAN(FI/FR)
           ELSE IF (FI.GT.0.0) THEN
              FA=SRD*(DATAN(FI/FR)+PI)
           ELSE IF (FI.LT.0.0) THEN
              FA=SRD*(DATAN(FI/FR)-PI)
           ENDIF
           XP=X*X+PI/4.0D0
           CS=DCOS(XP)
           SS=DSIN(XP)
           XQ2=1.0D0/DSQRT(PI)
           GR=XQ2*(FR*CS+FI0*SS)
           GI=(-1)**KS*XQ2*(FI0*CS-FR*SS)
           GM=DSQRT(GR*GR+GI*GI)
           IF (GR.GE.0.0) THEN
              GA=SRD*DATAN(GI/GR)
           ELSE IF (GI.GT.0.0) THEN
              GA=SRD*(DATAN(GI/GR)+PI)
           ELSE IF (GI.LT.0.0) THEN
              GA=SRD*(DATAN(GI/GR)-PI)
           ENDIF
           IF (X.LT.0.0D0) THEN
              FR=PP2-FR
              FI=(-1)**KS*PP2-FI
              FM=DSQRT(FR*FR+FI*FI)
              FA=SRD*DATAN(FI/FR)
              GR=DCOS(X*X)-GR
              GI=-(-1)**KS*DSIN(X*X)-GI
              GM=DSQRT(GR*GR+GI*GI)
              GA=SRD*DATAN(GI/GR)
           ENDIF
        ENDIF
        RETURN
        END
