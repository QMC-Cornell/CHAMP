        PROGRAM MPBDV

C       =========================================================
C       Purpose: This program computes the parabolic cylinder
C                functions Dv(x) and their derivatives using
C                subroutine PBDV
C       Input:   x --- Argument of Dv(x)
C                v --- Order of Dv(x)
C       Output:  DV(na) --- Dn+v0(x)
C                DP(na) --- Dn+v0'(x)
C                ( na = |n|, n = int(v), v0 = v-n, |v0| < 1
C                  n = 0,ס1,ס2,תתת, |n| ף 100 )
C                PDF --- Dv(x)
C                PDD --- Dv'(x)
C       Example: v = 5.5,  x =10.0,  v0 = 0.5,  n = 0,1,...,5

C                  n+v0      Dv(x)           Dv'(x)
C                ---------------------------------------
C                  0.5   .43971930D-10  -.21767183D-09
C                  1.5   .43753148D-09  -.21216995D-08
C                  2.5   .43093569D-08  -.20452956D-07
C                  3.5   .41999741D-07  -.19491595D-06
C                  4.5   .40491466D-06  -.18355745D-05
C                  5.5   .38601477D-05  -.17073708D-04

C                Dv(x)= .38601477D-05,  Dv'(x)=-.17073708D-04
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DV(0:100),DP(0:100)
        WRITE(*,*)'Please enter v and  x '
        READ(*,*)V,X
        WRITE(*,20)V,X
        NV=INT(V)
        V0=V-NV
        NA=ABS(NV)
        CALL PBDV(V,X,DV,DP,PDF,PDD)
        WRITE(*,*)
        WRITE(*,*)'   v       Dv(x)           Dv''(x)'
        WRITE(*,*)'---------------------------------------'
        DO 10 K=0,NA
           VK=K*ISIGN(1,NV)+V0
10         WRITE(*,30)VK,DV(K),DP(K)
        WRITE(*,*)
        WRITE(*,40)V,PDF,PDD
20      FORMAT(1X,'v =',F6.2,',    ','x =',F6.2)
30      FORMAT(1X,F5.1,2D16.8)
40      FORMAT(1X,'v =',F5.1,',  Dv(x)=',D14.8,',   Dv''(x)=',D14.8)
        END


        SUBROUTINE PBDV(V,X,DV,DP,PDF,PDD)

C       ====================================================
C       Purpose: Compute parabolic cylinder functions Dv(x)
C                and their derivatives
C       Input:   x --- Argument of Dv(x)
C                v --- Order of Dv(x)
C       Output:  DV(na) --- Dn+v0(x)
C                DP(na) --- Dn+v0'(x)
C                ( na = |n|, v0 = v-n, |v0| < 1,
C                  n = 0,ס1,ס2,תתת )
C                PDF --- Dv(x)
C                PDD --- Dv'(x)
C       Routines called:
C             (1) DVSA for computing Dv(x) for small |x|
C             (2) DVLA for computing Dv(x) for large |x|
C       ====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DV(0:*),DP(0:*)
        XA=DABS(X)
        VH=V
        V=V+DSIGN(1.0D0,V)
        NV=INT(V)
        V0=V-NV
        NA=ABS(NV)
        EP=DEXP(-.25D0*X*X)
        IF (NA.GE.1) JA=1
        IF (V.GE.0.0) THEN
           IF (V0.EQ.0.0) THEN
              PD0=EP
              PD1=X*EP
           ELSE
              DO 10 L=0,JA
                 V1=V0+L
                 IF (XA.LE.5.8) CALL DVSA(V1,X,PD1)
                 IF (XA.GT.5.8) CALL DVLA(V1,X,PD1)
                 IF (L.EQ.0) PD0=PD1
10            CONTINUE
           ENDIF
           DV(0)=PD0
           DV(1)=PD1
           DO 15 K=2,NA
              PDF=X*PD1-(K+V0-1.0D0)*PD0
              DV(K)=PDF
              PD0=PD1
15            PD1=PDF
        ELSE
           IF (X.LE.0.0) THEN
              IF (XA.LE.5.8D0)  THEN
                 CALL DVSA(V0,X,PD0)
                 V1=V0-1.0D0
                 CALL DVSA(V1,X,PD1)
              ELSE
                 CALL DVLA(V0,X,PD0)
                 V1=V0-1.0D0
                 CALL DVLA(V1,X,PD1)
              ENDIF
              DV(0)=PD0
              DV(1)=PD1
              DO 20 K=2,NA
                 PD=(-X*PD1+PD0)/(K-1.0D0-V0)
                 DV(K)=PD
                 PD0=PD1
20               PD1=PD
           ELSE IF (X.LE.2.0) THEN
              V2=NV+V0
              IF (NV.EQ.0) V2=V2-1.0D0
              NK=INT(-V2)
              CALL DVSA(V2,X,F1)
              V1=V2+1.0D0
              CALL DVSA(V1,X,F0)
              DV(NK)=F1
              DV(NK-1)=F0
              DO 25 K=NK-2,0,-1
                 F=X*F0+(K-V0+1.0D0)*F1
                 DV(K)=F
                 F1=F0
25               F0=F
           ELSE
              IF (XA.LE.5.8) CALL DVSA(V0,X,PD0)
              IF (XA.GT.5.8) CALL DVLA(V0,X,PD0)
              DV(0)=PD0
              M=100+NA
              F1=0.0D0
              F0=1.0D-30
              DO 30 K=M,0,-1
                 F=X*F0+(K-V0+1.0D0)*F1
                 IF (K.LE.NA) DV(K)=F
                 F1=F0
30               F0=F
              S0=PD0/F
              DO 35 K=0,NA
35               DV(K)=S0*DV(K)
           ENDIF
        ENDIF
        DO 40 K=0,NA-1
           V1=ABS(V0)+K
           IF (V.GE.0.0D0) THEN
              DP(K)=0.5D0*X*DV(K)-DV(K+1)
           ELSE
              DP(K)=-0.5D0*X*DV(K)-V1*DV(K+1)
           ENDIF
40      CONTINUE
        PDF=DV(NA-1)
        PDD=DP(NA-1)
        V=VH
        RETURN
        END


        SUBROUTINE DVSA(VA,X,PD)

C       ===================================================
C       Purpose: Compute parabolic cylinder function Dv(x)
C                for small argument
C       Input:   x  --- Argument
C                va --- Order
C       Output:  PD --- Dv(x)
C       Routine called: GAMMA for computing ג(x)
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        SQ2=DSQRT(2.0D0)
        EP=DEXP(-.25D0*X*X)
        VA0=0.5D0*(1.0D0-VA)
        IF (VA.EQ.0.0) THEN
           PD=EP
        ELSE
           IF (X.EQ.0.0) THEN
              IF (VA0.LE.0.0.AND.VA0.EQ.INT(VA0)) THEN
                 PD=0.0D0
              ELSE
                 CALL GAMMA(VA0,GA0)
                 PD=DSQRT(PI)/(2.0D0**(-.5D0*VA)*GA0)
              ENDIF
           ELSE
              CALL GAMMA(-VA,G1)
              A0=2.0D0**(-0.5D0*VA-1.0D0)*EP/G1
              VT=-.5D0*VA
              CALL GAMMA(VT,G0)
              PD=G0
              R=1.0D0
              DO 10 M=1,250
                 VM=.5D0*(M-VA)
                 CALL GAMMA(VM,GM)
                 R=-R*SQ2*X/M
                 R1=GM*R
                 PD=PD+R1
                 IF (DABS(R1).LT.DABS(PD)*EPS) GO TO 15
10            CONTINUE
15            PD=A0*PD
           ENDIF
        ENDIF
        RETURN
        END


        SUBROUTINE DVLA(VA,X,PD)

C       ====================================================
C       Purpose: Compute parabolic cylinder functions Dv(x)
C                for large argument
C       Input:   x  --- Argument
C                va --- Order
C       Output:  PD --- Dv(x)
C       Routines called:
C             (1) VVLA for computing Vv(x) for large |x|
C             (2) GAMMA for computing ג(x)
C       ====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EPS=1.0D-12
        EP=DEXP(-.25*X*X)
        A0=DABS(X)**VA*EP
        R=1.0D0
        PD=1.0D0
        DO 10 K=1,16
           R=-0.5D0*R*(2.0*K-VA-1.0)*(2.0*K-VA-2.0)/(K*X*X)
           PD=PD+R
           IF (DABS(R/PD).LT.EPS) GO TO 15
10      CONTINUE
15      PD=A0*PD
        IF (X.LT.0.0D0) THEN
            X1=-X
            CALL VVLA(VA,X1,VL)
            CALL GAMMA(-VA,GL)
            PD=PI*VL/GL+DCOS(PI*VA)*PD
        ENDIF
        RETURN
        END


        SUBROUTINE VVLA(VA,X,PV)

C       ===================================================
C       Purpose: Compute parabolic cylinder function Vv(x)
C                for large argument
C       Input:   x  --- Argument
C                va --- Order
C       Output:  PV --- Vv(x)
C       Routines called:
C             (1) DVLA for computing Dv(x) for large |x|
C             (2) GAMMA for computing ג(x)
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EPS=1.0D-12
        QE=DEXP(0.25*X*X)
        A0=DABS(X)**(-VA-1.0D0)*DSQRT(2.0D0/PI)*QE
        R=1.0D0
        PV=1.0D0
        DO 10 K=1,18
           R=0.5D0*R*(2.0*K+VA-1.0)*(2.0*K+VA)/(K*X*X)
           PV=PV+R
           IF (DABS(R/PV).LT.EPS) GO TO 15
10      CONTINUE
15      PV=A0*PV
        IF (X.LT.0.0D0) THEN
           X1=-X
           CALL DVLA(VA,X1,PDL)
           CALL GAMMA(-VA,GL)
           DSL=DSIN(PI*VA)*DSIN(PI*VA)
           PV=DSL*GL/PI*PDL-DCOS(PI*VA)*PV
        ENDIF
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
