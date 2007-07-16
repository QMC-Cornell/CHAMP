        PROGRAM MPBVV

C       ========================================================
C       Purpose: This program computes the parabolic cylinder
C                functions Vv(x) and Vv'(x) using subroutine
C                PBVV
C       Input:   x --- Argument of Vv(x)
C                v --- Order of Vv(x)
C       Output:  VV(na) --- Vv(x)
C                VP(na) --- Vv'(x)
C                ( na = |n|, v = n+v0, n = int(v), |v0| < 1
C                  n = 0,ס1,ס2,תתת, |n| ף 100 )
C                PVF --- Vv(x)
C                PVD --- Vv'(x)
C       Example: v = 5.5,  x =10.0,  v0 = 0.5,  n = 0,1,2,...,5

C                  n+v0      Vv(x)           Vv'(x)
C                ---------------------------------------
C                  0.5   .18522719D+10   .89761157D+10
C                  1.5   .19016268D+09   .90145854D+09
C                  2.5   .19741946D+08   .91452949D+08
C                  3.5   .20733667D+07   .93751130D+07
C                  4.5   .22038231D+06   .97145511D+06
C                  5.5   .23719356D+05   .10178553D+06

C                Vv(x)= .23719356D+05,  Vv'(x)= .10178553D+06
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION VV(0:100),VP(0:100)
        WRITE(*,*)'Please enter v and  x '
        READ(*,*)V,X
        WRITE(*,20)V,X
        NV=INT(V)
        V0=V-NV
        NA=ABS(NV)
        CALL PBVV(V,X,VV,VP,PVF,PVD)
        WRITE(*,*)
        WRITE(*,*)'   v       Vv(x)           Vv''(x)'
        WRITE(*,*)'---------------------------------------'
        DO 10 K=0,NA
           VK=K*ISIGN(1,NV)+V0
10         WRITE(*,30)VK,VV(K),VP(K)
        WRITE(*,*)
        WRITE(*,40)V,PVF,PVD
20      FORMAT(1X,'v =',F6.2,',    ','x =',F6.2)
30      FORMAT(1X,F5.1,2D16.8)
40      FORMAT(1X,'v =',F5.1,',  Vv(x)=',D14.8,',   Vv''(x)=',D14.8)
        END


        SUBROUTINE PBVV(V,X,VV,VP,PVF,PVD)

C       ===================================================
C       Purpose: Compute parabolic cylinder functions Vv(x)
C                and their derivatives
C       Input:   x --- Argument of Vv(x)
C                v --- Order of Vv(x)
C       Output:  VV(na) --- Vv(x)
C                VP(na) --- Vv'(x)
C                ( na = |n|, v = n+v0, |v0| < 1
C                  n = 0,ס1,ס2,תתת )
C                PVF --- Vv(x)
C                PVD --- Vv'(x)
C       Routines called:
C             (1) VVSA for computing Vv(x) for small |x|
C             (2) VVLA for computing Vv(x) for large |x|
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION VV(0:*),VP(0:*)
        PI=3.141592653589793D0
        XA=DABS(X)
        VH=V
        V=V+DSIGN(1.0D0,V)
        NV=INT(V)
        V0=V-NV
        NA=ABS(NV)
        QE=DEXP(0.25D0*X*X)
        Q2P=DSQRT(2.0D0/PI)
        IF (NA.GE.1) JA=1
        IF (V.LE.0.0) THEN
           IF (V0.EQ.0.0) THEN
              IF (XA.LE.7.5) CALL VVSA(V0,X,PV0)
              IF (XA.GT.7.5) CALL VVLA(V0,X,PV0)
              F0=Q2P*QE
              F1=X*F0
              VV(0)=PV0
              VV(1)=F0
              VV(2)=F1
           ELSE
              DO 10 L=0,JA
                 V1=V0-L
                 IF (XA.LE.7.5) CALL VVSA(V1,X,F1)
                 IF (XA.GT.7.5) CALL VVLA(V1,X,F1)
                 IF (L.EQ.0) F0=F1
10            CONTINUE
              VV(0)=F0
              VV(1)=F1
           ENDIF
           KV=2
           IF (V0.EQ.0.0) KV=3
           DO 15 K=KV,NA
              F=X*F1+(K-V0-2.0D0)*F0
              VV(K)=F
              F0=F1
15            F1=F
        ELSE
           IF (X.GE.0.0.AND.X.LE.7.5D0) THEN
              V2=V
              IF (V2.LT.1.0) V2=V2+1.0D0
              CALL VVSA(V2,X,F1)
              V1=V2-1.0D0
              KV=INT(V2)
              CALL VVSA(V1,X,F0)
              VV(KV)=F1
              VV(KV-1)=F0
              DO 20 K=KV-2,0,-1
                 F=X*F0-(K+V0+2.0D0)*F1
                 IF (K.LE.NA) VV(K)=F
                 F1=F0
20               F0=F
           ELSE IF (X.GT.7.5D0) THEN
              CALL VVLA(V0,X,PV0)
              M=100+ABS(NA)
              VV(1)=PV0
              F1=0.0D0
              F0=1.0D-40
              DO 25 K=M,0,-1
                 F=X*F0-(K+V0+2.0D0)*F1
                 IF (K.LE.NA) VV(K)=F
                 F1=F0
25               F0=F
              S0=PV0/F
              DO 30 K=0,NA
30               VV(K)=S0*VV(K)
           ELSE
              IF (XA.LE.7.5D0) THEN
                 CALL VVSA(V0,X,F0)
                 V1=V0+1.0
                 CALL VVSA(V1,X,F1)
              ELSE
                 CALL VVLA(V0,X,F0)
                 V1=V0+1.0D0
                 CALL VVLA(V1,X,F1)
              ENDIF
              VV(0)=F0
              VV(1)=F1
              DO 35 K=2,NA
                 F=(X*F1-F0)/(K+V0)
                 VV(K)=F
                 F0=F1
35               F1=F
           ENDIF
        ENDIF
        DO 40 K=0,NA-1
           V1=V0+K
           IF (V.GE.0.0D0) THEN
              VP(K)=0.5D0*X*VV(K)-(V1+1.0D0)*VV(K+1)
           ELSE
              VP(K)=-0.5D0*X*VV(K)+VV(K+1)
           ENDIF
40      CONTINUE
        PVF=VV(NA-1)
        PVD=VP(NA-1)
        V=VH
        RETURN
        END


        SUBROUTINE VVSA(VA,X,PV)

C       ===================================================
C       Purpose: Compute parabolic cylinder function Vv(x)
C                for small argument
C       Input:   x  --- Argument
C                va --- Order
C       Output:  PV --- Vv(x)
C       Routine called : GAMMA for computing ג(x)
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        EP=DEXP(-.25D0*X*X)
        VA0=1.0D0+0.5D0*VA
        IF (X.EQ.0.0) THEN
           IF (VA0.LE.0.0.AND.VA0.EQ.INT(VA0).OR.VA.EQ.0.0) THEN
              PV=0.0D0
           ELSE
              VB0=-0.5D0*VA
              SV0=DSIN(VA0*PI)
              CALL GAMMA(VA0,GA0)
              PV=2.0D0**VB0*SV0/GA0
           ENDIF
        ELSE
           SQ2=DSQRT(2.0D0)
           A0=2.0D0**(-.5D0*VA)*EP/(2.0D0*PI)
           SV=DSIN(-(VA+.5D0)*PI)
           V1=-.5D0*VA
           CALL GAMMA(V1,G1)
           PV=(SV+1.0D0)*G1
           R=1.0D0
           FAC=1.0D0
           DO 10 M=1,250
              VM=.5D0*(M-VA)
              CALL GAMMA(VM,GM)
              R=R*SQ2*X/M
              FAC=-FAC
              GW=FAC*SV+1.0D0
              R1=GW*R*GM
              PV=PV+R1
              IF (DABS(R1/PV).LT.EPS.AND.GW.NE.0.0) GO TO 15
10         CONTINUE
15         PV=A0*PV
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
        A0=X**VA*EP
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
