        PROGRAM MLPMV

C       =========================================================
C       Purpose: This program computes the associated Legendre
C                function Pmv(x) with an integer order and an
C                arbitrary nonnegative degree using subroutine
C                LPMV
C       Input :  x   --- Argument of Pm(x)  ( -1 ó x ó 1 )
C                m   --- Order of Pmv(x)
C                v   --- Degree of Pmv(x)
C       Output:  PMV --- Pmv(x)
C       Example:    m = 4,  x = 0.5
C                    v          Pmv(x)
C                 -----------------------
C                   1.5       .46218726
C                   1.6       .48103143
C                   1.7       .45031429
C                   1.8       .36216902
C                   1.9       .21206446
C                   2.0       .00000000
C                   2.5     -1.51996235
C       =========================================================

        IMPLICIT DOUBLE PRECISION (P,V,X)
        WRITE(*,*)'Please enter m,v,x = ?'
        READ(*,*) M,V,X
        WRITE(*,20)M,X
        WRITE(*,*)
        WRITE(*,*)'     v        Pmv(x)'
        WRITE(*,*)'  -----------------------'
        CALL LPMV(V,M,X,PMV)
        WRITE(*,10)V,PMV
10      FORMAT(3X,F5.1,E16.8)
20      FORMAT(3X,'m =',I2,',    ','x =',F6.2)
        END


        SUBROUTINE LPMV(V,M,X,PMV)

C       =======================================================
C       Purpose: Compute the associated Legendre function
C                Pmv(x) with an integer order and an arbitrary
C                nonnegative degree v
C       Input :  x   --- Argument of Pm(x)  ( -1 ó x ó 1 )
C                m   --- Order of Pmv(x)
C                v   --- Degree of Pmv(x)
C       Output:  PMV --- Pmv(x)
C       Routine called:  PSI for computing Psi function
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        EPS=1.0D-14
        NV=INT(V)
        V0=V-NV
        IF (X.EQ.-1.0D0.AND.V.NE.NV) THEN
           IF (M.EQ.0) PMV=-1.0D+300
           IF (M.NE.0) PMV=1.0D+300
           RETURN
        ENDIF
        C0=1.0D0
        IF (M.NE.0) THEN
           RG=V*(V+M)
           DO 10 J=1,M-1
10            RG=RG*(V*V-J*J)
           XQ=DSQRT(1.0D0-X*X)
           R0=1.0D0
           DO 15 J=1,M
15            R0=.5D0*R0*XQ/J
           C0=R0*RG
        ENDIF
        IF (V0.EQ.0.0D0) THEN
           PMV=1.0D0
           R=1.0D0
           DO 20 K=1,NV-M
              R=0.5D0*R*(-NV+M+K-1.0D0)*(NV+M+K)/(K*(K+M))
     &          *(1.0D0+X)
20            PMV=PMV+R
           PMV=(-1)**NV*C0*PMV
        ELSE
           IF (X.GE.-0.35D0) THEN
              PMV=1.0D0
              R=1.0D0
              DO 25 K=1,100
                 R=0.5D0*R*(-V+M+K-1.0D0)*(V+M+K)/(K*(M+K))*(1.0D0-X)
                 PMV=PMV+R
                 IF (K.GT.12.AND.DABS(R/PMV).LT.EPS) GO TO 30
25            CONTINUE
30            PMV=(-1)**M*C0*PMV
           ELSE
              VS=DSIN(V*PI)/PI
              PV0=0.0D0
              IF (M.NE.0) THEN
                 QR=DSQRT((1.0D0-X)/(1.0D0+X))
                 R2=1.0D0
                 DO 35 J=1,M
35                  R2=R2*QR*J
                 S0=1.0D0
                 R1=1.0D0
                 DO 40 K=1,M-1
                    R1=0.5D0*R1*(-V+K-1)*(V+K)/(K*(K-M))*(1.0D0+X)
40                  S0=S0+R1
                 PV0=-VS*R2/M*S0
              ENDIF
              CALL PSI(V,PSV)
              PA=2.0D0*(PSV+EL)+PI/DTAN(PI*V)+1.0D0/V
              S1=0.0D0
              DO 45 J=1,M
45               S1=S1+(J*J+V*V)/(J*(J*J-V*V))
              PMV=PA+S1-1.0D0/(M-V)+DLOG(0.5D0*(1.0D0+X))
              R=1.0D0
              DO 60 K=1,100
                 R=0.5D0*R*(-V+M+K-1.0D0)*(V+M+K)/(K*(K+M))*(1.0D0+X)
                 S=0.0D0
                 DO 50 J=1,M
50                  S=S+((K+J)**2+V*V)/((K+J)*((K+J)**2-V*V))
                 S2=0.0D0
                 DO 55 J=1,K
55                  S2=S2+1.0D0/(J*(J*J-V*V))
                 PSS=PA+S+2.0D0*V*V*S2-1.0D0/(M+K-V)
     &               +DLOG(0.5D0*(1.0D0+X))
                 R2=PSS*R
                 PMV=PMV+R2
                 IF (DABS(R2/PMV).LT.EPS) GO TO 65
60            CONTINUE
65            PMV=PV0+PMV*VS*C0
           ENDIF
        ENDIF
        RETURN
        END


        SUBROUTINE PSI(X,PS)

C       ======================================
C       Purpose: Compute psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+0.5.EQ.INT(XA+0.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
