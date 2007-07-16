        PROGRAM MSTVHV

C       ======================================================
C       Purpose:  This program computes Struve function Hv(x)
C                 for an arbitrary order using subroutine
C                 STVHV
C       Input :   v  --- Order of Hv(x)  ( -8.0 ף v ף 12.5 )
C                 x  --- Argument of Hv(x) ( x ע 0 )
C       Output:   HV --- Hv(x)
C       Example:  x = 10.0
C                   v           Hv(x)
C                 -----------------------
C                   .5     .46402212D+00
C                  1.5     .14452322D+01
C                  2.5     .31234632D+01
C                  3.5     .53730255D+01
C                  4.5     .72083122D+01
C                  5.5     .76851132D+01
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please v and x '
        READ(*,*)V,X
        WRITE(*,30)V,X
        WRITE(*,*)
        WRITE(*,*)'   v           Hv(x)'
        WRITE(*,*)' -----------------------'
        CAll STVHV(V,X,HV)
        WRITE(*,20)V,HV
20      FORMAT(1X,F5.1,D18.8)
30      FORMAT(1X,'v =',F5.1,6X,'x =',F5.1)
        END


        SUBROUTINE STVHV(V,X,HV)

C       =====================================================
C       Purpose: Compute Struve function Hv(x) with an
C                arbitrary order v
C       Input :  v  --- Order of Hv(x)  ( -8.0 ף v ף 12.5 )
C                x  --- Argument of Hv(x) ( x ע 0 )
C       Output:  HV --- Hv(x)
C       Routine called: GAMMA to compute the gamma function
C       =====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           IF (V.GT.-1.0.OR.INT(V)-V.EQ.0.5D0) THEN
              HV=0.0D0
           ELSE IF (V.LT.-1.0D0) THEN
              HV=(-1)**(INT(0.5D0-V)-1)*1.0D+300
           ELSE IF (V.EQ.-1.0D0) THEN
              HV=2.0D0/PI
           ENDIF
           RETURN
        ENDIF
        IF (X.LE.20.0D0) THEN
           V0=V+1.5D0
           CALL GAMMA(V0,GA)
           S=2.0D0/(DSQRT(PI)*GA)
           R1=1.0D0
           DO 10 K=1,100
              VA=K+1.5D0
              CALL GAMMA(VA,GA)
              VB=V+K+1.5D0
              CALL GAMMA(VB,GB)
              R1=-R1*(0.5D0*X)**2
              R2=R1/(GA*GB)
              S=S+R2
              IF (DABS(R2).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         HV=(0.5D0*X)**(V+1.0D0)*S
        ELSE
           SA=(0.5D0*X)**(V-1.0)/PI
           V0=V+0.5D0
           CALL GAMMA(V0,GA)
           S=DSQRT(PI)/GA
           R1=1.0D0
           DO 20 K=1,12
              VA=K+0.5D0
              CALL GAMMA(VA,GA)
              VB=-K+V+0.5D0
              CALL GAMMA(VB,GB)
              R1=R1/(0.5D0*X)**2
              S=S+R1*GA/GB
20         CONTINUE
           S0=SA*S
           U=DABS(V)
           N=INT(U)
           U0=U-N
           DO 35 L=0,1
              VT=4.0D0*(U0+L)**2
              R1=1.0D0
              PU1=1.0D0
              DO 25 K=1,12
                 R1=-0.0078125D0*R1*(VT-(4.0*K-3.0D0)**2)*
     &             (VT-(4.0D0*K-1.0)**2)/((2.0D0*K-1.0)*K*X*X)
                 PU1=PU1+R1
25            CONTINUE
              QU1=1.0D0
              R2=1.0D0
              DO 30 K=1,12
                 R2=-0.0078125D0*R2*(VT-(4.0D0*K-1.0)**2)*
     &             (VT-(4.0D0*K+1.0)**2)/((2.0D0*K+1.0)*K*X*X)
                 QU1=QU1+R2
30            CONTINUE
              QU1=0.125D0*(VT-1.0D0)/X*QU1
              IF (L.EQ.0) THEN
                 PU0=PU1
                 QU0=QU1
              ENDIF
35         CONTINUE
           T0=X-(0.5*U0+0.25D0)*PI
           T1=X-(0.5*U0+0.75D0)*PI
           SR=DSQRT(2.0D0/(PI*X))
           BY0=SR*(PU0*DSIN(T0)+QU0*DCOS(T0))
           BY1=SR*(PU1*DSIN(T1)+QU1*DCOS(T1))
           BF0=BY0
           BF1=BY1
           DO 40 K=2,N
              BF=2.0D0*(K-1.0+U0)/X*BF1-BF0
              BF0=BF1
40            BF1=BF
           IF (N.EQ.0) BYV=BY0
           IF (N.EQ.1) BYV=BY1
           IF (N.GT.1) BYV=BF
           HV=BYV+S0
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
