        PROGRAM MSTVLV

C       ======================================================
C       Purpose:  This program computes the modified Struve
C                 function Lv(x) for an arbitrary order v
C                 using subroutine STVLV
C       Input :   v   --- Order of Lv(x)  ( |v| ף 20 )
C                 x   --- Argument of Lv(x) ( x ע 0 )
C       Output:   SLV --- Lv(x)
C       Example:  x = 10.0
C                   v          Lv(x)
C                 ------------------------
C                  0.5     .27785323D+04
C                  1.5     .24996698D+04
C                  2.5     .20254774D+04
C                  3.5     .14816746D+04
C                  4.5     .98173460D+03
C                  5.5     .59154277D+03
C       =====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter v and x '
        READ(*,*)V,X
        WRITE(*,30)V,X
        WRITE(*,*)
        WRITE(*,*)'   v          Lv(x)'
        WRITE(*,*)' ------------------------'
        CAll STVLV(V,X,SLV)
        WRITE(*,20)V,SLV
20      FORMAT(1X,F5.1,D18.8)
30      FORMAT(1X,'v =',F5.1,6X,'x =',F5.1)
        END


        SUBROUTINE STVLV(V,X,SLV)

C       ======================================================
C       Purpose:  Compute modified Struve function Lv(x) with
C                 an arbitrary order v
C       Input :   v   --- Order of Lv(x)  ( |v| ף 20 )
C                 x   --- Argument of Lv(x) ( x ע 0 )
C       Output:   SLV --- Lv(x)
C       Routine called: GAMMA to compute the gamma function
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           IF (V.GT.-1.0.OR.INT(V)-V.EQ.0.5D0) THEN
              SLV=0.0D0
           ELSE IF (V.LT.-1.0D0) THEN
              SLV=(-1)**(INT(0.5D0-V)-1)*1.0D+300
           ELSE IF (V.EQ.-1.0D0) THEN
              SLV=2.0D0/PI
           ENDIF
           RETURN
        ENDIF
        IF (X.LE.40.0D0) THEN
           V0=V+1.5D0
           CALL GAMMA(V0,GA)
           S=2.0D0/(DSQRT(PI)*GA)
           R1=1.0D0
           DO 10 K=1,100
              VA=K+1.5D0
              CALL GAMMA(VA,GA)
              VB=V+K+1.5D0
              CALL GAMMA(VB,GB)
              R1=R1*(0.5D0*X)**2
              R2=R1/(GA*GB)
              S=S+R2
              IF (DABS(R2/S).LT.1.0D-12) GO TO 15
10         CONTINUE
15         SLV=(0.5D0*X)**(V+1.0D0)*S
        ELSE
           SA=-1.0D0/PI*(0.5D0*X)**(V-1.0)
           V0=V+0.5D0
           CALL GAMMA(V0,GA)
           S=-DSQRT(PI)/GA
           R1=-1.0D0
           DO 20 K=1,12
              VA=K+0.5D0
              CALL GAMMA(VA,GA)
              VB=-K+V+0.5D0
              CALL GAMMA(VB,GB)
              R1=-R1/(0.5D0*X)**2
              S=S+R1*GA/GB
20         CONTINUE
           S0=SA*S
           U=DABS(V)
           N=INT(U)
           U0=U-N
           DO 35 L=0,1
              VT=U0+L
              R=1.0D0
              BIV=1.0D0
              DO 25 K=1,16
                 R=-0.125*R*(4.0*VT*VT-(2.0*K-1.0D0)**2)/(K*X)
                 BIV=BIV+R
                 IF (DABS(R/BIV).LT.1.0D-12) GO TO 30
25            CONTINUE
30            IF (L.EQ.0) BIV0=BIV
35         CONTINUE
           BF0=BIV0
           BF1=BIV
           DO 40 K=2,N
              BF=-2.0D0*(K-1.0+U0)/X*BF1+BF0
              BF0=BF1
40            BF1=BF
           IF (N.EQ.0) BIV=BIV0
           IF (N.GT.1) BIV=BF
           SLV=DEXP(X)/DSQRT(2.0D0*PI*X)*BIV+S0
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
