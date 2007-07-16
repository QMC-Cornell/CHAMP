        PROGRAM MINCOG

C       ==========================================================
C       Purpose: This program computes the incomplete gamma
C                function r(a,x), ג(a,x) and P(a,x) using
C                subroutine INCOG
C       Input :  a   --- Parameter
C                x   --- Argument
C       Output:  GIN --- r(a,x)
C                GIM --- ג(a,x)
C                GIP --- P(a,x)
C       Example:
C            a     x      r(a,x)         ג(a,x)         P(a,x)
C           -------------------------------------------------------
C           3.0   5.0  .17506960D+01  .24930404D+00  .87534798D+00
C       ===========================================================


        DOUBLE PRECISION A,X,GIN,GIM,GIP
        WRITE(*,*)'Plese enter a and x'
        READ(*,*)A,X
        WRITE(*,*)
        WRITE(*,*)'   a     x      r(a,x)         ג(a,x)         P(a,x)'
        WRITE(*,*)' --------------------------------------------',
     &            '------------'
        CALL INCOG(A,X,GIN,GIM,GIP)
        WRITE(*,10)A,X,GIN,GIM,GIP
10      FORMAT(1X,F5.1,1X,F5.1,3D15.8)
        END


        SUBROUTINE INCOG(A,X,GIN,GIM,GIP)

C       ===================================================
C       Purpose: Compute the incomplete gamma function
C                r(a,x), ג(a,x) and P(a,x)
C       Input :  a   --- Parameter ( a ף 170 )
C                x   --- Argument
C       Output:  GIN --- r(a,x)
C                GIM --- ג(a,x)
C                GIP --- P(a,x)
C       Routine called: GAMMA for computing ג(x)
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'a and/or x too large'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           CALL GAMMA(A,GA)
           GIM=GA
           GIP=0.0
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           CALL GAMMA(A,GA)
           GIP=GIN/GA
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
           CALL GAMMA(A,GA)
           GIN=GA-GIM
           GIP=1.0D0-GIM/GA
        ENDIF
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
