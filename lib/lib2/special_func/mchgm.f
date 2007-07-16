        PROGRAM MCHGM

C       =======================================================
C       Purpose: This program computes the confluent
C                hypergeometric function M(a,b,x) using
C                subroutine CHGM
C       Input  : a  --- Parameter
C                b  --- Parameter ( b <> 0,-1,-2,... )
C                x  --- Argument
C       Output:  HG --- M(a,b,x)
C       Example:
C                   a       b       x          M(a,b,x)
C                 -----------------------------------------
C                  1.5     2.0    20.0     .1208527185D+09
C                  4.5     2.0    20.0     .1103561117D+12
C                 -1.5     2.0    20.0     .1004836854D+05
C                 -4.5     2.0    20.0    -.3936045244D+03
C                  1.5     2.0    50.0     .8231906643D+21
C                  4.5     2.0    50.0     .9310512715D+25
C                 -1.5     2.0    50.0     .2998660728D+16
C                 -4.5     2.0    50.0    -.1806547113D+13
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter a, b and x '
        READ(*,*)A,B,X
        WRITE(*,*)'   a       b       x          M(a,b,x)'
        WRITE(*,*)' -----------------------------------------'
        CALL CHGM(A,B,X,HG)
        WRITE(*,10)A,B,X,HG
10      FORMAT(1X,F5.1,3X,F5.1,3X,F5.1,D20.10)
        END


        SUBROUTINE CHGM(A,B,X,HG)

C       ===================================================
C       Purpose: Compute confluent hypergeometric function
C                M(a,b,x)
C       Input  : a  --- Parameter
C                b  --- Parameter ( b <> 0,-1,-2,... )
C                x  --- Argument
C       Output:  HG --- M(a,b,x)
C       Routine called: GAMMA for computing ג(x)
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        A0=A
        A1=A
        X0=X
        HG=0.0D0
        IF (B.EQ.0.0D0.OR.B.EQ.-ABS(INT(B))) THEN
           HG=1.0D+300
        ELSE IF (A.EQ.0.0D0.OR.X.EQ.0.0D0) THEN
           HG=1.0D0
        ELSE IF (A.EQ.-1.0D0) THEN
           HG=1.0D0-X/B
        ELSE IF (A.EQ.B) THEN
           HG=DEXP(X)
        ELSE IF (A-B.EQ.1.0D0) THEN
           HG=(1.0D0+X/B)*DEXP(X)
        ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
           HG=(DEXP(X)-1.0D0)/X
        ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
           M=INT(-A)
           R=1.0D0
           HG=1.0D0
           DO 10 K=1,M
              R=R*(A+K-1.0D0)/K/(B+K-1.0D0)*X
10            HG=HG+R
        ENDIF
        IF (HG.NE.0.0D0) RETURN
        IF (X.LT.0.0D0) THEN
           A=B-A
           A0=A
           X=DABS(X)
        ENDIF
        IF (A.LT.2.0D0) NL=0
        IF (A.GE.2.0D0) THEN
           NL=1
           LA=INT(A)
           A=A-LA-1.0D0
        ENDIF
        DO 30 N=0,NL
           IF (A0.GE.2.0D0) A=A+1.0D0
           IF (X.LE.30.0D0+DABS(B).OR.A.LT.0.0D0) THEN
              HG=1.0D0
              RG=1.0D0
              DO 15 J=1,500
                 RG=RG*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
                 HG=HG+RG
                 IF (DABS(RG/HG).LT.1.0D-15) GO TO 25
15            CONTINUE
           ELSE
              CALL GAMMA(A,TA)
              CALL GAMMA(B,TB)
              XG=B-A
              CALL GAMMA(XG,TBA)
              SUM1=1.0D0
              SUM2=1.0D0
              R1=1.0D0
              R2=1.0D0
              DO 20 I=1,8
                 R1=-R1*(A+I-1.0D0)*(A-B+I)/(X*I)
                 R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)
                 SUM1=SUM1+R1
20               SUM2=SUM2+R2
              HG1=TB/TBA*X**(-A)*DCOS(PI*A)*SUM1
              HG2=TB/TA*DEXP(X)*X**(A-B)*SUM2
              HG=HG1+HG2
           ENDIF
25         IF (N.EQ.0) Y0=HG
           IF (N.EQ.1) Y1=HG
30      CONTINUE
        IF (A0.GE.2.0D0) THEN
           DO 35 I=1,LA-1
              HG=((2.0D0*A-B+X)*Y1+(B-A)*Y0)/A
              Y0=Y1
              Y1=HG
35            A=A+1.0D0
        ENDIF
        IF (X0.LT.0.0D0) HG=HG*DEXP(X0)
        A=A1
        X=X0
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
