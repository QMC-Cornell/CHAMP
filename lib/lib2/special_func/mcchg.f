        PROGRAM MCCHG

C       ===========================================================
C       Purpose: This program computes confluent hypergeometric
C                function M(a,b,z) with real parameters a, b, and
C                a complex argument z using subroutine CCHG
C       Input :  a --- Parameter
C                b --- Parameter
C                z --- Complex argument
C       Output:  CHG --- M(a,b,z)
C       Examples:
C          a      b        z        Re[M(a,b,z)]   Im[M(a,b,z)]
C         -------------------------------------------------------
C         3.3   4.25    10 + 0i    .61677489D+04    0
C         3.3   4.25    25 + 0i    .95781835D+10  -.15738228D-03
C         3.3   4.25     3 -  i    .75828716D+01  -.86815474D+01
C         3.3   4.25    15 +10i   -.58313765D+06  -.48195426D+05
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A,B,X,Y)
        IMPLICIT COMPLEX *16 (C,Z)
        WRITE(*,*)'Please enter a, b, x and y (z=x+iy) '
        READ(*,*)A,B,X,Y
        WRITE(*,20)A,B,X,Y
        Z=CMPLX(X,Y)
        CALL CCHG(A,B,Z,CHG)
        WRITE(*,10)CHG
10      FORMAT(10X,'M(a,b,z) =',D18.8,' + i ',D18.8)
20      FORMAT(1X,'a =',F5.1,',  ','b =',F5.1,',  ','x =',F5.1,
     &         ',  ','y =',F5.1)
        END


        SUBROUTINE CCHG(A,B,Z,CHG)

C       ===================================================
C       Purpose: Compute confluent hypergeometric function
C                M(a,b,z) with real parameters a, b and a
C                complex argument z
C       Input :  a --- Parameter
C                b --- Parameter
C                z --- Complex argument
C       Output:  CHG --- M(a,b,z)
C       Routine called: GAMMA for computing gamma function
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX *16 (C,Z)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=A
        A1=A
        Z0=Z
        IF (B.EQ.0.0.OR.B.EQ.-INT(ABS(B))) THEN
           CHG=(1.0D+300,0.0D0)
        ELSE IF (A.EQ.0.0D0.OR.Z.EQ.0.0D0) THEN
           CHG=(1.0D0,0.0D0)
        ELSE IF (A.EQ.-1.0D0) THEN
           CHG=1.0D0-Z/B
        ELSE IF (A.EQ.B) THEN
           CHG=CDEXP(Z)
        ELSE IF (A-B.EQ.1.0D0) THEN
           CHG=(1.0D0+Z/B)*CDEXP(Z)
        ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
           CHG=(CDEXP(Z)-1.0D0)/Z
        ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
           M=INT(-A)
           CR=(1.0D0,0.0D0)
           CHG=(1.0D0,0.0D0)
           DO 10 K=1,M
              CR=CR*(A+K-1.0D0)/K/(B+K-1.0D0)*Z
10            CHG=CHG+CR
        ELSE
           X0=REAL(Z)
           IF (X0.LT.0.0D0) THEN
              A=B-A
              A0=A
              Z=-Z
           ENDIF
           IF (A.LT.2.0D0) NL=0
           IF (A.GE.2.0D0) THEN
              NL=1
              LA=INT(A)
              A=A-LA-1.0D0
           ENDIF
           DO 30 N=0,NL
              IF (A0.GE.2.0D0) A=A+1.0D0
              IF (CDABS(Z).LT.20.0D0+ABS(B).OR.A.LT.0.0D0) THEN
                 CHG=(1.0D0,0.0D0)
                 CRG=(1.0D0,0.0D0)
                 DO 15 J=1,500
                    CRG=CRG*(A+J-1.0D0)/(J*(B+J-1.0D0))*Z
                    CHG=CHG+CRG
                    IF (CDABS((CHG-CHW)/CHG).LT.1.D-15) GO TO 25
                    CHW=CHG
15               CONTINUE
              ELSE
                 CALL GAMMA(A,G1)
                 CALL GAMMA(B,G2)
                 BA=B-A
                 CALL GAMMA(BA,G3)
                 CS1=(1.0D0,0.0D0)
                 CS2=(1.0D0,0.0D0)
                 CR1=(1.0D0,0.0D0)
                 CR2=(1.0D0,0.0D0)
                 DO 20 I=1,8
                    CR1=-CR1*(A+I-1.0D0)*(A-B+I)/(Z*I)
                    CR2=CR2*(B-A+I-1.0D0)*(I-A)/(Z*I)
                    CS1=CS1+CR1
20                  CS2=CS2+CR2
                 X=REAL(Z)
                 Y=DIMAG(Z)
                 IF (X.EQ.0.0.AND.Y.GE.0.0) THEN
                    PHI=0.5D0*PI
                 ELSE IF (X.EQ.0.0.AND.Y.LE.0.0) THEN
                    PHI=-0.5D0*PI
                 ELSE
                    PHI=DATAN(Y/X)
                 ENDIF
                 IF (PHI.GT.-0.5*PI.AND.PHI.LT.1.5*PI) NS=1
                 IF (PHI.GT.-1.5*PI.AND.PHI.LE.-0.5*PI) NS=-1
                 CFAC=CDEXP(NS*CI*PI*A)
                 IF (Y.EQ.0.0D0) CFAC=DCOS(PI*A)
                 CHG1=G2/G3*Z**(-A)*CFAC*CS1
                 CHG2=G2/G1*CDEXP(Z)*Z**(A-B)*CS2
                 CHG=CHG1+CHG2
              ENDIF
25            IF (N.EQ.0) CY0=CHG
              IF (N.EQ.1) CY1=CHG
30         CONTINUE
           IF (A0.GE.2.0D0) THEN
              DO 35 I=1,LA-1
                 CHG=((2.0D0*A-B+Z)*CY1+(B-A)*CY0)/A
                 CY0=CY1
                 CY1=CHG
35               A=A+1.0D0
           ENDIF
           IF (X0.LT.0.0D0) CHG=CHG*CDEXP(-Z)
        ENDIF
        A=A1
        Z=Z0
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
