        PROGRAM MINCOB

C       =========================================================
C       Purpose: This program computes the incomplete beta
C                function Ix(a,b) using subroutine INCOB
C       Input :  a --- Parameter
C                b --- Parameter
C                x --- Argument ( 0 ף x ף 1 )
C       Output:  BIX --- Ix(a,b)
C       Example:
C                  a       b       x       Ix(a,b)
C                -----------------------------------
C                 1.0     3.0     .25     .57812500
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter a, b and x ( 0 ף x ף 1 )'
        READ(*,*)A,B,X
        WRITE(*,*)
        WRITE(*,*)'   a       b       x       Ix(a,b)'
        WRITE(*,*)' -----------------------------------'
        CALL INCOB(A,B,X,BIX)
        WRITE(*,10)A,B,X,BIX
10      FORMAT(1X,F5.1,3X,F5.1,3X,F5.2,F14.8)
        END


        SUBROUTINE INCOB(A,B,X,BIX)

C       ========================================================
C       Purpose: Compute the incomplete beta function Ix(a,b)
C       Input :  a --- Parameter
C                b --- Parameter
C                x --- Argument ( 0 ף x ף 1 )
C       Output:  BIX --- Ix(a,b)
C       Routine called: BETA for computing beta function B(p,q)
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DK(51),FK(51)
        S0=(A+1.0D0)/(A+B+2.0D0)
        CALL BETA(A,B,BT)
        IF (X.LE.S0) THEN
           DO 10 K=1,20
10            DK(2*K)=K*(B-K)*X/(A+2.0D0*K-1.0D0)/(A+2.0D0*K)
           DO 15 K=0,20
15            DK(2*K+1)=-(A+K)*(A+B+K)*X/(A+2.D0*K)/(A+2.0*K+1.0)
           T1=0.0D0
           DO 20 K=20,1,-1
20            T1=DK(K)/(1.0D0+T1)
           TA=1.0D0/(1.0D0+T1)
           BIX=X**A*(1.0D0-X)**B/(A*BT)*TA
        ELSE
           DO 25 K=1,20
25            FK(2*K)=K*(A-K)*(1.0D0-X)/(B+2.*K-1.0)/(B+2.0*K)
           DO 30 K=0,20
30            FK(2*K+1)=-(B+K)*(A+B+K)*(1.D0-X)/
     &                   (B+2.D0*K)/(B+2.D0*K+1.D0)
           T2=0.0D0
           DO 35 K=20,1,-1
35            T2=FK(K)/(1.0D0+T2)
           TB=1.0D0/(1.0D0+T2)
           BIX=1.0D0-X**A*(1.0D0-X)**B/(B*BT)*TB
        ENDIF
        RETURN
        END


        SUBROUTINE BETA(P,Q,BT)

C       ==========================================
C       Purpose: Compute the beta function B(p,q)
C       Input :  p --- Parameter  ( p > 0 )
C                q --- Parameter  ( q > 0 )
C       Output:  BT --- B(p,q)
C       Routine called: GAMMA for computing ג(x)
C       ==========================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CALL GAMMA(P,GP)
        CALL GAMMA(Q,GQ)
        PPQ=P+Q
        CALL GAMMA(PPQ,GPQ)
        BT=GP*GQ/GPQ
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
