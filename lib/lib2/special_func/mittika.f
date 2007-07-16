        PROGRAM MITTIKA

C       ============================================================
C       Purpose: This program computes the integral of [I0(t)-1]/t
C                with respect to t from 0 to x and K0(t)/t with
C                respect to t from x to ì using subroutine ITTIKA
C       Input :  x   --- Variable in the limits  ( x ò 0 )
C       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
C                TTK --- Integration of K0(t)/t from x to ì
C       Example:
C                   x     [1-I0(t)]/tdt     K0(t)/tdt
C                ---------------------------------------
C                  5.0   .71047763D+01   .58635626D-03
C                 10.0   .34081537D+03   .15629282D-05
C                 15.0   .25437619D+05   .59837472D-08
C                 20.0   .23673661D+07   .26790545D-10
C                 25.0   .24652751D+09   .13100706D-12
C       ============================================================

        DOUBLE PRECISION X,TTI,TTK
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x     [1-I0(t)]/tdt     K0(t)/tdt'
        WRITE(*,*)'---------------------------------------'
        CALL ITTIKA(X,TTI,TTK)
        WRITE(*,10)X,TTI,TTK
10      FORMAT(1X,F5.1,2D16.8)
        END


        SUBROUTINE ITTIKA(X,TTI,TTK)

C       =========================================================
C       Purpose: Integrate [I0(t)-1]/t with respect to t from 0
C                to x, and K0(t)/t with respect to t from x to ì
C       Input :  x   --- Variable in the limits  ( x ò 0 )
C       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
C                TTK --- Integration of K0(t)/t from x to ì
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION C(8)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        DATA C/1.625D0,4.1328125D0,
     &       1.45380859375D+1,6.553353881835D+1,
     &       3.6066157150269D+2,2.3448727161884D+3,
     &       1.7588273098916D+4,1.4950639538279D+5/
        IF (X.EQ.0.0D0) THEN
           TTI=0.0D0
           TTK=1.0D+300
           RETURN
        ENDIF
        IF (X.LT.40.0D0) THEN
           TTI=1.0D0
           R=1.0D0
           DO 10 K=2,50
              R=.25D0*R*(K-1.0D0)/(K*K*K)*X*X
              TTI=TTI+R
              IF (DABS(R/TTI).LT.1.0D-12) GO TO 15
10         CONTINUE
15         TTI=TTI*.125D0*X*X
        ELSE
           TTI=1.0D0
           R=1.0D0
           DO 20 K=1,8
              R=R/X
20            TTI=TTI+C(K)*R
           RC=X*DSQRT(2.0D0*PI*X)
           TTI=TTI*DEXP(X)/RC
        ENDIF
        IF (X.LE.12.0D0) THEN
           E0=(.5D0*DLOG(X/2.0D0)+EL)*DLOG(X/2.0D0)
     &        +PI*PI/24.0D0+.5D0*EL*EL
           B1=1.5D0-(EL+DLOG(X/2.0D0))
           RS=1.0D0
           R=1.0D0
           DO 25 K=2,50
              R=.25D0*R*(K-1.0D0)/(K*K*K)*X*X
              RS=RS+1.0D0/K
              R2=R*(RS+1.0D0/(2.0D0*K)-(EL+DLOG(X/2.0D0)))
              B1=B1+R2
              IF (DABS(R2/B1).LT.1.0D-12) GO TO 30
25         CONTINUE
30         TTK=E0-.125D0*X*X*B1
        ELSE
           TTK=1.0D0
           R=1.0D0
           DO 35 K=1,8
              R=-R/X
35            TTK=TTK+C(K)*R
           RC=X*DSQRT(2.0D0/PI*X)
           TTK=TTK*DEXP(-X)/RC
        ENDIF
        RETURN
        END
