        PROGRAM MITJYA

C       ===========================================================
C       Purpose: This program evaluates the integral of Bessel
C                functions J0(t) and Y0(t) with respect to t
C                from 0 to x using subroutine ITJYA
C       Input :  x  --- Upper limit of the integral ( x ò 0 )
C       Output:  TJ --- Integration of J0(t) from 0 to x
C                TY --- Integration of Y0(t) from 0 to x
C       Example:
C                   x         J0(t)dt          Y0(t)dt
C                ---------------------------------------
C                  5.0       .71531192       .19971938
C                 10.0      1.06701130       .24129032
C                 15.0      1.20516194       .00745772
C                 20.0      1.05837882      -.16821598
C                 25.0       .87101492      -.09360793
C                 30.0       .88424909       .08822971
C       ===========================================================

        DOUBLE PRECISION X,TJ,TY
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x         J0(t)dt          Y0(t)dt'
        WRITE(*,*)'---------------------------------------'
        CALL ITJYA(X,TJ,TY)
        WRITE(*,10)X,TJ,TY
10      FORMAT(1X,F5.1,2F16.8)
        END


        SUBROUTINE ITJYA(X,TJ,TY)

C       ==========================================================
C       Purpose: Integrate Bessel functions J0(t) & Y0(t) with
C                respect to t from 0 to x
C       Input :  x  --- Upper limit of the integral ( x ò 0 )
C       Output:  TJ --- Integration of J0(t) from 0 to x
C                TY --- Integration of Y0(t) from 0 to x
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(18)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        EPS=1.0D-12
        IF (X.EQ.0.0D0) THEN
           TJ=0.0D0
           TY=0.0D0
        ELSE IF (X.LE.20.0D0) THEN
           X2=X*X
           TJ=X
           R=X
           DO 10 K=1,60
              R=-.25D0*R*(2*K-1.0D0)/(2*K+1.0D0)/(K*K)*X2
              TJ=TJ+R
              IF (DABS(R).LT.DABS(TJ)*EPS) GO TO 15
10         CONTINUE
15         TY1=(EL+DLOG(X/2.0D0))*TJ
           RS=0.0D0
           TY2=1.0D0
           R=1.0D0
           DO 20 K=1,60
              R=-.25D0*R*(2*K-1.0D0)/(2*K+1.0D0)/(K*K)*X2
              RS=RS+1.0D0/K
              R2=R*(RS+1.0D0/(2.0D0*K+1.0D0))
              TY2=TY2+R2
              IF (DABS(R2).LT.DABS(TY2)*EPS) GO TO 25
20         CONTINUE
25         TY=(TY1-X*TY2)*2.0D0/PI
        ELSE
           A0=1.0D0
           A1=5.0D0/8.0D0
           A(1)=A1
           DO 30 K=1,16
              AF=((1.5D0*(K+.5D0)*(K+5.0D0/6.0D0)*A1-.5D0
     &           *(K+.5D0)*(K+.5D0)*(K-.5D0)*A0))/(K+1.0D0)
              A(K+1)=AF
              A0=A1
30            A1=AF
           BF=1.0D0
           R=1.0D0
           DO 35 K=1,8
              R=-R/(X*X)
35            BF=BF+A(2*K)*R
           BG=A(1)/X
           R=1.0D0/X
           DO 40 K=1,8
              R=-R/(X*X)
40            BG=BG+A(2*K+1)*R
           XP=X+.25D0*PI
           RC=DSQRT(2.0D0/(PI*X))
           TJ=1.0D0-RC*(BF*DCOS(XP)+BG*DSIN(XP))
           TY=RC*(BG*DCOS(XP)-BF*DSIN(XP))
        ENDIF
        RETURN
        END
