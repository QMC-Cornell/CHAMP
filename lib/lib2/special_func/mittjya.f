        PROGRAM MITTJYA

C       ===========================================================
C       Purpose: This program computes the integral of [1-J0(t)]/t
C                with respect to t from 0 to x and Y0(t)/t with
C                respect to t from x to ì using subroutine ITTJYA
C       Input :  x   --- Variable in the limits  ( x ò 0 )
C       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
C                TTY --- Integration of Y0(t)/t from x to ì
C       Example:
C                  x       [1-J0(t)]/tdt       Y0(t)/tdt
C               -------------------------------------------
C                 5.0     .15403472D+01    -.46322055D-01
C                10.0     .21778664D+01    -.22987934D-01
C                15.0     .25785507D+01     .38573574D-03
C                20.0     .28773106D+01     .85031527D-02
C                25.0     .31082313D+01     .35263393D-02
C       ===========================================================

        DOUBLE PRECISION X,TTJ,TTY
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x      [1-J0(t)]/tdt       Y0(t)/tdt'
        WRITE(*,*)'-------------------------------------------'
        CALL ITTJYA(X,TTJ,TTY)
        WRITE(*,10)X,TTJ,TTY
10      FORMAT(1X,F5.1,2D18.8)
        END


        SUBROUTINE ITTJYA(X,TTJ,TTY)

C       =========================================================
C       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
C                to x, and Y0(t)/t with respect to t from x to ì
C       Input :  x   --- Variable in the limits  ( x ò 0 )
C       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
C                TTY --- Integration of Y0(t)/t from x to ì
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IF (X.EQ.0.0D0) THEN
           TTJ=0.0D0
           TTY=-1.0D+300
        ELSE IF (X.LE.20.0D0) THEN
           TTJ=1.0D0
           R=1.0D0
           DO 10 K=2,100
              R=-.25D0*R*(K-1.0D0)/(K*K*K)*X*X
              TTJ=TTJ+R
              IF (DABS(R).LT.DABS(TTJ)*1.0D-12) GO TO 15
10         CONTINUE
15         TTJ=TTJ*.125D0*X*X
           E0=.5D0*(PI*PI/6.0D0-EL*EL)-(.5D0*DLOG(X/2.0D0)+EL)
     &        *DLOG(X/2.0D0)
           B1=EL+DLOG(X/2.0D0)-1.5D0
           RS=1.0D0
           R=-1.0D0
           DO 20 K=2,100
              R=-.25D0*R*(K-1.0D0)/(K*K*K)*X*X
              RS=RS+1.0D0/K
              R2=R*(RS+1.0D0/(2.0D0*K)-(EL+DLOG(X/2.0D0)))
              B1=B1+R2
              IF (DABS(R2).LT.DABS(B1)*1.0D-12) GO TO 25
20         CONTINUE
25         TTY=2.0D0/PI*(E0+.125D0*X*X*B1)
        ELSE
           A0=DSQRT(2.0D0/(PI*X))
           DO 50 L=0,1
              VT=4.0D0*L*L
              PX=1.0D0
              R=1.0D0
              DO 30 K=1,14
                 R=-.0078125D0*R*(VT-(4.0D0*K-3.0D0)**2)
     &             /(X*K)*(VT-(4.0D0*K-1.0D0)**2)
     &             /((2.0D0*K-1.0D0)*X)
                 PX=PX+R
                 IF (DABS(R).LT.DABS(PX)*1.0D-12) GO TO 35
30            CONTINUE
35            QX=1.0D0
              R=1.0D0
              DO 40 K=1,14
                 R=-.0078125D0*R*(VT-(4.0D0*K-1.0D0)**2)
     &             /(X*K)*(VT-(4.0D0*K+1.0D0)**2)
     &             /(2.0D0*K+1.0D0)/X
                 QX=QX+R
                 IF (DABS(R).LT.DABS(QX)*1.0D-12) GO TO 45
40            CONTINUE
45            QX=.125D0*(VT-1.0D0)/X*QX
              XK=X-(.25D0+.5D0*L)*PI
              BJ1=A0*(PX*DCOS(XK)-QX*DSIN(XK))
              BY1=A0*(PX*DSIN(XK)+QX*DCOS(XK))
              IF (L.EQ.0) THEN
                 BJ0=BJ1
                 BY0=BY1
              ENDIF
50         CONTINUE
           T=2.0D0/X
           G0=1.0D0
           R0=1.0D0
           DO 55 K=1,10
              R0=-K*K*T*T*R0
55            G0=G0+R0
           G1=1.0D0
           R1=1.0D0
           DO 60 K=1,10
              R1=-K*(K+1.0D0)*T*T*R1
60            G1=G1+R1
           TTJ=2.0D0*G1*BJ0/(X*X)-G0*BJ1/X+EL+DLOG(X/2.0D0)
           TTY=2.0D0*G1*BY0/(X*X)-G0*BY1/X
        ENDIF
        RETURN
        END
