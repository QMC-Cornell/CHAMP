        PROGRAM MITTJYB

C       ===========================================================
C       Purpose: This program computes the integral of [1-J0(t)]/t
C                with respect to t from 0 to x and Y0(t)/t with
C                respect to t from x to ì using subroutine ITTJYB
C       Input :  x   --- Variable in the limits  ( x ò 0 )
C       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
C                TTY --- Integration of Y0(t)/t from x to ì
C       Example:
C                  x      [1-J0(t)]/tdt       Y0(t)/tdt
C                ----------------------------------------
C                 5.0     .1540347D+01    -.4632208D-01
C                10.0     .2177866D+01    -.2298791D-01
C                15.0     .2578551D+01     .3857453D-03
C                20.0     .2877311D+01     .8503154D-02
C                25.0     .3108231D+01     .3526339D-02
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x     [1-J0(t)]/tdt      Y0(t)/tdt'
        WRITE(*,*)'----------------------------------------'
        CALL ITTJYB(X,TTJ,TTY)
        WRITE(*,10)X,TTJ,TTY
10      FORMAT(1X,F5.1,2D17.7)
        END


        SUBROUTINE ITTJYB(X,TTJ,TTY)

C       ==========================================================
C       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
C                to x, and Y0(t)/t with respect to t from x to ì
C       Input :  x   --- Variable in the limits  ( x ò 0 )
C       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
C                TTY --- Integration of Y0(t)/t from x to ì
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IF (X.EQ.0.0D0) THEN
           TTJ=0.0D0
           TTY=-1.0D+300
        ELSE IF (X.LE.4.0D0) THEN
           X1=X/4.0D0
           T=X1*X1
           TTJ=((((((.35817D-4*T-.639765D-3)*T+.7092535D-2)*T
     &         -.055544803D0)*T+.296292677D0)*T-.999999326D0)
     &         *T+1.999999936D0)*T
           TTY=(((((((-.3546D-5*T+.76217D-4)*T-.1059499D-2)*T
     &         +.010787555D0)*T-.07810271D0)*T+.377255736D0)
     &         *T-1.114084491D0)*T+1.909859297D0)*T
           E0=EL+DLOG(X/2.0D0)
           TTY=PI/6.0D0+E0/PI*(2.0D0*TTJ-E0)-TTY
        ELSE IF (X.LE.8.0D0) THEN
           XT=X+.25D0*PI
           T1=4.0D0/X
           T=T1*T1
           F0=(((((.0145369D0*T-.0666297D0)*T+.1341551D0)*T
     &        -.1647797D0)*T+.1608874D0)*T-.2021547D0)*T
     &        +.7977506D0
           G0=((((((.0160672D0*T-.0759339D0)*T+.1576116D0)*T
     &        -.1960154D0)*T+.1797457D0)*T-.1702778D0)*T
     &        +.3235819D0)*T1
           TTJ=(F0*DCOS(XT)+G0*DSIN(XT))/(DSQRT(X)*X)
           TTJ=TTJ+EL+DLOG(X/2.0D0)
           TTY=(F0*DSIN(XT)-G0*DCOS(XT))/(DSQRT(X)*X)
        ELSE
           T=8.0D0/X
           XT=X+.25D0*PI
           F0=(((((.18118D-2*T-.91909D-2)*T+.017033D0)*T
     &        -.9394D-3)*T-.051445D0)*T-.11D-5)*T+.7978846D0
           G0=(((((-.23731D-2*T+.59842D-2)*T+.24437D-2)*T
     &      -.0233178D0)*T+.595D-4)*T+.1620695D0)*T
           TTJ=(F0*DCOS(XT)+G0*DSIN(XT))/(DSQRT(X)*X)
     &         +EL+DLOG(X/2.0D0)
           TTY=(F0*DSIN(XT)-G0*DCOS(XT))/(DSQRT(X)*X)
        ENDIF
        RETURN
        END
