        PROGRAM MITJYB

C       ===========================================================
C       Purpose: This program evaluates the integral of Bessel
C                functions J0(t) and Y0(t) with respect to t
C                from 0 to x using subroutine ITJYB
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
        WRITE(*,*)'Pleas enter x '
        READ(*,*)X
        WRITE(*,*)'   x         J0(t)dt          Y0(t)dt'
        WRITE(*,*)'---------------------------------------'
        CALL ITJYB(X,TJ,TY)
        WRITE(*,10)X,TJ,TY
10      FORMAT(1X,F5.1,2F16.8)
        END


        SUBROUTINE ITJYB(X,TJ,TY)

C       =======================================================
C       Purpose: Integrate Bessel functions J0(t) and Y0(t)
C                with respect to t from 0 to x ( x ò 0 )
C       Input :  x  --- Upper limit of the integral
C       Output:  TJ --- Integration of J0(t) from 0 to x
C                TY --- Integration of Y0(t) from 0 to x
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           TJ=0.0D0
           TY=0.0D0
        ELSE IF (X.LE.4.0D0) THEN
           X1=X/4.0D0
           T=X1*X1
           TJ=(((((((-.133718D-3*T+.2362211D-2)*T
     &        -.025791036D0)*T+.197492634D0)*T-1.015860606D0)
     &        *T+3.199997842D0)*T-5.333333161D0)*T+4.0D0)*X1
           TY=((((((((.13351D-4*T-.235002D-3)*T+.3034322D-2)*
     &        T-.029600855D0)*T+.203380298D0)*T-.904755062D0)
     &        *T+2.287317974D0)*T-2.567250468D0)*T
     &        +1.076611469D0)*X1
           TY=2.0D0/PI*DLOG(X/2.0D0)*TJ-TY
        ELSE IF (X.LE.8.0D0) THEN
           XT=X-.25D0*PI
           T=16.0D0/(X*X)
           F0=((((((.1496119D-2*T-.739083D-2)*T+.016236617D0)
     &        *T-.022007499D0)*T+.023644978D0)
     &        *T-.031280848D0)*T+.124611058D0)*4.0D0/X
           G0=(((((.1076103D-2*T-.5434851D-2)*T+.01242264D0)
     &        *T-.018255209)*T+.023664841D0)*T-.049635633D0)
     &        *T+.79784879D0
           TJ=1.0D0-(F0*DCOS(XT)-G0*DSIN(XT))/DSQRT(X)
           TY=-(F0*DSIN(XT)+G0*DCOS(XT))/DSQRT(X)
        ELSE
           T=64.0D0/(X*X)
           XT=X-.25D0*PI
           F0=(((((((-.268482D-4*T+.1270039D-3)*T
     &        -.2755037D-3)*T+.3992825D-3)*T-.5366169D-3)*T
     &        +.10089872D-2)*T-.40403539D-2)*T+.0623347304D0)
     &        *8.0D0/X
           G0=((((((-.226238D-4*T+.1107299D-3)*T-.2543955D-3)
     &        *T+.4100676D-3)*T-.6740148D-3)*T+.17870944D-2)
     &        *T-.01256424405D0)*T+.79788456D0
           TJ=1.0D0-(F0*DCOS(XT)-G0*DSIN(XT))/DSQRT(X)
           TY=-(F0*DSIN(XT)+G0*DCOS(XT))/DSQRT(X)
        ENDIF
        RETURN
        END

