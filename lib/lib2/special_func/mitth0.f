        PROGRAM MITTH0

C       ===========================================================
C       Purpose: This program evaluates the integral of H0(t)/t
C                with respect to t from x to infinity using
C                subroutine ITTH0
C       Input :  x   --- Lower limit  ( x ò 0 )
C       Output:  TTH --- Integration of H0(t)/t from x to infinity
C       Example:
C                    x        H0(t)/t dt
C                 -----------------------
C                   0.0      1.57079633
C                   5.0       .07954575
C                  10.0       .04047175
C                  15.0       .04276558
C                  20.0       .04030796
C                  30.0       .01815256
C                  40.0       .01621331
C                  50.0       .01378661
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*) '   x       H0(t)/t dt'
        WRITE(*,*)'-----------------------'
        CALL ITTH0(X,TTH)
        WRITE(*,10)X,TTH
10      FORMAT(1X,F5.1,1X,E16.8)
        END


        SUBROUTINE ITTH0(X,TTH)

C       ===========================================================
C       Purpose: Evaluate the integral H0(t)/t with respect to t
C                from x to infinity
C       Input :  x   --- Lower limit  ( x ò 0 )
C       Output:  TTH --- Integration of H0(t)/t from x to infinity
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        S=1.0D0
        R=1.0D0
        IF (X.LT.24.5D0) THEN
           DO 10 K=1,60
              R=-R*X*X*(2.0*K-1.0D0)/(2.0*K+1.0D0)**3
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         TTH=PI/2.0D0-2.0D0/PI*X*S
        ELSE
           DO 20 K=1,10
              R=-R*(2.0*K-1.0D0)**3/((2.0*K+1.0D0)*X*X)
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 25
20            CONTINUE
25         TTH=2.0D0/(PI*X)*S
           T=8.0D0/X
           XT=X+.25D0*PI
           F0=(((((.18118D-2*T-.91909D-2)*T+.017033D0)*T
     &        -.9394D-3)*T-.051445D0)*T-.11D-5)*T+.7978846D0
           G0=(((((-.23731D-2*T+.59842D-2)*T+.24437D-2)*T
     &        -.0233178D0)*T+.595D-4)*T+.1620695D0)*T
           TTY=(F0*DSIN(XT)-G0*DCOS(XT))/(DSQRT(X)*X)
           TTH=TTH+TTY
        ENDIF
        RETURN
        END
