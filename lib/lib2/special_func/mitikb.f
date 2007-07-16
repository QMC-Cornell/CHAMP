        PROGRAM MITIKB

C       ============================================================
C       Purpose: This program evaluates the integral of modified
C                Bessel functions I0(t) and K0(t) with respect to t
C                from 0 to x using subroutine ITIKB
C       Input :  x  --- Upper limit of the integral  ( x ò 0 )
C       Output:  TI --- Integration of I0(t) from 0 to x
C                TK --- Integration of K0(t) from 0 to x
C       Example:
C                    x         I0(t)dt         K0(t)dt
C                 -------------------------------------
C                   5.0     .318487D+02       1.567387
C                  10.0     .299305D+04       1.570779
C                  15.0     .352619D+06       1.570796
C                  20.0     .447586D+08       1.570796
C                  25.0     .589919D+10       1.570796
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x         I0(t)dt         K0(t)dt'
        WRITE(*,*)'--------------------------------------'
        CALL ITIKB(X,TI,TK)
        WRITE(*,10)X,TI,TK
10      FORMAT(1X,F5.1,D16.6,F15.6)
        END


        SUBROUTINE ITIKB(X,TI,TK)

C       =======================================================
C       Purpose: Integrate Bessel functions I0(t) and K0(t)
C                with respect to t from 0 to x
C       Input :  x  --- Upper limit of the integral ( x ò 0 )
C       Output:  TI --- Integration of I0(t) from 0 to x
C                TK --- Integration of K0(t) from 0 to x
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           TI=0.0D0
        ELSE IF (X.LT.5.0D0) THEN
           T1=X/5.0D0
           T=T1*T1
           TI=((((((((.59434D-3*T+.4500642D-2)*T
     &        +.044686921D0)*T+.300704878D0)*T+1.471860153D0)
     &        *T+4.844024624D0)*T+9.765629849D0)*T
     &        +10.416666367D0)*T+5.0D0)*T1
        ELSE IF (X.GE.5.0.AND.X.LE.8.0D0) THEN
           T=5.0D0/X
           TI=(((-.015166D0*T-.0202292D0)*T+.1294122D0)*T
     &        -.0302912D0)*T+.4161224D0
           TI=TI*DEXP(X)/DSQRT(X)
        ELSE
           T=8.0D0/X
           TI=(((((-.0073995D0*T+.017744D0)*T-.0114858D0)*T
     &        +.55956D-2)*T+.59191D-2)*T+.0311734D0)*T
     &        +.3989423D0
           TI=TI*DEXP(X)/DSQRT(X)
        ENDIF
        IF (X.EQ.0.0D0) THEN
           TK=0.0D0
        ELSE IF (X.LE.2.0D0) THEN
           T1=X/2.0D0
           T=T1*T1
           TK=((((((.116D-5*T+.2069D-4)*T+.62664D-3)*T
     &        +.01110118D0)*T+.11227902D0)*T+.50407836D0)*T
     &        +.84556868D0)*T1
              TK=TK-DLOG(X/2.0D0)*TI
        ELSE IF (X.GT.2.0.AND.X.LE.4.0D0) THEN
           T=2.0D0/X
           TK=(((.0160395D0*T-.0781715D0)*T+.185984D0)*T
     &        -.3584641D0)*T+1.2494934D0
           TK=PI/2.0D0-TK*DEXP(-X)/DSQRT(X)
        ELSE IF (X.GT.4.0.AND.X.LE.7.0D0) THEN
           T=4.0D0/X
           TK=(((((.37128D-2*T-.0158449D0)*T+.0320504D0)*T
     &        -.0481455D0)*T+.0787284D0)*T-.1958273D0)*T
     &        +1.2533141D0
           TK=PI/2.0D0-TK*DEXP(-X)/DSQRT(X)
        ELSE
           T=7.0D0/X
           TK=(((((.33934D-3*T-.163271D-2)*T+.417454D-2)*T
     &        -.933944D-2)*T+.02576646D0)*T-.11190289D0)*T
     &        +1.25331414D0
           TK=PI/2.0D0-TK*DEXP(-X)/DSQRT(X)
        ENDIF
        RETURN
        END
