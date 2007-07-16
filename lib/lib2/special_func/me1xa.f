        PROGRAM ME1XA

C       =========================================================
C       Purpose: This program computes the exponential integral
C                E1(x) using subroutine E1XA
C       Input :  x  --- Argument of E1(x)  ( x > 0 )
C       Output:  E1 --- E1(x)
C       Example:
C                  x        E1(x)
C                ----------------------
C                 0.0     .1000000+301
C                 1.0     .2193839E+00
C                 2.0     .4890051E-01
C                 3.0     .1304838E-01
C                 4.0     .3779352E-02
C                 5.0     .1148296E-02
C       =========================================================

        DOUBLE PRECISION E1,X
        WRITE(*,*)'Please enter x '
        READ(*,*) X
        WRITE(*,*)'   x        E1(x)'
        WRITE(*,*)' ----------------------'
        CALL E1XA(X,E1)
        WRITE(*,10)X,E1
10      FORMAT(1X,F5.1,E17.7)
        END


        SUBROUTINE E1XA(X,E1)

C       ============================================
C       Purpose: Compute exponential integral E1(x)
C       Input :  x  --- Argument of E1(x)
C       Output:  E1 --- E1(x) ( x > 0 )
C       ============================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0) THEN
           E1=1.0D+300
        ELSE IF (X.LE.1.0) THEN
           E1=-DLOG(X)+((((1.07857D-3*X-9.76004D-3)*X+5.519968D-2)*X
     &        -0.24991055D0)*X+0.99999193D0)*X-0.57721566D0
        ELSE
           ES1=(((X+8.5733287401D0)*X+18.059016973D0)*X
     &         +8.6347608925D0)*X+0.2677737343D0
           ES2=(((X+9.5733223454D0)*X+25.6329561486D0)*X
     &         +21.0996530827D0)*X+3.9584969228D0
           E1=DEXP(-X)/X*ES1/ES2
        ENDIF
        RETURN
        END
