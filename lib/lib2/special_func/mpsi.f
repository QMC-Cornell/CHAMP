        PROGRAM MPSI

C       ==================================================
C       Purpose: This program computes the psi function
C                using subroutine PSI
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       Examples:
C                   x          Psi(x)
C                ------------------------
C                  .25      -4.227453533
C                  .50      -1.963510026
C                  .75      -1.085860880
C                 1.00       -.577215665
C                 1.25       -.227453533
C                 1.50        .036489974
C                 1.75        .247472454
C                 2.00        .422784335
C       ==================================================

        DOUBLE PRECISION X,PS
        WRITE(*,*)'Please enter x'
        READ(*,*)X
        WRITE(*,*)'    x          Psi(x)'
        WRITE(*,*)' ------------------------'
        CALL PSI(X,PS)
        WRITE(*,10)X,PS
10      FORMAT(1X,F6.2,F18.9)
        END


        SUBROUTINE PSI(X,PS)

C       ======================================
C       Purpose: Compute the psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
