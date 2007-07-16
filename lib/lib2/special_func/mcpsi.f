        PROGRAM MCPSI

C       =========================================================
C       Purpose: This program computes the psi function psi(z)
C                for a complex argument using subroutine CPSI
C       Input :  x   --- Real part of z
C                y   --- Imaginary part of z
C       Output:  PSR --- Real part of psi(z)
C                PSI --- Imaginary part of psi(z)
C       Examples:
C                   x       y      Re[psi(z)]     Im[psi(z)]
C                 -------------------------------------------
C                  3.0     2.0     1.16459152      .67080728
C                  3.0    -2.0     1.16459152     -.67080728
C                 -3.0     2.0     1.39536075     2.62465344
C                 -3.0    -2.0     1.39536075    -2.62465344
C       =========================================================

        DOUBLE PRECISION X,Y,PSR,PSI
        WRITE(*,*)'Please enter x and y ( z=x+iy)'
        READ(*,*)X,Y
        WRITE(*,20)X,Y
        WRITE(*,*)
        WRITE(*,*)'   x       y      Re[Psi(z)]      Im[Psi(z)]'
        WRITE(*,*)' ----------------------------------------------'
        CALL CPSI(X,Y,PSR,PSI)
        WRITE(*,10)X,Y,PSR,PSI
10      FORMAT(1X,F5.1,3X,F5.1,2X,2E16.8)
20      FORMAT(1X,2Hx=,F6.2,6X,2Hy=,F6.2)
        END


        SUBROUTINE CPSI(X,Y,PSR,PSI)

C       =============================================
C       Purpose: Compute the psi function for a
C                complex argument
C       Input :  x   --- Real part of z
C                y   --- Imaginary part of z
C       Output:  PSR --- Real part of psi(z)
C                PSI --- Imaginary part of psi(z)
C       =============================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(8)
        DATA A/-.8333333333333D-01,.83333333333333333D-02,
     &       -.39682539682539683D-02,.41666666666666667D-02,
     &       -.75757575757575758D-02,.21092796092796093D-01,
     &       -.83333333333333333D-01,.4432598039215686D0/
        PI=3.141592653589793D0
        IF (Y.EQ.0.0D0.AND.X.EQ.INT(X).AND.X.LE.0.0D0) THEN
           PSR=1.0D+300
           PSI=0.0D0
        ELSE
           IF (X.LT.0.0D0) THEN
              X1=X
              Y1=Y
              X=-X
              Y=-Y
           ENDIF
           X0=X
           IF (X.LT.8.0D0) THEN
              N=8-INT(X)
              X0=X+N
           ENDIF
           IF (X0.EQ.0.0D0.AND.Y.NE.0.0D0) TH=0.5D0*PI
           IF (X0.NE.0.0D0) TH=DATAN(Y/X0)
           Z2=X0*X0+Y*Y
           Z0=DSQRT(Z2)
           PSR=DLOG(Z0)-0.5D0*X0/Z2
           PSI=TH+0.5D0*Y/Z2
           DO 10 K=1,8
              PSR=PSR+A(K)*Z2**(-K)*DCOS(2.0D0*K*TH)
10            PSI=PSI-A(K)*Z2**(-K)*DSIN(2.0D0*K*TH)
           IF (X.LT.8.0D0) THEN
              RR=0.0D0
              RI=0.0D0
              DO 20 K=1,N
                 RR=RR+(X0-K)/((X0-K)**2.0D0+Y*Y)
20               RI=RI+Y/((X0-K)**2.0D0+Y*Y)
              PSR=PSR-RR
              PSI=PSI+RI
           ENDIF
           IF (X1.LT.0.0D0) THEN
              TN=DTAN(PI*X)
              TM=DTANH(PI*Y)
              CT2=TN*TN+TM*TM
              PSR=PSR+X/(X*X+Y*Y)+PI*(TN-TN*TM*TM)/CT2
              PSI=PSI-Y/(X*X+Y*Y)-PI*TM*(1.0D0+TN*TN)/CT2
              X=X1
              Y=Y1
           ENDIF
        ENDIF
        RETURN
        END
