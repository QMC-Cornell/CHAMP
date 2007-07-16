        PROGRAM ME1Z

C       =========================================================
C       Purpose: This program computes the complex exponential
C                integral E1(z) using subroutine E1Z
C       Example:
C                     z            Re[E1(z)]       Im[E1(z)]
C                -----------------------------------------------
C                 3.0    2.0    -.90959209D-02   -.69001793D-02
C                 3.0   -2.0    -.90959209D-02    .69001793D-02
C                -3.0    2.0    -.28074891D+01    .59603353D+01
C                -3.0   -2.0    -.28074891D+01   -.59603353D+01
C                25.0   10.0    -.29302080D-12    .40391222D-12
C                25.0  -10.0    -.29302080D-12   -.40391222D-12
C               -25.0   10.0     .27279957D+10   -.49430610D+09
C               -25.0  -10.0     .27279957D+10    .49430610D+09
C       =========================================================

        IMPLICIT COMPLEX*16 (C,Z)
        IMPLICIT DOUBLE PRECISION (D-H,O-Y)
        WRITE(*,*)'Please enter x and y ( z =x+iy ) '
        READ(*,*)X,Y
        Z=CMPLX(X,Y)
        CALL E1Z(Z,CE1)
        WRITE(*,*)
        WRITE(*,*)'       z           Re[E1(z)]        Im[E1(z)]'
        WRITE(*,*)' -----------------------------------------------'
        WRITE(*,10)X,Y,CE1
10      FORMAT(1X,F5.1,2X,F5.1,1X,2D17.8)
        END


        SUBROUTINE E1Z(Z,CE1)

C       ====================================================
C       Purpose: Compute complex exponential integral E1(z)
C       Input :  z   --- Argument of E1(z)
C       Output:  CE1 --- E1(z)
C       ====================================================

        IMPLICIT COMPLEX*16 (C,Z)
        IMPLICIT DOUBLE PRECISION (D-H,O-Y)
        PI=3.141592653589793D0
        EL=0.5772156649015328D0
        X=REAL(Z)
        A0=CDABS(Z)
        IF (A0.EQ.0.0D0) THEN
           CE1=(1.0D+300,0.0D0)
        ELSE IF (A0.LE.10.0.OR.X.LT.0.0.AND.A0.LT.20.0) THEN
           CE1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,150
              CR=-CR*K*Z/(K+1.0D0)**2
              CE1=CE1+CR
              IF (CDABS(CR).LE.CDABS(CE1)*1.0D-15) GO TO 15
10         CONTINUE
15         CE1=-EL-CDLOG(Z)+Z*CE1
        ELSE
           CT0=(0.0D0,0.0D0)
           DO 20 K=120,1,-1
              CT0=K/(1.0D0+K/(Z+CT0))
20         CONTINUE
           CT=1.0D0/(Z+CT0)
           CE1=CDEXP(-Z)*CT
           IF (X.LE.0.0.AND.DIMAG(Z).EQ.0.0) CE1=CE1-PI*(0.0D0,1.0D0)
        ENDIF
        RETURN
        END
