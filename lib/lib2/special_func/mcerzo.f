        PROGRAM MCERZO

C       ===============================================================
C       Purpose : This program evaluates the complex zeros of error
C                 function erf(z) using subroutine CERZO
C       Input:    NT --- Total number of zeros
C       Example:  NT = 10

C    n     complex zeros of erf(z)     n     complex zeros of erf(z)
C   -------------------------------------------------------------------
C    1   1.450616163 + i 1.880943000   6   4.158998400 + i 4.435571444
C    2   2.244659274 + i 2.616575141   7   4.516319400 + i 4.780447644
C    3   2.839741047 + i 3.175628100   8   4.847970309 + i 5.101588043
C    4   3.335460735 + i 3.646174376   9   5.158767908 + i 5.403332643
C    5   3.769005567 + i 4.060697234  10   5.452192201 + i 5.688837437
C       ===============================================================

        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION ZO(100)
        WRITE(*,*)'Please Enter NT '
        READ(*,*)NT
        WRITE(*,20)NT
        CALL CERZO(NT,ZO)
        WRITE(*,*)'  *****    Please Wait !    *****'
        WRITE(*,*)
        WRITE(*,*)'  n        Complex zeros of erf(z)'
        WRITE(*,*)'-------------------------------------'
        DO 10 I=1,NT
10         WRITE(*,30) I,ZO(I)
20      FORMAT(2X,'NT=',I3)
30      FORMAT(1X,I3,2X,F13.8,2X,2H+i,F13.8)
        END


        SUBROUTINE CERZO(NT,ZO)

C       ===============================================================
C       Purpose : Evaluate the complex zeros of error function erf(z)
C                 using the modified Newton's iteration method
C       Input :   NT --- Total number of zeros
C       Output:   ZO(L) --- L-th zero of erf(z), L=1,2,...,NT
C       Routine called: CERF for computing erf(z) and erf'(z)
C       ===============================================================

        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION ZO(NT)
        PI=3.141592653589793D0
        DO 35 NR=1,NT
           PU=DSQRT(PI*(4.0D0*NR-0.5D0))
           PV=PI*DSQRT(2.0D0*NR-0.25D0)
           PX=0.5*PU-0.5*DLOG(PV)/PU
           PY=0.5*PU+0.5*DLOG(PV)/PU
           Z=CMPLX(PX,PY)
           IT=0
15         IT=IT+1
           CALL CERF(Z,ZF,ZD)
           ZP=(1.0D0,0.0D0)
           DO 20 I=1,NR-1
20            ZP=ZP*(Z-ZO(I))
           ZFD=ZF/ZP
           ZQ=(0.0D0,0.0D0)
           DO 30 I=1,NR-1
              ZW=(1.0D0,0.0D0)
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 ZW=ZW*(Z-ZO(J))
25            CONTINUE
30            ZQ=ZQ+ZW
           ZGD=(ZD-ZQ*ZFD)/ZP
           Z=Z-ZFD/ZGD
           W0=W
           W=CDABS(Z)
           IF (IT.LE.50.AND.DABS((W-W0)/W).GT.1.0D-11) GO TO 15
35         ZO(NR)=Z
        RETURN
        END


        SUBROUTINE CERF(Z,CER,CDER)

C       ==========================================================
C       Purpose: Compute complex Error function erf(z) & erf'(z)
C       Input:   z   --- Complex argument of erf(z)
C                x   --- Real part of z
C                y   --- Imaginary part of z
C       Output:  CER --- erf(z)
C                CDER --- erf'(z)
C       ==========================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMPLEX *16 Z,CER,CDER
        EPS=1.0D-12
        PI=3.141592653589793D0
        X=REAL(Z)
        Y=DIMAG(Z)
        X2=X*X
        IF (X.LE.3.5D0) THEN
           ER=1.0D0
           R=1.0D0
           DO 10 K=1,100
              R=R*X2/(K+0.5D0)
              ER=ER+R
              IF (DABS(ER-W).LE.EPS*DABS(ER)) GO TO 15
10            W=ER
15         C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
           ER0=C0*ER
        ELSE
           ER=1.0D0
           R=1.0D0
           DO 20 K=1,12
              R=-R*(K-0.5D0)/X2
20            ER=ER+R
           C0=DEXP(-X2)/(X*DSQRT(PI))
           ER0=1.0D0-C0*ER
        ENDIF
        IF (Y.EQ.0.0D0) THEN
           ERR=ER0
           ERI=0.0D0
        ELSE
           CS=DCOS(2.0D0*X*Y)
           SS=DSIN(2.0D0*X*Y)
           ER1=DEXP(-X2)*(1.0D0-CS)/(2.0D0*PI*X)
           EI1=DEXP(-X2)*SS/(2.0D0*PI*X)
           ER2=0.0D0
           DO 25 N=1,100
              ER2=ER2+DEXP(-.25D0*N*N)/(N*N+4.0D0*X2)*(2.0D0*X
     &            -2.0D0*X*DCOSH(N*Y)*CS+N*DSINH(N*Y)*SS)
              IF (DABS((ER2-W1)/ER2).LT.EPS) GO TO 30
25            W1=ER2
30         C0=2.0D0*DEXP(-X2)/PI
           ERR=ER0+ER1+C0*ER2
           EI2=0.0D0
           DO 35 N=1,100
              EI2=EI2+DEXP(-.25D0*N*N)/(N*N+4.0D0*X2)*(2.0D0*X
     &            *DCOSH(N*Y)*SS+N*DSINH(N*Y)*CS)
              IF (DABS((EI2-W2)/EI2).LT.EPS) GO TO 40
35            W2=EI2
40         ERI=EI1+C0*EI2
        ENDIF
        CER=CMPLX(ERR,ERI)
        CDER=2.0D0/DSQRT(PI)*CDEXP(-Z*Z)
        RETURN
        END
