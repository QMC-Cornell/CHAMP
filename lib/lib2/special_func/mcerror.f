	PROGRAM MCERROR

C       ============================================================
C       Purpose: This program computes the error function erf(z)
C                for a complex argument using subroutine CERROR
C       Input :  x   --- Real part of z
C                y   --- Imaginary part of z  ( y ó 3.0 )
C       Output:  ERR --- Real part of erf(z)
C                ERI --- Imaginary part of erf(z)
C       Example:
C                   x       y       Re[erf(z)]      Im[erf(z)]
C                 ---------------------------------------------
C                  1.0     2.0      -.53664357     -5.04914370
C                  2.0     2.0      1.15131087       .12729163
C                  3.0     2.0       .99896328      -.00001155
C                  4.0     2.0      1.00000057      -.00000051
C                  5.0     2.0      1.00000000       .00000000
C       ============================================================

	IMPLICIT COMPLEX *16 (C,Z)
	DOUBLE PRECISION X,Y
	WRITE(*,*)'X,Y=?'
	READ(*,*)X,Y
	WRITE(*,*)'   x      y      Re[erf(z)]      Im[erf(z)]'
	WRITE(*,*)' ---------------------------------------------'
	Z=CMPLX(X,Y)
	CALL CERROR(Z,CER)
	WRITE(*,10) Z,CER
10      FORMAT(1X,F5.1,2X,F5.1,1X,2E16.8)
	END


	SUBROUTINE CERROR(Z,CER)

C       ====================================================
C       Purpose: Compute error function erf(z) for a complex
C                argument (z=x+iy)
C       Input :  z   --- Complex argument
C       Output:  CER --- erf(z)
C       ====================================================

	IMPLICIT COMPLEX *16 (C,Z)
	DOUBLE PRECISION A0,PI
	A0=CDABS(Z)
	C0=CDEXP(-Z*Z)
	PI=3.141592653589793D0
	Z1=Z
	IF (REAL(Z).LT.0.0) THEN
	   Z1=-Z
	ENDIF
	IF (A0.LE.5.8D0) THEN
	   CS=Z1
	   CR=Z1
	   DO 10 K=1,120
	      CR=CR*Z1*Z1/(K+0.5D0)
	      CS=CS+CR
	      IF (CDABS(CR/CS).LT.1.0D-15) GO TO 15
10         CONTINUE
15         CER=2.0D0*C0*CS/DSQRT(PI)
	ELSE
	   CL=1.0D0/Z1
	   CR=CL
	   DO 20 K=1,13
	      CR=-CR*(K-0.5D0)/(Z1*Z1)
	      CL=CL+CR
	      IF (CDABS(CR/CL).LT.1.0D-15) GO TO 25
20         CONTINUE
25         CER=1.0D0-C0*CL/DSQRT(PI)
	ENDIF
	IF (REAL(Z).LT.0.0) THEN
	   CER=-CER
	ENDIF
	RETURN
	END


