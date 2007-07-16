	PROGRAM MJYNB

C       ====================================================
C       Purpose: This program computes Bessel functions
C                Jn(x) and Yn(x), and their derivatives
C                using subroutine JYNB
C       Input :  x --- Argument of Jn(x) & Yn(x)  ( x ע 0 )
C                n --- Order of Jn(x) & Yn(x)
C                      ( n = 0,1,2,תתת, n ף 250 )
C       Output:  BJ(n) --- Jn(x)
C                DJ(n) --- Jn'(x)
C                BY(n) --- Yn(x)
C                DY(n) --- Yn'(x)
C       Example:
C                x = 10.0

C                n        Jn(x)           Jn'(x)
C              -------------------------------------
C                0    -.2459358D+00   -.4347275D-01
C               10     .2074861D+00    .8436958D-01
C               20     .1151337D-04    .2011954D-04
C               30     .1551096D-11    .4396479D-11

C                n        Yn(x)           Yn'(x)
C              -------------------------------------
C                0     .5567117D-01   -.2490154D+00
C               10    -.3598142D+00    .1605149D+00
C               20    -.1597484D+04    .2737803D+04
C               30    -.7256142D+10    .2047617D+11
C       ====================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION BJ(0:250),BY(0:250),DJ(0:250),DY(0:250)
	WRITE(*,*)'  Please enter n, x'
	READ(*,*)N,X
	WRITE(*,30)N,X
	IF (N.LE.8) THEN
	   NS=1
	ELSE
	   WRITE(*,*)'  Please enter order step Ns'
	   READ(*,*)NS
	ENDIF
	WRITE(*,*)
	CALL JYNB(N,X,NM,BJ,DJ,BY,DY)
	WRITE(*,*)'  n        Jn(x)           Jn''(x)'
	WRITE(*,*)'--------------------------------------'
	DO 10 K=0,NM,NS
10         WRITE(*,40)K,BJ(K),DJ(K)
	WRITE(*,*)
	WRITE(*,*)'  n        Yn(x)           Yn''(x)'
	WRITE(*,*)'--------------------------------------'
	DO 20 K=0,NM,NS
20         WRITE(*,40)K,BY(K),DY(K)
30      FORMAT(3X,6HNmax =,I3,',    ',3Hx =,F6.1)
40      FORMAT(1X,I3,1X,2D16.7)
	END


	SUBROUTINE JYNB(N,X,NM,BJ,DJ,BY,DY)

C       =====================================================
C       Purpose: Compute Bessel functions Jn(x), Yn(x) and
C                their derivatives
C       Input :  x --- Argument of Jn(x) and Yn(x) ( x ע 0 )
C                n --- Order of Jn(x) and Yn(x)
C       Output:  BJ(n) --- Jn(x)
C                DJ(n) --- Jn'(x)
C                BY(n) --- Yn(x)
C                DY(n) --- Yn'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 to calculate the starting
C                point for backward recurrence
C       =====================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION BJ(0:N),DJ(0:N),BY(0:N),DY(0:N),
     &            A(4),B(4),A1(4),B1(4)
	PI=3.141592653589793D0
	R2P=.63661977236758D0
	NM=N
	IF (X.LT.1.0D-100) THEN
	   DO 10 K=0,N
	      BJ(K)=0.0D0
	      DJ(K)=0.0D0
	      BY(K)=-1.0D+300
10            DY(K)=1.0D+300
	   BJ(0)=1.0D0
	   DJ(1)=0.5D0
	   RETURN
	ENDIF
	IF (X.LE.300.0.OR.N.GT.INT(0.9*X)) THEN
	   IF (N.EQ.0) NM=1
	   M=MSTA1(X,200)
	   IF (M.LT.NM) THEN
	      NM=M
	   ELSE
	      M=MSTA2(X,NM,15)
	   ENDIF
	   BS=0.0D0
	   SU=0.0D0
	   SV=0.0D0
	   F2=0.0D0
	   F1=1.0D-100
	   DO 15 K=M,0,-1
	      F=2.0D0*(K+1.0D0)/X*F1-F2
	      IF (K.LE.NM) BJ(K)=F
	      IF (K.EQ.2*INT(K/2).AND.K.NE.0) THEN
		 BS=BS+2.0D0*F
		 SU=SU+(-1)**(K/2)*F/K
	      ELSE IF (K.GT.1) THEN
		 SV=SV+(-1)**(K/2)*K/(K*K-1.0)*F
	      ENDIF
	      F2=F1
15            F1=F
	   S0=BS+F
	   DO 20 K=0,NM
20            BJ(K)=BJ(K)/S0
	   EC=DLOG(X/2.0D0)+0.5772156649015329D0
	   BY0=R2P*(EC*BJ(0)-4.0D0*SU/S0)
	   BY(0)=BY0
	   BY1=R2P*((EC-1.0D0)*BJ(1)-BJ(0)/X-4.0D0*SV/S0)
	   BY(1)=BY1
	ELSE
	   DATA A/-.7031250000000000D-01,.1121520996093750D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01/
	   DATA B/ .7324218750000000D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02/
	   DATA A1/.1171875000000000D+00,-.1441955566406250D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01/
	   DATA B1/-.1025390625000000D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02/
	   T1=X-0.25D0*PI
	   P0=1.0D0
	   Q0=-0.125D0/X
	   DO 25 K=1,4
	      P0=P0+A(K)*X**(-2*K)
25            Q0=Q0+B(K)*X**(-2*K-1)
	   CU=DSQRT(R2P/X)
	   BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
	   BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
	   BJ(0)=BJ0
	   BY(0)=BY0
	   T2=X-0.75D0*PI
	   P1=1.0D0
	   Q1=0.375D0/X
	   DO 30 K=1,4
	      P1=P1+A1(K)*X**(-2*K)
30            Q1=Q1+B1(K)*X**(-2*K-1)
	   BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
	   BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
	   BJ(1)=BJ1
	   BY(1)=BY1
	   DO 35 K=2,NM
	      BJK=2.0D0*(K-1.0D0)/X*BJ1-BJ0
	      BJ(K)=BJK
	      BJ0=BJ1
35            BJ1=BJK
	ENDIF
	DJ(0)=-BJ(1)
	DO 40 K=1,NM
40         DJ(K)=BJ(K-1)-K/X*BJ(K)
	DO 45 K=2,NM
	   BYK=2.0D0*(K-1.0D0)*BY1/X-BY0
	   BY(K)=BYK
	   BY0=BY1
45         BY1=BYK
	DY(0)=-BY(1)
	DO 50 K=1,NM
50         DY(K)=BY(K-1)-K*BY(K)/X
	RETURN
	END


	INTEGER FUNCTION MSTA1(X,MP)

C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that the magnitude of
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point
C       ===================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	A0=DABS(X)
	N0=INT(1.1*A0)+1
	F0=ENVJ(N0,A0)-MP
	N1=N0+5
	F1=ENVJ(N1,A0)-MP
	DO 10 IT=1,20
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)
	   F=ENVJ(NN,A0)-MP
	   IF(ABS(NN-N1).LT.1) GO TO 20
	   N0=N1
	   F0=F1
	   N1=NN
 10        F1=F
 20     MSTA1=NN
	RETURN
	END


	INTEGER FUNCTION MSTA2(X,N,MP)

C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	A0=DABS(X)
	HMP=0.5D0*MP
	EJN=ENVJ(N,A0)
	IF (EJN.LE.HMP) THEN
	   OBJ=MP
	   N0=INT(1.1*A0)
	ELSE
	   OBJ=HMP+EJN
	   N0=N
	ENDIF
	F0=ENVJ(N0,A0)-OBJ
	N1=N0+5
	F1=ENVJ(N1,A0)-OBJ
	DO 10 IT=1,20
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)
	   F=ENVJ(NN,A0)-OBJ
	   IF (ABS(NN-N1).LT.1) GO TO 20
	   N0=N1
	   F0=F1
	   N1=NN
10         F1=F
20      MSTA2=NN+10
	RETURN
	END

	REAL*8 FUNCTION ENVJ(N,X)
	DOUBLE PRECISION X
	ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
	RETURN
	END
