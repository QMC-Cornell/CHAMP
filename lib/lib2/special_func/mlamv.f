	PROGRAM MLAMV

C       =======================================================
C       Purpose: This program computes the lambda functions
C                for an arbitrary order, and their derivative
C                using subroutine LAMV
C       Input :  x --- Argument of lambda function
C                v --- Order of lambda function
C                      ( v = n+v0, 0 ó n ó 250, 0 ó v0 < 1 )
C       Output:  VL(n) --- Lambda function of order n+v0
C                DL(n) --- Derivative of lambda function
C       Example: x = 10.0

C                   v         Lambda(x)        Lambda'(x)
C                 ------------------------------------------
C                  0.25    -.12510515D+00    -.78558916D-01
C                  0.50    -.54402111D-01    -.78466942D-01
C                  0.75    -.13657787D-01    -.66234027D-01
C                  1.00     .86945492D-02    -.50926063D-01
C                  1.25     .19639729D-01    -.36186221D-01
C                  1.50     .23540083D-01    -.23382658D-01
C                  1.75     .23181910D-01    -.12893894D-01
C                  2.00     .20370425D-01    -.46703503D-02
C                  2.25     .16283799D-01     .15101684D-02
C                  2.50     .11691329D-01     .59243767D-02
C       =======================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION VL(0:250),DL(0:250)
	WRITE(*,*)'  Please enter v and x '
	READ(*,*)V,X
	WRITE(*,20)V,X
	IF (V.LE.8) THEN
	   NS=1
	ELSE
	   WRITE(*,*)'  Please enter order step Ns'
	   READ(*,*)NS
	ENDIF
	WRITE(*,*)
	WRITE(*,*) '   v         Lambda(x)        Lambda''(x)'
	WRITE(*,*)'-------------------------------------------'
	CALL LAMV(V,X,VM,VL,DL)
	NM=INT(VM)
	V0=VM-NM
	DO 10 K=0,NM,NS
	   VK=K+V0
10         WRITE(*,15)VK,VL(K),DL(K)
15      FORMAT(1X,F6.2,2D18.8)
20      FORMAT(1X,'v =',F6.2,'    ','x =',F8.2)
	END


	SUBROUTINE LAMV(V,X,VM,VL,DL)

C       =========================================================
C       Purpose: Compute lambda function with arbitrary order v,
C                and their derivative
C       Input :  x --- Argument of lambda function
C                v --- Order of lambda function
C       Output:  VL(n) --- Lambda function of order n+v0
C                DL(n) --- Derivative of lambda function
C                VM --- Highest order computed
C       Routines called:
C            (1) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C            (2) GAM0 for computing gamma function (|x| ó 1)
C       =========================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION VL(0:*),DL(0:*)
	PI=3.141592653589793D0
	RP2=0.63661977236758D0
	X=DABS(X)
	X2=X*X
	N=INT(V)
	V0=V-N
	VM=V
	IF (X.LE.12.0D0) THEN
	   DO 25 K=0,N
	      VK=V0+K
	      BK=1.0D0
	      R=1.0D0
	      DO 10 I=1,50
		 R=-0.25D0*R*X2/(I*(I+VK))
		 BK=BK+R
		 IF (DABS(R).LT.DABS(BK)*1.0D-15) GO TO 15
10            CONTINUE
15            VL(K)=BK
	      UK=1.0D0
	      R=1.0D0
	      DO 20 I=1,50
		 R=-0.25D0*R*X2/(I*(I+VK+1.0D0))
		 UK=UK+R
		 IF (DABS(R).LT.DABS(UK)*1.0D-15) GO TO 25
20            CONTINUE
25            DL(K)=-0.5D0*X/(VK+1.0D0)*UK
	   RETURN
	ENDIF
	K0=11
	IF (X.GE.35.0D0) K0=10
	IF (X.GE.50.0D0) K0=8
	DO 40 J=0,1
	   VV=4.0D0*(J+V0)*(J+V0)
	   PX=1.0D0
	   RP=1.0D0
	   DO 30 K=1,K0
	      RP=-0.78125D-2*RP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &            (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*X2)
30            PX=PX+RP
	   QX=1.0D0
	   RQ=1.0D0
	   DO 35 K=1,K0
	      RQ=-0.78125D-2*RQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &            (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*X2)
35            QX=QX+RQ
	   QX=0.125D0*(VV-1.0D0)*QX/X
	   XK=X-(0.5D0*(J+V0)+0.25D0)*PI
	   A0=DSQRT(RP2/X)
	   CK=DCOS(XK)
	   SK=DSIN(XK)
	   IF (J.EQ.0) BJV0=A0*(PX*CK-QX*SK)
	   IF (J.EQ.1) BJV1=A0*(PX*CK-QX*SK)
40      CONTINUE
	IF (V0.EQ.0.0D0) THEN
	   GA=1.0D0
	ELSE
	   CALL GAM0(V0,GA)
	   GA=V0*GA
	ENDIF
	FAC=(2.0D0/X)**V0*GA
	VL(0)=BJV0
	DL(0)=-BJV1+V0/X*BJV0
	VL(1)=BJV1
	DL(1)=BJV0-(1.0D0+V0)/X*BJV1
	R0=2.0D0*(1.0D0+V0)/X
	IF (N.LE.1) THEN
	   VL(0)=FAC*VL(0)
	   DL(0)=FAC*DL(0)-V0/X*VL(0)
	   VL(1)=FAC*R0*VL(1)
	   DL(1)=FAC*R0*DL(1)-(1.0D0+V0)/X*VL(1)
	   RETURN
	ENDIF
	IF (N.GE.2.AND.N.LE.INT(0.9*X)) THEN
	   F0=BJV0
	   F1=BJV1
	   DO 45 K=2,N
	      F=2.0D0*(K+V0-1.0D0)/X*F1-F0
	      F0=F1
	      F1=F
45            VL(K)=F
	ELSE IF (N.GE.2) THEN
	   M=MSTA1(X,200)
	   IF (M.LT.N) THEN
	      N=M
	   ELSE
	      M=MSTA2(X,N,15)
	   ENDIF
	   F2=0.0D0
	   F1=1.0D-100
	   DO 50 K=M,0,-1
	      F=2.0D0*(V0+K+1.0D0)/X*F1-F2
	      IF (K.LE.N) VL(K)=F
	      F2=F1
50            F1=F
	   IF (DABS(BJV0).GT.DABS(BJV1)) CS=BJV0/F
	   IF (DABS(BJV0).LE.DABS(BJV1)) CS=BJV1/F2
	   DO 55 K=0,N
55            VL(K)=CS*VL(K)
	ENDIF
	VL(0)=FAC*VL(0)
	DO 65 J=1,N
	   RC=FAC*R0
	   VL(J)=RC*VL(J)
	   DL(J-1)=-0.5D0*X/(J+V0)*VL(J)
65         R0=2.0D0*(J+V0+1)/X*R0
	DL(N)=2.0D0*(V0+N)*(VL(N-1)-VL(N))/X
	VM=N+V0
	RETURN
	END


	SUBROUTINE GAM0 (X,GA)

C       ================================================
C       Purpose: Compute gamma function â(x)
C       Input :  x  --- Argument of â(x)  ( |x| ó 1 )
C       Output:  GA --- â(x)
C       ================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION G(25)
	DATA G/1.0D0,0.5772156649015329D0,
     &       -0.6558780715202538D0, -0.420026350340952D-1,
     &        0.1665386113822915D0, -.421977345555443D-1,
     &        -.96219715278770D-2, .72189432466630D-2,
     &        -.11651675918591D-2, -.2152416741149D-3,
     &         .1280502823882D-3, -.201348547807D-4,
     &        -.12504934821D-5, .11330272320D-5,
     &        -.2056338417D-6, .61160950D-8,
     &         .50020075D-8, -.11812746D-8,
     &         .1043427D-9, .77823D-11,
     &        -.36968D-11, .51D-12,
     &        -.206D-13, -.54D-14, .14D-14/
	GR=(25)
	DO 20 K=24,1,-1
20         GR=GR*X+G(K)
	GA=1.0D0/(GR*X)
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
