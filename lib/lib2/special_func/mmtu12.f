	PROGRAM MMTU12

C       ===============================================================
C       Purpose: This program computes the modified Mathieu functions
C                of the first and second kinds, Mcm(1)(2)(x,q) and
C                Msm(1)(2)(x,q), and their derivatives using
C                subroutine MTU12
C       Input:   KF --- Function code
C                       KF=1 for computing Mcm(x,q)
C                       KF=2 for computing Msm(x,q)
C                KC --- Function Code
C                       KC=1 for computing Mcm(1)(x,q) and Mcm(1)'(x,q)
C                            or Msm(1)(x,q) and Msm(1)'(x,q)
C                       KC=2 for computing Mcm(2)(x,q) and Mcm(2)'(x,q)
C                            or Msm(2)(x,q) and Msm(2)'(x,q)
C                       KC=3 for both modified Mathieu functions of the
C                            first and second kinds, and their
C                            derivatives
C                m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                x  --- Argument of Mathieu functions
C       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
C                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
C                F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
C                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
C       ===============================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	WRITE(*,*)'Please enter KF, m, q and x '
	READ(*,*)KF,M,Q,X
	WRITE(*,10)KF,M,Q,X
	KC=3
	CALL MTU12(KF,KC,M,Q,X,F1R,D1R,F2R,D2R)
	WRITE(*,*)
	IF (KF.EQ.1) THEN
	   WRITE(*,*)'   x      Mcm(1)(x,q)    Mcm(1)''(x,q)',
     &               '    Mcm(2)(x,q)     Mcm(2)''(x,q)'
	ELSE
	   WRITE(*,*)'   x      Msm(1)(x,q)    Msm(1)''(x,q)',
     &               '    Msm(2)(x,q)     Msm(2)''(x,q)'
	ENDIF
	WRITE(*,*)' --------------------------------------',
     &            '-------------------------------'
	WRITE(*,20)X,F1R,D1R,F2R,D2R
	WRITE(*,*)
	WRITE(*,30)F1R*D2R-F2R*D1R,.63661977236758D0
	WRITE(*,40)
10      FORMAT(1X,4HKF =,I2,',  ',3Hm =,I3,',  ',
     &         3Hq =,F5.1,',  ',3Hx =,F5.1)
20      FORMAT(1X,F5.1,4D16.8)
30      FORMAT(1X,'WRONSKIAN=',E16.8,3X,'should equal   2/PI=',E16.8)
40      FORMAT(1X,/1X,'Caution: This check is not accurate if it ',
     &        'involves',/1X,'         the subtraction of two ',
     &        'similar numbers')
	END


	SUBROUTINE MTU12(KF,KC,M,Q,X,F1R,D1R,F2R,D2R)

C       ==============================================================
C       Purpose: Compute modified Mathieu functions of the first and
C                second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
C                and their derivatives
C       Input:   KF --- Function code
C                       KF=1 for computing Mcm(x,q)
C                       KF=2 for computing Msm(x,q)
C                KC --- Function Code
C                       KC=1 for computing the first kind
C                       KC=2 for computing the second kind
C                            or Msm(2)(x,q) and Msm(2)'(x,q)
C                       KC=3 for computing both the first
C                            and second kinds
C                m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions ( q ò 0 )
C                x  --- Argument of Mathieu functions
C       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
C                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
C                F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
C                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
C       Routines called:
C            (1) CVA2 for computing the characteristic values
C            (2) FCOEF for computing expansion coefficients
C            (3) JYNB for computing Jn(x), Yn(x) and their
C                derivatives
C       ==============================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION FG(251),BJ1(0:251),DJ1(0:251),BJ2(0:251),DJ2(0:251),
     &            BY1(0:251),DY1(0:251),BY2(0:251),DY2(0:251)
	EPS=1.0D-14
	IF (KF.EQ.1.AND.M.EQ.2*INT(M/2)) KD=1
	IF (KF.EQ.1.AND.M.NE.2*INT(M/2)) KD=2
	IF (KF.EQ.2.AND.M.NE.2*INT(M/2)) KD=3
	IF (KF.EQ.2.AND.M.EQ.2*INT(M/2)) KD=4
	CALL CVA2(KD,M,Q,A)
	IF (Q.LE.1.0D0) THEN
	   QM=7.5+56.1*SQRT(Q)-134.7*Q+90.7*SQRT(Q)*Q
	ELSE
	   QM=17.0+3.1*SQRT(Q)-.126*Q+.0037*SQRT(Q)*Q
	ENDIF
	KM=INT(QM+0.5*M)
	CALL FCOEF(KD,M,Q,A,FG)
	IC=INT(M/2)+1
	IF (KD.EQ.4) IC=M/2
	C1=DEXP(-X)
	C2=DEXP(X)
	U1=DSQRT(Q)*C1
	U2=DSQRT(Q)*C2
	CALL JYNB(KM,U1,NM,BJ1,DJ1,BY1,DY1)
	CALL JYNB(KM,U2,NM,BJ2,DJ2,BY2,DY2)
	IF (KC.EQ.2) GO TO 50
	F1R=0.0D0
	DO 30 K=1,KM
	   IF (KD.EQ.1) THEN
	      F1R=F1R+(-1)**(IC+K)*FG(K)*BJ1(K-1)*BJ2(K-1)
	   ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
	      F1R=F1R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BJ2(K)
     &            +(-1)**KD*BJ1(K)*BJ2(K-1))
	   ELSE
	      F1R=F1R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BJ2(K+1)
     &            -BJ1(K+1)*BJ2(K-1))
	   ENDIF
	   IF (K.GE.5.AND.DABS(F1R-W1).LT.DABS(F1R)*EPS) GO TO 35
30         W1=F1R
35      F1R=F1R/FG(1)
	D1R=0.0D0
	DO 40 K=1,KM
	   IF (KD.EQ.1) THEN
	      D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*BJ1(K-1)*DJ2(K-1)
     &            -C1*DJ1(K-1)*BJ2(K-1))
	   ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
	      D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DJ2(K)
     &            +(-1)**KD*BJ1(K)*DJ2(K-1))-C1*(DJ1(K-1)*BJ2(K)
     &            +(-1)**KD*DJ1(K)*BJ2(K-1)))
	   ELSE
	      D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DJ2(K+1)
     &            -BJ1(K+1)*DJ2(K-1))-C1*(DJ1(K-1)*BJ2(K+1)
     &            -DJ1(K+1)*BJ2(K-1)))
	   ENDIF
	   IF (K.GE.5.AND.DABS(D1R-W2).LT.DABS(D1R)*EPS) GO TO 45
40         W2=D1R
45      D1R=D1R*DSQRT(Q)/FG(1)
	IF (KC.EQ.1) RETURN
50      F2R=0.0D0
	DO 55 K=1,KM
	   IF (KD.EQ.1) THEN
	      F2R=F2R+(-1)**(IC+K)*FG(K)*BJ1(K-1)*BY2(K-1)
	   ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
	      F2R=F2R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BY2(K)
     &            +(-1)**KD*BJ1(K)*BY2(K-1))
	   ELSE
	      F2R=F2R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BY2(K+1)
     &            -BJ1(K+1)*BY2(K-1))
	   ENDIF
	   IF (K.GE.5.AND.DABS(F2R-W1).LT.DABS(F2R)*EPS) GO TO 60
55         W1=F2R
60      F2R=F2R/FG(1)
	D2R=0.0D0
	DO 65 K=1,KM
	   IF (KD.EQ.1) THEN
	      D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*BJ1(K-1)*DY2(K-1)
     &            -C1*DJ1(K-1)*BY2(K-1))
	   ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
	      D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DY2(K)
     &            +(-1)**KD*BJ1(K)*DY2(K-1))-C1*(DJ1(K-1)*BY2(K)
     &            +(-1)**KD*DJ1(K)*BY2(K-1)))
	   ELSE
	      D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DY2(K+1)
     &            -BJ1(K+1)*DY2(K-1))-C1*(DJ1(K-1)*BY2(K+1)
     &            -DJ1(K+1)*BY2(K-1)))
	   ENDIF
	   IF (K.GE.5.AND.DABS(D2R-W2).LT.DABS(D2R)*EPS) GO TO 70
65         W2=D2R
70         D2R=D2R*DSQRT(Q)/FG(1)
	RETURN
	END


	SUBROUTINE FCOEF(KD,M,Q,A,FC)

C       =====================================================
C       Purpose: Compute expansion coefficients for Mathieu
C                functions and modified Mathieu functions
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C                A  --- Characteristic value of Mathieu
C                       functions for given m and q
C       Output:  FC(k) --- Expansion coefficients of Mathieu
C                       functions ( k= 1,2,...,KM )
C                       FC(1),FC(2),FC(3),... correspond to
C                       A0,A2,A4,... for KD=1 case, A1,A3,
C                       A5,... for KD=2 case, B1,B3,B5,...
C                       for KD=3 case and B2,B4,B6,... for
C                       KD=4 case
C       =====================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION FC(251)
	IF (Q.LE.1.0D0) THEN
	   QM=7.5+56.1*SQRT(Q)-134.7*Q+90.7*SQRT(Q)*Q
	ELSE
	   QM=17.0+3.1*SQRT(Q)-.126*Q+.0037*SQRT(Q)*Q
	ENDIF
	KM=INT(QM+0.5*M)
	IF (Q.EQ.0.0D0) THEN
	   DO 10 K=1,KM
10            FC(K)=0.0D0
	   IF (KD.EQ.1) THEN
	      FC((M+2)/2)=1.0D0
	      IF (M.EQ.0) FC(1)=1.0D0/DSQRT(2.0D0)
	   ELSE IF (KD.EQ.4) THEN
	      FC(M/2)=1.0D0
	   ELSE
	      FC((M+1)/2)=1.0D0
	   ENDIF
	   RETURN
	ENDIF
	KB=0
	S=0.0D0
	F=1.0D-100
	U=0.0D0
	FC(KM)=0.0D0
	IF (KD.EQ.1) THEN
	   DO 25 K=KM,3,-1
	      V=U
	      U=F
	      F=(A-4.0D0*K*K)*U/Q-V
	      IF (DABS(F).LT.DABS(FC(K+1))) THEN
		 KB=K
		 FC(1)=1.0D-100
		 SP=0.0D0
		 F3=FC(K+1)
		 FC(2)=A/Q*FC(1)
		 FC(3)=(A-4.0D0)*FC(2)/Q-2.0D0*FC(1)
		 U=FC(2)
		 F1=FC(3)
		 DO 15 I=3,KB
		    V=U
		    U=F1
		    F1=(A-4.0D0*(I-1.0D0)**2)*U/Q-V
		    FC(I+1)=F1
		    IF (I.EQ.KB) F2=F1
		    IF (I.NE.KB) SP=SP+F1*F1
15               CONTINUE
		 SP=SP+2.0D0*FC(1)**2+FC(2)**2+FC(3)**2
		 SS=S+SP*(F3/F2)**2
		 S0=DSQRT(1.0D0/SS)
		 DO 20 J=1,KM
		    IF (J.LE.KB+1) THEN
		       FC(J)=S0*FC(J)*F3/F2
		    ELSE
		       FC(J)=S0*FC(J)
		    ENDIF
20               CONTINUE
		 GO TO 85
	      ELSE
		 FC(K)=F
		 S=S+F*F
	      ENDIF
25         CONTINUE
	   FC(2)=Q*FC(3)/(A-4.0D0-2.0D0*Q*Q/A)
	   FC(1)=Q/A*FC(2)
	   S=S+2.0D0*FC(1)**2+FC(2)**2
	   S0=DSQRT(1.0D0/S)
	   DO 30 K=1,KM
30            FC(K)=S0*FC(K)
	ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
	   DO 35 K=KM,3,-1
	      V=U
	      U=F
	      F=(A-(2.0D0*K-1)**2)*U/Q-V
	      IF (DABS(F).GE.DABS(FC(K))) THEN
		 FC(K-1)=F
		 S=S+F*F
	      ELSE
		 KB=K
		 F3=FC(K)
		 GO TO 45
	      ENDIF
35         CONTINUE
	   FC(1)=Q/(A-1.0D0-(-1)**KD*Q)*FC(2)
	   S=S+FC(1)*FC(1)
	   S0=DSQRT(1.0D0/S)
	   DO 40 K=1,KM
40            FC(K)=S0*FC(K)
	   GO TO 85
45         FC(1)=1.0D-100
	   FC(2)=(A-1.0D0-(-1)**KD*Q)/Q*FC(1)
	   SP=0.0D0
	   U=FC(1)
	   F1=FC(2)
	   DO 50 I=2,KB-1
	      V=U
	      U=F1
	      F1=(A-(2.0D0*I-1.0D0)**2)*U/Q-V
	      IF (I.NE.KB-1) THEN
		 FC(I+1)=F1
		 SP=SP+F1*F1
	      ELSE
		 F2=F1
	      ENDIF
50         CONTINUE
	   SP=SP+FC(1)**2+FC(2)**2
	   SS=S+SP*(F3/F2)**2
	   S0=1.0D0/DSQRT(SS)
	   DO 55 J=1,KM
	      IF (J.LT.KB) FC(J)=S0*FC(J)*F3/F2
	      IF (J.GE.KB) FC(J)=S0*FC(J)
55         CONTINUE
	ELSE IF (KD.EQ.4) THEN
	   DO 60 K=KM,3,-1
	      V=U
	      U=F
	      F=(A-4.0D0*K*K)*U/Q-V
	      IF (DABS(F).GE.DABS(FC(K))) THEN
		 FC(K-1)=F
		 S=S+F*F
	      ELSE
		 KB=K
		 F3=FC(K)
		 GO TO 70
	      ENDIF
60         CONTINUE
	   FC(1)=Q/(A-4.0D0)*FC(2)
	   S=S+FC(1)*FC(1)
	   S0=DSQRT(1.0D0/S)
	   DO 65 K=1,KM
65            FC(K)=S0*FC(K)
	   GO TO 85
70         FC(1)=1.0D-100
	   FC(2)=(A-4.0D0)/Q*FC(1)
	   SP=0.0D0
	   U=FC(1)
	   F1=FC(2)
	   DO 75 I=2,KB-1
	      V=U
	      U=F1
	      F1=(A-4.0D0*I*I)*U/Q-V
	      IF (I.NE.KB-1) THEN
		 FC(I+1)=F1
		 SP=SP+F1*F1
	      ELSE
		 F2=F1
	      ENDIF
75         CONTINUE
	   SP=SP+FC(1)**2+FC(2)**2
	   SS=S+SP*(F3/F2)**2
	   S0=1.0D0/DSQRT(SS)
	   DO 80 J=1,KM
	      IF (J.LT.KB) FC(J)=S0*FC(J)*F3/F2
	      IF (J.GE.KB) FC(J)=S0*FC(J)
80         CONTINUE
	ENDIF
85      IF (FC(1).LT.0.0D0) THEN
	   DO 90 J=1,KM
90            FC(J)=-FC(J)
	ENDIF
	RETURN
	END


	SUBROUTINE CVA2(KD,M,Q,A)

C       ======================================================
C       Purpose: Calculate a specific characteristic value of
C                Mathieu functions
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C       Output:  A  --- Characteristic value
C       Routines called:
C             (1) REFINE for finding accurate characteristic
C                 values using an iteration method
C             (2) CV0 for finding initial characteristic
C                 values using polynomial approximation
C             (3) CVQM for computing initial characteristic
C                 values for q ó 3*m
C             (3) CVQL for computing initial characteristic
C                 values for q ò m*m
C       ======================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IF (M.LE.12.OR.Q.LE.3.0*M.OR.Q.GT.M*M) THEN
	    CALL CV0(KD,M,Q,A)
	    IF (Q.NE.0.0D0) CALL REFINE(KD,M,Q,A,1)
	ELSE
	   NDIV=10
	   DELTA=(M-3.0)*M/NDIV
	   IF ((Q-3.0*M).LE.(M*M-Q)) THEN
5             NN=INT((Q-3.0*M)/DELTA)+1
	      DELTA=(Q-3.0*M)/NN
	      Q1=2.0*M
	      CALL CVQM(M,Q1,A1)
	      Q2=3.0*M
	      CALL CVQM(M,Q2,A2)
	      QQ=3.0*M
	      DO 10 I=1,NN
		 QQ=QQ+DELTA
		 A=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1)
		 IFLAG=1
		 IF (I.EQ.NN) IFLAG=-1
		 CALL REFINE(KD,M,QQ,A,IFLAG)
		 Q1=Q2
		 Q2=QQ
		 A1=A2
		 A2=A
10            CONTINUE
	      IF (IFLAG.EQ.-10) THEN
		 NDIV=NDIV*2
		 DELTA=(M-3.0)*M/NDIV
		 GO TO 5
	      ENDIF
	   ELSE
15            NN=INT((M*M-Q)/DELTA)+1
	      DELTA=(M*M-Q)/NN
	      Q1=M*(M-1.0)
	      CALL CVQL(KD,M,Q1,A1)
	      Q2=M*M
	      CALL CVQL(KD,M,Q2,A2)
	      QQ=M*M
	      DO 20 I=1,NN
		 QQ=QQ-DELTA
		 A=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1)
		 IFLAG=1
		 IF (I.EQ.NN) IFLAG=-1
		 CALL REFINE(KD,M,QQ,A,IFLAG)
		 Q1=Q2
		 Q2=QQ
		 A1=A2
		 A2=A
20            CONTINUE
	      IF (IFLAG.EQ.-10) THEN
		 NDIV=NDIV*2
		 DELTA=(M-3.0)*M/NDIV
		 GO TO 15
	      ENDIF
	   ENDIF
	ENDIF
	RETURN
	END


	SUBROUTINE REFINE(KD,M,Q,A,IFLAG)

C       =====================================================
C       Purpose: calculate the accurate characteristic value
C                by the secant method
C       Input :  m --- Order of Mathieu functions
C                q --- Parameter of Mathieu functions
C                A --- Initial characteristic value
C       Output:  A --- Refineed characteristic value
C       Routine called:  CVF for computing the value of F for
C                        characteristic equation
C       ========================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	EPS=1.0D-14
	MJ=10+M
	CA=A
	DELTA=0.0D0
	X0=A
	CALL CVF(KD,M,Q,X0,MJ,F0)
	X1=1.002*A
	CALL CVF(KD,M,Q,X1,MJ,F1)
5       DO 10 IT=1,100
	   MJ=MJ+1
	   X=X1-(X1-X0)/(1.0D0-F0/F1)
	   CALL CVF(KD,M,Q,X,MJ,F)
	   IF (ABS(1.0-X1/X).LT.EPS.OR.F.EQ.0.0) GO TO 15
	   X0=X1
	   F0=F1
	   X1=X
10         F1=F
15      A=X
	IF (DELTA.GT.0.05) THEN
	   A=CA
	   IF (IFLAG.LT.0) THEN
	      IFLAG=-10
	   ENDIF
	   RETURN
	ENDIF
	IF (ABS((A-CA)/CA).GT.0.05) THEN
	   X0=CA
	   DELTA=DELTA+0.005D0
	   CALL CVF(KD,M,Q,X0,MJ,F0)
	   X1=(1.0D0+DELTA)*CA
	   CALL CVF(KD,M,Q,X1,MJ,F1)
	   GO TO 5
	ENDIF
	RETURN
	END


	SUBROUTINE CVF(KD,M,Q,A,MJ,F)

C       ======================================================
C       Purpose: Compute the value of F for characteristic
C                equation of Mathieu functions
C       Input :  m --- Order of Mathieu functions
C                q --- Parameter of Mathieu functions
C                A --- Characteristic value
C       Output:  F --- Value of F for characteristic equation
C       ======================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	B=A
	IC=INT(M/2)
	L=0
	L0=0
	J0=2
	JF=IC
	IF (KD.EQ.1) L0=2
	IF (KD.EQ.1) J0=3
	IF (KD.EQ.2.OR.KD.EQ.3) L=1
	IF (KD.EQ.4) JF=IC-1
	T1=0.0D0
	DO 10 J=MJ,IC+1,-1
10         T1=-Q*Q/((2.0D0*J+L)**2-B+T1)
	IF (M.LE.2) THEN
	   T2=0.0D0
	   IF (KD.EQ.1.AND.M.EQ.0) T1=T1+T1
	   IF (KD.EQ.1.AND.M.EQ.2) T1=-2.0*Q*Q/(4.0-B+T1)-4.0
	   IF (KD.EQ.2.AND.M.EQ.1) T1=T1+Q
	   IF (KD.EQ.3.AND.M.EQ.1) T1=T1-Q
	ELSE
	   IF (KD.EQ.1) T0=4.0D0-B+2.0D0*Q*Q/B
	   IF (KD.EQ.2) T0=1.0D0-B+Q
	   IF (KD.EQ.3) T0=1.0D0-B-Q
	   IF (KD.EQ.4) T0=4.0D0-B
	   T2=-Q*Q/T0
	   DO 15 J=J0,JF
15            T2=-Q*Q/((2.0D0*J-L-L0)**2-B+T2)
	ENDIF
	F=(2.0D0*IC+L)**2+T1+T2-B
	RETURN
	END


	SUBROUTINE CV0(KD,M,Q,A0)

C       =====================================================
C       Purpose: Compute the initial characteristic value of
C                Mathieu functions for m ó 12  or q ó 300 or
C                q ò m*m
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C       Output:  A0 --- Characteristic value
C       Routines called:
C             (1) CVQM for computing initial characteristic
C                 value for q ó 3*m
C             (2) CVQL for computing initial characteristic
C                 value for q ò m*m
C       ====================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	Q2=Q*Q
	IF (M.EQ.0) THEN
	   IF (Q.LE.1.0) THEN
	      A0=(((.0036392*Q2-.0125868)*Q2+.0546875)*Q2-.5)*Q2
	   ELSE IF (Q.LE.10.0) THEN
	      A0=((3.999267D-3*Q-9.638957D-2)*Q-.88297)*Q
     &           +.5542818
	   ELSE
	      CALL CVQL(KD,M,Q,A0)
	   ENDIF
	ELSE IF (M.EQ.1) THEN
	   IF (Q.LE.1.0.AND.KD.EQ.2) THEN
	      A0=(((-6.51E-4*Q-.015625)*Q-.125)*Q+1.0)*Q+1.0
	   ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
	      A0=(((-6.51E-4*Q+.015625)*Q-.125)*Q-1.0)*Q+1.0
	   ELSE IF (Q.LE.10.0.AND. KD.EQ.2) THEN
	      A0=(((-4.94603D-4*Q+1.92917D-2)*Q-.3089229)
     &           *Q+1.33372)*Q+.811752
	   ELSE IF (Q.LE.10.0.AND.KD.EQ.3) THEN
	      A0=((1.971096D-3*Q-5.482465D-2)*Q-1.152218)
     &           *Q+1.10427
	   ELSE
	      CALL CVQL(KD,M,Q,A0)
	   ENDIF
	ELSE IF (M.EQ.2) THEN
	   IF (Q.LE.1.0.AND.KD.EQ.1) THEN
	      A0=(((-.0036391*Q2+.0125888)*Q2-.0551939)*Q2
     &           +.416667)*Q2+4.0
	   ELSE IF (Q.LE.1.0.AND.KD.EQ.4) THEN
	      A0=(.0003617*Q2-.0833333)*Q2+4.0
	   ELSE IF (Q.LE.15.AND.KD.EQ.1) THEN
	      A0=(((3.200972D-4*Q-8.667445D-3)*Q
     &           -1.829032D-4)*Q+.9919999)*Q+3.3290504
	   ELSE IF (Q.LE.10.0.AND.KD.EQ.4) THEN
	      A0=((2.38446D-3*Q-.08725329)*Q-4.732542D-3)
     &           *Q+4.00909
	   ELSE
	      CALL CVQL(KD,M,Q,A0)
	   ENDIF
	ELSE IF (M.EQ.3) THEN
	   IF (Q.LE.1.0.AND.KD.EQ.2) THEN
	      A0=((6.348E-4*Q+.015625)*Q+.0625)*Q2+9.0
	   ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
	      A0=((6.348E-4*Q-.015625)*Q+.0625)*Q2+9.0
	   ELSE IF (Q.LE.20.0.AND.KD.EQ.2) THEN
	      A0=(((3.035731D-4*Q-1.453021D-2)*Q
     &           +.19069602)*Q-.1039356)*Q+8.9449274
	   ELSE IF (Q.LE.15.0.AND.KD.EQ.3) THEN
	      A0=((9.369364D-5*Q-.03569325)*Q+.2689874)*Q
     &           +8.771735
	   ELSE
	      CALL CVQL(KD,M,Q,A0)
	   ENDIF
	ELSE IF (M.EQ.4) THEN
	   IF (Q.LE.1.0.AND.KD.EQ.1) THEN
	      A0=((-2.1E-6*Q2+5.012E-4)*Q2+.0333333)*Q2+16.0
	   ELSE IF (Q.LE.1.0.AND.KD.EQ.4) THEN
	      A0=((3.7E-6*Q2-3.669E-4)*Q2+.0333333)*Q2+16.0
	   ELSE IF (Q.LE.25.0.AND.KD.EQ.1) THEN
	      A0=(((1.076676D-4*Q-7.9684875D-3)*Q
     &           +.17344854)*Q-.5924058)*Q+16.620847
	   ELSE IF (Q.LE.20.0.AND.KD.EQ.4) THEN
	      A0=((-7.08719D-4*Q+3.8216144D-3)*Q
     &           +.1907493)*Q+15.744
	   ELSE
	      CALL CVQL(KD,M,Q,A0)
	   ENDIF
	ELSE IF (M.EQ.5) THEN
	   IF (Q.LE.1.0.AND.KD.EQ.2) THEN
	      A0=((6.8E-6*Q+1.42E-5)*Q2+.0208333)*Q2+25.0
	   ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
	      A0=((-6.8E-6*Q+1.42E-5)*Q2+.0208333)*Q2+25.0
	   ELSE IF (Q.LE.35.0.AND.KD.EQ.2) THEN
	      A0=(((2.238231D-5*Q-2.983416D-3)*Q
     &           +.10706975)*Q-.600205)*Q+25.93515
	   ELSE IF (Q.LE.25.0.AND.KD.EQ.3) THEN
	      A0=((-7.425364D-4*Q+2.18225D-2)*Q
     &           +4.16399D-2)*Q+24.897
	   ELSE
	      CALL CVQL(KD,M,Q,A0)
	   ENDIF
	ELSE IF (M.EQ.6) THEN
	   IF (Q.LE.1.0) THEN
	      A0=(.4D-6*Q2+.0142857)*Q2+36.0
	   ELSE IF (Q.LE.40.0.AND.KD.EQ.1) THEN
	      A0=(((-1.66846D-5*Q+4.80263D-4)*Q
     &           +2.53998D-2)*Q-.181233)*Q+36.423
	   ELSE IF (Q.LE.35.0.AND.KD.EQ.4) THEN
	      A0=((-4.57146D-4*Q+2.16609D-2)*Q-2.349616D-2)*Q
     &           +35.99251
	   ELSE
	      CALL CVQL(KD,M,Q,A0)
	   ENDIF
	ELSE IF (M.EQ.7) THEN
	   IF (Q.LE.10.0) THEN
	      CALL CVQM(M,Q,A0)
	   ELSE IF (Q.LE.50.0.AND.KD.EQ.2) THEN
	      A0=(((-1.411114D-5*Q+9.730514D-4)*Q
     &           -3.097887D-3)*Q+3.533597D-2)*Q+49.0547
	   ELSE IF (Q.LE.40.0.AND.KD.EQ.3) THEN
	      A0=((-3.043872D-4*Q+2.05511D-2)*Q
     &           -9.16292D-2)*Q+49.19035
	   ELSE
	      CALL CVQL(KD,M,Q,A0)
	   ENDIF
	ELSE IF (M.GE.8) THEN
	   IF (Q.LE.3.*M) THEN
	      CALL CVQM(M,Q,A0)
	   ELSE IF (Q.GT.M*M) THEN
	      CALL CVQL(KD,M,Q,A0)
	   ELSE
	      IF (M.EQ.8.AND.KD.EQ.1) THEN
		 A0=(((8.634308D-6*Q-2.100289D-3)*Q+.169072)*Q
     &              -4.64336)*Q+109.4211
	      ELSE IF (M.EQ.8.AND.KD.EQ.4) THEN
		 A0=((-6.7842D-5*Q+2.2057D-3)*Q+.48296)*Q+56.59
	      ELSE IF (M.EQ.9.AND.KD.EQ.2) THEN
		 A0=(((2.906435D-6*Q-1.019893D-3)*Q+.1101965)*Q
     &              -3.821851)*Q+127.6098
	      ELSE IF (M.EQ.9.AND.KD.EQ.3) THEN
		 A0=((-9.577289D-5*Q+.01043839)*Q+.06588934)*Q
     &              +78.0198
	      ELSE IF (M.EQ.10.AND.KD.EQ.1) THEN
		 A0=(((5.44927D-7*Q-3.926119D-4)*Q+.0612099)*Q
     &              -2.600805)*Q+138.1923
	      ELSE IF (M.EQ.10.AND.KD.EQ.4) THEN
		 A0=((-7.660143D-5*Q+.01132506)*Q-.09746023)*Q
     &              +99.29494
	      ELSE IF (M.EQ.11.AND.KD.EQ.2) THEN
		 A0=(((-5.67615D-7*Q+7.152722D-6)*Q+.01920291)*Q
     &              -1.081583)*Q+140.88
	      ELSE IF (M.EQ.11.AND.KD.EQ.3) THEN
		 A0=((-6.310551D-5*Q+.0119247)*Q-.2681195)*Q
     &              +123.667
	      ELSE IF (M.EQ.12.AND.KD.EQ.1) THEN
		 A0=(((-2.38351D-7*Q-2.90139D-5)*Q+.02023088)*Q
     &              -1.289)*Q+171.2723
	      ELSE IF (M.EQ.12.AND.KD.EQ.4) THEN
		 A0=(((3.08902D-7*Q-1.577869D-4)*Q+.0247911)*Q
     &              -1.05454)*Q+161.471
	      ENDIF
	   ENDIF
	ENDIF
	RETURN
	END


	SUBROUTINE CVQL(KD,M,Q,A0)

C       ========================================================
C       Purpose: Compute the characteristic value of Mathieu
C                functions  for q ò 3m
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C       Output:  A0 --- Initial characteristic value
C       ========================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IF (KD.EQ.1.OR.KD.EQ.2) W=2.0D0*M+1.0D0
	IF (KD.EQ.3.OR.KD.EQ.4) W=2.0D0*M-1.0D0
	W2=W*W
	W3=W*W2
	W4=W2*W2
	W6=W2*W4
	D1=5.0+34.0/W2+9.0/W4
	D2=(33.0+410.0/W2+405.0/W4)/W
	D3=(63.0+1260.0/W2+2943.0/W4+486.0/W6)/W2
	D4=(527.0+15617.0/W2+69001.0/W4+41607.0/W6)/W3
	C1=128.0
	P2=Q/W4
	P1=DSQRT(P2)
	CV1=-2.0*Q+2.0*W*DSQRT(Q)-(W2+1.0)/8.0
	CV2=(W+3.0/W)+D1/(32.0*P1)+D2/(8.0*C1*P2)
	CV2=CV2+D3/(64.0*C1*P1*P2)+D4/(16.0*C1*C1*P2*P2)
	A0=CV1-CV2/(C1*P1)
	RETURN
	END


	SUBROUTINE CVQM(M,Q,A0)

C       =====================================================
C       Purpose: Compute the characteristic value of Mathieu
C                functions for q ó m*m
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C       Output:  A0 --- Initial characteristic value
C       =====================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	HM1=.5*Q/(M*M-1.0)
	HM3=.25*HM1**3/(M*M-4.0)
	HM5=HM1*HM3*Q/((M*M-1.0)*(M*M-9.0))
	A0=M*M+Q*(HM1+(5.0*M*M+7.0)*HM3
     &     +(9.0*M**4+58.0*M*M+29.0)*HM5)
	RETURN
	END


	SUBROUTINE JYNB(N,X,NM,BJ,DJ,BY,DY)

C       =====================================================
C       Purpose: Compute Bessel functions Jn(x), Yn(x) and
C                their derivatives
C       Input :  x --- Argument of Jn(x) and Yn(x) ( x ò 0 )
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
