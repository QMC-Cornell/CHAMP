	PROGRAM MASWFA

C       ============================================================
C       Purpose: This program computes the prolate and oblate
C                spheroidal angular functions of the first
C                kind and their derivatives using subroutine ASWFA
C       Input :  m  --- Mode parameter,  m = 0,1,2,...
C                n  --- Mode parameter,  n = m,m+1,...
C                c  --- Spheroidal parameter
C                x  --- Argument of angular function, |x| < 1.0
C                KD --- Function code
C                       KD=1 for prolate;  KD=-1 for oblate
C                cv --- Characteristic value
C       Output:  S1F --- Angular function of the first kind
C                S1D --- Derivative of the angular function of
C                        the first kind
C       Examples:
C               KD = 1, m = 2, n = 3, c = 3.0 and cv = 14.8277782138
C                  x         Smn(c,x)            Smn'(c,x)
C                --------------------------------------------
C                 0.2      .28261309D+01       .12418631D+02
C                 0.5      .49938554D+01       .92761604D+00
C                 0.8      .31693975D+01      -.12646552D+02

C               KD =-1, m = 2, n = 3, c = 3.0 and cv = 8.8093939208
C                  x         Smn(-ic,x)         Smn'(-ic,x)
C                --------------------------------------------
C                 0.2      .29417848D+01       .14106305D+02
C                 0.5      .64138827D+01       .76007194D+01
C                 0.8      .60069873D+01      -.14387479D+02
C       ============================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION EG(200)
	WRITE(*,*)'Please enter KD, m, n and c '
	READ(*,*)KD,M,N,C
	WRITE(*,10)KD,M,N,C
	CALL SEGV(M,N,C,KD,CV,EG)
	WRITE(*,20)CV
	WRITE(*,*)
	IF (KD.EQ.1 ) THEN
	   WRITE(*,*)'    x         Smn(c,x)            Smn''(c,x)'
	ELSE IF (KD.EQ.-1) THEN
	   WRITE(*,*)'    x         Smn(-ic,x)         Smn''(-ic,x)'
	ENDIF
	WRITE(*,*)'  --------------------------------------------'
	DO 5 I=0,20
	   X=-1.0D0+0.1D0*I
	   CALL ASWFA(M,N,C,X,KD,CV,S1F,S1D)
5          WRITE(*,30)X,S1F,S1D
10      FORMAT(1X,'KD ='I2,', ','m =',I2,', ','n =',I2,', ','c =',F5.1)
20      FORMAT(1X,' cv =',F18.10)
30      FORMAT(1X,F5.1,2D20.8)
	END


	SUBROUTINE ASWFA(M,N,C,X,KD,CV,S1F,S1D)

C       ===========================================================
C       Purpose: Compute the prolate and oblate spheroidal angular
C                functions of the first kind and their derivatives
C       Input :  m  --- Mode parameter,  m = 0,1,2,...
C                n  --- Mode parameter,  n = m,m+1,...
C                c  --- Spheroidal parameter
C                x  --- Argument of angular function, |x| < 1.0
C                KD --- Function code
C                       KD=1 for prolate;  KD=-1 for oblate
C                cv --- Characteristic value
C       Output:  S1F --- Angular function of the first kind
C                S1D --- Derivative of the angular function of
C                        the first kind
C       Routine called:
C                SCKB for computing expansion coefficients ck
C       ===========================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION CK(200),DF(200)
	EPS=1.0D-14
	X0=X
	X=DABS(X)
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	NM=10+INT((N-M)/2+C)
	NM2=NM/2-2
	CALL SDMN(M,N,C,CV,KD,DF)
	CALL SCKB(M,N,C,DF,CK)
	X1=1.0D0-X*X
	IF (M.EQ.0.AND.X1.EQ.0.0D0) THEN
	   A0=1.0D0
	ELSE
	   A0=X1**(0.5D0*M)
	ENDIF
	SU1=CK(1)
	DO 10 K=1,NM2
	   R=CK(K+1)*X1**K
	   SU1=SU1+R
	   IF (K.GE.10.AND.DABS(R/SU1).LT.EPS) GO TO 15
10         CONTINUE
15      S1F=A0*X**IP*SU1
	IF (X.EQ.1.0D0) THEN
	   IF (M.EQ.0) S1D=IP*CK(1)-2.0D0*CK(2)
	   IF (M.EQ.1) S1D=-1.0D+100
	   IF (M.EQ.2) S1D=-2.0D0*CK(1)
	   IF (M.GE.3) S1D=0.0D0
	ELSE
	   D0=IP-M/X1*X**(IP+1.0D0)
	   D1=-2.0D0*A0*X**(IP+1.0D0)
	   SU2=CK(2)
	   DO 20 K=2,NM2
	      R=K*CK(K+1)*X1**(K-1.0D0)
	      SU2=SU2+R
	      IF (K.GE.10.AND.DABS(R/SU2).LT.EPS) GO TO 25
20            CONTINUE
25         S1D=D0*A0*SU1+D1*SU2
	ENDIF
	IF (X0.LT.0.0D0.AND.IP.EQ.0) S1D=-S1D
	IF (X0.LT.0.0D0.AND.IP.EQ.1) S1F=-S1F
	X=X0
	RETURN
	END


	SUBROUTINE SCKB(M,N,C,DF,CK)

C       ======================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                DF(k) --- Expansion coefficients dk
C       Output:  CK(k) --- Expansion coefficients ck;
C                          CK(1), CK(2), ... correspond to
C                          c0, c2, ...
C       ======================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION DF(200),CK(200)
	IF (C.LE.1.0D-10) C=1.0D-10
	NM=25+INT(0.5*(N-M)+C)
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	REG=1.0D0
	IF (M+NM.GT.80) REG=1.0D-200
	FAC=-0.5D0**M
	DO 35 K=0,NM-1
	   FAC=-FAC
	   I1=2*K+IP+1
	   R=REG
	   DO 10 I=I1,I1+2*M-1
10            R=R*I
	   I2=K+M+IP
	   DO 15 I=I2,I2+K-1
15            R=R*(I+0.5D0)
	   SUM=R*DF(K+1)
	   DO 20 I=K+1,NM
	      D1=2.0D0*I+IP
	      D2=2.0D0*M+D1
	      D3=I+M+IP-0.5D0
	      R=R*D2*(D2-1.0D0)*I*(D3+K)/(D1*(D1-1.0D0)*(I-K)*D3)
	      SUM=SUM+R*DF(I+1)
	      IF (DABS(SW-SUM).LT.DABS(SUM)*1.0D-14) GOTO 25
20            SW=SUM
25         R1=REG
	   DO 30 I=2,M+K
30            R1=R1*I
35         CK(K+1)=FAC*SUM/R1
	RETURN
	END


	SUBROUTINE SDMN(M,N,C,CV,KD,DF)

C       =====================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions, dk
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                cv --- Characteristic value
C                KD --- Function code
C                       KD=1 for prolate; KD=-1 for oblate
C       Output:  DF(k) --- Expansion coefficients dk;
C                          DF(1), DF(2), ... correspond to
C                          d0, d2, ... for even n-m and d1,
C                          d3, ... for odd n-m
C       =====================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION A(200),D(200),G(200),DF(200)
	NM=25+INT(0.5*(N-M)+C)
	IF (C.LT.1.0D-10) THEN
	   DO 5 I=1,NM
5             DF(I)=0D0
	   DF((N-M)/2+1)=1.0D0
	   RETURN
	ENDIF
	CS=C*C*KD
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	DO 10 I=1,NM+2
	   IF (IP.EQ.0) K=2*(I-1)
	   IF (IP.EQ.1) K=2*I-1
	   DK0=M+K
	   DK1=M+K+1
	   DK2=2*(M+K)
	   D2K=2*M+K
	   A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
	   D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &          *(DK2+3.0))*CS
	   G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
10      CONTINUE
	FS=1.0D0
	F1=0.0D0
	F0=1.0D-100
	KB=0
	DF(NM+1)=0.0D0
	DO 30 K=NM,1,-1
	   F=-((D(K+1)-CV)*F0+A(K+1)*F1)/G(K+1)
	   IF (DABS(F).GT.DABS(DF(K+1))) THEN
	      DF(K)=F
	      F1=F0
	      F0=F
	      IF (DABS(F).GT.1.0D+100) THEN
		 DO 12 K1=K,NM
12                  DF(K1)=DF(K1)*1.0D-100
		 F1=F1*1.0D-100
		 F0=F0*1.0D-100
	      ENDIF
	   ELSE
	      KB=K
	      FL=DF(K+1)
	      F1=1.0D-100
	      F2=-(D(1)-CV)/A(1)*F1
	      DF(1)=F1
	      IF (KB.EQ.1) THEN
		 FS=F2
	      ELSE IF (KB.EQ.2) THEN
		 DF(2)=F2
		 FS=-((D(2)-CV)*F2+G(2)*F1)/A(2)
	      ELSE
		 DF(2)=F2
		 DO 20 J=3,KB+1
		    F=-((D(J-1)-CV)*F2+G(J-1)*F1)/A(J-1)
		    IF (J.LE.KB) DF(J)=F
		    IF (DABS(F).GT.1.0D+100) THEN
		       DO 15 K1=1,J
15                        DF(K1)=DF(K1)*1.0D-100
		       F=F*1.0D-100
		       F2=F2*1.0D-100
		    ENDIF
		    F1=F2
20                  F2=F
		 FS=F
	      ENDIF
	      GO TO 35
	   ENDIF
30      CONTINUE
35      SU1=0.0D0
	R1=1.0D0
	DO 40 J=M+IP+1,2*(M+IP)
40         R1=R1*J
	SU1=DF(1)*R1
	DO 45 K=2,KB
	   R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
45           SU1=SU1+R1*DF(K)
	SU2=0.0D0
	DO 50 K=KB+1,NM
	   IF (K.NE.1) R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
	   SU2=SU2+R1*DF(K)
	   IF (DABS(SW-SU2).LT.DABS(SU2)*1.0D-14) GOTO 55
50         SW=SU2
55      R3=1.0D0
	DO 60 J=1,(M+N+IP)/2
60         R3=R3*(J+0.5D0*(N+M+IP))
	R4=1.0D0
	DO 65 J=1,(N-M-IP)/2
65         R4=-4.0D0*R4*J
	S0=R3/(FL*(SU1/FS)+SU2)/R4
	DO 70 K=1,KB
70         DF(K)=FL/FS*S0*DF(K)
	DO 75 K=KB+1,NM
75         DF(K)=S0*DF(K)
	RETURN
	END



	SUBROUTINE SEGV(M,N,C,KD,CV,EG)

C       =========================================================
C       Purpose: Compute the characteristic values of spheroidal
C                wave functions
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                KD --- Function code
C                       KD=1 for Prolate; KD=-1 for Oblate
C       Output:  CV --- Characteristic value for given m, n and c
C                EG(L) --- Characteristic value for mode m and n'
C                          ( L = n' - m + 1 )
C       =========================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION B(100),H(100),D(300),E(300),F(300),CV0(100),
     &            A(300),G(300),EG(200)
	IF (C.LT.1.0D-10) THEN
	   DO 5 I=1,N
5             EG(I)=(I+M)*(I+M-1.0D0)
	   GO TO 70
	ENDIF
	ICM=(N-M+2)/2
	NM=10+INT(0.5*(N-M)+C)
	CS=C*C*KD
	DO 60 L=0,1
	   DO 10 I=1,NM
	      IF (L.EQ.0) K=2*(I-1)
	      IF (L.EQ.1) K=2*I-1
	      DK0=M+K
	      DK1=M+K+1
	      DK2=2*(M+K)
	      D2K=2*M+K
	      A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
	      D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &             *(DK2+3.0))*CS
10            G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
	   DO 15 K=2,NM
	      E(K)=DSQRT(A(K-1)*G(K))
15            F(K)=E(K)*E(K)
	   F(1)=0.0D0
	   E(1)=0.0D0
	   XA=D(NM)+DABS(E(NM))
	   XB=D(NM)-DABS(E(NM))
	   NM1=NM-1
	   DO 20 I=1,NM1
	      T=DABS(E(I))+DABS(E(I+1))
	      T1=D(I)+T
	      IF (XA.LT.T1) XA=T1
	      T1=D(I)-T
	      IF (T1.LT.XB) XB=T1
20         CONTINUE
	   DO 25 I=1,ICM
	      B(I)=XA
25            H(I)=XB
	   DO 55 K=1,ICM
	      DO 30 K1=K,ICM
		 IF (B(K1).LT.B(K)) THEN
		    B(K)=B(K1)
		    GO TO 35
		 ENDIF
30            CONTINUE
35            IF (K.NE.1.AND.H(K).LT.H(K-1)) H(K)=H(K-1)
40            X1=(B(K)+H(K))/2.0D0
	      CV0(K)=X1
	      IF (DABS((B(K)-H(K))/X1).LT.1.0D-14) GO TO 50
	      J=0
	      S=1.0D0
	      DO 45 I=1,NM
		 IF (S.EQ.0.0D0) S=S+1.0D-30
		 T=F(I)/S
		 S=D(I)-T-X1
		 IF (S.LT.0.0D0) J=J+1
45            CONTINUE
	      IF (J.LT.K) THEN
		 H(K)=X1
	      ELSE
		 B(K)=X1
		 IF (J.GE.ICM) THEN
		    B(ICM)=X1
		 ELSE
		    IF (H(J+1).LT.X1) H(J+1)=X1
		    IF (X1.LT.B(J)) B(J)=X1
		 ENDIF
	      ENDIF
	      GO TO 40
50            CV0(K)=X1
	      IF (L.EQ.0) EG(2*K-1)=CV0(K)
	      IF (L.EQ.1) EG(2*K)=CV0(K)
55         CONTINUE
60      CONTINUE
70      CV=EG(N-M+1)
	RETURN
	END

