	PROGRAM MSCKB

C       ============================================================
C       Purpose: This program computes the expansion coefficients
C                of the prolate and oblate spheroidal functions,
C                c2k, using subroutine SCKB
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                cv --- Characteristic value
C                KD --- Function code
C                       KD=1 for prolate; KD=-1 for oblate
C       Output:  CK(k) --- Expansion coefficients ck;
C                          CK(1), CK(2), ... correspond to
C                          c0, c2, ...
C       Example: Compute the first 13 expansion coefficients C2k for
C                KD= 1, m=2, n=3, c=3.0 and cv=14.8277782138; and
C                KD=-1, m=2, n=3, c=3.0 and cv=8.80939392077

C                Coefficients of Prolate and oblate functions

C                  k         C2k(c)              C2k(-ic)
C                ---------------------------------------------
C                  0     .9173213327D+01     .2489664942D+02
C                  1     .4718258929D+01    -.1205287032D+02
C                  2     .9841212916D+00     .2410564082D+01
C                  3     .1151870224D+00    -.2735821590D+00
C                  4     .8733916403D-02     .2026057157D-01
C                  5     .4663888254D-03    -.1061946315D-02
C                  6     .1853910398D-04     .4158091152D-04
C                  7     .5708084895D-06    -.1264400411D-05
C                  8     .1402786472D-07     .3074963448D-07
C                  9     .2817194508D-09    -.6120579463D-09
C                 10     .4712094447D-11     .1015900041D-10
C                 11     .6667838485D-13    -.1427953361D-12
C                 12     .8087995432D-15     .1721924955D-14
C       ============================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION CK(200),DF(200),EG(200)
	WRITE(*,*)'Please KD, m, n and c '
	READ(*,*)KD,M,N,C
	CALL SEGV(M,N,C,KD,CV,EG)
	WRITE(*,30)KD,M,N,C,CV
	CALL SDMN(M,N,C,CV,KD,DF)
	CALL SCKB(M,N,C,DF,CK)
	WRITE(*,*)
	IF (KD.EQ.1) THEN
	   WRITE(*,*)'Coefficients of Prolate function'
	   WRITE(*,*)
	   WRITE(*,*)'   k            C2k(c)'
	ELSE
	   WRITE(*,*)'Coefficients of Oblate function'
	   WRITE(*,*)
	   WRITE(*,*)'   k           C2k(-ic)'
	ENDIF
	WRITE(*,*)'----------------------------'
	NM=25+INT((N-M)/2+C)
	DO 10 K=1,NM
10         WRITE(*,20)K-1,CK(K)
20      FORMAT(2X,I3,4X,D18.10)
30      FORMAT(1X,3HKD=,I3,',  ',2Hm=,I3,',  ',2Hn=,I3,
     &         ',  ',2Hc=,F5.1,',  ',4Hcv =,F18.10)
	END


	SUBROUTINE SCKB(M,N,C,DF,CK)

C       ======================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions, c2k
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
	FAC=-0.5D0**M
	REG=1.0D0
	IF (M+NM.GT.80) REG=1.0D-200
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
	      IF (DABS(SW-SUM).LT.DABS(SUM)*1.0D-14) GO TO 25
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
C                prolate and oblate spheroidal functions
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
