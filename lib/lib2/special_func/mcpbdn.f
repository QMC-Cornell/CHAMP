	PROGRAM MCPBDN

C       =============================================================
C       Purpose: This program computes parabolic cylinder functions
C                Dn(z) for an integer order and a complex argument
C                using subroutine CPBDN
C       Input :  x --- Real part of z
C                y --- Imaginary part of z
C                n --- Order of Dn(z)
C       Output:  CPB(|n|) --- Dn(z)
C                CPD(|n|) --- Dn'(z)
C       Example:
C                z = 5.0+ 5.0 i

C     n     Re[Dn(z)]      Im[Dn(z)]      Re[Dn'(z)]     Im[Dn'(z)]
C   -----------------------------------------------------------------
C     0   .99779828D+00  .66321897D-01 -.23286910D+01 -.26603004D+01
C     1   .46573819D+01  .53206009D+01  .26558457D+01 -.24878635D+02
C     2  -.43138931D+01  .49823592D+02  .14465848D+03 -.10313305D+03
C     3  -.28000219D+03  .21690729D+03  .12293320D+04  .30720802D+03
C     4  -.24716057D+04 -.46494526D+03  .38966424D+04  .82090067D+04
C    -1   .10813809D+00 -.90921592D-01 -.50014908D+00 -.23280660D-01
C    -2   .24998820D-02 -.19760577D-01 -.52486940D-01  .47769856D-01
C    -3  -.15821033D-02 -.23090595D-02 -.68249161D-03  .10032670D-01
C    -4  -.37829961D-03 -.10158757D-03  .89032322D-03  .11093416D-02
C       =============================================================

	IMPLICIT DOUBLE PRECISION (X,Y)
	IMPLICIT COMPLEX*16 (C,Z)
	DIMENSION CPB(0:100),CPD(0:100)
	WRITE(*,*)'Please enter n, x and y '
	READ(*,*)N,X,Y
	WRITE(*,20)N,X,Y
	Z=CMPLX(X,Y)
	N0=ABS(N)
	CALL CPBDN(N,Z,CPB,CPD)
	WRITE(*,*)
	IF (N.GE.0) THEN
	   WRITE(*,*)'  n     Re[Dn(z)]       Im[Dn(z)]       ',
     &               'Re[Dn''(z)]      Im[Dn''(z)]'
	ELSE
	   WRITE(*,*)' -n     Re[Dn(z)]       Im[Dn(z)]       ',
     &               'Re[Dn''(z)]      Im[Dn''(z)]'
	ENDIF
	WRITE(*,*)'-------------------------------------------',
     &            '-------------------------'
	DO 5 I=0,N0
5          WRITE(*,10)I,CPB(I),CPD(I)
10      FORMAT(1X,I3,4D16.8)
20      FORMAT(1X,'N =',I3,',   z =x+iy :',F6.2,'+',F6.2,' i')
	END


	SUBROUTINE CPBDN(N,Z,CPB,CPD)

C       ==================================================
C       Purpose: Compute the parabolic cylinder functions
C                 Dn(z) and Dn'(z) for a complex argument
C       Input:   z --- Complex argument of Dn(z)
C                n --- Order of Dn(z)  ( n=0,ס1,ס2,תתת )
C       Output:  CPB(|n|) --- Dn(z)
C                CPD(|n|) --- Dn'(z)
C       Routines called:
C            (1) CPDSA for computing Dn(z) for a small |z|
C            (2) CPDLA for computing Dn(z) for a large |z|
C       ==================================================

	IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Y)
	IMPLICIT COMPLEX*16 (C,Z)
	DIMENSION CPB(0:*),CPD(0:*)
	PI=3.141592653589793D0
	X=REAL(Z)
	A0=CDABS(Z)
	C0=(0.0D0,0.0D0)
	CA0=CDEXP(-0.25D0*Z*Z)
	IF (N.GE.0) THEN
	   CF0=CA0
	   CF1=Z*CA0
	   CPB(0)=CF0
	   CPB(1)=CF1
	   DO 10 K=2,N
	      CF=Z*CF1-(K-1.0D0)*CF0
	      CPB(K)=CF
	      CF0=CF1
10            CF1=CF
	ELSE
	   N0=-N
	   IF (X.LE.0.0.OR.CDABS(Z).EQ.0.0) THEN
	      CF0=CA0
	      CPB(0)=CF0
	      Z1=-Z
	      IF (A0.LE.7.0) THEN
		 CALL CPDSA(-1,Z1,CF1)
	      ELSE
		 CALL CPDLA(-1,Z1,CF1)
	      ENDIF
	      CF1=DSQRT(2.0D0*PI)/CA0-CF1
	      CPB(1)=CF1
	      DO 15 K=2,N0
		 CF=(-Z*CF1+CF0)/(K-1.0D0)
		 CPB(K)=CF
		 CF0=CF1
15               CF1=CF
	   ELSE
	      IF (A0.LE.3.0) THEN
		 CALL CPDSA(-N0,Z,CFA)
		 CPB(N0)=CFA
		 N1=N0+1
		 CALL CPDSA(-N1,Z,CFB)
		 CPB(N1)=CFB
		 NM1=N0-1
		 DO 20 K=NM1,0,-1
		    CF=Z*CFA+(K+1.0D0)*CFB
		    CPB(K)=CF
		    CFB=CFA
20                  CFA=CF
	      ELSE
		 M=100+ABS(N)
		 CFA=C0
		 CFB=(1.0D-30,0.0D0)
		 DO 25 K=M,0,-1
		    CF=Z*CFB+(K+1.0D0)*CFA
		    IF (K.LE.N0) CPB(K)=CF
		    CFA=CFB
25                  CFB=CF
		 CS0=CA0/CF
		 DO 30 K=0,N0
30                  CPB(K)=CS0*CPB(K)
	      ENDIF
	   ENDIF
	ENDIF
	CPD(0)=-0.5D0*Z*CPB(0)
	IF (N.GE.0) THEN
	   DO 35 K=1,N
35            CPD(K)=-0.5D0*Z*CPB(K)+K*CPB(K-1)
	ELSE
	   DO 40 K=1,N0
40            CPD(K)=0.5D0*Z*CPB(K)-CPB(K-1)
	ENDIF
	RETURN
	END


	SUBROUTINE CPDSA(N,Z,CDN)

C       ===========================================================
C       Purpose: Compute complex parabolic cylinder function Dn(z)
C                for small argument
C       Input:   z   --- complex argument of D(z)
C                n   --- Order of D(z) (n = 0,-1,-2,תתת)
C       Output:  CDN --- Dn(z)
C       Routine called: GAIH for computing ג(x), x=n/2 (n=1,2,...)
C       ===========================================================

	IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Y)
	IMPLICIT COMPLEX*16 (C,Z)
	EPS=1.0D-15
	PI=3.141592653589793D0
	SQ2=DSQRT(2.0D0)
	CA0=CDEXP(-.25D0*Z*Z)
	VA0=0.5D0*(1.0D0-N)
	IF (N.EQ.0.0) THEN
	   CDN=CA0
	ELSE
	   IF (CDABS(Z).EQ.0.0) THEN
	      IF (VA0.LE.0.0.AND.VA0.EQ.INT(VA0)) THEN
		 CDN=0.0D0
	      ELSE
		 CALL GAIH(VA0,GA0)
		 PD=DSQRT(PI)/(2.0D0**(-.5D0*N)*GA0)
		 CDN=CMPLX(PD,0.0D0)
	      ENDIF
	   ELSE
	      XN=-N
	      CALL GAIH(XN,G1)
	      CB0=2.0D0**(-0.5D0*N-1.0D0)*CA0/G1
	      VT=-.5D0*N
	      CALL GAIH(VT,G0)
	      CDN=CMPLX(G0,0.0D0)
	      CR=(1.0D0,0.0D0)
	      DO 10 M=1,250
		 VM=.5D0*(M-N)
		 CALL GAIH(VM,GM)
		 CR=-CR*SQ2*Z/M
		 CDW=GM*CR
		 CDN=CDN+CDW
		 IF (CDABS(CDW).LT.CDABS(CDN)*EPS) GO TO 20
10            CONTINUE
20            CDN=CB0*CDN
	   ENDIF
	ENDIF
	RETURN
	END


	SUBROUTINE CPDLA(N,Z,CDN)

C       ===========================================================
C       Purpose: Compute complex parabolic cylinder function Dn(z)
C                for large argument
C       Input:   z   --- Complex argument of Dn(z)
C                n   --- Order of Dn(z) (n = 0,ס1,ס2,תתת)
C       Output:  CDN --- Dn(z)
C       ===========================================================

	IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Y)
	IMPLICIT COMPLEX*16 (C,Z)
	CB0=Z**N*CDEXP(-.25D0*Z*Z)
	CR=(1.0D0,0.0D0)
	CDN=(1.0D0,0.0D0)
	DO 10 K=1,16
	   CR=-0.5D0*CR*(2.0*K-N-1.0)*(2.0*K-N-2.0)/(K*Z*Z)
	   CDN=CDN+CR
	   IF (CDABS(CR).LT.CDABS(CDN)*1.0D-12) GO TO 15
10      CONTINUE
15      CDN=CB0*CDN
	RETURN
	END


	SUBROUTINE GAIH(X,GA)

C       =====================================================
C       Purpose: Compute gamma function ג(x)
C       Input :  x  --- Argument of ג(x), x = n/2, n=1,2,תתת
C       Output:  GA --- ג(x)
C       =====================================================

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PI=3.141592653589793D0
	IF (X.EQ.INT(X).AND.X.GT.0.0) THEN
	   GA=1.0D0
	   M1=INT(X-1.0)
	   DO 10 K=2,M1
10            GA=GA*K
	ELSE IF (X+.5D0.EQ.INT(X+.5D0).AND.X.GT.0.0) THEN
	   M=INT(X)
	   GA=DSQRT(PI)
	   DO 15 K=1,M
15            GA=0.5D0*GA*(2.0D0*K-1.0D0)
	ENDIF
	RETURN
	END
