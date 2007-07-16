        PROGRAM MCIKNA

C       =============================================================
C       Purpose: This program computes the modified Bessel functions
C                In(z) and Kn(z), and their derivatives for a
C                complex argument using subroutine CIKNA
C       Input :  z --- Complex argument of In(z) and Kn(z)
C                n --- Order of In(z) and Kn(z)
C                      ( n = 0,1,תתת, n ף 250 )
C       Output:  CBI(n) --- In(z)
C                CDI(n) --- In'(z)
C                CBK(n) --- Kn(z)
C                CDK(n) --- Kn'(z)
C       Example: z = 4.0 + i 2.0 ,      Nmax = 5

C     n     Re[In(z)]      Im[In(z)]      Re[In'(z)]     Im[In'(z)]
C   -----------------------------------------------------------------
C     0  -.19056142D+01  .10403505D+02 -.23059657D+01  .92222463D+01
C     1  -.23059657D+01  .92222463D+01 -.23666457D+01  .83284588D+01
C     2  -.28276772D+01  .62534130D+01 -.24255774D+01  .61553456D+01
C     3  -.25451891D+01  .30884450D+01 -.22270972D+01  .36367893D+01
C     4  -.16265172D+01  .10201656D+01 -.16520416D+01  .16217056D+01
C     5  -.75889410D+00  .15496632D+00 -.94510625D+00  .48575220D+00

C     n     Re[Kn(z)]      Im[Kn(z)]      Re[Kn'(z)]     Im[Kn'(z)]
C   -----------------------------------------------------------------
C     0  -.64221754D-02 -.84393648D-02  .74307276D-02  .89585853D-02
C     1  -.74307276D-02 -.89585853D-02  .88041795D-02  .94880091D-02
C     2  -.11186184D-01 -.10536653D-01  .14012532D-01  .10936010D-01
C     3  -.20594336D-01 -.12913435D-01  .27416815D-01  .12106413D-01
C     4  -.43647447D-01 -.13676173D-01  .60982763D-01  .63953943D-02
C     5  -.10137119D+00  .12264588D-03  .14495731D+00 -.37132068D-01
C       =============================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        COMMON CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)
        WRITE(*,*)'  Please input n, x,y (z=x+iy)=?'
        READ(*,*)N,X,Y
        Z=CMPLX(X,Y)
        WRITE(*,40)X,Y,N
        IF (N.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)' Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL CIKNA(N,Z,NM,CBI,CDI,CBK,CDK)
        WRITE(*,*)
        WRITE(*,*)'   n      Re[In(z)]       Im[In(z)]',
     &            '       Re[In''(z)]      Im[In''(z)]'
        WRITE(*,*)' -----------------------------------',
     &            '----------------------------------'
        DO 10 K=0,NM,NS
           WRITE(*,30)K,CBI(K),CDI(K)
10      CONTINUE
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Kn(z)]       Im[Kn(z)]',
     &            '       Re[Kn''(z)]      Im[Kn''(z)]'
        WRITE(*,*)' -----------------------------------',
     &            '----------------------------------'
        DO 20 K=0,NM,NS
           WRITE(*,30)K,CBK(K),CDK(K)
20      CONTINUE
30      FORMAT(1X,I4,1X,4D16.8)
40      FORMAT(3X,3Hz =,F7.1,' + i',F7.1,' ,',6X,6HNmax =,I4)
        END


        SUBROUTINE CIKNA(N,Z,NM,CBI,CDI,CBK,CDK)

C       ========================================================
C       Purpose: Compute modified Bessel functions In(z), Kn(x)
C                and their derivatives for a complex argument
C       Input :  z --- Complex argument of In(z) and Kn(z)
C                n --- Order of In(z) and Kn(z)
C       Output:  CBI(n) --- In(z)
C                CDI(n) --- In'(z)
C                CBK(n) --- Kn(z)
C                CDK(n) --- Kn'(z)
C                NM --- Highest order computed
C       Routines called:
C             (1) CIK01 to compute I0(z), I1(z) K0(z) & K1(z)
C             (2) MSTA1 and MSTA2 to compute the starting
C                 point for backward recurrence
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A,B,P,W,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:N),CDI(0:N),CBK(0:N),CDK(0:N)
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBI(K)=(0.0D0,0.0D0)
              CDI(K)=(0.0D0,0.0D0)
              CBK(K)=-(1.0D+300,0.0D0)
10            CDK(K)=(1.0D+300,0.0D0)
           CBI(0)=(1.0D0,0.0D0)
           CDI(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        CALL CIK01(Z,CBI0,CDI0,CBI1,CDI1,CBK0,CDK0,CBK1,CDK1)
        CBI(0)=CBI0
        CBI(1)=CBI1
        CBK(0)=CBK0
        CBK(1)=CBK1
        CDI(0)=CDI0
        CDI(1)=CDI1
        CDK(0)=CDK0
        CDK(1)=CDK1
        IF (N.LE.1) RETURN
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           NM=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 45 K=M,0,-1
           CF=2.0D0*(K+1.0D0)/Z*CF1+CF2
           IF (K.LE.NM) CBI(K)=CF
           CF2=CF1
45         CF1=CF
        CS=CBI0/CF
        DO 50 K=0,NM
50         CBI(K)=CS*CBI(K)
        DO 60 K=2,NM
           IF (CDABS(CBI(K-1)).GT.CDABS(CBI(K-2))) THEN
              CKK=(1.0D0/Z-CBI(K)*CBK(K-1))/CBI(K-1)
           ELSE
              CKK=(CBI(K)*CBK(K-2)+2.0D0*(K-1.0D0)/(Z*Z))/CBI(K-2)
           ENDIF
60         CBK(K)=CKK
        DO 70 K=2,NM
           CDI(K)=CBI(K-1)-K/Z*CBI(K)
70         CDK(K)=-CBK(K-1)-K/Z*CBK(K)
        RETURN
        END


        SUBROUTINE CIK01(Z,CBI0,CDI0,CBI1,CDI1,CBK0,CDK0,CBK1,CDK1)

C       ==========================================================
C       Purpose: Compute modified complex Bessel functions I0(z),
C                I1(z), K0(z), K1(z), and their derivatives
C       Input :  z --- Complex argument
C       Output:  CBI0 --- I0(z)
C                CDI0 --- I0'(z)
C                CBI1 --- I1(z)
C                CDI1 --- I1'(z)
C                CBK0 --- K0(z)
C                CDK0 --- K0'(z)
C                CBK1 --- K1(z)
C                CDK1 --- K1'(z)
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION A(12),B(12),A1(10)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBI0=(1.0D0,0.0D0)
           CBI1=(0.0D0,0.0D0)
           CDI0=(0.0D0,0.0D0)
           CDI1=(0.5D0,0.0D0)
           CBK0=(1.0D+300,0.0D0)
           CBK1=(1.0D+300,0.0D0)
           CDK0=-(1.0D+300,0.0D0)
           CDK1=-(1.0D+300,0.0D0)
           RETURN
        ENDIF
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LE.18.0) THEN
           CBI0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,50
              CR=0.25D0*CR*Z2/(K*K)
              CBI0=CBI0+CR
              IF (CDABS(CR/CBI0).LT.1.0D-15) GO TO 15
10         CONTINUE
15         CBI1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,50
              CR=0.25D0*CR*Z2/(K*(K+1))
              CBI1=CBI1+CR
              IF (CDABS(CR/CBI1).LT.1.0D-15) GO TO 25
20         CONTINUE
25         CBI1=0.5D0*Z1*CBI1
        ELSE
           DATA A/0.125D0,7.03125D-2,
     &            7.32421875D-2,1.1215209960938D-1,
     &            2.2710800170898D-1,5.7250142097473D-1,
     &            1.7277275025845D0,6.0740420012735D0,
     &            2.4380529699556D01,1.1001714026925D02,
     &            5.5133589612202D02,3.0380905109224D03/
           DATA B/-0.375D0,-1.171875D-1,
     &            -1.025390625D-1,-1.4419555664063D-1,
     &            -2.7757644653320D-1,-6.7659258842468D-1,
     &            -1.9935317337513D0,-6.8839142681099D0,
     &            -2.7248827311269D01,-1.2159789187654D02,
     &            -6.0384407670507D02,-3.3022722944809D03/
           K0=12
           IF (A0.GE.35.0) K0=9
           IF (A0.GE.50.0) K0=7
           CA=CDEXP(Z1)/CDSQRT(2.0D0*PI*Z1)
           CBI0=(1.0D0,0.0D0)
           ZR=1.0D0/Z1
           DO 30 K=1,K0
30            CBI0=CBI0+A(K)*ZR**K
           CBI0=CA*CBI0
           CBI1=(1.0D0,0.0D0)
           DO 35 K=1,K0
35            CBI1=CBI1+B(K)*ZR**K
           CBI1=CA*CBI1
        ENDIF
        IF (A0.LE.9.0) THEN
           CS=(0.0D0,0.0D0)
           CT=-CDLOG(0.5D0*Z1)-0.5772156649015329D0
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           DO 40 K=1,50
              W0=W0+1.0D0/K
              CR=0.25D0*CR/(K*K)*Z2
              CS=CS+CR*(W0+CT)
              IF (CDABS((CS-CW)/CS).LT.1.0D-15) GO TO 45
40            CW=CS
45         CBK0=CT+CS
        ELSE
           DATA A1/0.125D0,0.2109375D0,
     &             1.0986328125D0,1.1775970458984D01,
     &             2.1461706161499D02,5.9511522710323D03,
     &             2.3347645606175D05,1.2312234987631D07,
     &             8.401390346421D08,7.2031420482627D10/
           CB=0.5D0/Z1
           ZR2=1.0D0/Z2
           CBK0=(1.0D0,0.0D0)
           DO 50 K=1,10
50            CBK0=CBK0+A1(K)*ZR2**K
           CBK0=CB*CBK0/CBI0
        ENDIF
        CBK1=(1.0D0/Z1-CBI1*CBK0)/CBI0
        IF (REAL(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBK0=CBK0+CI*PI*CBI0
           IF (DIMAG(Z).GT.0.0) CBK0=CBK0-CI*PI*CBI0
           IF (DIMAG(Z).LT.0.0) CBK1=-CBK1+CI*PI*CBI1
           IF (DIMAG(Z).GT.0.0) CBK1=-CBK1-CI*PI*CBI1
           CBI1=-CBI1
        ENDIF
        CDI0=CBI1
        CDI1=CBI0-1.0D0/Z*CBI1
        CDK0=-CBK1
        CDK1=-CBK0-1.0D0/Z*CBK1
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
