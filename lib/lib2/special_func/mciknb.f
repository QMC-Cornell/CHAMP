        PROGRAM MCIKNB

C       =============================================================
C       Purpose: This program computes the modified Bessel functions
C                In(z) and Kn(z), and their derivatives for a
C                complex argument using subroutine CIKNB
C       Input:   z --- Complex argument
C                n --- Order of In(z) and Kn(z)
C                      ( n = 0,1,תתת, n ף 250 )
C       Output:  CBI(n) --- In(z)
C                CDI(n) --- In'(z)
C                CBK(n) --- Kn(z)
C                CDK(n) --- Kn'(z)
C       Example: Nmax = 5,   z = 4.0 + i 2.0

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
        DIMENSION CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)
        WRITE(*,*)'  Please enter n, x and y (z = x+iy) '
        READ(*,*)N,X,Y
        WRITE(*,25)N,X,Y
        Z=CMPLX(X,Y)
        WRITE(*,*)
        IF (N.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL CIKNB(N,Z,NM,CBI,CDI,CBK,CDK)
        WRITE(*,*)'   n      Re[In(z)]       Im[In(z)]',
     &            '      Re[In''(z)]       Im[In''(z)]'
        WRITE(*,*)' ---------------------------------',
     &            '------------------------------------'
        DO 10 K=0,NM,NS
           WRITE(*,20)K,CBI(K),CDI(K)
10      CONTINUE
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Kn(z)]       Im[Kn(z)]',
     &            '      Re[Kn''(z)]       Im[Kn''(z)]'
        WRITE(*,*)' ---------------------------------',
     &            '------------------------------------'
        DO 15 K=0,NM,NS
           WRITE(*,20)K,CBK(K),CDK(K)
15      CONTINUE
20      FORMAT(1X,1X,I3,1X,4D16.8)
25      FORMAT(3X,6HNmaz =,I3,',    ','z =', F6.1,' + i',F6.1)
        END


        SUBROUTINE CIKNB(N,Z,NM,CBI,CDI,CBK,CDK)

C       ============================================================
C       Purpose: Compute modified Bessel functions In(z) and Kn(z),
C                and their derivatives for a complex argument
C       Input:   z --- Complex argument
C                n --- Order of In(z) and Kn(z)
C       Output:  CBI(n) --- In(z)
C                CDI(n) --- In'(z)
C                CBK(n) --- Kn(z)
C                CDK(n) --- Kn'(z)
C                NM --- Highest order computed
C       Routones called:
C                MSTA1 and MSTA2 to compute the starting point for
C                backward recurrence
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:N),CDI(0:N),CBK(0:N),CDK(0:N)
        PI=3.141592653589793D0
        EL=0.57721566490153D0
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBI(K)=(0.0D0,0.0D0)
              CBK(K)=(1.0D+300,0.0D0)
              CDI(K)=(0.0D0,0.0D0)
10            CDK(K)=-(1.0D+300,0.0D0)
           CBI(0)=(1.0D0,0.0D0)
           CDI(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        Z1=Z
        CI=(0.0D0,1.0D0)
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (N.EQ.0) NM=1
        M=MSTA1(A0,200)
        IF (M.LT.NM) THEN
           NM=M
        ELSE
           M=MSTA2(A0,NM,15)
        ENDIF
        CBS=0.0D0
        CSK0=0.0D0
        CF0=0.0D0
        CF1=1.0D-100
        DO 15 K=M,0,-1
           CF=2.0D0*(K+1.0D0)*CF1/Z1+CF0
           IF (K.LE.NM) CBI(K)=CF
           IF (K.NE.0.AND.K.EQ.2*INT(K/2)) CSK0=CSK0+4.0D0*CF/K
           CBS=CBS+2.0D0*CF
           CF0=CF1
15         CF1=CF
        CS0=CDEXP(Z1)/(CBS-CF)
        DO 20 K=0,NM
20         CBI(K)=CS0*CBI(K)
        IF (A0.LE.9.0) THEN
           CBK(0)=-(CDLOG(0.5D0*Z1)+EL)*CBI(0)+CS0*CSK0
           CBK(1)=(1.0D0/Z1-CBI(1)*CBK(0))/CBI(0)
        ELSE
           CA0=CDSQRT(PI/(2.0D0*Z1))*CDEXP(-Z1)
           K0=16
           IF (A0.GE.25.0) K0=10
           IF (A0.GE.80.0) K0=8
           IF (A0.GE.200.0) K0=6
           DO 30 L=0,1
              CBKL=1.0D0
              VT=4.0D0*L
              CR=(1.0D0,0.0D0)
              DO 25 K=1,K0
                 CR=0.125D0*CR*(VT-(2.0*K-1.0)**2)/(K*Z1)
25               CBKL=CBKL+CR
              CBK(L)=CA0*CBKL
30         CONTINUE
        ENDIF
        CG0=CBK(0)
        CG1=CBK(1)
        DO 35 K=2,NM
           CG=2.0D0*(K-1.0D0)/Z1*CG1+CG0
           CBK(K)=CG
           CG0=CG1
35         CG1=CG
        IF (REAL(Z).LT.0.0) THEN
           FAC=1.0D0
           DO 45 K=0,NM
              IF (DIMAG(Z).LT.0.0) THEN
                 CBK(K)=FAC*CBK(K)+CI*PI*CBI(K)
              ELSE
                 CBK(K)=FAC*CBK(K)-CI*PI*CBI(K)
              ENDIF
              CBI(K)=FAC*CBI(K)
              FAC=-FAC
45         CONTINUE
        ENDIF
        CDI(0)=CBI(1)
        CDK(0)=-CBK(1)
        DO 50 K=1,NM
           CDI(K)=CBI(K-1)-K/Z*CBI(K)
50         CDK(K)=-CBK(K-1)-K/Z*CBK(K)
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
