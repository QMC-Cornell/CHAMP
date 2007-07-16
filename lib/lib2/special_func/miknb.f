        PROGRAM MIKNB

C       =============================================================
C       Purpose: This program computes modified Bessel functions
C                In(x) and Kn(x), and their derivatives using
C                subroutine IKNB
C       Input:   x --- Argument of In(x) and Kn(x) ( 0 ó x ó 700 )
C                n --- Order of In(x) and Kn(x)
C                      ( n = 0,1,..., n ó 250 )
C       Output:  BI(n) --- In(x)
C                DI(n) --- In'(x)
C                BK(n) --- Kn(x)
C                DK(n) --- Kn'(x)
C       Example: Nmax = 5,    x = 10.0

C     n      In(x)          In'(x)         Kn(x)         Kn'(x)
C    ---------------------------------------------------------------
C     0   .2815717D+04   .2670988D+04   .1778006D-04  -.1864877D-04
C     1   .2670988D+04   .2548618D+04   .1864877D-04  -.1964494D-04
C     2   .2281519D+04   .2214685D+04   .2150982D-04  -.2295074D-04
C     3   .1758381D+04   .1754005D+04   .2725270D-04  -.2968563D-04
C     4   .1226491D+04   .1267785D+04   .3786144D-04  -.4239728D-04
C     5   .7771883D+03   .8378964D+03   .5754185D-04  -.6663236D-04
C       =============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:250),DI(0:250),BK(0:250),DK(0:250)
        WRITE(*,*)'  Please enter n, x '
        READ(*,*)N,X
        WRITE(*,25)N,X
        WRITE(*,*)
        IF (N.LE.10) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL IKNB(N,X,NM,BI,DI,BK,DK)
        WRITE(*,*)'  n      In(x)          In''(x) ',
     &            '        Kn(x)         Kn''(x) '
        WRITE(*,*)' -------------------------------',
     &            '--------------------------------'
        DO 10 K=0,NM,NS
           WRITE(*,20)K,BI(K),DI(K),BK(K),DK(K)
10      CONTINUE
20      FORMAT(1X,I3,4D15.7)
25      FORMAT(3X,6HNmax =,I3,',    ',3Hx =,F5.1)
        END


        SUBROUTINE IKNB(N,X,NM,BI,DI,BK,DK)

C       ============================================================
C       Purpose: Compute modified Bessel functions In(x) and Kn(x),
C                and their derivatives
C       Input:   x --- Argument of In(x) and Kn(x) ( 0 ó x ó 700 )
C                n --- Order of In(x) and Kn(x)
C       Output:  BI(n) --- In(x)
C                DI(n) --- In'(x)
C                BK(n) --- Kn(x)
C                DK(n) --- Kn'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting point
C                for backward recurrence
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:N),DI(0:N),BK(0:N),DK(0:N)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        NM=N
        IF (X.LE.1.0D-100) THEN
           DO 10 K=0,N
              BI(K)=0.0D0
              DI(K)=0.0D0
              BK(K)=1.0D+300
10            DK(K)=-1.0D+300
           BI(0)=1.0D0
           DI(1)=0.5D0
           RETURN
        ENDIF
        IF (N.EQ.0) NM=1
        M=MSTA1(X,200)
        IF (M.LT.NM) THEN
           NM=M
        ELSE
           M=MSTA2(X,NM,15)
        ENDIF
        BS=0.0D0
        SK0=0.0D0
        F0=0.0D0
        F1=1.0D-100
        DO 15 K=M,0,-1
           F=2.0D0*(K+1.0D0)/X*F1+F0
           IF (K.LE.NM) BI(K)=F
           IF (K.NE.0.AND.K.EQ.2*INT(K/2)) SK0=SK0+4.0D0*F/K
           BS=BS+2.0D0*F
           F0=F1
15         F1=F
        S0=DEXP(X)/(BS-F)
        DO 20 K=0,NM
20         BI(K)=S0*BI(K)
        IF (X.LE.8.0D0) THEN
           BK(0)=-(DLOG(0.5D0*X)+EL)*BI(0)+S0*SK0
           BK(1)=(1.0D0/X-BI(1)*BK(0))/BI(0)
        ELSE
           A0=DSQRT(PI/(2.0D0*X))*DEXP(-X)
           K0=16
           IF (X.GE.25.0) K0=10
           IF (X.GE.80.0) K0=8
           IF (X.GE.200.0) K0=6
           DO 30 L=0,1
              BKL=1.0D0
              VT=4.0D0*L
              R=1.0D0
              DO 25 K=1,K0
                 R=0.125D0*R*(VT-(2.0*K-1.0)**2)/(K*X)
25               BKL=BKL+R
              BK(L)=A0*BKL
30         CONTINUE
        ENDIF
        G0=BK(0)
        G1=BK(1)
        DO 35 K=2,NM
           G=2.0D0*(K-1.0D0)/X*G1+G0
           BK(K)=G
           G0=G1
35         G1=G
        DI(0)=BI(1)
        DK(0)=-BK(1)
        DO 40 K=1,NM
           DI(K)=BI(K-1)-K/X*BI(K)
40         DK(K)=-BK(K-1)-K/X*BK(K)
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
C       Input :  x     --- Argument of Jn(x)
C                n     --- Order of Jn(x)
C                MP    --- Significant digit
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
