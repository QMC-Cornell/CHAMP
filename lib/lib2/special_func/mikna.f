        PROGRAM MIKNA

C       =============================================================
C       Purpose: This program computes modified Bessel functions
C                In(x) and Kn(x), and their derivatives using
C                subroutine IKNA
C       Input:   x --- Argument of In(x) and Kn(x) ( x ע 0 )
C                n --- Order of In(x) and Kn(x)
C                      ( n = 0,1,תתת, n ף 250 )
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
        CALL IKNA(N,X,NM,BI,DI,BK,DK)
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


        SUBROUTINE IKNA(N,X,NM,BI,DI,BK,DK)

C       ========================================================
C       Purpose: Compute modified Bessel functions In(x) and
C                Kn(x), and their derivatives
C       Input:   x --- Argument of In(x) and Kn(x) ( x ע 0 )
C                n --- Order of In(x) and Kn(x)
C       Output:  BI(n) --- In(x)
C                DI(n) --- In'(x)
C                BK(n) --- Kn(x)
C                DK(n) --- Kn'(x)
C                NM --- Highest order computed
C       Routines called:
C            (1) IK01A for computing I0(x),I1(x),K0(x) & K1(x)
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:N),DI(0:N),BK(0:N),DK(0:N)
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
        CALL IK01A(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
        BI(0)=BI0
        BI(1)=BI1
        BK(0)=BK0
        BK(1)=BK1
        DI(0)=DI0
        DI(1)=DI1
        DK(0)=DK0
        DK(1)=DK1
        IF (N.LE.1) RETURN
        IF (X.GT.40.0.AND.N.LT.INT(0.25*X)) THEN
           H0=BI0
           H1=BI1
           DO 15 K=2,N
           H=-2.0D0*(K-1.0D0)/X*H1+H0
           BI(K)=H
           H0=H1
15         H1=H
        ELSE
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D-100
           DO 20 K=M,0,-1
              F=2.0D0*(K+1.0D0)*F1/X+F0
              IF (K.LE.NM) BI(K)=F
              F0=F1
20            F1=F
           S0=BI0/F
           DO 25 K=0,NM
25            BI(K)=S0*BI(K)
        ENDIF
        G0=BK0
        G1=BK1
        DO 30 K=2,NM
           G=2.0D0*(K-1.0D0)/X*G1+G0
           BK(K)=G
           G0=G1
30         G1=G
        DO 40 K=2,NM
           DI(K)=BI(K-1)-K/X*BI(K)
40         DK(K)=-BK(K-1)-K/X*BK(K)
        RETURN
        END


        SUBROUTINE IK01A(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)

C       =========================================================
C       Purpose: Compute modified Bessel functions I0(x), I1(1),
C                K0(x) and K1(x), and their derivatives
C       Input :  x   --- Argument ( x ע 0 )
C       Output:  BI0 --- I0(x)
C                DI0 --- I0'(x)
C                BI1 --- I1(x)
C                DI1 --- I1'(x)
C                BK0 --- K0(x)
C                DK0 --- K0'(x)
C                BK1 --- K1(x)
C                DK1 --- K1'(x)
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(12),B(12),A1(8)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           BI0=1.0D0
           BI1=0.0D0
           BK0=1.0D+300
           BK1=1.0D+300
           DI0=0.0D0
           DI1=0.5D0
           DK0=-1.0D+300
           DK1=-1.0D+300
           RETURN
        ELSE IF (X.LE.18.0) THEN
           BI0=1.0D0
           R=1.0D0
           DO 15 K=1,50
              R=0.25D0*R*X2/(K*K)
              BI0=BI0+R
              IF (DABS(R/BI0).LT.1.0D-15) GO TO 20
15         CONTINUE
20         BI1=1.0D0
           R=1.0D0
           DO 25 K=1,50
              R=0.25D0*R*X2/(K*(K+1))
              BI1=BI1+R
              IF (DABS(R/BI1).LT.1.0D-15) GO TO 30
25         CONTINUE
30         BI1=0.5D0*X*BI1
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
           IF (X.GE.35.0) K0=9
           IF (X.GE.50.0) K0=7
           CA=DEXP(X)/DSQRT(2.0D0*PI*X)
           BI0=1.0D0
           XR=1.0D0/X
           DO 35 K=1,K0
35            BI0=BI0+A(K)*XR**K
           BI0=CA*BI0
           BI1=1.0D0
           DO 40 K=1,K0
40            BI1=BI1+B(K)*XR**K
           BI1=CA*BI1
        ENDIF
        IF (X.LE.9.0D0) THEN
           CT=-(DLOG(X/2.0D0)+EL)
           BK0=0.0D0
           W0=0.0D0
           R=1.0D0
           DO 65 K=1,50
              W0=W0+1.0D0/K
              R=0.25D0*R/(K*K)*X2
              BK0=BK0+R*(W0+CT)
              IF (DABS((BK0-WW)/BK0).LT.1.0D-15) GO TO 70
65            WW=BK0
70         BK0=BK0+CT
        ELSE
           DATA A1/0.125D0,0.2109375D0,
     &             1.0986328125D0,1.1775970458984D01,
     &             2.1461706161499D02,5.9511522710323D03,
     &             2.3347645606175D05,1.2312234987631D07/
           CB=0.5D0/X
           XR2=1.0D0/X2
           BK0=1.0D0
           DO 75 K=1,8
75            BK0=BK0+A1(K)*XR2**K
           BK0=CB*BK0/BI0
        ENDIF
        BK1=(1.0D0/X-BI1*BK0)/BI0
        DI0=BI1
        DI1=BI0-BI1/X
        DK0=-BK1
        DK1=-BK0-BK1/X
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
