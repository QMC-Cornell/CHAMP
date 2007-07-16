        PROGRAM MCH12N

C       =====================================================
C       Purpose: This program computes Hankel functions of
C                the first and second kinds and their
C                derivatives for a complex argument using
C                subroutine CH12N
C       Input :  z --- Complex argument
C                n --- Order of Hn(1)(z) and Hn(2)(z)
C                      ( n = 0,1,תתת, n ף 250 )
C       Output:  CHF1(n) --- Hn(1)(z)
C                CHD1(n) --- Hn(1)'(z)
C                CHF2(n) --- Hn(2)(z)
C                CHD2(n) --- Hn(2)'(z)
C       =====================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CHF1(0:250),CHD1(0:250),CHF2(0:250),CHD2(0:250)
        WRITE(*,*)'  Please enter n, x and y (z=x+iy) '
        READ(*,*)N,X,Y
        WRITE(*,45)X,Y,N
        Z=CMPLX(X,Y)
        IF (N.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL CH12N(N,Z,NM,CHF1,CHD1,CHF2,CHD2)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Hn(1)(z)]     Im[Hn(1)(z)]',
     &            '      Re[Hn(1)''(z)]     Im[Hn(1)''(z)]'
        WRITE(*,*)' -------------------------------------',
     &               '---------------------------------------'
        DO 30 K=0,NM,NS
30         WRITE(*,40)K,CHF1(K),CHD1(K)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Hn(2)(z)]     Im[Hn(2)(z)]',
     &            '      Re[Hn(2)''(z)]     Im[Hn(2)''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '---------------------------------------'
        DO 35 K=0,NM,NS
35         WRITE(*,40)K,CHF2(K),CHD2(K)
40      FORMAT(1X,I4,4D18.10)
45      FORMAT(3X,3Hz =,F8.3,' + i ',F8.3,' ,',6X,6HNmax =,I4)
        END


        SUBROUTINE CH12N(N,Z,NM,CHF1,CHD1,CHF2,CHD2)

C       ====================================================
C       Purpose: Compute Hankel functions of the first and
C                second kinds and their derivatives for a
C                complex argument
C       Input :  z --- Complex argument
C                n --- Order of Hn(1)(z) and Hn(2)(z)
C       Output:  CHF1(n) --- Hn(1)(z)
C                CHD1(n) --- Hn(1)'(z)
C                CHF2(n) --- Hn(2)(z)
C                CHD2(n) --- Hn(2)'(z)
C                NM --- Highest order computed
C       Routines called:
C             (1) CJYNB for computing Jn(z) and Yn(z)
C             (2) CIKNB for computing In(z) and Kn(z)
C       ====================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:250),CDJ(0:250),CBY(0:250),CDY(0:250),
     &            CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)
        DIMENSION CHF1(0:N),CHD1(0:N),CHF2(0:N),CHD2(0:N)
        CI=(0.0D0,1.0D0)
        PI=3.141592653589793D0
        IF (DIMAG(Z).LT.0.0D0) THEN
           CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
           DO 10 K=0,NM
              CHF1(K)=CBJ(K)+CI*CBY(K)
10            CHD1(K)=CDJ(K)+CI*CDY(K)
           ZI=CI*Z
           CALL CIKNB(N,ZI,NM,CBI,CDI,CBK,CDK)
           CFAC=-2.0D0/(PI*CI)
           DO 15 K=0,NM
              CHF2(K)=CFAC*CBK(K)
              CHD2(K)=CFAC*CI*CDK(K)
15            CFAC=CFAC*CI
        ELSE IF (DIMAG(Z).GT.0.0D0) THEN
           ZI=-CI*Z
           CALL CIKNB(N,ZI,NM,CBI,CDI,CBK,CDK)
           CF1=-CI
           CFAC=2.0D0/(PI*CI)
           DO 20 K=0,NM
              CHF1(K)=CFAC*CBK(K)
              CHD1(K)=-CFAC*CI*CDK(K)
20            CFAC=CFAC*CF1
           CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
           DO 25 K=0,NM
              CHF2(K)=CBJ(K)-CI*CBY(K)
25            CHD2(K)=CDJ(K)-CI*CDY(K)
        ELSE
           CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
           DO 30 K=0,NM
              CHF1(K)=CBJ(K)+CI*CBY(K)
              CHD1(K)=CDJ(K)+CI*CDY(K)
              CHF2(K)=CBJ(K)-CI*CBY(K)
30            CHD2(K)=CDJ(K)-CI*CDY(K)
        ENDIF
        RETURN
        END


        SUBROUTINE CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)

C       =======================================================
C       Purpose: Compute Bessel functions Jn(z), Yn(z) and
C                their derivatives for a complex argument
C       Input :  z --- Complex argument of Jn(z) and Yn(z)
C                n --- Order of Jn(z) and Yn(z)
C       Output:  CBJ(n) --- Jn(z)
C                CDJ(n) --- Jn'(z)
C                CBY(n) --- Yn(z)
C                CDY(n) --- Yn'(z)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 to calculate the starting
C                point for backward recurrence
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:N),CDJ(0:N),CBY(0:N),CDY(0:N),
     &            A(4),B(4),A1(4),B1(4)
        EL=0.5772156649015329D0
        PI=3.141592653589793D0
        R2P=.63661977236758D0
        Y0=DABS(DIMAG(Z))
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           CBJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        IF (A0.LE.300.D0.OR.N.GT.80) THEN
           IF (N.EQ.0) NM=1
           M=MSTA1(A0,200)
           IF (M.LT.NM) THEN
              NM=M
           ELSE
              M=MSTA2(A0,NM,15)
           ENDIF
           CBS=(0.0D0,0.0D0)
           CSU=(0.0D0,0.0D0)
           CSV=(0.0D0,0.0D0)
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 15 K=M,0,-1
              CF=2.0D0*(K+1.0D0)/Z*CF1-CF2
              IF (K.LE.NM) CBJ(K)=CF
              IF (K.EQ.2*INT(K/2).AND.K.NE.0) THEN
                 IF (Y0.LE.1.0D0) THEN
                    CBS=CBS+2.0D0*CF
                 ELSE
                    CBS=CBS+(-1)**(K/2)*2.0D0*CF
                 ENDIF
                 CSU=CSU+(-1)**(K/2)*CF/K
              ELSE IF (K.GT.1) THEN
                 CSV=CSV+(-1)**(K/2)*K/(K*K-1.0D0)*CF
              ENDIF
              CF2=CF1
15            CF1=CF
           IF (Y0.LE.1.0D0) THEN
              CS0=CBS+CF
           ELSE
              CS0=(CBS+CF)/CDCOS(Z)
           ENDIF
           DO 20 K=0,NM
20            CBJ(K)=CBJ(K)/CS0
           CE=CDLOG(Z/2.0D0)+EL
           CBY(0)=R2P*(CE*CBJ(0)-4.0D0*CSU/CS0)
           CBY(1)=R2P*(-CBJ(0)/Z+(CE-1.0D0)*CBJ(1)-4.0D0*CSV/CS0)
        ELSE
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02/
           CT1=Z-0.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 25 K=1,4
25            CP0=CP0+A(K)*Z**(-2*K)
           CQ0=-0.125D0/Z
           DO 30 K=1,4
30            CQ0=CQ0+B(K)*Z**(-2*K-1)
           CU=CDSQRT(R2P/Z)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CBJ(0)=CBJ0
           CBY(0)=CBY0
           CT2=Z-0.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 35 K=1,4
35            CP1=CP1+A1(K)*Z**(-2*K)
           CQ1=0.375D0/Z
           DO 40 K=1,4
40            CQ1=CQ1+B1(K)*Z**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
           CBJ(1)=CBJ1
           CBY(1)=CBY1
           DO 45 K=2,NM
              CBJK=2.0D0*(K-1.0D0)/Z*CBJ1-CBJ0
              CBJ(K)=CBJK
              CBJ0=CBJ1
45            CBJ1=CBJK
        ENDIF
        CDJ(0)=-CBJ(1)
        DO 50 K=1,NM
50         CDJ(K)=CBJ(K-1)-K/Z*CBJ(K)
        IF (CDABS(CBJ(0)).GT.1.0D0) THEN
           CBY(1)=(CBJ(1)*CBY(0)-2.0D0/(PI*Z))/CBJ(0)
        ENDIF
        DO 55 K=2,NM
           IF (CDABS(CBJ(K-1)).GE.CDABS(CBJ(K-2))) THEN
              CYY=(CBJ(K)*CBY(K-1)-2.0D0/(PI*Z))/CBJ(K-1)
           ELSE
              CYY=(CBJ(K)*CBY(K-2)-4.0D0*(K-1.0D0)/(PI*Z*Z))/CBJ(K-2)
           ENDIF
           CBY(K)=CYY
55      CONTINUE
        CDY(0)=-CBY(1)
        DO 60 K=1,NM
60         CDY(K)=CBY(K-1)-K/Z*CBY(K)
        RETURN
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
