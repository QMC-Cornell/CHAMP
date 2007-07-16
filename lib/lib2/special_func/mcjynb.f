        PROGRAM MCJYNB

C       ================================================================
C       Purpose: This program computes Bessel functions Jn(z), Yn(z)
C                and their derivatives for a complex argument using
C                subroutine CJYNB
C       Input :  z --- Complex argument of Jn(z) and Yn(z)
C                n --- Order of Jn(z) and Yn(z)
C                      ( n = 0,1,תתת, n ף 250 )
C       Output:  CBJ(n) --- Jn(z)
C                CDJ(n) --- Jn'(z)
C                CBY(n) --- Yn(z)
C                CDY(n) --- Yn'(z)
C       Eaxmple: z = 4.0 + i 2.0

C     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
C    -------------------------------------------------------------------
C     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
C     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
C     2   .93050039D+00  -.77959350D+00  -.72363400D+00  -.72836666D+00
C     3   .93991546D+00   .23042918D+00   .29742236D+00  -.63587637D+00
C     4   .33565567D+00   .49215925D+00   .47452722D+00  -.29035945D-01
C     5  -.91389835D-02   .28850107D+00   .20054412D+00   .19908868D+00

C     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
C   --------------------------------------------------------------------
C     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
C     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
C     2   .79074211D+00   .86842120D+00   .78932897D+00  -.70142425D+00
C     3  -.29934789D+00   .89064431D+00   .70315755D+00   .24423024D+00
C     4  -.61557299D+00   .37996071D+00   .41126221D-01   .34044655D+00
C     5  -.38160033D+00   .20975121D+00  -.33884827D+00  -.20590670D-01

C                z = 20.0 + i  10.0 ,      Nmax =  5

C     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
C   --------------------------------------------------------------------
C     0   .15460268D+04  -.10391216D+04  -.10601232D+04  -.15098284D+04
C     1   .10601232D+04   .15098284D+04   .14734253D+04  -.10783122D+04
C     2  -.14008238D+04   .11175029D+04   .11274890D+04   .13643952D+04
C     3  -.11948548D+04  -.12189620D+04  -.11843035D+04   .11920871D+04
C     4   .96778325D+03  -.12666712D+04  -.12483664D+04  -.93887194D+03
C     5   .13018781D+04   .65878188D+03   .64152944D+03  -.12682398D+04

C     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
C   --------------------------------------------------------------------
C     0   .10391216D+04   .15460268D+04   .15098284D+04  -.10601232D+04
C     1  -.15098284D+04   .10601232D+04   .10783122D+04   .14734253D+04
C     2  -.11175029D+04  -.14008238D+04  -.13643952D+04   .11274890D+04
C     3   .12189620D+04  -.11948548D+04  -.11920871D+04  -.11843035D+04
C     4   .12666712D+04   .96778324D+03   .93887194D+03  -.12483664D+04
C     5  -.65878189D+03   .13018781D+04   .12682398D+04   .64152944D+03
C       ================================================================

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        COMMON CBJ(0:250),CDJ(0:250),CBY(0:250),CDY(0:250)
        WRITE(*,*)'  Please enter n, x,y (z=x+iy) '
        READ(*,*)N,X,Y
        Z=CMPLX(X,Y)
        WRITE(*,40)X,Y,N
        IF (N.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Jn(z)]       Im[Jn(z)]',
     &            '       Re[Jn''(z)]      Im[Jn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,30)K,CBJ(K),CDJ(K)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Yn(z)]       Im[Yn(z)]',
     &            '       Re[Yn''(z)]      Im[Yn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        DO 20 K=0,NM,NS
20         WRITE(*,30)K,CBY(K),CDY(K)
30      FORMAT(1X,I4,4D16.8)
40      FORMAT(3X,3Hz =,F5.1,' + i ',F5.1,' ,',6X,6HNmax =,I3)
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
        IF (A0.LE.300.D0.OR.N.GT.INT(0.25*A0)) THEN
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
