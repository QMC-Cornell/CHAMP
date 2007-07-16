        PROGRAM MCSPHIK

C       =============================================================
C       Purpose: This program computes the modified spherical Bessel
C                functions and their derivatives for a complex
C                argument using subroutine CSPHIK
C       Input :  z --- Complex argument
C                n --- Order of in(z) & kn(z) ( 0 ó n ó 250 )
C       Output:  CSI(n) --- in(z)
C                CDI(n) --- in'(z)
C                CSK(n) --- kn(z)
C                CDK(n) --- kn'(z)
C       Example: z =4.0+i 2.0

C     n     Re[in(z)]      Im[in(z)]     Re[in'(z)]     Im[in'(z)]
C    ---------------------------------------------------------------
C     0   .2118080D+00   .6101922D+01  -.4439356D+00   .4900150D+01
C     1  -.4439356D+00   .4900150D+01  -.5906477D+00   .4053075D+01
C     2  -.9918756D+00   .3028652D+01  -.7574058D+00   .2785396D+01
C     3  -.9663859D+00   .1375561D+01  -.7689911D+00   .1541649D+01
C     4  -.6018277D+00   .4263967D+00  -.5777565D+00   .6482500D+00
C     5  -.2668530D+00   .6640148D-01  -.3214450D+00   .1866032D+00

C     n     Re[kn(z)]      Im[kn(z)]     Re[kn'(z)]     Im[kn'(z)]
C    ---------------------------------------------------------------
C     0  -.5010582D-02  -.4034862D-02   .6416184D-02   .4340777D-02
C     1  -.6416184D-02  -.4340777D-02   .8445211D-02   .4487936D-02
C     2  -.1016253D-01  -.4714473D-02   .1392804D-01   .4120703D-02
C     3  -.1893595D-01  -.3973987D-02   .2690088D-01   .3192843D-03
C     4  -.3945464D-01   .2977107D-02   .5690203D-01  -.1873044D-01
C     5  -.8727490D-01   .3689398D-01   .1220481D+00  -.9961483D-01
C       =============================================================

        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION X,Y
        DIMENSION CSI(0:250),CDI(0:250),CSK(0:250),CDK(0:250)
        WRITE(*,*)'Please enter n,x,y (z=x+iy) '
        READ(*,*)N,X,Y
        WRITE(*,30)N,X,Y
        Z=CMPLX(X,Y)
        IF (N.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)'Please enter order step Ns '
           READ(*,*)NS
        ENDIF
        CALL CSPHIK(N,Z,NM,CSI,CDI,CSK,CDK)
        WRITE(*,*)
        WRITE(*,*)'  n      Re[in(z)]        Im[in(z)]',
     &  '        Re[in''(z)]       Im[in''(z)]'
        WRITE(*,*)'--------------------------------------------',
     &  '----------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,CSI(K),CDI(K)
        WRITE(*,*)
        WRITE(*,*)'  n      Re[kn(z)]        Im[kn(z)]',
     &  '        Re[kn''(z)]       Im[kn''(z)]'
        WRITE(*,*)'--------------------------------------------',
     &  '----------------------------'
        DO 15 K=0,NM,NS
15         WRITE(*,20)K,CSK(K),CDK(K)
20      FORMAT(1X,I3,4D17.8)
30      FORMAT(3X,'Nmaz =',I3,',     ','z = ',F8.1,'+ i',F8.1)
        END


        SUBROUTINE CSPHIK(N,Z,NM,CSI,CDI,CSK,CDK)

C       =======================================================
C       Purpose: Compute modified spherical Bessel functions
C                and their derivatives for a complex argument
C       Input :  z --- Complex argument
C                n --- Order of in(z) & kn(z) ( n = 0,1,2,... )
C       Output:  CSI(n) --- in(z)
C                CDI(n) --- in'(z)
C                CSK(n) --- kn(z)
C                CDK(n) --- kn'(z)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       =======================================================

        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION A0,PI
        DIMENSION CSI(0:N),CDI(0:N),CSK(0:N),CDK(0:N)
        PI=3.141592653589793D0
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-60) THEN
           DO 10 K=0,N
              CSI(K)=0.0D0
              CDI(K)=0.0D0
              CSK(K)=1.0D+300
10            CDK(K)=-1.0D+300
           CSI(0)=1.0D0
           CDI(1)=0.3333333333333333D0
           RETURN
        ENDIF
        CI=CMPLX(0.0D0,1.0D0)
        CSINH=CDSIN(CI*Z)/CI
        CCOSH=CDCOS(CI*Z)
        CSI0=CSINH/Z
        CSI1=(-CSINH/Z+CCOSH)/Z
        CSI(0)=CSI0
        CSI(1)=CSI1
        IF (N.GE.2) THEN
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF0=0.0D0
           CF1=1.0D0-100
           DO 15 K=M,0,-1
              CF=(2.0D0*K+3.0D0)*CF1/Z+CF0
              IF (K.LE.NM) CSI(K)=CF
              CF0=CF1
15            CF1=CF
           IF (CDABS(CSI0).GT.CDABS(CSI1)) CS=CSI0/CF
           IF (CDABS(CSI0).LE.CDABS(CSI1)) CS=CSI1/CF0
           DO 20 K=0,NM
20            CSI(K)=CS*CSI(K)
        ENDIF
        CDI(0)=CSI(1)
        DO 25 K=1,NM
25         CDI(K)=CSI(K-1)-(K+1.0D0)*CSI(K)/Z
        CSK(0)=0.5D0*PI/Z*CDEXP(-Z)
        CSK(1)=CSK(0)*(1.0D0+1.0D0/Z)
        DO 30 K=2,NM
           IF (CDABS(CSI(K-1)).GT.CDABS(CSI(K-2))) THEN
              CSK(K)=(0.5D0*PI/(Z*Z)-CSI(K)*CSK(K-1))/CSI(K-1)
           ELSE
              CSK(K)=(CSI(K)*CSK(K-2)+(K-0.5D0)*PI/Z**3)/CSI(K-2)
           ENDIF
30      CONTINUE
        CDK(0)=-CSK(1)
        DO 35 K=1,NM
35         CDK(K)=-CSK(K-1)-(K+1.0D0)*CSK(K)/Z
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
