        PROGRAM MRCTJ

C       =======================================================
C       Purpose: This program computes the Riccati-Bessel
C                functions of the first kind, and their
C                derivatives using subroutine RCTJ
C       Input:   x --- Argument of Riccati-Bessel function
C                n --- Order of jn(x)  ( 0 ó n ó 250 )
C       Output:  RJ(n) --- xújn(x)
C                DJ(n) --- [xújn(x)]'
C       Example: x = 10.0
C                  n        xújn(x)             [xújn(x)]'
C                --------------------------------------------
C                  0    -.5440211109D+00    -.8390715291D+00
C                  1     .7846694180D+00    -.6224880527D+00
C                  2     .7794219363D+00     .6287850307D+00
C                  3    -.3949584498D+00     .8979094712D+00
C                  4    -.1055892851D+01     .2739869063D-01
C                  5    -.5553451162D+00    -.7782202931D+00
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RJ(0:250),DJ(0:250)
        WRITE(*,*)'  Please enter n and x '
        READ(*,*)N,X
        WRITE(*,30)N,X
        IF (N.LE.10) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        WRITE(*,*)
        CALL RCTJ(N,X,NM,RJ,DJ)
        WRITE(*,*)
        WRITE(*,*)'  n        xújn(x)             [xújn(x)]'''
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,RJ(K),DJ(K)
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,6HNmax =,I3,',    ',3Hx =,F7.2)
        END


        SUBROUTINE RCTJ(N,X,NM,RJ,DJ)

C       ========================================================
C       Purpose: Compute Riccati-Bessel functions of the first
C                kind and their derivatives
C       Input:   x --- Argument of Riccati-Bessel function
C                n --- Order of jn(x)  ( n = 0,1,2,... )
C       Output:  RJ(n) --- xújn(x)
C                DJ(n) --- [xújn(x)]'
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RJ(0:N),DJ(0:N)
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              RJ(K)=0.0D0
10            DJ(K)=0.0D0
           DJ(0)=1.0D0
           RETURN
        ENDIF
        RJ(0)=DSIN(X)
        RJ(1)=RJ(0)/X-DCOS(X)
        RJ0=RJ(0)
        RJ1=RJ(1)
        IF (N.GE.2) THEN
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X-F0
              IF (K.LE.NM) RJ(K)=F
              F0=F1
15            F1=F
           IF (DABS(RJ0).GT.DABS(RJ1)) CS=RJ0/F
           IF (DABS(RJ0).LE.DABS(RJ1)) CS=RJ1/F0
           DO 20 K=0,NM
20            RJ(K)=CS*RJ(K)
        ENDIF
        DJ(0)=DCOS(X)
        DO 25 K=1,NM
25         DJ(K)=-K*RJ(K)/X+RJ(K-1)
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
