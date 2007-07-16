        PROGRAM MSPHJ

C       =======================================================
C       Purpose: This program computes the spherical Bessel
C                functions jn(x) and jn'(x) using subroutine
C                SPHJ
C       Input :  x --- Argument of jn(x)
C                n --- Order of jn(x)  ( n = 0,1,תתת,ף 250 )
C       Output:  SJ(n) --- jn(x)
C                DJ(n) --- jn'(x)
C       Example:   x =10.0
C                  n          jn(x)               jn(x)
C                --------------------------------------------
C                  0    -.5440211109D-01    -.7846694180D-01
C                  1     .7846694180D-01    -.7009549945D-01
C                  2     .7794219363D-01     .5508428371D-01
C                  3    -.3949584498D-01     .9374053162D-01
C                  4    -.1055892851D+00     .1329879757D-01
C                  5    -.5553451162D-01    -.7226857814D-01
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SJ(0:250),DJ(0:250)
        WRITE(*,*)'Please enter n and x '
        READ(*,*)N,X
        WRITE(*,30)N,X
        IF (N.LE.10) THEN
           NS=1
        ELSE
           WRITE(*,*)'Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL SPHJ(N,X,NM,SJ,DJ)
        WRITE(*,*)
        WRITE(*,*)'  n          jn(x)               jn''(x)'
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,SJ(K),DJ(K)
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,6HNmax =,I3,',     ',2Hx=,F5.1)
        END


        SUBROUTINE SPHJ(N,X,NM,SJ,DJ)

C       =======================================================
C       Purpose: Compute spherical Bessel functions jn(x) and
C                their derivatives
C       Input :  x --- Argument of jn(x)
C                n --- Order of jn(x)  ( n = 0,1,תתת )
C       Output:  SJ(n) --- jn(x)
C                DJ(n) --- jn'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SJ(0:N),DJ(0:N)
        NM=N
        IF (DABS(X).EQ.1.0D-100) THEN
           DO 10 K=0,N
              SJ(K)=0.0D0
10            DJ(K)=0.0D0
           SJ(0)=1.0D0
           DJ(1)=.3333333333333333D0
           RETURN
        ENDIF
        SJ(0)=DSIN(X)/X
        SJ(1)=(SJ(0)-DCOS(X))/X
        IF (N.GE.2) THEN
           SA=SJ(0)
           SB=SJ(1)
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D0-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X-F0
              IF (K.LE.NM) SJ(K)=F
              F0=F1
15            F1=F
           IF (DABS(SA).GT.DABS(SB)) CS=SA/F
           IF (DABS(SA).LE.DABS(SB)) CS=SB/F0
           DO 20 K=0,NM
20            SJ(K)=CS*SJ(K)
        ENDIF
        DJ(0)=(DCOS(X)-DSIN(X)/X)/X
        DO 25 K=1,NM
25         DJ(K)=SJ(K-1)-(K+1.0D0)*SJ(K)/X
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
