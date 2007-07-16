        PROGRAM MLAMN

C       ====================================================
C       Purpose: This program computes the lambda functions
C                and their derivatives using subroutine
C                LAMN
C       Input:   x --- Argument of lambda function
C                n --- Order of lambda function
C                      ( n = 0,1,..., n ó 250 )
C       Output:  BL(n) --- Lambda function of order n
C                DL(n) --- Derivative of lambda function
C       Example: Nmax = 5,  x = 10.00

C                 n       lambda(x)        lambda'(x)
C                ---------------------------------------
C                 0    -.24593576D+00    -.43472746D-01
C                 1     .86945492D-02    -.50926063D-01
C                 2     .20370425D-01    -.46703503D-02
C                 3     .28022102D-02     .10540929D-01
C                 4    -.84327431D-02     .89879627D-02
C                 5    -.89879627D-02     .55521954D-03
C       ====================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BL(0:250),DL(0:250)
        WRITE(*,*)'  Please enter n,x = ?'
        READ(*,*)N,X
        WRITE(*,15)N,X
        IF (N.LE.10) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL LAMN(N,X,NM,BL,DL)
        WRITE(*,*)
        WRITE(*,*) '  n       lambda(x)        lambda''(x)'
        WRITE(*,*)' ---------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,BL(K),DL(K)
15      FORMAT(1X,3HN =,I4,6X,3Hx =,F8.2)
20      FORMAT(1X,I3,2D18.8)
        END


        SUBROUTINE LAMN(N,X,NM,BL,DL)

C       =========================================================
C       Purpose: Compute lambda functions and their derivatives
C       Input:   x --- Argument of lambda function
C                n --- Order of lambda function
C       Output:  BL(n) --- Lambda function of order n
C                DL(n) --- Derivative of lambda function
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the start
C                point for backward recurrence
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BL(0:N),DL(0:N)
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              BL(K)=0.0D0
10            DL(K)=0.0D0
           BL(0)=1.0D0
           DL(1)=0.5D0
           RETURN
        ENDIF
        IF (X.LE.12.0D0) THEN
           X2=X*X
           DO 25 K=0,N
              BK=1.0D0
              R=1.0D0
              DO 15 I=1,50
                 R=-0.25D0*R*X2/(I*(I+K))
                 BK=BK+R
                 IF (DABS(R).LT.DABS(BK)*1.0D-15) GO TO 20
15            CONTINUE
20            BL(K)=BK
25            IF (K.GE.1) DL(K-1)=-0.5D0*X/K*BK
           UK=1.0D0
           R=1.0D0
           DO 30 I=1,50
              R=-0.25D0*R*X2/(I*(I+N+1.0D0))
              UK=UK+R
              IF (DABS(R).LT.DABS(UK)*1.0D-15) GO TO 35
30            CONTINUE
35         DL(N)=-0.5D0*X/(N+1.0D0)*UK
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
        F0=0.0D0
        F1=1.0D-100
        DO 40 K=M,0,-1
           F=2.0D0*(K+1.0D0)*F1/X-F0
           IF (K.LE.NM) BL(K)=F
           IF (K.EQ.2*INT(K/2)) BS=BS+2.0D0*F
           F0=F1
40         F1=F
        BG=BS-F
        DO 45 K=0,NM
45         BL(K)=BL(K)/BG
        R0=1.0D0
        DO 50 K=1,NM
           R0=2.0D0*R0*K/X
50         BL(K)=R0*BL(K)
        DL(0)=-0.5D0*X*BL(1)
        DO 55 K=1,NM
55         DL(K)=2.0D0*K/X*(BL(K-1)-BL(K))
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
