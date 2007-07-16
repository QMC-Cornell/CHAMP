        PROGRAM MLQNB

C       ===============================================================
C       Purpose: This program computes the Legendre functions Qn(x)
C                and Qn'(x) using subroutine LQNB
C       Input :  x  --- Argument of Qn(x)
C                n  --- Degree of Qn(x)  ( n = 0,1,תתת)
C       Output:  QN(n) --- Qn(x)
C                QD(n) --- Qn'(x)
C       Examples:     x1 = 0.50,    x2 = 2.50

C       n      Qn(x1)        Qn'(x1)       Qn(x2)          Qn'(x2)
C     ----------------------------------------------------------------
C       0     .54930614    1.33333333   .42364893D+00  -.19047619D+00
C       1    -.72534693    1.21597281   .59122325D-01  -.52541546D-01
C       2    -.81866327    -.84270745   .98842555D-02  -.13109214D-01
C       3    -.19865477   -2.87734353   .17695141D-02  -.31202687D-02
C       4     .44017453   -2.23329085   .32843271D-03  -.72261513D-03
C       5     .55508089    1.08422720   .62335892D-04  -.16437427D-03
C       ===============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION QN(0:100),QD(0:100)
        WRITE(*,*)'Please enter Nmax and x '
        READ(*,*)N,X
        WRITE(*,40)X
        WRITE(*,*)
        WRITE(*,*)'  n          Qn(x)           Qn''(x)'
        WRITE(*,*)'--------------------------------------'
        CALL LQNB(N,X,QN,QD)
        DO 10 K=0,N
           IF (X.LE.1.0) THEN
              WRITE(*,20)K,QN(K),QD(K)
           ELSE
              WRITE(*,30)K,QN(K),QD(K)
           ENDIF
10      CONTINUE
20      FORMAT(1X,I3,2F17.8)
30      FORMAT(1X,I3,2D17.8)
40      FORMAT(3X,'x =',F5.2)
        END


        SUBROUTINE LQNB(N,X,QN,QD)

C       ====================================================
C       Purpose: Compute Legendre functions Qn(x) & Qn'(x)
C       Input :  x  --- Argument of Qn(x)
C                n  --- Degree of Qn(x)  ( n = 0,1,2,תתת)
C       Output:  QN(n) --- Qn(x)
C                QD(n) --- Qn'(x)
C       ====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION QN(0:N),QD(0:N)
        EPS=1.0D-14
        IF (DABS(X).EQ.1.0D0) THEN
           DO 10 K=0,N
              QN(K)=1.0D+300
10            QD(K)=1.0D+300
           RETURN
        ENDIF
        IF (X.LE.1.021D0) THEN
           X2=DABS((1.0D0+X)/(1.0D0-X))
           Q0=0.5D0*DLOG(X2)
           Q1=X*Q0-1.0D0
           QN(0)=Q0
           QN(1)=Q1
           QD(0)=1.0D0/(1.0D0-X*X)
           QD(1)=QN(0)+X*QD(0)
           DO 15 K=2,N
              QF=((2.0D0*K-1.0D0)*X*Q1-(K-1.0D0)*Q0)/K
              QN(K)=QF
              QD(K)=(QN(K-1)-X*QF)*K/(1.0D0-X*X)
              Q0=Q1
15            Q1=QF
        ELSE
           QC2=1.0D0/X
           DO 20 J=1,N
              QC2=QC2*J/((2.0*J+1.0D0)*X)
              IF (J.EQ.N-1) QC1=QC2
20         CONTINUE
           DO 35 L=0,1
              NL=N+L
              QF=1.0D0
              QR=1.0D0
              DO 25 K=1,500
                 QR=QR*(0.5D0*NL+K-1.0D0)*(0.5D0*(NL-1)+K)
     &              /((NL+K-0.5D0)*K*X*X)
                 QF=QF+QR
                 IF (DABS(QR/QF).LT.EPS) GO TO 30
25            CONTINUE
30            IF (L.EQ.0) THEN
                 QN(N-1)=QF*QC1
              ELSE
                 QN(N)=QF*QC2
              ENDIF
35         CONTINUE
           QF2=QN(N)
           QF1=QN(N-1)
           DO 40 K=N,2,-1
              QF0=((2*K-1.0D0)*X*QF1-K*QF2)/(K-1.0D0)
              QN(K-2)=QF0
              QF2=QF1
40            QF1=QF0
           QD(0)=1.0D0/(1.0D0-X*X)
           DO 45 K=1,N
45            QD(K)=K*(QN(K-1)-X*QN(K))/(1.0D0-X*X)
        ENDIF
        RETURN
        END
