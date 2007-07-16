        PROGRAM MLQNA

C       ======================================================
C       Purpose: This program computes the Legendre functions
C                Qn(x) and Qn'(x) using subroutine LQNA
C       Input :  x  --- Argument of Qn(x)  ( -1 ף x ף 1 )
C                n  --- Degree of Qn(x)  ( n = 0,1,תתת )
C       Output:  QN(n) --- Qn(x)
C                QD(n) --- Qn'(x)
C       Example:  x = 0.50
C                 n        Qn(x)         Qn'(x)
C                ---------------------------------
C                 0      .54930614     1.33333333
C                 1     -.72534693     1.21597281
C                 2     -.81866327     -.84270745
C                 3     -.19865477    -2.87734353
C                 4      .44017453    -2.23329085
C                 5      .55508089     1.08422720
C       ======================================================

        DOUBLE PRECISION QN,QD,X
        DIMENSION QN(0:100),QD(0:100)
        WRITE(*,*)'  Please enter Nmax and x'
        READ(*,*)N,X
        WRITE(*,30) X
        WRITE(*,*)
        CALL LQNA(N,X,QN,QD)
        WRITE(*,*)'  n        Qn(x)         Qn''(x)'
        WRITE(*,*)' ---------------------------------'
        DO 10 K=0,N
10         WRITE(*,20)K,QN(K),QD(K)
20      FORMAT(1X,I3,2F15.8)
30      FORMAT(3X,'x =',F5.2)
        END


        SUBROUTINE LQNA(N,X,QN,QD)

C       =====================================================
C       Purpose: Compute Legendre functions Qn(x) and Qn'(x)
C       Input :  x  --- Argument of Qn(x) ( -1 ף x ף 1 )
C                n  --- Degree of Qn(x) ( n = 0,1,2,תתת )
C       Output:  QN(n) --- Qn(x)
C                QD(n) --- Qn'(x)
C                ( 1.0D+300 stands for infinity )
C       =====================================================

        IMPLICIT DOUBLE PRECISION (Q,X)
        DIMENSION QN(0:N),QD(0:N)
        IF (DABS(X).EQ.1.0D0) THEN
           DO 10 K=0,N
              QN(K)=1.0D+300
              QD(K)=-1.0D+300
10         CONTINUE
        ELSE IF (DABS(X).LT.1.0D0) THEN
           Q0=0.5D0*DLOG((1.0D0+X)/(1.0D0-X))
           Q1=X*Q0-1.0D0
           QN(0)=Q0
           QN(1)=Q1
           QD(0)=1.0D0/(1.0D0-X*X)
           QD(1)=QN(0)+X*QD(0)
           DO 15 K=2,N
              QF=((2*K-1)*X*Q1-(K-1)*Q0)/K
              QN(K)=QF
              QD(K)=(QN(K-1)-X*QF)*K/(1.0D0-X*X)
              Q0=Q1
15            Q1=QF
        ENDIF
        RETURN
        END
