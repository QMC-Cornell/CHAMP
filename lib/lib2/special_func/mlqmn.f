        PROGRAM MLQMN

C       ===============================================================
C       Purpose: This program computes the associated Legendre
C                functions Qmn(x) and their derivatives Qmn'(x) using
C                subroutine LQMN
C       Input :  x --- Argument of Qmn(x)
C                m --- Order of Qmn(x)  ( m = 0,1,2,תתת )
C                n --- Degree of Qmn(x) ( n = 0,1,2,תתת )
C       Output:  QM(m,n) --- Qmn(x)
C                QD(m,n) --- Qmn'(x)
C       Examples:

C       Qmn(x):  x = 0.5
C       n\m      0           1           2           3           4
C       ---------------------------------------------------------------
C        0     .549306   -1.154701    1.333333   -5.388603   26.666667
C        1    -.725347   -1.053063    2.666667   -6.158403   32.000000
C        2    -.818663     .729806    4.069272  -12.316806   42.666667
C        3    -.198655    2.491853    -.493486  -23.778868   85.333333
C        4     .440175    1.934087  -11.036781   -9.325204  186.818394

C       Qmn'(x): x = 0.5
C       n\m      0           1           2           3           4
C       ---------------------------------------------------------------
C        0    1.333333    -.769800    4.444444  -20.014809  145.777778
C        1    1.215973   -2.377159    3.555556  -24.633611  156.444444
C        2    -.842707   -5.185328    8.796526  -24.633611  199.111111
C        3   -2.877344   -1.091406   28.115454  -50.976710  227.555556
C        4   -2.233291   11.454786   25.483527 -197.068892  412.039838

C       Qmn(x): x = 2.0
C       n\m      0           1           2           3           4
C       ---------------------------------------------------------------
C        0     .549306    -.577350    1.333333   -5.003702   26.666667
C        1     .098612    -.203274     .666667   -3.079201   18.666667
C        2     .021184    -.064946     .277089   -1.539601   10.666667
C        3     .004871    -.019817     .104220    -.679543    5.333333
C        4     .001161    -.005887     .036816    -.276005    2.427640

C       Qmn'(x): x = 2.0
C       n\m      0           1           2           3           4
C       ---------------------------------------------------------------
C        0    -.333333     .384900   -1.111111    5.388603  -36.444444
C        1    -.117361     .249384    -.888889    4.618802  -32.000000
C        2    -.037496     .116680    -.519437    3.079201  -23.111111
C        3    -.011442     .046960    -.253375    1.720114  -14.222222
C        4    -.003399     .017331    -.110263     .849589   -7.748516
C       ===============================================================

        IMPLICIT DOUBLE PRECISION (Q,X)
        DIMENSION QM(0:100,0:100),QD(0:100,0:100)
        WRITE(*,*)'  Please enter m, n and x'
        READ(*,*) M,N,X
        WRITE(*,*)
        WRITE(*,*)'  m     n      x          Qmn(x)         Qmn''(x)'
        WRITE(*,*)' ---------------------------------------------------'
        CALL LQMN(100,M,N,X,QM,QD)
        DO 15 J=0,N
           WRITE(*,10)M,J,X,QM(M,J),QD(M,J)
15      CONTINUE
10      FORMAT(1X,I3,3X,I3,3X,F5.1,2D17.8)
        END


        SUBROUTINE LQMN(MM,M,N,X,QM,QD)

C       ==========================================================
C       Purpose: Compute the associated Legendre functions of the
C                second kind, Qmn(x) and Qmn'(x)
C       Input :  x  --- Argument of Qmn(x)
C                m  --- Order of Qmn(x)  ( m = 0,1,2,תתת )
C                n  --- Degree of Qmn(x) ( n = 0,1,2,תתת )
C                mm --- Physical dimension of QM and QD
C       Output:  QM(m,n) --- Qmn(x)
C                QD(m,n) --- Qmn'(x)
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (Q,X)
        DIMENSION QM(0:MM,0:N),QD(0:MM,0:N)
        IF (DABS(X).EQ.1.0D0) THEN
           DO 10 I=0,M
           DO 10 J=0,N
              QM(I,J)=1.0D+300
              QD(I,J)=1.0D+300
10         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XS=LS*(1.0D0-X*X)
        XQ=DSQRT(XS)
        Q0=0.5D0*DLOG(DABS((X+1.0D0)/(X-1.0D0)))
        IF (DABS(X).LT.1.0001D0) THEN
           QM(0,0)=Q0
           QM(0,1)=X*Q0-1.0D0
           QM(1,0)=-1.0D0/XQ
           QM(1,1)=-XQ*(Q0+X/(1.0D0-X*X))
           DO 15 I=0,1
           DO 15 J=2,N
              QM(I,J)=((2.0D0*J-1.0D0)*X*QM(I,J-1)
     &               -(J+I-1.0D0)*QM(I,J-2))/(J-I)
15         CONTINUE
           DO 20 J=0,N
           DO 20 I=2,M
              QM(I,J)=-2.0D0*(I-1.0D0)*X/XQ*QM(I-1,J)-LS*
     &                (J+I-1.0D0)*(J-I+2.0D0)*QM(I-2,J)
20         CONTINUE
        ELSE
           IF (DABS(X).GT.1.1D0) THEN
              KM=40+M+N
           ELSE
              KM=(40+M+N)*INT(-1.0-1.8*LOG(X-1.0))
           ENDIF
           QF2=0.0D0
           QF1=1.0D0
           DO 25 K=KM,0,-1
              QF0=((2*K+3.0D0)*X*QF1-(K+2.0D0)*QF2)/(K+1.0D0)
              IF (K.LE.N) QM(0,K)=QF0
              QF2=QF1
25            QF1=QF0
           DO 30 K=0,N
30            QM(0,K)=Q0*QM(0,K)/QF0
           QF2=0.0D0
           QF1=1.0D0
           DO 35 K=KM,0,-1
              QF0=((2*K+3.0D0)*X*QF1-(K+1.0D0)*QF2)/(K+2.0D0)
              IF (K.LE.N) QM(1,K)=QF0
              QF2=QF1
35            QF1=QF0
           Q10=-1.0D0/XQ
           DO 40 K=0,N
40            QM(1,K)=Q10*QM(1,K)/QF0
           DO 45 J=0,N
              Q0=QM(0,J)
              Q1=QM(1,J)
              DO 45 I=0,M-2
                 QF=-2.0D0*(I+1)*X/XQ*Q1+(J-I)*(J+I+1.0D0)*Q0
                 QM(I+2,J)=QF
                 Q0=Q1
                 Q1=QF
45         CONTINUE
        ENDIF
        QD(0,0)=LS/XS
        DO 50 J=1,N
50         QD(0,J)=LS*J*(QM(0,J-1)-X*QM(0,J))/XS
        DO 55 J=0,N
        DO 55 I=1,M
           QD(I,J)=LS*I*X/XS*QM(I,J)+(I+J)*(J-I+1.0D0)/XQ*QM(I-1,J)
55      CONTINUE
        RETURN
        END
