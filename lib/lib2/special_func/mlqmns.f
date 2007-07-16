        PROGRAM MLQMNS

C       =========================================================
C       Purpose: This program computes the associated Legendre
C                functions Qmn(x) and their derivatives Qmn'(x)
C                for a given order using subroutine LQMNS
C       Input :  x --- Argument of Qmn(x)
C                m --- Order of Qmn(x),  m = 0,1,2,...
C                n --- Degree of Qmn(x), n = 0,1,2,...
C       Output:  QM(n) --- Qmn(x)
C                QD(n) --- Qmn'(x)
C       Examples:
C                m = 1,  N = 5,  x = .5
C                n        Qmn(x)           Qmn'(x)
C               -------------------------------------
C                0    .11547005D+01    .76980036D+00
C                1    .10530633D+01    .23771592D+01
C                2   -.72980606D+00    .51853281D+01
C                3   -.24918526D+01    .10914062D+01
C                4   -.19340866D+01   -.11454786D+02
C                5    .93896830D+00   -.18602587D+02

C                m = 2,  N = 5,  x = 2.5
C                n        Qmn(x)           Qmn'(x)
C               -------------------------------------
C                0    .95238095D+00   -.52607710D+00
C                1    .38095238D+00   -.36281179D+00
C                2    .12485160D+00   -.17134314D+00
C                3    .36835513D-01   -.66284127D-01
C                4    .10181730D-01   -.22703958D-01
C                5    .26919481D-02   -.71662396D-02
C       =========================================================

        IMPLICIT DOUBLE PRECISION (Q,X,Y)
        DIMENSION QM(0:200),QD(0:200)
        WRITE(*,*)'Please enter m, N, and x '
        READ(*,*)M,N,X
        WRITE(*,30)M,N,X
        CALL LQMNS(M,N,X,QM,QD)
        WRITE(*,*)
        WRITE(*,*)'  n        Qmn(x)           Qmn''(x)'
        WRITE(*,*)' -------------------------------------'
        DO 10 J=0,N
        WRITE(*,20)J,QM(J),QD(J)
10      CONTINUE
20      FORMAT(1X,I3,2D17.8)
30      FORMAT(1X,'m =',I2,',  ','n =',I2,',  ','x =',F5.1)
        END


        SUBROUTINE LQMNS(M,N,X,QM,QD)

C       ========================================================
C       Purpose: Compute associated Legendre functions Qmn(x)
C                and Qmn'(x) for a given order
C       Input :  x --- Argument of Qmn(x)
C                m --- Order of Qmn(x),  m = 0,1,2,...
C                n --- Degree of Qmn(x), n = 0,1,2,...
C       Output:  QM(n) --- Qmn(x)
C                QD(n) --- Qmn'(x)
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION QM(0:N),QD(0:N)
        DO 10 K=0,N
           QM(K)=0.0D0
10         QD(K)=0.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 K=0,N
              QM(K)=1.0D+300
15            QD(K)=1.0D+300
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        Q0=0.5D0*DLOG(DABS((X+1.0)/(X-1.0)))
        Q00=Q0
        Q10=-1.0D0/XQ
        Q01=X*Q0-1.0D0
        Q11=-LS*XQ*(Q0+X/(1.0D0-X*X))
        QF0=Q00
        QF1=Q10
        DO 20 K=2,M
           QM0=-2.0D0*(K-1.0)/XQ*X*QF1-LS*(K-1.0)*(2.0-K)*QF0
           QF0=QF1
20         QF1=QM0
        IF (M.EQ.0) QM0=Q00
        IF (M.EQ.1) QM0=Q10
        QM(0)=QM0
        IF (DABS(X).LT.1.0001D0) THEN
           IF (M.EQ.0.AND.N.GT.0) THEN
              QF0=Q00
              QF1=Q01
              DO 25 K=2,N
                 QF2=((2.0*K-1.0D0)*X*QF1-(K-1.0)*QF0)/K
                 QM(K)=QF2
                 QF0=QF1
25               QF1=QF2
           ENDIF
           QG0=Q01
           QG1=Q11
           DO 30 K=2,M
              QM1=-2.0D0*(K-1.0)/XQ*X*QG1-LS*K*(3.0-K)*QG0
              QG0=QG1
30            QG1=QM1
           IF (M.EQ.0) QM1=Q01
           IF (M.EQ.1) QM1=Q11
           QM(1)=QM1
           IF (M.EQ.1.AND.N.GT.1) THEN
              QH0=Q10
              QH1=Q11
              DO 35 K=2,N
                 QH2=((2.0*K-1.0D0)*X*QH1-K*QH0)/(K-1.0)
                 QM(K)=QH2
                 QH0=QH1
35               QH1=QH2
           ELSE IF (M.GE.2) THEN
              QG0=Q00
              QG1=Q01
              QH0=Q10
              QH1=Q11
              DO 45 L=2,N
                 Q0L=((2.0D0*L-1.0D0)*X*QG1-(L-1.0D0)*QG0)/L
                 Q1L=((2.0*L-1.0D0)*X*QH1-L*QH0)/(L-1.0D0)
                 QF0=Q0L
                 QF1=Q1L
                 DO 40 K=2,M
                    QMK=-2.0D0*(K-1.0)/XQ*X*QF1-LS*(K+L-1.0)*
     &                  (L+2.0-K)*QF0
                    QF0=QF1
40                  QF1=QMK
                 QM(L)=QMK
                 QG0=QG1
                 QG1=Q0L
                 QH0=QH1
45               QH1=Q1L
           ENDIF
        ELSE
           IF (DABS(X).GT.1.1) THEN
              KM=40+M+N
           ELSE
              KM=(40+M+N)*INT(-1.0-1.8*LOG(X-1.0))
           ENDIF
           QF2=0.0D0
           QF1=1.0D0
           DO 50 K=KM,0,-1
              QF0=((2.0*K+3.0D0)*X*QF1-(K+2.0-M)*QF2)/(K+M+1.0)
              IF (K.LE.N) QM(K)=QF0
              QF2=QF1
50            QF1=QF0
           DO 55 K=0,N
55            QM(K)=QM(K)*QM0/QF0
        ENDIF
        IF (DABS(X).LT.1.0D0) THEN
           DO 60 K=0,N
60            QM(K)=(-1)**M*QM(K)
        ENDIF
        QD(0)=((1.0D0-M)*QM(1)-X*QM(0))/(X*X-1.0)
        DO 65 K=1,N
65         QD(K)=(K*X*QM(K)-(K+M)*QM(K-1))/(X*X-1.0)
        RETURN
        END
