        PROGRAM MCVA1

C       ============================================================
C       Purpose: This program computes a sequence of characteristic
C                values of Mathieu functions using subroutine CVA1
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C       Output:  CV(I) --- Characteristic values; I = 1,2,3,...
C                For KD=1, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of cem for m = 0,2,4,...
C                For KD=2, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of cem for m = 1,3,5,...
C                For KD=3, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of sem for m = 1,3,5,...
C                For KD=4, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of sem for m = 0,2,4,...

C       Example: Mmax = 12,    q = 25.00

C                Characteristic values of Mathieu functions

C                  m            a                  b
C                ------------------------------------------
C                  0      -40.256779547
C                  1      -21.314899691      -40.256778985
C                  2       -3.522164727      -21.314860622
C                  3       12.964079444       -3.520941527
C                  4       27.805240581       12.986489953
C                  5       40.050190986       28.062765899
C                  6       48.975786716       41.801071292
C                  7       57.534689001       55.002957151
C                  8       69.524065166       69.057988351
C                  9       85.076999882       85.023356505
C                 10      103.230204804      103.225680042
C                 11      123.643012376      123.642713667
C                 12      146.207690643      146.207674647
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION CV1(200),CV2(200),CVE(200),CVS(200)
        WRITE(*,*)'Please enter Mmax,q =?'
        READ(*,*)MMAX,Q
        WRITE(*,25)MMAX,Q
        WRITE(*,*)
        CALL CVA1(1,MMAX,Q,CV1)
        CALL CVA1(2,MMAX,Q,CV2)
        DO 10 J=1,MMAX/2+1
           CVE(2*J-1)=CV1(J)
10         CVE(2*J)=CV2(J)
        CALL CVA1(3,MMAX,Q,CV1)
        CALL CVA1(4,MMAX,Q,CV2)
        DO 15 J=1,MMAX/2+1
           CVS(2*J)=CV1(J)
15         CVS(2*J+1)=CV2(J)
        WRITE(*,35)
        WRITE(*,*)
        WRITE(*,*)'  m            a                  b'
        WRITE(*,*)'------------------------------------------'
        DO 20 J=0,MMAX
           IF (J.EQ.0) WRITE(*,30)J,CVE(J+1)
           IF (J.NE.0) WRITE(*,30)J,CVE(J+1),CVS(J+1)
20      CONTINUE
25      FORMAT(3X,6HMmax =,I3,',    ',3Hq =,F6.2)
30      FORMAT(1X,I3,2F19.9)
35      FORMAT(1X,'Characteristic values of Mathieu functions')
        END


        SUBROUTINE CVA1(KD,M,Q,CV)

C       ============================================================
C       Purpose: Compute a sequence of characteristic values of
C                Mathieu functions
C       Input :  M  --- Maximum order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,תתת )
C                       KD=2 for cem(x,q)  ( m = 1,3,5,תתת )
C                       KD=3 for sem(x,q)  ( m = 1,3,5,תתת )
C                       KD=4 for sem(x,q)  ( m = 2,4,6,תתת )
C       Output:  CV(I) --- Characteristic values; I = 1,2,3,...
C                For KD=1, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of cem for m = 0,2,4,...
C                For KD=2, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of cem for m = 1,3,5,...
C                For KD=3, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of sem for m = 1,3,5,...
C                For KD=4, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of sem for m = 0,2,4,...
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(200),H(200),D(500),E(500),F(500),CV(200)
        EPS=1.0D-14
        ICM=INT(M/2)+1
        IF (KD.EQ.4) ICM=M/2
        IF (Q.EQ.0.0D0) THEN
           IF (KD.EQ.1) THEN
              DO 10 IC=1,ICM
10               CV(IC)=4.0D0*(IC-1.0D0)**2
           ELSE IF (KD.NE.4) THEN
              DO 15 IC=1,ICM
15               CV(IC)=(2.0D0*IC-1.0D0)**2
           ELSE
              DO 20 IC=1,ICM
20               CV(IC)=4.0D0*IC*IC
           ENDIF
        ELSE
           NM=INT(10+1.5*M+0.5*Q)
           E(1)=0.0D0
           F(1)=0.0D0
           IF (KD.EQ.1) THEN
              D(1)=0.0D0
              DO 25 I=2,NM
                 D(I)=4.0D0*(I-1.0D0)**2
                 E(I)=Q
25               F(I)=Q*Q
              E(2)=DSQRT(2.0D0)*Q
              F(2)=2.0D0*Q*Q
           ELSE IF (KD.NE.4) THEN
              D(1)=1.0D0+(-1)**KD*Q
              DO 30 I=2,NM
                 D(I)=(2.0D0*I-1.0D0)**2
                 E(I)=Q
30               F(I)=Q*Q
           ELSE
              D(1)=4.0D0
              DO 35 I=2,NM
                 D(I)=4.0D0*I*I
                 E(I)=Q
35               F(I)=Q*Q
           ENDIF
           XA=D(NM)+DABS(E(NM))
           XB=D(NM)-DABS(E(NM))
           NM1=NM-1
           DO 40 I=1,NM1
              T=DABS(E(I))+DABS(E(I+1))
              T1=D(I)+T
              IF (XA.LT.T1) XA=T1
              T1=D(I)-T
              IF (T1.LT.XB) XB=T1
40         CONTINUE
           DO 45 I=1,ICM
              G(I)=XA
45            H(I)=XB
           DO 75 K=1,ICM
              DO 50 K1=K,ICM
                 IF (G(K1).LT.G(K)) THEN
                    G(K)=G(K1)
                    GO TO 55
                 ENDIF
50            CONTINUE
55            IF (K.NE.1.AND.H(K).LT.H(K-1)) H(K)=H(K-1)
60            X1=(G(K)+H(K))/2.0D0
              CV(K)=X1
              IF (DABS((G(K)-H(K))/X1).LT.EPS) GO TO 70
              J=0
              S=1.0D0
              DO 65 I=1,NM
                 IF (S.EQ.0.0D0) S=S+1.0D-30
                 T=F(I)/S
                 S=D(I)-T-X1
                 IF (S.LT.0.0) J=J+1
65            CONTINUE
              IF (J.LT.K) THEN
                 H(K)=X1
              ELSE
                 G(K)=X1
                 IF (J.GE.ICM) THEN
                    G(ICM)=X1
                 ELSE
                    IF (H(J+1).LT.X1) H(J+1)=X1
                    IF (X1.LT.G(J)) G(J)=X1
                 ENDIF
              ENDIF
              GO TO 60
70            CV(K)=X1
75         CONTINUE
        ENDIF
        RETURN
        END
