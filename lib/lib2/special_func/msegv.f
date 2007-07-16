        PROGRAM MSEGV

C       ============================================================
C       Purpose: This program computes a sequence of characteristic
C                values of spheroidal prolate and oblate functions
C                using subroutine SEGV
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                KD --- Function code
C                       KD=1 for Prolate; KD=-1 for Oblate
C       Output:  CV --- Characteristic value for given m, n and c
C                EG(L) --- Characteristic value for mode m and n'
C                          ( L = n' - m + 1 )
C       Examples:
C                Prolate: ( KD = 1 , m = 1, n = 5, c = 5.0 )

C                m      n       c        Lambda mn(c)
C              ---------------------------------------
C                1      1      5.0        5.35042230
C                1      2      5.0       14.64295624
C                1      3      5.0       23.39761312
C                1      4      5.0       32.42194359
C                1      5      5.0       42.65818215

C                Oblate: ( KD = -1 , m = 1, n = 5, c = 5.0 )

C                m      n       c      Lambda mn(-ic)
C               --------------------------------------
C                1      1      5.0       -7.49338828
C                1      2      5.0       -7.12783752
C                1      3      5.0        2.75036721
C                1      4      5.0        8.69495925
C                1      5      5.0       18.43931577
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EG(100)
        WRITE(*,*)'Please enter KD, m, n and c '
        READ(*,*)KD,M,N,C
        WRITE(*,15)KD,M,N,C
        WRITE(*,*)
        CALL SEGV(M,N,C,KD,CV,EG)
        IF (KD.EQ.1) THEN
           WRITE(*,*)'  m      n       c       Lambda mn(c)'
        ELSE IF (KD.EQ.-1) THEN
           WRITE(*,*)'  m      n       c      Lambda mn(-ic)'
        ENDIF
        WRITE(*,*)'---------------------------------------'
        DO 10 L=1,N-M+1
           N1=M+L-1
10         WRITE(*,20)M,N1,C,EG(L)
15      FORMAT(1X,'KD =',I2,',  m =',I3,',   n =',I3,',  c =',F5.1)
20      FORMAT(1X,I3,4X,I3,4X,F5.1,F18.8)
        END


        SUBROUTINE SEGV(M,N,C,KD,CV,EG)

C       =========================================================
C       Purpose: Compute the characteristic values of spheroidal
C                wave functions
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                KD --- Function code
C                       KD=1 for Prolate; KD=-1 for Oblate
C       Output:  CV --- Characteristic value for given m, n and c
C                EG(L) --- Characteristic value for mode m and n'
C                          ( L = n' - m + 1 )
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION B(100),H(100),D(300),E(300),F(300),CV0(100),
     &            A(300),G(300),EG(200)
        IF (C.LT.1.0D-10) THEN
           DO 5 I=1,N
5             EG(I)=(I+M)*(I+M-1.0D0)
           GO TO 70
        ENDIF
        ICM=(N-M+2)/2
        NM=10+INT(0.5*(N-M)+C)
        CS=C*C*KD
        DO 60 L=0,1
           DO 10 I=1,NM
              IF (L.EQ.0) K=2*(I-1)
              IF (L.EQ.1) K=2*I-1
              DK0=M+K
              DK1=M+K+1
              DK2=2*(M+K)
              D2K=2*M+K
              A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
              D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &             *(DK2+3.0))*CS
10            G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
           DO 15 K=2,NM
              E(K)=DSQRT(A(K-1)*G(K))
15            F(K)=E(K)*E(K)
           F(1)=0.0D0
           E(1)=0.0D0
           XA=D(NM)+DABS(E(NM))
           XB=D(NM)-DABS(E(NM))
           NM1=NM-1
           DO 20 I=1,NM1
              T=DABS(E(I))+DABS(E(I+1))
              T1=D(I)+T
              IF (XA.LT.T1) XA=T1
              T1=D(I)-T
              IF (T1.LT.XB) XB=T1
20         CONTINUE
           DO 25 I=1,ICM
              B(I)=XA
25            H(I)=XB
           DO 55 K=1,ICM
              DO 30 K1=K,ICM
                 IF (B(K1).LT.B(K)) THEN
                    B(K)=B(K1)
                    GO TO 35
                 ENDIF
30            CONTINUE
35            IF (K.NE.1.AND.H(K).LT.H(K-1)) H(K)=H(K-1)
40            X1=(B(K)+H(K))/2.0D0
              CV0(K)=X1
              IF (DABS((B(K)-H(K))/X1).LT.1.0D-14) GO TO 50
              J=0
              S=1.0D0
              DO 45 I=1,NM
                 IF (S.EQ.0.0D0) S=S+1.0D-30
                 T=F(I)/S
                 S=D(I)-T-X1
                 IF (S.LT.0.0D0) J=J+1
45            CONTINUE
              IF (J.LT.K) THEN
                 H(K)=X1
              ELSE
                 B(K)=X1
                 IF (J.GE.ICM) THEN
                    B(ICM)=X1
                 ELSE
                    IF (H(J+1).LT.X1) H(J+1)=X1
                    IF (X1.LT.B(J)) B(J)=X1
                 ENDIF
              ENDIF
              GO TO 40
50            CV0(K)=X1
              IF (L.EQ.0) EG(2*K-1)=CV0(K)
              IF (L.EQ.1) EG(2*K)=CV0(K)
55         CONTINUE
60      CONTINUE
70      CV=EG(N-M+1)
        RETURN
        END
