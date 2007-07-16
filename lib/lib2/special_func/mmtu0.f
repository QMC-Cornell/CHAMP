        PROGRAM MMTU0

C       ============================================================
C       Purpose: This program computes Mathieu functions cem(x,q),
C                sem(x,q) and their derivatives using subroutine
C                MTU0 ( q ò 0 )
C       Input :  KF  --- Function code
C                        KF=1 for computing cem(x,q) and cem'(x,q)
C                        KF=2 for computing sem(x,q) and sem'(x,q)
C                m   --- Order of Mathieu functions
C                q   --- Parameter of Mathieu functions
C                x   --- Argument of Mathieu functions (in degrees)
C       Output:  CSF --- cem(x,q) or sem(x,q)
C                CSD --- cem'x,q) or sem'x,q)
C       Example: x = 40
C           m     q    cem(x,q)   cem'(x,q)    sem(x,q)  sem'(x,q)
C          --------------------------------------------------------
C           0    5.0   .3025683    .9470247
C           1    5.0   .7669652   1.2873097    .2988052   .9606824
C           2    5.0   .9102723   -.3463855    .7549264  1.4743128
C           5    5.0  -.9810931   -.6328576    .1694850 -4.8676455
C           0   25.0   .0515371    .3823737
C           1   25.0   .2074402   1.2646301    .0515365   .3823777
C           2   25.0  -.5297051  -2.4292679    .2074275  1.2646996
C           5   25.0   .7507159  -3.9047012   1.1881232   .3258081
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter KF, m, q and x (in degrees)'
        READ(*,*)KF,M,Q,X
        WRITE(*,10)KF,M,Q,X
        WRITE(*,*)
        IF (KF.EQ.1) WRITE(*,*)' x(degs.)    cem(x,q)       cem''(x,q)'
        IF (KF.EQ.2) WRITE(*,*)' x(degs.)    sem(x,q)       sem''(x,q)'
        WRITE(*,*)' --------------------------------------'
        CALL MTU0(KF,M,Q,X,CSF,CSD)
        WRITE(*,20)X,CSF,CSD
10      FORMAT(1X,4HKF =,I2,',  ',3Hm =,I2,',  ',
     &         3Hq =,F5.1,',  ',3Hx =,F5.1)
20      FORMAT(2X,F5.1,2F16.9)
        END


        SUBROUTINE MTU0(KF,M,Q,X,CSF,CSD)

C       ===============================================================
C       Purpose: Compute Mathieu functions cem(x,q) and sem(x,q)
C                and their derivatives ( q ò 0 )
C       Input :  KF  --- Function code
C                        KF=1 for computing cem(x,q) and cem'(x,q)
C                        KF=2 for computing sem(x,q) and sem'(x,q)
C                m   --- Order of Mathieu functions
C                q   --- Parameter of Mathieu functions
C                x   --- Argument of Mathieu functions (in degrees)
C       Output:  CSF --- cem(x,q) or sem(x,q)
C                CSD --- cem'x,q) or sem'x,q)
C       Routines called:
C            (1) CVA2 for computing the characteristic values
C            (2) FCOEF for computing the expansion coefficients
C       ===============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION FG(251)
        EPS=1.0D-14
        IF (KF.EQ.1.AND.M.EQ.2*INT(M/2)) KD=1
        IF (KF.EQ.1.AND.M.NE.2*INT(M/2)) KD=2
        IF (KF.EQ.2.AND.M.NE.2*INT(M/2)) KD=3
        IF (KF.EQ.2.AND.M.EQ.2*INT(M/2)) KD=4
        CALL CVA2(KD,M,Q,A)
        IF (Q.LE.1.0D0) THEN
           QM=7.5+56.1*SQRT(Q)-134.7*Q+90.7*SQRT(Q)*Q
        ELSE
           QM=17.0+3.1*SQRT(Q)-.126*Q+.0037*SQRT(Q)*Q
        ENDIF
        KM=INT(QM+0.5*M)
        CALL FCOEF(KD,M,Q,A,FG)
        IC=INT(M/2)+1
        RD=1.74532925199433D-2
        XR=X*RD
        CSF=0.0D0
        DO 10 K=1,KM
           IF (KD.EQ.1) THEN
              CSF=CSF+FG(K)*DCOS((2*K-2)*XR)
           ELSE IF (KD.EQ.2) THEN
              CSF=CSF+FG(K)*DCOS((2*K-1)*XR)
           ELSE IF (KD.EQ.3) THEN
              CSF=CSF+FG(K)*DSIN((2*K-1)*XR)
           ELSE IF (KD.EQ.4) THEN
              CSF=CSF+FG(K)*DSIN(2*K*XR)
           ENDIF
           IF (K.GE.IC.AND.DABS(FG(K)).LT.DABS(CSF)*EPS) GO TO 15
10         CONTINUE
15      CSD=0.0D0
        DO 20 K=1,KM
           IF (KD.EQ.1) THEN
              CSD=CSD-(2*K-2)*FG(K)*DSIN((2*K-2)*XR)
           ELSE IF (KD.EQ.2) THEN
              CSD=CSD-(2*K-1)*FG(K)*DSIN((2*K-1)*XR)
           ELSE IF (KD.EQ.3) THEN
              CSD=CSD+(2*K-1)*FG(K)*DCOS((2*K-1)*XR)
           ELSE IF (KD.EQ.4) THEN
              CSD=CSD+2.0D0*K*FG(K)*DCOS(2*K*XR)
           ENDIF
           IF (K.GE.IC.AND.DABS(FG(K)).LT.DABS(CSD)*EPS) GO TO 25
20         CONTINUE
25      RETURN
        END


        SUBROUTINE FCOEF(KD,M,Q,A,FC)

C       =====================================================
C       Purpose: Compute expansion coefficients for Mathieu
C                functions and modified Mathieu functions
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C                A  --- Characteristic value of Mathieu
C                       functions for given m and q
C       Output:  FC(k) --- Expansion coefficients of Mathieu
C                       functions ( k= 1,2,...,KM )
C                       FC(1),FC(2),FC(3),... correspond to
C                       A0,A2,A4,... for KD=1 case, A1,A3,
C                       A5,... for KD=2 case, B1,B3,B5,...
C                       for KD=3 case and B2,B4,B6,... for
C                       KD=4 case
C       =====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION FC(251)
        IF (Q.LE.1.0D0) THEN
           QM=7.5+56.1*SQRT(Q)-134.7*Q+90.7*SQRT(Q)*Q
        ELSE
           QM=17.0+3.1*SQRT(Q)-.126*Q+.0037*SQRT(Q)*Q
        ENDIF
        KM=INT(QM+0.5*M)
        IF (Q.EQ.0.0D0) THEN
           DO 10 K=1,KM
10            FC(K)=0.0D0
           IF (KD.EQ.1) THEN
              FC((M+2)/2)=1.0D0
              IF (M.EQ.0) FC(1)=1.0D0/DSQRT(2.0D0)
           ELSE IF (KD.EQ.4) THEN
              FC(M/2)=1.0D0
           ELSE
              FC((M+1)/2)=1.0D0
           ENDIF
           RETURN
        ENDIF
        KB=0
        S=0.0D0
        F=1.0D-100
        U=0.0D0
        FC(KM)=0.0D0
        IF (KD.EQ.1) THEN
           DO 25 K=KM,3,-1
              V=U
              U=F
              F=(A-4.0D0*K*K)*U/Q-V
              IF (DABS(F).LT.DABS(FC(K+1))) THEN
                 KB=K
                 FC(1)=1.0D-100
                 SP=0.0D0
                 F3=FC(K+1)
                 FC(2)=A/Q*FC(1)
                 FC(3)=(A-4.0D0)*FC(2)/Q-2.0D0*FC(1)
                 U=FC(2)
                 F1=FC(3)
                 DO 15 I=3,KB
                    V=U
                    U=F1
                    F1=(A-4.0D0*(I-1.0D0)**2)*U/Q-V
                    FC(I+1)=F1
                    IF (I.EQ.KB) F2=F1
                    IF (I.NE.KB) SP=SP+F1*F1
15               CONTINUE
                 SP=SP+2.0D0*FC(1)**2+FC(2)**2+FC(3)**2
                 SS=S+SP*(F3/F2)**2
                 S0=DSQRT(1.0D0/SS)
                 DO 20 J=1,KM
                    IF (J.LE.KB+1) THEN
                       FC(J)=S0*FC(J)*F3/F2
                    ELSE
                       FC(J)=S0*FC(J)
                    ENDIF
20               CONTINUE
                 GO TO 85
              ELSE
                 FC(K)=F
                 S=S+F*F
              ENDIF
25         CONTINUE
           FC(2)=Q*FC(3)/(A-4.0D0-2.0D0*Q*Q/A)
           FC(1)=Q/A*FC(2)
           S=S+2.0D0*FC(1)**2+FC(2)**2
           S0=DSQRT(1.0D0/S)
           DO 30 K=1,KM
30            FC(K)=S0*FC(K)
        ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
           DO 35 K=KM,3,-1
              V=U
              U=F
              F=(A-(2.0D0*K-1)**2)*U/Q-V
              IF (DABS(F).GE.DABS(FC(K))) THEN
                 FC(K-1)=F
                 S=S+F*F
              ELSE
                 KB=K
                 F3=FC(K)
                 GO TO 45
              ENDIF
35         CONTINUE
           FC(1)=Q/(A-1.0D0-(-1)**KD*Q)*FC(2)
           S=S+FC(1)*FC(1)
           S0=DSQRT(1.0D0/S)
           DO 40 K=1,KM
40            FC(K)=S0*FC(K)
           GO TO 85
45         FC(1)=1.0D-100
           FC(2)=(A-1.0D0-(-1)**KD*Q)/Q*FC(1)
           SP=0.0D0
           U=FC(1)
           F1=FC(2)
           DO 50 I=2,KB-1
              V=U
              U=F1
              F1=(A-(2.0D0*I-1.0D0)**2)*U/Q-V
              IF (I.NE.KB-1) THEN
                 FC(I+1)=F1
                 SP=SP+F1*F1
              ELSE
                 F2=F1
              ENDIF
50         CONTINUE
           SP=SP+FC(1)**2+FC(2)**2
           SS=S+SP*(F3/F2)**2
           S0=1.0D0/DSQRT(SS)
           DO 55 J=1,KM
              IF (J.LT.KB) FC(J)=S0*FC(J)*F3/F2
              IF (J.GE.KB) FC(J)=S0*FC(J)
55         CONTINUE
        ELSE IF (KD.EQ.4) THEN
           DO 60 K=KM,3,-1
              V=U
              U=F
              F=(A-4.0D0*K*K)*U/Q-V
              IF (DABS(F).GE.DABS(FC(K))) THEN
                 FC(K-1)=F
                 S=S+F*F
              ELSE
                 KB=K
                 F3=FC(K)
                 GO TO 70
              ENDIF
60         CONTINUE
           FC(1)=Q/(A-4.0D0)*FC(2)
           S=S+FC(1)*FC(1)
           S0=DSQRT(1.0D0/S)
           DO 65 K=1,KM
65            FC(K)=S0*FC(K)
           GO TO 85
70         FC(1)=1.0D-100
           FC(2)=(A-4.0D0)/Q*FC(1)
           SP=0.0D0
           U=FC(1)
           F1=FC(2)
           DO 75 I=2,KB-1
              V=U
              U=F1
              F1=(A-4.0D0*I*I)*U/Q-V
              IF (I.NE.KB-1) THEN
                 FC(I+1)=F1
                 SP=SP+F1*F1
              ELSE
                 F2=F1
              ENDIF
75         CONTINUE
           SP=SP+FC(1)**2+FC(2)**2
           SS=S+SP*(F3/F2)**2
           S0=1.0D0/DSQRT(SS)
           DO 80 J=1,KM
              IF (J.LT.KB) FC(J)=S0*FC(J)*F3/F2
              IF (J.GE.KB) FC(J)=S0*FC(J)
80         CONTINUE
        ENDIF
85      IF (FC(1).LT.0.0D0) THEN
           DO 90 J=1,KM
90            FC(J)=-FC(J)
        ENDIF
        RETURN
        END


        SUBROUTINE CVA2(KD,M,Q,A)

C       ======================================================
C       Purpose: Calculate a specific characteristic value of
C                Mathieu functions
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C       Output:  A  --- Characteristic value
C       Routines called:
C             (1) REFINE for finding accurate characteristic
C                 values using an iteration method
C             (2) CV0 for finding initial characteristic
C                 values using polynomial approximation
C             (3) CVQM for computing initial characteristic
C                 values for q ó 3*m
C             (3) CVQL for computing initial characteristic
C                 values for q ò m*m
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (M.LE.12.OR.Q.LE.3.0*M.OR.Q.GT.M*M) THEN
            CALL CV0(KD,M,Q,A)
            IF (Q.NE.0.0D0) CALL REFINE(KD,M,Q,A,1)
        ELSE
           NDIV=10
           DELTA=(M-3.0)*M/NDIV
           IF ((Q-3.0*M).LE.(M*M-Q)) THEN
5             NN=INT((Q-3.0*M)/DELTA)+1
              DELTA=(Q-3.0*M)/NN
              Q1=2.0*M
              CALL CVQM(M,Q1,A1)
              Q2=3.0*M
              CALL CVQM(M,Q2,A2)
              QQ=3.0*M
              DO 10 I=1,NN
                 QQ=QQ+DELTA
                 A=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1)
                 IFLAG=1
                 IF (I.EQ.NN) IFLAG=-1
                 CALL REFINE(KD,M,QQ,A,IFLAG)
                 Q1=Q2
                 Q2=QQ
                 A1=A2
                 A2=A
10            CONTINUE
              IF (IFLAG.EQ.-10) THEN
                 NDIV=NDIV*2
                 DELTA=(M-3.0)*M/NDIV
                 GO TO 5
              ENDIF
           ELSE
15            NN=INT((M*M-Q)/DELTA)+1
              DELTA=(M*M-Q)/NN
              Q1=M*(M-1.0)
              CALL CVQL(KD,M,Q1,A1)
              Q2=M*M
              CALL CVQL(KD,M,Q2,A2)
              QQ=M*M
              DO 20 I=1,NN
                 QQ=QQ-DELTA
                 A=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1)
                 IFLAG=1
                 IF (I.EQ.NN) IFLAG=-1
                 CALL REFINE(KD,M,QQ,A,IFLAG)
                 Q1=Q2
                 Q2=QQ
                 A1=A2
                 A2=A
20            CONTINUE
              IF (IFLAG.EQ.-10) THEN
                 NDIV=NDIV*2
                 DELTA=(M-3.0)*M/NDIV
                 GO TO 15
              ENDIF
           ENDIF
        ENDIF
        RETURN
        END


        SUBROUTINE REFINE(KD,M,Q,A,IFLAG)

C       =====================================================
C       Purpose: calculate the accurate characteristic value
C                by the secant method
C       Input :  m --- Order of Mathieu functions
C                q --- Parameter of Mathieu functions
C                A --- Initial characteristic value
C       Output:  A --- Refineed characteristic value
C       Routine called:  CVF for computing the value of F for
C                        characteristic equation
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-14
        MJ=10+M
        CA=A
        DELTA=0.0D0
        X0=A
        CALL CVF(KD,M,Q,X0,MJ,F0)
        X1=1.002*A
        CALL CVF(KD,M,Q,X1,MJ,F1)
5       DO 10 IT=1,100
           MJ=MJ+1
           X=X1-(X1-X0)/(1.0D0-F0/F1)
           CALL CVF(KD,M,Q,X,MJ,F)
           IF (ABS(1.0-X1/X).LT.EPS.OR.F.EQ.0.0) GO TO 15
           X0=X1
           F0=F1
           X1=X
10         F1=F
15      A=X
        IF (DELTA.GT.0.05) THEN
           A=CA
           IF (IFLAG.LT.0) THEN
              IFLAG=-10
           ENDIF
           RETURN
        ENDIF
        IF (ABS((A-CA)/CA).GT.0.05) THEN
           X0=CA
           DELTA=DELTA+0.005D0
           CALL CVF(KD,M,Q,X0,MJ,F0)
           X1=(1.0D0+DELTA)*CA
           CALL CVF(KD,M,Q,X1,MJ,F1)
           GO TO 5
        ENDIF
        RETURN
        END


        SUBROUTINE CVF(KD,M,Q,A,MJ,F)

C       ======================================================
C       Purpose: Compute the value of F for characteristic
C                equation of Mathieu functions
C       Input :  m --- Order of Mathieu functions
C                q --- Parameter of Mathieu functions
C                A --- Characteristic value
C       Output:  F --- Value of F for characteristic equation
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        B=A
        IC=INT(M/2)
        L=0
        L0=0
        J0=2
        JF=IC
        IF (KD.EQ.1) L0=2
        IF (KD.EQ.1) J0=3
        IF (KD.EQ.2.OR.KD.EQ.3) L=1
        IF (KD.EQ.4) JF=IC-1
        T1=0.0D0
        DO 10 J=MJ,IC+1,-1
10         T1=-Q*Q/((2.0D0*J+L)**2-B+T1)
        IF (M.LE.2) THEN
           T2=0.0D0
           IF (KD.EQ.1.AND.M.EQ.0) T1=T1+T1
           IF (KD.EQ.1.AND.M.EQ.2) T1=-2.0*Q*Q/(4.0-B+T1)-4.0
           IF (KD.EQ.2.AND.M.EQ.1) T1=T1+Q
           IF (KD.EQ.3.AND.M.EQ.1) T1=T1-Q
        ELSE
           IF (KD.EQ.1) T0=4.0D0-B+2.0D0*Q*Q/B
           IF (KD.EQ.2) T0=1.0D0-B+Q
           IF (KD.EQ.3) T0=1.0D0-B-Q
           IF (KD.EQ.4) T0=4.0D0-B
           T2=-Q*Q/T0
           DO 15 J=J0,JF
15            T2=-Q*Q/((2.0D0*J-L-L0)**2-B+T2)
        ENDIF
        F=(2.0D0*IC+L)**2+T1+T2-B
        RETURN
        END


        SUBROUTINE CV0(KD,M,Q,A0)

C       =====================================================
C       Purpose: Compute the initial characteristic value of
C                Mathieu functions for m ó 12  or q ó 300 or
C                q ò m*m
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C       Output:  A0 --- Characteristic value
C       Routines called:
C             (1) CVQM for computing initial characteristic
C                 value for q ó 3*m
C             (2) CVQL for computing initial characteristic
C                 value for q ò m*m
C       ====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        Q2=Q*Q
        IF (M.EQ.0) THEN
           IF (Q.LE.1.0) THEN
              A0=(((.0036392*Q2-.0125868)*Q2+.0546875)*Q2-.5)*Q2
           ELSE IF (Q.LE.10.0) THEN
              A0=((3.999267D-3*Q-9.638957D-2)*Q-.88297)*Q
     &           +.5542818
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.1) THEN
           IF (Q.LE.1.0.AND.KD.EQ.2) THEN
              A0=(((-6.51E-4*Q-.015625)*Q-.125)*Q+1.0)*Q+1.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
              A0=(((-6.51E-4*Q+.015625)*Q-.125)*Q-1.0)*Q+1.0
           ELSE IF (Q.LE.10.0.AND. KD.EQ.2) THEN
              A0=(((-4.94603D-4*Q+1.92917D-2)*Q-.3089229)
     &           *Q+1.33372)*Q+.811752
           ELSE IF (Q.LE.10.0.AND.KD.EQ.3) THEN
              A0=((1.971096D-3*Q-5.482465D-2)*Q-1.152218)
     &           *Q+1.10427
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.2) THEN
           IF (Q.LE.1.0.AND.KD.EQ.1) THEN
              A0=(((-.0036391*Q2+.0125888)*Q2-.0551939)*Q2
     &           +.416667)*Q2+4.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.4) THEN
              A0=(.0003617*Q2-.0833333)*Q2+4.0
           ELSE IF (Q.LE.15.AND.KD.EQ.1) THEN
              A0=(((3.200972D-4*Q-8.667445D-3)*Q
     &           -1.829032D-4)*Q+.9919999)*Q+3.3290504
           ELSE IF (Q.LE.10.0.AND.KD.EQ.4) THEN
              A0=((2.38446D-3*Q-.08725329)*Q-4.732542D-3)
     &           *Q+4.00909
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.3) THEN
           IF (Q.LE.1.0.AND.KD.EQ.2) THEN
              A0=((6.348E-4*Q+.015625)*Q+.0625)*Q2+9.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
              A0=((6.348E-4*Q-.015625)*Q+.0625)*Q2+9.0
           ELSE IF (Q.LE.20.0.AND.KD.EQ.2) THEN
              A0=(((3.035731D-4*Q-1.453021D-2)*Q
     &           +.19069602)*Q-.1039356)*Q+8.9449274
           ELSE IF (Q.LE.15.0.AND.KD.EQ.3) THEN
              A0=((9.369364D-5*Q-.03569325)*Q+.2689874)*Q
     &           +8.771735
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.4) THEN
           IF (Q.LE.1.0.AND.KD.EQ.1) THEN
              A0=((-2.1E-6*Q2+5.012E-4)*Q2+.0333333)*Q2+16.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.4) THEN
              A0=((3.7E-6*Q2-3.669E-4)*Q2+.0333333)*Q2+16.0
           ELSE IF (Q.LE.25.0.AND.KD.EQ.1) THEN
              A0=(((1.076676D-4*Q-7.9684875D-3)*Q
     &           +.17344854)*Q-.5924058)*Q+16.620847
           ELSE IF (Q.LE.20.0.AND.KD.EQ.4) THEN
              A0=((-7.08719D-4*Q+3.8216144D-3)*Q
     &           +.1907493)*Q+15.744
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.5) THEN
           IF (Q.LE.1.0.AND.KD.EQ.2) THEN
              A0=((6.8E-6*Q+1.42E-5)*Q2+.0208333)*Q2+25.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
              A0=((-6.8E-6*Q+1.42E-5)*Q2+.0208333)*Q2+25.0
           ELSE IF (Q.LE.35.0.AND.KD.EQ.2) THEN
              A0=(((2.238231D-5*Q-2.983416D-3)*Q
     &           +.10706975)*Q-.600205)*Q+25.93515
           ELSE IF (Q.LE.25.0.AND.KD.EQ.3) THEN
              A0=((-7.425364D-4*Q+2.18225D-2)*Q
     &           +4.16399D-2)*Q+24.897
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.6) THEN
           IF (Q.LE.1.0) THEN
              A0=(.4D-6*Q2+.0142857)*Q2+36.0
           ELSE IF (Q.LE.40.0.AND.KD.EQ.1) THEN
              A0=(((-1.66846D-5*Q+4.80263D-4)*Q
     &           +2.53998D-2)*Q-.181233)*Q+36.423
           ELSE IF (Q.LE.35.0.AND.KD.EQ.4) THEN
              A0=((-4.57146D-4*Q+2.16609D-2)*Q-2.349616D-2)*Q
     &           +35.99251
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.7) THEN
           IF (Q.LE.10.0) THEN
              CALL CVQM(M,Q,A0)
           ELSE IF (Q.LE.50.0.AND.KD.EQ.2) THEN
              A0=(((-1.411114D-5*Q+9.730514D-4)*Q
     &           -3.097887D-3)*Q+3.533597D-2)*Q+49.0547
           ELSE IF (Q.LE.40.0.AND.KD.EQ.3) THEN
              A0=((-3.043872D-4*Q+2.05511D-2)*Q
     &           -9.16292D-2)*Q+49.19035
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.GE.8) THEN
           IF (Q.LE.3.*M) THEN
              CALL CVQM(M,Q,A0)
           ELSE IF (Q.GT.M*M) THEN
              CALL CVQL(KD,M,Q,A0)
           ELSE
              IF (M.EQ.8.AND.KD.EQ.1) THEN
                 A0=(((8.634308D-6*Q-2.100289D-3)*Q+.169072)*Q
     &              -4.64336)*Q+109.4211
              ELSE IF (M.EQ.8.AND.KD.EQ.4) THEN
                 A0=((-6.7842D-5*Q+2.2057D-3)*Q+.48296)*Q+56.59
              ELSE IF (M.EQ.9.AND.KD.EQ.2) THEN
                 A0=(((2.906435D-6*Q-1.019893D-3)*Q+.1101965)*Q
     &              -3.821851)*Q+127.6098
              ELSE IF (M.EQ.9.AND.KD.EQ.3) THEN
                 A0=((-9.577289D-5*Q+.01043839)*Q+.06588934)*Q
     &              +78.0198
              ELSE IF (M.EQ.10.AND.KD.EQ.1) THEN
                 A0=(((5.44927D-7*Q-3.926119D-4)*Q+.0612099)*Q
     &              -2.600805)*Q+138.1923
              ELSE IF (M.EQ.10.AND.KD.EQ.4) THEN
                 A0=((-7.660143D-5*Q+.01132506)*Q-.09746023)*Q
     &              +99.29494
              ELSE IF (M.EQ.11.AND.KD.EQ.2) THEN
                 A0=(((-5.67615D-7*Q+7.152722D-6)*Q+.01920291)*Q
     &              -1.081583)*Q+140.88
              ELSE IF (M.EQ.11.AND.KD.EQ.3) THEN
                 A0=((-6.310551D-5*Q+.0119247)*Q-.2681195)*Q
     &              +123.667
              ELSE IF (M.EQ.12.AND.KD.EQ.1) THEN
                 A0=(((-2.38351D-7*Q-2.90139D-5)*Q+.02023088)*Q
     &              -1.289)*Q+171.2723
              ELSE IF (M.EQ.12.AND.KD.EQ.4) THEN
                 A0=(((3.08902D-7*Q-1.577869D-4)*Q+.0247911)*Q
     &              -1.05454)*Q+161.471
              ENDIF
           ENDIF
        ENDIF
        RETURN
        END


        SUBROUTINE CVQL(KD,M,Q,A0)

C       ========================================================
C       Purpose: Compute the characteristic value of Mathieu
C                functions  for q ò 3m
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C       Output:  A0 --- Initial characteristic value
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (KD.EQ.1.OR.KD.EQ.2) W=2.0D0*M+1.0D0
        IF (KD.EQ.3.OR.KD.EQ.4) W=2.0D0*M-1.0D0
        W2=W*W
        W3=W*W2
        W4=W2*W2
        W6=W2*W4
        D1=5.0+34.0/W2+9.0/W4
        D2=(33.0+410.0/W2+405.0/W4)/W
        D3=(63.0+1260.0/W2+2943.0/W4+486.0/W6)/W2
        D4=(527.0+15617.0/W2+69001.0/W4+41607.0/W6)/W3
        C1=128.0
        P2=Q/W4
        P1=DSQRT(P2)
        CV1=-2.0*Q+2.0*W*DSQRT(Q)-(W2+1.0)/8.0
        CV2=(W+3.0/W)+D1/(32.0*P1)+D2/(8.0*C1*P2)
        CV2=CV2+D3/(64.0*C1*P1*P2)+D4/(16.0*C1*C1*P2*P2)
        A0=CV1-CV2/(C1*P1)
        RETURN
        END


        SUBROUTINE CVQM(M,Q,A0)

C       =====================================================
C       Purpose: Compute the characteristic value of Mathieu
C                functions for q ó m*m
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C       Output:  A0 --- Initial characteristic value
C       =====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        HM1=.5*Q/(M*M-1.0)
        HM3=.25*HM1**3/(M*M-4.0)
        HM5=HM1*HM3*Q/((M*M-1.0)*(M*M-9.0))
        A0=M*M+Q*(HM1+(5.0*M*M+7.0)*HM3
     &     +(9.0*M**4+58.0*M*M+29.0)*HM5)
        RETURN
        END


