        PROGRAM MCVA2

C       =============================================================
C       Purpose: This program calculates a specific characteristic
C                value of Mathieu functions using subroutine CVA2
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C       Output:  A  --- Characteristic value
C       Example: q = 25.0, m = 0,1,2,...,12

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
C       =============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter KD, m and q '
        READ(*,*)KD,M,Q
        CALL CVA2(KD,M,Q,A)
        IF (KD.LE.2) THEN
           WRITE(*,20)
        ELSE
           WRITE(*,30)
        ENDIF
        WRITE(*,*)
        WRITE(*,*)'  m         a (or b)'
        WRITE(*,*)'------------------------'
        WRITE(*,10)M,A
10      FORMAT(1X,I3,F18.8)
20      FORMAT(1X,'Characteristic value of even Mathieu function, a')
30      FORMAT(1X,'Characteristic value of odd Mathieu function, b')
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
C                 value using an iteration method
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
