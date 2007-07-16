        PROGRAM MJDZO

C       =============================================================
C       Purpose: This program computes the zeros of Bessel functions
C                Jn(x) and Jn'(x), and arranges them in the order
C                of their values
C       Input :  NT    --- Number of total zeros ( NT ף 1200 )
C       Output:  ZO(L) --- Value of the L-th zero of Jn(x) and
C                          Jn'(x)
C                N(L)  --- n, order of Jn(x) or Jn'(x) associated
C                          with the L-th zero
C                M(L)  --- m, serial number of the zeros of Jn(x)
C                          or Jn'(x) associated with the L-th zero
C                          ( L is the serial number of all the
C                            zeros of Jn(x) and Jn'(x) )
C                P(L)  --- TM or TE, a code for designating the
C                          zeros of Jn(x) or Jn'(x)
C                          In the waveguide applications, the zeros
C                          of Jn(x) correspond to TM modes and those
C                          of Jn'(x) correspond to TE modes.
C       =============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CHARACTER P(1400)*4
        DIMENSION N(1400),M(1400),ZO(1400)
        WRITE(*,*)'NT=?'
        READ(*,*)NT
        WRITE(*,60)NT
        WRITE(*,70)
        CALL JDZO(NT,N,M,P,ZO)
        WRITE(*,*)
        KS=NT/101+1
        DO 55 K0=1,KS
           WRITE(*,*)' Table           Zeros of Bessel',
     &               ' functions Jn(x) and Jn''(x)'
           WRITE(*,*)
           WRITE(*,*)' ----------------------------------',
     &               '----------------------------------'
           DO 50 K=1,50
              J1=100*(K0-1)+K+1
              J2=J1+50
              IF (J1.LE.NT+1.AND.J2.LE.NT+1) THEN
                  WRITE(*,65)J1-1,P(J1),N(J1),M(J1),ZO(J1),
     &                 J2-1,P(J2),N(J2),M(J2),ZO(J2)
              ELSE IF (J1.LE.NT+1.AND.J2.GT.NT+1) THEN
                  WRITE(*,65)J1-1,P(J1),N(J1),M(J1),ZO(J1)
              ENDIF
50         CONTINUE
           WRITE(*,*)' ----------------------------------',
     &               '----------------------------------'
           WRITE(*,75)
55      CONTINUE
60      FORMAT(1X,'Total number of the zeros:',I5)
65      FORMAT(1X,I4,3X,A2,I4,2H -,I2,F14.8,3X,1H|,2X,I4,
     &         3X,A2,I4,2H -,I2,F14.8)
70      FORMAT(15X,'***  Please wait !  The program is running  ***')
75      FORMAT(/)
        END


        SUBROUTINE JDZO(NT,N,M,P,ZO)

C       ===========================================================
C       Purpose: Compute the zeros of Bessel functions Jn(x) and
C                Jn'(x), and arrange them in the order of their
C                magnitudes
C       Input :  NT    --- Number of total zeros ( NT ף 1200 )
C       Output:  ZO(L) --- Value of the L-th zero of Jn(x)
C                          and Jn'(x)
C                N(L)  --- n, order of Jn(x) or Jn'(x) associated
C                          with the L-th zero
C                M(L)  --- m, serial number of the zeros of Jn(x)
C                          or Jn'(x) associated with the L-th zero
C                          ( L is the serial number of all the
C                            zeros of Jn(x) and Jn'(x) )
C                P(L)  --- TM or TE, a code for designating the
C                          zeros of Jn(x)  or Jn'(x).
C                          In the waveguide applications, the zeros
C                          of Jn(x) correspond to TM modes and
C                          those of Jn'(x) correspond to TE modes
C       Routine called:    BJNDD for computing Jn(x), Jn'(x) and
C                          Jn''(x)
C       =============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CHARACTER P(1400)*4,P1(70)*4
        DIMENSION N(1400),M(1400),ZO(1400),N1(70),M1(70),
     &            ZOC(70),BJ(101),DJ(101),FJ(101)
        IF (NT.LT.600) THEN
           XM=-1.0+2.248485*NT**0.5-.0159382*NT+3.208775E-4
     &        *NT**1.5
           NM=INT(14.5+.05875*NT)
           MM=INT(.02*NT)+6
        ELSE
           XM=5.0+1.445389*NT**.5+.01889876*NT-2.147763E-4
     &        *NT**1.5
           NM=INT(27.8+.0327*NT)
           MM=INT(.01088*NT)+10
        ENDIF
        L0=0
        DO 45 I=1,NM
           X1=.407658+.4795504*(I-1)**.5+.983618*(I-1)
           X2=1.99535+.8333883*(I-1)**.5+.984584*(I-1)
           L1=0
           DO 30 J=1,MM
              IF (I.EQ.1.AND.J.EQ.1) GO TO 15
              X=X1
10            CALL BJNDD(I,X,BJ,DJ,FJ)
              X0=X
              X=X-DJ(I)/FJ(I)
              IF (X1.GT.XM) GO TO 20
              IF (DABS(X-X0).GT.1.0D-10) GO TO 10
15            L1=L1+1
              N1(L1)=I-1
              M1(L1)=J
              IF (I.EQ.1) M1(L1)=J-1
              P1(L1)='TE'
              ZOC(L1)=X
              IF (I.LE.15) THEN
                 X1=X+3.057+.0122*(I-1)+(1.555+.41575*(I-1))/(J+1)**2
              ELSE
                 X1=X+2.918+.01924*(I-1)+(6.26+.13205*(I-1))/(J+1)**2
              ENDIF
20            X=X2
25            CALL BJNDD(I,X,BJ,DJ,FJ)
              X0=X
              X=X-BJ(I)/DJ(I)
              IF (X.GT.XM) GO TO 30
              IF (DABS(X-X0).GT.1.0D-10) GO TO 25
              L1=L1+1
              N1(L1)=I-1
              M1(L1)=J
              P1(L1)='TM'
              ZOC(L1)=X
              IF (I.LE.15) THEN
                 X2=X+3.11+.0138*(I-1)+(.04832+.2804*(I-1))/(J+1)**2
              ELSE
                 X2=X+3.001+.0105*(I-1)+(11.52+.48525*(I-1))/(J+3)**2
              ENDIF
30         CONTINUE
           L=L0+L1
           L2=L
35         IF (L0.EQ.0) THEN
              DO 40 K=1,L
                 ZO(K)=ZOC(K)
                 N(K)=N1(K)
                 M(K)=M1(K)
40               P(K)=P1(K)
              L1=0
           ELSE IF (L0.NE.0) THEN
              IF (ZO(L0).GE.ZOC(L1)) THEN
                 ZO(L0+L1)=ZO(L0)
                 N(L0+L1)=N(L0)
                 M(L0+L1)=M(L0)
                 P(L0+L1)=P(L0)
                 L0=L0-1
              ELSE
                 ZO(L0+L1)=ZOC(L1)
                 N(L0+L1)=N1(L1)
                 M(L0+L1)=M1(L1)
                 P(L0+L1)=P1(L1)
                 L1=L1-1
              ENDIF
           ENDIF
           IF (L1.NE.0) GO TO 35
45         L0=L2
        RETURN
        END


        SUBROUTINE BJNDD(N,X,BJ,DJ,FJ)

C       =====================================================
C       Purpose: Compute Bessel functions Jn(x) and their
C                first and second derivatives ( n= 0,1,תתת )
C       Input:   x ---  Argument of Jn(x)  ( x ע 0 )
C                n ---  Order of Jn(x)
C       Output:  BJ(n+1) ---  Jn(x)
C                DJ(n+1) ---  Jn'(x)
C                FJ(n+1) ---  Jn"(x)
C       =====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(101),DJ(101),FJ(101)
        DO 10 NT=1,900
           MT=INT(0.5*LOG10(6.28*NT)-NT*LOG10(1.36*DABS(X)/NT))
           IF (MT.GT.20) GO TO 15
10      CONTINUE
15      M=NT
        BS=0.0D0
        F0=0.0D0
        F1=1.0D-35
        DO 20 K=M,0,-1
           F=2.0D0*(K+1.0D0)*F1/X-F0
           IF (K.LE.N) BJ(K+1)=F
           IF (K.EQ.2*INT(K/2)) BS=BS+2.0D0*F
           F0=F1
20         F1=F
        DO 25 K=0,N
25         BJ(K+1)=BJ(K+1)/(BS-F)
        DJ(1)=-BJ(2)
        FJ(1)=-1.0D0*BJ(1)-DJ(1)/X
        DO 30 K=1,N
           DJ(K+1)=BJ(K)-K*BJ(K+1)/X
30         FJ(K+1)=(K*K/(X*X)-1.0D0)*BJ(K+1)-DJ(K+1)/X
        RETURN
        END
