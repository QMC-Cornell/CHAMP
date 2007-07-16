        PROGRAM MCHGU

C       =======================================================
C       Purpose: This program computes the confluent
C                hypergeometric function U(a,b,x) using
C                subroutine CHGU
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument  ( x ע 0 )
C       Output:  HU --- U(a,b,x)
C                MD --- Method code
C       Example:
C                a       b       x        U(a,b,x)
C             --------------------------------------
C              -2.5     2.5     5.0     -9.02812446
C              -1.5     2.5     5.0      2.15780560
C               -.5     2.5     5.0      1.76649370
C                .0     2.5     5.0      1.00000000
C                .5     2.5     5.0       .49193496
C               1.5     2.5     5.0       .08944272
C               2.5     2.5     5.0       .01239387

C                a       b       x        U(a,b,x)
C             --------------------------------------
C              -2.5     5.0    10.0     -2.31982196
C              -1.5     5.0    10.0      8.65747115
C               -.5     5.0    10.0      2.37997143
C                .0     5.0    10.0      1.00000000
C                .5     5.0    10.0       .38329536
C               1.5     5.0    10.0       .04582817
C               2.5     5.0    10.0       .00444535
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Pleas enter a, b and x '
        READ(*,*)A,B,X
        WRITE(*,*)'   a       b       x        U(a,b,x)'
        WRITE(*,*)'--------------------------------------'
        CALL CHGU(A,B,X,HU,MD)
        WRITE(*,10)A,B,X,HU
10      FORMAT(1X,F5.1,3X,F5.1,3X,F5.1,E15.8)
        END


        SUBROUTINE CHGU(A,B,X,HU,MD)

C       =======================================================
C       Purpose: Compute the confluent hypergeometric function
C                U(a,b,x)
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument  ( x > 0 )
C       Output:  HU --- U(a,b,x)
C                MD --- Method code
C       Routines called:
C            (1) CHGUS for small x ( MD=1 )
C            (2) CHGUL for large x ( MD=2 )
C            (3) CHGUBI for integer b ( MD=3 )
C            (4) CHGUIT for numerical integration ( MD=4 )
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL IL1,IL2,IL3,BL1,BL2,BL3,BN
        AA=A-B+1.0D0
        IL1=A.EQ.INT(A).AND.A.LE.0.0
        IL2=AA.EQ.INT(AA).AND.AA.LE.0.0
        IL3=ABS(A*(A-B+1.0))/X.LE.2.0
        BL1=X.LE.5.0.OR.(X.LE.10.0.AND.A.LE.2.0)
        BL2=(X.GT.5.0.AND.X.LE.12.5).AND.(A.GE.1.0.AND.B.GE.A+4.0)
        BL3=X.GT.12.5.AND.A.GE.5.0.AND.B.GE.A+5.0
        BN=B.EQ.INT(B).AND.B.NE.0.0
        ID1=-100
        IF (B.NE.INT(B)) THEN
           CALL CHGUS(A,B,X,HU,ID1)
           MD=1
           IF (ID1.GE.6) RETURN
           HU1=HU
        ENDIF
        IF (IL1.OR.IL2.OR.IL3) THEN
           CALL CHGUL(A,B,X,HU,ID)
           MD=2
           IF (ID.GE.6) RETURN
           IF (ID1.GT.ID) THEN
              MD=1
              ID=ID1
              HU=HU1
           ENDIF
        ENDIF
        IF (A.GE.0.0) THEN
           IF (BN.AND.(BL1.OR.BL2.OR.BL3)) THEN
              CALL CHGUBI(A,B,X,HU,ID)
              MD=3
           ELSE
              CALL CHGUIT(A,B,X,HU,ID)
              MD=4
           ENDIF
        ELSE
           IF (B.LE.A) THEN
              A00=A
              B00=B
              A=A-B+1.0D0
              B=2.0D0-B
              CALL CHGUIT(A,B,X,HU,ID)
              HU=X**(1.0D0-B00)*HU
              A=A00
              B=B00
              MD=4
           ELSE IF (BN.AND.(.NOT.IL1)) THEN
              CALL CHGUBI(A,B,X,HU,ID)
              MD=3
           ENDIF
        ENDIF
        IF (ID.LT.6) WRITE(*,*)'No accurate result obtained'
        RETURN
        END


        SUBROUTINE CHGUS(A,B,X,HU,ID)

C       ======================================================
C       Purpose: Compute confluent hypergeometric function
C                U(a,b,x) for small argument x
C       Input  : a  --- Parameter
C                b  --- Parameter ( b <> 0,-1,-2,...)
C                x  --- Argument
C       Output:  HU --- U(a,b,x)
C                ID --- Estimated number of significant digits
C       Routine called: GAMMA for computing gamma function
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        ID=-100
        PI=3.141592653589793D0
        CALL GAMMA(A,GA)
        CALL GAMMA(B,GB)
        XG1=1.0D0+A-B
        CALL GAMMA(XG1,GAB)
        XG2=2.0D0-B
        CALL GAMMA(XG2,GB2)
        HU0=PI/DSIN(PI*B)
        R1=HU0/(GAB*GB)
        R2=HU0*X**(1.0D0-B)/(GA*GB2)
        HU=R1-R2
        HMAX=0.0D0
        HMIN=1.0D+300
        DO 10 J=1,150
           R1=R1*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
           R2=R2*(A-B+J)/(J*(1.0D0-B+J))*X
           HU=HU+R1-R2
           HUA=DABS(HU)
           IF (HUA.GT.HMAX) HMAX=HUA
           IF (HUA.LT.HMIN) HMIN=HUA
           IF (DABS(HU-H0).LT.DABS(HU)*1.0D-15) GO TO 15
10         H0=HU
15      D1=LOG10(HMAX)
        IF (HMIN.NE.0.0) D2=LOG10(HMIN)
        ID=15-ABS(D1-D2)
        RETURN
        END


        SUBROUTINE CHGUL(A,B,X,HU,ID)

C       =======================================================
C       Purpose: Compute the confluent hypergeometric function
C                U(a,b,x) for large argument x
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument
C       Output:  HU --- U(a,b,x)
C                ID --- Estimated number of significant digits
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL IL1,IL2
        ID=-100
        AA=A-B+1.0D0
        IL1=A.EQ.INT(A).AND.A.LE.0.0
        IL2=AA.EQ.INT(AA).AND.AA.LE.0.0
        IF (IL1) NM=ABS(A)
        IF (IL2) NM=ABS(AA)
        IF (IL1.OR.IL2) THEN
           HU=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=-R*(A+K-1.0D0)*(A-B+K)/(K*X)
              HU=HU+R
10         CONTINUE
           HU=X**(-A)*HU
           ID=10
        ELSE
           HU=1.0D0
           R=1.0D0
           DO 15 K=1,25
              R=-R*(A+K-1.0D0)*(A-B+K)/(K*X)
              RA=DABS(R)
              IF (K.GT.5.AND.RA.GE.R0.OR.RA.LT.1.0D-15) GO TO 20
              R0=RA
15            HU=HU+R
20         ID=ABS(LOG10(RA))
           HU=X**(-A)*HU
        ENDIF
        RETURN
        END


        SUBROUTINE CHGUBI(A,B,X,HU,ID)

C       ======================================================
C       Purpose: Compute confluent hypergeometric function
C                U(a,b,x) with integer b ( b = ס1,ס2,... )
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument
C       Output:  HU --- U(a,b,x)
C                ID --- Estimated number of significant digits
C       Routines called:
C            (1) GAMMA for computing gamma function ג(x)
C            (2) PSI for computing psi function
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        ID=-100
        EL=0.5772156649015329D0
        N=ABS(B-1)
        RN1=1.0D0
        RN=1.0D0
        DO 10 J=1,N
           RN=RN*J
           IF (J.EQ.N-1) RN1=RN
10      CONTINUE
        CALL PSI(A,PS)
        CALL GAMMA(A,GA)
        IF (B.GT.0.0) THEN
           A0=A
           A1=A-N
           A2=A1
           CALL GAMMA(A1,GA1)
           UA=(-1)**(N-1)/(RN*GA1)
           UB=RN1/GA*X**(-N)
        ELSE
           A0=A+N
           A1=A0
           A2=A
           CALL GAMMA(A1,GA1)
           UA=(-1)**(N-1)/(RN*GA)*X**N
           UB=RN1/GA1
        ENDIF
        HM1=1.0D0
        R=1.0D0
        HMAX=0.0D0
        HMIN=1.0D+300
        DO 15 K=1,150
           R=R*(A0+K-1.0D0)*X/((N+K)*K)
           HM1=HM1+R
           HU1=DABS(HM1)
           IF (HU1.GT.HMAX) HMAX=HU1
           IF (HU1.LT.HMIN) HMIN=HU1
           IF (DABS(HM1-H0).LT.DABS(HM1)*1.0D-15) GO TO 20
15         H0=HM1
20      DA1=LOG10(HMAX)
        IF (HMIN.NE.0.0) DA2=LOG10(HMIN)
        ID=15-ABS(DA1-DA2)
        HM1=HM1*DLOG(X)
        S0=0.0D0
        DO 25 M=1,N
           IF (B.GE.0.0) S0=S0-1.0D0/M
25         IF (B.LT.0.0) S0=S0+(1.0D0-A)/(M*(A+M-1.0D0))
        HM2=PS+2.0D0*EL+S0
        R=1.0D0
        HMAX=0.0D0
        HMIN=1.0D+300
        DO 50 K=1,150
           S1=0.0D0
           S2=0.0D0
           IF (B.GT.0.0) THEN
              DO 30 M=1,K
30               S1=S1-(M+2.0D0*A-2.0D0)/(M*(M+A-1.0D0))
              DO 35 M=1,N
35               S2=S2+1.0D0/(K+M)
           ELSE
              DO 40 M=1,K+N
40               S1=S1+(1.0D0-A)/(M*(M+A-1.0D0))
              DO 45 M=1,K
45               S2=S2+1.0D0/M
           ENDIF
           HW=2.0D0*EL+PS+S1-S2
           R=R*(A0+K-1.0D0)*X/((N+K)*K)
           HM2=HM2+R*HW
           HU2=DABS(HM2)
           IF (HU2.GT.HMAX) HMAX=HU2
           IF (HU2.LT.HMIN) HMIN=HU2
           IF (DABS((HM2-H0)/HM2).LT.1.0D-15) GO TO 55
50         H0=HM2
55      DB1=LOG10(HMAX)
        IF (HMIN.NE.0.0) DB2=LOG10(HMIN)
        ID1=15-ABS(DB1-DB2)
        IF (ID1.LT.ID) ID=ID1
        HM3=1.0D0
        IF (N.EQ.0) HM3=0.0D0
        R=1.0D0
        DO 60 K=1,N-1
           R=R*(A2+K-1.0D0)/((K-N)*K)*X
60         HM3=HM3+R
        SA=UA*(HM1+HM2)
        SB=UB*HM3
        HU=SA+SB
        IF (SA.NE.0.0) ID1=INT(LOG10(ABS(SA)))
        IF (HU.NE.0.0) ID2=INT(LOG10(ABS(HU)))
        IF (SA*SB.LT.0.0) ID=ID-ABS(ID1-ID2)
        RETURN
        END


        SUBROUTINE CHGUIT(A,B,X,HU,ID)

C       ======================================================
C       Purpose: Compute hypergeometric function U(a,b,x) by
C                using Gaussian-Legendre integration (n=60)
C       Input  : a  --- Parameter ( a > 0 )
C                b  --- Parameter
C                x  --- Argument ( x > 0 )
C       Output:  HU --- U(a,b,z)
C                ID --- Estimated number of significant digits
C       Routine called: GAMMA for computing ג(x)
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION T(30),W(30)
        DATA T/ .259597723012478D-01, .778093339495366D-01,
     &          .129449135396945D+00, .180739964873425D+00,
     &          .231543551376029D+00, .281722937423262D+00,
     &          .331142848268448D+00, .379670056576798D+00,
     &          .427173741583078D+00, .473525841761707D+00,
     &          .518601400058570D+00, .562278900753945D+00,
     &          .604440597048510D+00, .644972828489477D+00,
     &          .683766327381356D+00, .720716513355730D+00,
     &          .755723775306586D+00, .788693739932264D+00,
     &          .819537526162146D+00, .848171984785930D+00,
     &          .874519922646898D+00, .898510310810046D+00,
     &          .920078476177628D+00, .939166276116423D+00,
     &          .955722255839996D+00, .969701788765053D+00,
     &          .981067201752598D+00, .989787895222222D+00,
     &          .995840525118838D+00, .999210123227436D+00/
        DATA W/ .519078776312206D-01, .517679431749102D-01,
     &          .514884515009810D-01, .510701560698557D-01,
     &          .505141845325094D-01, .498220356905502D-01,
     &          .489955754557568D-01, .480370318199712D-01,
     &          .469489888489122D-01, .457343797161145D-01,
     &          .443964787957872D-01, .429388928359356D-01,
     &          .413655512355848D-01, .396806954523808D-01,
     &          .378888675692434D-01, .359948980510845D-01,
     &          .340038927249464D-01, .319212190192963D-01,
     &          .297524915007890D-01, .275035567499248D-01,
     &          .251804776215213D-01, .227895169439978D-01,
     &          .203371207294572D-01, .178299010142074D-01,
     &          .152746185967848D-01, .126781664768159D-01,
     &          .100475571822880D-01, .738993116334531D-02,
     &          .471272992695363D-02, .202681196887362D-02/
        ID=7
        A1=A-1.0D0
        B1=B-A-1.0D0
        C=12.0/X
        DO 20 M=10,100,5
           HU1=0.0D0
           G=0.5D0*C/M
           D=G
           DO 15 J=1,M
              S=0.0D0
              DO 10 K=1,30
                 T1=D+G*T(K)
                 T2=D-G*T(K)
                 F1=DEXP(-X*T1)*T1**A1*(1.0D0+T1)**B1
                 F2=DEXP(-X*T2)*T2**A1*(1.0D0+T2)**B1
                 S=S+W(K)*(F1+F2)
10            CONTINUE
              HU1=HU1+S*G
              D=D+2.0D0*G
15         CONTINUE
           IF (DABS(1.0D0-HU0/HU1).LT.1.0D-7) GO TO 25
           HU0=HU1
20      CONTINUE
25      CALL GAMMA(A,GA)
        HU1=HU1/GA
        DO 40 M=2,10,2
           HU2=0.0D0
           G=0.5D0/M
           D=G
           DO 35 J=1,M
              S=0.0D0
              DO 30 K=1,30
                 T1=D+G*T(K)
                 T2=D-G*T(K)
                 T3=C/(1.0D0-T1)
                 T4=C/(1.0D0-T2)
                 F1=T3*T3/C*DEXP(-X*T3)*T3**A1*(1.0D0+T3)**B1
                 F2=T4*T4/C*DEXP(-X*T4)*T4**A1*(1.0D0+T4)**B1
                 S=S+W(K)*(F1+F2)
30            CONTINUE
              HU2=HU2+S*G
              D=D+2.0D0*G
35         CONTINUE
           IF (DABS(1.0D0-HU0/HU2).LT.1.0D-7) GO TO 45
           HU0=HU2
40      CONTINUE
45      CALL GAMMA(A,GA)
        HU2=HU2/GA
        HU=HU1+HU2
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)

C       ==================================================
C       Purpose: Compute gamma function ג(x)
C       Input :  x  --- Argument of ג(x)
C                       ( x is not equal to 0,-1,-2,תתת)
C       Output:  GA --- ג(x)
C       ==================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        SUBROUTINE PSI(X,PS)

C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
