        PROGRAM MJYV

C       ============================================================
C       Purpose: This program computes Bessel functions Jv(x) and
C                Yv(x) and their derivatives using subroutine JYV
C       Input :  x --- Argument of Jv(x) and Yv(x)
C                v --- Order of Jv(x) and Yv(x)
C                      ( v = n+v0,  0 ף n ף 250, 0 ף v0 < 1 )
C       Output:  BJ(n) --- Jn+v0(x)
C                DJ(n) --- Jn+v0'(x)
C                BY(n) --- Yn+v0(x)
C                DY(n) --- Yn+v0'(x)
C       Example: Compute Jv(x) and Yv(x) and their derivatives
C                for v = 0.25(1.0)5.25 and x = 10.0
C                Computation results:

C                v =  5.25,      x = 10.00

C        v        Jv(x)         Jv'(x)        Yv(x)         Yv'(x)
C       ------------------------------------------------------------
C        .25   -.20639379    -.13476340     .14493044    -.21381777
C       1.25    .12960355    -.22259423     .21744103     .11775031
C       2.25    .23879467     .07587475    -.09057018     .23781932
C       3.25   -.02214595     .24599211    -.25819761    -.00665596
C       4.25   -.25318954     .08545961    -.07725827    -.22536285
C       5.25   -.19306516    -.15183033     .19252809    -.17833551
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMMON BJ(0:250),DJ(0:250),BY(0:250),DY(0:250)
        WRITE(*,*)'  Please enter v, x '
        READ(*,*)V,X
        WRITE(*,20)V,X
        IF (V.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL JYV(V,X,VM,BJ,DJ,BY,DY)
        NM=INT(VM)
        V0=VM-NM
        WRITE(*,*)
        WRITE(*,*)'   v         Jv(x)           Jv''(x)',
     &            '          Yv(x)           Yv''(x)'
        WRITE(*,*)' ---------------------------------------------',
     &            '------------------------'
        DO 10 K=0,NM,NS
           VK=K+V0
10         WRITE(*,15)VK,BJ(K),DJ(K),BY(K),DY(K)
15      FORMAT(1X,F6.2,4D16.8)
20      FORMAT(8X,3Hv =,F6.2,',    ',3Hx =,F6.2)
        END


        SUBROUTINE JYV(V,X,VM,BJ,DJ,BY,DY)

C       =======================================================
C       Purpose: Compute Bessel functions Jv(x) and Yv(x)
C                and their derivatives
C       Input :  x --- Argument of Jv(x) and Yv(x)
C                v --- Order of Jv(x) and Yv(x)
C                      ( v = n+v0, 0 ף v0 < 1, n = 0,1,2,... )
C       Output:  BJ(n) --- Jn+v0(x)
C                DJ(n) --- Jn+v0'(x)
C                BY(n) --- Yn+v0(x)
C                DY(n) --- Yn+v0'(x)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(0:*),DJ(0:*),BY(0:*),DY(0:*)
        EL=.5772156649015329D0
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        X2=X*X
        N=INT(V)
        V0=V-N
        IF (X.LT.1.0D-100) THEN
           DO 10 K=0,N
              BJ(K)=0.0D0
              DJ(K)=0.0D0
              BY(K)=-1.0D+300
10            DY(K)=1.0D+300
           IF (V0.EQ.0.0) THEN
              BJ(0)=1.0D0
              DJ(1)=0.5D0
           ELSE
              DJ(0)=1.0D+300
           ENDIF
           VM=V
           RETURN
        ENDIF
        IF (X.LE.12.0) THEN
           DO 25 L=0,1
              VL=V0+L
              BJVL=1.0D0
              R=1.0D0
              DO 15 K=1,40
                 R=-0.25D0*R*X2/(K*(K+VL))
                 BJVL=BJVL+R
                 IF (DABS(R).LT.DABS(BJVL)*1.0D-15) GO TO 20
15            CONTINUE
20            VG=1.0D0+VL
              CALL GAMMA(VG,GA)
              A=(0.5D0*X)**VL/GA
              IF (L.EQ.0) BJV0=BJVL*A
              IF (L.EQ.1) BJV1=BJVL*A
25         CONTINUE
        ELSE
           K0=11
           IF (X.GE.35.0) K0=10
           IF (X.GE.50.0) K0=8
           DO 40 J=0,1
              VV=4.0D0*(J+V0)*(J+V0)
              PX=1.0D0
              RP=1.0D0
              DO 30 K=1,K0
                 RP=-0.78125D-2*RP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &              (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*X2)
30               PX=PX+RP
              QX=1.0D0
              RQ=1.0D0
              DO 35 K=1,K0
                 RQ=-0.78125D-2*RQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &              (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*X2)
35               QX=QX+RQ
              QX=0.125D0*(VV-1.0)*QX/X
              XK=X-(0.5D0*(J+V0)+0.25D0)*PI
              A0=DSQRT(RP2/X)
              CK=DCOS(XK)
              SK=DSIN(XK)
              IF (J.EQ.0) THEN
                 BJV0=A0*(PX*CK-QX*SK)
                 BYV0=A0*(PX*SK+QX*CK)
              ELSE IF (J.EQ.1) THEN
                 BJV1=A0*(PX*CK-QX*SK)
                 BYV1=A0*(PX*SK+QX*CK)
              ENDIF
40         CONTINUE
        ENDIF
        BJ(0)=BJV0
        BJ(1)=BJV1
        DJ(0)=V0/X*BJ(0)-BJ(1)
        DJ(1)=-(1.0D0+V0)/X*BJ(1)+BJ(0)
        IF (N.GE.2.AND.N.LE.INT(0.9*X)) THEN
           F0=BJV0
           F1=BJV1
           DO 45 K=2,N
              F=2.0D0*(K+V0-1.0D0)/X*F1-F0
              BJ(K)=F
              F0=F1
45            F1=F
        ELSE IF (N.GE.2) THEN
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              N=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F2=0.0D0
           F1=1.0D-100
           DO 50 K=M,0,-1
              F=2.0D0*(V0+K+1.0D0)/X*F1-F2
              IF (K.LE.N) BJ(K)=F
              F2=F1
50            F1=F
           IF (DABS(BJV0).GT.DABS(BJV1)) THEN
               CS=BJV0/F
           ELSE
               CS=BJV1/F2
           ENDIF
           DO 55 K=0,N
55            BJ(K)=CS*BJ(K)
        ENDIF
        DO 60 K=2,N
60         DJ(K)=-(K+V0)/X*BJ(K)+BJ(K-1)
        IF (X.LE.12.0D0) THEN
           IF (V0.NE.0.0) THEN
              DO 75 L=0,1
                 VL=V0+L
                 BJVL=1.0D0
                 R=1.0D0
                 DO 65 K=1,40
                    R=-0.25D0*R*X2/(K*(K-VL))
                    BJVL=BJVL+R
                    IF (DABS(R).LT.DABS(BJVL)*1.0D-15) GO TO 70
65               CONTINUE
70               VG=1.0D0-VL
                 CALL GAMMA(VG,GB)
                 B=(2.0D0/X)**VL/GB
                 IF (L.EQ.0) BJU0=BJVL*B
                 IF (L.EQ.1) BJU1=BJVL*B
75            CONTINUE
              PV0=PI*V0
              PV1=PI*(1.0D0+V0)
              BYV0=(BJV0*DCOS(PV0)-BJU0)/DSIN(PV0)
              BYV1=(BJV1*DCOS(PV1)-BJU1)/DSIN(PV1)
           ELSE
              EC=DLOG(X/2.0D0)+EL
              CS0=0.0D0
              W0=0.0D0
              R0=1.0D0
              DO 80 K=1,30
                 W0=W0+1.0D0/K
                 R0=-0.25D0*R0/(K*K)*X2
80               CS0=CS0+R0*W0
              BYV0=RP2*(EC*BJV0-CS0)
              CS1=1.0D0
              W1=0.0D0
              R1=1.0D0
              DO 85 K=1,30
                 W1=W1+1.0D0/K
                 R1=-0.25D0*R1/(K*(K+1))*X2
85               CS1=CS1+R1*(2.0D0*W1+1.0D0/(K+1.0D0))
              BYV1=RP2*(EC*BJV1-1.0D0/X-0.25D0*X*CS1)
           ENDIF
        ENDIF
        BY(0)=BYV0
        BY(1)=BYV1
        DO 90 K=2,N
           BYVK=2.0D0*(V0+K-1.0D0)/X*BYV1-BYV0
           BY(K)=BYVK
           BYV0=BYV1
90         BYV1=BYVK
        DY(0)=V0/X*BY(0)-BY(1)
        DO 95 K=1,N
95         DY(K)=-(K+V0)/X*BY(K)+BY(K-1)
        VM=N+V0
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


        INTEGER FUNCTION MSTA1(X,MP)

C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that the magnitude of
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)

C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
