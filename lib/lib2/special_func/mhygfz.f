        PROGRAM MHYGFZ

C       ============================================================
C       Purpose: This program computes hypergeometric function for
C                a complex argument, F(a,b,c,z), using subroutine
C                HYGFZ
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter,  c <> 0,-1,-2,...
C                x --- Real part of complex argument z
C                y --- Imaginary part of complex argument z
C                      ( z = x+iy )
C       Output:  ZHF --- F(a,b,c,z)
C       Examples:
C     a     b     c    z = x+ iy             F(a,b,c,z)
C   --------------------------------------------------------------
C    3.2   1.8   6.7   1.0+0.0 i    .54689992D+01+.00000000D+00 i
C    3.2  -1.8   6.7   1.0+0.0 i    .33750635D+00+.00000000D+00 i
C   -5.0   3.3   6.7   5.2+4.8 i    .11682745D+03+.60389104D+03 i
C    3.3  -6.0   3.7   5.2-4.8 i    .17620425D+05+.38293812D+05 i
C   -7.0   3.3  -3.7   5.2-4.8 i   -.11772779D+11-.14382286D+11 i
C    4.3  -8.0  -3.7   5.2+4.8 i    .13161188D+13-.10129870D+12 i
C    3.3   5.8   6.7   0.2+0.1 i    .17330557D+01+.63401030D+00 i
C    3.5  -2.4   6.7   0.2+0.5 i    .64762241D+00-.52110507D+00 i
C    3.3   4.3   6.7   0.8+0.3 i   -.14830086D+01+.83744258D+01 i
C    7.0   5.0   4.1   3.0-1.0 i   -.40376095D-02-.29566326D-02 i
C    5.0   7.0   4.1   3.0-1.0 i   -.40376095D-02-.29566326D-02 i
C    3.5   1.2   9.7   0.6+0.9 i    .10343044D+01+.54473814D+00 i
C    2.1   5.4   9.7   0.5+0.7 i    .68850442D+00+.12274187D+01 i
C    8.7   3.2   6.7   0.5+0.7 i   -.90046505D+00-.11198900D+01 i
C    8.7   2.7   6.7   0.6+0.9 i   -.46083890D+00-.54575701D+00 i
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Y)
        IMPLICIT COMPLEX *16 (Z)
        WRITE(*,*)'Please enter a,b,c,x and y '
        READ(*,*)A,B,C,X,Y
        Z=CMPLX(X,Y)
        CALL HYGFZ(A,B,C,Z,ZHF)
        WRITE(*,*)
        WRITE(*,*)'     a      b      c      x      y',
     &            '          Re[F]           Im[F]'
        WRITE(*,*)'   --------------------------------',
     &            '-----------------------------------'
        WRITE(*,10)A,B,C,X,Y,ZHF
10      FORMAT(1X,5F7.1,2X,2D16.8)
        END


        SUBROUTINE HYGFZ(A,B,C,Z,ZHF)

C       ======================================================
C       Purpose: Compute the hypergeometric function for a
C                complex argument, F(a,b,c,z)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter,  c <> 0,-1,-2,...
C                z --- Complex argument
C       Output:  ZHF --- F(a,b,c,z)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI for computing psi function
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Y)
        IMPLICIT COMPLEX *16 (Z)
        LOGICAL L0,L1,L2,L3,L4,L5,L6
        X=REAL(Z)
        Y=DIMAG(Z)
        EPS=1.0D-15
        L0=C.EQ.INT(C).AND.C.LT.0.0D0
        L1=DABS(1.0D0-X).LT.EPS.AND.Y.EQ.0.0D0.AND.C-A-B.LE.0.0D0
        L2=CDABS(Z+1.0D0).LT.EPS.AND.DABS(C-A+B-1.0D0).LT.EPS
        L3=A.EQ.INT(A).AND.A.LT.0.0D0
        L4=B.EQ.INT(B).AND.B.LT.0.0D0
        L5=C-A.EQ.INT(C-A).AND.C-A.LE.0.0D0
        L6=C-B.EQ.INT(C-B).AND.C-B.LE.0.0D0
        AA=A
        BB=B
        A0=CDABS(Z)
        IF (A0.GT.0.95D0) EPS=1.0D-8
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        IF (A0.EQ.0.0D0.OR.A.EQ.0.0D0.OR.B.EQ.0.0D0) THEN
           ZHF=(1.0D0,0.0D0)
        ELSE IF (Z.EQ.1.0D0.AND.C-A-B.GT.0.0D0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           ZHF=GC*GCAB/(GCA*GCB)
        ELSE IF (L2) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0D0+A/2.0D0-B,G2)
           CALL GAMMA(0.5D0+0.5D0*A,G3)
           ZHF=G0*G1/(G2*G3)
        ELSE IF (L3.OR.L4) THEN
           IF (L3) NM=INT(ABS(A))
           IF (L4) NM=INT(ABS(B))
           ZHF=(1.0D0,0.0D0)
           ZR=(1.0D0,0.0D0)
           DO 10 K=1,NM
              ZR=ZR*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*Z
10            ZHF=ZHF+ZR
        ELSE IF (L5.OR.L6) THEN
           IF (L5) NM=INT(ABS(C-A))
           IF (L6) NM=INT(ABS(C-B))
           ZHF=(1.0D0,0.0D0)
           ZR=(1.0D0,0.0D0)
           DO 15 K=1,NM
              ZR=ZR*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*Z
15            ZHF=ZHF+ZR
           ZHF=(1.0D0-Z)**(C-A-B)*ZHF
        ELSE IF (A0.LE.1.0D0) THEN
           IF (X.LT.0.0D0) THEN
              Z1=Z/(Z-1.0D0)
              IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
                 A=BB
                 B=AA
              ENDIF
              ZC0=1.0D0/((1.0D0-Z)**A)
              ZHF=(1.0D0,0.0D0)
              ZR0=(1.0D0,0.0D0)
              DO 20 K=1,500
                 ZR0=ZR0*(A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*Z1
                 ZHF=ZHF+ZR0
                 IF (CDABS(ZHF-ZW).LT.CDABS(ZHF)*EPS) GO TO 25
20               ZW=ZHF
25            ZHF=ZC0*ZHF
           ELSE IF (A0.GE.0.90D0) THEN
              GM=0.0D0
              MCAB=INT(C-A-B+EPS*DSIGN(1.0D0,C-A-B))
              IF (DABS(C-A-B-MCAB).LT.EPS) THEN
                 M=INT(C-A-B)
                 CALL GAMMA(A,GA)
                 CALL GAMMA(B,GB)
                 CALL GAMMA(C,GC)
                 CALL GAMMA(A+M,GAM)
                 CALL GAMMA(B+M,GBM)
                 CALL PSI(A,PA)
                 CALL PSI(B,PB)
                 IF (M.NE.0) GM=1.0D0
                 DO 30 J=1,ABS(M)-1
30                  GM=GM*J
                 RM=1.0D0
                 DO 35 J=1,ABS(M)
35                  RM=RM*J
                 ZF0=(1.0D0,0.0D0)
                 ZR0=(1.0D0,0.0D0)
                 ZR1=(1.0D0,0.0D0)
                 SP0=0.D0
                 SP=0.0D0
                 IF (M.GE.0) THEN
                    ZC0=GM*GC/(GAM*GBM)
                    ZC1=-GC*(Z-1.0D0)**M/(GA*GB*RM)
                    DO 40 K=1,M-1
                       ZR0=ZR0*(A+K-1.D0)*(B+K-1.D0)/(K*(K-M))*(1.D0-Z)
40                     ZF0=ZF0+ZR0
                    DO 45 K=1,M
45                     SP0=SP0+1.0D0/(A+K-1.0D0)+1.0/(B+K-1.0D0)-1.D0/K
                    ZF1=PA+PB+SP0+2.0D0*EL+CDLOG(1.0D0-Z)
                    DO 55 K=1,500
                       SP=SP+(1.0D0-A)/(K*(A+K-1.0D0))+(1.0D0-B)/
     &                    (K*(B+K-1.0D0))
                       SM=0.0D0
                       DO 50 J=1,M
                          SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0D0))
     &                       +1.0D0/(B+J+K-1.0D0)
50                     CONTINUE
                       ZP=PA+PB+2.0D0*EL+SP+SM+CDLOG(1.0D0-Z)
                       ZR1=ZR1*(A+M+K-1.0D0)*(B+M+K-1.0D0)/(K*(M+K))
     &                     *(1.0D0-Z)
                       ZF1=ZF1+ZR1*ZP
                       IF (CDABS(ZF1-ZW).LT.CDABS(ZF1)*EPS) GO TO 60
55                     ZW=ZF1
60                  ZHF=ZF0*ZC0+ZF1*ZC1
                 ELSE IF (M.LT.0) THEN
                    M=-M
                    ZC0=GM*GC/(GA*GB*(1.0D0-Z)**M)
                    ZC1=-(-1)**M*GC/(GAM*GBM*RM)
                    DO 65 K=1,M-1
                       ZR0=ZR0*(A-M+K-1.0D0)*(B-M+K-1.0D0)/(K*(K-M))
     &                     *(1.0D0-Z)
65                     ZF0=ZF0+ZR0
                    DO 70 K=1,M
70                     SP0=SP0+1.0D0/K
                    ZF1=PA+PB-SP0+2.0D0*EL+CDLOG(1.0D0-Z)
                    DO 80 K=1,500
                       SP=SP+(1.0D0-A)/(K*(A+K-1.0D0))+(1.0D0-B)/(K*
     &                    (B+K-1.0D0))
                       SM=0.0D0
                       DO 75 J=1,M
75                        SM=SM+1.0D0/(J+K)
                       ZP=PA+PB+2.0D0*EL+SP-SM+CDLOG(1.0D0-Z)
                       ZR1=ZR1*(A+K-1.D0)*(B+K-1.D0)/(K*(M+K))*(1.D0-Z)
                       ZF1=ZF1+ZR1*ZP
                       IF (CDABS(ZF1-ZW).LT.CDABS(ZF1)*EPS) GO TO 85
80                     ZW=ZF1
85                  ZHF=ZF0*ZC0+ZF1*ZC1
                 ENDIF
              ELSE
                 CALL GAMMA(A,GA)
                 CALL GAMMA(B,GB)
                 CALL GAMMA(C,GC)
                 CALL GAMMA(C-A,GCA)
                 CALL GAMMA(C-B,GCB)
                 CALL GAMMA(C-A-B,GCAB)
                 CALL GAMMA(A+B-C,GABC)
                 ZC0=GC*GCAB/(GCA*GCB)
                 ZC1=GC*GABC/(GA*GB)*(1.0D0-Z)**(C-A-B)
                 ZHF=(0.0D0,0.0D0)
                 ZR0=ZC0
                 ZR1=ZC1
                 DO 90 K=1,500
                    ZR0=ZR0*(A+K-1.D0)*(B+K-1.D0)/(K*(A+B-C+K))*(1.D0-Z)
                    ZR1=ZR1*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C-A-B+K))
     &                  *(1.0D0-Z)
                    ZHF=ZHF+ZR0+ZR1
                    IF (CDABS(ZHF-ZW).LT.CDABS(ZHF)*EPS) GO TO 95
90                  ZW=ZHF
95               ZHF=ZHF+ZC0+ZC1
              ENDIF
           ELSE
              Z00=(1.0D0,0.0D0)
              IF (C-A.LT.A.AND.C-B.LT.B) THEN
                  Z00=(1.0D0-Z)**(C-A-B)
                  A=C-A
                  B=C-B
              ENDIF
              ZHF=(1.0D0,0.D0)
              ZR=(1.0D0,0.0D0)
              DO 100 K=1,1500
                 ZR=ZR*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*Z
                 ZHF=ZHF+ZR
                 IF (CDABS(ZHF-ZW).LE.CDABS(ZHF)*EPS) GO TO 105
100              ZW=ZHF
105           ZHF=Z00*ZHF
           ENDIF
        ELSE IF (A0.GT.1.0D0) THEN
           MAB=INT(A-B+EPS*DSIGN(1.0D0,A-B))
           IF (DABS(A-B-MAB).LT.EPS.AND.A0.LE.1.1D0) B=B+EPS
           IF (DABS(A-B-MAB).GT.EPS) THEN
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A-B,GAB)
              CALL GAMMA(B-A,GBA)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              ZC0=GC*GBA/(GCA*GB*(-Z)**A)
              ZC1=GC*GAB/(GCB*GA*(-Z)**B)
              ZR0=ZC0
              ZR1=ZC1
              ZHF=(0.0D0,0.0D0)
              DO 110 K=1,500
                 ZR0=ZR0*(A+K-1.0D0)*(A-C+K)/((A-B+K)*K*Z)
                 ZR1=ZR1*(B+K-1.0D0)*(B-C+K)/((B-A+K)*K*Z)
                 ZHF=ZHF+ZR0+ZR1
                 IF (CDABS((ZHF-ZW)/ZHF).LE.EPS) GO TO 115
110              ZW=ZHF
115           ZHF=ZHF+ZC0+ZC1
           ELSE
              IF (A-B.LT.0.0D0) THEN
                 A=BB
                 B=AA
              ENDIF
              CA=C-A
              CB=C-B
              NCA=INT(CA+EPS*DSIGN(1.0D0,CA))
              NCB=INT(CB+EPS*DSIGN(1.0D0,CB))
              IF (DABS(CA-NCA).LT.EPS.OR.DABS(CB-NCB).LT.EPS) C=C+EPS
              CALL GAMMA(A,GA)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-B,GCB)
              CALL PSI(A,PA)
              CALL PSI(C-A,PCA)
              CALL PSI(A-C,PAC)
              MAB=INT(A-B+EPS)
              ZC0=GC/(GA*(-Z)**B)
              CALL GAMMA(A-B,GM)
              ZF0=GM/GCB*ZC0
              ZR=ZC0
              DO 120 K=1,MAB-1
                 ZR=ZR*(B+K-1.0D0)/(K*Z)
                 T0=A-B-K
                 CALL GAMMA(T0,G0)
                 CALL GAMMA(C-B-K,GCBK)
120              ZF0=ZF0+ZR*G0/GCBK
              IF (MAB.EQ.0) ZF0=(0.0D0,0.0D0)
              ZC1=GC/(GA*GCB*(-Z)**A)
              SP=-2.0D0*EL-PA-PCA
              DO 125 J=1,MAB
125              SP=SP+1.0D0/J
              ZP0=SP+CDLOG(-Z)
              SQ=1.0D0
              DO 130 J=1,MAB
130              SQ=SQ*(B+J-1.0D0)*(B-C+J)/J
              ZF1=(SQ*ZP0)*ZC1
              ZR=ZC1
              RK1=1.0D0
              SJ1=0.0D0
              DO 145 K=1,10000
                 ZR=ZR/Z
                 RK1=RK1*(B+K-1.0D0)*(B-C+K)/(K*K)
                 RK2=RK1
                 DO 135 J=K+1,K+MAB
135                 RK2=RK2*(B+J-1.0D0)*(B-C+J)/J
                 SJ1=SJ1+(A-1.0D0)/(K*(A+K-1.0D0))+(A-C-1.0D0)/
     &               (K*(A-C+K-1.0D0))
                 SJ2=SJ1
                 DO 140 J=K+1,K+MAB
140                 SJ2=SJ2+1.0D0/J
                 ZP=-2.0D0*EL-PA-PAC+SJ2-1.0D0/(K+A-C)
     &              -PI/DTAN(PI*(K+A-C))+CDLOG(-Z)
                 ZF1=ZF1+RK2*ZR*ZP
                 WS=CDABS(ZF1)
                 IF (DABS((WS-W0)/WS).LT.EPS) GO TO 150
145              W0=WS
150           ZHF=ZF0+ZF1
           ENDIF
        ENDIF
155     A=AA
        B=BB
        IF (K.GT.150) WRITE(*,160)
160     FORMAT(1X,'Warning! You should check the accuracy')
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
