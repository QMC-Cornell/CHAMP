        PROGRAM MKLVNA

C       =======================================================
C       Purpose: This program computes Kelvin functions ber x,
C                bei x, ker x and kei x, and their derivatives
C                using subroutine KLVNA
C       Input :  x   --- Argument of Kelvin functions
C       Output:  BER --- ber x
C                BEI --- bei x
C                GER --- ker x
C                GEI --- kei x
C                DER --- ber'x
C                DEI --- bei'x
C                HER --- ker'x
C                HEI --- kei'x
C       Example:

C      x       ber x          bei x          ker x          kei x
C    -----------------------------------------------------------------
C      0    .1000000D+01    0              ì            -.7853982D+00
C      5   -.6230082D+01   .1160344D+00  -.1151173D-01   .1118759D-01
C     10    .1388405D+03   .5637046D+02   .1294663D-03  -.3075246D-03
C     15   -.2967255D+04  -.2952708D+04  -.1514347D-07   .7962894D-05
C     20    .4748937D+05   .1147752D+06  -.7715233D-07  -.1858942D-06

C      x       ber'x          bei'x          ker'x          kei'x
C    -----------------------------------------------------------------
C      0     0              0            - ì              0
C      5   -.3845339D+01  -.4354141D+01   .1719340D-01  -.8199865D-03
C     10    .5119526D+02   .1353093D+03  -.3155969D-03   .1409138D-03
C     15    .9105533D+02  -.4087755D+04   .5644678D-05  -.5882223D-05
C     20   -.4880320D+05   .1118550D+06  -.7501859D-07   .1906243D-06
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        CALL KLVNA(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
        WRITE(*,*)'   x        ber x           bei x',
     &            '           ker x           kei x'
        WRITE(*,*)'--------------------------------',
     &            '--------------------------------------'
        WRITE(*,10)X,BER,BEI,GER,GEI
        WRITE(*,*)
        WRITE(*,*)'   x        ber''x           bei''x',
     &            '           ker''x           kei''x'
        WRITE(*,*)'--------------------------------',
     &            '--------------------------------------'
        WRITE(*,10)X,DER,DEI,HER,HEI
10      FORMAT(1X,F5.1,4D16.8)
        END


        SUBROUTINE KLVNA(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)

C       ======================================================
C       Purpose: Compute Kelvin functions ber x, bei x, ker x
C                and kei x, and their derivatives  ( x > 0 )
C       Input :  x   --- Argument of Kelvin functions
C       Output:  BER --- ber x
C                BEI --- bei x
C                GER --- ker x
C                GEI --- kei x
C                DER --- ber'x
C                DEI --- bei'x
C                HER --- ker'x
C                HEI --- kei'x
C       ================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        EPS=1.0D-15
        IF (X.EQ.0.0D0) THEN
           BER=1.0D0
           BEI=0.0D0
           GER=1.0D+300
           GEI=-0.25D0*PI
           DER=0.0D0
           DEI=0.0D0
           HER=-1.0D+300
           HEI=0.0D0
           RETURN
        ENDIF
        X2=0.25D0*X*X
        X4=X2*X2
        IF (DABS(X).LT.10.0D0) THEN
           BER=1.0D0
           R=1.0D0
           DO 10 M=1,60
              R=-0.25D0*R/(M*M)/(2.0D0*M-1.0D0)**2*X4
              BER=BER+R
              IF (DABS(R).LT.DABS(BER)*EPS) GO TO 15
10         CONTINUE
15         BEI=X2
           R=X2
           DO 20 M=1,60
              R=-0.25D0*R/(M*M)/(2.0D0*M+1.0D0)**2*X4
              BEI=BEI+R
              IF (DABS(R).LT.DABS(BEI)*EPS) GO TO 25
20         CONTINUE
25         GER=-(DLOG(X/2.0D0)+EL)*BER+0.25D0*PI*BEI
           R=1.0D0
           GS=0.0D0
           DO 30 M=1,60
              R=-0.25D0*R/(M*M)/(2.0D0*M-1.0D0)**2*X4
              GS=GS+1.0D0/(2.0D0*M-1.0D0)+1.0D0/(2.0D0*M)
              GER=GER+R*GS
              IF (DABS(R*GS).LT.DABS(GER)*EPS) GO TO 35
30         CONTINUE
35         GEI=X2-(DLOG(X/2.0D0)+EL)*BEI-0.25D0*PI*BER
           R=X2
           GS=1.0D0
           DO 40 M=1,60
              R=-0.25D0*R/(M*M)/(2.0D0*M+1.0D0)**2*X4
              GS=GS+1.0D0/(2.0D0*M)+1.0D0/(2.0D0*M+1.0D0)
              GEI=GEI+R*GS
              IF (DABS(R*GS).LT.DABS(GEI)*EPS) GO TO 45
40         CONTINUE
45         DER=-0.25D0*X*X2
           R=DER
           DO 50 M=1,60
              R=-0.25D0*R/M/(M+1.0D0)/(2.0D0*M+1.0D0)**2*X4
              DER=DER+R
              IF (DABS(R).LT.DABS(DER)*EPS) GO TO 55
50         CONTINUE
55         DEI=0.5D0*X
           R=DEI
           DO 60 M=1,60
              R=-0.25D0*R/(M*M)/(2.D0*M-1.D0)/(2.D0*M+1.D0)*X4
              DEI=DEI+R
              IF (DABS(R).LT.DABS(DEI)*EPS) GO TO 65
60            CONTINUE
65         R=-0.25D0*X*X2
           GS=1.5D0
           HER=1.5D0*R-BER/X-(DLOG(X/2.D0)+EL)*DER+0.25*PI*DEI
           DO 70 M=1,60
              R=-0.25D0*R/M/(M+1.0D0)/(2.0D0*M+1.0D0)**2*X4
              GS=GS+1.0D0/(2*M+1.0D0)+1.0D0/(2*M+2.0D0)
              HER=HER+R*GS
              IF (DABS(R*GS).LT.DABS(HER)*EPS) GO TO 75
70         CONTINUE
75         R=0.5D0*X
           GS=1.0D0
           HEI=0.5D0*X-BEI/X-(DLOG(X/2.D0)+EL)*DEI-0.25*PI*DER
           DO 80 M=1,60
              R=-0.25D0*R/(M*M)/(2*M-1.0D0)/(2*M+1.0D0)*X4
              GS=GS+1.0D0/(2.0D0*M)+1.0D0/(2*M+1.0D0)
              HEI=HEI+R*GS
              IF (DABS(R*GS).LT.DABS(HEI)*EPS) RETURN
80         CONTINUE
        ELSE
           PP0=1.0D0
           PN0=1.0D0
           QP0=0.0D0
           QN0=0.0D0
           R0=1.0D0
           KM=18
           IF (DABS(X).GE.40.0) KM=10
           FAC=1.0D0
           DO 85 K=1,KM
              FAC=-FAC
              XT=0.25D0*K*PI-INT(0.125D0*K)*2.0D0*PI
              CS=COS(XT)
              SS=SIN(XT)
              R0=0.125D0*R0*(2.0D0*K-1.0D0)**2/K/X
              RC=R0*CS
              RS=R0*SS
              PP0=PP0+RC
              PN0=PN0+FAC*RC
              QP0=QP0+RS
85            QN0=QN0+FAC*RS
           XD=X/DSQRT(2.0D0)
           XE1=DEXP(XD)
           XE2=DEXP(-XD)
           XC1=1.D0/DSQRT(2.0D0*PI*X)
           XC2=DSQRT(.5D0*PI/X)
           CP0=DCOS(XD+0.125D0*PI)
           CN0=DCOS(XD-0.125D0*PI)
           SP0=DSIN(XD+0.125D0*PI)
           SN0=DSIN(XD-0.125D0*PI)
           GER=XC2*XE2*(PN0*CP0-QN0*SP0)
           GEI=XC2*XE2*(-PN0*SP0-QN0*CP0)
           BER=XC1*XE1*(PP0*CN0+QP0*SN0)-GEI/PI
           BEI=XC1*XE1*(PP0*SN0-QP0*CN0)+GER/PI
           PP1=1.0D0
           PN1=1.0D0
           QP1=0.0D0
           QN1=0.0D0
           R1=1.0D0
           FAC=1.0D0
           DO 90 K=1,KM
              FAC=-FAC
              XT=0.25D0*K*PI-INT(0.125D0*K)*2.0D0*PI
              CS=DCOS(XT)
              SS=DSIN(XT)
              R1=0.125D0*R1*(4.D0-(2.0D0*K-1.0D0)**2)/K/X
              RC=R1*CS
              RS=R1*SS
              PP1=PP1+FAC*RC
              PN1=PN1+RC
              QP1=QP1+FAC*RS
              QN1=QN1+RS
90         CONTINUE
           HER=XC2*XE2*(-PN1*CN0+QN1*SN0)
           HEI=XC2*XE2*(PN1*SN0+QN1*CN0)
           DER=XC1*XE1*(PP1*CP0+QP1*SP0)-HEI/PI
           DEI=XC1*XE1*(PP1*SP0-QP1*CP0)+HER/PI
        ENDIF
        RETURN
        END
