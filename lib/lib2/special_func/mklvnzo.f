        PROGRAM MKLVNZO

C       ==============================================================
C       Purpose: This program computes the first NT zeros of Kelvin
C                functions and their derivatives using subroutine
C                KLVNZO
C       Input :  NT --- Total number of zeros
C       Example: NT = 5

C           Zeros of Kelvin functions ber x, bei x, ker x and kei x

C        m       ber x          bei x          ker x          kei x
C       ---------------------------------------------------------------
C        1     2.84891782     5.02622395     1.71854296     3.91466761
C        2     7.23882945     9.45540630     6.12727913     8.34422506
C        3    11.67396355    13.89348785    10.56294271    12.78255715
C        4    16.11356383    18.33398346    15.00268812    17.22314372
C        5    20.55463158    22.77543929    19.44381663    21.66464214

C          Zeros of Kelvin Functions ber'x, bei'x, ker'x and kei'x

C        m       ber'x          bei'x          ker'x          kei'x
C       ---------------------------------------------------------------
C        1     6.03871081     3.77267330     2.66583979     4.93181194
C        2    10.51364251     8.28098785     7.17212212     9.40405458
C        3    14.96844542    12.74214752    11.63218639    13.85826916
C        4    19.41757493    17.19343175    16.08312025    18.30717294
C        5    23.86430432    21.64114394    20.53067845    22.75379258
C       ==============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION R1(50),R2(50),R3(50),R4(50),R5(50),R6(50),
     &            R7(50),R8(50)
        WRITE(*,35)
        WRITE(*,*)'Please enter NT '
        READ(*,*)NT
        WRITE(*,25)
        WRITE(*,*)
        WRITE(*,*)'   m       ber x          bei x          ker x',
     &            '          kei x'
        WRITE(*,*)' ------------------------------------------------',
     &            '---------------'
        CALL KLVNZO(NT,1,R1)
        CALL KLVNZO(NT,2,R2)
        CALL KLVNZO(NT,3,R3)
        CALL KLVNZO(NT,4,R4)
        DO 10 L=1,NT
10         WRITE(*,20)L,R1(L),R2(L),R3(L),R4(L)
        CALL KLVNZO(NT,5,R5)
        CALL KLVNZO(NT,6,R6)
        CALL KLVNZO(NT,7,R7)
        CALL KLVNZO(NT,8,R8)
        WRITE(*,*)
        WRITE(*,30)
        WRITE(*,*)
        WRITE(*,*)'   m       ber''x          bei''x          ker''x',
     &            '          kei''x'
        WRITE(*,*)' ------------------------------------------------',
     &            '---------------'
        DO 15 L=1,NT
15         WRITE(*,20)L,R5(L),R6(L),R7(L),R8(L)
        CLOSE (01)
20      FORMAT(1X,I3,1X,F14.8,1X,F14.8,1X,F14.8,1X,F14.8)
25      FORMAT(4X,'Zeros of Kelvin functions ber x, bei x,'
     &        ,' ker x and kei x')
30      FORMAT(4X,'Zeros of Kelvin functions ber''x, bei''x,'
     &        ,' ker''x and kei''x')
35      FORMAT(1X,'NT is the number of the zeros')
        END


        SUBROUTINE KLVNZO(NT,KD,ZO)

C       ====================================================
C       Purpose: Compute the zeros of Kelvin functions
C       Input :  NT  --- Total number of zeros
C                KD  --- Function code
C                KD=1 to 8 for ber x, bei x, ker x, kei x,
C                          ber'x, bei'x, ker'x and kei'x,
C                          respectively.
C       Output:  ZO(M) --- the M-th zero of Kelvin function
C                          for code KD
C       Routine called:
C                KLVNA for computing Kelvin functions and
C                their derivatives
C       ====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION ZO(NT),RT0(8)
        RT0(1)=2.84891
        RT0(2)=5.02622
        RT0(3)=1.71854
        RT0(4)=3.91467
        RT0(5)=6.03871
        RT0(6)=3.77268
        RT0(7)=2.66584
        RT0(8)=4.93181
        RT=RT0(KD)
        DO 15 M=1,NT
10         CALL KLVNA(RT,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
           IF (KD.EQ.1) THEN
              RT=RT-BER/DER
           ELSE IF (KD.EQ.2) THEN
              RT=RT-BEI/DEI
           ELSE IF (KD.EQ.3) THEN
              RT=RT-GER/HER
           ELSE IF (KD.EQ.4) THEN
              RT=RT-GEI/HEI
           ELSE IF (KD.EQ.5) THEN
              DDR=-BEI-DER/RT
              RT=RT-DER/DDR
           ELSE IF (KD.EQ.6) THEN
              DDI=BER-DEI/RT
              RT=RT-DEI/DDI
           ELSE IF (KD.EQ.7) THEN
              GDR=-GEI-HER/RT
              RT=RT-HER/GDR
           ELSE
              GDI=GER-HEI/RT
              RT=RT-HEI/GDI
           ENDIF
           IF (DABS(RT-RT0(KD)).GT.5.0D-10) THEN
              RT0(KD)=RT
              GO TO 10
           ENDIF
           ZO(M)=RT
15         RT=RT+4.44D0
        RETURN
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
           GEI=-.25D0*PI
           DER=0.0D0
           DEI=0.0D0
           HER=-1.0D+300
           HEI=0.0D0
           RETURN
        ENDIF
        X2=.25D0*X*X
        X4=X2*X2
        IF (DABS(X).LT.10.0D0) THEN
           BER=1.0D0
           R=1.0D0
           DO 10 M=1,60
              R=-.25D0*R/(M*M)/(2.0D0*M-1.0D0)**2*X4
              BER=BER+R
              IF (DABS(R/BER).LT.EPS) GO TO 15
10         CONTINUE
15         BEI=X2
           R=X2
           DO 20 M=1,60
              R=-.25D0*R/(M*M)/(2.0D0*M+1.0D0)**2*X4
              BEI=BEI+R
              IF (DABS(R/BEI).LT.EPS) GO TO 25
20         CONTINUE
25         GER=-(DLOG(X/2.0D0)+EL)*BER+.25D0*PI*BEI
           R=1.0D0
           GS=0.0D0
           DO 30 M=1,60
              R=-.25D0*R/(M*M)/(2.0D0*M-1.0D0)**2*X4
              GS=GS+1.0D0/(2.0D0*M-1.0D0)+1.0D0/(2.0D0*M)
              GER=GER+R*GS
              IF (DABS(R*GS/GER).LT.EPS) GO TO 35
30         CONTINUE
35         GEI=X2-(DLOG(X/2.0D0)+EL)*BEI-.25D0*PI*BER
           R=X2
           GS=1.0D0
           DO 40 M=1,60
              R=-.25D0*R/(M*M)/(2.0D0*M+1.0D0)**2*X4
              GS=GS+1.0D0/(2.0D0*M)+1.0D0/(2.0D0*M+1.0D0)
              GEI=GEI+R*GS
              IF (DABS(R*GS/GEI).LT.EPS) GO TO 45
40         CONTINUE
45         DER=-.25D0*X*X2
           R=DER
           DO 50 M=1,60
              R=-.25D0*R/M/(M+1.0D0)/(2.0D0*M+1.0D0)**2*X4
              DER=DER+R
              IF (DABS(R/DER).LT.EPS) GO TO 55
50         CONTINUE
55         DEI=.5D0*X
           R=DEI
           DO 60 M=1,60
              R=-.25D0*R/(M*M)/(2.D0*M-1.D0)/(2.D0*M+1.D0)*X4
              DEI=DEI+R
              IF (DABS(R/DEI).LT.EPS) GO TO 65
60            CONTINUE
65         R=-.25D0*X*X2
           GS=1.5D0
           HER=1.5D0*R-BER/X-(DLOG(X/2.D0)+EL)*DER+.25*PI*DEI
           DO 70 M=1,60
              R=-.25D0*R/M/(M+1.0D0)/(2.0D0*M+1.0D0)**2*X4
              GS=GS+1.0D0/(2*M+1.0D0)+1.0D0/(2*M+2.0D0)
              HER=HER+R*GS
              IF (DABS(R*GS/HER).LT.EPS) GO TO 75
70         CONTINUE
75         R=.5D0*X
           GS=1.0D0
           HEI=.5D0*X-BEI/X-(DLOG(X/2.D0)+EL)*DEI-.25*PI*DER
           DO 80 M=1,60
              R=-.25D0*R/(M*M)/(2*M-1.0D0)/(2*M+1.0D0)*X4
              GS=GS+1.0D0/(2.0D0*M)+1.0D0/(2*M+1.0D0)
              HEI=HEI+R*GS
              IF (DABS(R*GS/HEI).LT.EPS) RETURN
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
              XT=.25D0*K*PI-INT(.125D0*K)*2.0D0*PI
              CS=COS(XT)
              SS=SIN(XT)
              R0=.125D0*R0*(2.0D0*K-1.0D0)**2/K/X
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
           CP0=DCOS(XD+.125D0*PI)
           CN0=DCOS(XD-.125D0*PI)
           SP0=DSIN(XD+.125D0*PI)
           SN0=DSIN(XD-.125D0*PI)
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
              XT=.25D0*K*PI-INT(.125D0*K)*2.0D0*PI
              CS=DCOS(XT)
              SS=DSIN(XT)
              R1=.125D0*R1*(4.D0-(2.0D0*K-1.0D0)**2)/K/X
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
