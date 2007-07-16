        PROGRAM MAIRYZO

C       =========================================================
C       Purpose: This program computes the first NT zeros of Airy
C                functions Ai(x) and Ai'(x), and the associated
C                values of Ai(a') and Ai'(a), and the first NT
C                zeros of Airy functions Bi(x) and Bi'(x), and
C                the associated values of Bi(b') and Bi'(b) using
C                subroutine AIRYZO
C       Input :  NT    --- Total number of zeros
C                KF    --- Function code
C                          KF=1 for Ai(x) and Ai'(x)
C                          KF=2 for Bi(x) and Bi'(x)
C       Output:  XA(m) --- a, the m-th zero of Ai(x) or
C                          b, the m-th zero of Bi(x)
C                XB(m) --- a', the m-th zero of Ai'(x) or
C                          b', the m-th zero of Bi'(x)
C                XC(m) --- Ai(a') or Bi(b')
C                XD(m) --- Ai'(a) or Bi'(b)
C                          ( m --- Serial number of zeros )
C       Example: NT=5

C       m         a            Ai'(a)         a'          Ai(a')
C      -----------------------------------------------------------
C       1    -2.33810741     .70121082   -1.01879297    .53565666
C       2    -4.08794944    -.80311137   -3.24819758   -.41901548
C       3    -5.52055983     .86520403   -4.82009921    .38040647
C       4    -6.78670809    -.91085074   -6.16330736   -.35790794
C       5    -7.94413359     .94733571   -7.37217726    .34230124

C       m         b            Bi'(b)         b'          Bi(b')
C      -----------------------------------------------------------
C       1    -1.17371322     .60195789   -2.29443968   -.45494438
C       2    -3.27109330    -.76031014   -4.07315509    .39652284
C       3    -4.83073784     .83699101   -5.51239573   -.36796916
C       4    -6.16985213    -.88947990   -6.78129445    .34949912
C       5    -7.37676208     .92998364   -7.94017869   -.33602624
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XA(50),XB(50),XC(50),XD(50)
        WRITE(*,35)
        WRITE(*,40)
        WRITE(*,*)'Please enter KF,NT '
        READ(*,*)KF,NT
        WRITE(*,30)KF,NT
        IF (KF.EQ.1) THEN
           WRITE(*,*)'  m        a             Ai''(a)        a''',
     &               '           Ai(a'')'
        ELSE IF (KF.EQ.2) THEN
           WRITE(*,*)'  m        b             Bi''(b)        b''',
     &               '           Bi(b'')'
        ENDIF
        WRITE(*,*)'---------------------------------',
     &
     &            '---------------------------'
        CALL AIRYZO(NT,KF,XA,XB,XC,XD)
        DO 10 K=1,NT
10         WRITE(*,20)K,XA(K),XD(K),XB(K),XC(K)
20      FORMAT(1X,I3,1X,3F14.8,F13.8)
30      FORMAT(1X,3HKF=,I2,',     ',3HNT=,I3)
35      FORMAT(10X,'KF=1 for Ai(x) and Ai''(x); KF=2 for Bi(x)',
     &          ' and Bi''(x)')
40      FORMAT(10X,'NT is the number of the zeros')
        END


        SUBROUTINE AIRYZO(NT,KF,XA,XB,XC,XD)

C       ========================================================
C       Purpose: Compute the first NT zeros of Airy functions
C                Ai(x) and Ai'(x), a and a', and the associated
C                values of Ai(a') and Ai'(a); and the first NT
C                zeros of Airy functions Bi(x) and Bi'(x), b and
C                b', and the associated values of Bi(b') and
C                Bi'(b)
C       Input :  NT    --- Total number of zeros
C                KF    --- Function code
C                          KF=1 for Ai(x) and Ai'(x)
C                          KF=2 for Bi(x) and Bi'(x)
C       Output:  XA(m) --- a, the m-th zero of Ai(x) or
C                          b, the m-th zero of Bi(x)
C                XB(m) --- a', the m-th zero of Ai'(x) or
C                          b', the m-th zero of Bi'(x)
C                XC(m) --- Ai(a') or Bi(b')
C                XD(m) --- Ai'(a) or Bi'(b)
C                          ( m --- Serial number of zeros )
C       Routine called: AIRYB for computing Airy functions and
C                       their derivatives
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XA(NT),XB(NT),XC(NT),XD(NT)
        PI=3.141592653589793D0
        DO 15 I=1,NT
           IF (KF.EQ.1) THEN
              U=3.0*PI*(4.0*I-1)/8.0D0
              U1=1/(U*U)
              RT0=-(U*U)**(1.0/3.0)*((((-15.5902*U1+.929844)*U1
     &            -.138889)*U1+.10416667D0)*U1+1.0D0)
           ELSE IF (KF.EQ.2) THEN
              IF (I.EQ.1) THEN
                 RT0=-1.17371
              ELSE
                 U=3.0*PI*(4.0*I-3.0)/8.0
                 U1=1.0D0/(U*U)
                 RT0=-(U*U)**(1.0/3.0)*((((-15.5902*U1+.929844)*U1
     &               -.138889)*U1+.10416667)*U1+1.0)
              ENDIF
           ENDIF
10         X=RT0
           CALL AIRYB(X,AI,BI,AD,BD)
           IF (KF.EQ.1) RT=RT0-AI/AD
           IF (KF.EQ.2) RT=RT0-BI/BD
           IF (DABS((RT-RT0)/RT).GT.1.D-9) THEN
              RT0=RT
              GOTO 10
           ELSE
              XA(I)=RT
              IF (KF.EQ.1) XD(I)=AD
              IF (KF.EQ.2) XD(I)=BD
           ENDIF
15      CONTINUE
        DO 25 I=1,NT
           IF (KF.EQ.1) THEN
              IF (I.EQ.1) THEN
                 RT0=-1.01879
              ELSE
                 U=3.0*PI*(4.0*I-3.0)/8.0
                 U1=1/(U*U)
                 RT0=-(U*U)**(1.0/3.0)*((((15.0168*U1-.873954)
     &            *U1+.121528)*U1-.145833D0)*U1+1.0D0)
              ENDIF
           ELSE IF (KF.EQ.2) THEN
              IF (I.EQ.1) THEN
                 RT0=-2.29444
              ELSE
                 U=3.0*PI*(4.0*I-1.0)/8.0
                 U1=1.0/(U*U)
                 RT0=-(U*U)**(1.0/3.0)*((((15.0168*U1-.873954)
     &               *U1+.121528)*U1-.145833)*U1+1.0)
              ENDIF
           ENDIF
20         X=RT0
           CALL AIRYB(X,AI,BI,AD,BD)
           IF (KF.EQ.1) RT=RT0-AD/(AI*X)
           IF (KF.EQ.2) RT=RT0-BD/(BI*X)
           IF (DABS((RT-RT0)/RT).GT.1.0D-9) THEN
              RT0=RT
              GOTO 20
           ELSE
              XB(I)=RT
              IF (KF.EQ.1) XC(I)=AI
              IF (KF.EQ.2) XC(I)=BI
           ENDIF
25      CONTINUE
        RETURN
        END


        SUBROUTINE AIRYB(X,AI,BI,AD,BD)

C       =======================================================
C       Purpose: Compute Airy functions and their derivatives
C       Input:   x  --- Argument of Airy function
C       Output:  AI --- Ai(x)
C                BI --- Bi(x)
C                AD --- Ai'(x)
C                BD --- Bi'(x)
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION CK(41),DK(41)
        EPS=1.0D-15
        PI=3.141592653589793D0
        C1=0.355028053887817D0
        C2=0.258819403792807D0
        SR3=1.732050807568877D0
        XA=DABS(X)
        XQ=DSQRT(XA)
        IF (X.GT.0.0D0) XM=5.0
        IF (X.LE.0.0D0) XM=8.0
        IF (X.EQ.0.0D0) THEN
           AI=C1
           BI=SR3*C1
           AD=-C2
           BD=SR3*C2
           RETURN
        ENDIF
        IF (XA.LE.XM) THEN
           FX=1.0D0
           R=1.0D0
           DO 10 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K-1.0D0)*X
              FX=FX+R
              IF (DABS(R/FX).LT.EPS) GO TO 15
10         CONTINUE
15         GX=X
           R=X
           DO 20 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K+1.0D0)*X
              GX=GX+R
              IF (DABS(R/GX).LT.EPS) GO TO 25
20         CONTINUE
25         AI=C1*FX-C2*GX
           BI=SR3*(C1*FX+C2*GX)
           DF=.5D0*X*X
           R=DF
           DO 30 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K+2.0D0)*X
              DF=DF+R
              IF (DABS(R/DF).LT.EPS) GO TO 35
30         CONTINUE
35         DG=1.0D0
           R=1.0D0
           DO 40 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K-2.0D0)*X
              DG=DG+R
              IF (DABS(R/DG).LT.EPS) GO TO 45
40         CONTINUE
45         AD=C1*DF-C2*DG
           BD=SR3*(C1*DF+C2*DG)
        ELSE
           XE=XA*XQ/1.5D0
           XR1=1.0D0/XE
           XAR=1.0D0/XQ
           XF=DSQRT(XAR)
           RP=.5641895835477563D0
           R=1.0D0
           DO 50 K=1,40
              R=R*(6.0D0*K-1.0D0)/216.0D0*(6.0D0*K-3.0D0)
     &          /K*(6.0D0*K-5.0D0)/(2.0D0*K-1.0D0)
              CK(K)=R
50            DK(K)=-(6.0D0*K+1.0D0)/(6.0D0*K-1.0D0)*CK(K)
           KM=INT(24.5-XA)
           IF (XA.LT.6.0) KM=14
           IF (XA.GT.15.0) KM=10
           IF (X.GT.0.0D0) THEN
              SAI=1.0D0
              SAD=1.0D0
              R=1.0D0
              DO 55 K=1,KM
                 R=-R*XR1
                 SAI=SAI+CK(K)*R
55               SAD=SAD+DK(K)*R
              SBI=1.0D0
              SBD=1.0D0
              R=1.0D0
              DO 60 K=1,KM
                 R=R*XR1
                 SBI=SBI+CK(K)*R
60               SBD=SBD+DK(K)*R
              XP1=DEXP(-XE)
              AI=.5D0*RP*XF*XP1*SAI
              BI=RP*XF/XP1*SBI
              AD=-.5D0*RP/XF*XP1*SAD
              BD=RP/XF/XP1*SBD
           ELSE
              XCS=DCOS(XE+PI/4.0D0)
              XSS=DSIN(XE+PI/4.0D0)
              SSA=1.0D0
              SDA=1.0D0
              R=1.0D0
              XR2=1.0D0/(XE*XE)
              DO 65 K=1,KM
                 R=-R*XR2
                 SSA=SSA+CK(2*K)*R
65               SDA=SDA+DK(2*K)*R
              SSB=CK(1)*XR1
              SDB=DK(1)*XR1
              R=XR1
              DO 70 K=1,KM
                 R=-R*XR2
                 SSB=SSB+CK(2*K+1)*R
70               SDB=SDB+DK(2*K+1)*R
              AI=RP*XF*(XSS*SSA-XCS*SSB)
              BI=RP*XF*(XCS*SSA+XSS*SSB)
              AD=-RP/XF*(XCS*SDA+XSS*SDB)
              BD=RP/XF*(XSS*SDA-XCS*SDB)
           ENDIF
        ENDIF
        RETURN
        END
