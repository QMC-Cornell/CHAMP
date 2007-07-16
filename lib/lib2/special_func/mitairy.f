        PROGRAM MITAIRY

C       ===========================================================
C       Purpose: This program computes the integrals of Airy
C                functions using subroutine ITAIRY
C       Input  : x   --- Upper limit of the integral
C       Output : APT --- Integration of Ai(t) from 0 and x
C                BPT --- Integration of Bi(t) from 0 and x
C                ANT --- Integration of Ai(-t) from 0 and x
C                BNT --- Integration of Bi(-t) from 0 and x
C       Example:

C         x      Ai(t)dt       Bi(t)dt       Ai(-t)dt     Bi(-t)dt
C        ----------------------------------------------------------
C         5    .33328759   .32147832D+03    .71788220    .15873094
C        10    .33333333   .14780980D+09    .76569840    .01504043
C        15    .33333333   .49673090D+16    .68358063    .07202621
C        20    .33333333   .47447423D+25    .71173925   -.03906173
C        25    .33333333   .78920820D+35    .70489539    .03293190
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,20)
        WRITE(*,30)
        CALL ITAIRY(X,APT,BPT,ANT,BNT)
        WRITE(*,10)X,APT,BPT,ANT,BNT
10      FORMAT(1X,F5.1,F14.8,2X,D15.8,2F14.8)
20      FORMAT(3X,'x',8X,'Ai(t)dt',7X,'Bi(t)dt',9X,
     &        'Ai(-t)dt',6X,'Bi(-t)dt')
30      FORMAT(2X,'----------------------------------',
     &     '------------------------------')
        END


        SUBROUTINE ITAIRY(X,APT,BPT,ANT,BNT)

C       ======================================================
C       Purpose: Compute the integrals of Airy fnctions with
C                respect to t from 0 and x ( x ò 0 )
C       Input  : x   --- Upper limit of the integral
C       Output : APT --- Integration of Ai(t) from 0 and x
C                BPT --- Integration of Bi(t) from 0 and x
C                ANT --- Integration of Ai(-t) from 0 and x
C                BNT --- Integration of Bi(-t) from 0 and x
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(16)
        EPS=1.0D-15
        PI=3.141592653589793D0
        C1=.355028053887817D0
        C2=.258819403792807D0
        SR3=1.732050807568877D0
        IF (X.EQ.0.0D0) THEN
           APT=0.0D0
           BPT=0.0D0
           ANT=0.0D0
           BNT=0.0D0
        ELSE
           IF (DABS(X).LE.9.25D0) THEN
              DO 30 L=0,1
                 X=(-1)**L*X
                 FX=X
                 R=X
                 DO 10 K=1,40
                    R=R*(3.0*K-2.0D0)/(3.0*K+1.0D0)*X/(3.0*K)
     &                *X/(3.0*K-1.0D0)*X
                    FX=FX+R
                    IF (DABS(R).LT.DABS(FX)*EPS) GO TO 15
10               CONTINUE
15               GX=.5D0*X*X
                 R=GX
                 DO 20 K=1,40
                    R=R*(3.0*K-1.0D0)/(3.0*K+2.0D0)*X/(3.0*K)
     &                *X/(3.0*K+1.0D0)*X
                    GX=GX+R
                    IF (DABS(R).LT.DABS(GX)*EPS) GO TO 25
20               CONTINUE
25               ANT=C1*FX-C2*GX
                 BNT=SR3*(C1*FX+C2*GX)
                 IF (L.EQ.0) THEN
                    APT=ANT
                    BPT=BNT
                 ELSE
                    ANT=-ANT
                    BNT=-BNT
                    X=-X
                 ENDIF
30            CONTINUE
           ELSE
              DATA A/.569444444444444D0,.891300154320988D0,
     &             .226624344493027D+01,.798950124766861D+01,
     &             .360688546785343D+02,.198670292131169D+03,
     &             .129223456582211D+04,.969483869669600D+04,
     &             .824184704952483D+05,.783031092490225D+06,
     &             .822210493622814D+07,.945557399360556D+08,
     &             .118195595640730D+10,.159564653040121D+11,
     &             .231369166433050D+12,.358622522796969D+13/
              Q2=1.414213562373095D0
              Q0=.3333333333333333D0
              Q1=.6666666666666667D0
              XE=X*DSQRT(X)/1.5D0
              XP6=1.0D0/DSQRT(6.0D0*PI*XE)
              SU1=1.0D0
              R=1.0D0
              XR1=1.0D0/XE
              DO 35 K=1,16
                 R=-R*XR1
35               SU1=SU1+A(K)*R
              SU2=1.0D0
              R=1.0D0
              DO 40 K=1,16
                 R=R*XR1
40               SU2=SU2+A(K)*R
              APT=Q0-DEXP(-XE)*XP6*SU1
              BPT=2.0D0*DEXP(XE)*XP6*SU2
              SU3=1.0D0
              R=1.0D0
              XR2=1.0D0/(XE*XE)
              DO 45 K=1,8
                 R=-R*XR2
45               SU3=SU3+A(2*K)*R
              SU4=A(1)*XR1
              R=XR1
              DO 50 K=1,7
                 R=-R*XR2
50               SU4=SU4+A(2*K+1)*R
              SU5=SU3+SU4
              SU6=SU3-SU4
              ANT=Q1-Q2*XP6*(SU5*DCOS(XE)-SU6*DSIN(XE))
              BNT=Q2*XP6*(SU5*DSIN(XE)+SU6*DCOS(XE))
           ENDIF
        ENDIF
        RETURN
        END
