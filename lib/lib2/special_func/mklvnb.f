        PROGRAM MKLVNB

C       ========================================================
C       Purpose: This program computes Kelvin functions ber x,
C                bei x, ker x and kei x, and their derivatives
C                using subroutine KLVNB
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
C     x       ber x         bei x         ker x         kei x
C   -------------------------------------------------------------
C     0    .100000D+01   .000000D+00    ì           -.785398D+00
C     5   -.623008D+01   .116034D+00  -.115117D-01   .111876D-01
C    10    .138840D+03   .563705D+02   .129466D-03  -.307525D-03
C    15   -.296725D+04  -.295271D+04  -.151433D-07   .796289D-05
C    20    .474894D+05   .114775D+06  -.771523D-07  -.185894D-06

C     x       ber'x         bei'x         ker'x         kei'x
C   -------------------------------------------------------------
C     0    .000000D+00   .000000D+00  - ì            .000000D+00
C     5   -.384534D+01  -.435414D+01   .171934D-01  -.819979D-03
C    10    .511952D+02   .135309D+03  -.315597D-03   .140914D-03
C    15    .910555D+02  -.408776D+04   .564468D-05  -.588222D-05
C    20   -.488032D+05   .111855D+06  -.750186D-07   .190624D-06
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x        ber x           bei x',
     &            '           ker x           kei x'
        WRITE(*,*)'--------------------------------',
     &            '--------------------------------------'
        CALL KLVNB(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
        WRITE(*,10)X,BER,BEI,GER,GEI
        WRITE(*,*)
        WRITE(*,*)'   x        ber''x           bei''x',
     &            '           ker''x           kei''x'
        WRITE(*,*)'--------------------------------',
     &            '--------------------------------------'
        WRITE(*,10)X,DER,DEI,HER,HEI
10      FORMAT(1X,F5.1,4D16.6)
        END


        SUBROUTINE KLVNB(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)

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
        IF (X.EQ.0.0D0) THEN
           BER=1.0D0
           BEI=0.0D0
           GER=1.0D+300
           GEI=-.25D0*PI
           DER=0.0D0
           DEI=0.0D0
           HER=-1.0D+300
           HEI=0.0D0
        ELSE IF (X.LT.8.0D0) THEN
           T=X/8.0D0
           T2=T*T
           U=T2*T2
           BER=((((((-.901D-5*U+.122552D-2)*U-.08349609D0)*U
     &         +2.64191397D0)*U-32.36345652D0)*U
     &         +113.77777774D0)*U-64.0D0)*U+1.0D0
           BEI=T*T*((((((.11346D-3*U-.01103667D0)*U
     &         +.52185615D0)*U-10.56765779D0)*U
     &         +72.81777742D0)*U-113.77777774D0)*U+16.0D0)
           GER=((((((-.2458D-4*U+.309699D-2)*U-.19636347D0)
     &         *U+5.65539121D0)*U-60.60977451D0)*U+
     &         171.36272133D0)*U-59.05819744D0)*U-.57721566D0
           GER=GER-DLOG(.5D0*X)*BER+.25D0*PI*BEI
           GEI=T2*((((((.29532D-3*U-.02695875D0)*U
     &         +1.17509064D0)*U-21.30060904D0)*U
     &         +124.2356965D0)*U-142.91827687D0)*U
     &         +6.76454936D0)
           GEI=GEI-DLOG(.5D0*X)*BEI-.25D0*PI*BER
           DER=X*T2*((((((-.394D-5*U+.45957D-3)*U
     &         -.02609253D0)*U+.66047849D0)*U-6.0681481D0)*U
     &         +14.22222222D0)*U-4.0D0)
           DEI=X*((((((.4609D-4*U-.379386D-2)*U+.14677204D0)
     &         *U-2.31167514D0)*U+11.37777772D0)*U
     &         -10.66666666D0)*U+.5D0)
           HER=X*T2*((((((-.1075D-4*U+.116137D-2)*U
     &         -.06136358D0)*U+1.4138478D0)*U-11.36433272D0)
     &         *U+21.42034017D0)*U-3.69113734D0)
           HER=HER-DLOG(.5D0*X)*DER-BER/X+.25D0*PI*DEI
           HEI=X*((((((.11997D-3*U-.926707D-2)*U
     &         +.33049424D0)*U-4.65950823D0)*U+19.41182758D0)
     &         *U-13.39858846D0)*U+.21139217D0)
           HEI=HEI-DLOG(.5D0*X)*DEI-BEI/X-.25D0*PI*DER
        ELSE
           T=8.0D0/X
           DO 10 L=1,2
              V=(-1)**L*T
              TPR=((((.6D-6*V-.34D-5)*V-.252D-4)*V-.906D-4)
     &            *V*V+.0110486D0)*V
              TPI=((((.19D-5*V+.51D-5)*V*V-.901D-4)*V
     &            -.9765D-3)*V-.0110485D0)*V-.3926991D0
              IF (L.EQ.1) THEN
                 TNR=TPR
                 TNI=TPI
              ENDIF
10         CONTINUE
           YD=X/DSQRT(2.0D0)
           YE1=DEXP(YD+TPR)
           YE2=DEXP(-YD+TNR)
           YC1=1.0D0/DSQRT(2.0D0*PI*X)
           YC2=DSQRT(PI/(2.0D0*X))
           CSP=DCOS(YD+TPI)
           SSP=DSIN(YD+TPI)
           CSN=DCOS(-YD+TNI)
           SSN=DSIN(-YD+TNI)
           GER=YC2*YE2*CSN
           GEI=YC2*YE2*SSN
           FXR=YC1*YE1*CSP
           FXI=YC1*YE1*SSP
           BER=FXR-GEI/PI
           BEI=FXI+GER/PI
           DO 15 L=1,2
              V=(-1)**L*T
              PPR=(((((.16D-5*V+.117D-4)*V+.346D-4)*V+.5D-6)
     &            *V-.13813D-2)*V-.0625001D0)*V+.7071068D0
              PPI=(((((-.32D-5*V-.24D-5)*V+.338D-4)*V+
     &           .2452D-3)*V+.13811D-2)*V-.1D-6)*V+.7071068D0
              IF (L.EQ.1) THEN
                 PNR=PPR
                 PNI=PPI
              ENDIF
15         CONTINUE
           HER=GEI*PNI-GER*PNR
           HEI=-(GEI*PNR+GER*PNI)
           DER=FXR*PPR-FXI*PPI-HEI/PI
           DEI=FXI*PPR+FXR*PPI+HER/PI
        ENDIF
        RETURN
        END
