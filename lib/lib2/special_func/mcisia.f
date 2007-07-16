        PROGRAM MCISIA

C       ========================================================
C       Purpose: This program computes the cosine and sine
C                integrals using subroutine CISIA
C       Input :  x  --- Argument of Ci(x) and Si(x)
C       Output:  CI --- Ci(x)
C                SI --- Si(x)
C       Example:
C                    x         Ci(x)          Si(x)
C                 ------------------------------------
C                   0.0     - ì             .00000000
C                   5.0     -.19002975     1.54993124
C                  10.0     -.04545643     1.65834759
C                  20.0      .04441982     1.54824170
C                  30.0     -.03303242     1.56675654
C                  40.0      .01902001     1.58698512
C       ========================================================

        DOUBLE PRECISION CI,SI,X
        WRITE(*,*)'Please enter x '
        READ(*,*) X
        WRITE(*,*)'   x         Ci(x)          Si(x)'
        WRITE(*,*)'------------------------------------'
        CALL CISIA(X,CI,SI)
        IF (X.NE.0.0D0) WRITE(*,10)X,CI,SI
        IF (X.EQ.0.0D0) WRITE(*,20)
10      FORMAT(1X,F5.1,2F15.8)
20      FORMAT(3X,' .0',4X,' - ì',13X,'.00000000')
        END



        SUBROUTINE CISIA(X,CI,SI)

C       =============================================
C       Purpose: Compute cosine and sine integrals
C                Si(x) and Ci(x)  ( x ò 0 )
C       Input :  x  --- Argument of Ci(x) and Si(x)
C       Output:  CI --- Ci(x)
C                SI --- Si(x)
C       =============================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(101)
        P2=1.570796326794897D0
        EL=.5772156649015329D0
        EPS=1.0D-15
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           CI=-1.0D+300
           SI=0.0D0
        ELSE IF (X.LE.16.0D0) THEN
           XR=-.25D0*X2
           CI=EL+DLOG(X)+XR
           DO 10 K=2,40
              XR=-.5D0*XR*(K-1)/(K*K*(2*K-1))*X2
              CI=CI+XR
              IF (DABS(XR).LT.DABS(CI)*EPS) GO TO 15
10         CONTINUE
15         XR=X
           SI=X
           DO 20 K=1,40
              XR=-.5D0*XR*(2*K-1)/K/(4*K*K+4*K+1)*X2
              SI=SI+XR
              IF (DABS(XR).LT.DABS(SI)*EPS) RETURN
20         CONTINUE
        ELSE IF (X.LE.32.0D0) THEN
           M=INT(47.2+.82*X)
           XA1=0.0D0
           XA0=1.0D-100
           DO 25 K=M,1,-1
              XA=4.0D0*K*XA0/X-XA1
              BJ(K)=XA
              XA1=XA0
25            XA0=XA
           XS=BJ(1)
           DO 30 K=3,M,2
30            XS=XS+2.0D0*BJ(K)
           BJ(1)=BJ(1)/XS
           DO 35 K=2,M
35            BJ(K)=BJ(K)/XS
           XR=1.0D0
           XG1=BJ(1)
           DO 40 K=2,M
              XR=.25D0*XR*(2.0*K-3.0)**2/((K-1.0)*(2.0*K-1.0)**2)*X
40            XG1=XG1+BJ(K)*XR
           XR=1.0D0
           XG2=BJ(1)
           DO 45 K=2,M
              XR=.25D0*XR*(2.0*K-5.0)**2/((K-1.0)*(2.0*K-3.0)**2)*X
45            XG2=XG2+BJ(K)*XR
           XCS=DCOS(X/2.0D0)
           XSS=DSIN(X/2.0D0)
           CI=EL+DLOG(X)-X*XSS*XG1+2*XCS*XG2-2*XCS*XCS
           SI=X*XCS*XG1+2*XSS*XG2-DSIN(X)
        ELSE
           XR=1.0D0
           XF=1.0D0
           DO 50 K=1,9
              XR=-2.0D0*XR*K*(2*K-1)/X2
50            XF=XF+XR
           XR=1.0D0/X
           XG=XR
           DO 55 K=1,8
              XR=-2.0D0*XR*(2*K+1)*K/X2
55            XG=XG+XR
           CI=XF*DSIN(X)/X-XG*DCOS(X)/X
           SI=P2-XF*DCOS(X)/X-XG*DSIN(X)/X
        ENDIF
        RETURN
        END

