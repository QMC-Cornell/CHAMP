        PROGRAM MFCS

C       =======================================================
C       Purpose: This program computes the Fresnel integrals
C                C(x) and S(x) using subroutine FCS
C       Input :  x --- Argument of C(x) and S(x)
C       Output:  C --- C(x)
C                S --- S(x)
C       Example:
C                  x          C(x)          S(x)
C                -----------------------------------
C                 0.0      .00000000      .00000000
C                 0.5      .49234423      .06473243
C                 1.0      .77989340      .43825915
C                 1.5      .44526118      .69750496
C                 2.0      .48825341      .34341568
C                 2.5      .45741301      .61918176
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*) X
        WRITE(*,*)'   x          C(x)          S(x)'
        WRITE(*,*)' -----------------------------------'
        CALL FCS(X,C,S)
        WRITE(*,10)X,C,S
10      FORMAT(1X,F5.1,2F15.8)
        END


        SUBROUTINE FCS(X,C,S)

C       =================================================
C       Purpose: Compute Fresnel integrals C(x) and S(x)
C       Input :  x --- Argument of C(x) and S(x)
C       Output:  C --- C(x)
C                S --- S(x)
C       =================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        XA=DABS(X)
        PX=PI*XA
        T=.5D0*PX*XA
        T2=T*T
        IF (XA.EQ.0.0) THEN
           C=0.0D0
           S=0.0D0
        ELSE IF (XA.LT.2.5D0) THEN
           R=XA
           C=R
           DO 10 K=1,50
              R=-.5D0*R*(4.0D0*K-3.0D0)/K/(2.0D0*K-1.0D0)
     &          /(4.0D0*K+1.0D0)*T2
              C=C+R
              IF (DABS(R).LT.DABS(C)*EPS) GO TO 15
10         CONTINUE
15         S=XA*T/3.0D0
           R=S
           DO 20 K=1,50
              R=-.5D0*R*(4.0D0*K-1.0D0)/K/(2.0D0*K+1.0D0)
     &          /(4.0D0*K+3.0D0)*T2
              S=S+R
              IF (DABS(R).LT.DABS(S)*EPS) GO TO 40
20         CONTINUE
        ELSE IF (XA.LT.4.5D0) THEN
           M=INT(42.0+1.75*T)
           SU=0.0D0
           C=0.0D0
           S=0.0D0
           F1=0.0D0
           F0=1.0D-100
           DO 25 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F0/T-F1
              IF (K.EQ.INT(K/2)*2) THEN
                 C=C+F
              ELSE
                 S=S+F
              ENDIF
              SU=SU+(2.0D0*K+1.0D0)*F*F
              F1=F0
25            F0=F
           Q=DSQRT(SU)
           C=C*XA/Q
           S=S*XA/Q
        ELSE
           R=1.0D0
           F=1.0D0
           DO 30 K=1,20
              R=-.25D0*R*(4.0D0*K-1.0D0)*(4.0D0*K-3.0D0)/T2
30            F=F+R
           R=1.0D0/(PX*XA)
           G=R
           DO 35 K=1,12
              R=-.25D0*R*(4.0D0*K+1.0D0)*(4.0D0*K-1.0D0)/T2
35            G=G+R
           T0=T-INT(T/(2.0D0*PI))*2.0D0*PI
           C=.5D0+(F*DSIN(T0)-G*DCOS(T0))/PX
           S=.5D0-(F*DCOS(T0)+G*DSIN(T0))/PX
        ENDIF
40      IF (X.LT.0.0D0) THEN
           C=-C
           S=-S
        ENDIF
        RETURN
        END
