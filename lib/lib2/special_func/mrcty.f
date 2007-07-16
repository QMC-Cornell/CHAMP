        PROGRAM MRCTY

C       =======================================================
C       Purpose: This program computes the Riccati-Bessel
C                functions of the second kind and their
C                derivatives using subroutine RCTY
C       Input:   x --- Argument of Riccati-Bessel function
C                n --- Order of yn(x)
C       Output:  RY(n) --- xúyn(x)
C                DY(n) --- [xúyn(x)]'
C       Example: x = 10.0
C                  n        xúyn(x)             [xúyn(x)]'
C                --------------------------------------------
C                  0     .8390715291D+00    -.5440211109D+00
C                  1     .6279282638D+00     .7762787027D+00
C                  2    -.6506930499D+00     .7580668738D+00
C                  3    -.9532747888D+00    -.3647106133D+00
C                  4    -.1659930220D-01    -.9466350679D+00
C                  5     .9383354168D+00    -.4857670106D+00
C       =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RY(0:250),DY(0:250)
        WRITE(*,*)'  Please enter n and x '
        READ(*,*)N,X
        WRITE(*,30)N,X
        IF (N.LE.10) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        WRITE(*,*)
        CALL RCTY(N,X,NM,RY,DY)
        WRITE(*,*)
        WRITE(*,*)'  n        xúyn(x)             [xúyn(x)]'''
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
           WRITE(*,20)K,RY(K),DY(K)
10      CONTINUE
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,6HNmax =,I3,',    ',3Hx =,F6.2)
        END


        SUBROUTINE RCTY(N,X,NM,RY,DY)

C       ========================================================
C       Purpose: Compute Riccati-Bessel functions of the second
C                kind and their derivatives
C       Input:   x --- Argument of Riccati-Bessel function
C                n --- Order of yn(x)
C       Output:  RY(n) --- xúyn(x)
C                DY(n) --- [xúyn(x)]'
C                NM --- Highest order computed
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RY(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              RY(K)=-1.0D+300
10            DY(K)=1.0D+300
           RY(0)=-1.0D0
           DY(0)=0.0D0
           RETURN
        ENDIF
        RY(0)=-DCOS(X)
        RY(1)=RY(0)/X-DSIN(X)
        RF0=RY(0)
        RF1=RY(1)
        DO 15 K=2,N
           RF2=(2.0D0*K-1.0D0)*RF1/X-RF0
           IF (DABS(RF2).GT.1.0D+300) GO TO 20
           RY(K)=RF2
           RF0=RF1
15         RF1=RF2
20      NM=K-1
        DY(0)=DSIN(X)
        DO 25 K=1,NM
25         DY(K)=-K*RY(K)/X+RY(K-1)
        RETURN
        END
