        PROGRAM MLPMNS

C       ========================================================
C       Purpose: This program computes the associated Legendre
C                functions Pmn(x) and their derivatives Pmn'(x)
C                for a given order using subroutine LPMNS
C       Input :  x --- Argument of Pmn(x)
C                m --- Order of Pmn(x),  m = 0,1,2,...,n
C                n --- Degree of Pmn(x), n = 0,1,2,...,N
C       Output:  PM(n) --- Pmn(x)
C                PD(n) --- Pmn'(x)
C       Examples:
C                m = 1,  N = 5,  x = .5
C                n        Pmn(x)           Pmn'(x)
C               -------------------------------------
C                0    .00000000D+00    .00000000D+00
C                1    .86602540D+00   -.57735027D+00
C                2    .12990381D+01    .17320508D+01
C                3    .32475953D+00    .62786842D+01
C                4   -.13531647D+01    .57735027D+01
C                5   -.19282597D+01   -.43977853D+01

C                m = 2,  N = 6,  x = 2.5
C                n        Pmn(x)           Pmn'(x)
C               -------------------------------------
C                0    .00000000D+00    .00000000D+00
C                1    .00000000D+00    .00000000D+00
C                2    .15750000D+02    .15000000D+02
C                3    .19687500D+03    .26625000D+03
C                4    .16832813D+04    .29812500D+04
C                5    .12230859D+05    .26876719D+05
C                6    .81141416D+05    .21319512D+06
C       =======================================================

        IMPLICIT DOUBLE PRECISION (P,X,Y)
        DIMENSION PM(0:200),PD(0:200)
        WRITE(*,*)'Please enter m, N, and x '
        READ(*,*)M,N,X
        WRITE(*,30)M,N,X
        CALL LPMNS(M,N,X,PM,PD)
        WRITE(*,*)
        WRITE(*,*)'  n        Pmn(x)           Pmn''(x)    '
        WRITE(*,*)' -------------------------------------'
        DO 10 J=0,N
        WRITE(*,20)J,PM(J),PD(J)
10      CONTINUE
20      FORMAT(1X,I3,2D17.8)
30      FORMAT(1X,'m =',I2,',  ','n =',I2,',  ','x =',F5.1)
        END


        SUBROUTINE LPMNS(M,N,X,PM,PD)

C       ========================================================
C       Purpose: Compute associated Legendre functions Pmn(x)
C                and Pmn'(x) for a given order
C       Input :  x --- Argument of Pmn(x)
C                m --- Order of Pmn(x),  m = 0,1,2,...,n
C                n --- Degree of Pmn(x), n = 0,1,2,...,N
C       Output:  PM(n) --- Pmn(x)
C                PD(n) --- Pmn'(x)
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION PM(0:N),PD(0:N)
        DO 10 K=0,N
           PM(K)=0.0D0
10         PD(K)=0.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 K=0,N
              IF (M.EQ.0) THEN
                 PM(K)=1.0D0
                 PD(K)=0.5D0*K*(K+1.0)
                 IF (X.LT.0.0) THEN
                    PM(K)=(-1)**K*PM(K)
                    PD(K)=(-1)**(K+1)*PD(K)
                 ENDIF
              ELSE IF (M.EQ.1) THEN
                 PD(K)=1.0D+300
              ELSE IF (M.EQ.2) THEN
                 PD(K)=-0.25D0*(K+2.0)*(K+1.0)*K*(K-1.0)
                 IF (X.LT.0.0) PD(K)=(-1)**(K+1)*PD(K)
              ENDIF
15         CONTINUE
           RETURN
        ENDIF
        X0=DABS(1.0D0-X*X)
        PM0=1.0D0
        PMK=PM0
        DO 20 K=1,M
           PMK=(2.0D0*K-1.0D0)*DSQRT(X0)*PM0
20         PM0=PMK
        PM1=(2.0D0*M+1.0D0)*X*PM0
        PM(M)=PMK
        PM(M+1)=PM1
        DO 25 K=M+2,N
           PM2=((2.0D0*K-1.0D0)*X*PM1-(K+M-1.0D0)*PMK)/(K-M)
           PM(K)=PM2
           PMK=PM1
25         PM1=PM2
        PD(0)=((1.0D0-M)*PM(1)-X*PM(0))/(X*X-1.0)
        DO 30 K=1,N
30          PD(K)=(K*X*PM(K)-(K+M)*PM(K-1))/(X*X-1.0D0)
        RETURN
        END
