        PROGRAM MSPHY

C       ========================================================
C       Purpose: This program computes the spherical Bessel
C                functions yn(x) and yn'(x) using subroutine
C                SPHY
C       Input :  x --- Argument of yn(x) ( x ע 0 )
C                n --- Order of yn(x) ( n = 0,1,תתת, ף 250 )
C       Output:  SY(n) --- yn(x)
C                DY(n) --- yn'(x)
C       Example:   x = 10.0
C                  n          yn(x)               yn'(x)
C                --------------------------------------------
C                  0     .8390715291D-01    -.6279282638D-01
C                  1     .6279282638D-01     .7134858763D-01
C                  2    -.6506930499D-01     .8231361788D-01
C                  3    -.9532747888D-01    -.2693831344D-01
C                  4    -.1659930220D-02    -.9449751377D-01
C                  5     .9383354168D-01    -.5796005523D-01
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SY(0:250),DY(0:250)
        WRITE(*,*)'Please enter n and x '
        READ(*,*)N,X
        WRITE(*,30)N,X
        IF (N.LE.10) THEN
           NS=1
        ELSE
           WRITE(*,*)'Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL SPHY(N,X,NM,SY,DY)
        WRITE(*,*)
        WRITE(*,*)'  n          yn(x)               yn''(x)'
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,SY(K),DY(K)
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,6HNmax =,I3,',     ',2Hx=,F6.1)
        END


        SUBROUTINE SPHY(N,X,NM,SY,DY)

C       ======================================================
C       Purpose: Compute spherical Bessel functions yn(x) and
C                their derivatives
C       Input :  x --- Argument of yn(x) ( x ע 0 )
C                n --- Order of yn(x) ( n = 0,1,תתת )
C       Output:  SY(n) --- yn(x)
C                DY(n) --- yn'(x)
C                NM --- Highest order computed
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SY(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              SY(K)=-1.0D+300
10            DY(K)=1.0D+300
           RETURN
        ENDIF
        SY(0)=-DCOS(X)/X
        SY(1)=(SY(0)-DSIN(X))/X
        F0=SY(0)
        F1=SY(1)
        DO 15 K=2,N
           F=(2.0D0*K-1.0D0)*F1/X-F0
           SY(K)=F
           IF (DABS(F).GE.1.0D+300) GO TO 20
           F0=F1
15         F1=F
20      NM=K-1
           DY(0)=(DSIN(X)+DCOS(X)/X)/X
           DO 25 K=1,NM
25            DY(K)=SY(K-1)-(K+1.0D0)*SY(K)/X
        RETURN
        END
