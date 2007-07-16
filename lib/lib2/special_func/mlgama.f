        PROGRAM MLGAMA

C       ===================================================
C       Purpose: This program computes the gamma function
C                â(x) for x > 0 using subroutine LGAMA
C       Examples:
C                  x           â(x)
C                -------------------------
C                 0.5     .1772453851D+01
C                 2.5     .1329340388D+01
C                 5.0     .2400000000D+02
C                 7.5     .1871254306D+04
C                10.0     .3628800000D+06
C       ===================================================

        IMPLICIT DOUBLE PRECISION (G,X)
        WRITE(*,*)'   x           â(x)'
        WRITE(*,*)' -------------------------'
        DO 10 L=0,20,5
           X=0.5D0*L
           IF (L.EQ.0) X=0.5
           CALL LGAMA(1,X,GL)
           WRITE(*,20)X,GL
10      CONTINUE
        WRITE(*,*) 'Please enter x:'
        READ(*,*) X
        CALL LGAMA(1,X,GL)
        WRITE(*,20)X,GL
20      FORMAT(1X,F5.1,D20.10)
        END


        SUBROUTINE LGAMA(KF,X,GL)

C       ==================================================
C       Purpose: Compute gamma function â(x) or ln[â(x)]
C       Input:   x  --- Argument of â(x) ( x > 0 )
C                KF --- Function code
C                       KF=1 for â(x); KF=0 for ln[â(x)]
C       Output:  GL --- â(x) or ln[â(x)]
C       ==================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/
        X0=X
        IF (X.EQ.1.0.OR.X.EQ.2.0) THEN
           GL=0.0D0
           GO TO 20
        ELSE IF (X.LE.7.0) THEN
           N=INT(7-X)
           X0=X+N
        ENDIF
        X2=1.0D0/(X0*X0)
        XP=6.283185307179586477D0
        GL0=A(10)
        DO 10 K=9,1,-1
10         GL0=GL0*X2+A(K)
        GL=GL0/X0+0.5D0*DLOG(XP)+(X0-.5D0)*DLOG(X0)-X0
        IF (X.LE.7.0) THEN
           DO 15 K=1,N
              GL=GL-DLOG(X0-1.0D0)
15            X0=X0-1.0D0
        ENDIF
20      IF (KF.EQ.1) GL=DEXP(GL)
        RETURN
        END
