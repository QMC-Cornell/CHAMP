        PROGRAM MENXA

C       =========================================================
C       Purpose: This program computes the exponential integral
C                En(x) using subroutine ENXA
C       Example: x = 10.0
C                   n         En(x)
C                 ----------------------
C                   0     .45399930D-05
C                   1     .41569689D-05
C                   2     .38302405D-05
C                   3     .35487626D-05
C                   4     .33041014D-05
C                   5     .30897289D-05
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:100)
        WRITE(*,*)'Please enter n and x '
        READ(*,*)N,X
        WRITE(*,20)N,X
        WRITE(*,*)
        WRITE(*,*)'   n         En(x)'
        WRITE(*,*)' ----------------------'
        CALL ENXA(N,X,EN)
        DO 10 K=0,N
           WRITE(*,30)K,EN(K)
10      CONTINUE
20      FORMAT(5X,I3,',   ','x=',F5.1)
30      FORMAT(2X,I3,D18.8)
        END


        SUBROUTINE ENXA(N,X,EN)

C       ============================================
C       Purpose: Compute exponential integral En(x)
C       Input :  x --- Argument of En(x) ( x ó 20 )
C                n --- Order of En(x)
C       Output:  EN(n) --- En(x)
C       Routine called: E1XB for computing E1(x)
C       ============================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:N)
        EN(0)=DEXP(-X)/X
        CALL E1XB(X,E1)
        EN(1)=E1
        DO 10 K=2,N
           EK=(DEXP(-X)-X*E1)/(K-1.0D0)
           EN(K)=EK
10         E1=EK
        RETURN
        END


        SUBROUTINE E1XB(X,E1)

C       ============================================
C       Purpose: Compute exponential integral E1(x)
C       Input :  x  --- Argument of E1(x)
C       Output:  E1 --- E1(x)
C       ============================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0) THEN
           E1=1.0D+300
        ELSE IF (X.LE.1.0) THEN
           E1=1.0D0
           R=1.0D0
           DO 10 K=1,25
              R=-R*K*X/(K+1.0D0)**2
              E1=E1+R
              IF (DABS(R).LE.DABS(E1)*1.0D-15) GO TO 15
10         CONTINUE
15         GA=0.5772156649015328D0
           E1=-GA-DLOG(X)+X*E1
        ELSE
           M=20+INT(80.0/X)
           T0=0.0D0
           DO 20 K=M,1,-1
              T0=K/(1.0D0+K/(X+T0))
20         CONTINUE
           T=1.0D0/(X+T0)
           E1=DEXP(-X)*T
        ENDIF
        RETURN
        END
