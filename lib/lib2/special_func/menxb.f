        PROGRAM MENXB

C       =========================================================
C       Purpose: This program computes the exponential integral
C                En(x) using subroutine ENXB
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
        CALL ENXB(N,X,EN)
        DO 10 K=0,N
           WRITE(*,30)K,EN(K)
10      CONTINUE
20      FORMAT(5X,I3,',   ','x=',F5.1)
30      FORMAT(2X,I3,D18.8)
        END


        SUBROUTINE ENXB(N,X,EN)

C       ===============================================
C       Purpose: Compute exponential integral En(x)
C       Input :  x --- Argument of En(x)
C                n --- Order of En(x)  (n = 0,1,2,...)
C       Output:  EN(n) --- En(x)
C       ===============================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:N)
        IF (X.EQ.0.0) THEN
           EN(0)=1.0D+300
           EN(1)=1.0D+300
           DO 10 K=2,N
10            EN(K)=1.0D0/(K-1.0)
           RETURN
        ELSE IF (X.LE.1.0) THEN
           EN(0)=DEXP(-X)/X
           DO 40 L=1,N
              RP=1.0D0
              DO 15 J=1,L-1
15               RP=-RP*X/J
              PS=-0.5772156649015328D0
              DO 20 M=1,L-1
20               PS=PS+1.0D0/M
              ENS=RP*(-DLOG(X)+PS)
              S=0.0D0
              DO 30 M=0,20
                 IF (M.EQ.L-1) GO TO 30
                 R=1.0D0
                 DO 25 J=1,M
25                  R=-R*X/J
                 S=S+R/(M-L+1.0D0)
                 IF (DABS(S-S0).LT.DABS(S)*1.0D-15) GO TO 35
                 S0=S
30            CONTINUE
35            EN(L)=ENS-S
40         CONTINUE
        ELSE
           EN(0)=DEXP(-X)/X
           M=15+INT(100.0/X)
           DO 50 L=1,N
              T0=0.0D0
              DO 45 K=M,1,-1
45               T0=(L+K-1.0D0)/(1.0D0+K/(X+T0))
              T=1.0D0/(X+T0)
50            EN(L)=DEXP(-X)*T
        ENDIF
        END
