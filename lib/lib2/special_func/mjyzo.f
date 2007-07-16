        PROGRAM MJYZO

C       ==========================================================
C       Purpose: This program computes the zeros of Bessel
C                functions Jn(x), Yn(x), and their derivatives
C                using subroutine JYZO
C       Input :  n --- Order of Bessel functions ( n ó 100 )
C                NT --- Number of zeros
C       Output:  RJ0(m) --- m-th zero of Jn(x),  m=1,2,...,NT
C                RJ1(m) --- m-th zero of Jn'(x), m=1,2,...,NT
C                RY0(m) --- m-th zero of Yn(x),  m=1,2,...,NT
C                RY1(m) --- m-th zero of Yn'(x), m=1,2,...,NT
C       Example: n = 1, NT =5

C      Zeros of Bessel funcions Jn(x), Yn(x) and their derivatives
C                                 ( n = 1 )
C       m       jnm           j'nm          ynm           y'nm
C      -----------------------------------------------------------
C       1     3.8317060     1.8411838     2.1971413     3.6830229
C       2     7.0155867     5.3314428     5.4296810     6.9415000
C       3    10.1734681     8.5363164     8.5960059    10.1234047
C       4    13.3236919    11.7060049    11.7491548    13.2857582
C       5    16.4706301    14.8635886    14.8974421    16.4400580
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RJ0(101),RJ1(101),RY0(101),RY1(101)
        WRITE(*,*)'Please enter n and NT '
        READ(*,*)N,NT
        WRITE(*,*)
        CALL JYZO(N,NT,RJ0,RJ1,RY0,RY1)
        WRITE(*,30)
        WRITE(*,40)N
        WRITE(*,*)'  m       jnm           j''nm          ynm',
     &            '           y''nm'
        WRITE(*,*)' ----------------------------------------',
     &            '-------------------'
        DO 10 M=1,NT
10         WRITE(*,50)M,RJ0(M),RJ1(M),RY0(M),RY1(M)
30      FORMAT(2X,'Zeros of Bessel funcions Jn(x), Yn(x)',
     &         ' and their derivatives')
40      FORMAT(30X,'( n =',I2,' )')
50      FORMAT(1X,I3,4F14.7)
        END


        SUBROUTINE JYZO(N,NT,RJ0,RJ1,RY0,RY1)

C       ======================================================
C       Purpose: Compute the zeros of Bessel functions Jn(x),
C                Yn(x), and their derivatives
C       Input :  n  --- Order of Bessel functions ( n ó 101 )
C                NT --- Number of zeros (roots)
C       Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
C                RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
C                RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
C                RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
C       Routine called: JYNDD for computing Jn(x), Yn(x), and
C                       their first and second derivatives
C       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RJ0(NT),RJ1(NT),RY0(NT),RY1(NT)
        IF (N.LE.20) THEN
           X=2.82141+1.15859*N
        ELSE
           X=N+1.85576*N**0.33333+1.03315/N**0.33333
        ENDIF
        L=0
10      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-BJN/DJN
        IF (DABS(X-X0).GT.1.0D-9) GO TO 10
        L=L+1
        RJ0(L)=X
        X=X+3.1416+(0.0972+0.0679*N-0.000354*N**2)/L
        IF (L.LT.NT) GO TO 10
        IF (N.LE.20) THEN
           X=0.961587+1.07703*N
        ELSE
           X=N+0.80861*N**0.33333+0.07249/N**0.33333
        ENDIF
        IF (N.EQ.0) X=3.8317
        L=0
15      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-DJN/FJN
        IF (DABS(X-X0).GT.1.0D-9) GO TO 15
        L=L+1
        RJ1(L)=X
        X=X+3.1416+(0.4955+0.0915*N-0.000435*N**2)/L
        IF (L.LT.NT) GO TO 15
        IF (N.LE.20) THEN
           X=1.19477+1.08933*N
        ELSE
           X=N+0.93158*N**0.33333+0.26035/N**0.33333
        ENDIF
        L=0
20      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-BYN/DYN
        IF (DABS(X-X0).GT.1.0D-9) GO TO 20
        L=L+1
        RY0(L)=X
        X=X+3.1416+(0.312+0.0852*N-0.000403*N**2)/L
        IF (L.LT.NT) GO TO 20
        IF (N.LE.20) THEN
           X=2.67257+1.16099*N
        ELSE
           X=N+1.8211*N**0.33333+0.94001/N**0.33333
        ENDIF
        L=0
25      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-DYN/FYN
        IF (DABS(X-X0).GT.1.0D-9) GO TO 25
        L=L+1
        RY1(L)=X
        X=X+3.1416+(0.197+0.0643*N-0.000286*N**2)/L
        IF (L.LT.NT) GO TO 25
        RETURN
        END


        SUBROUTINE JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)

C       ===========================================================
C       Purpose: Compute Bessel functions Jn(x) and Yn(x), and
C                their first and second derivatives
C       Input:   x   ---  Argument of Jn(x) and Yn(x) ( x > 0 )
C                n   ---  Order of Jn(x) and Yn(x)
C       Output:  BJN ---  Jn(x)
C                DJN ---  Jn'(x)
C                FJN ---  Jn"(x)
C                BYN ---  Yn(x)
C                DYN ---  Yn'(x)
C                FYN ---  Yn"(x)
C       ===========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(102),BY(102)
        DO 10 NT=1,900
           MT=INT(0.5*LOG10(6.28*NT)-NT*LOG10(1.36*DABS(X)/NT))
           IF (MT.GT.20) GO TO 15
10      CONTINUE
15      M=NT
        BS=0.0D0
        F0=0.0D0
        F1=1.0D-35
        SU=0.0D0
        DO 20 K=M,0,-1
           F=2.0D0*(K+1.0D0)*F1/X-F0
           IF (K.LE.N+1) BJ(K+1)=F
           IF (K.EQ.2*INT(K/2)) THEN
              BS=BS+2.0D0*F
              IF (K.NE.0) SU=SU+(-1)**(K/2)*F/K
           ENDIF
           F0=F1
20         F1=F
        DO 25 K=0,N+1
25         BJ(K+1)=BJ(K+1)/(BS-F)
        BJN=BJ(N+1)
        EC=0.5772156649015329D0
        E0=0.3183098861837907D0
        S1=2.0D0*E0*(DLOG(X/2.0D0)+EC)*BJ(1)
        F0=S1-8.0D0*E0*SU/(BS-F)
        F1=(BJ(2)*F0-2.0D0*E0/X)/BJ(1)
        BY(1)=F0
        BY(2)=F1
        DO 30 K=2,N+1
           F=2.0D0*(K-1.0D0)*F1/X-F0
           BY(K+1)=F
           F0=F1
30         F1=F
        BYN=BY(N+1)
        DJN=-BJ(N+2)+N*BJ(N+1)/X
        DYN=-BY(N+2)+N*BY(N+1)/X
        FJN=(N*N/(X*X)-1.0D0)*BJN-DJN/X
        FYN=(N*N/(X*X)-1.0D0)*BYN-DYN/X
        RETURN
        END



