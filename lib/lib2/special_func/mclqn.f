        PROGRAM MCLQN

C       ==========================================================
C       Purpose: This program computes the Legendre polynomials
C                Qn(z) and Qn'(z) for a complex argument using
C                subroutine CLQN
C       Input :  x --- Real part of z
C                y --- Imaginary part of z
C                n --- Degree of Qn(z), n = 0,1,...
C       Output:  CQN(n) --- Qn(z)
C                CQD(n) --- Qn'(z)
C       Examples:

C       z = 0.5 + 0.5 i
C       n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn'(z)]    Im[Qn'(z)]
C      -----------------------------------------------------------
C       0   .402359D+00   .553574D+00   .800000D+00   .400000D+00
C       1  -.107561D+01   .477967D+00   .602359D+00   .115357D+01
C       2  -.136636D+01  -.725018D+00  -.242682D+01   .183390D+01
C       3   .182619D+00  -.206146D+01  -.622944D+01  -.247151D+01
C       4   .298834D+01  -.110022D+01  -.114849D+01  -.125963D+02
C       5   .353361D+01   .334847D+01   .206656D+02  -.123735D+02

C       z = 3.0 + 2.0 i
C       n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn'(z)]    Im[Qn'(z)]
C      -----------------------------------------------------------
C       0   .229073D+00  -.160875D+00  -.250000D-01   .750000D-01
C       1   .896860D-02  -.244805D-01   .407268D-02   .141247D-01
C       2  -.736230D-03  -.281865D-02   .190581D-02   .155860D-02
C       3  -.264727D-03  -.227023D-03   .391535D-03   .314880D-04
C       4  -.430648D-04  -.443187D-05   .527190D-04  -.305592D-04
C       5  -.481362D-05   .265297D-05   .395108D-05  -.839883D-05
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CQN(0:100),CQD(0:100)
        WRITE(*,*)'  Please enter Nmax, x and y (z=x+iy)'
        READ(*,*)N,X,Y
        WRITE(*,30)X,Y
        WRITE(*,*)
        CALL CLQN(N,X,Y,CQN,CQD)
        WRITE(*,*)'  n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn''(z)]',
     &            '    Im[Qn''(z)]'
        WRITE(*,*)' ---------------------------------------------',
     &            '--------------'
        DO 10 K=0,N
10         WRITE(*,20)K,CQN(K),CQD(K)
20      FORMAT(1X,I3,4D14.6)
30      FORMAT(3X,'x =',F5.1,',  ','y =',F5.1)
        END


        SUBROUTINE CLQN(N,X,Y,CQN,CQD)

C       ==================================================
C       Purpose: Compute the Legendre functions Qn(z) and
C                their derivatives Qn'(z) for a complex
C                argument
C       Input :  x --- Real part of z
C                y --- Imaginary part of z
C                n --- Degree of Qn(z), n = 0,1,2,...
C       Output:  CQN(n) --- Qn(z)
C                CQD(n) --- Qn'(z)
C       ==================================================

        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CQN(0:N),CQD(0:N)
        Z=CMPLX(X,Y)
        IF (Z.EQ.1.0D0) THEN
           DO 10 K=0,N
              CQN(K)=(1.0D+300,0.0D0)
10            CQD(K)=(1.0D+300,0.0D0)
           RETURN
        ENDIF
        LS=1
        IF (CDABS(Z).GT.1.0D0) LS=-1
        CQ0=0.5D0*CDLOG(LS*(1.0D0+Z)/(1.0D0-Z))
        CQ1=Z*CQ0-1.0D0
        CQN(0)=CQ0
        CQN(1)=CQ1
        IF (CDABS(Z).LT.1.0001D0) THEN
           CQF0=CQ0
           CQF1=CQ1
           DO 15 K=2,N
              CQF2=((2.0D0*K-1.0D0)*Z*CQF1-(K-1.0D0)*CQF0)/K
              CQN(K)=CQF2
              CQF0=CQF1
15            CQF1=CQF2
        ELSE
           IF (CDABS(Z).GT.1.1D0) THEN
              KM=40+N
           ELSE
              KM=(40+N)*INT(-1.0-1.8*LOG(CDABS(Z-1.0)))
           ENDIF
           CQF2=0.0D0
           CQF1=1.0D0
           DO 20 K=KM,0,-1
              CQF0=((2*K+3.0D0)*Z*CQF1-(K+2.0D0)*CQF2)/(K+1.0D0)
              IF (K.LE.N) CQN(K)=CQF0
              CQF2=CQF1
20            CQF1=CQF0
           DO 25 K=0,N
25            CQN(K)=CQN(K)*CQ0/CQF0
        ENDIF
        CQD(0)=(CQN(1)-Z*CQN(0))/(Z*Z-1.0D0)
        DO 30 K=1,N
30         CQD(K)=(K*Z*CQN(K)-K*CQN(K-1))/(Z*Z-1.0D0)
        RETURN
        END
