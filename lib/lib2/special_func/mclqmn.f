        PROGRAM MCLQMN

C       ============================================================
C       Purpose: This program computes the associated Legendre
C                functions Qmn(z) and their derivatives Qmn'(z) for
C                a complex argument using subroutine CLQMN
C       Definition: Qmn(z)=(-1)**m*(1-z*z)**(m/2)*dm/dzm[Qn(z)]
C                   Q0(z)=1/2*LOG[(1+z)/(1-z)]     ( for |z|<1 )
C                   Qmn(z)=(z*z-1)**(m/2)*dm/dzm[Qn(z)]
C                   Q0(z)=1/2*LOG[(z+1)/(z-1)]     ( for |z|>1 )
C       Input :  x --- Real part of z
C                y --- Imaginary part of z
C                m --- Order of Qmn(z)  ( m = 0,1,2,תתת )
C                n --- Degree of Qmn(z) ( n = 0,1,2,תתת )
C       Output:  CQM(m,n) --- Qmn(z)
C                CQD(m,n) --- Qmn'(z)
C       Examples:
C                n = 5, x = 0.5, y = 0.2

C       m     Re[Qmn(z)]    Im[Qmn(z)]    Re[Qmn'(z)]   Im[Qmn'(z)]
C      -------------------------------------------------------------
C       0    .987156D+00   .354345D+00    .324023D+01  -.447297D+01
C       1   -.240328D+01   .436861D+01    .281158D+02   .171437D+02
C       2   -.245853D+02  -.138072D+02   -.106283D+03   .913792D+02
C       3    .102723D+03  -.651233D+02   -.362578D+03  -.429802D+03
C       4    .155510D+03   .357712D+03    .196975D+04  -.287414D+02
C       5   -.167357D+04  -.680954D+03   -.193093D+04  -.925757D+03

C                n = 5, x = 2.5, y = 1.0

C       m     Re[Qmn(z)]    Im[Qmn(z)]    Re[Qmn'(z)]   Im[Qmn'(z)]
C      -------------------------------------------------------------
C       0   -.274023D-04  -.227141D-04    .809834D-04   .210884D-04
C       1    .165620D-03   .136108D-03   -.489095D-03  -.124400D-03
C       2   -.118481D-02  -.948832D-03    .349090D-02   .825057D-03
C       3    .982179D-02   .753264D-02   -.288271D-01  -.596384D-02
C       4   -.927915D-01  -.669521D-01    .270840D+00   .451376D-01
C       5    .985601D+00   .656737D+00   -.285567D+01  -.332533D+00
C       ============================================================

        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CQM(0:40,0:40),CQD(0:40,0:40)
        WRITE(*,*)'  Please enter m, n, x and y '
        READ(*,*) M,N,X,Y
        WRITE(*,30)M,N,X,Y
        CALL CLQMN(40,M,N,X,Y,CQM,CQD)
        WRITE(*,*)'   m   n   Re[Qmn(z)]    Im[Qmn(z)]    ',
     &            'Re[Qmn''(z)]   Im[Qmn''(z)]'
        WRITE(*,*)' -----------------------------------',
     &            '------------------------------'
        DO 10 J=0,N
10         WRITE(*,20)M,J,CQM(M,J),CQD(M,J)
20      FORMAT(1X,2I4,2D14.6,1X,2D14.6)
30      FORMAT(1X,'m =',I2,', ','n =',I2,', ','x =',F4.1,
     &         ', ','y =',F4.1)
        END


        SUBROUTINE CLQMN(MM,M,N,X,Y,CQM,CQD)

C       =======================================================
C       Purpose: Compute the associated Legendre functions of
C                the second kind, Qmn(z) and Qmn'(z), for a
C                complex argument
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                m  --- Order of Qmn(z)  ( m = 0,1,2,תתת )
C                n  --- Degree of Qmn(z) ( n = 0,1,2,תתת )
C                mm --- Physical dimension of CQM and CQD
C       Output:  CQM(m,n) --- Qmn(z)
C                CQD(m,n) --- Qmn'(z)
C       =======================================================

        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CQM(0:MM,0:N),CQD(0:MM,0:N)
        Z=CMPLX(X,Y)
        IF (DABS(X).EQ.1.0D0.AND.Y.EQ.0.0D0) THEN
           DO 10 I=0,M
           DO 10 J=0,N
              CQM(I,J)=(1.0D+300,0.0D0)
              CQD(I,J)=(1.0D+300,0.0D0)
10         CONTINUE
           RETURN
        ENDIF
        XC=CDABS(Z)
        IF (DIMAG(Z).EQ.0.0D0.OR.XC.LT.1.0D0) LS=1
        IF (XC.GT.1.0D0) LS=-1
        ZQ=CDSQRT(LS*(1.0D0-Z*Z))
        ZS=LS*(1.0D0-Z*Z)
        CQ0=0.5D0*CDLOG(LS*(1.0D0+Z)/(1.0D0-Z))
        IF (XC.LT.1.0001D0) THEN
           CQM(0,0)=CQ0
           CQM(0,1)=Z*CQ0-1.0D0
           CQM(1,0)=-1.0D0/ZQ
           CQM(1,1)=-ZQ*(CQ0+Z/(1.0D0-Z*Z))
           DO 15 I=0,1
           DO 15 J=2,N
              CQM(I,J)=((2.0D0*J-1.0D0)*Z*CQM(I,J-1)
     &                -(J+I-1.0D0)*CQM(I,J-2))/(J-I)
15         CONTINUE
           DO 20 J=0,N
           DO 20 I=2,M
              CQM(I,J)=-2.0D0*(I-1.0D0)*Z/ZQ*CQM(I-1,J)-LS*
     &                 (J+I-1.0D0)*(J-I+2.0D0)*CQM(I-2,J)
20         CONTINUE
        ELSE
           IF (XC.GT.1.1) THEN
              KM=40+M+N
           ELSE
              KM=(40+M+N)*INT(-1.0-1.8*LOG(XC-1.0))
           ENDIF
           CQF2=(0.0D0,0.0D0)
           CQF1=(1.0D0,0.0D0)
           DO 25 K=KM,0,-1
              CQF0=((2*K+3.0D0)*Z*CQF1-(K+2.0D0)*CQF2)/(K+1.0D0)
              IF (K.LE.N) CQM(0,K)=CQF0
              CQF2=CQF1
25            CQF1=CQF0
           DO 30 K=0,N
30            CQM(0,K)=CQ0*CQM(0,K)/CQF0
           CQF2=0.0D0
           CQF1=1.0D0
           DO 35 K=KM,0,-1
              CQF0=((2*K+3.0D0)*Z*CQF1-(K+1.0D0)*CQF2)/(K+2.0D0)
              IF (K.LE.N) CQM(1,K)=CQF0
              CQF2=CQF1
35            CQF1=CQF0
           CQ10=-1.0D0/ZQ
           DO 40 K=0,N
40            CQM(1,K)=CQ10*CQM(1,K)/CQF0
           DO 45 J=0,N
              CQ0=CQM(0,J)
              CQ1=CQM(1,J)
              DO 45 I=0,M-2
                 CQF=-2.0D0*(I+1)*Z/ZQ*CQ1+(J-I)*(J+I+1.0D0)*CQ0
                 CQM(I+2,J)=CQF
                 CQ0=CQ1
                 CQ1=CQF
45         CONTINUE
        ENDIF
        CQD(0,0)=LS/ZS
        DO 50 J=1,N
50         CQD(0,J)=LS*J*(CQM(0,J-1)-Z*CQM(0,J))/ZS
        DO 55 J=0,N
        DO 55 I=1,M
           CQD(I,J)=LS*I*Z/ZS*CQM(I,J)+(I+J)*(J-I+1.0D0)
     &              /ZQ*CQM(I-1,J)
55      CONTINUE
        RETURN
        END


