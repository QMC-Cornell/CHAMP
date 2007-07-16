        PROGRAM MCLPN

C       ==========================================================
C       Purpose: This program computes the Legendre polynomials
C                Pn(z) and Pn'(z) for a complex argument using
C                subroutine CLPN
C       Input :  x --- Real part of z
C                y --- Imaginary part of z
C                n --- Degree of Pn(z), n = 0,1,...,N
C       Output:  CPN(n) --- Pn(z)
C                CPD(n) --- Pn'(z)
C       Example: z = 3.0 +2.0 i

C       n    Re[Pn(z)]     Im[Pn(z)]     Re[Pn'(z)]   Im[Pn'(z)]
C      -----------------------------------------------------------
C       0   .100000D+01   .000000D+00   .000000D+00   .000000D+00
C       1   .300000D+01   .200000D+01   .100000D+01   .000000D+00
C       2   .700000D+01   .180000D+02   .900000D+01   .600000D+01
C       3  -.270000D+02   .112000D+03   .360000D+02   .900000D+02
C       4  -.539000D+03   .480000D+03  -.180000D+03   .790000D+03
C       5  -.461700D+04   .562000D+03  -.481500D+04   .441000D+04
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION CPN(0:100),CPD(0:100)
        WRITE(*,*)'  Please enter Nmax, x and y (z=x+iy)'
        READ(*,*)N,X,Y
        WRITE(*,30)X,Y
        WRITE(*,*)
        CALL CLPN(N,X,Y,CPN,CPD)
        WRITE(*,*)'  n    Re[Pn(z)]     Im[Pn(z)]     Re[Pn''(z)]',
     &            '   Im[Pn''(z)]'
        WRITE(*,*)' ---------------------------------------------',
     &            '--------------'
        DO 10 K=0,N
10         WRITE(*,20)K,CPN(K),CPD(K)
20      FORMAT(1X,I3,4D14.6)
30      FORMAT(3X,'x =',F5.1,',  ','y =',F5.1)
        END


        SUBROUTINE CLPN(N,X,Y,CPN,CPD)

C       ==================================================
C       Purpose: Compute Legendre polynomials Pn(z) and
C                their derivatives Pn'(z) for a complex
C                argument
C       Input :  x --- Real part of z
C                y --- Imaginary part of z
C                n --- Degree of Pn(z), n = 0,1,2,...
C       Output:  CPN(n) --- Pn(z)
C                CPD(n) --- Pn'(z)
C       ==================================================

        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION CPN(0:N),CPD(0:N)
        Z=CMPLX(X,Y)
        CPN(0)=(1.0D0,0.0D0)
        CPN(1)=Z
        CPD(0)=(0.0D0,0.0D0)
        CPD(1)=(1.0D0,0.0D0)
        CP0=(1.0D0,0.0D0)
        CP1=Z
        DO 10 K=2,N
           CPF=(2.0D0*K-1.0D0)/K*Z*CP1-(K-1.0D0)/K*CP0
           CPN(K)=CPF
           IF (DABS(X).EQ.1.0D0.AND.Y.EQ.0.0D0) THEN
              CPD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              CPD(K)=K*(CP1-Z*CPF)/(1.0D0-Z*Z)
           ENDIF
           CP0=CP1
10         CP1=CPF
        RETURN
        END
