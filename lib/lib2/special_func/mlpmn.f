        PROGRAM MLPMN

C       ==========================================================
C       Purpose: This program computes the associated Legendre
C                functions Pmn(x) and their derivatives Pmn'(x)
C                using subroutine LPMN
C       Input :  x --- Argument of Pmn(x)
C                m --- Order of Pmn(x),  m = 0,1,2,...,n
C                n --- Degree of Pmn(x), n = 0,1,2,...,N
C       Output:  PM(m,n) --- Pmn(x)
C                PD(m,n) --- Pmn'(x)
C       Example: x = 0.50
C          Pmn(x):
C          m\n        1            2            3            4
C         --------------------------------------------------------
C           0      .500000     -.125000     -.437500     -.289063
C           1     -.866025    -1.299038     -.324760     1.353165
C           2      .000000     2.250000     5.625000     4.218750
C           3      .000000      .000000    -9.742786   -34.099750
C           4      .000000      .000000      .000000    59.062500

C          Pmn'(x):
C          m\n        1            2            3            4
C         --------------------------------------------------------
C           0     1.000000     1.500000      .375000    -1.562500
C           1      .577350    -1.732051    -6.278684    -5.773503
C           2      .000000    -3.000000     3.750000    33.750000
C           3      .000000      .000000    19.485572      .000000
C           4      .000000      .000000      .000000  -157.500000
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PM(0:100,0:100),PD(0:100,0:100)
        WRITE(*,*)'  Please enter m, n and x'
        READ(*,*) M,N,X
        WRITE(*,*)
        WRITE(*,*)'  m     n      x          Pmn(x)         Pmn''(x)'
        WRITE(*,*)' ---------------------------------------------------'
        CALL LPMN(100,M,N,X,PM,PD)
        DO 15 J=0,N
           WRITE(*,10)M,J,X,PM(M,J),PD(M,J)
15      CONTINUE
10      FORMAT(1X,I3,3X,I3,3X,F5.1,2E17.8)
        END


        SUBROUTINE LPMN(MM,M,N,X,PM,PD)

C       =====================================================
C       Purpose: Compute the associated Legendre functions
C                Pmn(x) and their derivatives Pmn'(x)
C       Input :  x  --- Argument of Pmn(x)
C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
C                mm --- Physical dimension of PM and PD
C       Output:  PM(m,n) --- Pmn(x)
C                PD(m,n) --- Pmn'(x)
C       =====================================================

        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PM(0:MM,0:N),PD(0:MM,0:N)
        DO 10 I=0,N
        DO 10 J=0,M
           PM(J,I)=0.0D0
10         PD(J,I)=0.0D0
        PM(0,0)=1.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
              PM(0,I)=X**I
15            PD(0,I)=0.5D0*I*(I+1.0D0)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 PD(I,J)=1.0D+300
              ELSE IF (I.EQ.2) THEN
                 PD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1)
        DO 35 I=0,M
35         PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)
        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)-
     &             (I+J-1.0D0)*PM(I,J-2))/(J-I)
40      CONTINUE
        PD(0,0)=0.0D0
        DO 45 J=1,N
45         PD(0,J)=LS*J*(PM(0,J-1)-X*PM(0,J))/XS
        DO 50 I=1,M
        DO 50 J=I,N
           PD(I,J)=LS*I*X*PM(I,J)/XS+(J+I)
     &             *(J-I+1.0D0)/XQ*PM(I-1,J)
50      CONTINUE
        RETURN
        END
