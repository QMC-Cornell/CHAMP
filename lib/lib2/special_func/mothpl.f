        PROGRAM MOTHPL

C       =========================================================
C       Purpose: This program computes orthogonal polynomials:
C                Tn(x) or Un(x) or Ln(x) or Hn(x), and their
C                derivatives using subroutine OTHPL
C       Input :  KF --- Function code
C                       KF=1 for Chebyshev polynomial Tn(x)
C                       KF=2 for Chebyshev polynomial Un(x)
C                       KF=3 for Laguerre polynomial Ln(x)
C                       KF=4 for Hermite polynomial Hn(x)
C                n ---  Order of orthogonal polynomials
C                x ---  Argument
C       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
C                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
C                          n = 0,1,2,...,N ( N ó 100 )
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION PL(0:100),DPL(0:100)
        WRITE(*,*)'KF,N,x = ?'
        READ(*,*)KF,N,X
        WRITE(*,20)KF,N,X
        WRITE(*,*)
        CALL OTHPL(KF,N,X,PL,DPL)
        IF (KF.EQ.1) WRITE(*,*)'  n          Tn(x)            Tn''(x)'
        IF (KF.EQ.2) WRITE(*,*)'  n          Un(x)            Un''(x)'
        IF (KF.EQ.3) WRITE(*,*)'  n          Ln(x)            Ln''(x)'
        IF (KF.EQ.4) WRITE(*,*)'  n          Hn(x)            Hn''(x)'
        WRITE(*,*)'-----------------------------------------'
        DO 10 K=0,N
10         WRITE(*,30)K,PL(K),DPL(K)
20      FORMAT(1X,3HKF=,I3,5X,5HNmax=,I3,5X,2Hx=,F6.3)
30      FORMAT(1X,I3,2D18.8)
        END


        SUBROUTINE OTHPL(KF,N,X,PL,DPL)

C       ==========================================================
C       Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
C                or Ln(x) or Hn(x), and their derivatives
C       Input :  KF --- Function code
C                       KF=1 for Chebyshev polynomial Tn(x)
C                       KF=2 for Chebyshev polynomial Un(x)
C                       KF=3 for Laguerre polynomial Ln(x)
C                       KF=4 for Hermite polynomial Hn(x)
C                n ---  Order of orthogonal polynomials
C                x ---  Argument of orthogonal polynomials
C       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
C                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION PL(0:N),DPL(0:N)
        A=2.0D0
        B=0.0D0
        C=1.0D0
        Y0=1.0D0
        Y1=2.0D0*X
        DY0=0.0D0
        DY1=2.0D0
        PL(0)=1.0D0
        PL(1)=2.0D0*X
        DPL(0)=0.0D0
        DPL(1)=2.0D0
        IF (KF.EQ.1) THEN
           Y1=X
           DY1=1.0D0
           PL(1)=X
           DPL(1)=1.0D0
        ELSE IF (KF.EQ.3) THEN
           Y1=1.0D0-X
           DY1=-1.0D0
           PL(1)=1.0D0-X
           DPL(1)=-1.0D0
        ENDIF
        DO 10 K=2,N
           IF (KF.EQ.3) THEN
              A=-1.0D0/K
              B=2.0D0+A
              C=1.0D0+A
           ELSE IF (KF.EQ.4) THEN
              C=2.0D0*(K-1.0D0)
           ENDIF
           YN=(A*X+B)*Y1-C*Y0
           DYN=A*Y1+(A*X+B)*DY1-C*DY0
           PL(K)=YN
           DPL(K)=DYN
           Y0=Y1
           Y1=YN
           DY0=DY1
10         DY1=DYN
        RETURN
        END
