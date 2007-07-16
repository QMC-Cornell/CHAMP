        PROGRAM MCLPMN

C       ============================================================
C       Purpose: This program computes the associated Legendre
C                functions Pmn(z) and their derivatives Pmn'(z) for
C                a complex argument using subroutine CLPMN
C       Input :  x --- Real part of z
C                y --- Imaginary part of z
C                m --- Order of Pmn(z),  m = 0,1,2,...,n
C                n --- Degree of Pmn(z), n = 0,1,2,...,N
C       Output:  CPM(m,n) --- Pmn(z)
C                CPD(m,n) --- Pmn'(z)
C       Examples:
C                n = 5, x = 0.5, y = 0.2

C       m     Re[Pmn(z)]    Im[Pmn(z)]    Re[Pmn'(z)]   Im[Pmn'(z)]
C      -------------------------------------------------------------
C       0    .252594D+00  -.530293D+00   -.347606D+01  -.194250D+01
C       1    .333071D+01   .135206D+01    .117643D+02  -.144329D+02
C       2   -.102769D+02   .125759D+02    .765713D+02   .598500D+02
C       3   -.572879D+02  -.522744D+02   -.343414D+03   .147389D+03
C       4    .335711D+03  -.389151D+02   -.226328D+03  -.737100D+03
C       5   -.461125D+03   .329122D+03    .187180D+04   .160494D+02

C                n = 5, x = 2.5, y = 1.0

C       m     Re[Pmn(z)]    Im[Pmn(z)]    Re[Pmn'(z)]   Im[Pmn'(z)]
C      -------------------------------------------------------------
C       0   -.429395D+03   .900336D+03   -.350391D+02   .193594D+04
C       1   -.216303D+04   .446358D+04   -.208935D+03   .964685D+04
C       2   -.883477D+04   .174005D+05   -.123703D+04   .381938D+05
C       3   -.273211D+05   .499684D+05   -.568080D+04   .112614D+06
C       4   -.565523D+05   .938503D+05   -.167147D+05   .219713D+06
C       5   -.584268D+05   .863328D+05   -.233002D+05   .212595D+06
C       ============================================================

        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CPM(0:40,0:40),CPD(0:40,0:40)
        WRITE(*,*)'  Please enter m, n, x and y '
        READ(*,*) M,N,X,Y
        WRITE(*,30) M,N,X,Y
        CALL CLPMN(40,M,N,X,Y,CPM,CPD)
        WRITE(*,*)'   m   n    Re[Pmn(z)]    Im[Pmn(z)]    ',
     &            'Re[Pmn''(z)]   Im[Pmn''(z)]'
        WRITE(*,*)' -----------------------------------',
     &            '-------------------------------'
        DO 10 J=0,N
10         WRITE(*,20)M,J,CPM(M,J),CPD(M,J)
20      FORMAT(1X,2I4,1X,2D14.6,1X,2D14.6)
30      FORMAT(1X,'m =',I2,', ','n =',I2,', ','x =',F5.1,
     &         ', ','y =',F5.1)
        END


        SUBROUTINE CLPMN(MM,M,N,X,Y,CPM,CPD)

C       =========================================================
C       Purpose: Compute the associated Legendre functions Pmn(z)
C                and their derivatives Pmn'(z) for a complex
C                argument
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                m  --- Order of Pmn(z),  m = 0,1,2,...,n
C                n  --- Degree of Pmn(z), n = 0,1,2,...,N
C                mm --- Physical dimension of CPM and CPD
C       Output:  CPM(m,n) --- Pmn(z)
C                CPD(m,n) --- Pmn'(z)
C       =========================================================

        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CPM(0:MM,0:N),CPD(0:MM,0:N)
        Z=CMPLX(X,Y)
        DO 10 I=0,N
        DO 10 J=0,M
           CPM(J,I)=(0.0D0,0.0D0)
10         CPD(J,I)=(0.0D0,0.0D0)
        CPM(0,0)=(1.0D0,0.0D0)
        IF (DABS(X).EQ.1.0D0.AND.Y.EQ.0.0D0) THEN
           DO 15 I=1,N
              CPM(0,I)=X**I
15            CPD(0,I)=0.5D0*I*(I+1)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 CPD(I,J)=(1.0D+300,0.0D0)
              ELSE IF (I.EQ.2) THEN
                 CPD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (CDABS(Z).GT.1.0D0) LS=-1
        ZQ=CDSQRT(LS*(1.0D0-Z*Z))
        ZS=LS*(1.0D0-Z*Z)
        DO 25 I=1,M
25         CPM(I,I)=-LS*(2.0D0*I-1.0D0)*ZQ*CPM(I-1,I-1)
        DO 30 I=0,M
30         CPM(I,I+1)=(2.0D0*I+1.0D0)*Z*CPM(I,I)
        DO 35 I=0,M
        DO 35 J=I+2,N
           CPM(I,J)=((2.0D0*J-1.0D0)*Z*CPM(I,J-1)-(I+J-
     &              1.0D0)*CPM(I,J-2))/(J-I)
35      CONTINUE
        CPD(0,0)=(0.0D0,0.0D0)
        DO 40 J=1,N
40         CPD(0,J)=LS*J*(CPM(0,J-1)-Z*CPM(0,J))/ZS
        DO 45 I=1,M
        DO 45 J=I,N
           CPD(I,J)=LS*I*Z*CPM(I,J)/ZS+(J+I)*(J-I+1.0D0)
     &              /ZQ*CPM(I-1,J)
45      CONTINUE
        RETURN
        END
