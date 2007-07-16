        PROGRAM MELIT

C       ==========================================================
C       Purpose: This program computes complete and incomplete
C                elliptic integrals F(k,phi) and E(k,phi) using
C                subroutine ELIT
C       Input  : HK  --- Modulus k ( 0 ó k ó 1 )
C                Phi --- Argument ( in degrees )
C       Output : FE  --- F(k,phi)
C                EE  --- E(k,phi)
C       Example:
C                k = .5

C                 phi     F(k,phi)       E(k,phi)
C                -----------------------------------
C                   0      .00000000      .00000000
C                  15      .26254249      .26106005
C                  30      .52942863      .51788193
C                  45      .80436610      .76719599
C                  60     1.08955067     1.00755556
C                  75     1.38457455     1.23988858
C                  90     1.68575035     1.46746221
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter k and phi (in degs.) '
        READ(*,*)HK,PHI
        WRITE(*,*)
        WRITE(*,*)'  phi     F(k,phi)       E(k,phi)'
        WRITE(*,*)' -----------------------------------'
        CALL ELIT(HK,PHI,FE,EE)
        WRITE(*,10)PHI,FE,EE
10      FORMAT(1X,I4,2F15.8)
        END


        SUBROUTINE ELIT(HK,PHI,FE,EE)

C       ==================================================
C       Purpose: Compute complete and incomplete elliptic
C                integrals F(k,phi) and E(k,phi)
C       Input  : HK  --- Modulus k ( 0 ó k ó 1 )
C                Phi --- Argument ( in degrees )
C       Output : FE  --- F(k,phi)
C                EE  --- E(k,phi)
C       ==================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        G=0.0D0
        PI=3.14159265358979D0
        A0=1.0D0
        B0=DSQRT(1.0D0-HK*HK)
        D0=(PI/180.0D0)*PHI
        R=HK*HK
        IF (HK.EQ.1.0D0.AND.PHI.EQ.90.0D0) THEN
           FE=1.0D+300
           EE=1.0D0
        ELSE IF (HK.EQ.1.0D0) THEN
           FE=DLOG((1.0D0+DSIN(D0))/DCOS(D0))
           EE=DSIN(D0)
        ELSE
           FAC=1.0D0
           DO 10 N=1,40
              A=(A0+B0)/2.0D0
              B=DSQRT(A0*B0)
              C=(A0-B0)/2.0D0
              FAC=2.0D0*FAC
              R=R+FAC*C*C
              IF (PHI.NE.90.0D0) THEN
                 D=D0+DATAN((B0/A0)*DTAN(D0))
                 G=G+C*DSIN(D)
                 D0=D+PI*INT(D/PI+.5D0)
              ENDIF
              A0=A
              B0=B
              IF (C.LT.1.0D-7) GO TO 15
10         CONTINUE
15         CK=PI/(2.0D0*A)
           CE=PI*(2.0D0-R)/(4.0D0*A)
           IF (PHI.EQ.90.0D0) THEN
              FE=CK
              EE=CE
           ELSE
              FE=D/(FAC*A)
              EE=FE*CE/CK+G
           ENDIF
        ENDIF
        RETURN
        END
