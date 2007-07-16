        PROGRAM MCOMELP

C       ===================================================
C       Purpose: This program computes complete elliptic
C                integrals K(k) and E(k) using subroutine
C                COMELP
C       Input  : K  --- Modulus k ( 0 ó k ó 1 )
C       Output : CK --- K(k)
C                CE --- E(k)
C       Example:
C                  k         K(k)          E(K)
C                ---------------------------------
C                 .00      1.570796      1.570796
C                 .25      1.596242      1.545957
C                 .50      1.685750      1.467462
C                 .75      1.910990      1.318472
C                1.00       ì            1.000000
C       ===================================================

        DOUBLE PRECISION HK,CK,CE
        WRITE(*,*)'Please enter the modulus k '
        READ(*,*) HK
        WRITE(*,*)'    k         K(k)          E(K)'
        WRITE(*,*)'  ---------------------------------'
        CALL COMELP(HK,CK,CE)
        IF (HK.NE.1.0) WRITE(*,10) HK,CK,CE
        IF (HK.EQ.1.0) WRITE(*,20) HK,CE
10      FORMAT(2X,F5.2,2F14.6)
20      FORMAT(2X,F5.2,7X,'ì',6X,F14.6)
        END


        SUBROUTINE COMELP(HK,CK,CE)

C       ==================================================
C       Purpose: Compute complete elliptic integrals K(k)
C                and E(k)
C       Input  : K  --- Modulus k ( 0 ó k ó 1 )
C       Output : CK --- K(k)
C                CE --- E(k)
C       ==================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PK=1.0D0-HK*HK
        IF (HK.EQ.1.0) THEN
           CK=1.0D+300
           CE=1.0D0
        ELSE
           AK=(((.01451196212D0*PK+.03742563713D0)*PK
     &        +.03590092383D0)*PK+.09666344259D0)*PK+
     &        1.38629436112D0
           BK=(((.00441787012D0*PK+.03328355346D0)*PK+
     &        .06880248576D0)*PK+.12498593597D0)*PK+.5D0
           CK=AK-BK*DLOG(PK)
           AE=(((.01736506451D0*PK+.04757383546D0)*PK+
     &        .0626060122D0)*PK+.44325141463D0)*PK+1.0D0
           BE=(((.00526449639D0*PK+.04069697526D0)*PK+
     &        .09200180037D0)*PK+.2499836831D0)*PK
           CE=AE-BE*DLOG(PK)
        ENDIF
        RETURN
        END
