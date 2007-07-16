        PROGRAM MELIT3

C       ==========================================================
C       Purpose: This program computes the elliptic integral of
C                the third kind using subroutine ELIT3
C       Input :  Phi --- Argument ( in degrees )
C                 k  --- Modulus   ( 0 ó k ó 1 )
C                 c  --- Parameter ( 0 ó c ó 1 )
C       Output:  EL3 ÄÄÄ Value of the elliptic integral of the
C                        third kind
C       ==========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter phi, k and c '
        READ(*,*)PHI,HK,C
        CALL ELIT3(PHI,HK,C,EL3)
        WRITE(*,10)EL3
10      FORMAT(1X,'EL3=',F12.8)
        END


        SUBROUTINE ELIT3(PHI,HK,C,EL3)

C       =========================================================
C       Purpose: Compute the elliptic integral of the third kind
C                using Gauss-Legendre quadrature
C       Input :  Phi --- Argument ( in degrees )
C                 k  --- Modulus   ( 0 ó k ó 1.0 )
C                 c  --- Parameter ( 0 ó c ó 1.0 )
C       Output:  EL3 --- Value of the elliptic integral of the
C                        third kind
C       =========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION T(10),W(10)
        LOGICAL LB1,LB2
        DATA T/.9931285991850949,.9639719272779138,
     &         .9122344282513259,.8391169718222188,
     &         .7463319064601508,.6360536807265150,
     &         .5108670019508271,.3737060887154195,
     &         .2277858511416451,.7652652113349734D-1/
        DATA W/.1761400713915212D-1,.4060142980038694D-1,
     &         .6267204833410907D-1,.8327674157670475D-1,
     &         .1019301198172404,.1181945319615184,
     &         .1316886384491766,.1420961093183820,
     &         .1491729864726037,.1527533871307258/
        LB1=HK.EQ.1.0D0.AND.DABS(PHI-90.0).LE.1.0D-8
        LB2=C.EQ.1.0D0.AND.DABS(PHI-90.0).LE.1.0D-8
        IF (LB1.OR.LB2) THEN
            EL3=1.0D+300
            RETURN
        ENDIF
        C1=0.87266462599716D-2*PHI
        C2=C1
        EL3=0.0D0
        DO 10 I=1,10
           C0=C2*T(I)
           T1=C1+C0
           T2=C1-C0
           F1=1.0D0/((1.0D0-C*DSIN(T1)*DSIN(T1))*
     &              DSQRT(1.0D0-HK*HK*DSIN(T1)*DSIN(T1)))
           F2=1.0D0/((1.0D0-C*DSIN(T2)*DSIN(T2))*
     &              DSQRT(1.0D0-HK*HK*DSIN(T2)*DSIN(T2)))
10         EL3=EL3+W(I)*(F1+F2)
        EL3=C1*EL3
        RETURN
        END
