        PROGRAM MJELP

C       ============================================================
C       Purpose: This program computes Jacobian elliptic functions
C                sn u, cn u and dn u using subroutine JELP
C       Input  : u   --- Argument of Jacobian elliptic fuctions
C                Hk  --- Modulus k ( 0 ó k ó 1 )
C       Output : ESN --- sn u
C                ECN --- cn u
C                EDN --- dn u
C                EPH --- phi ( in degrees )
C       Example:
C                k = .5, ( K(k) = 1.68575035 ), and u = u0*K

C                u0       phi       sn u        cn u        dn u
C              ----------------------------------------------------
C               0.0      .0000    .0000000   1.0000000   1.0000000
C               0.5    47.0586    .7320508    .6812500    .9306049
C               1.0    90.0000   1.0000000    .0000000    .8660254
C               1.5   132.9414    .7320508   -.6812500    .9306049
C               2.0   180.0000    .0000000  -1.0000000   1.0000000
C       ============================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter k and u '
        READ(*,*)HK,U
        WRITE(*,*)
        WRITE(*,*)'   k        u          phi        sn u',
     &            '        cn u        dn u'
        WRITE(*,*)' -------------------------------------',
     &            '---------------------------'
        CALL JELP(U,HK,ESN,ECN,EDN,EPH)
        WRITE(*,10)HK,U,EPH,ESN,ECN,EDN
10      FORMAT(1X,F5.3,F12.7,2X,F9.5,3F12.7)
        END


        SUBROUTINE JELP(U,HK,ESN,ECN,EDN,EPH)

C       ========================================================
C       Purpose: Compute Jacobian elliptic functions sn u, cn u
C                and dn u
C       Input  : u   --- Argument of Jacobian elliptic fuctions
C                Hk  --- Modulus k ( 0 ó k ó 1 )
C       Output : ESN --- sn u
C                ECN --- cn u
C                EDN --- dn u
C                EPH --- phi ( in degrees )
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION R(40)
        PI=3.14159265358979D0
        A0=1.0D0
        B0=DSQRT(1.0D0-HK*HK)
        DO 10 N=1,40
           A=(A0+B0)/2.0D0
           B=DSQRT(A0*B0)
           C=(A0-B0)/2.0D0
           R(N)=C/A
           IF (C.LT.1.0D-7) GO TO 15
           A0=A
10         B0=B
15      DN=2.0D0**N*A*U
        DO 20 J=N,1,-1
           T=R(J)*DSIN(DN)
           SA=DATAN(T/DSQRT(DABS(1.0D0-T*T)))
           D=.5D0*(DN+SA)
20         DN=D
        EPH=D*180.0D0/PI
        ESN=DSIN(D)
        ECN=DCOS(D)
        EDN=DSQRT(1.0D0-HK*HK*ESN*ESN)
        RETURN
        END
