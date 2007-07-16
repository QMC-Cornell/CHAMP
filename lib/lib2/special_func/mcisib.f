        PROGRAM MCISIB

C       ========================================================
C       Purpose: This program computes the cosine and sine
C                integrals using subroutine CISIB
C       Input :  x  --- Argument of Ci(x) and Si(x)
C       Output:  CI --- Ci(x)
C                SI --- Si(x)
C       Example:

C                   x        Ci(x)           Si(x)
C                ------------------------------------
C                  0.0    - ì                 0
C                  5.0    -.190030D+00      1.549931
C                 10.0    -.454563D-01      1.658348
C                 20.0     .444201D-01      1.548241
C                 30.0    -.330326D-01      1.566757
C                 40.0     .190201D-01      1.586985
C       ========================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*) X
        WRITE(*,*)'   x        Ci(x)           Si(x)'
        WRITE(*,*)'------------------------------------'
        CALL CISIB(X,CI,SI)
        IF (X.NE.0.0D0) WRITE(*,10)X,CI,SI
        IF (X.EQ.0.0D0) WRITE(*,20)
10      FORMAT(1X,F5.1,D16.6,F14.6)
20      FORMAT(3X,' .0',3X,' - ì',17X,'0')
        END


        SUBROUTINE CISIB(X,CI,SI)

C       =============================================
C       Purpose: Compute cosine and sine integrals
C                Si(x) and Ci(x) ( x ò 0 )
C       Input :  x  --- Argument of Ci(x) and Si(x)
C       Output:  CI --- Ci(x)
C                SI --- Si(x)
C       =============================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        X2=X*X
        IF (X.EQ.0.0) THEN
           CI=-1.0D+300
           SI=0.0D0
        ELSE IF (X.LE.1.0D0) THEN
           CI=((((-3.0D-8*X2+3.10D-6)*X2-2.3148D-4)
     &        *X2+1.041667D-2)*X2-0.25)*X2+0.577215665D0+LOG(X)
           SI=((((3.1D-7*X2-2.834D-5)*X2+1.66667D-003)
     &        *X2-5.555556D-002)*X2+1.0)*X
        ELSE
           FX=((((X2+38.027264D0)*X2+265.187033D0)*X2
     &        +335.67732D0)*X2+38.102495D0)/((((X2
     &        +40.021433D0)*X2+322.624911D0)*X2
     &        +570.23628D0)*X2+157.105423D0)
           GX=((((X2+42.242855D0)*X2+302.757865D0)*X2
     &        +352.018498D0)*X2+21.821899D0)/((((X2
     &        +48.196927D0)*X2+482.485984D0)*X2
     &        +1114.978885D0)*X2+449.690326D0)/X
           CI=FX*SIN(X)/X-GX*COS(X)/X
           SI=1.570796327D0-FX*COS(X)/X-GX*SIN(X)/X
        ENDIF
        RETURN
        END
