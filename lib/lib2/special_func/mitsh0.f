        PROGRAM MITSH0

C       ====================================================
C       Purpose: This program evaluates the integral of
C                Struve function H0(t) with respect to t
C                from 0 and x using subroutine ITSH0
C       Input :  x   --- Upper limit  ( x ò 0 )
C       Output:  TH0 --- Integration of H0(t) from 0 and x
C       Example:
C                    x        H0(t)dt
C                 ----------------------
C                   0.0       .0000000
C                   5.0      2.0442437
C                  10.0      2.5189577
C                  15.0      2.5415824
C                  20.0      2.5484517
C                  30.0      3.0625848
C                  40.0      3.1484123
C                  50.0      3.2445168
C       ====================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        WRITE(*,*)'Please enter x '
        READ(*,*)X
        WRITE(*,*)'   x         H0(t)dt'
        WRITE(*,*)'----------------------Ä'
        CALL ITSH0(X,TH0)
        WRITE(*,10)X,TH0
10      FORMAT(1X,F5.1,E16.7)
        END


        SUBROUTINE ITSH0(X,TH0)

C       ===================================================
C       Purpose: Evaluate the integral of Struve function
C                H0(t) with respect to t from 0 and x
C       Input :  x   --- Upper limit  ( x ò 0 )
C       Output:  TH0 --- Integration of H0(t) from 0 and x
C       ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(25)
        PI=3.141592653589793D0
        R=1.0D0
        IF (X.LE.30.0) THEN
           S=0.5D0
           DO 10 K=1,100
              RD=1.0D0
              IF (K.EQ.1) RD=0.5D0
              R=-R*RD*K/(K+1.0D0)*(X/(2.0D0*K+1.0D0))**2
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         TH0=2.0D0/PI*X*X*S
        ELSE
           S=1.0D0
           DO 20 K=1,12
              R=-R*K/(K+1.0D0)*((2.0D0*K+1.0D0)/X)**2
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 25
20         CONTINUE
25         EL=.57721566490153D0
           S0=S/(PI*X*X)+2.0D0/PI*(DLOG(2.0D0*X)+EL)
           A0=1.0D0
           A1=5.0D0/8.0D0
           A(1)=A1
           DO 30 K=1,20
              AF=((1.5D0*(K+.5D0)*(K+5.0D0/6.0D0)*A1-.5D0
     &           *(K+.5D0)*(K+.5D0)*(K-.5D0)*A0))/(K+1.0D0)
              A(K+1)=AF
              A0=A1
30            A1=AF
           BF=1.0D0
           R=1.0D0
           DO 35 K=1,10
              R=-R/(X*X)
35            BF=BF+A(2*K)*R
           BG=A(1)/X
           R=1.0D0/X
           DO 40 K=1,10
              R=-R/(X*X)
40            BG=BG+A(2*K+1)*R
           XP=X+.25D0*PI
           TY=DSQRT(2.0D0/(PI*X))*(BG*DCOS(XP)-BF*DSIN(XP))
           TH0=TY+S0
        ENDIF
        RETURN
        END
