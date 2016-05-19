!        SUBROUTINE DGELG

!        PURPOSE
!           TO SOLVE A GENERAL SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS.

!        USAGE
!           CALL DGELG(R,A,M,N,EPS,IER)

!        DESCRIPTION OF PARAMETERS
!           R      - DOUBLE PRECISION M BY N RIGHT HAND SIDE MATRIX
!                    (DESTROYED). ON RETURN R CONTAINS THE SOLUTIONS
!                    OF THE EQUATIONS.
!           A      - DOUBLE PRECISION M BY M COEFFICIENT MATRIX
!                    (DESTROYED).
!           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.
!           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.
!           EPS    - SINGLE PRECISION INPUT CONSTANT WHICH IS USED AS
!                    RELATIVE TOLERANCE FOR TEST ON LOSS OF
!                    SIGNIFICANCE.
!           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS
!                    IER=0  - NO ERROR,
!                    IER=-1 - NO RESULT BECAUSE OF M LESS THAN 1 OR
!                             PIVOT ELEMENT AT ANY ELIMINATION STEP
!                             EQUAL TO 0,
!                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-
!                             CANCE INDICATED AT ELIMINATION STEP K+1,
!                             WHERE PIVOT ELEMENT WAS LESS THAN OR
!                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES
!                             ABSOLUTELY GREATEST ELEMENT OF MATRIX A.

!        REMARKS
!           INPUT MATRICES R AND A ARE ASSUMED TO BE STORED COLUMNWISE
!           IN M*N RESP. M*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN
!           SOLUTION MATRIX R IS STORED COLUMNWISE TOO.
!           THE PROCEDURE GIVES RESULTS IF THE NUMBER OF EQUATIONS M IS
!           GREATER THAN 0 AND PIVOT ELEMENTS AT ALL ELIMINATION STEPS
!           ARE DIFFERENT FROM 0. HOWEVER WARNING IER=K - IF GIVEN -
!           INDICATES POSSIBLE LOSS OF SIGNIFICANCE. IN CASE OF A WELL
!           SCALED MATRIX A AND APPROPRIATE TOLERANCE EPS, IER=K MAY BE
!           INTERPRETED THAT MATRIX A HAS THE RANK K. NO WARNING IS
!           GIVEN IN CASE M=1.

!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!           NONE

!        METHOD
!           SOLUTION IS DONE BY MEANS OF GAUSS-ELIMINATION WITH
!           COMPLETE PIVOTING.

!     ..................................................................

      SUBROUTINE DGELG(R,A,M,N,EPS,IER)


      DIMENSION A(*),R(*)
      REAL*8 R,A,EPS,PIV,TB,TOL,PIVI
      IF(M)23,23,1

!     SEARCH FOR GREATEST ELEMENT IN MATRIX A
    1 IER=0
      PIV=0.D0
      MM=M*M
      NM=N*M
      DO 3 L=1,MM
      TB=DABS(A(L))
      IF(TB-PIV)3,3,2
    2 PIV=TB
      I=L
    3 CONTINUE
      TOL=EPS*PIV
!     A(I) IS PIVOT ELEMENT. PIV CONTAINS THE ABSOLUTE VALUE OF A(I).


!     START ELIMINATION LOOP
      LST=1
      DO 17 K=1,M

!     TEST ON SINGULARITY
      IF(PIV)23,23,4
    4 IF(IER)7,5,7
    5 IF(PIV-TOL)6,6,7
    6 IER=K-1
    7 PIVI=1.D0/A(I)
      J=(I-1)/M
      I=I-J*M-K
      J=J+1-K
!     I+K IS ROW-INDEX, J+K COLUMN-INDEX OF PIVOT ELEMENT

!     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R
      DO 8 L=K,NM,M
      LL=L+I
      TB=PIVI*R(LL)
      R(LL)=R(L)
    8 R(L)=TB

!     IS ELIMINATION TERMINATED
      IF(K-M)9,18,18

!     COLUMN INTERCHANGE IN MATRIX A
    9 LEND=LST+M-K
      IF(J)12,12,10
   10 II=J*M
      DO 11 L=LST,LEND
      TB=A(L)
      LL=L+II
      A(L)=A(LL)
   11 A(LL)=TB

!     ROW INTERCHANGE AND PIVOT ROW REDUCTION IN MATRIX A
   12 DO 13 L=LST,MM,M
      LL=L+I
      TB=PIVI*A(LL)
      A(LL)=A(L)
   13 A(L)=TB

!     SAVE COLUMN INTERCHANGE INFORMATION
      A(LST)=J

!     ELEMENT REDUCTION AND NEXT PIVOT SEARCH
      PIV=0.D0
      LST=LST+1
      J=0
      DO 16 II=LST,LEND
      PIVI=-A(II)
      IST=II+M
      J=J+1
      DO 15 L=IST,MM,M
      LL=L-J
      A(L)=A(L)+PIVI*A(LL)
      TB=DABS(A(L))
      IF(TB-PIV)15,15,14
   14 PIV=TB
      I=L
   15 CONTINUE
      DO 16 L=K,NM,M
      LL=L+J
   16 R(LL)=R(LL)+PIVI*R(L)
   17 LST=LST+M
!     END OF ELIMINATION LOOP


!     BACK SUBSTITUTION AND BACK INTERCHANGE
   18 IF(M-1)23,22,19
   19 IST=MM+M
      LST=M+1
      DO 21 I=2,M
      II=LST-I
      IST=IST-LST
      L=IST-M
      L=int(A(L)+.5D0)
      DO 21 J=II,NM,M
      TB=R(J)
      LL=J
      DO 20 K=IST,MM,M
      LL=LL+1
   20 TB=TB-A(K)*R(LL)
      K=J+L
      R(J)=R(K)
   21 R(K)=TB
   22 RETURN


!     ERROR RETURN
   23 IER=-1
      RETURN
      END
!----------------------------------------------------------CBO----------

!     SUBROUTINE CBO(JP,CR,X1,CI,X2,BR,X3,BI,X4,EPS1,ALFR,ALFI,
!    *BETA,MATV,VR,X5,VI,X6,ITER,IFAIL)

!-----------------------------------------------------------------------
!c    CBO EST UNE INTERFACE ENTRE V ET CBORIS
!-----------------------------------------------------------------------
!c    INTEGER JP,X1,X2,X3,X4,X5,X6,ITER,IFAIL,I,IP,J,K
!     INTEGER JP,X1,IFAIL,I,IP,J,K
!     REAL*8 CR(X1,X1),CI(X1,X1),BR(X1,X1),BI(X1,X1),VR(X1,X1),VI(X1,X1)
!c    REAL*8 EPS1,ALFR(X1),ALFI(X1),BETA(X1)
!     REAL*8 ALFR(X1),ALFI(X1),BETA(X1)
!c    LOGICAL MATV
!     DO 20 I=1,JP
!     DO 20 J=1,JP
!  20 VR(I,J)=CR(I,J)
!     J=JP-1
!     DO 30 I=1,J
!     IP=I+1
!     DO 30 K=IP,JP
!     VR(I,K)=CI(K,I)
!  30 BR(I,K)=BI(K,I)
!     CALL CBORIS(JP,X1,VR,BR,VI,ALFR,ALFI,BETA,IFAIL)
!     DO 100 I=1,X1
! 100 BETA(I)=1.0D0
!     RETURN
!     END
!-----------------------------------------------------------------------

      SUBROUTINE CBORIS(N,ND,A,B,C,D,E,F,FAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      INTEGER*4 N,ND
      DIMENSION A(ND,ND),B(ND,ND),C(ND,ND),D(ND),E(ND),F(ND)
      INTEGER*4 FAIL



! SUBROUTINE CBORIS SOLVES THE COMPLEX EIGENVALUE
! PROBLEM

!      A*X=B*X*LAMBDA

! WHERE

! A     --- GIVEN MATRIX OF ORDER N
! B     --- GIVEN POSITIVE DEFINITE MATRIX OF ORDER N
! X     --- EIGENVECTORS, COLUMN BY COLUMN,
!           NORMALIZED TO (BX,X)=I
! LAMBDA --- EIGENVALUES IN INCREASING ORDER

! MATRICES A AND B SHOULD BE GIVEN IN
! ARRAYS A AND B RESPECTIVELY IN THE FOLLOWING WAY

! DIAGONAL ELEMENTS IN THE DIAGONAL
! REAL PARTS OF THE LOWER TRIANGLE IN THE LOWER
! TRIANGLE
! IMAG PARTS OF THE LOWER TRIANGLE IN THE UPPER
! TRIANGLE

! REAL PARTS OF EIGENVECTORS ARE STORED IN A,
! IMAG PARTS IN ARRAY C.  ARRAY B IS DESTROYED DURING
! COMPUTATION.  IT HOLDS THE LOWER TRIANGLE OF THE
! CHOLESKI DECOMPOSITION OF MATRIX B (SEE
! SUBROUTINE CCHOL).

! ALL MATRICES ARE OF ORDER N, DECLARED IN THE
! CALLING PROGRAM WITH DIMENSION ND WHICH NEED
! NOT TO BE EQUAL TO N.

! EIGENVALUES ARE STORED IN ARRAY D.  ARRAYS E
! AND F ARE USED FOR INTERMEDIATE RESULTS.

! FAIL GETS THE FOLLOWING VALUES

!  0 --- COMPUTATION FINISHED SUCCESFULLY
!  1 --- B NOT POSITIVE DEFINITE
!  2 --- QR ALGORITHM DOES NOT CONVERGE

! SUBROUTINES USED

!     CCHOL
!     CTRED2
!     CTQL2


! PROGRAMMED BY E. ZAKRAJSEK
! JUNE 21,1974


! DECOMPOSE MATRIX B

      FAIL=1
      CALL CCHOL(N,ND,B,LF)
      IF(LF.NE.0) RETURN

! MOVE MATRIX A

      DO 11 I=1,N
      C(I,I)=0.D0
      IF(I.EQ.1) GO TO 11
      IA=I-1
      DO 10 J=1,IA
      C(I,J)=A(J,I)
      C(J,I)=-A(J,I)
   10 A(J,I)=A(I,J)
   11 CONTINUE

!  COMPUTE (L(-1)*A)

      DO 22 J=1,N
      DO 22 I=1,N
      IF(I.EQ.1) GO TO 21
      IA=I-1
      DO 20 K=1,IA
      A(I,J)=A(I,J)-A(K,J)*B(I,K)+C(K,J)*B(K,I)
   20 C(I,J)=C(I,J)-A(K,J)*B(K,I)-C(K,J)*B(I,K)
   21 A(I,J)=A(I,J)/B(I,I)
   22 C(I,J)=C(I,J)/B(I,I)

!  COMPUTE  A*L(-H)

      DO 32 I=1,N
      DO 32 J=1,I
      IF(J.EQ.1) GO TO 31
      JA=J-1
      DO 30 K=1,JA
      A(I,J)=A(I,J)-A(I,K)*B(J,K)-C(I,K)*B(K,J)
   30 C(I,J)=C(I,J)+A(I,K)*B(K,J)-C(I,K)*B(J,K)
   31 A(I,J)=A(I,J)/B(J,J)
   32 C(I,J)=C(I,J)/B(J,J)

!     PUT MATRIX TOGETHER INTO A

      DO 41 I=1,N
      IF(I.EQ.N) GO TO 41
      IA=I+1
      DO 40 J=IA,N
   40 A(I,J)=C(J,I)
   41 CONTINUE

!     DIAGONALIZE A

      FAIL=2
      CALL CTRED2(N,ND,A,C,D,E,F)
      CALL CTQL2(N,ND,D,E,F,A,C,LF)
      IF(LF.NE.0) RETURN

!     COMPUTE L(-H)*A

      DO 52 J=1,N
      DO 52 II=1,N
      I=N-II+1
      IF(I.EQ.N) GO TO 51
      IA=I+1
      DO 50 K=IA,N
      A(I,J)=A(I,J)-A(K,J)*B(K,I)-C(K,J)*B(I,K)
   50 C(I,J)=C(I,J)+A(K,J)*B(I,K)-C(K,J)*B(K,I)
   51 A(I,J)=A(I,J)/B(I,I)
      C(I,J)=C(I,J)/B(I,I)
   52 CONTINUE
      FAIL=0
      RETURN
      END
!-----------------------------------------------------------------------

      SUBROUTINE CCHOL(N,ND,A,FAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      INTEGER*4 N,ND
      DIMENSION A(ND,ND)
      INTEGER*4 FAIL

! SUBROUTINE CCHOL COMPUTES CHOLESKI
! DECOMPOSITION OF GIVEN COMPLEX POSITIVE DEFINITE
! MATRIX A.

! INPUT DATA

!  N --- ORDER OF MATRIX
!  ND -- DIMENSION OF ARRAY A (IT CAN BE
!        GREATER THAN OR EQUAL TO N)
!  A --- GIVEN MATRIX
!        IT IS SUPPOSED TO BE STORED IN THE
!        FOLLOWING WAY.  DIAGONAL ELEMENTS,
!        BEING REAL, ARE STORED ON THE DIAGONAL,
!        REAL PARTS OF OFF-DIAGONAL ELEMENTS
!        ARE STORED IN THE LOWER TRIANGLE OF A,
!        IMAG PARTS OF THE LOWER TRIANGLE ARE
!        STORED IN THE UPPER TRIANGLE OF A.

!     EXIT INFORMATION

!  A --- COMPLEX ELEMENTS OF MATRIX L, DEFINED BY
!     A=L*L(H)
!        ARE STORED IN THE SAME WAY AS ORIGINAL
!        ELEMENTS OF A, THAT IS REAL PARTS OF THE
!        LOWER TRIANGLE OF L IN THE LOWER TRIANGLE
!        OF A AND THE CORRESPONDING IMAG PARTS IN
!        THE UPPER TRIANGLE OF A.
!  FAIL --- IS SET TO ZERO IF THE DECOMPOSITION WAS
!           SUCCESSFUL AND TO NONZERO IF
!           THE MATRIX WAS NOT POSITIVE DEFINITE.


!    PROGRAMMED BY E. ZAKRAJSEK
!     JUNE 20, 1974



!     SUPPOSE DECOMPOSITION WILL FAIL

      FAIL=1

      DO 13 I=1,N

!     TEST FOR POSITIVE DEFINITENESS

      IF(A(I,I).LE.0.0D0) RETURN

!      COMPUTE COLUMN I

      A(I,I)=DSQRT(A(I,I))
      IF(I.EQ.N) GO TO 13
      IA=I+1
      DO 10 J=IA,N
      A(J,I)=A(J,I)/A(I,I)
   10 A(I,J)=A(I,J)/A(I,I)

!     REDUCE REMAINING COLUMNS

      DO 12 K=IA,N
      A(K,K)=A(K,K)-A(K,I)**2-A(I,K)**2
      IF(K.EQ.N) GO TO 12
      KA=K+1
      DO 11 J=KA,N
      A(J,K)=A(J,K)-A(J,I)*A(K,I)-A(I,J)*A(I,K)
   11 A(K,J)=A(K,J)-A(I,J)*A(K,I)+A(J,I)*A(I,K)
   12 CONTINUE
   13 CONTINUE
      FAIL=0
      RETURN
      END
!-----------------------------------------------------------------------

      SUBROUTINE CTRED2(N,ND,A,B,D,E,F)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      INTEGER*4 N,ND
      DIMENSION A(ND,ND),B(ND,ND),D(ND),E(ND),F(ND)


!     SUBROUTINE CTRED2 REDUCES GIVEN COMPLEX
!     HERMITIAN MATRIX TO A TRIDIAGONAL FORM

!     PARAMETERS:

!     N    --- ORDER OF THE MATRIX
!     ND   --- DIMENSION OF ARRAYS A AND B
!     A    --- GIVEN MATRIX, REPLACED BY REAL PART
!              OF THE TRANSFORMATION MATRIX
!     B    --- IMAG PART OF TRANSFORMATION MATRIX
!     D    --- DIAGONAL PART OF THE TRIADIAGONAL MATRIX
!     E    --- REAL PART OF THE CODIAGONAL OF THE
!              TRIDIAGONAL MATRIX
!              (LAST N-1 LOCATIONS)
!     F    --- IMAG PARTS OF THE LOWER CODIAGONAL.

!     THE GIVEN MATRIX SHOULD BE STORED IN THE
!     FOLLOWING WAY:

!          --- DIAGONAL ELEMENTS IN THE DIAGONAL
!          --- REAL PART OF THE LOWER TRIANGLE IN THE
!              LOWER TRIANGLE
!          --- IMAG PARTS OF THE LOWER TRIANGLE
!              IN THE UPPER TRIANGLE


!     PROGRAMMED BY E. ZAKRAJSEK
!     JUNE 20,1974



!     CHEP=2.0D0**(-56)
      D(1)=A(1,1)
      IF(N.EQ.1) GO TO 31

! MAIN K LOOP

      DO 30 K=2,N
      L=K-1

!     COMPUTE NORM

      ALL=0.D0
      DO 10 I=K,N
   10 ALL=ALL+A(I,L)**2 + A(L,I)**2
      ALL=DSQRT(ALL)

!     COMPUTE CONSTANTS

      C=1.0D0
      S=0.D0
      R = DSQRT(A(K,L)**2 + A(L,K)**2)
      IF(DABS(R).LT.1.0D-19) R=0.D0
      IF(R.EQ.0.0D0) GOTO 11
      C=A(K,L)/R
      S=A(L,K)/R
   11 ALR=ALL*C
      ALI=ALL*S
      A(L,L)=0.0D0

!     TEST FOR SUPERFLUOUS TRANSFORMATION

      SM=ALL*(ALL+R)
      IF(DABS(SM).LT.1.0D-19) SM=0.D0
      IF(SM.EQ.0.D0) GO TO 20
      G=1.0D0/SM
      A(L,L)=G

      A(K,L)=A(K,L)+ALR
      A(L,K) =A(L,K)+ALI

!     NOW COMPUTE U=A*W
!     AND STORE INTO (E,F)

      T=0.0D0
      DO 16 I=K,N
      C=A(I,I)*A(I,L)
      S=A(I,I)*A(L,I)
      IF(I.EQ.K)GOTO 13
      IA=I-1
      DO 12 J=K,IA
      C=C+A(I,J)*A(J,L)-A(J,I)*A(L,J)
   12 S=S+A(I,J)*A(L,J)+A(J,I)*A(J,L)
   13 IF(I.EQ.N) GOTO 15
      IA=I+1
      DO 14 J=IA,N
      C=C+A(J,I)*A(J,L)+A(I,J)*A(L,J)
   14 S=S+A(J,I)*A(L,J)-A(I,J)*A(J,L)
   15 E(I)=G*C
      F(I)=G*S
   16 T=T+A(I,L)*C+A(L,I) *S
      T=T*G**2

!    TRANSFORM  MATRIX

      DO 18 I=K,N
      A(I,I)=A(I,I)-2.0D0*(A(I,L)*E(I)+A(L,I)*F(I))+
     *      T*(A(I,L)**2+A(L,I)**2)
      IF(I.EQ.K) GOTO 18
      IA=I-1
      DO 17 J=K,IA
      A(I,J)=A(I,J)-A(I,L)*E(J)-A(L,I)*F(J)
     *      -A(J,L)*E(I)-A(L,J)*F(I)
     *      +T*(A(I,L)*A(J,L)+A(L,I)*A(L,J))
      A(J,I)=A(J,I)-A(L,I)*E(J)+A(I,L)*F(J)
     *      +A(L,J)*E(I)-A(J,L)*F(I)
     *      +T*(A(L,I)*A(J,L)-A(I,L)*A(L,J))
   17 CONTINUE
   18 CONTINUE

!     STORE DIAGONAL AND CODIAGONAL ELEMENTS

   20 D(K)=A(K,K)
      E(K)=-ALR
      F(K)=-ALI
   30 CONTINUE

!     NOW ACCUMULATE TRANSFORMATIONS

   31 A(N,N)=1.D0
      B(N,N)=0.D0
      IF(N.EQ.1) RETURN
      DO 40 KK=2,N
      K=N-KK+2
      L=K-1

!     SKIP TRANSFORMATION IF UNIT

      IF(DABS(A(L,L)).LT.1.0D-19) A(L,L)=0.D0
      IF(A(L,L).EQ.0.0D0) GO TO 36

!     COMPUTE PRODUCT

      DO 35 J=K,N
      C=0.D0
      S=0.D0
      DO 33 I=K,N
      C=C+A(I,L)*A(I,J)+A(L,I)*B(I,J)
   33 S=S+A(I,L)*B(I,J)-A(L,I)*A(I,J)
      C=C*A(L,L)
      S=S*A(L,L)
      DO 34 I=K,N
      A(I,J)=A(I,J)-C*A(I,L)+S*A(L,I)
   34 B(I,J)=B(I,J)-C*A(L,I)-S*A(I,L)
   35 CONTINUE

!     MAKE NEW LINE

   36 DO 37 I=K,N
      A(I,L)=0.D0
      A(L,I)=0.D0
      B(I,L)=0.D0
   37 B(L,I)=0.D0
      A(L,L)=1.D0
      B(L,L)=0.D0
   40 CONTINUE
      RETURN
      END
!-----------------------------------------------------------------------

      SUBROUTINE CTQL2(N,ND,D,E,F,A,B,FAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      INTEGER*4 N,ND
      DIMENSION A(ND,ND), B(ND,ND), D(ND), E(ND), F(ND)
      INTEGER*4 FAIL


!     SUBROUTINE CTQL2 COMPUTES THE EIGENVALUES AND
!     EIGENVECTORS OF A COMPLEX HERMITIAN TRIDIAGONAL
!     MATRIX

!     PARAMETERS:

!     N    --- ORDER OF MATRIX
!     ND   --- DIMENSION OF A AND B
!     D    --- DIAGONAL GIVEN
!     E    --- REAL PART OF CODIAGONAL GIVEN
!              (LAST N-1 LOCATIONS)
!     F    --- IMAG PART OF THE LOWER CODIAGONAL
!     A    --- REAL PART OF EIGENVECTORS
!     B    --- IMAG PART OF EIGENVECTORS
!     FAIL --- RECEIVES VALUE OF 1 INSTEAD OF ZERO
!              IF SOME EIGENVALUE TAKES MORE THAN 30
!              ITERATIONS.


!     EIGENVALUES ARE OBTAINED IN INCREASING OF
!     MAGNITUDE IN VECTOR D, EIGENVECTORS ARE STORED
!     BY COLUMNS.  ARRAYS A AND B SHOULD BE PRESET TO
!     SOME UNITARY MATRIX SUCH AS THE IDENTITY MATRIX
!     OR THE MATRIX PRODUCED BY CTRED2.


!     PROGRAMMED BY E.  ZAKRAJSEK
!     JUNE 21, 1974




!     ***************************************
!     *                                     *
!     * NEXT LINE OF PROGRAM DEFINES        *
!     * MACHINE DEPENDENT CONSTANT CHEP     *
!     * DEFINED AS THE SMALLEST REAL        *
!     * NUMBER FOR WHICH                    *
!     *                                     *
!     *        1.0+CHEP .GT. 1.0            *
!     *                                     *
!     ***************************************

      CHEP=2.0D0**(-56)

!     FIRST MAKE REAL CODIAGONAL MOVED DOWN
!     TO FIRST LOCATION

      IF(N.EQ.1) GO TO 12
      DO 11 K=2,N
      R=DSQRT(E(K)**2+F(K)**2)
      IF(DABS(R).LT.1.0D-19) R=0.D0
      IF(R.EQ.0.0D0) GO TO 11
      C=E(K)/R
      S=F(K)/R

!     ACCUMULATE ROTATION

      DO 10 I=1,N
      P=A(I,K)*C-B(I,K)*S
      B(I,K)=A(I,K)*S+B(I,K)*C
   10 A(I,K)=P

!     TRANSFORM NEXT E

      IF(K.EQ.N) GO TO 11
      L=K+1
      P=E(L)*C-F(L)*S
      F(L)=E(L)*S+F(L)*C
      E(L)=P
   11 E(K-1)=R
   12 E(N)=0.D0

!     INITIALIZE

      BB=0.D0
      FF=0.D0
      FAIL=1

!     MAIN LOOP

      DO 32 L=1,N
      J=0
      H=CHEP*(DABS(D(L))+DABS(E(L)))
      IF(BB.LT.H) BB=H

!     LOOK FOR SMALL SUBDIAGONAL ELEMENT

      DO 20 M=L,N
      IF(DABS(E(M)).LE.BB) GO TO 21
   20 CONTINUE
   21 IF(M.EQ.L) GO TO 31

!     NEXT ITERATION

   24 IF(J.EQ.30) RETURN
      J=J+1

!     FORM SHIFT

      P=(D(L+1)-D(L))/(2.0D0*E(L))
      R=DSQRT(1.0D0+P**2)
      H=D(L)-E(L)/(P+DSIGN(R,P))

      DO 25 I=L,N
   25 D(I)=D(I)-H
      FF=FF+H

!     QL TRANSFORMATION

      P=D(M)
      C=1.0D0
      S=0.0D0
      MA=M-1

      DO 30 IA=L,MA
      I=MA-IA+L
      I1=I+1
      G=C*E(I)
      H=C*P
      IF(DABS(P).LT.DABS(E(I))) GO TO 26

      C=E(I)/P
      R=DSQRT(C**2+1.0D0)
      E(I1 )=S*P*R
      S=C/R
      C=1.0D0/R
      GO TO 27

   26 C=P/E(I)
      R=DSQRT(C**2+1.D0)
      E(I1 )=S*E(I)*R
      S=1.0D0/R
      C=C/R

   27 P=C*D(I)-S*G
      D(I1 )=H+S*(C*G+S*D(I))

!     FORM VECTOR

      DO 28 K=1,N
      HR=A(K,I1 )
      HI=B(K,I1 )
      A(K,I1 )=S*A(K,I)+C*HR
      B(K,I1 )=S*B(K,I)+C*HI
      A(K,I)=C*A(K,I)-S*HR
   28 B(K,I)=C*B(K,I)-S*HI

   30 CONTINUE

      E(L)=S*P
      D(L)=C*P
      IF(DABS(E(L)).GT.BB) GO TO 24

!     ROOT FOUND

   31 D(L)=D(L)+FF
   32 CONTINUE

!     ORDER EIGENVALUES AND EIGENVECTORS

      DO 42 I=1,N
      K=I
      DO 40 J=I,N
      IF(D(J).LT.D(K)) K=J
   40 CONTINUE


      IF(K.EQ.I) GO TO 42
      P=D(I)
      D(I)=D(K)
      D(K)=P

      DO 41 J=  1,N
      P=A(J,I)
      A(J,I)=A(J,K)
      A(J,K)=P

      P=B(J,I)
      B(J,I)=B(J,K)
   41 B(J,K)=P
   42 CONTINUE
      FAIL=0
      RETURN
      END
