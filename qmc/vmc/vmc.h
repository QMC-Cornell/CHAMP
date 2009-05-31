! Used to set dimensions of arrays.
! MORB           orbitals
! MORB_OCC       occupied orbitals

! MORB     must be at least as large as the largest iworb in the input
! MORB_OCC must be at least as large as the number of occupied orbitals
! Since MORB_OCC <= MORB, so far I have just made sure that in the large arrays
! I use MORB_OCC if that is all that is needed, but there are probably some smaller
! arrays that could have their dimensions changed from MORB to MORB_OCC.

! DBLMIN, DBLMAX should be the smallest and largest double precision numbers
! that can be represented, but approximate values are OK.
!JT   parameter(DBLMIN=1.d-300,DBLMAX=1.d+300)
      parameter(DBLMIN=1.d-30,DBLMAX=1.d+30)

! Cyrus:
!      parameter(MORB_OCC=24, MORB=MAX(64,MORB_OCC))
!      parameter(ML_BAS=4)
       parameter(MORB_OCC=42, MORB=MAX(80,MORB_OCC))
!      parameter(ML_BAS=4)

! Julien:
!      parameter(MORB_OCC=10, MORB=MAX(68,MORB_OCC))
       parameter(ML_BAS=4)

! Ryo, Masayoshi:
!     parameter(MORB_OCC=90, MORB=MAX(90,MORB_OCC))
!     parameter(ML_BAS=6)

! Wissam:
!     parameter(MORB_OCC=34, MORB=MAX(50,MORB_OCC))
!     parameter(ML_BAS=4)

      parameter(MORDJ=9, MPARMJ=70, NSPLIN=1001, NRAD=1001, RADMAX=8.d0)
      parameter(NAX=50,NAK1=40,NAK2=1)
      parameter(MORDJ1=MORDJ+1,delri=(NRAD-1)/RADMAX)
