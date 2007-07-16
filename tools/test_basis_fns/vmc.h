! DBLMIN, DBLMAX should be the smallest and largest double precision numbers
! that can be represented, but approximate values are OK.
!JT      parameter(DBLMIN=1.d-300,DBLMAX=1.d+300)
      parameter(DBLMIN=1.d-30,DBLMAX=1.d+30)

! MORB     must be at least as large as the largest iworb in the input
! MORB_OCC must be at least as large as the number of occupied orbitals
! Since MORB_OCC <= MORB, so far I have just made sure that in the large arrays
! I use MORB_OCC if that is all that is needed, but there are probably some smaller
! arrays that could have their dimensions changed from MORB to MORB_OCC.

!     parameter(MELEC=128,MORB=80,MBASIS=36,MBASIS_CTYPE=36,MDET=1,MCENT=16,MCTYPE=1,
!     parameter(MELEC=162,MORB=82,MBASIS=900,MBASIS_CTYPE=36,MDET=1,MCENT=63,MCTYPE=3,
!     parameter(MELEC=64,MORB=40,MBASIS=246,MBASIS_CTYPE=36,MDET=45,MCENT=23,MCTYPE=2,
!     parameter(MELEC=128,MORB=80,MBASIS=246,MBASIS_CTYPE=36,MDET=9,MCENT=23,MCTYPE=2,
!     parameter(MELEC=128,MORB=80,MBASIS=246,MBASIS_CTYPE=36,MDET=12,MCENT=23,MCTYPE=2,
!     parameter(MELEC=64,MORB=40,MBASIS=100,MBASIS_CTYPE=36,MDET=46,MCENT=16,MCTYPE=4,
!     parameter(MELEC=80,MORB=48,MBASIS=100,MBASIS_CTYPE=36,MDET=46,MCENT=16,MCTYPE=4,
!     parameter(MELEC=64,MORB=40,MBASIS=246,MBASIS_CTYPE=36,MDET=46,MCENT=23,MCTYPE=4,
!     parameter(MELEC=216,MORB=134,MBASIS=246,MBASIS_CTYPE=36,MDET=1,MCENT=23,MCTYPE=4,
!     parameter(MELEC=256,MORB=128,MBASIS=20,MBASIS_CTYPE=20,MDET=18,MCENT=64,MCTYPE=2,
!     parameter(MELEC=256,MORB=136,MBASIS=4,MBASIS_CTYPE=4,MDET=2,MCENT=64,MCTYPE=2,
!     parameter(MELEC=512,MORB=256,MBASIS=2,MBASIS_CTYPE=2,MDET=1,MCENT=128,MCTYPE=2,
!     parameter(MELEC=256,MORB=136,MORB_OCC=128,MBASIS=4,MBASIS_CTYPE=4,MDET=2,MCENT=64,MCTYPE=2,
!     parameter(MELEC=30,MORB=14,MORB_OCC=14,MBASIS=14,MBASIS_CTYPE=14,MDET=8,MCENT=2,MCTYPE=2,
!     parameter(MELEC=216,MORB=140,MORB_OCC=140,MBASIS=92,MBASIS_CTYPE=40,MDET=36,MCSF=20,MDET_CSF=8,MCENT=3,MCTYPE=2,
!JT      parameter(MELEC=30,MORB=18,MORB_OCC=18,MBASIS=92,MBASIS_CTYPE=48,MDET=68,MCSF=20,MDET_CSF=20,MCENT=3,MCTYPE=2,
      parameter(MELEC=30,MORB=230,MORB_OCC=18,MBASIS=228,MBASIS_CTYPE=150)
      parameter(MDET=660,MCSF=200,MDET_CSF=50,MCENT=12,MCTYPE=2)
      parameter(NSPLIN=1001,MORDJ=7,MPARMJ=70,ML_BAS=4,radmax=8.d0,nrad=801)
      parameter(MMAT_DIM=(MELEC*MELEC)/4,MMAT_DIM2=(MELEC*(MELEC-1))/2)
      parameter(MORDJ1=MORDJ+1,delri=(nrad-1)/radmax)
