! Used to set dimensions of arrays.
! MCENT          centers (nuclei)
! MCTYPE         center types
! MELEC          electrons
! MELECUD        up or dn electrons
! MORB           orbitals
! MORB_OCC       occupied orbitals
! MBASIS         basis functions for finite systems
! MBASIS_CTYPE   basis function on a single center type
! MDET           determinants
! MDETUD         up or dn determinants
! MCSF           configuration state functions
! MDET_CSF       determinants in a CSF

! MORB     must be at least as large as the largest iworb in the input
! MORB_OCC must be at least as large as the number of occupied orbitals
! Since MORB_OCC <= MORB, so far I have just made sure that in the large arrays
! I use MORB_OCC if that is all that is needed, but there are probably some smaller
! arrays that could have their dimensions changed from MORB to MORB_OCC.

! DBLMIN, DBLMAX should be the smallest and largest double precision numbers
! that can be represented, but approximate values are OK.
!JT   parameter(DBLMIN=1.d-300,DBLMAX=1.d+300)
      parameter(DBLMIN=1.d-30,DBLMAX=1.d+30)

! Standard:
!     parameter(MCENT=66,MCTYPE=2,MELEC=66,MELECUD=33,MORB=66,MORB_OCC=33)
!     parameter(MBASIS=66,MBASIS_CTYPE=40)
!     parameter(MDET=100,MDETUD=10,MCSF=30,MDET_CSF=20)
!     parameter(MOTYPE=5,MPARMD=170)
!     parameter(NSPLIN=1001,MORDJ=7,MPARMJ=70,ML_BAS=4,NRAD=1001,RADMAX=8.d0)

! Julien:
      parameter(MCENT=66,MCTYPE=2,MELEC=66,MELECUD=33,MORB=66,MORB_OCC=33)
      parameter(MBASIS=66,MBASIS_CTYPE=40)
      parameter(MDET=100,MDETUD=10,MCSF=30,MDET_CSF=20)
      parameter(MOTYPE=5,MPARMD=170)
      parameter(NSPLIN=1001,MORDJ=7,MPARMJ=70,ML_BAS=4,NRAD=1001,RADMAX=8.d0)

! Cyrus:
!      parameter(MCENT=8,MCTYPE=2,MELEC=32,MELECUD=20,MORB=102,MORB_OCC=18)
!      parameter(MBASIS=272,MBASIS_CTYPE=102)
!      parameter(MDET=660,MDETUD=190,MCSF=165,MDET_CSF=80)
!      parameter(MOTYPE=5,MPARMD=770)
!      parameter(NSPLIN=1001,MORDJ=9,MPARMJ=70,ML_BAS=4,NRAD=1001,RADMAX=8.d0)

! Ryo, Masayoshi:
!     parameter(MCENT=1,MCTYPE=1,MELEC=180,MELECUD=90,MORB=90,MORB_OCC=90)
!     parameter(MBASIS=180,MBASIS_CTYPE=180)
!     parameter(MDET=100,MDETUD=10,MCSF=30,MDET_CSF=20)
!     parameter(MOTYPE=5,MPARMD=170)
!     parameter(NSPLIN=1001,MORDJ=7,MPARMJ=70,ML_BAS=6,NRAD=1001,RADMAX=8.d0)


      parameter(NAX=25,NAK1=40,NAK2=1)
      parameter(MMAT_DIM=(MELECUD*MELECUD),MMAT_DIM2=(MELEC*(MELEC-1))/2)
      parameter(MORDJ1=MORDJ+1,delri=(NRAD-1)/RADMAX)
