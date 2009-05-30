! Used to set dimensions of arrays.
! MORB           orbitals
! MORB_OCC       occupied orbitals
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

! Cyrus:
! I have set MPARMD>= MCSF but if this is not an optimization run it could be less.
!      parameter(MORB_OCC=24, MORB=MAX(64,MORB_OCC))
!!     parameter(MCSF=165, MDET_CSF=140, MDETUD=190, MDET=MAX(660,MDET_CSF,MDETUD))
!      parameter(MCSF=60, MDET_CSF=132, MDETUD=162, MDET=MAX(350,MDET_CSF,MDETUD))
!      parameter(MOTYPE=5, MPARMD=MAX(770,MCSF), ML_BAS=4)
       parameter(MORB_OCC=42, MORB=MAX(80,MORB_OCC))
!      parameter(MCSF=165, MDET_CSF=140, MDETUD=190, MDET=MAX(660,MDET_CSF,MDETUD))
!      parameter(MCSF=60, MDET_CSF=132, MDETUD=162, MDET=MAX(350,MDET_CSF,MDETUD))
!      parameter(MCSF=20, MDET_CSF=32, MDETUD=10, MDET=MAX(50,MDET_CSF,MDETUD))
!      parameter(MCSF=10, MDET_CSF=10, MDETUD=10, MDET=MAX(20,MDET_CSF,MDETUD))
!      parameter(MOTYPE=5, MPARMD=MAX(100,MCSF), ML_BAS=4)

! Julien:
!      parameter(MORB_OCC=10, MORB=MAX(68,MORB_OCC))
       parameter(MCSF=165, MDET_CSF=80, MDETUD=100, MDET=MAX(660,MDET_CSF,MDETUD))
       parameter(MOTYPE=5, MPARMD=MAX(170,MCSF), ML_BAS=4)

! Ryo, Masayoshi:
!     parameter(MORB_OCC=90, MORB=MAX(90,MORB_OCC))
!     parameter(MCSF=30, MDET_CSF=20, MDETUD=10, MDET=MAX(100,MDET_CSF,MDETUD))
!     parameter(MOTYPE=5, MPARMD=MAX(170,MCSF), ML_BAS=6)

! Wissam:
!     parameter(MORB_OCC=34, MORB=MAX(50,MORB_OCC))
!     parameter(MCSF=2, MDET_CSF=2, MDETUD=30, MDET=MAX(2,MDET_CSF,MDETUD))
!     parameter(MOTYPE=5, MPARMD=MAX(170,MCSF), ML_BAS=4)

      parameter(MORDJ=9, MPARMJ=70, NSPLIN=1001, NRAD=1001, RADMAX=8.d0)
      parameter(NAX=50,NAK1=40,NAK2=1)
      parameter(MORDJ1=MORDJ+1,delri=(NRAD-1)/RADMAX)
