! DBLMIN, DBLMAX should be the smallest and largest double precision numbers
! that can be represented, but approximate values are OK.
!JT   parameter(DBLMIN=1.d-300,DBLMAX=1.d+300)
      parameter(DBLMIN=1.d-30,DBLMAX=1.d+30)

      parameter(MORDJ=9, NRAD=1001, RADMAX=8.d0)
      parameter(NAX=50,NAK1=40,NAK2=1)
      parameter(MORDJ1=MORDJ+1,delri=(NRAD-1)/RADMAX)
