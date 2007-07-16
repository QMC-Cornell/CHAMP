c Stefan Goedecker's routine to get atom positions for supercells.
c computes positions of the atoms in diamond structure
c and writes them in a format that Jose-Luis Martin's PW code can read.
      implicit real*8(a-h,o-z)
      parameter (MCENT=500)
      dimension rxyz(3,MCENT),alat(3)

      write(6,*) 'input ncell1,ncell2,ncell3'
      read(5,*) ncell1,ncell2,ncell3
      write(6,*) 'input the number of centers "ncent"'
      read(5,*) ncent
      call diamond(alat,ncell1,ncell2,ncell3,ncent,rxyz)

      do 10 i=1,ncent
   10   write(6,'(3(3x,f12.6))') (rxyz(k,i),k=1,3)

      stop
      end
c-----------------------------------------------------------------------

      subroutine diamond(alat,ncell1,ncell2,ncell3,nat,rxyz)
c Stefan Goedecker's routine to get atom positions for supercells.
c computes positions of the atoms in diamond structure
      implicit real*8(a-h,o-z)

c Silicon simple cubic lattice constant in Angstroem
c     parameter(acell=5.42981d0)
      parameter(acell=.5d0)

      dimension rxyz(3,nat),alat(3)

      if (nat.ne.8*ncell1*ncell2*ncell3) stop 'wrong nat'
      alat(1)=acell*ncell1
      alat(2)=acell*ncell2
      alat(3)=acell*ncell3

      acell1=acell
      acell2=acell
      acell3=acell
c perfect lattice position routine
      do 30 j3=0,ncell3-1
       do 30 j2=0,ncell2-1
        do 30 j1=0,ncell1-1
           jj=8*(j1+ncell1*j2+ncell2*ncell1*j3)

           rxyz(1,jj+1)=acell1*j1
           rxyz(2,jj+1)=acell2*j2
           rxyz(3,jj+1)=acell3*j3

           rxyz(1,jj+2)=acell1*j1 +.5d0*acell1
           rxyz(2,jj+2)=acell2*j2 +.5d0*acell2
           rxyz(3,jj+2)=acell3*j3

           rxyz(1,jj+3)=acell1*j1 +.5d0*acell1
           rxyz(2,jj+3)=acell2*j2
           rxyz(3,jj+3)=acell3*j3 +.5d0*acell3

           rxyz(1,jj+4)=acell1*j1
           rxyz(2,jj+4)=acell2*j2 +.5d0*acell2
           rxyz(3,jj+4)=acell3*j3 +.5d0*acell3

           rxyz(1,jj+5)=acell1*j1 + .25d0*acell1
           rxyz(2,jj+5)=acell2*j2 + .25d0*acell2
           rxyz(3,jj+5)=acell3*j3 + .25d0*acell3

           rxyz(1,jj+6)=acell1*j1 + .25d0*acell1 +.5d0*acell1
           rxyz(2,jj+6)=acell2*j2 + .25d0*acell2 +.5d0*acell2
           rxyz(3,jj+6)=acell3*j3 + .25d0*acell3

           rxyz(1,jj+7)=acell1*j1 + .25d0*acell1 +.5d0*acell1
           rxyz(2,jj+7)=acell2*j2 + .25d0*acell2
           rxyz(3,jj+7)=acell3*j3 + .25d0*acell3 +.5d0*acell3

           rxyz(1,jj+8)=acell1*j1 + .25d0*acell1
           rxyz(2,jj+8)=acell2*j2 + .25d0*acell2 +.5d0*acell2
           rxyz(3,jj+8)=acell3*j3 + .25d0*acell3 +.5d0*acell3
30    continue

      return
      end
