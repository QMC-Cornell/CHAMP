c     implicit real*8(a-h,o-z)
c     parameter(MGRID=200)
c     dimension f(0:MGRID-1,0:MGRID-1,0:MGRID-1)

c     pi=4*datan(1.d0)
c     read(5,*) ngrid,rmax
c     do 10 i=0,ngrid-1
c     do 10 j=0,ngrid-1
c     do 10 k=0,ngrid-1
c  10   f(i,j,k)=sin(2*pi*i/ngrid)*sin(2*pi*j/ngrid)*sin(2*pi*k/ngrid)

ccFirst check on the grid pts themselves
c     do 15 i=0,ngrid-1
c     do 15 j=0,ngrid-1
c     do 15 k=0,ngrid-1
c       x=i/dfloat(ngrid)
c       y=j/dfloat(ngrid)
c       z=k/dfloat(ngrid)
c       xi=ngrid*x
c       yi=ngrid*y
c       zi=ngrid*z
c       fi=interpol_orb(MGRID,MGRID,MGRID,ngrid,ngrid,ngrid,f,xi,yi,zi)
c  15   write(6,'(''xi,yi,zi,fi,f'',3f8.4,9f12.8)') xi,yi,zi,fi,
c    &  sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)

ccNow check on other pts.
c     do 20 i=1,2
c     do 20 j=1,2
c     do 20 k=1,2
c     x=mod(i*sqrt(17.d0),10.d0*rmax)-5*rmax
c     y=mod(j*sqrt(19.d0),10.d0*rmax)-5*rmax
c     z=mod(k*sqrt(21.d0),10.d0*rmax)-5*rmax
c     if(x.ge.0) then
c       xi=ngrid*x
c      else
c       xi=(x-int(x)+1)*ngrid
c     endif
c     if(y.ge.0) then
c       yi=ngrid*y
c      else
c       yi=(y-int(y)+1)*ngrid
c     endif
c     if(z.ge.0) then
c       zi=ngrid*z
c      else
c       zi=(z-int(z)+1)*ngrid
c     endif
c     fi=interpol_orb(MGRID,MGRID,MGRID,ngrid,ngrid,ngrid,f,xi,yi,zi)
c  20 write(6,'(''xi,yi,zi,fi,f'',3f8.4,9f12.8)') xi,yi,zi,fi,
c    &  sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
c     stop
c     end
c-----------------------------------------------------------------------

      subroutine interpol_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
c Written by Cyrus Umrigar

      implicit real*8(a-h,o-z)

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      dimension orb(*),dorb(3,*),ddorb(*)

      if(abs(inum_orb).eq.1) then
        call interpol1_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
       elseif(abs(inum_orb).eq.4) then
        call interpol4_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine interpol1_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c Evaluate orbitals, their gradient and Laplacian by looking up closest pt.
c This routine is never used in practice.
c The mesh pts. on which the function values, f, are given, are assumed
c to be at 0,1,2,...nx-1, and similarly for y and z.
c So, pt. nx is the same as pt. 0.
c It is easy to write a routine for general order and unequally-spaced
c grid, but this routine is specialized to make it fast.
c Warning: This routine assumes fn. is periodic and tabulated on entire period.
c If function is tabulated on half-period then this routine needs modification.

      use orbital_grid_mod
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'numorb.h'

      common /dim/ ndim
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ngrid_orbx,ngrid_orby,ngrid_orbz

      dimension orb(*),dorb(3,*),ddorb(*)

c     if(xi.ge.0) then
c       x=xi-int(xi)
c      else
c       x=xi-int(xi)+1
c     endif
      ix=int(mod(xi,dfloat(nx)))
      if(ix.lt.0) ix=ix+nx

c     if(yi.ge.0) then
c       y=yi-int(yi)
c      else
c       y=yi-int(yi)+1
c     endif
      iy=int(mod(yi,dfloat(ny)))
      if(iy.lt.0) iy=iy+ny

c     if(zi.ge.0) then
c       z=zi-int(zi)
c      else
c       z=zi-int(zi)+1
c     endif
      iz=int(mod(zi,dfloat(nz)))
      if(iz.lt.0) iz=iz+nz

      do 10 iorb=1,norb
        orb(iorb)=orb_num(iorb,ix,iy,iz)
        ddorb(iorb)=ddorb_num(iorb,ix,iy,iz)
        do 10 k=1,ndim
   10     dorb(k,iorb)=dorb_num(k,iorb,ix,iy,iz)

      return
      end
c-----------------------------------------------------------------------

      subroutine interpol4_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c Evaluate orbitals, their gradient and Laplacian by a cubic (4 mesh pts) Lagrange
c interpolation on an equally-spaced 3D grid with periodic boundary conditions.
c The mesh pts. on which the function values, f, are given, are assumed
c to be at 0,1,2,...nx-1, and similarly for y and z.
c So, pt. nx is the same as pt. 0.
c It is easy to write a routine for general order and unequally-spaced
c grid, but this routine is specialized to make it fast.
c If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
c If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'numorb.h'
      parameter(half=0.5d0,sixth=1.d0/6.d0)

      common /dim/ ndim
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /periodic2/ rkvec_shift_latt(3)

      dimension orb(*),dorb(3,*),ddorb(*),
     &orb1(MORB),dorb1(3,MORB),ddorb1(MORB),term(4),iz(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

c input coordinates are positive
c     if(zi.ge.0) then
        z=zi-int(zi)
c      else
c       z=zi-int(zi)+1
c     endif
c     iz(2)=int(mod(zi,dfloat(nz)))
c     if(iz(2).lt.0) iz(2)=iz(2)+nz
      iz(2)=int(zi)
      iz(1)=iz(2)-1
      if(iz(1).lt.0) then
        iz(1)=iz(1)+nz
        if(rkvec_shift_latt(3).ne.0.d0) isgn(1)=-1
      endif
c     iz(3)=mod(iz(2)+1,nz)
c     iz(4)=mod(iz(2)+2,nz)
      iz(3)=iz(2)+1
      if(iz(3).ge.nz) then
        iz(3)=iz(3)-nz
        if(rkvec_shift_latt(3).ne.0.d0) isgn(3)=-1
      endif
      iz(4)=iz(2)+2
      if(iz(4).ge.nz) then
        iz(4)=iz(4)-nz
        if(rkvec_shift_latt(3).ne.0.d0) isgn(4)=-1
      endif

c     write(6,'(''iz'',4i2)') (iz(j),j=1,4)

      term(1)=-sixth*z*(z-1)*(z-2)*isgn(1)
      term(2)=half*(z+1)*(z-1)*(z-2)
      term(3)=-half*(z+1)*z*(z-2)*isgn(3)
      term(4)=sixth*(z+1)*z*(z-1)*isgn(4)

      do 5 iorb=1,norb
        orb(iorb)=0
        ddorb(iorb)=0
        do 5 k=1,ndim
    5     dorb(k,iorb)=0

      do 10 j=1,4
        call interpol4_y(nx,ny,xi,yi,iz(j),orb1,dorb1,ddorb1)
        do 10 iorb=1,norb
          orb(iorb)=orb(iorb)+term(j)*orb1(iorb)
          ddorb(iorb)=ddorb(iorb)+term(j)*ddorb1(iorb)
          do 10 k=1,ndim
   10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb1(k,iorb)

      return
      end
c-----------------------------------------------------------------------

      subroutine interpol4_y(nx,ny,xi,yi,iz,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
c If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'numorb.h'
      parameter(half=0.5d0,sixth=1.d0/6.d0)

      common /dim/ ndim
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /periodic2/ rkvec_shift_latt(3)

      dimension orb(*),dorb(3,*),ddorb(*),
     &orb1(MORB),dorb1(3,MORB),ddorb1(MORB),term(4),iy(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

c input coordinates are positive
c     if(yi.ge.0) then
        y=yi-int(yi)
c      else
c       y=yi-int(yi)+1
c     endif
c     iy(2)=int(mod(yi,dfloat(ny)))
c     if(iy(2).lt.0) iy(2)=iy(2)+ny
      iy(2)=int(yi)
      iy(1)=iy(2)-1
      if(iy(1).lt.0) then
        iy(1)=iy(1)+ny
        if(rkvec_shift_latt(2).ne.0.d0) isgn(1)=-1
      endif
c     iy(3)=mod(iy(2)+1,ny)
c     iy(4)=mod(iy(2)+2,ny)
      iy(3)=iy(2)+1
      if(iy(3).ge.ny) then
        iy(3)=iy(3)-ny
        if(rkvec_shift_latt(2).ne.0.d0) isgn(3)=-1
      endif
      iy(4)=iy(2)+2
      if(iy(4).ge.ny) then
        iy(4)=iy(4)-ny
        if(rkvec_shift_latt(2).ne.0.d0) isgn(4)=-1
      endif

c     write(6,'(''iy'',4i2)') (iy(j),j=1,4)

      term(1)=-sixth*y*(y-1)*(y-2)*isgn(1)
      term(2)=half*(y+1)*(y-1)*(y-2)
      term(3)=-half*(y+1)*y*(y-2)*isgn(3)
      term(4)=sixth*(y+1)*y*(y-1)*isgn(4)

      do 5 iorb=1,norb
        orb(iorb)=0
        ddorb(iorb)=0
        do 5 k=1,ndim
    5     dorb(k,iorb)=0

      do 10 j=1,4
        call interpol4_x(nx,xi,iy(j),iz,orb1,dorb1,ddorb1)
        do 10 iorb=1,norb
          orb(iorb)=orb(iorb)+term(j)*orb1(iorb)
          ddorb(iorb)=ddorb(iorb)+term(j)*ddorb1(iorb)
          do 10 k=1,ndim
   10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb1(k,iorb)

      return
      end
c-----------------------------------------------------------------------

      subroutine interpol4_x(nx,xi,iy,iz,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
c If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      use orbital_grid_mod
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'numorb.h'
      parameter(half=0.5d0,sixth=1.d0/6.d0)

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ngrid_orbx,ngrid_orby,ngrid_orbz
      common /periodic2/ rkvec_shift_latt(3)

      dimension orb(*),dorb(3,*),ddorb(*)
     &,term(4),ix(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

c input coordinates are positive
c     if(xi.ge.0) then
        x=xi-int(xi)
c      else
c       x=xi-int(xi)+1
c     endif
c     ix(2)=int(mod(xi,dfloat(nx)))
c     if(ix(2).lt.0) ix(2)=ix(2)+nx
      ix(2)=int(xi)
      ix(1)=ix(2)-1
      if(ix(1).lt.0) then
        ix(1)=ix(1)+nx
        if(rkvec_shift_latt(1).ne.0.d0) isgn(1)=-1
      endif
c     ix(3)=mod(ix(2)+1,nx)
c     ix(4)=mod(ix(2)+2,nx)
      ix(3)=ix(2)+1
      if(ix(3).ge.nx) then
        ix(3)=ix(3)-nx
        if(rkvec_shift_latt(1).ne.0.d0) isgn(3)=-1
      endif
      ix(4)=ix(2)+2
      if(ix(4).ge.nx) then
        ix(4)=ix(4)-nx
        if(rkvec_shift_latt(1).ne.0.d0) isgn(4)=-1
      endif

c     write(6,'(''ix'',4i2)') (ix(j),j=1,4)

      term(1)=-sixth*x*(x-1)*(x-2)*isgn(1)
      term(2)=half*(x+1)*(x-1)*(x-2)
      term(3)=-half*(x+1)*x*(x-2)*isgn(3)
      term(4)=sixth*(x+1)*x*(x-1)*isgn(4)

c Timing tests done with g77-2.95 O(3) on newton (file server)
c Various possible codings are identified with 1), 1a), 2), etc.

c 1) This is the straightforward reference 7.90:
c     do 10 iorb=1,norb
c       orb(iorb)=0
c       ddorb(iorb)=0
c       do 5 k=1,ndim
c   5     dorb(k,iorb)=0
c       do 10 j=1,4
c         orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)
c         ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)
c         do 10 k=1,ndim
c  10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


c 1a) this is 20% faster 6.76
c     do 10 iorb=1,norb
c       orb(iorb)=term(1)*orb_num(iorb,ix(1),iy,iz)
c       ddorb(iorb)=term(1)*ddorb_num(iorb,ix(1),iy,iz)
c       do 5 k=1,ndim
c   5     dorb(k,iorb)=term(1)*dorb_num(k,iorb,ix(1),iy,iz)
c       do 10 j=2,4
c         orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)
c         ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)
c         do 10 k=1,ndim
c  10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


c 2) this is 15% faster 6.99
c     do 5 iorb=1,norb
c       orb(iorb)=0
c       ddorb(iorb)=0
c       do 5 k=1,ndim
c   5     dorb(k,iorb)=0

c     do 10 j=1,4
c       do 10 iorb=1,norb
c         orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)
c         ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)
c         do 10 k=1,ndim
c  10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


c 2a) this is 20% faster 6.43
c     do 5 iorb=1,norb
c       orb(iorb)=term(1)*orb_num(iorb,ix(1),iy,iz)
c       ddorb(iorb)=term(1)*ddorb_num(iorb,ix(1),iy,iz)
c       do 5 k=1,ndim
c   5     dorb(k,iorb)=term(1)*dorb_num(k,iorb,ix(1),iy,iz)

c     do 10 j=2,4
c       do 10 iorb=1,norb
c         orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)
c         ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)
c         do 10 k=1,ndim
c  10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


c 2b) this is 25% faster 5.68
      do 5 iorb=1,norb
        orb(iorb)=term(1)*orb_num(iorb,ix(1),iy,iz)
        ddorb(iorb)=term(1)*ddorb_num(iorb,ix(1),iy,iz)
          dorb(1,iorb)=term(1)*dorb_num(1,iorb,ix(1),iy,iz)
          dorb(2,iorb)=term(1)*dorb_num(2,iorb,ix(1),iy,iz)
    5     dorb(3,iorb)=term(1)*dorb_num(3,iorb,ix(1),iy,iz)

      do 10 j=2,4
        do 10 iorb=1,norb
          orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)
          ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)
            dorb(1,iorb)=dorb(1,iorb)+term(j)*dorb_num(1,iorb,ix(j),iy,iz)
            dorb(2,iorb)=dorb(2,iorb)+term(j)*dorb_num(2,iorb,ix(j),iy,iz)
   10       dorb(3,iorb)=dorb(3,iorb)+term(j)*dorb_num(3,iorb,ix(j),iy,iz)


c 3) this is about 20% slower (this routine only) 9.12
c     do 5 iorb=1,norb
c   5   orb(iorb)=0

c     do 6 iorb=1,norb
c   6   ddorb(iorb)=0

c     do 7 iorb=1,norb
c       do 7 k=1,ndim
c   7     dorb(k,iorb)=0

c     do 10 j=1,4
c       do 10 iorb=1,norb
c  10     orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)

c     do 20 j=1,4
c       do 20 iorb=1,norb
c  20     ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)

c     do 30 j=1,4
c       do 30 iorb=1,norb
c         do 30 k=1,ndim
c  30       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


c 4) 10% faster 7.31
c     do 10 iorb=1,norb
c       orb(iorb)=0
c       do 10 j=1,4
c  10     orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)

c     do 20 iorb=1,norb
c       ddorb(iorb)=0
c       do 20 j=1,4
c  20     ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)

c     do 40 iorb=1,norb
c       do 30 k=1,ndim
c  30     dorb(k,iorb)=0
c       do 40 j=1,4
c         do 40 k=1,ndim
c  40       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


c 4a) this is 20% faster 6.43
c     do 10 iorb=1,norb
c       orb(iorb)=term(1)*orb_num(iorb,ix(1),iy,iz)
c       do 10 j=2,4
c  10     orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)

c     do 20 iorb=1,norb
c       ddorb(iorb)=term(1)*ddorb_num(iorb,ix(1),iy,iz)
c       do 20 j=2,4
c  20     ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)

c     do 40 iorb=1,norb
c       do 30 k=1,ndim
c  30     dorb(k,iorb)=term(1)*dorb_num(k,iorb,ix(1),iy,iz)
c       do 40 j=2,4
c         do 40 k=1,ndim
c  40       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)

      return
      end
c-----------------------------------------------------------------------

      subroutine interpol_orbe(nx,ny,nz,xi,yi,zi,orb)
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)

      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt

      dimension orb(*)

      if(abs(inum_orb).eq.1) then
        call interpol1_orbe(nx,ny,nz,xi,yi,zi,orb)
       elseif(abs(inum_orb).eq.4) then
        call interpol4_orbe(nx,ny,nz,xi,yi,zi,orb)
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine interpol1_orbe(nx,ny,nz,xi,yi,zi,orb)
c Written by Cyrus Umrigar
c Evaluate orbitals by looking up closest pt.
c This routine is never used in practice.
c The mesh pts. on which the function values, f, are given, are assumed
c to be at 0,1,2,...nx-1, and similarly for y and z.
c So, pt. nx is the same as pt. 0.
c It is easy to write a routine for general order and unequally-spaced
c grid, but this routine is specialized to make it fast.
c Warning: This routine assumes fn. is periodic and tabulated on entire period.
c If function is tabulated on half-period then this routine needs modification.

      use orbital_grid_mod
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'numorb.h'

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ngrid_orbx,ngrid_orby,ngrid_orbz

      dimension orb(*)

c     if(xi.ge.0) then
c       x=xi-int(xi)
c      else
c       x=xi-int(xi)+1
c     endif
      ix=int(mod(xi,dfloat(nx)))
      if(ix.lt.0) ix=ix+nx

c     if(yi.ge.0) then
c       y=yi-int(yi)
c      else
c       y=yi-int(yi)+1
c     endif
      iy=int(mod(yi,dfloat(ny)))
      if(iy.lt.0) iy=iy+ny

c     if(zi.ge.0) then
c       z=zi-int(zi)
c      else
c       z=zi-int(zi)+1
c     endif
      iz=int(mod(zi,dfloat(nz)))
      if(iz.lt.0) iz=iz+nz

      do 10 iorb=1,norb
   10   orb(iorb)=orb_num(iorb,ix,iy,iz)

      return
      end
c-----------------------------------------------------------------------
      subroutine interpol4_orbe(nx,ny,nz,xi,yi,zi,orb)
c Written by Cyrus Umrigar
c Evaluate orbitals by a cubic (4 mesh pts) Lagrange
c interpolation on an equally-spaced 3D grid with periodic boundary conditions.
c The mesh pts. on which the function values, f, are given, are assumed
c to be at 0,1,2,...nx-1, and similarly for y and z.
c So, pt. nx is the same as pt. 0.
c It is easy to write a routine for general order and unequally-spaced
c grid, but this routine is specialized to make it fast.
c If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
c If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'numorb.h'
      parameter(half=0.5d0,sixth=1.d0/6.d0)

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /periodic2/ rkvec_shift_latt(3)

      dimension orb(*),
     &orb1(MORB),term(4),iz(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

c input coordinates are positive
c     if(zi.ge.0) then
        z=zi-int(zi)
c      else
c       z=zi-int(zi)+1
c     endif
c     iz(2)=int(mod(zi,dfloat(nz)))
c     if(iz(2).lt.0) iz(2)=iz(2)+nz
      iz(2)=int(zi)
      iz(1)=iz(2)-1
      if(iz(1).lt.0) then
        iz(1)=iz(1)+nz
        if(rkvec_shift_latt(3).ne.0.d0) isgn(1)=-1
      endif
c     iz(3)=mod(iz(2)+1,nz)
c     iz(4)=mod(iz(2)+2,nz)
      iz(3)=iz(2)+1
      if(iz(3).ge.nz) then
        iz(3)=iz(3)-nz
        if(rkvec_shift_latt(3).ne.0.d0) isgn(3)=-1
      endif
      iz(4)=iz(2)+2
      if(iz(4).ge.nz) then
        iz(4)=iz(4)-nz
        if(rkvec_shift_latt(3).ne.0.d0) isgn(4)=-1
      endif

c     write(6,'(''iz'',4i2)') (iz(j),j=1,4)

      term(1)=-sixth*z*(z-1)*(z-2)*isgn(1)
      term(2)=half*(z+1)*(z-1)*(z-2)
      term(3)=-half*(z+1)*z*(z-2)*isgn(3)
      term(4)=sixth*(z+1)*z*(z-1)*isgn(4)

      do 5 iorb=1,norb
    5   orb(iorb)=0

      do 10 j=1,4
        call interpol4_ye(nx,ny,xi,yi,iz(j),orb1)
        do 10 iorb=1,norb
   10     orb(iorb)=orb(iorb)+term(j)*orb1(iorb)

      return
      end
c-----------------------------------------------------------------------

      subroutine interpol4_ye(nx,ny,xi,yi,iz,orb)
c Written by Cyrus Umrigar
c If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
c If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'numorb.h'
      parameter(half=0.5d0,sixth=1.d0/6.d0)

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /periodic2/ rkvec_shift_latt(3)

      dimension orb(*),
     &orb1(MORB),term(4),iy(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

c input coordinates are positive
c     if(yi.ge.0) then
        y=yi-int(yi)
c      else
c       y=yi-int(yi)+1
c     endif
c     iy(2)=int(mod(yi,dfloat(ny)))
c     if(iy(2).lt.0) iy(2)=iy(2)+ny
      iy(2)=int(yi)
      iy(1)=iy(2)-1
      if(iy(1).lt.0) then
        iy(1)=iy(1)+ny
        if(rkvec_shift_latt(2).ne.0.d0) isgn(1)=-1
      endif
c     iy(3)=mod(iy(2)+1,ny)
c     iy(4)=mod(iy(2)+2,ny)
      iy(3)=iy(2)+1
      if(iy(3).ge.ny) then
        iy(3)=iy(3)-ny
        if(rkvec_shift_latt(2).ne.0.d0) isgn(3)=-1
      endif
      iy(4)=iy(2)+2
      if(iy(4).ge.ny) then
        iy(4)=iy(4)-ny
        if(rkvec_shift_latt(2).ne.0.d0) isgn(4)=-1
      endif

c     write(6,'(''iy'',4i2)') (iy(j),j=1,4)

      term(1)=-sixth*y*(y-1)*(y-2)*isgn(1)
      term(2)=half*(y+1)*(y-1)*(y-2)
      term(3)=-half*(y+1)*y*(y-2)*isgn(3)
      term(4)=sixth*(y+1)*y*(y-1)*isgn(4)

      do 5 iorb=1,norb
    5   orb(iorb)=0

      do 10 j=1,4
        call interpol4_xe(nx,xi,iy(j),iz,orb1)
        do 10 iorb=1,norb
   10     orb(iorb)=orb(iorb)+term(j)*orb1(iorb)

      return
      end
c-----------------------------------------------------------------------

      subroutine interpol4_xe(nx,xi,iy,iz,orb)
c Written by Cyrus Umrigar
c If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
c If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      use orbital_grid_mod
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'numorb.h'
      parameter(half=0.5d0,sixth=1.d0/6.d0)

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ngrid_orbx,ngrid_orby,ngrid_orbz
      common /periodic2/ rkvec_shift_latt(3)

      dimension orb(*)
     &,term(4),ix(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

c input coordinates are positive
c     if(xi.ge.0) then
        x=xi-int(xi)
c      else
c       x=xi-int(xi)+1
c     endif
c     ix(2)=int(mod(xi,dfloat(nx)))
c     if(ix(2).lt.0) ix(2)=ix(2)+nx
      ix(2)=int(xi)
      ix(1)=ix(2)-1
      if(ix(1).lt.0) then
        ix(1)=ix(1)+nx
        if(rkvec_shift_latt(1).ne.0.d0) isgn(1)=-1
      endif
c     ix(3)=mod(ix(2)+1,nx)
c     ix(4)=mod(ix(2)+2,nx)
      ix(3)=ix(2)+1
      if(ix(3).ge.nx) then
        ix(3)=ix(3)-nx
        if(rkvec_shift_latt(1).ne.0.d0) isgn(3)=-1
      endif
      ix(4)=ix(2)+2
      if(ix(4).ge.nx) then
        ix(4)=ix(4)-nx
        if(rkvec_shift_latt(1).ne.0.d0) isgn(4)=-1
      endif

c     write(6,'(''ix'',4i2)') (ix(j),j=1,4)

      term(1)=-sixth*x*(x-1)*(x-2)*isgn(1)
      term(2)=half*(x+1)*(x-1)*(x-2)
      term(3)=-half*(x+1)*x*(x-2)*isgn(3)
      term(4)=sixth*(x+1)*x*(x-1)*isgn(4)

      do 10 iorb=1,norb
        orb(iorb)=term(1)*orb_num(iorb,ix(1),iy,iz)
        do 10 j=2,4
   10     orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)

      return
      end
