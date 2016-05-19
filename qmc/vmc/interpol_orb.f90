!     implicit real*8(a-h,o-z)
!     parameter(MGRID=200)
!     dimension f(0:MGRID-1,0:MGRID-1,0:MGRID-1)

!     pi=4*datan(1.d0)
!     read(5,*) ngrid,rmax
!     do 10 i=0,ngrid-1
!     do 10 j=0,ngrid-1
!     do 10 k=0,ngrid-1
!  10   f(i,j,k)=sin(2*pi*i/ngrid)*sin(2*pi*j/ngrid)*sin(2*pi*k/ngrid)

!cFirst check on the grid pts themselves
!     do 15 i=0,ngrid-1
!     do 15 j=0,ngrid-1
!     do 15 k=0,ngrid-1
!       x=i/dfloat(ngrid)
!       y=j/dfloat(ngrid)
!       z=k/dfloat(ngrid)
!       xi=ngrid*x
!       yi=ngrid*y
!       zi=ngrid*z
!       fi=interpol_orb(MGRID,MGRID,MGRID,ngrid,ngrid,ngrid,f,xi,yi,zi)
!  15   write(6,'(''xi,yi,zi,fi,f'',3f8.4,9f12.8)') xi,yi,zi,fi,
!    &  sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)

!cNow check on other pts.
!     do 20 i=1,2
!     do 20 j=1,2
!     do 20 k=1,2
!     x=mod(i*sqrt(17.d0),10.d0*rmax)-5*rmax
!     y=mod(j*sqrt(19.d0),10.d0*rmax)-5*rmax
!     z=mod(k*sqrt(21.d0),10.d0*rmax)-5*rmax
!     if(x.ge.0) then
!       xi=ngrid*x
!      else
!       xi=(x-int(x)+1)*ngrid
!     endif
!     if(y.ge.0) then
!       yi=ngrid*y
!      else
!       yi=(y-int(y)+1)*ngrid
!     endif
!     if(z.ge.0) then
!       zi=ngrid*z
!      else
!       zi=(z-int(z)+1)*ngrid
!     endif
!     fi=interpol_orb(MGRID,MGRID,MGRID,ngrid,ngrid,ngrid,f,xi,yi,zi)
!  20 write(6,'(''xi,yi,zi,fi,f'',3f8.4,9f12.8)') xi,yi,zi,fi,
!    &  sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
!     stop
!     end
!-----------------------------------------------------------------------

      subroutine interpol_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
! Written by Cyrus Umrigar

      use contr2_mod
      implicit real*8(a-h,o-z)


      dimension orb(*),dorb(3,*),ddorb(*)

      if(abs(inum_orb).eq.1) then
        call interpol1_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
       elseif(abs(inum_orb).eq.4) then
        call interpol4_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine interpol1_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
! Written by Cyrus Umrigar
! Evaluate orbitals, their gradient and Laplacian by looking up closest pt.
! This routine is never used in practice.
! The mesh pts. on which the function values, f, are given, are assumed
! to be at 0,1,2,...nx-1, and similarly for y and z.
! So, pt. nx is the same as pt. 0.
! It is easy to write a routine for general order and unequally-spaced
! grid, but this routine is specialized to make it fast.
! Warning: This routine assumes fn. is periodic and tabulated on entire period.
! If function is tabulated on half-period then this routine needs modification.

      use orbital_grid_mod
      use coefs_mod
      use dim_mod
      implicit real*8(a-h,o-z)

      dimension orb(*),dorb(3,*),ddorb(*)

!     if(xi.ge.0) then
!       x=xi-int(xi)
!      else
!       x=xi-int(xi)+1
!     endif
      ix=int(mod(xi,dfloat(nx)))
      if(ix.lt.0) ix=ix+nx

!     if(yi.ge.0) then
!       y=yi-int(yi)
!      else
!       y=yi-int(yi)+1
!     endif
      iy=int(mod(yi,dfloat(ny)))
      if(iy.lt.0) iy=iy+ny

!     if(zi.ge.0) then
!       z=zi-int(zi)
!      else
!       z=zi-int(zi)+1
!     endif
      iz=int(mod(zi,dfloat(nz)))
      if(iz.lt.0) iz=iz+nz

      do 10 iorb=1,norb
        orb(iorb)=orb_num(iorb,ix,iy,iz)
        ddorb(iorb)=ddorb_num(iorb,ix,iy,iz)
        do 10 k=1,ndim
   10     dorb(k,iorb)=dorb_num(k,iorb,ix,iy,iz)

      return
      end
!-----------------------------------------------------------------------

      subroutine interpol4_orb(nx,ny,nz,xi,yi,zi,orb,dorb,ddorb)
! Written by Cyrus Umrigar
! Evaluate orbitals, their gradient and Laplacian by a cubic (4 mesh pts) Lagrange
! interpolation on an equally-spaced 3D grid with periodic boundary conditions.
! The mesh pts. on which the function values, f, are given, are assumed
! to be at 0,1,2,...nx-1, and similarly for y and z.
! So, pt. nx is the same as pt. 0.
! It is easy to write a routine for general order and unequally-spaced
! grid, but this routine is specialized to make it fast.
! If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
! If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      use constants_mod
      use coefs_mod
      use dim_mod
      use periodic2_mod
      implicit real*8(a-h,o-z)

      dimension orb(*),dorb(3,*),ddorb(*)
      dimension orb1(norb),dorb1(3,norb),ddorb1(norb),term(4),iz(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

! input coordinates are positive
!     if(zi.ge.0) then
        z=zi-int(zi)
!      else
!       z=zi-int(zi)+1
!     endif
!     iz(2)=int(mod(zi,dfloat(nz)))
!     if(iz(2).lt.0) iz(2)=iz(2)+nz
      iz(2)=int(zi)
      iz(1)=iz(2)-1
      if(iz(1).lt.0) then
        iz(1)=iz(1)+nz
        if(rkvec_shift_latt(3).ne.0.d0) isgn(1)=-1
      endif
!     iz(3)=mod(iz(2)+1,nz)
!     iz(4)=mod(iz(2)+2,nz)
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

!     write(6,'(''iz'',4i2)') (iz(j),j=1,4)

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
!-----------------------------------------------------------------------

      subroutine interpol4_y(nx,ny,xi,yi,iz,orb,dorb,ddorb)
! Written by Cyrus Umrigar
! If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
! If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      use constants_mod
      use coefs_mod
      use dim_mod
      use periodic2_mod
      implicit real*8(a-h,o-z)

      dimension orb(*),dorb(3,*),ddorb(*)
      dimension orb1(norb),dorb1(3,norb),ddorb1(norb),term(4),iy(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

! input coordinates are positive
!     if(yi.ge.0) then
        y=yi-int(yi)
!      else
!       y=yi-int(yi)+1
!     endif
!     iy(2)=int(mod(yi,dfloat(ny)))
!     if(iy(2).lt.0) iy(2)=iy(2)+ny
      iy(2)=int(yi)
      iy(1)=iy(2)-1
      if(iy(1).lt.0) then
        iy(1)=iy(1)+ny
        if(rkvec_shift_latt(2).ne.0.d0) isgn(1)=-1
      endif
!     iy(3)=mod(iy(2)+1,ny)
!     iy(4)=mod(iy(2)+2,ny)
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

!     write(6,'(''iy'',4i2)') (iy(j),j=1,4)

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
!-----------------------------------------------------------------------

      subroutine interpol4_x(nx,xi,iy,iz,orb,dorb,ddorb)
! Written by Cyrus Umrigar
! If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
! If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      use constants_mod
      use orbital_grid_mod
      use coefs_mod
      use periodic2_mod
      implicit real*8(a-h,o-z)

      dimension orb(*),dorb(3,*),ddorb(*) &
     &,term(4),ix(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

! input coordinates are positive
!     if(xi.ge.0) then
        x=xi-int(xi)
!      else
!       x=xi-int(xi)+1
!     endif
!     ix(2)=int(mod(xi,dfloat(nx)))
!     if(ix(2).lt.0) ix(2)=ix(2)+nx
      ix(2)=int(xi)
      ix(1)=ix(2)-1
      if(ix(1).lt.0) then
        ix(1)=ix(1)+nx
        if(rkvec_shift_latt(1).ne.0.d0) isgn(1)=-1
      endif
!     ix(3)=mod(ix(2)+1,nx)
!     ix(4)=mod(ix(2)+2,nx)
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

!     write(6,'(''ix'',4i2)') (ix(j),j=1,4)

      term(1)=-sixth*x*(x-1)*(x-2)*isgn(1)
      term(2)=half*(x+1)*(x-1)*(x-2)
      term(3)=-half*(x+1)*x*(x-2)*isgn(3)
      term(4)=sixth*(x+1)*x*(x-1)*isgn(4)

! Timing tests done with g77-2.95 O(3) on newton (file server)
! Various possible codings are identified with 1), 1a), 2), etc.

! 1) This is the straightforward reference 7.90:
!     do 10 iorb=1,norb
!       orb(iorb)=0
!       ddorb(iorb)=0
!       do 5 k=1,ndim
!   5     dorb(k,iorb)=0
!       do 10 j=1,4
!         orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)
!         ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)
!         do 10 k=1,ndim
!  10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


! 1a) this is 20% faster 6.76
!     do 10 iorb=1,norb
!       orb(iorb)=term(1)*orb_num(iorb,ix(1),iy,iz)
!       ddorb(iorb)=term(1)*ddorb_num(iorb,ix(1),iy,iz)
!       do 5 k=1,ndim
!   5     dorb(k,iorb)=term(1)*dorb_num(k,iorb,ix(1),iy,iz)
!       do 10 j=2,4
!         orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)
!         ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)
!         do 10 k=1,ndim
!  10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


! 2) this is 15% faster 6.99
!     do 5 iorb=1,norb
!       orb(iorb)=0
!       ddorb(iorb)=0
!       do 5 k=1,ndim
!   5     dorb(k,iorb)=0

!     do 10 j=1,4
!       do 10 iorb=1,norb
!         orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)
!         ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)
!         do 10 k=1,ndim
!  10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


! 2a) this is 20% faster 6.43
!     do 5 iorb=1,norb
!       orb(iorb)=term(1)*orb_num(iorb,ix(1),iy,iz)
!       ddorb(iorb)=term(1)*ddorb_num(iorb,ix(1),iy,iz)
!       do 5 k=1,ndim
!   5     dorb(k,iorb)=term(1)*dorb_num(k,iorb,ix(1),iy,iz)

!     do 10 j=2,4
!       do 10 iorb=1,norb
!         orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)
!         ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)
!         do 10 k=1,ndim
!  10       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


! 2b) this is 25% faster 5.68
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


! 3) this is about 20% slower (this routine only) 9.12
!     do 5 iorb=1,norb
!   5   orb(iorb)=0

!     do 6 iorb=1,norb
!   6   ddorb(iorb)=0

!     do 7 iorb=1,norb
!       do 7 k=1,ndim
!   7     dorb(k,iorb)=0

!     do 10 j=1,4
!       do 10 iorb=1,norb
!  10     orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)

!     do 20 j=1,4
!       do 20 iorb=1,norb
!  20     ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)

!     do 30 j=1,4
!       do 30 iorb=1,norb
!         do 30 k=1,ndim
!  30       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


! 4) 10% faster 7.31
!     do 10 iorb=1,norb
!       orb(iorb)=0
!       do 10 j=1,4
!  10     orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)

!     do 20 iorb=1,norb
!       ddorb(iorb)=0
!       do 20 j=1,4
!  20     ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)

!     do 40 iorb=1,norb
!       do 30 k=1,ndim
!  30     dorb(k,iorb)=0
!       do 40 j=1,4
!         do 40 k=1,ndim
!  40       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)


! 4a) this is 20% faster 6.43
!     do 10 iorb=1,norb
!       orb(iorb)=term(1)*orb_num(iorb,ix(1),iy,iz)
!       do 10 j=2,4
!  10     orb(iorb)=orb(iorb)+term(j)*orb_num(iorb,ix(j),iy,iz)

!     do 20 iorb=1,norb
!       ddorb(iorb)=term(1)*ddorb_num(iorb,ix(1),iy,iz)
!       do 20 j=2,4
!  20     ddorb(iorb)=ddorb(iorb)+term(j)*ddorb_num(iorb,ix(j),iy,iz)

!     do 40 iorb=1,norb
!       do 30 k=1,ndim
!  30     dorb(k,iorb)=term(1)*dorb_num(k,iorb,ix(1),iy,iz)
!       do 40 j=2,4
!         do 40 k=1,ndim
!  40       dorb(k,iorb)=dorb(k,iorb)+term(j)*dorb_num(k,iorb,ix(j),iy,iz)

      return
      end
!-----------------------------------------------------------------------

      subroutine interpol_orbe(nx,ny,nz,xi,yi,zi,orb)
! Written by Cyrus Umrigar
      use contr2_mod
      implicit real*8(a-h,o-z)


      dimension orb(*)

      if(abs(inum_orb).eq.1) then
        call interpol1_orbe(nx,ny,nz,xi,yi,zi,orb)
       elseif(abs(inum_orb).eq.4) then
        call interpol4_orbe(nx,ny,nz,xi,yi,zi,orb)
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine interpol1_orbe(nx,ny,nz,xi,yi,zi,orb)
! Written by Cyrus Umrigar
! Evaluate orbitals by looking up closest pt.
! This routine is never used in practice.
! The mesh pts. on which the function values, f, are given, are assumed
! to be at 0,1,2,...nx-1, and similarly for y and z.
! So, pt. nx is the same as pt. 0.
! It is easy to write a routine for general order and unequally-spaced
! grid, but this routine is specialized to make it fast.
! Warning: This routine assumes fn. is periodic and tabulated on entire period.
! If function is tabulated on half-period then this routine needs modification.

      use orbital_grid_mod
      use coefs_mod
      implicit real*8(a-h,o-z)

      dimension orb(*)

!     if(xi.ge.0) then
!       x=xi-int(xi)
!      else
!       x=xi-int(xi)+1
!     endif
      ix=int(mod(xi,dfloat(nx)))
      if(ix.lt.0) ix=ix+nx

!     if(yi.ge.0) then
!       y=yi-int(yi)
!      else
!       y=yi-int(yi)+1
!     endif
      iy=int(mod(yi,dfloat(ny)))
      if(iy.lt.0) iy=iy+ny

!     if(zi.ge.0) then
!       z=zi-int(zi)
!      else
!       z=zi-int(zi)+1
!     endif
      iz=int(mod(zi,dfloat(nz)))
      if(iz.lt.0) iz=iz+nz

      do 10 iorb=1,norb
   10   orb(iorb)=orb_num(iorb,ix,iy,iz)

      return
      end
!-----------------------------------------------------------------------
      subroutine interpol4_orbe(nx,ny,nz,xi,yi,zi,orb)
! Written by Cyrus Umrigar
! Evaluate orbitals by a cubic (4 mesh pts) Lagrange
! interpolation on an equally-spaced 3D grid with periodic boundary conditions.
! The mesh pts. on which the function values, f, are given, are assumed
! to be at 0,1,2,...nx-1, and similarly for y and z.
! So, pt. nx is the same as pt. 0.
! It is easy to write a routine for general order and unequally-spaced
! grid, but this routine is specialized to make it fast.
! If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
! If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      use constants_mod
      use coefs_mod
      use periodic2_mod
      implicit real*8(a-h,o-z)

      dimension orb(*),orb1(norb),term(4),iz(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

! input coordinates are positive
!     if(zi.ge.0) then
        z=zi-int(zi)
!      else
!       z=zi-int(zi)+1
!     endif
!     iz(2)=int(mod(zi,dfloat(nz)))
!     if(iz(2).lt.0) iz(2)=iz(2)+nz
      iz(2)=int(zi)
      iz(1)=iz(2)-1
      if(iz(1).lt.0) then
        iz(1)=iz(1)+nz
        if(rkvec_shift_latt(3).ne.0.d0) isgn(1)=-1
      endif
!     iz(3)=mod(iz(2)+1,nz)
!     iz(4)=mod(iz(2)+2,nz)
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

!     write(6,'(''iz'',4i2)') (iz(j),j=1,4)

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
!-----------------------------------------------------------------------

      subroutine interpol4_ye(nx,ny,xi,yi,iz,orb)
! Written by Cyrus Umrigar
! If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
! If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      use constants_mod
      use coefs_mod
      use periodic2_mod
      implicit real*8(a-h,o-z)

      dimension orb(*),orb1(norb),term(4),iy(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

! input coordinates are positive
!     if(yi.ge.0) then
        y=yi-int(yi)
!      else
!       y=yi-int(yi)+1
!     endif
!     iy(2)=int(mod(yi,dfloat(ny)))
!     if(iy(2).lt.0) iy(2)=iy(2)+ny
      iy(2)=int(yi)
      iy(1)=iy(2)-1
      if(iy(1).lt.0) then
        iy(1)=iy(1)+ny
        if(rkvec_shift_latt(2).ne.0.d0) isgn(1)=-1
      endif
!     iy(3)=mod(iy(2)+1,ny)
!     iy(4)=mod(iy(2)+2,ny)
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

!     write(6,'(''iy'',4i2)') (iy(j),j=1,4)

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
!-----------------------------------------------------------------------

      subroutine interpol4_xe(nx,xi,iy,iz,orb)
! Written by Cyrus Umrigar
! If rkvec_shift_latt(k) =0 function is tabulated on entire period in this dimension
! If rkvec_shift_latt(k)!=0 function is tabulated on half period in this dimension

      use constants_mod
      use orbital_grid_mod
      use coefs_mod
      use periodic2_mod
      implicit real*8(a-h,o-z)

      dimension orb(*) &
     &,term(4),ix(4),isgn(4)

      do 3 i=1,4
    3   isgn(i)=1

! input coordinates are positive
!     if(xi.ge.0) then
        x=xi-int(xi)
!      else
!       x=xi-int(xi)+1
!     endif
!     ix(2)=int(mod(xi,dfloat(nx)))
!     if(ix(2).lt.0) ix(2)=ix(2)+nx
      ix(2)=int(xi)
      ix(1)=ix(2)-1
      if(ix(1).lt.0) then
        ix(1)=ix(1)+nx
        if(rkvec_shift_latt(1).ne.0.d0) isgn(1)=-1
      endif
!     ix(3)=mod(ix(2)+1,nx)
!     ix(4)=mod(ix(2)+2,nx)
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

!     write(6,'(''ix'',4i2)') (ix(j),j=1,4)

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
