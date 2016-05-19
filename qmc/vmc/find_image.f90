      subroutine check_lattice(rlatt,cutr,isim_cell)
! Written by Cyrus Umrigar
! Checks to see if the lattice vectors specified are the smallest
! ones possible.  This is necessary for the simple heuristic algorithm
! Mathew Foulkes suggested to find the image particle (if any) that lies in the
! inscribing sphere of the nearest Wigner Seitz cell.  Also set cutjas_en/cutjas_ee
! to 1/2 the shortest primitive/simulation cell lattice vector (inscribing sphere radius.
! If the input cutjas is smaller, it will be reset to the smaller value in read_input.
! Also calculate rlenmin to set cutr to 1/2 the shortest lattice vector.
! If it finds one or more shorter lattice vector it stops.  The user should then
! replace in the input the longest lattice vector by the shortest found.  If it finds
! more than one shorter lattice vector the user should do the replacement one at a time
! otherwise one can end up with a linearly dependent set and a zero volume.

      use dim_mod
      use jaspar6_mod
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-12)


      dimension rlatt(3,3),imax(3)

! nimax is the number of longest lattice vectors
      nimax=0
      rlenmax=0
      rlenmin=9.d99
      do 20 i=1,ndim
        rlen=0
        do 10 k=1,ndim
   10     rlen=rlen+rlatt(k,i)**2
        if(abs(rlen-rlenmax).lt.eps) then
          rlenmax=max(rlen,rlenmax)
          nimax=nimax+1
          imax(nimax)=i
         elseif(rlen.gt.rlenmax-eps) then
          rlenmax=max(rlen,rlenmax)
          nimax=1
          imax(nimax)=i
        endif
        if(rlen.lt.rlenmin) then
          rlenmin=min(rlen,rlenmin)
          imin=i
        endif
   20 continue
      rlenmax=sqrt(rlenmax)
      rlenmin=sqrt(rlenmin)
      cutr=rlenmin/2

      if(isim_cell.eq.0) then
        write(6,'(''number of longest primitive cell lattice vectors='',i3)') nimax
        write(6,'(''primitive  cell lattice vector'',i3,'' is longest ; length='',f8.3)') &
     &  imax(nimax),rlenmax
        write(6,'(''primitive  cell lattice vector'',i3,'' is shortest; length='',f8.3)') &
     &  imin,rlenmin
        cutjas_en=rlenmin/2
       else
        write(6,'(''number of longest simulation cell lattice vectors='',i3)') nimax
        write(6,'(''simulation cell lattice vector'',i3,'' is longest ; length='',f8.3)') &
     &  imax(nimax),rlenmax
        write(6,'(''simulation cell lattice vector'',i3,'' is shortest; length='',f8.3)') &
     &  imin,rlenmin
        cutjas_ee=rlenmin/2
      endif

      do 40 i1=-1,1
        do 40 i2=-1,1
          do 40 i3=-1,1
            do 40 i=1,nimax
            if((imax(i).eq.1.and.i1.ne.0).or.(imax(i).eq.2.and.i2.ne.0) &
     &      .or.(imax(i).eq.3.and.i3.ne.0)) then
              rlen=0
              do 30 k=1,ndim
   30           rlen=rlen+(i1*rlatt(k,1)+i2*rlatt(k,2)+i3*rlatt(k,3))**2
              rlen=sqrt(rlen)
!             write(6,'(''imax(i),i1,i2,i3,rlenmax,rlen'',4i3,9f9.5)') imax(i),i1,i2,i3,rlenmax,rlen
              if(rlen.lt.rlenmax-eps) then
                if(isim_cell.eq.0) then
                  write(6,*) 'Warning: found shorter primitive cell lattice vector'
                 else
                  write(6,*) 'Warning: found shorter simulation cell lattice vector'
                endif
                write(6,'(''i1,i2,i3,rlen='',3i3,f8.3)') i1,i2,i3,rlen
                write(6,'(''new rlatt='',3f8.3)') &
     &          i1*rlatt(1,1)+i2*rlatt(1,2)+i3*rlatt(1,3), &
     &          i1*rlatt(2,1)+i2*rlatt(2,2)+i3*rlatt(2,3), &
     &          i1*rlatt(3,1)+i2*rlatt(3,2)+i3*rlatt(3,3)
! Warning: For the moment to run AF MnO comment out the stop.
!               if(isim_cell.eq.0) then
!                 stop 'found shorter primitive cell lattice vector in check_lattice'
!                else
!                 stop 'found shorter simulation cell lattice vector in check_lattice'
!               endif
              endif
            endif
   40 continue

      return
      end
!-----------------------------------------------------------------------

!     subroutine reduce_sim_cell(r,rlatt_sim,rlatt_sim_inv)
      subroutine reduce_sim_cell(r)
! Written by Cyrus Umrigar
! For any electron position, replace it by the equivalent
! position in the simulation cell centered at the origin.
! r       = position in cartesian coords
! r_basis = position in lattice coords
! rlatt   = lattice vectors
! r       = rlatt * r_basis
! r_basis = rlatt_inv * r
! Note we add eps before taking nint so that when points are on the cell
! boundary (as they will be when removing half simul. lattice vectors from rshift),
! then one consistently gets the positive boundary and not the negative one.
! (Half simul. lattice vectors, because primitive lattice vectors may be linear
! combinations of half. simul. lattice vectors.)
! Warning: I have modified it so it removes only even multiples of simulation
! cell vectors. This could cause problems for e-e-n terms in J.

      use dim_mod
      use periodic_mod
      use contrl_per_mod
      use periodic_1d_mod
      implicit real*8(a-h,o-z)
      parameter(eps=1.d-9)


!     dimension r(3),r_basis(3),rlatt_sim(3,3),rlatt_sim_inv(3,3)
      dimension r(3),r_basis(3)

! ACM, June 1, 2010: add case for systems periodic in only 1D so DMC will work
      if (iperiodic.eq.1) then
         r(1) = modulo(r(1),alattice) ! returns a number between 0 and alattice
         if (r(1).ge.(alattice/2.)) r(1) = r(1) - alattice
         return
      endif

! Find vector in basis coordinates
      do 20 k=1,ndim
        r_basis(k)=0
        do 10 i=1,ndim
   10     r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*r(i)
!  20   r_basis(k)=r_basis(k)-nint(r_basis(k)-eps)
   20   r_basis(k)=r_basis(k)-2*(nint(r_basis(k)-eps)/2)

!     write(6,'(''r_basis'',9f9.4)') r_basis

! Convert back to cartesian coodinates
      do 30 k=1,ndim
        r(k)=0
        do 30 i=1,ndim
   30     r(k)=r(k)+rlatt_sim(k,i)*r_basis(i)

      return
      end
!-----------------------------------------------------------------------

      subroutine find_sim_cell(r,rlatt_inv,r_basis,i_basis)
! Written by Cyrus Umrigar
! For any electron position, find its lattice coordinates
! r       = position in cartesian coords
! r_basis = position in lattice coords
! i_basis = which simulation cell it is in
! rlatt   = lattice vectors
! r       = rlatt * r_basis
! r_basis = rlatt_inv * r

      use dim_mod
      implicit real*8(a-h,o-z)

      dimension r(3),r_basis(3),i_basis(3),rlatt_inv(3,3)

! Find vector in basis coordinates
      do 20 k=1,ndim
        r_basis(k)=0
        do 10 i=1,ndim
   10     r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
   20   i_basis(k)=nint(r_basis(k))

      return
      end
!-----------------------------------------------------------------------

      subroutine find_image(r,rlatt,rlatt_inv)
! Written by Cyrus Umrigar
! For any vector (from one particle to another) it finds the
! image that is closest.

      use dim_mod
      implicit real*8(a-h,o-z)

      dimension r(3),r_basis(3),rlatt(3,3),rlatt_inv(3,3) &
     &,r1_try(3),r2_try(3),r3_try(3),i_sav(3),isign(3)

! Starting from a vector, which is a diff. of 2 vectors, each of which
! have been reduced to the central lattice cell, calculate
! a) its length
! b) sign along each of lattice directions

      r2=0
      do 20 k=1,ndim
        r2=r2+r(k)**2
        r_basis(k)=0
        do 10 i=1,ndim
   10     r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        if(abs(r_basis(k)).gt.1.d0) write(6,'(''**Warning, abs(r_basis)>1'')')
   20   isign(k)=nint(sign(1.d0,r_basis(k)))

      do 25 k=1,ndim
   25   i_sav(k)=0

! Check just 8, rather than 27, trapezoids
      do 60 i1=0,isign(1),isign(1)
        do 30 k=1,ndim
   30     r1_try(k)=r(k)-i1*rlatt(k,1)
        do 60 i2=0,isign(2),isign(2)
          do 40 k=1,ndim
   40       r2_try(k)=r1_try(k)-i2*rlatt(k,2)
          do 60 i3=0,isign(3),isign(3)
            r_try2=0
            do 50 k=1,ndim
              r3_try(k)=r2_try(k)-i3*rlatt(k,3)
   50         r_try2=r_try2+r3_try(k)**2
          if(r_try2.lt.r2) then
            i_sav(1)=i1
            i_sav(2)=i2
            i_sav(3)=i3
            r2=r_try2
          endif
   60 continue

! Replace r by its shortest image
      do 70 i=1,ndim
        do 70 k=1,ndim
   70     r(k)=r(k)-i_sav(i)*rlatt(k,i)

!     write(6,'(''rnew'',9f10.5)') (r(k),k=1,ndim),sqrt(r2)

! debug
!     r2_tmp=0
!     do 80 k=1,ndim
!  80  r2_tmp=r2_tmp+r(k)**2
!     if(r2_tmp.ne.r2) write(6,'(''r2,r2_tmp'',3d12.4)') r2,r2_tmp,r2-r2_tmp

      return
      end
!-----------------------------------------------------------------------

      subroutine find_image2(r,rlatt,r_basis1,r_basis2,i_basis1,i_basis2)
! Written by Cyrus Umrigar
! For any electron positions in lattice coordinates, it finds the
! image that is closest.
! Needs precomputed r_basis1,r_basis2,i_basis1,i_basis2.

      use dim_mod
      implicit real*8(a-h,o-z)

      dimension r(3),rlatt(3,3) &
     &,r_basis1(3),r_basis2(3),i_basis1(3),i_basis2(3) &
     &,r1_try(3),r2_try(3),r3_try(3),i_sav(3),isign(3)

! Find length of original vector and sign along each of lattice directions
      r2=0
      do 20 k=1,ndim
      r2=r2+r(k)**2
   20   isign(k)=int(sign(1.d0,r_basis2(k)-r_basis1(k)-i_basis2(k)+i_basis1(k)))

      do 25 k=1,ndim
   25   i_sav(k)=0

! Check just 8, rather than 27, trapezoids (not needed for orthorhombic lattice)
      do 60 i1=0,isign(1),isign(1)
        do 30 k=1,ndim
   30     r1_try(k)=r(k)-rlatt(k,1)*(i1+i_basis2(1)-i_basis1(1))
        do 60 i2=0,isign(2),isign(2)
          do 40 k=1,ndim
   40       r2_try(k)=r1_try(k)-rlatt(k,2)*(i2+i_basis2(2)-i_basis1(2))
          do 60 i3=0,isign(3),isign(3)
            r_try2=0
            do 50 k=1,ndim
              r3_try(k)=r2_try(k)-rlatt(k,3)*(i3+i_basis2(3)-i_basis1(3))
   50         r_try2=r_try2+r3_try(k)**2
          if(r_try2.lt.r2) then
            i_sav(1)=i1+i_basis2(1)-i_basis1(1)
            i_sav(2)=i2+i_basis2(2)-i_basis1(2)
            i_sav(3)=i3+i_basis2(3)-i_basis1(3)
            r2=r_try2
          endif
   60 continue

! Replace r by its shortest image
      do 70 i=1,ndim
        do 70 k=1,ndim
   70     r(k)=r(k)-rlatt(k,i)*i_sav(i)

      return
      end
!-----------------------------------------------------------------------

      subroutine find_image3(r,rnorm,rlatt,rlatt_inv)
! Written by Cyrus Umrigar
! For any vector r (from one particle to another) it replaces the vector
! by its closest image and finds its norm

      use dim_mod
      implicit real*8(a-h,o-z)

      dimension r(3),r_basis(3),rlatt(3,3),rlatt_inv(3,3) &
     &,r1_try(3),r2_try(3),r3_try(3),i_sav(3),isign(3)

! Warning: tempor
!     dimension rsav(3)
!     do 5 k=1,ndim
!   5   rsav(k)=r(k)

! a) reduce vector to central cell by expressing vector in lattice coordinates and
!    removing nint of it in each direction
! b) sign along each of lattice directions of vector reduced to central cell
      do 20 k=1,ndim
        r_basis(k)=0
        do 10 i=1,ndim
   10     r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        r_basis(k)=r_basis(k)-nint(r_basis(k))
   20   isign(k)=nint(sign(1.d0,r_basis(k)))

! Convert back to cartesian coodinates and find squared length
      r2=0
      do 23 k=1,ndim
        r(k)=0
        do 22 i=1,ndim
   22     r(k)=r(k)+rlatt(k,i)*r_basis(i)
   23   r2=r2+r(k)**2

      do 25 k=1,ndim
   25   i_sav(k)=0

! Check just 8, rather than 27, trapezoids (not needed for orthorhombic lattice)
      do 60 i1=0,isign(1),isign(1)
!     do 60 i1=-1,1,1
        do 30 k=1,ndim
   30     r1_try(k)=r(k)-i1*rlatt(k,1)
        do 60 i2=0,isign(2),isign(2)
!       do 60 i2=-1,1,1
          do 40 k=1,ndim
   40       r2_try(k)=r1_try(k)-i2*rlatt(k,2)
          do 60 i3=0,isign(3),isign(3)
!         do 60 i3=-1,1,1
            r_try2=0
            do 50 k=1,ndim
              r3_try(k)=r2_try(k)-i3*rlatt(k,3)
   50         r_try2=r_try2+r3_try(k)**2
          if(r_try2.lt.r2) then
            i_sav(1)=i1
            i_sav(2)=i2
            i_sav(3)=i3
            r2=r_try2
          endif
   60 continue

! Replace r by its shortest image
      rnorm=0
      do 80 k=1,ndim
        do 70 i=1,ndim
   70     r(k)=r(k)-rlatt(k,i)*i_sav(i)
   80   rnorm=rnorm+r(k)**2
      rnorm=sqrt(rnorm)

!     if(rnorm.gt.5.d0) write(6,'(''long'',6i2,10f8.4)')
!    &(isign(k),k=1,ndim),(i_sav(k),k=1,ndim),rnorm,(r(k),k=1,ndim),(rsav(k),k=1,ndim),(r_basis(k),k=1,ndim)

      return
      end
!-----------------------------------------------------------------------

      subroutine find_image4(rshift,r,rnorm,rlatt,rlatt_inv)
! Written by Cyrus Umrigar
! For any vector r (from one particle to another) it replaces the vector
! by its closest image and finds its norm and the shift needed.
! The shift is modulo simulation lattice vectors.  So if the simulation
! cell is the primitive cell, then rshift is always zero.  The shift is
! used to make sure that two electrons are close to the same nucleus in
! the simulation cell and not just to the same nucleus in the primitive cell.

      use dim_mod
      implicit real*8(a-h,o-z)

      dimension r(3),r_basis(3),rshift(3),rlatt(3,3),rlatt_inv(3,3) &
     &,r1_try(3),r2_try(3),r3_try(3),i_sav(3),isign(3)

! a) reduce vector to central cell by expressing vector in lattice coordinates and
!    removing nint of it in each direction
! b) sign along each of lattice directions of vector reduced to central cell
! Note: rhift is just a work array here; calculated for real only at end.
      do 20 k=1,ndim
        r_basis(k)=0
        do 10 i=1,ndim
   10     r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        rshift(k)=r_basis(k)-nint(r_basis(k))
   20   isign(k)=nint(sign(1.d0,rshift(k)))

! Convert back to cartesian coodinates and find squared length
      r2=0
      do 23 k=1,ndim
        r(k)=0
        do 22 i=1,ndim
   22     r(k)=r(k)+rlatt(k,i)*rshift(i)
   23   r2=r2+r(k)**2

      do 25 k=1,ndim
   25   i_sav(k)=0

! Check just 8, rather than 27, trapezoids (not needed for orthorhombic lattice)
      do 60 i1=0,isign(1),isign(1)
        do 30 k=1,ndim
   30     r1_try(k)=r(k)-i1*rlatt(k,1)
        do 60 i2=0,isign(2),isign(2)
          do 40 k=1,ndim
   40       r2_try(k)=r1_try(k)-i2*rlatt(k,2)
          do 60 i3=0,isign(3),isign(3)
            r_try2=0
            do 50 k=1,ndim
              r3_try(k)=r2_try(k)-i3*rlatt(k,3)
   50         r_try2=r_try2+r3_try(k)**2
          if(r_try2.lt.r2) then
            i_sav(1)=i1
            i_sav(2)=i2
            i_sav(3)=i3
            r2=r_try2
          endif
   60 continue

! Replace r by its shortest image and calculate rshift
      rnorm=0
      do 80 k=1,ndim
        rshift(k)=0
        do 70 i=1,ndim
          rshift(k)=rshift(k)+rlatt(k,i)*(nint(r_basis(i))+i_sav(i))
   70     r(k)=r(k)-rlatt(k,i)*i_sav(i)
   80   rnorm=rnorm+r(k)**2
      rnorm=sqrt(rnorm)

! Reduce shift to central simulation cell
!     write(6,'(''rshift_bef'',9f9.4)') rshift
! Warning: calling reduce_sim_cell when rkvec_shift !=0 messes up the
! calculation of the nonlocal pseudopotential.  I do not know why.
      call reduce_sim_cell(rshift)
!     write(6,'(''rshift_aft'',9f9.4)') rshift

      return
      end


!-----------------------------------------------------------------------

      subroutine find_image_1d(r,rnorm)
! Written by Abhijit Mehta
!  (Based on find_image4)
! For a system with 1D periodic BC's, this takes a vector r (from one particle
!   to another) in and replaces the vector by its closest image.
! It also calculates the norm and returns that in rnorm.
! The shift is modulo the simulation lattice vector 'alattice'
! Our convention is that -alattice/2 <= r < alattice/2
! Since the simulation cell is the primitive cell for our 1D periodic systems,
! we don't need to calculate an 'rshift' since it is always zero.

!  We previously did the modulo math explicitly in distances.f and ewald.f
!  Even though this is very little code, we make it a separate subroutine
!   for two reasons:
!  - This way, the structure for 1D periodic BC's parallels that for 3D periodic
!  - If we decide later that we want r to run from 0 to alattice instead of
!      from -alattice/2 to alattice/2, it's easier to just make the change here

      use dim_mod
      use periodic_1d_mod
      implicit real*8(a-h,o-z)

      dimension r(3)

!    Note that "modulo(a,b)" is a Fortran 90 function which returns
!     a mod b, so that the sign of the answer always matches the sign of b
!    This is different from the function "MOD(a,b)", which returns an
!     answer that has the same sign as a
      r(1) = modulo(r(1),alattice) ! returns a number between 0 and alattice
      if (r(1).ge.(alattice/2.)) r(1) = r(1) - alattice
      rnorm = 0.
      do k = 1,ndim
         rnorm = rnorm + r(k)**2
      enddo
      rnorm = dsqrt(rnorm)

      return
      end
