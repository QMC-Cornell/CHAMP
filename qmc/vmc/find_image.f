      subroutine check_lattice(rlatt,cutr,isim_cell)
c Written by Cyrus Umrigar
c Checks to see if the lattice vectors specified are the smallest
c ones possible.  This is necessary for the simple heuristic algorithm
c Mathew Foulkes suggested to find the image particle (if any) that lies in the
c inscribing sphere of the nearest Wigner Seitz cell.  Also set cutjas_en/cutjas_ee
c to 1/2 the shortest primitive/simulation cell lattice vector (inscribing sphere radius.
c If the input cutjas is smaller, it will be reset to the smaller value in read_input.
c Also calculate rlenmin to set cutr to 1/2 the shortest lattice vector.
c If it finds one or more shorter lattice vector it stops.  The user should then
c replace in the input the longest lattice vector by the shortest found.  If it finds
c more than one shorter lattice vector the user should do the replacement one at a time
c otherwise one can end up with a linearly dependent set and a zero volume.

      use dim_mod
      use jaspar6_mod
      implicit real*8(a-h,o-z)

!JT      include 'vmc.h'
!JT      include 'force.h'
      parameter (eps=1.d-12)

!JT      common /dim/ ndim
!JT      common /jaspar6/ asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF)
!JT     &,dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF)
!JT     &,asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF)
!JT     &,asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF)
!JT     &,cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF)
!JT     &,cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)

      dimension rlatt(3,3),imax(3)

c nimax is the number of longest lattice vectors
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
        write(6,'(''primitive  cell lattice vector'',i3,'' is longest ; length='',f8.3)')
     &  imax(nimax),rlenmax
        write(6,'(''primitive  cell lattice vector'',i3,'' is shortest; length='',f8.3)')
     &  imin,rlenmin
        cutjas_en=rlenmin/2
       else
        write(6,'(''number of longest simulation cell lattice vectors='',i3)') nimax
        write(6,'(''simulation cell lattice vector'',i3,'' is longest ; length='',f8.3)')
     &  imax(nimax),rlenmax
        write(6,'(''simulation cell lattice vector'',i3,'' is shortest; length='',f8.3)')
     &  imin,rlenmin
        cutjas_ee=rlenmin/2
      endif

      do 40 i1=-1,1
        do 40 i2=-1,1
          do 40 i3=-1,1
            do 40 i=1,nimax
            if((imax(i).eq.1.and.i1.ne.0).or.(imax(i).eq.2.and.i2.ne.0)
     &      .or.(imax(i).eq.3.and.i3.ne.0)) then
              rlen=0
              do 30 k=1,ndim
   30           rlen=rlen+(i1*rlatt(k,1)+i2*rlatt(k,2)+i3*rlatt(k,3))**2
              rlen=sqrt(rlen)
c             write(6,'(''imax(i),i1,i2,i3,rlenmax,rlen'',4i3,9f9.5)') imax(i),i1,i2,i3,rlenmax,rlen
              if(rlen.lt.rlenmax-eps) then
                if(isim_cell.eq.0) then
                  write(6,*) 'Warning: found shorter primitive cell lattice vector'
                 else
                  write(6,*) 'Warning: found shorter simulation cell lattice vector'
                endif
                write(6,'(''i1,i2,i3,rlen='',3i3,f8.3)') i1,i2,i3,rlen
                write(6,'(''new rlatt='',3f8.3)')
     &          i1*rlatt(1,1)+i2*rlatt(1,2)+i3*rlatt(1,3),
     &          i1*rlatt(2,1)+i2*rlatt(2,2)+i3*rlatt(2,3),
     &          i1*rlatt(3,1)+i2*rlatt(3,2)+i3*rlatt(3,3)
c Warning: For the moment to run AF MnO comment out the stop.
c               if(isim_cell.eq.0) then
c                 stop 'found shorter primitive cell lattice vector in check_lattice'
c                else
c                 stop 'found shorter simulation cell lattice vector in check_lattice'
c               endif
              endif
            endif
   40 continue

      return
      end
c-----------------------------------------------------------------------

c     subroutine reduce_sim_cell(r,rlatt_sim,rlatt_sim_inv)
      subroutine reduce_sim_cell(r)
c Written by Cyrus Umrigar
c For any electron position, replace it by the equivalent
c position in the simulation cell centered at the origin.
c r       = position in cartesian coords
c r_basis = position in lattice coords
c rlatt   = lattice vectors
c r       = rlatt * r_basis
c r_basis = rlatt_inv * r
c Note we add eps before taking nint so that when points are on the cell
c boundary (as they will be when removing half simul. lattice vectors from rshift),
c then one consistently gets the positive boundary and not the negative one.
c (Half simul. lattice vectors, because primitive lattice vectors may be linear
c combinations of half. simul. lattice vectors.)
c Warning: I have modified it so it removes only even multiples of simulation
c cell vectors. This could cause problems for e-e-n terms in J.

      use dim_mod
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'ewald.h'
      parameter(eps=1.d-9)

!JT      common /dim/ ndim
      common /periodic/ rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
     &,rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
     &,cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
     &,igvec(3,NGVEC_BIGX),gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX),igmult(NGNORM_BIGX)
     &,igvec_sim(3,NGVEC_SIM_BIGX),gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX),igmult_sim(NGNORM_SIM_BIGX)
     &,rkvec_shift(3),kvec(3,MKPTS),rkvec(3,MKPTS),rknorm(MKPTS)
     &,k_inv(MKPTS),nband(MKPTS),ireal_imag(MORB)
     &,znuc_sum,znuc2_sum,vcell,vcell_sim
     &,ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
     &,ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
     &,ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange

c     dimension r(3),r_basis(3),rlatt_sim(3,3),rlatt_sim_inv(3,3)
      dimension r(3),r_basis(3)

c Find vector in basis coordinates
      do 20 k=1,ndim
        r_basis(k)=0
        do 10 i=1,ndim
   10     r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*r(i)
c  20   r_basis(k)=r_basis(k)-nint(r_basis(k)-eps)
   20   r_basis(k)=r_basis(k)-2*(nint(r_basis(k)-eps)/2)

c     write(6,'(''r_basis'',9f9.4)') r_basis

c Convert back to cartesian coodinates
      do 30 k=1,ndim
        r(k)=0
        do 30 i=1,ndim
   30     r(k)=r(k)+rlatt_sim(k,i)*r_basis(i)

      return
      end
c-----------------------------------------------------------------------

      subroutine find_sim_cell(r,rlatt_inv,r_basis,i_basis)
c Written by Cyrus Umrigar
c For any electron position, find its lattice coordinates
c r       = position in cartesian coords
c r_basis = position in lattice coords
c i_basis = which simulation cell it is in
c rlatt   = lattice vectors
c r       = rlatt * r_basis
c r_basis = rlatt_inv * r

      use dim_mod
      implicit real*8(a-h,o-z)

!JT      common /dim/ ndim
      dimension r(3),r_basis(3),i_basis(3),rlatt_inv(3,3)

c Find vector in basis coordinates
      do 20 k=1,ndim
        r_basis(k)=0
        do 10 i=1,ndim
   10     r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
   20   i_basis(k)=nint(r_basis(k))

      return
      end
c-----------------------------------------------------------------------

      subroutine find_image(r,rlatt,rlatt_inv)
c Written by Cyrus Umrigar
c For any vector (from one particle to another) it finds the
c image that is closest.

      use dim_mod
      implicit real*8(a-h,o-z)

!JT      common /dim/ ndim
      dimension r(3),r_basis(3),rlatt(3,3),rlatt_inv(3,3)
     &,r1_try(3),r2_try(3),r3_try(3),i_sav(3),isign(3)

c Starting from a vector, which is a diff. of 2 vectors, each of which
c have been reduced to the central lattice cell, calculate
c a) its length
c b) sign along each of lattice directions

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

c Check just 8, rather than 27, trapezoids
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

c Replace r by its shortest image
      do 70 i=1,ndim
        do 70 k=1,ndim
   70     r(k)=r(k)-i_sav(i)*rlatt(k,i)

c     write(6,'(''rnew'',9f10.5)') (r(k),k=1,ndim),sqrt(r2)

c debug
c     r2_tmp=0
c     do 80 k=1,ndim
c  80  r2_tmp=r2_tmp+r(k)**2
c     if(r2_tmp.ne.r2) write(6,'(''r2,r2_tmp'',3d12.4)') r2,r2_tmp,r2-r2_tmp

      return
      end
c-----------------------------------------------------------------------

      subroutine find_image2(r,rlatt,r_basis1,r_basis2,i_basis1,i_basis2)
c Written by Cyrus Umrigar
c For any electron positions in lattice coordinates, it finds the
c image that is closest.
c Needs precomputed r_basis1,r_basis2,i_basis1,i_basis2.

      use dim_mod
      implicit real*8(a-h,o-z)

!JT      common /dim/ ndim
      dimension r(3),rlatt(3,3)
     &,r_basis1(3),r_basis2(3),i_basis1(3),i_basis2(3)
     &,r1_try(3),r2_try(3),r3_try(3),i_sav(3),isign(3)

c Find length of original vector and sign along each of lattice directions
      r2=0
      do 20 k=1,ndim
      r2=r2+r(k)**2
   20   isign(k)=int(sign(1.d0,r_basis2(k)-r_basis1(k)-i_basis2(k)+i_basis1(k)))

      do 25 k=1,ndim
   25   i_sav(k)=0

c Check just 8, rather than 27, trapezoids (not needed for orthorhombic lattice)
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

c Replace r by its shortest image
      do 70 i=1,ndim
        do 70 k=1,ndim
   70     r(k)=r(k)-rlatt(k,i)*i_sav(i)

      return
      end
c-----------------------------------------------------------------------

      subroutine find_image3(r,rnorm,rlatt,rlatt_inv)
c Written by Cyrus Umrigar
c For any vector r (from one particle to another) it replaces the vector
c by its closest image and finds its norm

      use dim_mod
      implicit real*8(a-h,o-z)

!JT      common /dim/ ndim
      dimension r(3),r_basis(3),rlatt(3,3),rlatt_inv(3,3)
     &,r1_try(3),r2_try(3),r3_try(3),i_sav(3),isign(3)

c Warning: tempor
c     dimension rsav(3)
c     do 5 k=1,ndim
c   5   rsav(k)=r(k)

c a) reduce vector to central cell by expressing vector in lattice coordinates and
c    removing nint of it in each direction
c b) sign along each of lattice directions of vector reduced to central cell
      do 20 k=1,ndim
        r_basis(k)=0
        do 10 i=1,ndim
   10     r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        r_basis(k)=r_basis(k)-nint(r_basis(k))
   20   isign(k)=nint(sign(1.d0,r_basis(k)))

c Convert back to cartesian coodinates and find squared length
      r2=0
      do 23 k=1,ndim
        r(k)=0
        do 22 i=1,ndim
   22     r(k)=r(k)+rlatt(k,i)*r_basis(i)
   23   r2=r2+r(k)**2

      do 25 k=1,ndim
   25   i_sav(k)=0

c Check just 8, rather than 27, trapezoids (not needed for orthorhombic lattice)
      do 60 i1=0,isign(1),isign(1)
c     do 60 i1=-1,1,1
        do 30 k=1,ndim
   30     r1_try(k)=r(k)-i1*rlatt(k,1)
        do 60 i2=0,isign(2),isign(2)
c       do 60 i2=-1,1,1
          do 40 k=1,ndim
   40       r2_try(k)=r1_try(k)-i2*rlatt(k,2)
          do 60 i3=0,isign(3),isign(3)
c         do 60 i3=-1,1,1
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

c Replace r by its shortest image
      rnorm=0
      do 80 k=1,ndim
        do 70 i=1,ndim
   70     r(k)=r(k)-rlatt(k,i)*i_sav(i)
   80   rnorm=rnorm+r(k)**2
      rnorm=sqrt(rnorm)

c     if(rnorm.gt.5.d0) write(6,'(''long'',6i2,10f8.4)')
c    &(isign(k),k=1,ndim),(i_sav(k),k=1,ndim),rnorm,(r(k),k=1,ndim),(rsav(k),k=1,ndim),(r_basis(k),k=1,ndim)

      return
      end
c-----------------------------------------------------------------------

      subroutine find_image4(rshift,r,rnorm,rlatt,rlatt_inv)
c Written by Cyrus Umrigar
c For any vector r (from one particle to another) it replaces the vector
c by its closest image and finds its norm and the shift needed.
c The shift is modulo simulation lattice vectors.  So if the simulation
c cell is the primitive cell, then rshift is always zero.  The shift is
c used to make sure that two electrons are close to the same nucleus in
c the simulation cell and not just to the same nucleus in the primitive cell.

      use dim_mod
      implicit real*8(a-h,o-z)

!JT      common /dim/ ndim
      dimension r(3),r_basis(3),rshift(3),rlatt(3,3),rlatt_inv(3,3)
     &,r1_try(3),r2_try(3),r3_try(3),i_sav(3),isign(3)

c a) reduce vector to central cell by expressing vector in lattice coordinates and
c    removing nint of it in each direction
c b) sign along each of lattice directions of vector reduced to central cell
c Note: rhift is just a work array here; calculated for real only at end.
      do 20 k=1,ndim
        r_basis(k)=0
        do 10 i=1,ndim
   10     r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        rshift(k)=r_basis(k)-nint(r_basis(k))
   20   isign(k)=nint(sign(1.d0,rshift(k)))

c Convert back to cartesian coodinates and find squared length
      r2=0
      do 23 k=1,ndim
        r(k)=0
        do 22 i=1,ndim
   22     r(k)=r(k)+rlatt(k,i)*rshift(i)
   23   r2=r2+r(k)**2

      do 25 k=1,ndim
   25   i_sav(k)=0

c Check just 8, rather than 27, trapezoids (not needed for orthorhombic lattice)
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

c Replace r by its shortest image and calculate rshift
      rnorm=0
      do 80 k=1,ndim
        rshift(k)=0
        do 70 i=1,ndim
          rshift(k)=rshift(k)+rlatt(k,i)*(nint(r_basis(i))+i_sav(i))
   70     r(k)=r(k)-rlatt(k,i)*i_sav(i)
   80   rnorm=rnorm+r(k)**2
      rnorm=sqrt(rnorm)

c Reduce shift to central simulation cell
c     write(6,'(''rshift_bef'',9f9.4)') rshift
c Warning: calling reduce_sim_cell when rkvec_shift !=0 messes up the
c calculation of the nonlocal pseudopotential.  I do not know why.
      call reduce_sim_cell(rshift)
c     write(6,'(''rshift_aft'',9f9.4)') rshift

      return
      end
