      subroutine orbitals_period_num(x,orb,dorb,ddorb)
c Written by Cyrus Umrigar
c Calculate pw orbitals, gradient and laplacian
c by interpolating on a grid.
      use bwfdet_mod

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'ewald.h'
      include 'numorb.h'

      common /dim/ ndim
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,ngrid_orbx,ngrid_orby,ngrid_orbz
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
!    &,orb_splines(8,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,
!    &0:MGRID_ORB_PER-1,MORB_OCC)
!    &,grid_orbx(0:MGRID_ORB_PER-1)
!    &,grid_orby(0:MGRID_ORB_PER-1),grid_orbz(0:MGRID_ORB_PER-1)
c    &,orb_splines_explicit(4,4,4,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:
c    &MGRID_ORB_PER-1,MORB_OCC)

      common /periodic2/ rkvec_shift_latt(3)

      dimension x(3,*),r_basis(3)
      dimension orb(MELEC,*),dorb(3,MELEC,*),ddorb(MELEC,*)
      dimension orb_tmp(MORB),dorb_tmp(3,MORB),ddorb_tmp(MORB)
      dimension ict(10),orb_splines_tmp(10),ddorb_splines_tmp(3)

      dimension orb_blip_tmp(MORB_OCC,MDET),dorb_blip_tmp(3,MORB_OCC,MDET),
     &     ddorb_blip_tmp(MORB_OCC,MDET)

      integer ix,iy,iz

cwparker Set up the ict array to tell PSPLINE to return the value
c        of the function, its first derivative in every direction
c        and its second derivative in every direction and the mixed
c        derivatives

      do i=1,10
        ict(i)=1
      enddo

      do iel=1,nelec

c Determine position in lattice coordinates in line 10, same modulo 1 in line 20, and whether phase factor, isgn, is +1 or -1.
c Note that modulo and mod are different functions in fortran.
c Note we add 1 to r_basis because interpol_orb expects a positive input
c The 'if' after 10 does the foll.  Since the k-vectors can be half sim-cell reciprocal lattice vectors
c and we reduce r_basis(k) to be between 0 and 1 in line 20, we need to multiply by -1 for
c each rbasis(k) in line 10 for which floor(r_basis) is odd.
c r_basis in the line just before 20 is between  0 and 1 if r_basis in line 10 is > 0
c                                   but between -1 and 0 if r_basis in line 10 is < 0
c r_basis in line 20 is always between 0 and 1.
        isgn=1
        do 20 k=1,ndim
          r_basis(k)=0
          do 10 i=1,ndim
   10       r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*x(i,iel)
          if(rkvec_shift_latt(k).ne.0.d0) then
            if(r_basis(k).ge.0.d0) then
              isgn=isgn*(-1)**int(r_basis(k))
             else
              isgn=isgn*(-1)**(int(r_basis(k))+1)
            endif
          endif
          r_basis(k)=r_basis(k)-int(r_basis(k))
   20     if(r_basis(k).lt.0.d0) r_basis(k)=r_basis(k)+1

c Warning: I could replace the above with the foll., but since I do not have the time to test it I am
c for the moment just putting in the foll. commented lines:
c Determine position in lattice coordinates in line 10, same modulo 1 in line 20, and whether phase factor, isgn, is +1 or -1.
c Note that modulo and mod are different functions in fortran.
c Since the k-vectors can be half sim-cell reciprocal lattice vectors
c and we reduce r_basis(k) to be between 0 and 1 in line 20, we need to multiply by -1 for 
c each rbasis(k) in line 10 for which floor(r_basis) is odd. 
c       isgn=1
c       do 20 k=1,ndim
c         r_basis(k)=0
c         do 10 i=1,ndim
c  10       r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*x(i,iel)
c         isgn=isgn*(-1)**floor(r_basis(k))
c  20     r_basis(k)=r_basis(k)-floor(r_basis(k))


        if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
          xi=r_basis(1)*ngrid_orbx
          yi=r_basis(2)*ngrid_orby
          zi=r_basis(3)*ngrid_orbz

c       write(6,'(''r_basis'',9f9.4)') r_basis,xi,yi,zi

c Warning:  If we changed the order of iel and iorb then we would not need _tmp arrays
c       call interpol_orb(ngrid_orbx,ngrid_orby,ngrid_orbz,xi,yi,zi,orb(1,iel),dorb(1,1,iel),ddorb(1,iel))
        call interpol_orb(ngrid_orbx,ngrid_orby,ngrid_orbz,xi,yi,zi,orb_tmp,dorb_tmp,ddorb_tmp)

          do iorb=1,norb
            orb(iel,iorb)=orb_tmp(iorb)*isgn
            do k=1,ndim
              dorb(k,iel,iorb)=dorb_tmp(k,iorb)*isgn
            enddo
            ddorb(iel,iorb)=ddorb_tmp(iorb)*isgn
c            if((iorb.eq.1).and.(iel.eq.1)) then
c              write(6,'(''r_basis='',3f9.4)')r_basis(1),r_basis(2),
c     &                                       r_basis(3)
c             write(6,*)'ddorb_lagrange=',ddorb(iel,iorb)
c            endif
          enddo
        endif

!        if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
!          do iorb=1,norb
!
!cwparker for a given point, r_basis(k), calculate the value of the
!c	 splines at the point as well as selected derivatives
!cwparker Compact representation of splines
!             call r8evtricub(r_basis(1),r_basis(2),r_basis(3),
!     &                       grid_orbx,ngrid_orbx,grid_orby,
!     &                       ngrid_orby,grid_orbz,ngrid_orbz,1,1,1,
!     &                       orb_splines(1,0,0,0,iorb),ngrid_orby,
!     &                       ngrid_orbz,ict,orb_splines_tmp,ier)
!
!             if(ier.ne.0) stop 'error in r8evtricub'
!cwparker Explicit splines
!c                     call r8tcspeval(r_basis(1),r_basis(2),r_basis(3),
!c    &                                ict,orb_splines_tmp,grid_orbx,
!c    &                                ngrid_orbx,grid_orby,ngrid_orby,
!c    &                                grid_orbz,ngrid_orbz,1,1,1,
!c    &                                orb_splines_explicit(1,1,1,0,0,0,
!c    &                                                            iorb),
!c    &                                ngrid_orbx,ngrid_orby,ier)
!
!c            if(ier.ne.0) stop 'error in r8tcspeval'
!
!cwparker the first index of the output array, orb_splines_tmp(1), returns
!c        the value of the splines at the input point
!
!             orb_tmp(iorb)    = orb_splines_tmp(1)
!
!cwparker the second, third and fourth indices of the output array,
!c	 orb_splines_tmp(2), orb_splines_tmp(3), etc. return the
!c	 x, y and z derivatives, respectively, of the splines at the
!c	 input point; they must be scaled by the reciprocal simulation
!c        lattice vectors
!
!             dorb_tmp(1,iorb) = rlatt_sim_inv(1,1)*orb_splines_tmp(2)+
!     &                          rlatt_sim_inv(1,2)*orb_splines_tmp(3)+
!     &                          rlatt_sim_inv(1,3)*orb_splines_tmp(4)
!
!             dorb_tmp(2,iorb) = rlatt_sim_inv(2,1)*orb_splines_tmp(2)+
!     &                          rlatt_sim_inv(2,2)*orb_splines_tmp(3)+
!     &                          rlatt_sim_inv(2,3)*orb_splines_tmp(4)
!
!             dorb_tmp(3,iorb) = rlatt_sim_inv(3,1)*orb_splines_tmp(2)+
!     &                          rlatt_sim_inv(3,2)*orb_splines_tmp(3)+
!     &                          rlatt_sim_inv(3,3)*orb_splines_tmp(4)
!
!
!cwparker the fifth, sixth and seventh indices of the output array,
!c	 orb_splines_tmp(5), etc. return the second x, y and z
!c	 derivatives, respectively of the splines at the input point
!c	 the eighth, ninth and tenth indices of the output array
!c	 return the mixed derivatives, d2f/dxdy, d2f/dxdz, d2f/dydz
!c	 respectively, that we need to transform to Cartesian coordinates
!
!             ddorb_splines_tmp(1)=rlatt_sim_inv(1,1)**2*
!     &                      orb_splines_tmp(5)+
!     &                      2*rlatt_sim_inv(1,1)*rlatt_sim_inv(1,2)*
!     &                      orb_splines_tmp(8)+rlatt_sim_inv(1,2)**2*
!     &                      orb_splines_tmp(6)+2*rlatt_sim_inv(1,2)*
!     &                      rlatt_sim_inv(1,3)*orb_splines_tmp(10)+
!     &                      rlatt_sim_inv(1,3)**2*orb_splines_tmp(7)+
!     &                      2*rlatt_sim_inv(1,3)*rlatt_sim_inv(1,1)*
!     &                      orb_splines_tmp(9)
!
!             ddorb_splines_tmp(2)=rlatt_sim_inv(2,1)**2*
!     &                      orb_splines_tmp(5)+
!     &                      2*rlatt_sim_inv(2,1)*rlatt_sim_inv(2,2)*
!     &                      orb_splines_tmp(8)+rlatt_sim_inv(2,2)**2*
!     &                      orb_splines_tmp(6)+2*rlatt_sim_inv(2,2)*
!     &                      rlatt_sim_inv(2,3)*orb_splines_tmp(10)+
!     &                      rlatt_sim_inv(2,3)**2*orb_splines_tmp(7)+
!     &                      2*rlatt_sim_inv(2,3)*rlatt_sim_inv(2,1)*
!     &                      orb_splines_tmp(9)
!
!             ddorb_splines_tmp(3)=rlatt_sim_inv(3,1)**2*
!     &                      orb_splines_tmp(5)+
!     &                      2*rlatt_sim_inv(3,1)*rlatt_sim_inv(3,2)*
!     &                      orb_splines_tmp(8)+rlatt_sim_inv(3,2)**2*
!     &                      orb_splines_tmp(6)+2*rlatt_sim_inv(3,2)*
!     &                      rlatt_sim_inv(3,3)*orb_splines_tmp(10)+
!     &                      rlatt_sim_inv(3,3)**2*orb_splines_tmp(7)+
!     &                      2*rlatt_sim_inv(3,3)*rlatt_sim_inv(3,1)*
!     &                      orb_splines_tmp(9)
!
!             ddorb_tmp(iorb) =  ddorb_splines_tmp(1)+
!     &                          ddorb_splines_tmp(2)+
!     &                          ddorb_splines_tmp(3)
!
!             orb(iel,iorb)=orb_tmp(iorb)*isgn
!             do k=1,ndim
!                dorb(k,iel,iorb)=dorb_tmp(k,iorb)*isgn
!             enddo
!             ddorb(iel,iorb)=ddorb_tmp(iorb)*isgn
!
!cwparker     if((iorb.eq.1).and.(iel.eq.1)) then
!c             write(6,'(''r_basis='',3f9.4)')r_basis(1),r_basis(2),
!c    &                                       r_basis(3)
!c             write(6,*)'ddorb_splines=',ddorb(iel,iorb)
!c            endif
!
!
!cwparker End of do loop over the orbitals
!          enddo
!
!cwparker End of if statement for splines
!        endif

cwparker If statment for blips
        if(inum_orb.eq.6 .or. inum_orb.eq.-6) then

         if(ndet.gt.1) stop 'blip implementation can only handle one det
     &erminant'

cwparker Use x since bwfdet_main converts to crystal lattice units for us
c	 Pass 1 for iw because we need the wavefunction here
c	 Pass 1 for igl because we need the gradient and Laplacian here
c	 Pass 1 for spin because we don't differentiate between spins yet

cwparker Loop over orbitals and set temporary arrays equal to permanent ones

             call bwfdet_main(x(:,iel),1,1,1,orb_blip_tmp,dorb_blip_tmp
     &                   ,ddorb_blip_tmp)

          do iorb=1,norb

             orb(iel,iorb)=orb_blip_tmp(iorb,1)
             dorb(:,iel,iorb)=dorb_blip_tmp(:,iorb,1)
             ddorb(iel,iorb)=ddorb_blip_tmp(iorb,1)

cwparker End of loop over orbitals
          enddo

cwparker End of if statement for blips
        endif

cwparker End of loop over the electrons
      enddo

cwp   write(6,*)'End of orbitals_period_num'

      return
      end
c-----------------------------------------------------------------------
      subroutine orbitals_period_nume(x,orb)
c Written by Cyrus Umrigar
c Calculate pw orbitals for electron iel by interpolating on a grid.
      use bwfdet_mod

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'ewald.h'
      include 'numorb.h'

      common /dim/ ndim
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
     &,ngrid_orbx,ngrid_orby,ngrid_orbz
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
!    &,orb_splines(8,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,
!    &0:MGRID_ORB_PER-1,MORB_OCC)
!    &,grid_orbx(0:MGRID_ORB_PER-1)
!    &,grid_orby(0:MGRID_ORB_PER-1),grid_orbz(0:MGRID_ORB_PER-1)
c    &,orb_splines_explicit(4,4,4,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:
c    &MGRID_ORB_PER-1,MORB_OCC)
      common /periodic2/ rkvec_shift_latt(3)

      dimension x(3),orb(*)
      dimension r_basis(3),ict(10),orb_splines_tmp(10)

      dimension orb_blip_tmp(MORB_OCC,MDET),dorb_blip_tmp(3,MORB_OCC,MDET),
     &     ddorb_blip_tmp(MORB_OCC,MDET)

c Determine position in lattice coordinates
c Find vector in basis coordinates
c Note we add 1 to r_basis because interpol_orb expects a positive input
c r_basis in the line after 10 is between  0 and 1 if r_basis in line 10 is > 0
c                             but between -1 and 0 if r_basis in line 10 is < 0
c r_basis in line 20 is always between 0 and 1
        isgn=1
        do 20 k=1,ndim
          r_basis(k)=0
          do 10 i=1,ndim
   10       r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*x(i)
          if(rkvec_shift_latt(k).ne.0.d0) then
            if(r_basis(k).ge.0.d0) then
              isgn=isgn*(-1)**int(r_basis(k))
             else
              isgn=isgn*(-1)**(int(r_basis(k))+1)
            endif
          endif
          r_basis(k)=r_basis(k)-int(r_basis(k))
   20     if(r_basis(k).lt.0.d0) r_basis(k)=r_basis(k)+1

      xi=r_basis(1)*ngrid_orbx
      yi=r_basis(2)*ngrid_orby
      zi=r_basis(3)*ngrid_orbz

c     write(6,'(''r_basis'',9f9.4)') r_basis,xi,yi,zi

      if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
         call interpol_orbe(ngrid_orbx,ngrid_orby,ngrid_orbz,xi,yi,zi,orb)

         do 30 iorb=1,norb
   30       orb(iorb)=orb(iorb)*isgn

      endif

!      if(inum_orb.eq.5 .or. inum_orb.eq.-5) then
!
!         ict(1)=1
!
!         do i=2,10
!            ict(i)=0
!         enddo
!
!c        ngrid_orbx_plus1=ngrid_orbx+1
!c        ngrid_orby_plus1=ngrid_orby+1
!c        ngrid_orbz_plus1=ngrid_orbz+1
!
!         do iorb=1,norb
!
!cwparker Compact representation of splines
!             call r8evtricub(r_basis(1),r_basis(2),r_basis(3),
!     &                       grid_orbx,ngrid_orbx,grid_orby,
!     &                       ngrid_orby,grid_orbz,ngrid_orbz
!     &                       ,1,1,1,orb_splines(1,0,0,0,iorb),
!     &                       ngrid_orby,ngrid_orbz,
!     &                       ict,orb_splines_tmp,ier)
!
!             if(ier.ne.0) stop 'error in r8evtricub'
!cwparker Explicit representation of splines
!c                     call r8tcspeval(r_basis(1),r_basis(2),r_basis(3),
!c    &                                ict,orb_splines_tmp,grid_orbx,
!c    &                                ngrid_orbx,grid_orby,ngrid_orby,
!c    &                                grid_orbz,ngrid_orbz,1,1,1,
!c    &                                orb_splines_explicit(1,1,1,0,0,0,
!c    &                                                            iorb),
!c    &                                ngrid_orbx,ngrid_orby,ier)
!c
!c            if(ier.ne.0) stop 'error in r8tcspeval'
!
!
!             orb(iorb)=orb_splines_tmp(1)*isgn
!
!        enddo
!
!      endif

      if(inum_orb.eq.6 .or. inum_orb.eq.-6) then

         if(ndet.gt.1) stop 'can only handle one determinant'

         call bwfdet_main(x,1,0,1,orb_blip_tmp,dorb_blip_tmp,
     &                    ddorb_blip_tmp)

         do iorb=1,norb
            orb(iorb)=orb_blip_tmp(iorb,1)
         enddo

      endif

      return
      end
