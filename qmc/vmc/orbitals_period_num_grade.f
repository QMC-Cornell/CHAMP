      subroutine orbitals_period_num_grade(x,orb,dorb,ddorb)
c Written by Cyrus Umrigar. Modified by William Parker to interface to Einspline.
c Calculate pw orbitals, gradient and laplacian for electron at x
c by interpolating on a grid.
      use bwfdet_mod
      use bsplines_mod

      use orbital_grid_mod
      implicit real*8(a-h,o-z)

!     include 'vmc.h'
!     include 'force.h'
!     include 'ewald.h'
!     include 'numorb.h'

      common /dim/ ndim
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
c     common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
c    &,ngrid_orbx,ngrid_orby,ngrid_orbz
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
      common /periodic2/ rkvec_shift_latt(3)

      dimension x(3),r_basis(3),orb(*),dorb(3,*),ddorb(*)

      dimension orb_blip_tmp(MORB_OCC,MDET),dorb_blip_tmp(3,MORB_OCC,MDET),
     &     ddorb_blip_tmp(MORB_OCC,MDET)

c Determine position in lattice coordinates
c Find vector in basis coordinates
c Note we add 1 to r_basis because interpol_orb expects a positive input


cwparker This was what was originally here
c     do 20 k=1,ndim
c       r_basis(k)=0
c       do 10 i=1,ndim
c  10     r_basis(k)=r_basis(k)+rlatt_sim_inv(k,i)*x(i)
c  20   r_basis(k)=r_basis(k)-nint(r_basis(k))+1

cwparker Here is what came from the later version of this subroutine
cr_basis in the line after 10 is between  0 and 1 if r_basis in line 10 is > 0
c                            but between -1 and 0 if r_basis in line 10 is < 0
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
cwparker Get the values from Lagrange interpolation routine

      if(inum_orb.eq.4 .or. inum_orb.eq.-4) then
         call interpol_orb(ngrid_orbx,ngrid_orby,ngrid_orbz,xi,yi,zi,
     &                     orb,dorb,ddorb)

cwparker Added from the later version of this subroutine
               do 40 iorb=1,norb
                  orb(iorb)=orb(iorb)*isgn
                  do 30 k=1,ndim
   30                dorb(k,iorb)=dorb(k,iorb)*isgn
   40             ddorb(iorb)=ddorb(iorb)*isgn

      endif

      if(inum_orb.eq.6 .or. inum_orb.eq.-6) then
         if(ndet.gt.1) stop 'can only handle one determinant with blips'
         call bwfdet_main(x,1,1,1,orb_blip_tmp,dorb_blip_tmp,
     &                  ddorb_blip_tmp)
         do iorb=1,norb
              orb(iorb)=orb_blip_tmp(iorb,1)
              dorb(:,iorb)=dorb_blip_tmp(:,iorb,1)
              ddorb(iorb)=ddorb_blip_tmp(iorb,1)
         enddo
      endif

      if(inum_orb.eq.8 .or. inum_orb.eq.-8) then
         call evaluate_bsplines_with_derivatives(x,orb,dorb,ddorb)
      endif

      return
      end
