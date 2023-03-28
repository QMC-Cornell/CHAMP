      subroutine nonloc(x,rshift,rvec_en,r_en)
! This routine is called by vmc and dmc. l_do_tmoves=true only if doing dmc and tmoves=true in input.
! Warning: I may still need to update the new tmove parts for the derivatives put in by Julien
! Warning: At present it uses the same points for the reverse T-move to save time, but this may create a tiny bias.
! Warning: Assumes for the reverse T-move probability that the nonlocal ranges of the atoms are non overlapping.
! npotd(ict) is the number of channels, local and nonlocal, for nucleus type ict. In most of the code npotd=npotu and only npotd is used.
! lpotp1=L+1, where L is the local channel.  It would have been more clear to start some arrays at 0, in which case the +1 would not be needed.
! exp(value) is the ratio of the Jastrows at the quadrature points and the initial point.
! w_psi_quad is the wt times the ratio of the wavefns at the quadrature points and the initial point
! Written by Claudia Filippi, modified by Cyrus Umrigar

      use constants_mod
      use control_mod
      use deriv_orb_mod
      use periodic_jastrow_mod !WAS
      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use pseudo_mod
      use contrl_per_mod
      use periodic_mod
      use qua_mod
      use contrldmc_mod, only : tau,taunow,tmoves
      use config_dmc_mod, only : xoldw, voldw, psidow, psijow
      use average_mod, only : current_walker
      use variables_mod, only: l_mode_dmc
      use pjase_mod, only: ido_pjas
      use gamma_mod    !TA
      use orbe_mod     !TA
      use orbitals_mod !TA
! Temporary:
      use div_v_dmc_mod
      use delocc_mod

      implicit real*8(a-h,o-z)
      integer, parameter :: max_lpotp1=5

      dimension x(3,*),rshift(3,nelec,*),rvec_en(3,nelec,*),r_en(nelec,*)
      dimension rr_en(nelec,ncent),rr_en2(nelec,ncent),rr_en_sav(ncent),rr_en2_sav(ncent) &
     &,xsav(3),rshift_sav(3,ncent),rvec_en_sav(3,ncent),r_en_sav(ncent)
      dimension xnew(3)               ! local variables for tmoves
      dimension x_tmove_sav(3,nquad*ncent), gpsp_tmove(nquad*ncent), gpsp_tmove_heatbath(nquad*ncent), w_psi_quad(nquad*ncent), pl_costh(max_lpotp1,nquad,nquad)
      integer iwhich_cent_tmove(nquad*ncent)

      if(maxval(npotd).gt.5) stop 'maxval(npotd).gt.max_lpotp1'
      do iq=1,nquad
        do jq=1,iq
          costh=min(1.d0,max(-1.d0,xq(iq)*xq(jq)+yq(iq)*yq(jq)+zq(iq)*zq(jq)))
          do l=1,maxval(npotd)
            pl_costh(l,iq,jq)=yl0(l,costh)
            pl_costh(l,jq,iq)=pl_costh(l,iq,jq)
          enddo
        enddo
      enddo

      do 10 i=1,nelec
        do 10 ic=1,ncent
          call scale_dist(r_en(i,ic),rr_en(i,ic),1)
   10     call scale_dist(r_en(i,ic),rr_en2(i,ic),3)

      psp_nonloc_orb(:,:) = 0d0 !TA

      do 150 i=1,nelec

! Save position of ith electron and its distances etc. from all nuclei
        do 11 k=1,ndim
   11     xsav(k)=x(k,i)
        do 12 jc=1,ncent
          r_en_sav(jc)=r_en(i,jc)
          rr_en_sav(jc)=rr_en(i,jc)
          rr_en2_sav(jc)=rr_en2(i,jc)
          do 12 k=1,ndim
            rshift_sav(k,jc)=rshift(k,i,jc)
   12       rvec_en_sav(k,jc)=rvec_en(k,i,jc)

! Note if nonlocal parts of psp. are overlapping, ntmove_pts can be multiple of nquad
        ntmove_pts=0
        do 100 ic=1,ncent
          ict=iwctype(ic)

! vps was calculated by calling getvps_xx
          iskip=1
          do 15 l=1,npotd(ict)
   15       if(l.ne.lpotp1(ict) .and. dabs(vps(i,ic,l)).gt.1.d-4) iskip=0

          if(iskip.eq.0) then ! If it is within the radius of any of the nonlocal channels

            if(l_mode_dmc) call rotqua ! If we are doing t-moves then the grid needs to be rotated more often to prevent electrons from being along rays

            ri=one/r_en(i,ic)

            do 60 iq=1,nquad
              ntmove_pts=ntmove_pts+1
              iwhich_cent_tmove(ntmove_pts)=ic
              costh=rvec_en_sav(1,ic)*xq(iq)+rvec_en_sav(2,ic)*yq(iq)+rvec_en_sav(3,ic)*zq(iq)
              costh=costh*ri

              if (l_mode_dmc) then !TA
                quadx(1,ntmove_pts,i,current_walker) = xq(iq)
                quadx(2,ntmove_pts,i,current_walker) = yq(iq)
                quadx(3,ntmove_pts,i,current_walker) = zq(iq)
              endif

              if(iperiodic.eq.0) then
                x(1,i)=r_en(i,ic)*xq(iq)+cent(1,ic)
                x(2,i)=r_en(i,ic)*yq(iq)+cent(2,ic)
                x(3,i)=r_en(i,ic)*zq(iq)+cent(3,ic)
              else
                x(1,i)=r_en(i,ic)*xq(iq)+cent(1,ic)+rshift(1,i,ic)
                x(2,i)=r_en(i,ic)*yq(iq)+cent(2,ic)+rshift(2,i,ic)
                x(3,i)=r_en(i,ic)*zq(iq)+cent(3,ic)+rshift(3,i,ic)
              endif

! Since we are rotating on sphere around nucleus ic, that elec-nucl distance does not change but distances to other nuclei do
              do 40 jc=1,ncent
                do 38 k=1,ndim
   38             rvec_en(k,i,jc)=x(k,i)-cent(k,jc)
                if(jc.eq.ic) cycle

                if(jc.ne.ic) then
                  if(iperiodic.eq.0) then
                    r_en(i,jc)=0
                    do 39 k=1,ndim
   39                 r_en(i,jc)=r_en(i,jc)+rvec_en(k,i,jc)**2
                    r_en(i,jc)=dsqrt(r_en(i,jc))
                  else
                    call find_image4(rshift(1,i,jc),rvec_en(1,i,jc),r_en(i,jc),rlatt,rlatt_inv)
                  endif

                  call scale_dist(r_en(i,jc),rr_en(i,jc),1)
                  call scale_dist(r_en(i,jc),rr_en2(i,jc),3)
                endif
   40         continue

              iel=i

              electron = iel !JT
              call object_modified_by_index (electron_index) !JT

!              call nonlocd(iel,x(1,i),rvec_en,r_en,deter) TA
              !update orbe (was called in nonlocd)
              if(iperiodic.eq.0) then
                if(inum_orb.eq.0) then
                  call orbitals_loc_anae(iel,rvec_en,r_en,orbe)
                 else
                  call orbitals_loc_nume(x(1,i),orbe)
                endif
               else
                if(inum_orb.eq.0) then
                  call orbitals_pwe(iel,x(1,i),orbe)
                 else
                  call orbitals_period_nume(x(1,i),orbe)
                endif
              endif
              call object_modified_by_index (orbe_index) !JT

              call nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,value)
              if(ido_pjas.eq.1) then ! periodic Jastrow implemented by WAS
                call nonloc_pjas (iel, x(:,1:nelec), value)
              endif

              if (l_mode_dmc) then
!                call nonlocd(iel,x(1,i),rvec_en,r_en,deter)
                if (iel.le.nup) then
                  call object_provide_by_index(gup_index)
                  psid_ratio = sum(gup(:,i)*orbe(occup))
                else
                  call object_provide_by_index(gdn_index)
                  psid_ratio = sum(gdn(:,i-nup)*orbe(occdn))
                endif
                quadr(ntmove_pts,i,current_walker) = psid_ratio*exp(value)
!                quadr(ntmove_pts,i,current_walker) = &
!                deter*exp(value)/psidow(current_walker,1)
              endif

              if(ipr.ge.4) then
                write(6,'(''rr_en,rr_en2'',2d14.6)') rr_en(1,1),rr_en2(1,1)
                write(6,'(''ic,i,iq,value'',3i3,2d14.6)') ic,i,iq,value
              endif

              vpsp_tmove=0d0
              do 50 l=1,npotd(ict)
                if(l.ne.lpotp1(ict)) then
                  vpsp_tmove=vpsp_tmove + vps(i,ic,l)*yl0(l,costh)
!                  if(l_do_tmoves) then
!                    if(ipr.ge.1) write(6,'(''vps(i,ic,l),wq(iq),yl0(l,costh),deter,psidow(current_walker,1),exp(value),tau'',9d12.4)') &
!     &                vps(i,ic,l),wq(iq),yl0(l,costh),deter,psidow(current_walker,1),exp(value),tau
!!                   gpsp_tmove(ntmove_pts)=gpsp_tmove(ntmove_pts) + vps(i,ic,l)*yl0(l,costh)                  ! linear approx in Casula
!                    gpsp_tmove(ntmove_pts)=gpsp_tmove(ntmove_pts) + (exp(-tau*vps(i,ic,l))-1)*yl0(l,costh) ! exact
!                    call systemflush(6)
!                  endif
!                  if(ipr.ge.1) write(6,'(''nonloc: l,yl0(l,costh),deter,exp(value),yl0(l,costh)*deter*exp(value)'',i3,9f20.15)') &
!     &            l,yl0(l,costh),deter,exp(value),yl0(l,costh)*deter*exp(value)
                  if(ipr.ge.1) write(6,'(''nonloc: l,yl0(l,costh),exp(value)'',i3,9f20.15)') &
     &            l,yl0(l,costh),exp(value)
                endif
   50         continue ! npotd(ict)
              factor=wq(iq)*exp(value)
!              vpsp=vpsp+factor*deter*vpsp_tmove

              do iorb=1,orb_tot_nb !TA
                psp_nonloc_orb(iel,iorb) = &
                psp_nonloc_orb(iel,iorb) + factor*orbe(iorb)*vpsp_tmove
              enddo
   60       continue ! nquad

! Restore electron i and its distances to value it had before doing angular integration
            do 68 k=1,ndim
   68         x(k,i)=xsav(k)
            do 70 jc=1,ncent
              r_en(i,jc)=r_en_sav(jc)
              rr_en(i,jc)=rr_en_sav(jc)
              rr_en2(i,jc)=rr_en2_sav(jc)
              do 70 k=1,ndim
                rshift(k,i,jc)=rshift_sav(k,jc)
   70           rvec_en(k,i,jc)=rvec_en_sav(k,jc)

          endif  ! iskip, i.e. electron within r_c for this nucleus
  100   continue ! ncent

! Warning: For the moment correlated sampling for forces does not work with t-moves.
! Decide which if any tmove to perform for electron i
        call systemflush(6)
  150 continue ! nelec

! Note that above we called hpsiedmc but not hpsi.  So, the velocity is updated after the tmove but not the energy.
! So, in the reweighting we are using the energies from before and after the drift-dfus-acc/rej step and not after tmove.
! This is slightly different from what is presented on pg. 5 of the Casula et al. paper, if I take it literally.

!      call object_modified_by_index (vpsp_ex_index) ! JT

      call object_modified_by_index (psp_nonloc_orb_index) !TA

      return
      end
!-----------------------------------------------------------------------

      function yl0(l,costh)
! (2L+1)*P_L(costh)
! (2L+1)*P_L(costh)
! This is not quite Y_L0 but sqrt[4pi*(2L+1)]Y_L0 = (2L+1)*P_L(costh)
! Note that the associated P_L^m and the unassociated P_L Legendre polynomials are the same for m=0.
! l is actually L+1.

      implicit real*8(a-h,o-z)

      if(l.eq.1) then
        yl0=1.d0
       elseif(l.eq.2) then
        yl0=3.d0*costh
       elseif(l.eq.3) then
        yl0=2.5d0*(3*costh*costh-1)
       elseif(l.eq.4) then
        yl0=3.5d0*costh*(5*costh*costh-3)
       elseif(l.eq.5) then
        yl0=1.125d0*(35*costh**4-30*costh**2+3)
       else
        stop 'yl0 implemented to l=4 only (Warning: l is l+1)'
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine nonlocd(iel,x,rvec_en,r_en,determ)
! Written by Claudia Filippi, modified by Cyrus Umrigar
      use all_tools_mod
      use control_mod
      use eloc_mod
      use dorb_mod
      use slatn_mod
      use orbe_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use contr2_mod
      use contrl_opt2_mod
      use wfsec_mod
      use contrl_per_mod
      use contr3_mod
      use phifun_mod
      use const_mod
      use slatn2_mod
      use dete_mod
      implicit real*8(a-h,o-z)

      dimension x(3),rvec_en(3,nelec,*),r_en(nelec,*)

      real(dp), allocatable :: ui(:), vi(:)
      real(dp), pointer :: det, inv(:) 
      integer , pointer :: ref(:), exc(:)
      real(dp) :: cdet, dchi
      integer  :: iup, idn, sgn, idet

      if (iel.le.nup) then
        allocate(ui(nup), vi(nup))
      else
        allocate(ui(ndn), vi(ndn))
      endif

!     determ=0

! get orbitals for all electron iel
! Note x is 3-dim but rvec_en,r_en have info about all electrons.
! So, iel is passed to select elements of rvec_en,r_en,phin and for IO.
      if(iperiodic.eq.0) then

        if(inum_orb.eq.0) then
          call orbitals_loc_anae(iel,rvec_en,r_en,orbe)
         else
          call orbitals_loc_nume(x,orbe)
        endif

       else

        if(inum_orb.eq.0) then
          call orbitals_pwe(iel,x,orbe)
         else
          call orbitals_period_nume(x,orbe)
        endif

      endif

      call object_modified_by_index (orbe_index) !JT

      call dete_init(iel)

      if (iel.le.nup) then
        rat = dot_product(orbe(ref_orbs_up), aiup(:,iel))
        ui  = -aiup(:,iel)/rat
        vi  = matmul(orbe(ref_orbs_up), tup(:,occup)) - orbe(occup)

        !update reference det w/ Sherman-Morrison
        deta_upn = deta_upn*rat

        !update tupn w/ Sherman-Morrison 
!        do i=1,orb_tot_nb 
        do i=1,noccup
          do j=1,nup
            tupn(j,i) = tupn(j,i) + ui(j)*vi(i) 
          enddo 
        enddo 

        !update each alpha matrix
        do i=1,ndetup
          k   =   ordup(i)
          ref =>  refup(1:k,i)
          exc =>  excup(1:k,i)
          inv => invupn(1:k*k,i)
          det => detupn(i)

          do ii=1,k
            do jj=1,k
              inv(jj+k*(ii-1)) = tupn(ref(jj),exc(ii))
            enddo
          enddo

          call matinv(inv, k, det)
!          call matinve(inv, k, det, ui(ref), vi(exc))
!          call det_sherman_morrison(det, inv, k, ui(ref), vi(exc))
        enddo
      else
        rat = dot_product(orbe(ref_orbs_dn), aidn(:,iel-nup))
        ui  = -aidn(:,iel-nup)/rat
        vi  = matmul(orbe(ref_orbs_dn), tdn(:,occdn)) - orbe(occdn)

        deta_dnn = deta_dnn*rat

        !update tdnn w/ Sherman-Morrison 
!        do i=1,orb_tot_nb 
        do i=1,noccdn
          do j=1,ndn
            tdnn(j,i) = tdnn(j,i) + ui(j)*vi(i) 
          enddo 
        enddo 

        do i=1,ndetdn
          k   =   orddn(i)
          ref =>  refdn(1:k,i)
          exc =>  excdn(1:k,i)
          inv => invdnn(1:k*k,i)
          det => detdnn(i)

          do ii=1,k
            do jj=1,k
              inv(jj+k*(ii-1)) = tdnn(ref(jj),exc(ii))
            enddo
          enddo

          call matinv(inv, k, det)
!          call matinve(inv, k, det, ui(ref), vi(exc))
!          call det_sherman_morrison(det, inv, k, ui(ref), vi(exc))
        enddo
      endif

      chin = 0d0
      do idet=1,ndet
        cdet = cdet_in_wf(idet, iwf)
        iup  = iwdetup(idet)
        idn  = iwdetdn(idet)
        sgn  = sgnup(iup)*sgndn(idn)
        if (iel.le.nup) then
            dchi = sgn*cdet*detupn(iup)*detdn(idn)
        else
            dchi = sgn*cdet*detup(iup)*detdnn(idn)
        endif
        chin = chin + dchi
      enddo

      if (iel.le.nup) then
        determ = chin*deta_upn*deta_dn
      else
        determ = chin*deta_dnn*deta_up
      endif

      return
      end
!-----------------------------------------------------------------------

!WAS  subroutine nonlocj(iel,x,rshift,rr_en,rr_en2,value)
      subroutine nonlocj(iel,x,rshift,r_en,rr_en,rr_en2,value)
! fso are pair contributions in the Jastrow exponent before electron iel moves
! fsn are pair contributions in the Jastrow exponent after  electron iel moves
! value=sum(fsn-fso), i.e., the change in the Jastrow exponent when electron iel moves
! Written by Claudia Filippi, modified by Cyrus Umrigar

      use control_mod
      use atom_mod
      use dets_mod
      use const_mod
      use dim_mod
      use contrl_per_mod
      use jaspar_mod
      use bparm_mod
      use periodic_mod
      use jaso_mod
      implicit real*8(a-h,o-z)

      dimension x(3,*),rshift(3,nelec,ncent),rr_en(nelec,ncent),rr_en2(nelec,ncent),fsn(nelec,nelec),dx(3)

      dimension r_en(nelec,ncent)

      fsumn=0

      if(nelec.lt.2) goto 47

      do 45 jj=1,nelec

        if(jj.eq.iel) goto 45
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif

        sspinn=1
        ipar=0
        isb=1
        if(i.le.nup .or. j.gt.nup) then
          if(nspin2b.eq.2) then
            isb=2
           elseif(nocuspb.eq.0) then
            sspinn=half
          endif
          ipar=1
        endif

        do 10 k=1,ndim
   10     dx(k)=x(k,jj)-x(k,iel)

        if(iperiodic.eq.0) then
          rij=0
          do 20 k=1,ndim
   20       rij=rij+dx(k)**2
          rij=dsqrt(rij)
         else
          call find_image3(dx,rij,rlatt_sim,rlatt_sim_inv)
        endif

! e-e terms
        call scale_dist(rij,u,2)

        fsn(i,j)=psibnl(u,isb,ipar)

! e-e-n terms
! The scaling is switched in psinl, so do not do it here.
!     if(isc.ge.12) call scale_dist(rij,u,3)
      call scale_dist(rij,u,4)

        do 40 ic=1,ncent
          it=iwctype(ic)
!WAS 40   fsn(i,j)=fsn(i,j) + psinl(u,rshift(1,i,ic),rshift(1,j,ic),rr_en2(i,ic),rr_en2(j,ic),it)
   40     fsn(i,j)=fsn(i,j) + psinl(u,rshift(1,i,ic),rshift(1,j,ic),r_en(i,ic),r_en(j,ic),rr_en2(i,ic),rr_en2(j,ic),it)
        fsumn=fsumn+fsn(i,j)-fso(i,j)

   45 continue

! e-n terms
   47 fsn(iel,iel)=0
      do 50 ic=1,ncent
        it=iwctype(ic)
   50   fsn(iel,iel)=fsn(iel,iel)+psianl(rr_en(iel,ic),it)

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)

      value=fsumn

      return
      end
