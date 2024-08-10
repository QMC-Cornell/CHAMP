subroutine tmove
! Written by Tyler Anderson by modifying nonloc.f90
! This routine is called by dmc. 
! Warning: Assumes for the reverse T-move probability that the nonlocal ranges of the atoms are non overlapping.
! npotd(ict) is the number of channels, local and nonlocal, for nucleus type ict. In most of the code npotd=npotu and only npotd is used.
! lpotp1=L+1, where L is the local channel.  It would have been more clear to start some arrays at 0, in which case the +1 would not be needed.
  use types_mod, only: dp
  use average_mod, only: current_walker
  use distance_mod, only: r_en, rvec_en, rshift
  use const_mod, only: nelec, ipr
  use config_dmc_mod, only: xoldw, voldw, psijow, psidow
  use atom_mod, only: ncent, cent, iwctype
  use qua_mod, only: nquad, xq, yq, zq, wq, iaccept_tmove, quadr, quadx
  use contrldmc_mod, only: tau
  use pseudo_mod, only: npotd, lpotp1, vps
  use contrl_per_mod, only: iperiodic, ibasis
  use control_mod, only: tmoves_type
  use dete_mod, only: dete_save
  use debug_shared_mod, only: is, iw
  implicit none

!  integer  :: iw, i, ic, ict, iq, jq, ntmove_pts, iwhich_tmove, l, ioffset
  integer  :: i, ic, ict, iq, jq, ntmove_pts, iwhich_tmove, l, ioffset
  integer  :: iwhich_cent_tmove(nquad*ncent)
  integer  :: discrete_rand
  real(dp) :: psidn, psijn, p, ratio, costh
  real(dp) :: gpsp, sum_gpsp_fwd, sum_gpsp_bwd, gpsp_tmove(0:nquad*ncent)
  real(dp) :: psi_ratio(nquad*ncent), xnew(3), vnew(3,nelec)
  real(dp) :: xtmove(3), x_tmove_sav(3,nquad*ncent), x_en(3)
!  real(dp) :: wfwd, w_tmove_sav(nquad*ncent)
  real(dp) :: rr_en(nelec, ncent), rr_en2(nelec, ncent)
  real(dp) :: rannyu, yl0
  logical  :: skip, load_quad

  do i=1,nelec
    do ic=1,ncent
      call scale_dist(r_en(i,ic),rr_en(i,ic),1)
      call scale_dist(r_en(i,ic),rr_en2(i,ic),3)
    enddo
  enddo

  if (iperiodic.ne.0) stop 'tmove not working on periodic systems.'

!  iw=current_walker
  current_walker=iw !needed in quad_psi subroutine
  iaccept_tmove=0
  do i=1,nelec
    !Save old distances in r_en_sav, rvec_en_sav, r_ee_sav, rvec_ee_sav
    call distancese(i,xoldw(1,1,iw,1))

    !Calculate tmove probs
    gpsp_tmove(0) = 1d0
    ntmove_pts=0
    do ic=1,ncent
      ict=iwctype(ic)

      skip=.true.
      do l=1,npotd(ict)
        if(l.ne.lpotp1(ict) .and. dabs(vps(i,ic,l)).gt.1.d-4) skip=.false.
      enddo
      if (skip) cycle

      if (iaccept_tmove.eq.1) call rotqua

      do iq=1,nquad
        ntmove_pts=ntmove_pts+1

        if (iaccept_tmove.ne.1) then
          xq(iq) = quadx(1,ntmove_pts,i)
          yq(iq) = quadx(2,ntmove_pts,i)
          zq(iq) = quadx(3,ntmove_pts,i)
        endif

        costh=(rvec_en(1,i,ic)*xq(iq) &
              +rvec_en(2,i,ic)*yq(iq) &
              +rvec_en(3,i,ic)*zq(iq))/r_en(i,ic)

        xtmove(1)=r_en(i,ic)*xq(iq)+cent(1,ic)
        xtmove(2)=r_en(i,ic)*yq(iq)+cent(2,ic)
        xtmove(3)=r_en(i,ic)*zq(iq)+cent(3,ic)
        if (iperiodic.ne.0) xtmove(:) = xtmove(:) + rshift(:,i,ic)

        if (iaccept_tmove.ne.1) then
          ratio = quadr(ntmove_pts,i)
        else
          call quad_psi(i, ic, rr_en, rr_en2, xtmove, ratio)
        endif

        gpsp=0d0
        do l=1,npotd(ict)
          if (l.eq.lpotp1(ict)) cycle
          gpsp=gpsp+wq(iq)*(exp(-tau*vps(i,ic,l))-1)*yl0(l,costh)*ratio
!         gpsp=gpsp-wq(iq)*tau*vps(i,ic,l)*yl0(l,costh)*ratio                  ! linear approx in Casula
        enddo

        x_tmove_sav(:,ntmove_pts)=xtmove(:)
        psi_ratio(ntmove_pts)=ratio
        gpsp_tmove(ntmove_pts)=max(0d0, gpsp)
        iwhich_cent_tmove(ntmove_pts)=ic
      enddo !iq=1,nquad
    enddo !ic=1,ncent

    !propose a tmove
    iwhich_tmove = discrete_rand(gpsp_tmove(0:ntmove_pts), ntmove_pts+1)-1
    !no tmove
    if (iwhich_tmove.eq.0) cycle
    !do tmove
    ic = iwhich_cent_tmove(iwhich_tmove)
    ioffset = int((iwhich_tmove-1)/nquad)*nquad
    iq = 1+mod( iwhich_tmove-1 ,nquad)
    ict= iwctype(ic)
    xnew(:)=x_tmove_sav(:,iwhich_tmove)
    x_en(:)=xnew(:)-cent(:,ic)

    if (.not.(tmoves_type.eq.'casula1')) then
      sum_gpsp_fwd=sum(gpsp_tmove(0:ntmove_pts))
      !resample backward move
      sum_gpsp_bwd=1d0
      call rotqua
      do jq=1,nquad
        !warning: r_en only valid for this ic.
        costh=(x_en(1)*xq(jq) &
              +x_en(2)*yq(jq) &
              +x_en(3)*zq(jq))/r_en(i,ic)
  
        xtmove(1)=r_en(i,ic)*xq(jq)+cent(1,ic)
        xtmove(2)=r_en(i,ic)*yq(jq)+cent(2,ic)
        xtmove(3)=r_en(i,ic)*zq(jq)+cent(3,ic)
        if (iperiodic.ne.0) xtmove(:) = xtmove(:) + rshift(:,i,ic)
  
        call quad_psi(i, ic, rr_en, rr_en2, xtmove, ratio)
        !after line above, ratio = psi(xtmove)/psi(xoldw(:,:,iw,1))
        ratio=ratio/psi_ratio(iwhich_tmove)
        !after line above, ratio = psi(xtmove)/psi(iwhich_tmove)
  
        !warning: vps only correct for this ic.
        gpsp=0d0
        do l=1,npotd(ict)
          if (l.eq.lpotp1(ict)) cycle
          gpsp=gpsp+wq(jq)*(exp(-tau*vps(i,ic,l))-1)*yl0(l,costh)*ratio
!         gpsp=gpsp-wq(jq)*tau*vps(i,ic,l)*yl0(l,costh)*ratio                  ! linear approx in Casula
        enddo
        sum_gpsp_bwd=sum_gpsp_bwd+max(0d0,gpsp)
      enddo !jq=1,nquad
      p = min(1d0, sum_gpsp_fwd/sum_gpsp_bwd)
    endif !.not.(tmoves_type.eq.'casula1')
    !accept/reject tmove

    if((rannyu(0).lt.p).or.(tmoves_type.eq.'casula1')) then
      iaccept_tmove=1
      call hpsiedmc(i,iw,xnew,psidn,psijn,vnew)
      xoldw(:,i,iw,1)=xnew(:)
      voldw(:,:,iw,1)=vnew(:,:)
      psidow(iw,1)=psidn
      psijow(iw,1)=psijn
      call jassav(i)
      if(ibasis.eq.3) then
        call cdetsav(i)
      else
        call detsav(i)
        call dete_save(i)
      endif
    endif !rannyu(0).lt.p
  enddo !i=1,nelec
end subroutine

function discrete_rand(probs, n)
!Returns a random int i between 1 to n inclusive
!with probability probs(i)
  use types_mod, only: dp
  implicit none
  integer :: discrete_rand

  integer , intent(in) :: n
  real(dp), intent(in) :: probs(1:n)

  real(dp) :: psum, rand, rannyu
  integer  :: i

  rand=sum(probs(1:n))*rannyu(0)
  psum=0d0
  do i=1,n
    psum=psum+probs(i)
    if (psum.gt.rand) then
      discrete_rand=i
      return
    endif
  enddo
  return
end function

subroutine quad_psi(i, ic, rr_en, rr_en2, xtmove, ratio)
! exp(lnjas) is the ratio of the Jastrows at the quadrature points and the initial point.
!Calculates the ratio of the wavefunction at a quadrature point located at
!xtmove and the wavefunction at the original point (xoldw(:,:,iw,1))
!This function doesn't change any global variables
!only 'ratio' calculated and returned (rr_en, rr_en2 return to original value) 
  use types_mod, only: dp
  use distance_mod, only: r_en, rvec_en, rshift
  use atom_mod, only: ncent, cent
  use periodic_mod, only: rlatt, rlatt_inv
  use dim_mod, only: ndim
  use average_mod, only: current_walker
  use const_mod, only: nelec
  use contrl_per_mod, only: iperiodic
  use config_dmc_mod, only: xoldw, psidow, psijow
  use pjase_mod, only: ido_pjas
  use periodic_jastrow_mod, only: nonloc_pjas
  implicit none

  integer , intent(in)    :: i, ic
  real(dp), intent(in)    :: xtmove(3)
  real(dp), intent(inout) :: rr_en(nelec, ncent), rr_en2(nelec, ncent)
  real(dp), intent(out)   :: ratio

  integer  :: jc, k
  real(dp) :: det, lnjas, xnew(3,nelec), rr_en_sav(ncent), rr_en2_sav(ncent)

  rr_en_sav(:) = rr_en(i,:)
  rr_en2_sav(:) = rr_en2(i,:)
  do jc=1,ncent
    do k=1,ndim
      rvec_en(k,i,jc)=xtmove(k)-cent(k,jc)
    enddo

    if(jc.eq.ic) cycle

    if(iperiodic.eq.0) then
      r_en(i,jc)=0
      do k=1,ndim
        r_en(i,jc)=r_en(i,jc)+rvec_en(k,i,jc)**2
      enddo
      r_en(i,jc)=dsqrt(r_en(i,jc))
    else
      call find_image4(rshift(1,i,jc),rvec_en(1,i,jc),r_en(i,jc),rlatt,rlatt_inv)
    endif

    call scale_dist(r_en(i,jc),rr_en(i,jc),1)
    call scale_dist(r_en(i,jc),rr_en2(i,jc),3)
  enddo

  xnew(:,:) = xoldw(:,:,current_walker,1)
  xnew(:,i) = xtmove(:)
  call nonlocd(i, xtmove, rvec_en, r_en, det)
  call nonlocj(i, xnew, rshift, r_en, rr_en, rr_en2, lnjas)
  if (ido_pjas.eq.1) then
    call nonloc_pjas(i, xnew(:,:), lnjas)
  endif

  ratio = det*exp(lnjas)/psidow(current_walker,1)

  call distancese_restore(i) 
  rr_en(i,:) = rr_en_sav(:)
  rr_en2(i,:) = rr_en2_sav(:)
end subroutine quad_psi
