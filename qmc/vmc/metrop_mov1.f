      subroutine metrop_mov1(ipass)
!     subroutine metrop(ipass)
! Written by Cyrus Umrigar and Claudia Filippi
! routine to move configuration by one step using a
! velocity bias type Monte Carlo.
! This routine has been changed to move one electron at a time and
! to keep track of the acceptance prob. as a function of delta
! Minor mods by A.D.Guclu to include pair-density function calculation
      use all_tools_mod
      use constants_mod
      use control_mod
      use montecarlo_mod
      use eloc_mod
      use opt_lin_mod
      use opt_ptb_mod
      use config_mod
      use dets_mod
      use const_mod
      use dim_mod
      use forcepar_mod
!      use doefp_mod
      use contrl_per_mod
      use delocc_mod
      use div_v_mod
      use denupdn_mod
      use stepv_mod
      use kinet_mod
      use pairden_mod
      use fourier_mod
      use forcewt_mod
      use estsig_mod
      use estsum_mod
      use determinants_mod
      use distance_mod
      use zigzag_mod, only: izigzag
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-10)

      logical vgreater

! Warning: program has a bug, but it does not matter much since
! we never use it.

!     common /contr3/ mode
      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys

      dimension dx(3,nelec),idist(nelec),ixo(3),ixn(3)
      dimension xstrech(3,nelec)

!     write(6,'(''entering metrop_mov1'')')

!     mode='vmc_mov1    '

! Sample transition probability from current state to new state
! and record value of probability in fxop.
! The transition probabilty is given by statement label 17
! for the reverse transition.
! The transition probability is an approximately psi(new)/psi(old)
! and is positive definite.

      deltfi=two*fbias*deltai
      do 300 i=1,nelec
        if(i.ne.ifixe) then

          fxop=one
          do 10 k=1,ndim
            fx=dmin1(dabs(vold(k,i)),deltfi)
            ftest=fx*delta*half
            fx=sign(fx,vold(k,i))
            if(ftest.le.rannyu(0)) then
!           sample from 1 piece
              dx(k,i)=(rannyu(0)-half)*delta
             else
!             sample from the (xnew-xold) piece assuming velocity positive
              dx(k,i)=(dmax1(rannyu(0),rannyu(0))-half)*delta
!             if velocity was negative then dx has opposite sign
              if(vold(k,i).lt.zero) dx(k,i)=-dx(k,i)
            endif
            fxop=fxop*(one+fx*dx(k,i))
   10       xnew(k,i)=xold(k,i)+dx(k,i)

! calculate psi etc. at new configuration
          iel=i
!         call hpsi(xnew,psidn,psijn,vnew,div_vn,d2,pen,enew(1),denergy,1)
!         write(6,'(''calling hpsie from metrop_mov1'')')
          call hpsie(iel,xnew,psidn,psijn,vnew)
          psi2n(1)=2*(dlog(dabs(psidn))+psijn)

! save electrostatic potential at new configuration
!   distances.f is not called from from hpsie, so we add a call to distances here,
!   since we need it to compute pot_ee(i) (ACM)

!          call distances(xnew,pe_dummy,pei_dummy)
          pot_ee_new = pot_ee  ! this is an array assignment

! If error is large then save config. to use in optimizing routine

!         if(dabs((enew(1)-etrial)/etrial).gt.1.d-3) then
!c        if(dabs((enew(1)-etrial)/etrial).gt.1.d+1) then
!         write(8,'(i6,f8.2,2d10.2,(8f8.4))') ipass,enew(1)-etrial,psin,
!    &    (enew(1)-etrial)*psin,((xnew(k,jj),k=1,ndim),jj=1,nelec)
!         endif
!         if(dabs(psin*(enew(1)-etrial)/etrial).gt.1.d-7) then
!           write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass,enew(1)-etrial,psin,
!    &      (enew(1)-etrial)*psin,((xnew(k,jj),k=1,ndim),jj=1,nelec)
!         endif

! calculate probability for reverse transition
          fxnp=one
          do 20 k=1,ndim
   20       fxnp=fxnp*(one-dx(k,i)
     &      *sign(dmin1(dabs(vnew(k,i)),deltfi),vnew(k,i)))

! form the Jackson Feenberg kinetic energy expression

!         tjfn=d2
!         tjfn=-tjfn*half*hb

! p is the probability of accepting new move

          p=exp(psi2n(1)-psi2o(1))*fxnp/fxop
!         write(6,*) ((xold(kr,ir),kr=1,ndim),ir=1,nelec)
!         write(6,*) ((xnew(kr,ir),kr=1,ndim),ir=1,nelec)
!         write(6,'(''psi2n,psi2o,fxnp,fxop,psi2n/psi2o,fxnp/fxop,p'',9f12.6
!    &    )')psi2n(1),psi2o(1),fxnp,fxop,psi2n(1)/psi2o(1),fxnp/fxop,p
!         write(6,'(''psido, psidn, psijo, psijn'',4f12.6)') psido, psidn, psijo, psijn
          p=dmin1(one,p)
          q=one-p

! form expected values of e, pe, etc.

!     esum1=          p*enew(1)+q*eold(1)
!     esum(1)=esum(1)+p*enew(1)+q*eold(1)
!     pesum=pesum+p*pen+q*peo
!     tpbsum=tpbsum+p*(enew(1)-pen)+q*(eold(1)-peo)
!     tjfsum=tjfsum+p*tjfn+q*tjfo
!     do 25 k=1,ndim
!  25     r2sum=r2sum+p*xnew(k,i)**2+q*xold(k,i)**2
!     if(nefp.gt.0) then
!       call sample_efp(0,xold,eold(1),q)
!       call sample_efp(1,xnew,enew(1),p)
!     endif
!     call acues1

! Calculate as a function of the distance to the nucleus
! 1) acceptance,  2) force-bias truncation probability,
! 3) kinetic energy and it's fluctuation
! The K.E. is not quite correct, since we should use p times new
! and q times old, and keep track of which bin the old was in
! The reason why I changed
! itryo=min(int(delri*rold)+1,NRAD)  to
! itryo=int(min(delri*rold+1,dfloat(NRAD))+eps)
! is that 2147483647 is the largest 32-bit integer and 1 more than that gives -2147483648.
          rold=0
          rnew=0
          do 100 idim=1,ndim
            rold=rold+xold(idim,i)**2
  100       rnew=rnew+xnew(idim,i)**2
          rold=dsqrt(rold)
          rnew=dsqrt(rnew)
          itryo=int(min(delri*rold+1,dfloat(NRAD))+eps)
          itryn=int(min(delri*rnew+1,dfloat(NRAD))+eps)
          try(itryo)=try(itryo)+1
          suc(itryo)=suc(itryo)+p
! this was for ndim=3
!      if((dabs(vold(1,i)).gt.deltfi).or.(dabs(vold(2,i)).gt.deltfi)
!     &.or.(dabs(vold(3,i)).gt.deltfi)) trunfb(itryo)=trunfb(itryo)+1
! for general ndim
          vgreater=.false.
          do 150 idim=1,ndim
  150       if(dabs(vold(idim,i)).gt.deltfi) vgreater=.true.
          if(vgreater) trunfb(itryo)=trunfb(itryo)+1
          if(i.le.nup) then
            rprobup(itryo)=rprobup(itryo)+q
            rprobup(itryn)=rprobup(itryn)+p
          else
            rprobdn(itryo)=rprobdn(itryo)+q
            rprobdn(itryn)=rprobdn(itryn)+p
          endif
          rprob(itryo)=rprob(itryo)+q
          rprob(itryn)=rprob(itryn)+p
          ekin(itryo)=ekin(itryo)+q*ekineo(i)
          ekin2(itryo)=ekin2(itryo)+q*ekineo(i)**2
          ekin(itryn)=ekin(itryn)+p*ekinen(i)
          ekin2(itryn)=ekin2(itryn)+p*ekinen(i)**2

! calculate 2d-density related functions
          if(iperiodic.eq.1 .or. iperiodic.eq.3) then  ! In 1D and 3D reduce position to simulation cell
            call reduce_sim_cell(xold(:,i))
            call reduce_sim_cell(xnew(:,i))
          endif
          if(izigzag.gt.0) call zigzag2d(p,q,xold,xnew,i)
          if(ifixe.le.-2) call pairden2d(p,q,xold,xnew)
          if(ifourier.eq.1 .or. ifourier.eq.3) call fourierrk(p,q,xold,xnew)
          if(ifourier.eq.2 .or. ifourier.eq.3) call fourierkk(p,q,xold,xnew)
!          if(ifourier.eq.1) call fourier2d(1.d0,0.d0,xold,xnew)
          if(ifixe.ne.0 .and. ifixe.ne.-2) then   !  2d density (ifixe=-1,-3) or fixed electron pair density (ifixe>0)
            if(icoosys.eq.1) then
              do 170 idim=1,ndim
! note that ix can be negative or positive. nint is a better choice.
                ixo(idim)=nint(delxi(idim)*xold(idim,i))
  170           ixn(idim)=nint(delxi(idim)*xnew(idim,i))
             else
! same trick adapted to circular coordinates
               ixo(1)=nint(delradi*(rold-rmean))
               ixn(1)=nint(delradi*(rnew-rmean))
               ixo(2)=nint(delti*(datan2(xold(2,i),xold(1,i))))
               ixn(2)=nint(delti*(datan2(xnew(2,i),xnew(1,i))))
            endif

            if(abs(ixo(1)).le.NAX .and. abs(ixo(2)).le.NAX) then
              den2d_t(ixo(1),ixo(2))=den2d_t(ixo(1),ixo(2))+q
              pot_ee2d_t(ixo(1),ixo(2))=pot_ee2d_t(ixo(1),ixo(2))+q*pot_ee_old(i)
              if(i.le.nup) then
                den2d_u(ixo(1),ixo(2))=den2d_u(ixo(1),ixo(2))+q
                pot_ee2d_u(ixo(1),ixo(2))=pot_ee2d_u(ixo(1),ixo(2))+q*pot_ee_old(i)
               else
                den2d_d(ixo(1),ixo(2))=den2d_d(ixo(1),ixo(2))+q
                pot_ee2d_d(ixo(1),ixo(2))=pot_ee2d_d(ixo(1),ixo(2))+q*pot_ee_old(i)
              endif
            endif
            if(abs(ixn(1)).le.NAX .and. abs(ixn(2)).le.NAX) then
              den2d_t(ixn(1),ixn(2))=den2d_t(ixn(1),ixn(2))+p
              pot_ee2d_t(ixn(1),ixn(2))=pot_ee2d_t(ixn(1),ixn(2))+p*pot_ee_new(i)
              if(i.le.nup) then
                den2d_u(ixn(1),ixn(2))=den2d_u(ixn(1),ixn(2))+p
                pot_ee2d_u(ixn(1),ixn(2))=pot_ee2d_u(ixn(1),ixn(2))+p*pot_ee_new(i)
               else
                den2d_d(ixn(1),ixn(2))=den2d_d(ixn(1),ixn(2))+p
                pot_ee2d_d(ixn(1),ixn(2))=pot_ee2d_d(ixn(1),ixn(2))+p*pot_ee_new(i)
              endif
            endif
          endif

          r2sumeo = 0.
          r2sumen = 0.
          do k=1,ndim
            r2sumeo = r2sumeo + xold(k,i)**2
            r2sumen = r2sumen + xnew(k,i)**2
          enddo
          r1sum=r1sum+p*dsqrt(r2sumen) + q*dsqrt(r2sumeo)
          r2sum=r2sum+ p*r2sumen + q*r2sumeo
          r3sum=r3sum+p*r2sumen*dsqrt(r2sumen) + q*r2sumeo*dsqrt(r2sumeo)
          r4sum=r4sum+ p*r2sumen*r2sumen + q*r2sumeo*r2sumeo

!     eksum=zero
!     do 26 ii=1,nelec
!  26   eksum=eksum+ekine(ii)
!     write(6,'(''ke='',9d13.5)') eksum,enew(1)-pen


! accept new move with probability p
! Note when one electron moves the velocity on all electrons change.
          if(rannyu(0).lt.p) then
            idist(i)=itryn
            psi2o(1)=psi2n(1)
            pot_ee_old = pot_ee_new
            do 240 k=1,ndim
              xold(k,i)=xnew(k,i)
              rvmino(k,i)=rvminn(k,i)
              do 240 ii=1,nelec
  240           vold(k,ii)=vnew(k,ii)
            psido=psidn
            psijo=psijn
            accsum=accsum+one
            call jassav(i)
            if(ibasis.eq.3) then
              call cdetsav(i)         ! complex case
             else
              call detsav(i)
            endif

           else
            idist(i)=itryo
            do 250 k=1,ndim
  250         xnew(k,i)=xold(k,i)
            call distancese_restore(i)
          endif

        endif
  300 continue

      call object_modified_by_index (xold_index)  !JT

! loop over secondary configurations
      do 350 ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
!     write(6,'(''calling hpsi from metrop_mov1,ifr,nforce'',9i5)') ifr,nforce
        call hpsi(xstrech,psido,psijo,vold,div_vo,d2,peo,peio,eold(ifr),denergy,ifr)
  350   psi2o(ifr)=2*(dlog(dabs(psido))+psijo)+dlog(ajacob)

! primary configuration
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido,psijo,vold,div_vo,d2o,peo,peio,eold(1),denergy,1)
      psi2o(1)=2*(dlog(dabs(psido))+psijo)

      eloc = eold(1)                 !JT
      psi_det = psido                !JT
      call object_modified_by_index (eold_index)  !JT
      call object_modified_by_index (eloc_index)  !JT
      call object_modified_by_index (psi_det_index) !JT
      call object_modified_by_index (psijo_index)  !JT
      call object_modified_by_index (denergy_index) !JT
      call object_modified_by_index (vold_index) !JT
      call object_modified_by_index (div_vo_index) !JT

      tjfo=-d2o*half*hb
! form expected values of e, pe, etc.
      esum1=eold(1)
      esum(1)=esum(1)+eold(1)
      pesum=pesum+peo
      peisum=peisum+peio
      tpbsum=tpbsum+(eold(1)-peo)
      tjfsum=tjfsum+tjfo
!      if(nefp.gt.0) then
!        call sample_efp(1,xold,eold(1),1.d0)
!        call efpsav
!      endif

      call grad_hess_jas_sum(1.d0,0.d0,eold(1),eold(1),1.d0,0.d0)

      call acues1

      do 380 ifr=1,nforce
        if(ifr.eq.1) then
          esum1s(ifr)=eold(ifr)
         else
          wstro=max(min(exp(psi2o(ifr)-psi2o(1)),1.d99),1.d-99)
          wsum1s(ifr)=wstro
          esum(ifr)=esum(ifr)+eold(ifr)*wstro
          wsum(ifr)=wsum(ifr)+wstro
          esum1s(ifr)=eold(ifr)*wstro
! psi2o has the log of the wavefn.
!         if(abs(psi2o(ifr)-psi2o(1)).gt.50.d0)
!    &    write(6,'(''ifr,esum(ifr),wstro'',i3,100d12.4)')
!    &    ifr,esum(ifr),wstro,psi2o(ifr),psi2o(1)
        endif
  380 continue

      call acusig

      do 390 i=1,nelec
        if(i.ne.ifixe) then
          ekineo(i)=ekinen(i)
          ekin(idist(i))=ekin(idist(i))+ekineo(i)
          ekin2(idist(i))=ekin2(idist(i))+ekineo(i)**2
        endif
  390 continue

      return
      end
