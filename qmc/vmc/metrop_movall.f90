      subroutine metrop_movall(ipass)
!     subroutine metrop(ipass)
! Written by Cyrus Umrigar
! routine to move configuration by one step using a
! force-bias type Monte Carlo.
! Minor mods added by A.D.Guclu to include correlated sampling.
      use all_tools_mod
      use constants_mod
      use config_mod
      use dets_mod
      use const_mod
      use dim_mod
      use forcepar_mod
!      use doefp_mod
      use delocc_mod
      use denupdn_mod
      use stepv_mod
      use pairden_mod
      use fourier_mod
      use forcewt_mod
      use estsig_mod
      use estsum_mod
      use eloc_mod
      use contrl_per_mod
      use determinants_mod
      use distance_mod
      use zigzag_mod, only: izigzag

      implicit real*8(a-h,o-z)

      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys

      dimension dx(3,nelec),div_vn(nelec),ixo(3),ixn(3)
      dimension xstrech(3,nelec)
      parameter(eps=1.d-10)

! Sample transition probability from current state to new state
! and record value of probability in fxop.
! The transition probabilty is given by statement label 17
! for the reverse transition.
! The transition probability is an approximately psi(new)/psi(old)
! and is positive definite.


      deltfi=two*fbias*deltai
      fxop=one
      do 10 i=1,nelec
! worry about the fixe electron only in this loop
        if(i.ne.ifixe) then
          do 5 k=1,ndim
            fx=dmin1(dabs(vold(k,i)),deltfi)
            ftest=fx*delta*half
            fx=sign(fx,vold(k,i))
            if(ftest.le.rannyu(0)) then
!           sample from 1 piece
              dx(k,i)=(rannyu(0)-half)*delta
             else
!           sample from the (xnew-xold) piece assuming velocity positive
              dx(k,i)=(dmax1(rannyu(0),rannyu(0))-half)*delta
!           if velocity was negative then dx has opposite sign
              if(vold(k,i).lt.zero) dx(k,i)=-dx(k,i)
            endif
            fxop=fxop*(one+fx*dx(k,i))
            xnew(k,i)=xold(k,i)+dx(k,i)
 5       continue
        else
          do 7 k=1,ndim
            dx(k,i)=0.d0
            xnew(k,i)=xold(k,i)
 7        continue
        endif
 10   continue

! calculate psi etc. at new configuration
! loop over secondary configurations

      do 15 ifr=2,nforce
        call strech(xnew,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psidn,psijn,vnew,div_vn,d2,pen,pein,enew(ifr),denergy,ifr)
 15     psi2n(ifr)=2*(dlog(dabs(psidn))+psijn)+dlog(ajacob)

! primary configuration
      if(nforce.gt.1) call strech(xnew,xstrech,ajacob,1,0)
      call hpsi(xnew,psidn,psijn,vnew,div_vn,d2,pen,pein,enew(1),denergy,1)
      psi2n(1)=2*(dlog(dabs(psidn))+psijn)

! save electrostatic potential at new configuration
!   Note that we already called distances.f from hpsi, otherwise
!   we should add a call to distances here (ACM)

          pot_ee_new = pot_ee  ! this is an array assignment

! If error is large then save config. to use in optimizing routine

!     if(dabs((enew(1)-etrial)/etrial).gt.1.d0.and.iwrit8.le.1000) then
!        iwrit8=iwrit8+1
!        write(8,'(i6,f8.2,2d10.2,(8f8.4))') ipass,enew(1)-etrial,psin,
!    &   (enew(1)-etrial)*psin,((xnew(k,jj),k=1,ndim),jj=1,nelec)
!c       write(8,'(10f8.4)') ((xnew(k,jj),k=1,ndim),jj=1,nelec)
!c    endif
!     if(dabs(psin*(enew(1)-etrial)/etrial).gt.1.d-7) then
!        write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass,enew(1)-etrial,psin,
!    &   (enew(1)-etrial)*psin,((xnew(k,jj),k=1,ndim),jj=1,nelec)
!     endif

! calculate probability for reverse transition

      fxnp=one
      do 20 k=1,ndim
      do 20 i=1,nelec
   20 fxnp=fxnp*(one-dx(k,i) &
     &   *sign(dmin1(dabs(vnew(k,i)),deltfi),vnew(k,i)))

! form the Jackson Feenberg kinetic energy expression

      tjfn=d2
      tjfn=-tjfn*half*hb

! p is the probability of accepting new move

      p=exp(psi2n(1)-psi2o(1))*fxnp/fxop
      p=dmin1(one,p)
      q=one-p

! form expected values of e, pe, etc.

      esum1=          p*enew(1)+q*eold(1)
      esum(1)=esum(1)+p*enew(1)+q*eold(1)
      pesum=pesum+p*pen+q*peo
      peisum=peisum+p*pein+q*peio
      tpbsum=tpbsum+p*(enew(1)-pen)+q*(eold(1)-peo)
      tjfsum=tjfsum+p*tjfn+q*tjfo
      do 25 k=1,ndim
        do 25 i=1,nelec
   25     r2sum=r2sum+p*xnew(k,i)**2+q*xold(k,i)**2
!      if(nefp.gt.0) then
!        call sample_efp(0,xold,eold(1),q)
!        call sample_efp(1,xnew,enew(1),p)
!      endif

      call grad_hess_jas_sum(p,q,enew(1),eold(1),1.d0,0.d0)

      call acues1

      do 26 ifr=2,nforce
        wstro=max(min(exp(psi2o(ifr)-psi2o(1)),1.d99),1.d-99)
        wstrn=max(min(exp(psi2n(ifr)-psi2n(1)),1.d99),1.d-99)
        esum(ifr)=esum(ifr)+p*enew(ifr)*wstrn+q*eold(ifr)*wstro
 26     wsum(ifr)=wsum(ifr)+p*wstrn+q*wstro

! The reason why I changed
! itryo=min(int(delri*rold)+1,NRAD)  to
! itryo=int(min(delri*rold+1,dfloat(NRAD))+eps)
! is that 2147483647 is the largest 32-bit integer and 1 more than that gives -2147483648.
      do 28 i=1,nelec
        rold=dsqrt(xold(1,i)**2+xold(2,i)**2+xold(3,i)**2)
        rnew=dsqrt(xnew(1,i)**2+xnew(2,i)**2+xnew(3,i)**2)
!       itryo=min(int(delri*rold)+1,NRAD)
!       itryn=min(int(delri*rnew)+1,NRAD)
        itryo=int(min(delri*rold+1,dfloat(NRAD))+eps)
        itryn=int(min(delri*rnew+1,dfloat(NRAD))+eps)
        if(i.le.nup) then
          rprobup(itryo)=rprobup(itryo)+q
          rprobup(itryn)=rprobup(itryn)+p
         else
          rprobdn(itryo)=rprobdn(itryo)+q
          rprobdn(itryn)=rprobdn(itryn)+p
        endif
        rprob(itryo)=rprob(itryo)+q
        rprob(itryn)=rprob(itryn)+p
! calculate 2d density related functions:
        if(iperiodic.eq.1 .or. iperiodic.eq.3) then  ! In 1D and 3D reduce position to simulation cell
          call reduce_sim_cell(xold(:,i))
          call reduce_sim_cell(xnew(:,i))
        endif
        if(ifixe.ne.0 .and. ifixe.ne.i .and. ifixe.ne.-2) then     ! fixed electron density or 2d density
            if(icoosys.eq.1) then
              do 27 idim=1,ndim
! note that ix can be negative or positive. nint is a better choice.
                ixo(idim)=nint(delxi(idim)*xold(idim,i))
  27            ixn(idim)=nint(delxi(idim)*xnew(idim,i))
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

   28 continue
      if(izigzag.gt.0) call zigzag2d(p,q,xold,xnew,0)
      if(ifixe.le.-2) call pairden2d(p,q,xold,xnew)  ! full pair-density
      if(ifourier.eq.1 .or. ifourier.eq.3) call fourierrk(p,q,xold,xnew)
      if(ifourier.eq.2 .or. ifourier.eq.3) call fourierkk(p,q,xold,xnew)

! accept new move with probability p

      if(rannyu(0).le.p) then

! move is accepted so update positions etc.
        pot_ee_old = pot_ee_new

        do 30 k=1,ndim
          do 30 i=1,nelec
            xold(k,i)=xnew(k,i)
   30       vold(k,i)=vnew(k,i)
        do 40 ifr=1,nforce
          eold(ifr)=enew(ifr)
   40     psi2o(ifr)=psi2n(ifr)
        accsum=accsum+one
        peo=pen
        peio=pein
        tjfo=tjfn
        psido=psidn
        psijo=psijn

!        if(nefp.gt.0) call efpsav

        call grad_hess_jas_save

      endif

      eloc = eold(1)
      call object_modified_by_index (eloc_index)
      call object_modified_by_index (xold_index)  !JT

      do 380 ifr=1,nforce
        if(ifr.eq.1) then
          esum1s(ifr)=eold(ifr)
         else
          wstro=exp(psi2o(ifr)-psi2o(1))
          wsum1s(ifr)=wstro
          esum1s(ifr)=eold(ifr)*wstro
        endif
  380 continue

      call acusig

      return
      end
