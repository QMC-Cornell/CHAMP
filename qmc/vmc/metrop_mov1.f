      subroutine metrop_mov1(ipass)
c     subroutine metrop(ipass)
c Written by Cyrus Umrigar and Claudia Filippi
c routine to move configuration by one step using a
c velocity bias type Monte Carlo.
c This routine has been changed to move one electron at a time and
c to keep track of the acceptance prob. as a function of delta
c Minor mods by A.D.Guclu to include pair-density function calculation
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
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-10)

      logical vgreater

c Warning: program has a bug, but it does not matter much since
c we never use it.

c     common /contr3/ mode
      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys

      dimension dx(3,nelec),idist(nelec),ixo(3),ixn(3)
      dimension xstrech(3,nelec)

c     write(6,'(''entering metrop_mov1'')')

c     mode='vmc_mov1    '

c Sample transition probability from current state to new state
c and record value of probability in fxop.
c The transition probabilty is given by statement label 17
c for the reverse transition.
c The transition probability is an approximately psi(new)/psi(old)
c and is positive definite.

      deltfi=two*fbias*deltai
      do 300 i=1,nelec
        if(i.ne.ifixe) then

          fxop=one
          do 10 k=1,ndim
            fx=dmin1(dabs(vold(k,i)),deltfi)
            ftest=fx*delta*half
            fx=sign(fx,vold(k,i))
            if(ftest.le.rannyu(0)) then
c           sample from 1 piece
              dx(k,i)=(rannyu(0)-half)*delta
             else
c             sample from the (xnew-xold) piece assuming velocity positive
              dx(k,i)=(dmax1(rannyu(0),rannyu(0))-half)*delta
c             if velocity was negative then dx has opposite sign
              if(vold(k,i).lt.zero) dx(k,i)=-dx(k,i)
            endif
            fxop=fxop*(one+fx*dx(k,i))
   10       xnew(k,i)=xold(k,i)+dx(k,i)

c calculate psi etc. at new configuration
          iel=i
c         call hpsi(xnew,psidn,psijn,vnew,div_vn,d2,pen,enew(1),denergy,1)
c         write(6,'(''calling hpsie from metrop_mov1'')')
          call hpsie(iel,xnew,psidn,psijn,vnew)
          psi2n(1)=2*(dlog(dabs(psidn))+psijn)

c save electrostatic potential at new configuration
c   Note that we already called distances.f from hpsi, otherwise
c   we should add a call to distances here (ACM)

          pot_ee_new = pot_ee  ! this is an array assignment

c If error is large then save config. to use in optimizing routine

c         if(dabs((enew(1)-etrial)/etrial).gt.1.d-3) then
cc        if(dabs((enew(1)-etrial)/etrial).gt.1.d+1) then
c         write(8,'(i6,f8.2,2d10.2,(8f8.4))') ipass,enew(1)-etrial,psin,
c    &    (enew(1)-etrial)*psin,((xnew(k,jj),k=1,ndim),jj=1,nelec)
c         endif
c         if(dabs(psin*(enew(1)-etrial)/etrial).gt.1.d-7) then
c           write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass,enew(1)-etrial,psin,
c    &      (enew(1)-etrial)*psin,((xnew(k,jj),k=1,ndim),jj=1,nelec)
c         endif

c calculate probability for reverse transition
          fxnp=one
          do 20 k=1,ndim
   20       fxnp=fxnp*(one-dx(k,i)
     &      *sign(dmin1(dabs(vnew(k,i)),deltfi),vnew(k,i)))

c form the Jackson Feenberg kinetic energy expression

c         tjfn=d2
c         tjfn=-tjfn*half*hb

c p is the probability of accepting new move

          p=exp(psi2n(1)-psi2o(1))*fxnp/fxop
cr        write(6,'(''psi2n,psi2o,fxnp,fxop,psi2n/psi2o,fxnp/fxop,p'',9f12.6
cr   &    )')psi2n(1),psi2o(1),fxnp,fxop,psi2n(1)/psi2o(1),fxnp/fxop,p
          p=dmin1(one,p)
          q=one-p

c form expected values of e, pe, etc.

c     esum1=          p*enew(1)+q*eold(1)
c     esum(1)=esum(1)+p*enew(1)+q*eold(1)
c     pesum=pesum+p*pen+q*peo
c     tpbsum=tpbsum+p*(enew(1)-pen)+q*(eold(1)-peo)
c     tjfsum=tjfsum+p*tjfn+q*tjfo
c     do 25 k=1,ndim
c  25     r2sum=r2sum+p*xnew(k,i)**2+q*xold(k,i)**2
c     if(nefp.gt.0) then
c       call sample_efp(0,xold,eold(1),q)
c       call sample_efp(1,xnew,enew(1),p)
c     endif
c     call acues1

c Calculate as a function of the distance to the nucleus
c 1) acceptance,  2) force-bias truncation probability,
c 3) kinetic energy and it's fluctuation
c The K.E. is not quite correct, since we should use p times new
c and q times old, and keep track of which bin the old was in
c The reason why I changed
c itryo=min(int(delri*rold)+1,NRAD)  to
c itryo=int(min(delri*rold+1,dfloat(NRAD))+eps)
c is that 2147483647 is the largest 32-bit integer and 1 more than that gives -2147483648.
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
c this was for ndim=3
c      if((dabs(vold(1,i)).gt.deltfi).or.(dabs(vold(2,i)).gt.deltfi)
c     &.or.(dabs(vold(3,i)).gt.deltfi)) trunfb(itryo)=trunfb(itryo)+1
c for general ndim
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

c calculate 2d-density related functions
          if(iperiodic.eq.1) then  ! 1D periodic bc's, so make sure x-posn between -a/2 and a/2
            call reduce_sim_cell(xold(:,i))
            call reduce_sim_cell(xnew(:,i))
          endif
          if(ifixe.le.-2) call pairden2d(p,q,xold,xnew)
          if(ifourier.eq.1 .or. ifourier.eq.3) call fourierrk(p,q,xold,xnew)
          if(ifourier.eq.2 .or. ifourier.eq.3) call fourierkk(p,q,xold,xnew)
c          if(ifourier.eq.1) call fourier2d(1.d0,0.d0,xold,xnew)
          if(ifixe.ne.0 .and. ifixe.ne.-2) then   !  2d density (ifixe=-1,-3) or fixed electron pair density (ifixe>0)
            if(icoosys.eq.1) then 
              do 170 idim=1,ndim
c note that ix can be negative or positive. nint is a better choice.
                ixo(idim)=nint(delxi*xold(idim,i))
  170           ixn(idim)=nint(delxi*xnew(idim,i))
             else
c same trick adapted to circular coordinates
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

          do 210 k=1,ndim
  210       r2sum=r2sum+p*xnew(k,i)**2+q*xold(k,i)**2

c     eksum=zero
c     do 26 ii=1,nelec
c  26   eksum=eksum+ekine(ii)
c     write(6,'(''ke='',9d13.5)') eksum,enew(1)-pen


c accept new move with probability p
c Note when one electron moves the velocity on all electrons change.
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

c loop over secondary configurations
      do 350 ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
c     write(6,'(''calling hpsi from metrop_mov1,ifr,nforce'',9i5)') ifr,nforce
        call hpsi(xstrech,psido,psijo,vold,div_vo,d2,peo,peio,eold(ifr),denergy,ifr)
  350   psi2o(ifr)=2*(dlog(dabs(psido))+psijo)+dlog(ajacob)

c primary configuration
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
c form expected values of e, pe, etc.
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
          wstro=exp(psi2o(ifr)-psi2o(1))
          wsum1s(ifr)=wstro
          esum(ifr)=esum(ifr)+eold(ifr)*wstro
          wsum(ifr)=wsum(ifr)+wstro
          esum1s(ifr)=eold(ifr)*wstro
c psi2o has the log of the wavefn.
c         if(abs(psi2o(ifr)-psi2o(1)).gt.50.d0)
c    &    write(6,'(''ifr,esum(ifr),wstro'',i3,100d12.4)')
c    &    ifr,esum(ifr),wstro,psi2o(ifr),psi2o(1)
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
