      subroutine dmc_good_ps_movall
! Written by Cyrus Umrigar and Claudia Filippi starting from dmc_good
! Uses the diffusion Monte Carlo algorithm described in:
! 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
!    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993).
! with portions related to nuclear cusps removed.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Control variables are:
! idmc         < 0     VMC
!              > 0     DMC
! abs(idmc)    = 1 **  simple kernel using dmc.brock.f
!              = 2     good kernel using dmc_good or dmc_good_inhom
! ipq         <= 0 *   do not use expected averages
!             >= 1     use expected averages (mostly in all-electron move algorithm)
! itau_eff    <=-1 *   always use tau in branching (not implemented)
!              = 0     use 0 / tau for acc /acc_int moves in branching
!             >= 1     use tau_eff (calcul in equilibration runs) for all moves
! iacc_rej    <=-1 **  accept all moves (except possibly node crossings)
!              = 0 **  use weights rather than accept/reject
!             >= 1     use accept/reject
! icross      <=-1 **  kill walkers that cross nodes (not implemented)
!              = 0     reject walkers that cross nodes
!             >= 1     allow walkers to cross nodes
!                      (OK since crossing prob. goes as tau^(3/2))
! icuspg      <= 0     approximate cusp in Green function
!             >= 1     impose correct cusp condition on Green function
! idiv_v      <= 0     do not use div_v to modify diffusion in good kernel
!              = 1     use homog. div_v to modify diffusion in good kernel
!             >= 2     use inhomog. div_v to modify diffusion in good kernel
! icut_br     <= 0     do not limit branching
!             >= 1 *   use smooth formulae to limit branching to (1/2,2)
!                      (bad because it makes energies depend on E_trial)
! icut_e      <= 0     do not limit energy
!             >= 1 *   use smooth formulae to limit energy (not implemented)

! *  => bad option, modest deterioration in efficiency or time-step error
! ** => very bad option, big deterioration in efficiency or time-step error
! So, idmc=6,66 correspond to the foll. two:
! 21 1 10 0 00 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
! 21 0 11 0 00 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
! Another reasonable choice is:
! 21 0 11 1 10 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      use constants_mod
      use atom_mod
      use dets_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use forcepar_mod
      use delocc_mod
      use force_dmc_mod
      use iterat_mod
      use jacobsave_mod
      use denupdn_mod
      use stepv_mod
      use config_dmc_mod
      use branch_mod
      use estsum_dmc_mod
      use estcum_dmc_mod
      use div_v_dmc_mod
      use contrldmc_mod
      use stats_mod
      use age_mod
      use pairden_mod
      use fourier_mod
      use pop_control_mod, only : ffn
      use distance_mod, only: pot_ee
      use config_mod, only: pot_ee_new, pot_ee_old
      use contrl_per_mod
      use zigzag_mod, only: izigzag
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-10, adrift=0.5d0)

      common /tmp/ eacc,enacc,macc,mnacc
      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys
      dimension xnc(3,nelec),xoc(3,nelec)
      dimension ixo(3),ixn(3)

      dimension xnew(3,nelec,nforce),vnew(3,nelec,nforce),psidn(nforce), &
     &psijn(nforce),enew(nforce),pen(nforce),pein(nforce),d2n(nforce),div_vn(nelec)
      dimension xbac(3),ajacnew(nforce)

      gauss()=dcos(two*pi*rannyu(0))*sqrt(-two*dlog(rannyu(0)))

!     term=(sqrt(two*pi*tau))**3/pi

      ajacnew(1)=one

! Undo products
      ipmod=mod(ipass,nfprod)
      ipmod2=mod(ipass+1,nfprod)
      if(idmc.gt.0) then
        ginv=min(1.d0,tau)
        ffn=eigv*(wdsumo/nconf_global)**ginv
        ffi=one/ffn
        fprod=fprod*ffn/ff(ipmod)
        ff(ipmod)=ffn
        fprodd=fprod/ff(ipmod2)
       else
        ginv=1
        ffn=1
        ffi=1
        fprod=1
        ff(ipmod)=1
        expon=1
        dwt=1
        fprodd=1
      endif

! Undo weights
      iwmod=mod(ipass,nwprod)

      ioldest=0
      do 300 iw=1,nwalk
! Loop over primary+secondary paths to compute forces
        do 300 ifr=1,nforce

! Stretch nuclei and electrons if calculating forces
          if(nforce.gt.1) then
            if(ifr.eq.1.or.istrech.eq.0) then
! Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
              call strech(xnew(1,1,1),xnew(1,1,ifr),ajacob,ifr,0)
              ajacnew(ifr)=one
             else
              call strech(xnew(1,1,1),xnew(1,1,ifr),ajacob,ifr,1)
              ajacnew(ifr)=ajacob
            endif
          endif

! Sample Green function for forward move
          dr2=zero
          dfus2o=zero
          vav2sumo=zero
          v2sumo=zero
          fnormo=one
          do 100 i=1,nelec
! Tau secondary in drift in order to compute tau_s in second equil. block
! Note: vavvo=vav2sumo/v2sumo appears in the branching
            if(ifr.gt.1.and.itausec.eq.1) then
              if(itau_eff.ge.1) then
                tauu=tau*taueff(ifr)/taueff(1)
               else
                tauu=taueff(ifr)
              endif
             else
              tauu=tau
            endif

! Use more accurate formula for the drift
            v2old=0
            do 5 k=1,ndim
    5         v2old=v2old+voldw(k,i,iw,ifr)**2
            volda=sqrt(v2old)
            vav=(dsqrt(1+2*adrift*v2old*tauu)-1)/(adrift*volda*tauu)
            vav2sumo=vav2sumo+vav**2
            v2sumo=v2sumo+v2old

! Calculate drifted position
            do 60 k=1,ndim
              xbac(k)=xoldw(k,i,iw,ifr)+voldw(k,i,iw,ifr)*tauu*vav/volda
              if(ifr.gt.1) then
                dfus=xnew(k,i,ifr)-xbac(k)
                dfus2o=dfus2o+dfus**2
              endif
   60       continue

! Do the diffusion.  Actually, this diffusion contains part of the drift too
! if idiv_v.ge.1
! Use div_v_hom to set tau_hom
            div_v_hom=div_vow(i,iw)/3

            if(ifr.eq.1) then
              if(idiv_v.ge.1) then
                if(div_v_hom.lt.0.d0) then
                  tau_hom=3*tau-1/div_v_hom &
     &            -((tau-1/div_v_hom)*sqrt(-div_v_hom) &
     &            *(1-derfc(1/sqrt(-2*div_v_hom*tau))) &
     &            +sqrt(2*tau/pi)*exp(1/(2*div_v_hom*tau)))**2
                 else
                  tau_hom=tau*(1+2*div_v_hom*tau)/(1+div_v_hom*tau)
                endif
               else
                tau_hom=tau
              endif
              rttau_hom=dsqrt(tau_hom)

! Sample gaussian
              dfus2a=0
              do 80 k=1,ndim
                drift=xbac(k)-xoldw(k,i,iw,ifr)
                dfus=gauss()*rttau_hom
                dx=drift+dfus
                dr2=dr2+dx**2
                dfus2a=dfus2a+dfus**2
                dfus2o=dfus2o+dfus**2
   80           xnew(k,i,ifr)=xbac(k)+dfus
            endif
            gaus_norm=1/sqrt(two*pi*tau_hom)**3
  100       fnormo=fnormo*gaus_norm*exp(-half*dfus2a/tau_hom)
          vavvo=sqrt(vav2sumo/v2sumo)

! calculate psi etc. at new configuration
          call hpsi(xnew(1,1,ifr),psidn(ifr),psijn(ifr),vnew(1,1,ifr),div_vn,d2n(ifr),pen(ifr),pein(ifr),enew(ifr),denergy,ifr)

! Check for node crossings
          if(psidn(ifr)*psidow(iw,ifr).le.zero.and.ifr.eq.1) then
            nodecr=nodecr+1
            if(icross.le.0) then
              p=zero
              goto 210
            endif
          endif

! Calculate Green function for the reverse move
          dfus2n=zero
          vav2sumn=zero
          v2sumn=zero
          fnormn=one
          do 200 i=1,nelec
            v2new=0
            do 151 k=1,ndim
  151         v2new=v2new+vnew(k,i,ifr)**2
            vnewa=sqrt(v2new)
            vav=(dsqrt(1+2*adrift*v2new*tauu)-1)/(adrift*vnewa*tauu)
            vav2sumn=vav2sumn+vav**2
            v2sumn=v2sumn+v2new

! Calculate drifted position
            dfus2a=0
            do 160 k=1,ndim
              xbac(k)=xnew(k,i,ifr)+vnew(k,i,ifr)*tauu*vav/vnewa
  160         dfus2a=dfus2a+(xoldw(k,i,iw,ifr)-xbac(k))**2

! Use div_v_hom to set tau_hom
            div_v_hom=div_vn(i)/3

            if(idiv_v.ge.1) then
              if(div_v_hom.lt.0.d0) then
                tau_hom=3*tau-1/div_v_hom &
     &          -((tau-1/div_v_hom)*sqrt(-div_v_hom) &
     &          *(1-derfc(1/sqrt(-2*div_v_hom*tau))) &
     &          +sqrt(2*tau/pi)*exp(1/(2*div_v_hom*tau)))**2
               else
                tau_hom=tau*(1+2*div_v_hom*tau)/(1+div_v_hom*tau)
              endif
             else
              tau_hom=tau
            endif

            gaus_norm=1/sqrt(two*pi*tau_hom)**3
            fnormn=fnormn*gaus_norm*exp(-half*dfus2a/tau_hom)

            if(ipr.ge.1) then
              write(6,'(''xold'',9f10.6)')(xoldw(k,i,iw,ifr),k=1,ndim), &
     &        (xnew(k,i,ifr),k=1,ndim), (xbac(k),k=1,ndim)
              write(6,'(''iel,fnormn,fnormn/fnormo'',i3,f10.6,d12.4)')i,fnormn,fnormn/fnormo
            endif
  200     continue
          vavvn=sqrt(vav2sumn/v2sumn)

          p=(psidn(ifr)/psidow(iw,ifr))**2*exp(2*(psijn(ifr)- &
     &    psijow(iw,ifr)))*fnormn/fnormo

          if(ipr.ge.1) write(6,'(''p'',9f10.6)') &
     &    p,(psidn(ifr)/psidow(iw,ifr))**2*exp(2*(psijn(ifr)- &
     &    psijow(iw,ifr))),exp((dfus2o-dfus2n)/(two*tau)),psidn(ifr), &
     &    psidow(iw,ifr),psijn(ifr),psijow(iw,ifr),dfus2o,dfus2n

! The following is one reasonable way to cure persistent configurations
! Not needed if itau_eff <=0 and in practice we have never needed it even
! otherwise
          if(iage(iw).gt.50) p=p*1.1d0**(iage(iw)-50)

          pp=p
          p=dmin1(one,p)
  210     q=one-p

          dfus2unf(ifr)=dfus2unf(ifr)+dfus2o

! Set taunow for branching.  Note that for primary walk taueff(1) always
! is tau*dfus2ac/dfus2unf so if itau_eff.le.0 then we reset taunow to tau
! On the other hand, secondary walk only has dfus2ac/dfus2unf if
! itau_eff.ge.1 so there is no need to reset taunow.
          if(ifr.eq.1) then
            acc=acc+p
            try_int=try_int+1
            dfus2ac=dfus2ac+p*dfus2o
            dr2ac=dr2ac+p*dr2
            dr2un=dr2un+dr2
            tautot=tautot+tau*dfus2ac/dfus2unf(1)

! If we are using weights rather than accept/reject
            if(iacc_rej.le.0) then
              p=one
              q=zero
            endif

            taunow=taueff(ifr)
            if(rannyu(0).lt.p) then
              iaccept=1
              acc_int=acc_int+1
              if(itau_eff.le.0) taunow=tau
              if(ipq.le.0) p=one
              eacc=eacc+enew(ifr)
              macc=macc+1
             else
              iaccept=0
              if(itau_eff.le.0) taunow=zero
              if(ipq.le.0) p=zero
              enacc=enacc+eoldw(iw,ifr)
              mnacc=mnacc+1
            endif
            q=one-p

            psav=p
            qsav=q

           elseif(itausec.eq.1) then

            taunow=taueff(ifr)
            if(itau_eff.le.0.and.iaccept.eq.0) taunow=zero

          endif
          if(ipr.ge.1)write(6,'(''wt'',9f10.5)') wt(iw),etrial,eest

! Warning: Change UNR93 reweighting factor because it gives large time-step error at small tau for pseudo systems
!         ewto=eest-(eest-eoldw(iw,ifr))*vavvo                                                        ! UNR93
!         ewtn=eest-(eest   -enew(ifr))*vavvn                                                         ! UNR93
!         ewto=eoldw(iw,ifr)                                                                          ! no_ene_int
!         ewtn=enew(ifr)                                                                              ! no_ene_int
!         ewto=eest-(eest-eoldw(iw,ifr))*vavvo**2*(3-2*vavvo)                                         ! new_ene_int
!         ewtn=eest-(eest-enew(ifr))*vavvn**2*(3-2*vavvn)                                             ! new_ene_int
!         ewto=eest-(eest-eoldw(iw,ifr))*vavvo**3*(10-15*vavvo+6*vavvo**2)                            ! new_ene_int2
!         ewtn=eest-(eest-enew(ifr))*vavvn**3*(10-15*vavvn+6*vavvn**2)                                ! new_ene_int2
          ecut=0.2*sqrt(nelec/tau)                                                                    ! Alfe
          ewto=etrial+min(ecut,max(-ecut,eoldw(iw,ifr)-eest))                                         ! Alfe
          ewtn=etrial+min(ecut,max(-ecut,enew(ifr)-eest))                                             ! Alfe
! Note in the 2 lines above we use etrial instead of eest for the first term because we do population control with a fixed etrial
! so that we keep track of the fluctuating factor and undo the population control

          if(idmc.gt.0) then
            expon=(etrial-half*((one+qsav)*ewto+psav*ewtn))*taunow
            if(icut_br.le.0) then
              dwt=dexp(expon)
             else
              dwt=0.5d0+1/(1+exp(-4*expon))
!             if(expon.gt.0) then
!               dwt=(1+2*expon)/(1+expon)
!              else
!               dwt=(1-expon)/(1-2*expon)
!             endif
            endif
          endif

! If we are using weights rather than accept/reject
          if(iacc_rej.eq.0) dwt=dwt*pp

! Exercise population control if dmc or vmc with weights
          if(idmc.gt.0.or.iacc_rej.eq.0) dwt=dwt*ffi

! Set weights and product of weights over last nwprod steps
          pwt(iw,ifr)=pwt(iw,ifr)*dwt/wthist(iw,iwmod,ifr)
          wthist(iw,iwmod,ifr)=dwt
          if(ifr.eq.1) then
            wt(iw)=wt(iw)*dwt
            wtnow=wt(iw)
           else
            wtnow=wt(iw)*pwt(iw,ifr)/pwt(iw,1)
          endif

          if(ipr.ge.1)write(6,'(''eold,enew,wt'',9f10.5)') &
     &    eoldw(iw,ifr),enew(ifr),wtnow

          wtg=wtnow*fprod

          enowo=eoldw(iw,ifr)
          enown=enew(ifr)

          ! used to plot interaction potential (ACM)
          pot_ee_new = pot_ee ! array assignment
          pot_ee_old = pot_ee_oldw(:,iw,ifr) ! array assignment

          if(ifr.eq.1) then
            psi2savo=2*(psijow(iw,1)+dlog(dabs(psidow(iw,1))))
            psi2savn=2*(psijn(1)+dlog(dabs(psidn(1))))

            wsum1(ifr)=wsum1(ifr)+wtnow
            esum1(ifr)=esum1(ifr)+wtnow*(q*enowo+p*enown)
            pesum(ifr)=pesum(ifr)+  wtg*(q*peow(iw,ifr) +p*pen(ifr))
            peisum(ifr)=peisum(ifr)+  wtg*(q*peiow(iw,ifr) +p*pein(ifr))
            tpbsum(ifr)=tpbsum(ifr)+wtg*(q*(enowo-peow(iw,ifr))+ &
     &                                   p*(enown-pen(ifr)))
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*(q*d2ow(iw,ifr)+ &
     &                                           p*d2n(ifr))

           else

            ro=ajacold(iw,ifr)*psidow(iw,ifr)**2* &
     &         exp(2*psijow(iw,ifr)-psi2savo)
            rn=ajacnew(ifr)*psidn(ifr)**2* &
     &         exp(2*psijn(ifr)-psi2savn)

            wsum1(ifr)=wsum1(ifr)+wtnow*(qsav*ro+psav*rn)
            esum1(ifr)=esum1(ifr)+wtnow*(qsav*enowo*ro+psav*enown*rn)
            pesum(ifr)=pesum(ifr)+  wtg*(qsav*peow(iw,ifr)*ro+ &
     &                                   psav*pen(ifr)*rn)
            peisum(ifr)=peisum(ifr)+  wtg*(qsav*peiow(iw,ifr)*ro+ &
     &                                   psav*pein(ifr)*rn)
            tpbsum(ifr)=tpbsum(ifr)+wtg*(qsav*(enowo-peow(iw,ifr))*ro+ &
     &                                   psav*(enown-pen(ifr))*rn)
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*(qsav*d2ow(iw,ifr)*ro+ &
     &                                           psav*d2n(ifr)*rn)
          endif

! Collect density only for primary walk
          if(ifr.eq.1) then
            wtgp=wtg*p
            wtgq=wtg*q
            do 250 i=1,nelec
              r2o=zero
              r2n=zero
              do 240 k=1,ndim
                r2n=r2n+xnew(k,i,ifr)**2
  240           r2o=r2o+xoldw(k,i,iw,ifr)**2
              rmino=zero
              rminn=zero
              do 245 k=1,ndim
                rmino=rmino+(xoldw(k,i,iw,ifr)-cent(k,1))**2
  245           rminn=rminn+(xnew(k,i,ifr)-cent(k,1))**2
              rmino=sqrt(rmino)
              rminn=sqrt(rminn)
!             itryo=min(int(delri*rmino)+1,NRAD)
!             itryn=min(int(delri*rminn)+1,NRAD)
              itryo=int(min(delri*rmino+1,dfloat(NRAD))+eps)
              itryn=int(min(delri*rminn+1,dfloat(NRAD))+eps)
              if(i.le.nup) then
                rprobup(itryo)=rprobup(itryo)+wtgq
                rprobup(itryn)=rprobup(itryn)+wtgp
               else
                rprobdn(itryo)=rprobdn(itryo)+wtgq
                rprobdn(itryn)=rprobdn(itryn)+wtgp
              endif
              rprob(itryo)=rprob(itryo)+wtgq
              rprob(itryn)=rprob(itryn)+wtgp
              r1sum=r1sum+wtg*(q*dsqrt(r2o)+p*dsqrt(r2n))
              r2sum=r2sum+wtg*(q*r2o+p*r2n)
              r3sum=r3sum+wtg*(q*r2o*dsqrt(r2o)+p*r2n*dsqrt(r2n))
              r4sum=r4sum+wtg*(q*r2o*r2o+p*r2n*r2n)
              risum=risum+wtg*(q/dsqrt(r2o)+p/dsqrt(r2n))

! calculate 2d density related functions:
              if(iperiodic.eq.1 .or. iperiodic.eq.3) then  ! In 1D and 3D reduce position to simulation cell
                call reduce_sim_cell(xoldw(:,i,iw,ifr))
                call reduce_sim_cell(xnew(:,i,ifr))
              endif
              if(izigzag.gt.0) call zigzag2d(wtgp,wtgq,xoldw(:,:,iw,ifr),xnew(:,:,ifr), i)
              if(ifixe.eq.-1 .or. ifixe.eq.-3) then

                if(icoosys.eq.1) then
                  do 247 idim=1,ndim
! note that ix can be negative or positive. nint is a better choice.
                    ixo(idim)=nint(delxi(idim)*xoldw(idim,i,iw,ifr))
  247               ixn(idim)=nint(delxi(idim)*xnew(idim,i,ifr))
                  else
! same trick adapted to circular coordinates
                    ixo(1)=nint(delradi*(rmino-rmean))
                    ixn(1)=nint(delradi*(rminn-rmean))
                    ixo(2)=nint(delti*(datan2(xoldw(2,i,iw,ifr),xoldw(1,i,iw,ifr))))
                    ixn(2)=nint(delti*(datan2(xnew(2,i,ifr),xnew(1,i,ifr))))
                endif

                if(abs(ixo(1)).le.NAX .and. abs(ixo(2)).le.NAX) then
                  den2d_t(ixo(1),ixo(2))=den2d_t(ixo(1),ixo(2))+wtgq
                  pot_ee2d_t(ixo(1),ixo(2))=pot_ee2d_t(ixo(1),ixo(2))+wtgq*pot_ee_old(i)
                  if(i.le.nup) then
                    den2d_u(ixo(1),ixo(2))=den2d_u(ixo(1),ixo(2))+wtgq
                    pot_ee2d_u(ixo(1),ixo(2))=pot_ee2d_u(ixo(1),ixo(2))+wtgq*pot_ee_old(i)
                  else
                    den2d_d(ixo(1),ixo(2))=den2d_d(ixo(1),ixo(2))+wtgq
                    pot_ee2d_d(ixo(1),ixo(2))=pot_ee2d_d(ixo(1),ixo(2))+wtgq*pot_ee_old(i)
                  endif
                endif
                if(abs(ixn(1)).le.NAX .and. abs(ixn(2)).le.NAX) then
                  den2d_t(ixn(1),ixn(2))=den2d_t(ixn(1),ixn(2))+wtgp
                  pot_ee2d_t(ixn(1),ixn(2))=pot_ee2d_t(ixn(1),ixn(2))+wtgp*pot_ee_new(i)
                  if(i.le.nup) then
                    den2d_u(ixn(1),ixn(2))=den2d_u(ixn(1),ixn(2))+wtgp
                    pot_ee2d_u(ixn(1),ixn(2))=pot_ee2d_u(ixn(1),ixn(2))+wtgp*pot_ee_new(i)
                  else
                    den2d_d(ixn(1),ixn(2))=den2d_d(ixn(1),ixn(2))+wtgp
                    pot_ee2d_d(ixn(1),ixn(2))=pot_ee2d_d(ixn(1),ixn(2))+wtgp*pot_ee_new(i)
                  endif
                endif
              endif

  250       continue

            if(ifixe.le.-2 .or. ifourier.ne.0) then     ! full pair-density
              do 255 k=1,ndim
                do 255 i=1,nelec
                  xnc(k,i)=xnew(k,i,ifr)
  255             xoc(k,i)=xoldw(k,i,iw,ifr)
              if(ifixe.le.-2) call pairden2d(wtgp,wtgq,xoc,xnc)
              if(ifourier.eq.1 .or. ifourier.eq.3) call fourierrk(wtgp,wtgq,xoc,xnc)
              if(ifourier.eq.2 .or. ifourier.eq.3) call fourierkk(wtgp,wtgq,xoc,xnc)
!             if(ifourier.eq.1) call fourier2d(wtgp,wtgq,xoc,xnc)
            endif

          endif

          if(iaccept.eq.1) then

!           if(dabs((enew(ifr)-etrial)/etrial).gt.0.2d+0) then
            if(dabs((enew(ifr)-etrial)/etrial).gt.5.0d+0) then
               write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass, &
     &         enew(ifr)-etrial,psidn(ifr),psijn(ifr), &
     &         ((xnew(k,jj,ifr),k=1,ndim),jj=1,nelec)
            endif

            if(wt(iw).gt.3) write(18,'(i6,i4,3f8.2,30f8.4)') ipass,iw, &
     &      wt(iw),enew(ifr)-etrial,eoldw(iw,ifr)-etrial, &
     &      ((xnew(k,jj,ifr),k=1,ndim),jj=1,nelec)

            psidow(iw,ifr)=psidn(ifr)
            psijow(iw,ifr)=psijn(ifr)
            eoldw(iw,ifr)=enew(ifr)
            peow(iw,ifr)=pen(ifr)
            peiow(iw,ifr)=pein(ifr)
            d2ow(iw,ifr)=d2n(ifr)
            pot_ee_oldw(:,iw,ifr) = pot_ee_new ! array assignment
            do 260 i=1,nelec
              div_vow(i,iw)=div_vn(i)
              do 260 k=1,ndim
                xoldw(k,i,iw,ifr)=xnew(k,i,ifr)
  260           voldw(k,i,iw,ifr)=vnew(k,i,ifr)
            iage(iw)=0
            ajacold(iw,ifr)=ajacnew(ifr)
           else
            if(ifr.eq.1) then
              iage(iw)=iage(iw)+1
              ioldest=max(ioldest,iage(iw))
              ioldestmx=max(ioldestmx,iage(iw))
            endif

          endif

  300 continue

      if(wsum1(1).gt.1.1*nconf_global) write(18,'(i6,9d12.4)') ipass,ffn,fprod, &
     &fprodd,wsum1(1),wgdsumo

      wdsumn=wsum1(1)
      wdsum1=wdsumo
!     wgdsum1=wgdsumo
      if(idmc.gt.0.or.iacc_rej.eq.0) then
        wfsum1=wsum1(1)*ffn
        wgdsumn=wsum1(1)*fprodd
        efsum1=esum1(1)*ffn
       else
        wfsum1=wsum1(1)
        wgdsumn=wsum1(1)
        efsum1=esum1(1)
      endif

      wsum=wsum+wsum1(1)
      wfsum=wfsum+wfsum1
      wdsum=wdsum+wdsumo
      wgdsum=wgdsum+wgdsumo
      esum=esum+esum1(1)
      efsum=efsum+efsum1
      eisum=eisum+wfsum1/wdsumo

      do 305 ifr=1,nforce
        if(idmc.gt.0.or.iacc_rej.eq.0) then
          wgsum1(ifr)=wsum1(ifr)*fprod
          egsum1(ifr)=esum1(ifr)*fprod
         else
          wgsum1(ifr)=wsum1(ifr)
          egsum1(ifr)=esum1(ifr)
        endif

        wgsum(ifr)=wgsum(ifr)+wgsum1(ifr)
  305   egsum(ifr)=egsum(ifr)+egsum1(ifr)

! Do split-join
!JT      call splitj ! moved outside the routine

! Estimate eigenvalue of G from the energy
      eest=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
      if(itau_eff.ge.1) then
        eigv=dexp((etrial-eest)*taueff(1))
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest, &
     &  egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
       else
        accavn=acc_int/try_int
        if(icut_br.le.0) then
          eigv=accavn*exp((etrial-eest)*tau)+(one-accavn)
         else
          eigv=one-accavn+accavn*(0.5d0+1/(1+exp(-4*(etrial-eest)*tau)))
        endif
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest,accavn, &
     &  egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wdsumn
      wgdsumo=wgdsumn
      wtgen(ipmod)=wdsumn

      return
      end
