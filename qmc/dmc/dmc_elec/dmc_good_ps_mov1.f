      subroutine dmc_good_ps_mov1
c Written by Cyrus Umrigar and Claudia Filippi
c Uses the diffusion Monte Carlo algorithm described in:
c 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
c    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993)
c modified to do accept/reject after single-electron moves and to
c remove portions related to nuclear cusps.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Control variables are:
c idmc         < 0     VMC
c              > 0     DMC
c abs(idmc)    = 1 **  simple kernel using dmc.brock.f
c              = 2     good kernel using dmc_good or dmc_good_inhom
c ipq         <= 0 *   do not use expected averages
c             >= 1     use expected averages (mostly in all-electron move algorithm)
c itau_eff    <=-1 *   always use tau in branching (not implemented)
c              = 0     use 0 / tau for acc /acc_int moves in branching
c             >= 1     use tau_eff (calcul in equilibration runs) for all moves
c iacc_rej    <=-1 **  accept all moves (except possibly node crossings)
c              = 0 **  use weights rather than accept/reject
c             >= 1     use accept/reject
c icross      <=-1 **  kill walkers that cross nodes (not implemented)
c              = 0     reject walkers that cross nodes
c             >= 1     allow walkers to cross nodes
c                      (OK since crossing prob. goes as tau^(3/2))
c icuspg      <= 0     approximate cusp in Green function
c             >= 1     impose correct cusp condition on Green function
c idiv_v      <= 0     do not use div_v to modify diffusion in good kernel
c              = 1     use homog. div_v to modify diffusion in good kernel
c             >= 2     use inhomog. div_v to modify diffusion in good kernel
c icut_br     <= 0     do not limit branching
c              = 1 *   use first 2 terms in Taylor expansion of Exp if exponent>0
c             >= 2 *   use smooth formulae to limit branching to (1/2,2)
c                      (bad because it makes energies depend on E_trial)
c icut_e      <= 0     do not limit energy
c             >= 1 *   use smooth formulae to limit energy (not implemented)

c *  => bad option, modest deterioration in efficiency or time-step error
c ** => very bad option, big deterioration in efficiency or time-step error
c So, idmc=6,66 correspond to the foll. two:
c 2 1 1 1 0 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
c 2 1 0 1 1 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
c Another reasonable choice is:
c 2 1 0 1 1 1 1 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      use all_tools_mod
      use constants_mod
      use control_mod
      use average_mod
      use atom_mod
      use dets_mod
      use jastrow_mod, only: psi_jas
      use basis1_mod
      use contrl_mod
      use contr3_mod, only : mode
      use const_mod
      use dim_mod
      use forcepar_mod
      use contrl_per_mod
      use contrl_opt_mod
      use delocc_mod
      use force_dmc_mod
      use iterat_mod
      use jacobsave_mod
      use forcest_dmc_mod
      use denupdn_mod
      use stepv_mod
      use config_dmc_mod
      use branch_mod
      use estsum_dmc_mod
      use estcum_dmc_mod
      use estcm2_mod
      use div_v_dmc_mod
      use contrldmc_mod
      use stats_mod
      use age_mod
      use branch_dmc_opt_mod
      use pairden_mod
      use fourier_mod
      use velratio_mod
      use pop_control_mod, only : ffn
      use eloc_mod
!     use optimization_mod, only : l_opt_ovlp_fn
      use opt_ovlp_fn_mod, only : wt_lambda
      use distance_mod, only: pot_ee
      use config_mod, only: pot_ee_new, pot_ee_old
      use zigzag_mod, only: izigzag
      implicit real*8(a-h,o-z)

      parameter (adrift=0.5d0)

      common /circularmesh/ rmin,rmax,rmean,delradi,delti,nmeshr,nmesht,icoosys

      dimension xstrech(3,nelec)
      dimension xnew(3),vnew(3,nelec)
      dimension xbac(3)
      dimension itryo(nelec),itryn(nelec),unacp(nelec)
      dimension xnc(3,nelec),xoc(3,nelec),xnci(3,nelec,nelec),xoci(3,nelec,nelec)
      dimension ixo(3),ixn(3)
      dimension dewto(nparm),dewtn(nparm),dexponent(nparm)

      data ncall,ipr_sav /0,0/
      save ipr_sav

c     gauss()=dcos(two*pi*rannyu(0))*sqrt(-two*dlog(rannyu(0)))
      sigma(x,x2,w)=sqrt(max((x2/w-(x/w)**2)*nconf,0.d0))

c     rn_eff(w,w2)=w**2/w2
c     error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
c     errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

      if(idiv_v.ge.1) stop 'div_v not implemented yet in 1-electron move algorithm'

c     term=(sqrt(two*pi*tau))**3/pi

c Temporarily create wt_lambda_tau here, but it should be done at start of program
      wt_lambda_tau=wt_lambda**tau

c Undo products
      ipmod=mod(ipass,nfprod)
      ipmod2=mod(ipass+1,nfprod)
      if(idmc.gt.0) then
c     if(idmc.gt.0 .and. .not. l_opt_ovlp_fn) then
        ginv=min(1.d0,tau)
c       ginv=min(1.d0,10*tau)
        ffn=eigv*(wdsumo/nconf_global)**ginv
        ffi=one/ffn
        fprod=fprod*ffn/ff(ipmod)
        ff(ipmod)=ffn
       else
        ginv=1
        ffn=1
        ffi=1
        fprod=1
        ff(ipmod)=1
        expon=1
        dwt=1
      endif

c Undo weights
      iwmod=mod(ipass,nwprod)

c Store (well behaved velocity/velocity)
      if(ncall.eq.0.and.irstar.eq.0) then
        do 7 iw=1,nwalk
          do 7 ifr=1,nforce

            vav2sumo=zero
            v2sumo=zero
            do 6 i=1,nelec

c Tau secondary in drift is one (first time around)
              tratio=one

              v2old=0
              do 5 k=1,ndim
    5           v2old=v2old+voldw(k,i,iw,ifr)**2
              vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/
     &            (adrift*v2old)
              vavvo=vavvt/(tau*tratio)
              vav2sumo=vav2sumo+vavvo*vavvo*v2old
              v2sumo=v2sumo+v2old

    6         continue
            fratio(iw,ifr)=dsqrt(vav2sumo/v2sumo)
    7     continue
        ncall=ncall+1
      endif

      ioldest=0
      do 300 iw=1,nwalk
c Loop over primary walker
        current_walker = iw !JT
        call object_modified_by_index (current_walker_index) !JT

c Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
        if(nforce.gt.1)
     &  call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)

        if(ibasis.eq.3) then          ! complex basis set
          call cwalkstrdet(iw)
        else
          call walkstrdet(iw)
        endif
        call walkstrjas(iw)

c Sample Green function for forward move
        r1sume=zero
        r2sume=zero
        risume=zero
        dfus2ac=zero
        dfus2un=zero
        dr2ac=zero
        dr2un=zero
        drifdif=zero
        iaccept=0
        do 200 i=1,nelec
c Use more accurate formula for the drift
          v2old=0
          do 75 k=1,ndim
   75       v2old=v2old+voldw(k,i,iw,1)**2
c Tau primary -> tratio=one
          vavvt=(dsqrt(one+two*adrift*v2old*tau)-one)/(adrift*v2old)

          dr2=zero
          dfus2o=zero
          do 80 k=1,ndim
            drift=vavvt*voldw(k,i,iw,1)
            dfus=gauss()*rttau
            dx=drift+dfus
            dr2=dr2+dx**2
            dfus2o=dfus2o+dfus**2
   80       xnew(k)=xoldw(k,i,iw,1)+dx

          if(ipr.ge.1)
     &    write(6,'(''xnewdr'',2i4,9f8.5)') iw,i,(xnew(k),k=1,ndim)
c calculate psi and velocity at new configuration
          call hpsiedmc(i,iw,xnew,psidn,psijn,vnew)

c Check for node crossings
          if(psidn*psidow(iw,1).le.zero) then
            nodecr=nodecr+1
            if(icross.le.0) then
              p=zero
              goto 160
            endif
          endif

c Calculate Green function for the reverse move

          v2new=0
            do 149 k=1,ndim
  149         v2new=v2new+vnew(k,i)**2
          vavvt=(dsqrt(one+two*adrift*v2new*tau)-one)/(adrift*v2new)


          dfus2n=zero
          do 150 k=1,ndim
            drift=vavvt*vnew(k,i)
            xbac(k)=xnew(k)+drift
            dfus=xbac(k)-xoldw(k,i,iw,1)
  150       dfus2n=dfus2n+dfus**2

          if(ipr.ge.1) then
            write(6,'(''xoldw'',9f10.6)')(xoldw(k,i,iw,1),k=1,ndim),
     &      (xnew(k),k=1,ndim), (xbac(k),k=1,ndim)
            write(6,'(''dfus2o'',9f10.6)')dfus2o,dfus2n,
     &      psidow(iw,1),psidn,psijow(iw,1),psijn
          endif

          p=(psidn/psidow(iw,1))**2*exp(2*(psijn-psijow(iw,1)))*
     &    exp((dfus2o-dfus2n)/(two*tau))

          if(ipr.ge.1) write(6,'(''p'',11f10.6)')
     &    p,(psidn/psidow(iw,1))**2*exp(2*(psijn-psijow(iw,1))),
     &    exp((dfus2o-dfus2n)/(two*tau)),psidn,psidow(iw,1),
     &    psijn,psijow(iw,1),dfus2o,dfus2n

c The following is one reasonable way to cure persistent configurations
c Not needed if itau_eff <=0 and in practice we have never needed it even
c otherwise
          if(iage(iw).gt.50) p=p*1.1d0**(iage(iw)-50)

c         pp=p
          p=dmin1(one,p)
  160     q=one-p

          acc=acc+p
          try_int=try_int+1
          dfus2ac=dfus2ac+p*dfus2o
          dfus2un=dfus2un+dfus2o
          dr2ac=dr2ac+p*dr2
          dr2un=dr2un+dr2

c Calculate density and moments of r for primary walk
          r2o=zero
          r2n=zero
          rmino=zero
          rminn=zero
          do 165 k=1,ndim
            r2o=r2o+xoldw(k,i,iw,1)**2
            r2n=r2n+xnew(k)**2
            rmino=rmino+(xoldw(k,i,iw,1)-cent(k,1))**2
  165       rminn=rminn+(xnew(k)-cent(k,1))**2
          rmino=sqrt(rmino)
          rminn=sqrt(rminn)
          itryo(i)=min(int(delri*rmino)+1,NRAD)
          itryn(i)=min(int(delri*rminn)+1,NRAD)

c for pair-density calculation we will need full old/new positions:
          if(ifixe.lt.0 .or. ifourier.ne.0 .or. izigzag.gt.0) then
            do 167 k=1,ndim
              do 167 j=1,nelec
                xoci(k,j,i)=xoldw(k,j,iw,1)
                if(i.eq.j) then
                  xnci(k,j,i)=xnew(k)
                else
                  xnci(k,j,i)=xoldw(k,j,iw,1)
                endif
  167       continue
          endif

c If we are using weights rather than accept/reject
          if(iacc_rej.le.0) then
            p=one
            q=zero
          endif

          if(rannyu(0).lt.p) then
            iaccept=1
            acc_int=acc_int+1
            if(ipq.le.0) p=one

            iage(iw)=0
            do 170 k=1,ndim
              drifdif=drifdif+(xoldw(k,i,iw,1)-xnew(k))**2
              xoldw(k,i,iw,1)=xnew(k)
              do 170 l=1,nelec
                voldw(k,l,iw,1)=vnew(k,l)
!ACM            write (6,*) 'voldw:'
!ACM            write (6,*) i, voldw(1,i,iw,1)
  170           continue
                
            psidow(iw,1)=psidn
            psijow(iw,1)=psijn
            call jassav(i)
            if(ibasis.eq.3) then                        ! complex calculations
                call cdetsav(i)
             else
                call detsav(i)
            endif

           else
            if(ipq.le.0) p=zero
            call distancese_restore(i)
          endif
          q=one-p

c Calculate moments of r and save rejection probability for primary walk
          r1sume=r1sume+(q*dsqrt(r2o)+p*dsqrt(r2n))
          r2sume=r2sume+(q*r2o+p*r2n)
          risume=risume+(q/dsqrt(r2o)+p/dsqrt(r2n))
          unacp(i)=q

  200   continue ! nelec

c Effective tau for branching
        tauprim=tau*dfus2ac/dfus2un

        do 280 ifr=1,nforce

          if(iaccept.eq.1) then   ! If any of the 1-electron moves were accepted then
            if(ifr.eq.1) then
c Primary configuration
              drifdifr=one
              if(nforce.gt.1)
     &        call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)
              call hpsi(xoldw(1,1,iw,1),psidn,psijn,voldw(1,1,iw,1),div_vow(1,iw),d2n,pen,pein,enew,denergy,1)

              psi_det = psidn                             !JT
              psi_jas = exp(psijn)
              call object_modified_by_index (voldw_index) !JT
              call object_modified_by_index (psi_det_index) !JT
              call object_modified_by_index (psi_jas_index) !JT
              call object_modified_by_index (div_vow_index) !JT

              if(ibasis.eq.3) then                     !complex basis set
                call cwalksav_det(iw)
              else
                call walksav_det(iw)
              endif
              call walksav_jas(iw)
             else
c Secondary configuration
              if(istrech.eq.0) then
                drifdifr=one
c No streched positions for electrons
                do 210 i=1,nelec
                  do 210 k=1,ndim
  210               xoldw(k,i,iw,ifr)=xoldw(k,i,iw,1)
                ajacold(iw,ifr)=one
               else
c Compute streched electronic positions for all nucleus displacement
                call strech(xoldw(1,1,iw,1),xstrech,ajacob,ifr,1)
                drifdifs=zero
                do 220 i=1,nelec
                  do 220 k=1,ndim
                    drifdifs=drifdifs+(xstrech(k,i)-xoldw(k,i,iw,ifr))**2
  220               xoldw(k,i,iw,ifr)=xstrech(k,i)
                ajacold(iw,ifr)=ajacob
                drifdifr=drifdifs/drifdif
              endif
              call hpsi(xoldw(1,1,iw,ifr),psidn,psijn,voldw(1,1,iw,ifr),div_vow(1,iw),d2n,pen,pein,enew,denergy,ifr)
c Temporary test to see how well correlated sampling works for excited states
c Limit outliers such that tau->0 limit is unchanged
c Warning: I should calculate sigma_est
c             sigma_est=0.2d0
c             if(wgcum(1)+wgsum(1).ne.0.d0) then
c               eest_pri=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
c               eest_sec=(egcum(ifr)+egsum(ifr))/(wgcum(ifr)+wgsum(ifr))
c               eest_dif=eest_sec-eest_pri
c               write(6,'(''enew'',9f9.5)')
c    & enew,eest_pri,eest_sec,eest_dif,max(eoldw(iw,1)+eest_dif-3*sigma_est,min(eoldw(iw,1)+eest_dif+3*sigma_est,enew))
ccc             enew=max(eoldw(iw,1)+eest_dif-1/tau,min(eoldw(iw,1)+eest_dif+1/tau,enew))
c               enew=max(eoldw(iw,1)+eest_dif-3*sigma_est,min(eoldw(iw,1)+eest_dif+3*sigma_est,enew))
c              else
c               enew=max(eoldw(iw,1)-sigma_est,min(eoldw(iw,1)+sigma_est,enew))
c             endif
            endif ! ifr.eq.1

            vav2sumn=zero
            v2sumn=zero
            do 260 i=1,nelec

c Use more accurate formula for the drift and tau secondary in drift
              tratio=one
              if(ifr.gt.1.and.itausec.eq.1) tratio=drifdifr

              v2old=0
              do 251 k=1,ndim
  251           v2old=v2old+voldw(k,i,iw,ifr)**2
              vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/
     &            (adrift*v2old)
              vavvn=vavvt/(tau*tratio)

              vav2sumn=vav2sumn+vavvn**2*v2old
  260         v2sumn=v2sumn+v2old
            fration=dsqrt(vav2sumn/v2sumn)
           else
            drifdifr=one
            fration=fratio(iw,ifr)
            enew=eoldw(iw,ifr)
          endif ! accept

          taunow=tauprim*drifdifr

          if(ipr.ge.1)write(6,'(''wt'',9f10.5)') wt(iw),etrial,eest

          ewto=eest-(eest-eoldw(iw,ifr))*fratio(iw,ifr)
          ewtn=eest-(eest-enew)*fration

          do 262 iparm=1,nparm
            dewto(iparm)=denergy_old_dmc(iparm,iw)*fratio(iw,ifr)
  262       dewtn(iparm)=denergy(iparm)*fration

          if(idmc.gt.0) then
            expon=(etrial-half*(ewto+ewtn))*taunow
c Warning we are temporarily ignoring the term that comes from the derivative of (V_av/V) because
c it should be small compared to the term that we keep.
            do 264 iparm=1,nparm
  264         dexponent(iparm)=-half*(dewto(iparm)+dewtn(iparm))*taunow
            if(icut_br.le.0) then
              dwt=dexp(expon)
             elseif(icut_br.eq.1) then
              if(expon.le.0.d0) then
                dwt=dexp(expon)
               else
c Warning: tmp
c               dwt=1+expon+0.5d0*expon**2
                dwt=1+expon/(1+expon)
              endif
             else
              dwt=0.5d0+1/(1+exp(-4*expon))
            endif
          endif

c Warning: These lines were added to reduce the probability of population explosions.
c These occur mostly for nonlocal psps.
c A better solution would be to employ a better way of treating nonlocal psps. in DMC similar to Casula.
c We truncate wts that come from energies that are too low by more than 10*energy_sigma.
c This gives a DMC energy that is too high even in the tau->0 limit, but by a really negligible amount.
c For mpi1 runs a different energy_sigma is calculated on each processor because I did not want to add new MPI calls.
c For mpi2/3 runs a mpi_allreduce is done in the acues1 routine so that it is summed over the processors.
c So, multiply energy_sigma by sqrt(float(nproc)).
c It is more stable to use the energy_sigma with the population control bias than the one with the bias removed.
c         if(iblk.ge.2. or. (iblk.ge.1 .and. nstep.ge.2)) then
          if(ipass-nstep*2*nblkeq .gt. 5) then
            energy_sigma=sigma(ecum1,ecm21,wcum1)
            if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            if(dwt.gt.1+10*energy_sigma*tau) then
              ipr_sav=ipr_sav+1
              if(ipr_sav.le.3) then
                write(6,'(''Warning: dwt>1+10*energy_sigma*tau: nwalk,energy_sigma,dwt,ewto,ewtn,fratio(iw,ifr),fration='',i5,9d12.4
     &          )') nwalk,energy_sigma,dwt,ewto,ewtn,fratio(iw,ifr),fration
                if(ipr_sav.eq.1) write(6,'(''This should add a totally negligible positive bias to the energy'')')
               elseif(ipr_sav.eq.4) then
                write(6,'(''Warning: Additional warning msgs. of dwt>1+10*energy_sigma*tau suppressed'')')
              endif
              dwt=1+10*energy_sigma*tau
            endif
          endif

c ffi has already been raised to wt_lambda.  Do the same for dwt.  We do this even for the current move so that wt_lambda can serve to limit size of move.
          dwt=dwt**wt_lambda

c Exercise population control if dmc or vmc with weights
          if(idmc.gt.0.or.iacc_rej.eq.0) dwt=dwt*ffi

c Set weights and product of weights over last nwprod steps
          if(ifr.eq.1) then

c Raise product of previous generation wts to power wt_lambda_tau to keep product under control if branching is turned off
            wt(iw)=(wt(iw)**wt_lambda_tau)*dwt
c           write(6,'(''wt_lambda_tau='',9d20.12)') wt_lambda_tau, wt(iw)
            do 266 iparm=1,nparm
  266         wi_w(iparm,iw)=wt_lambda_tau*wi_w(iparm,iw)+dexponent(iparm)
            wtnow=wt(iw)
            pwt(iw,ifr)=pwt(iw,ifr)*dwt/wthist(iw,iwmod,ifr)
            wthist(iw,iwmod,ifr)=dwt

           elseif(ifr.gt.1) then

            pwt(iw,ifr)=pwt(iw,ifr)*dwt/wthist(iw,iwmod,ifr)
            wthist(iw,iwmod,ifr)=dwt
            wtnow=wt(iw)*pwt(iw,ifr)/pwt(iw,1)

          endif

          if(ipr.ge.1)write(6,'(''eoldw,enew,wt'',9f10.5)')
     &    eoldw(iw,ifr),enew,wtnow

          wtg=wtnow*fprod
          current_walker_weight = wt(iw) * fprod !JT
          call object_modified_by_index (current_walker_weight_index) !JT
          
          ! used to plot interaction potential (ACM)
          pot_ee_new = pot_ee ! array assignment
          pot_ee_old = pot_ee_oldw(:,iw,ifr) ! array assignment
          
          if(ifr.eq.1) then

            r1sum=r1sum+wtg*r1sume
            r2sum=r2sum+wtg*r2sume
            risum=risum+wtg*risume
            do 270 i=1,nelec
              wtgq=wtg*unacp(i)
              wtgp=wtg-wtgq
              if(i.le.nup) then
                rprobup(itryo(i))=rprobup(itryo(i))+wtgq
                rprobup(itryn(i))=rprobup(itryn(i))+wtgp
               else
                rprobdn(itryo(i))=rprobdn(itryo(i))+wtgq
                rprobdn(itryn(i))=rprobdn(itryn(i))+wtgp
              endif
              rprob(itryo(i))=rprob(itryo(i))+wtgq
              rprob(itryn(i))=rprob(itryn(i))+wtgp
c             r2sum=r2sum+wtg*(unacp(i)*r2o+(one-unacp(i)*r2n)
c 270         risum=risum+wtg*(unacp(i)/dsqrt(r2o)+(one-unacp(i)/dsqrt(r2n))

              if(ifixe.le.-2 .or. ifourier.ne.0 .or. izigzag.gt.0) then
                do j=1,nelec
                  do idim=1,ndim
c note that xoci and xnci represent the old/new positions of all electrons-j when an
c electron-i is being moved
                    xoc(idim,j)=xoci(idim,j,i)
                    xnc(idim,j)=xnci(idim,j,i)
                  enddo
                  if(iperiodic.eq.1) then  ! 1D periodic bc's, so make sure x-posn between -a/2 and a/2
                    call reduce_sim_cell(xoc(:,j))
                    call reduce_sim_cell(xnc(:,j))
                  endif
                enddo
c               write (6,*) 'in dmc_good_ps_mov1:'
c               write (6,*) wtgp, wtgq, i
c               write (6,*) (xoc(1,iii),iii=1,nelec)
c               write (6,*) (xnc(1,iii),iii=1,nelec)
                if(izigzag.gt.0) call zigzag2d(Wtgp,wtgq,xoc,xnc,i)
                if(ifixe.le.-2) call pairden2d(wtgp,wtgq,xoc,xnc)
                if(ifourier.eq.1 .or. ifourier.eq.3) call fourierrk(wtgp,wtgq,xoc,xnc)
                if(ifourier.eq.2 .or. ifourier.eq.3) call fourierkk(wtgp,wtgq,xoc,xnc)
              endif
              
              if(ifixe.eq.-1 .or. ifixe.eq.-3) then
                if(icoosys.eq.1) then 
                  if(iperiodic.eq.1) then  ! 1D periodic bc's, so make sure x-posn between -a/2 and a/2
                    call reduce_sim_cell(xoci(:,i,i))
                    call reduce_sim_cell(xnci(:,i,i))
                  endif
                  do 265 idim=1,ndim
c note that ix can be negative or positive. nint is a better choice.
                    ixo(idim)=nint(delxi(idim)*xoci(idim,i,i))
  265               ixn(idim)=nint(delxi(idim)*xnci(idim,i,i))
                 else
c same trick adapted to circular coordinates
                    ixo(1)=nint(delradi*(dsqrt(xoci(1,i,i)*xoci(1,i,i)+xoci(2,i,i)*xoci(2,i,i))-rmean))
                    ixn(1)=nint(delradi*(dsqrt(xnci(1,i,i)*xnci(1,i,i)+xnci(2,i,i)*xnci(2,i,i))-rmean))
                    ixo(2)=nint(delti*(datan2(xoci(2,i,i),xoci(1,i,i))))
                    ixn(2)=nint(delti*(datan2(xnci(2,i,i),xnci(1,i,i))))
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

  270       continue

          endif
          tausum(ifr)=tausum(ifr)+wtg*taunow

          if(iaccept.eq.1) then   ! If any of the 1-electron moves were accepted then

!JT            if(dabs((enew-etrial)/etrial).gt.1.0d+0) then
!JT             write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass,
!JT            endif

!JT            if(wt(iw).gt.3) write(18,'(i6,i4,3f8.2,30f8.4)') ipass,iw,

            eoldw(iw,ifr)=enew
            do 272 iparm=1,nparm
  272         denergy_old_dmc(iparm,iw)=denergy(iparm)
            peow(iw,ifr)=pen
            peiow(iw,ifr)=pein
            d2ow(iw,ifr)=d2n
            psidow(iw,ifr)=psidn
            psijow(iw,ifr)=psijn
            fratio(iw,ifr)=fration
            pot_ee_oldw(:,iw,ifr) = pot_ee_new ! array assignment
           else
            if(ifr.eq.1) then
              iage(iw)=iage(iw)+1
              ioldest=max(ioldest,iage(iw))
              ioldestmx=max(ioldestmx,iage(iw))
            endif
          endif

          if(ifr.eq.1) then
            psi2savo=2*(dlog(dabs(psidow(iw,1)))+psijow(iw,1))

            wsum1(ifr)=wsum1(ifr)+wtnow
            esum1(ifr)=esum1(ifr)+wtnow*eoldw(iw,ifr)
            pesum(ifr)=pesum(ifr)+wtg*peow(iw,ifr)
            peisum(ifr)=peisum(ifr)+wtg*peiow(iw,ifr)
            tpbsum(ifr)=tpbsum(ifr)+wtg*(eoldw(iw,ifr)-peow(iw,ifr))
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*d2ow(iw,ifr)

!           local energy for current walker
            eloc = eoldw(iw,1)
            call object_modified_by_index (eloc_index)

            call grad_hess_jas_sum(1.d0,0.d0,eoldw(iw,1),eoldw(iw,1))
           else
            ro=ajacold(iw,ifr)*psidow(iw,ifr)**2*
     &         exp(2*psijow(iw,ifr)-psi2savo)
c Warning: Limit weight ratio
c           write(6,'(''eest_pri,eest_sec,eest_dif,ro='',9f9.5)') eest_pri,eest_sec,eest_dif,ro

            wsum1(ifr)=wsum1(ifr)+wtnow*ro
            esum1(ifr)=esum1(ifr)+wtnow*eoldw(iw,ifr)*ro
            pesum(ifr)=pesum(ifr)+wtg*peow(iw,ifr)*ro
            peisum(ifr)=peisum(ifr)+wtg*peiow(iw,ifr)*ro
            tpbsum(ifr)=tpbsum(ifr)+wtg*(eoldw(iw,ifr)-peow(iw,ifr))*ro
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*d2ow(iw,ifr)*ro
          endif

  280   continue
c Call to grad_hess_jas_sum() used to be for optimizing Jastrow for periodic systems.
        call grad_hess_jas_sum(1.d0,0.d0,eoldw(iw,1),eoldw(iw,1),wt(iw)*fprod,wi_w(:,iw))

        call compute_averages_step !JT
  300 continue ! nwalk

!JT      if(wsum1(1).gt.1.1d0*nconf_global) write(18,'(i6,9d12.4)') ipass,ffn,fprod,

      if(idmc.gt.0.or.iacc_rej.eq.0) then
        wfsum1=wsum1(1)*ffn
        efsum1=esum1(1)*ffn
      endif
      do 305 ifr=1,nforce
        if(idmc.gt.0.or.iacc_rej.eq.0) then
          wgsum1(ifr)=wsum1(ifr)*fprod
          egsum1(ifr)=esum1(ifr)*fprod
         else
          wgsum1(ifr)=wsum1(ifr)
          egsum1(ifr)=esum1(ifr)
        endif
  305 continue

      call object_modified_by_index (eoldw_index) !JT
      call object_modified_by_index (wt_index) !JT
      call object_modified_by_index (fprod_index) !JT

!JT      call splitj ! moved outside the routine

      return
      end
