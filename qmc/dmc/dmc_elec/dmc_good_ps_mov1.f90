      subroutine dmc_good_ps_mov1
! Written by Cyrus Umrigar and Claudia Filippi
! Uses the diffusion Monte Carlo algorithm described in:
! 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
!    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993)
! modified to do accept/reject after single-electron moves and to
! remove portions related to nuclear cusps.
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
!              = 1 *   use first 2 terms in Taylor expansion of Exp if exponent>0
!             >= 2 *   use smooth formulae to limit branching to (1/2,2)
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
      use distance_mod, only: pot_ee, rshift, rvec_en, r_en
      use config_mod, only: pot_ee_new, pot_ee_old
      use zigzag_mod, only: izigzag
      use qua_mod, only: l_do_tmoves, iaccept_tmove
      use pseudo_mod, only: nloc
      use dete_mod, only: eval_grad, dete_save                 !TA
      use jaso_mod, only: fjo                                  !TA
      use deriv_fast_mod, only: aiup, aidn, tup, tdn, yup, ydn !TA
      use orb_mod, only: dorb                                  !TA
      use pseudo_mod, only: nloc

      implicit real*8(a-h,o-z)

      parameter (eps=1.d-10)

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

!     gauss()=dcos(two*pi*rannyu(0))*sqrt(-two*dlog(rannyu(0)))
      e_sigma(x,x2,w)=sqrt(max((x2/w-(x/w)**2)*nconf,0.d0))

!     rn_eff(w,w2)=w**2/w2
!     error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
!     errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

      if(idiv_v.ge.1) stop 'div_v not implemented yet in 1-electron move algorithm'

!     term=(sqrt(two*pi*tau))**3/pi

! Temporarily create wt_lambda_tau here, but it should be done at start of program
      wt_lambda_tau=wt_lambda**tau

! Undo products
      ipmod=mod(ipass,nfprod)
      ipmod2=mod(ipass+1,nfprod)
      if(idmc.gt.0) then
!     if(idmc.gt.0 .and. .not. l_opt_ovlp_fn) then
        ginv=min(1.d0,tau)
!       ginv=min(1.d0,10*tau)
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
        dwt=1
      endif

! Undo weights
      iwmod=mod(ipass,nwprod)

! Store fratio(iw,ifr) (well behaved velocity/velocity)
      if(ncall.eq.0.and.irstar.eq.0) then
        do 7 iw=1,nwalk
          do 7 ifr=1,nforce

            vav2sumo=zero
            v2sumo=zero
            do 6 i=1,nelec

! Tau secondary in drift is one (first time around)
              tratio=one

              v2old=0
              do 5 k=1,ndim
    5           v2old=v2old+voldw(k,i,iw,ifr)**2

              if(drift_type=='unr93' .or. adrift <= 1.d0) then
                vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/ (adrift*v2old)
              elseif(drift_type=='quadratic' .and. adrift > 1.d0) then
                vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/ (adrift*v2old +(adrift-1)*sqrt(2*adrift*v2old*tau*tratio))
              endif
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
! Loop over primary walker
        current_walker = iw !JT
        call object_modified_by_index (current_walker_index) !JT

! Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
        if(nforce.gt.1) call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)

        if(ibasis.eq.3) then          ! complex basis set
          call cwalkstrdet(iw)
        else
          call walkstrdet(iw)
        endif
        call walkstrjas(iw)

! This was put in when working on tmoves.  It is commented out because it appears to be unnecessary if distances is called from hpsie.f
!       call distances(xoldw(1,1,iw,1),pe,pei)

! Sample Green function for forward move
        r1sume=zero
        r2sume=zero
        r3sume=zero
        r4sume=zero
        risume=zero
        dfus2ac=zero
        dfus2un=zero
!       dr2ac=zero
!       dr2un=zero
        drifdif=zero
        iaccept=0
        l_do_tmoves=.false.
        psum=0

        ! calculate v2sumo for ifr=1 TA
        v2sumo=0d0
        do i=1,nelec
          do k=1,ndim
            v2sumo=v2sumo+voldw(k,i,iw,1)**2
          enddo
        enddo

        call distances(xoldw(1,1,iw,1),pe,pei)

        do i=1,nelec
          if(nloc.eq.1) then
            call getvps_fahy(r_en,i)
           elseif(nloc.ge.2 .and. nloc.le.5) then
            call getvps_champ(r_en,i)
           elseif(nloc.eq.6) then
            call getvps_gauss(r_en,i)
           else
            stop 'nloc < 1 or 6 < nloc'
          endif
        enddo

        if(tmoves) then
          call tmove
          if(iaccept_tmove.eq.1) iaccept=1
        endif

        do 200 i=1,nelec

! Save distances (not really needed for locality approx. but is needed for tmoves)
          call distancese(i,xoldw(1,1,iw,1))

          if (i.le.nup) call eval_grad(i, dorb(:,i,:), aiup, tup, yup, voldw(:,i,iw,1)) !TA
          if (i.gt.nup) call eval_grad(i, dorb(:,i,:), aidn, tdn, ydn, voldw(:,i,iw,1))
          voldw(1,i,iw,1) = voldw(1,i,iw,1) + fjo(1,i)
          voldw(2,i,iw,1) = voldw(2,i,iw,1) + fjo(2,i)
          voldw(3,i,iw,1) = voldw(3,i,iw,1) + fjo(3,i)

! Use more accurate formula for the drift
          v2old=0
          do 75 k=1,ndim
   75       v2old=v2old+voldw(k,i,iw,1)**2
! Tau primary -> tratio=one

          if(drift_type=='unr93' .or. adrift <= 1.d0) then
            vavvt=(dsqrt(one+two*adrift*v2old*tau)-one)/ (adrift*v2old)
          elseif(drift_type=='quadratic' .and. adrift > 1.d0) then
            vavvt=(dsqrt(one+two*adrift*v2old*tau)-one)/ (adrift*v2old +(adrift-1)*sqrt(2*adrift*v2old*tau))
          endif

          dr2=zero
          dfus2o=zero
          do 80 k=1,ndim
            drift=vavvt*voldw(k,i,iw,1)
            dfus=gauss()*rttau
            dx=drift+dfus
            dr2=dr2+dx**2
            dfus2o=dfus2o+dfus**2
   80       xnew(k)=xoldw(k,i,iw,1)+dx

!      if(v2old.lt.-1.d99 .or. v2old.gt.1.d99) write(6,'(''v2old'',d12.4)') v2old
!      if(abs(v2old).gt.1.d99) write(6,'(''v2old'',d12.4)') v2old
!      write(6,'(''i, iw, adrift, v2old, (voldw(:,:,iw,1)'',2i4,99d12.4)')
!    &  i, iw, adrift, v2old, ((voldw(k,ii,iw,1),k=1,3),ii=1,nelec)
!      write(6,'(''(xoldw(:,:,iw,1)'',99d22.14)')
!    &  ((xoldw(k,ii,iw,1),k=1,3),ii=1,nelec)
!      call systemflush(6)

          if(ipr.ge.1) write(6,'(''xnewdr'',2i4,9f8.5)') iw,i,(xnew(k),k=1,ndim)
! calculate psi and velocity at new configuration
          call hpsiedmc(i,iw,xnew,psidn,psijn,vnew)

! Check for node crossings
          if(psidn*psidow(iw,1).le.zero) then
            nodecr=nodecr+1
            if(icross.le.0) then
              p=zero
              goto 160
            endif
          endif

! Calculate Green function for the reverse move

          v2new=0
          do 149 k=1,ndim
  149       v2new=v2new+vnew(k,i)**2

          if(drift_type=='unr93' .or. adrift <= 1.d0) then
            vavvt=(dsqrt(one+two*adrift*v2new*tau)-one)/ (adrift*v2new)
          elseif(drift_type=='quadratic' .and. adrift > 1.d0) then
            vavvt=(dsqrt(one+two*adrift*v2new*tau)-one)/ (adrift*v2new +(adrift-1)*sqrt(2*adrift*v2new*tau))
          endif

          dfus2n=zero
          do 150 k=1,ndim
            drift=vavvt*vnew(k,i)
            xbac(k)=xnew(k)+drift
            dfus=xbac(k)-xoldw(k,i,iw,1)
  150       dfus2n=dfus2n+dfus**2

          if(ipr.ge.1) then
            write(6,'(''xoldw'',9f10.6)')(xoldw(k,i,iw,1),k=1,ndim), &
     &      (xnew(k),k=1,ndim), (xbac(k),k=1,ndim)
            write(6,'(''dfus2o'',9f10.6)')dfus2o,dfus2n, &
     &      psidow(iw,1),psidn,psijow(iw,1),psijn
          endif

          p=(psidn/psidow(iw,1))**2*exp(2*(psijn-psijow(iw,1)))* &
     &    exp((dfus2o-dfus2n)/(two*tau))

          if(ipr.ge.1) write(6,'(''p'',11f10.6)') &
     &    p,(psidn/psidow(iw,1))**2*exp(2*(psijn-psijow(iw,1))), &
     &    exp((dfus2o-dfus2n)/(two*tau)),psidn,psidow(iw,1), &
     &    psijn,psijow(iw,1),dfus2o,dfus2n

! The following is one reasonable way to cure persistent configurations
! Not needed if itau_eff <=0 and in practice we have never needed it even
! otherwise
          if(idmc > 0 .and. iage(iw).gt.50) p=p*1.1d0**(iage(iw)-50)

!         pp=p
          p=dmin1(one,p)
  160     q=one-p
          psum=psum+p

          acc=acc+p
          try_int=try_int+1
          dfus2ac=dfus2ac+p*dfus2o
          dfus2un=dfus2un+dfus2o
          dr2ac=dr2ac+p*dr2
          dr2un=dr2un+dr2

! Calculate density and moments of r for primary walk
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
!         itryo(i)=min(int(delri*rmino)+1,NRAD)
!         itryn(i)=min(int(delri*rminn)+1,NRAD)
          itryo(i)=int(min(delri*rmino+1,dfloat(NRAD))+eps)
          itryn(i)=int(min(delri*rminn+1,dfloat(NRAD))+eps)

! for pair-density calculation we will need full old/new positions:
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

! If we are using weights rather than accept/reject
          if(iacc_rej.le.0) then
            p=one
            q=zero
          endif

          if(rannyu(0).lt.p) then
            iaccept=1
            acc_int=acc_int+1
            if(ipq.le.0) p=one

            iage(iw)=0

!           Since move is accepted, copy from new to old
            do 170 k=1,ndim
              drifdif=drifdif+(xoldw(k,i,iw,1)-xnew(k))**2
              xoldw(k,i,iw,1)=xnew(k)
              do 170 l=1,nelec
                voldw(k,l,iw,1)=vnew(k,l)
  170           continue
            psidow(iw,1)=psidn
            psijow(iw,1)=psijn
            call jassav(i)
            if(ibasis.eq.3) then                        ! complex calculations
                call cdetsav(i)
            else
                call detsav(i)
                call dete_save(i) !TA
            endif

          else
            if(ipq.le.0) p=zero
            call distancese_restore(i)
          endif
          q=one-p

! Calculate moments of r and save rejection probability for primary walk
          r1sume=r1sume+(q*dsqrt(r2o)+p*dsqrt(r2n))
          r2sume=r2sume+(q*r2o+p*r2n)
          r3sume=r3sume+(q*r2o*dsqrt(r2o)+p*r2n*dsqrt(r2n))
          r4sume=r4sume+(q*r2o*r2o+p*r2n*r2n)
          risume=risume+(q/dsqrt(r2o)+p/dsqrt(r2n))
          unacp(i)=q

  200   continue ! nelec

! Effective tau for branching
        tauprim=tau*dfus2ac/dfus2un  ! taunow is used in nonloc, called from hpsi.  taunow=tauprim for ifr=1, tauprim*drifdifr for ifr>1

        do 280 ifr=1,nforce

          if(ifr.gt.1) then !Calculate v2sumo for ifr != 1 - TA
            v2sumo=0d0
            do i=1,nelec
              do k=1,ndim
                v2sumo=v2sumo+voldw(k,i,iw,ifr)**2
              enddo
            enddo
          endif

          if(iaccept.eq.1) then   ! If any of the 1-electron moves were accepted then
            if(ifr.eq.1) then
! Primary configuration
              drifdifr=one
              if(nforce.gt.1) call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)
              taunow=tauprim
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

            else                                       ! ifr>1, i.e. secondary walk
! Secondary configuration

              if(istrech.eq.0) then
                drifdifr=one
! No streched positions for electrons
                do 210 i=1,nelec
                  do 210 k=1,ndim
  210               xoldw(k,i,iw,ifr)=xoldw(k,i,iw,1)
                ajacold(iw,ifr)=one
              else
! Compute streched electronic positions for all nucleus displacement
                call strech(xoldw(1,1,iw,1),xstrech,ajacob,ifr,1)
                drifdifs=zero
                do 220 i=1,nelec
                  do 220 k=1,ndim
                    drifdifs=drifdifs+(xstrech(k,i)-xoldw(k,i,iw,ifr))**2
  220               xoldw(k,i,iw,ifr)=xstrech(k,i)
                ajacold(iw,ifr)=ajacob
                drifdifr=drifdifs/drifdif
              endif
              taunow=tauprim*drifdifr
              call hpsi(xoldw(1,1,iw,ifr),psidn,psijn,voldw(1,1,iw,ifr),div_vow(1,iw),d2n,pen,pein,enew,denergy,ifr)
! Temporary test to see how well correlated sampling works for excited states
! Limit outliers such that tau->0 limit is unchanged
! Warning: I should calculate sigma_est
!             sigma_est=0.2d0
!             if(wgcum(1)+wgsum(1).ne.0.d0) then
!               eest_pri=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
!               eest_sec=(egcum(ifr)+egsum(ifr))/(wgcum(ifr)+wgsum(ifr))
!               eest_dif=eest_sec-eest_pri
!               write(6,'(''enew'',9f9.5)')
!    & enew,eest_pri,eest_sec,eest_dif,max(eoldw(iw,1)+eest_dif-3*sigma_est,min(eoldw(iw,1)+eest_dif+3*sigma_est,enew))
!cc             enew=max(eoldw(iw,1)+eest_dif-1/tau,min(eoldw(iw,1)+eest_dif+1/tau,enew))
!               enew=max(eoldw(iw,1)+eest_dif-3*sigma_est,min(eoldw(iw,1)+eest_dif+3*sigma_est,enew))
!             else
!               enew=max(eoldw(iw,1)-sigma_est,min(eoldw(iw,1)+sigma_est,enew))
!             endif
            endif ! ifr.eq.1

            vav2sumn=zero
            v2sumn=zero
            do 260 i=1,nelec

! Use more accurate formula for the drift and tau secondary in drift
              tratio=one
              if(ifr.gt.1.and.itausec.eq.1) tratio=drifdifr

              v2old=0
              do 251 k=1,ndim
  251           v2old=v2old+voldw(k,i,iw,ifr)**2

              if(drift_type=='unr93' .or. adrift <= 1.d0) then
                vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/ (adrift*v2old)
              elseif(drift_type=='quadratic' .and. adrift > 1.d0) then
                vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/ (adrift*v2old +(adrift-1)*sqrt(2*adrift*v2old*tau*tratio))
              endif
              vavvn=vavvt/(tau*tratio)

              vav2sumn=vav2sumn+vavvn**2*v2old
  260         v2sumn=v2sumn+v2old
            fration=dsqrt(vav2sumn/v2sumn)
          else  ! none of the 1-electron drift-diffusion moves were accepted
            drifdifr=one
            fration=fratio(iw,ifr)
            enew=eoldw(iw,ifr)
          endif ! accept

!         taunow=tauprim*drifdifr

          if(ipr.ge.1)write(6,'(''wt'',9f10.5)') wt(iw),etrial,eest

!There are 2 kinds of limits we put on the reweight factor.
!For all-electron calculations the purpose is to reduce the time-step error.
!For locality-approximation pseudpotential calculations one of both of these are necessary to avoid large negative spikes in the energy.
!First we use e=eest-(eest-e)*f(fratio) where f is a function of fratio=vel_av/vel that is 0 at 0, 1 at 1 and deviates from 1 at 1 as some power, e.g.
!UNR93:        f(x)=x                           deviates from 1 and 0 linearly.
!new_ene_int:  f(x)=x**2*(3-2*x)                deviates from 1 and 0 quadratically.
!new_ene_int2: f(x)=x**3*(10-15*x+6*x**2)       deviates from 1 and 0 cubically.
!new_ene_int3: f(x)=1-(1-x)**3 = x*(3-3*x+x**2) deviates from 1 cubically and from 0 linearly.
!no_ene_int:   f(x)=1                           does not deviate from 1     
!All these will give the same energy at tau=0, but progressively lower energies, going from UNR93 to no_ene_int, at finite values of tau
!The second limit is described and imposed about 50 lines down.

! Warning: Change UNR93 reweighting factor because it gives large time-step error at small tau for pseudo systems as pointed out by Alfe
          if(ene_int=='unr93') then
            ewto=eest-(eest-eoldw(iw,ifr))*fratio(iw,ifr)                                               ! UNR93
            ewtn=eest-(eest-enew)*fration                                                               ! UNR93
          elseif(ene_int=='new_ene_int') then
            ewto=eest-(eest-eoldw(iw,ifr))*fratio(iw,ifr)**2*(3-2*fratio(iw,ifr))                       ! new_ene_int
            ewtn=eest-(eest-enew)*fration**2*(3-2*fration)                                              ! new_ene_int
          elseif(ene_int=='new_ene_int2') then
            ewto=eest-(eest-eoldw(iw,ifr))*fratio(iw,ifr)**3*(10-15*fratio(iw,ifr)+6*fratio(iw,ifr)**2) ! new_ene_int2
            ewtn=eest-(eest-enew)*fration**3*(10-15*fration+6*fration**2)                               ! new_ene_int2
          elseif(ene_int=='new_ene_int3') then
            ewto=eest-(eest-eoldw(iw,ifr))*(1-(1-fratio(iw,ifr))**3)                                    ! new_ene_int3
            ewtn=eest-(eest-enew)*(1-(1-fration)**3)                                                    ! new_ene_int3
          elseif(ene_int=='new_ene_int4') then
            if(ipass .gt. nstep*2*nblkeq + max(2,nint(1.d0/tau))) then
              energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
              if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            else
              energy_sigma=0.2d0*sqrt(real(nelec))
            endif
            ecut=10*energy_sigma/(1+sqrt(tau))                                                          ! new_ene_int4
            ewto=eest+min(ecut,max(-ecut,eoldw(iw,ifr)-eest))                                           ! new_ene_int4
            ewtn=eest+min(ecut,max(-ecut,enew-eest))                                                    ! new_ene_int4
          elseif(ene_int=='new_ene_int5') then
            if(ipass .gt. nstep*2*nblkeq + max(2,nint(1.d0/tau))) then
              energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
              if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            else
              energy_sigma=0.2d0*sqrt(real(nelec))
            endif
            ecut=10*energy_sigma/(1+tau)                                                                ! new_ene_int5
            ewto=eest+min(ecut,max(-ecut,eoldw(iw,ifr)-eest))                                           ! new_ene_int5
            ewtn=eest+min(ecut,max(-ecut,enew-eest))                                                    ! new_ene_int5
          elseif(ene_int=='new_ene_int7') then
            if(ipass .gt. nstep*2*nblkeq + max(2,nint(1.d0/tau))) then
              energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
              if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            else
              energy_sigma=0.2d0*sqrt(real(nelec))
            endif
            ecut=10*energy_sigma/(1+v2sumo*tau/nelec)                                                   ! new_ene_int7
            ewto=eest+min(ecut,max(-ecut,eoldw(iw,ifr)-eest))                                           ! new_ene_int7
            ecut=10*energy_sigma/(1+v2sumn*tau/nelec)                                                   ! new_ene_int7
            ewtn=eest+min(ecut,max(-ecut,enew-eest))                                                    ! new_ene_int7
          elseif(ene_int=='new_ene_int8') then
            if(ipass .gt. nstep*2*nblkeq + max(2,nint(1.d0/tau))) then
              energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
              if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            else
              energy_sigma=0.2d0*sqrt(real(nelec))
            endif
!           if(ipass.le.10000) write(6,'(''ipass,etrial,eest,energy_sigma'',i7,9f10.4)') ipass,etrial,eest,energy_sigma
            ecut=min(abs(eest-eoldw(iw,ifr)),10*energy_sigma)                                           ! new_ene_int8
            ewto=eest-sign(ecut,eest-eoldw(iw,ifr))/(1+(v2sumo*tau/nelec)**2)                           ! new_ene_int8
            ecut=min(abs(eest-enew),10*energy_sigma)                                                    ! new_ene_int8
            ewtn=eest-sign(ecut,eest-enew)/(1+(v2sumn*tau/nelec)**2)                                    ! new_ene_int8
          elseif(ene_int=='new_ene_int9') then
            if(ipass .gt. nstep*2*nblkeq + max(2,nint(1.d0/tau))) then
              energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
              if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            else
              energy_sigma=0.2d0*sqrt(real(nelec))
            endif
            ecut=min(abs(eest-eoldw(iw,ifr)),10*energy_sigma)                                           ! new_ene_int9
            ewto=eest-sign(ecut,eest-eoldw(iw,ifr))/(1+(v2sumo*tau/nelec))                              ! new_ene_int9
            ecut=min(abs(eest-enew),10*energy_sigma)                                                    ! new_ene_int9
            ewtn=eest-sign(ecut,eest-enew)/(1+(v2sumn*tau/nelec))                                       ! new_ene_int9
!           write(6,'(''v2sumo,v2sumn,ewto,ewtn,ecut'',9es12.4)') v2sumo,v2sumn,ewto,ewtn,ecut
          elseif(ene_int=='new_ene_int10') then
            if(ipass .gt. nstep*2*nblkeq + max(2,nint(1.d0/tau))) then
              energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
              if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            else
              energy_sigma=0.2d0*sqrt(real(nelec))
            endif
            ecut=min(abs(eest-eoldw(iw,ifr)),10*energy_sigma)                                           ! new_ene_int10
           !ewto=eest-sign(ecut,eest-eoldw(iw,ifr))/(1+(v2sumo*tau/(4*nelec)))                          ! new_ene_int10
            ewto=eest-sign(ecut,eest-eoldw(iw,ifr))/(1+(sqrt(v2sumo)*tau/nelec))                        ! new_ene_int10
            ecut=min(abs(eest-enew),10*energy_sigma)                                                    ! new_ene_int10
           !ewtn=eest-sign(ecut,eest-enew)/(1+(v2sumn*tau/(4*nelec)))                                   ! new_ene_int10
            ewtn=eest-sign(ecut,eest-enew)/(1+(sqrt(v2sumn)*tau/nelec))                                 ! new_ene_int10
          elseif(ene_int=='new_ene_int11') then
            if(ipass .gt. nstep*2*nblkeq + max(2,nint(1.d0/tau))) then
              energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
              if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            else
              energy_sigma=0.2d0*sqrt(real(nelec))
            endif
            ecut=min(abs(eest-eoldw(iw,ifr)),10*energy_sigma)                                           ! new_ene_int11
           !ewto=eest-sign(ecut,eest-eoldw(iw,ifr))/(1+(v2sumo*tau/(4*nelec)))                          ! new_ene_int11
            ewto=eest-sign(ecut,eest-eoldw(iw,ifr))/(1+2*(sqrt(v2sumo/nelec)*tau**2))                   ! new_ene_int11
            ecut=min(abs(eest-enew),10*energy_sigma)                                                    ! new_ene_int11
           !ewtn=eest-sign(ecut,eest-enew)/(1+(v2sumn*tau/(4*nelec)))                                   ! new_ene_int11
            ewtn=eest-sign(ecut,eest-enew)/(1+2*(sqrt(v2sumn/nelec)*tau**2))                            ! new_ene_int11
          elseif(ene_int=='new_ene_int12') then
            if(ipass .gt. nstep*2*nblkeq + max(2,nint(1.d0/tau))) then
              energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
              if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            else
              energy_sigma=0.2d0*sqrt(real(nelec))
            endif
            ecut=min(abs(eest-eoldw(iw,ifr)),10*energy_sigma)                                           ! new_ene_int12
           !ewto=eest-sign(ecut,eest-eoldw(iw,ifr))/(1+abs(eest-eoldw(iw,ifr))*tau)                     ! new_ene_int12
            ewto=eest-sign(ecut,eest-eoldw(iw,ifr))/(1+ecut*tau)                                        ! new_ene_int12
            ecut=min(abs(eest-enew),10*energy_sigma)                                                    ! new_ene_int12
           !ewtn=eest-sign(ecut,eest-enew)/(1+abs(eest-enew)*tau)                                       ! new_ene_int12
            ewtn=eest-sign(ecut,eest-enew)/(1+ecut*tau)                                                 ! new_ene_int12
          elseif(ene_int=='new_ene_int13') then                                                         ! cut only the low energy side
            if(ipass .gt. nstep*2*nblkeq + max(2,nint(1.d0/tau))) then
              energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
              if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            else
              energy_sigma=0.2d0*sqrt(real(nelec))
            endif
            ecut=min(eest-eoldw(iw,ifr),10*energy_sigma)                                                ! new_ene_int13
            ewto=eest-sign(ecut,eest-eoldw(iw,ifr))/(1+abs(ecut)*tau)                                   ! new_ene_int13
            ecut=min(eest-enew,10*energy_sigma)                                                         ! new_ene_int13
            ewtn=eest-sign(ecut,eest-enew)/(1+abs(ecut)*tau)                                            ! new_ene_int13
          elseif(ene_int=='no_ene_int') then
            ewto=eoldw(iw,ifr)                                                                          ! no_ene_int
            ewtn=enew                                                                                   ! no_ene_int
          elseif(ene_int=='alfe') then
            ecut=0.2*sqrt(nelec/tau)                                                                    ! Alfe
            ewto=eest+min(ecut,max(-ecut,eoldw(iw,ifr)-eest))                                           ! Alfe
            ewtn=eest+min(ecut,max(-ecut,enew-eest))                                                    ! Alfe
          endif

          do 262 iparm=1,nparm
            dewto(iparm)=denergy_old_dmc(iparm,iw)*fratio(iw,ifr)
  262       dewtn(iparm)=denergy(iparm)*fration

          psum=psum/nelec
          qsum=1-psum

          if(idmc.gt.0) then
! rewt_type=='pq' gives a higher energy than rewt_type=='sym' at all tau.
! So, it gives smaller negative time-step error than the usual rewt_type=='sym' at large tau, but does it give an improvement at small tau?
            if(rewt_type=='sym') then
              expon=(etrial-half*(ewto+ewtn))*taunow
            elseif(rewt_type=='pq') then
              expon=(etrial-((1-half*psum)*ewto+half*psum*ewtn))*taunow
            endif
! Warning we are temporarily ignoring the term that comes from the derivative of (V_av/V) or other modification of reweighting because
! it should be small compared to the term that we keep.
            do 264 iparm=1,nparm
  264         dexponent(iparm)=-half*(dewto(iparm)+dewtn(iparm))*taunow
            if(icut_br.le.0) then
              dwt=dexp(expon)
            elseif(icut_br.eq.1) then
              if(expon.le.0.d0) then
                dwt=dexp(expon)
              else
                dwt=1+expon/(1+expon)
              endif
            else
              dwt=0.5d0+1/(1+exp(-4*expon))
            endif
          endif

! Warning: These lines were added to reduce the probability of population explosions.
! These occur for nonlocal psps. with the locality approx, and are cured by our slightly modified version of Casula et al'
! size-consistent tmoves version 1 and by our version of tmoves that has an additional accept-reject step.  For tmoves, these lines are no longer needed.
! If we truncate wts that are larger than 1+limit_rewt_dmc*energy_sigma then there would be no effect as tau->0.
! If we truncate wts that are larger than 1+limit_rewt_dmc*energy_sigma*sqrt(tau) then also there would be no effect as tau->0.
! If we truncate wts that are larger than 1+limit_rewt_dmc*energy_sigma*tau then there is a bias, but it is tiny if energy_sigma > 10.
! However, for reasons I need to understand, when using the locality approx. we get a too low energy as tau->0, if limit_rewt_dmc > 10.
! (see e.g. /data/cyrus/qmc_runs/atoms_ps/c_dolg_2zeta/csf01/test_timestep_new_ene_int2_20sigma.)
! For mpi1 runs a different energy_sigma is calculated on each processor because I did not want to add new MPI calls.
! For mpi2/3 runs an mpi_allreduce is done in the acues1 routine so that it is summed over the processors.
! So, multiply energy_sigma by sqrt(float(nproc)).
! It is more stable to use the energy_sigma with the population control bias than the one with the bias removed.
!         if(ipass-nstep*2*nblkeq .gt. 5) then
          if(ipass .gt. nstep*2*nblkeq + max(10,nint(10.d0/tau))) then
!           energy_sigma=e_sigma(ecum1,ecm21,wcum1)
            energy_sigma=e_sigma(egcum1(1),egcm21(1),wgcum1(1))
            if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
          else
            energy_sigma=0.2d0*sqrt(real(nelec))
          endif
!         if(iw==1) write(6,'(''ipass,e_sigma(ecum1,ecm21,wcum1),e_sigma(egcum(1),egcm2(1),wgcum(1)),energy_sigma,dwt,1+limit_rewt_dmc*energy_sigma*tau'',i6,9f10.6)') ipass,e_sigma(ecum1,ecm21,wcum1),e_sigma(egcum1(1),egcm21(1),wgcum1(1)),e_sigma(egcum(1),egcm2(1),wgcum(1)),energy_sigma,dwt,1+limit_rewt_dmc*energy_sigma*tau
!         if(dwt.gt.1+limit_rewt_dmc*energy_sigma*tau) then
          if(dwt.gt.exp((etrial-eest+limit_rewt_dmc*energy_sigma)*tau)) then
            ipr_sav=ipr_sav+1
            if(ipr_sav.le.3) then
              write(6,'(''Warning: dwt>exp((etrial-eest+limit_rewt_dmc*energy_sigma)*tau): nwalk,energy_sigma,dwt,ewto,ewtn,fratio(iw,ifr),fration='',i5,9d12.4)') &
     &        nwalk,energy_sigma,dwt,ewto,ewtn,fratio(iw,ifr),fration
              if(ipr_sav.eq.1) write(6,'(''This should add a totally negligible positive bias to the energy'')')
            elseif(ipr_sav.eq.4) then
              write(6,'(''Warning: Additional warning msgs. of dwt>1+limit_rewt_dmc*energy_sigma*tau suppressed'')')
            endif
!           dwt=1+limit_rewt_dmc*energy_sigma*tau
            dwt=exp((etrial-eest+limit_rewt_dmc*energy_sigma)*tau)
          endif

! ffi has already been raised to wt_lambda.  Do the same for dwt.  We do this even for the current move so that wt_lambda can serve to limit size of move.
          dwt=dwt**wt_lambda

! Exercise population control if dmc or vmc with weights
          if(idmc.gt.0.or.iacc_rej.eq.0) dwt=dwt*ffi

! Set weights and product of weights over last nwprod steps
          if(ifr.eq.1) then

! Raise product of previous generation wts to power wt_lambda_tau to keep product under control if branching is turned off
            wt(iw)=(wt(iw)**wt_lambda_tau)*dwt
!           write(6,'(''wt_lambda_tau='',9d20.12)') wt_lambda_tau, wt(iw)
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

          if(ipr.ge.1)write(6,'(''eoldw,enew,wt'',9f10.5)') eoldw(iw,ifr),enew,wtnow

          wtg=wtnow*fprod
          current_walker_weight = wt(iw) * fprod !JT
          call object_modified_by_index (current_walker_weight_index) !JT

          ! used to plot interaction potential (ACM)
          pot_ee_new = pot_ee ! array assignment
          pot_ee_old = pot_ee_oldw(:,iw,ifr) ! array assignment

          if(ifr.eq.1) then

            r1sum=r1sum+wtg*r1sume
            r2sum=r2sum+wtg*r2sume
            r3sum=r3sum+wtg*r3sume
            r4sum=r4sum+wtg*r4sume
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
!             r2sum=r2sum+wtg*(unacp(i)*r2o+(one-unacp(i)*r2n)
! 270         risum=risum+wtg*(unacp(i)/dsqrt(r2o)+(one-unacp(i)/dsqrt(r2n))

              if(ifixe.le.-2 .or. ifourier.ne.0 .or. izigzag.gt.0) then
                do j=1,nelec
                  do idim=1,ndim
! note that xoci and xnci represent the old/new positions of all electrons-j when an
! electron-i is being moved
                    xoc(idim,j)=xoci(idim,j,i)
                    xnc(idim,j)=xnci(idim,j,i)
                  enddo
                  if(iperiodic.eq.1) then  ! 1D periodic bc's, so make sure x-posn between -a/2 and a/2
                    call reduce_sim_cell(xoc(:,j))
                    call reduce_sim_cell(xnc(:,j))
                  endif
                enddo
!               write (6,*) 'in dmc_good_ps_mov1:'
!               write (6,*) wtgp, wtgq, i
!               write (6,*) (xoc(1,iii),iii=1,nelec)
!               write (6,*) (xnc(1,iii),iii=1,nelec)
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
! note that ix can be negative or positive. nint is a better choice.
                    ixo(idim)=nint(delxi(idim)*xoci(idim,i,i))
  265               ixn(idim)=nint(delxi(idim)*xnci(idim,i,i))
                else
! same trick adapted to circular coordinates
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
            eloc_dmc = eoldw(iw,1)
            call object_modified_by_index (eloc_dmc_index)

            call grad_hess_jas_sum(1.d0,0.d0,eoldw(iw,1),eoldw(iw,1))
          else
            ro=ajacold(iw,ifr)*psidow(iw,ifr)**2*exp(2*psijow(iw,ifr)-psi2savo)
! Warning: Limit weight ratio
!           write(6,'(''eest_pri,eest_sec,eest_dif,ro='',9f9.5)') eest_pri,eest_sec,eest_dif,ro

            wsum1(ifr)=wsum1(ifr)+wtnow*ro
            esum1(ifr)=esum1(ifr)+wtnow*eoldw(iw,ifr)*ro
            pesum(ifr)=pesum(ifr)+wtg*peow(iw,ifr)*ro
            peisum(ifr)=peisum(ifr)+wtg*peiow(iw,ifr)*ro
            tpbsum(ifr)=tpbsum(ifr)+wtg*(eoldw(iw,ifr)-peow(iw,ifr))*ro
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*d2ow(iw,ifr)*ro
          endif

  280   continue ! ifr=1,nforce
! Call to grad_hess_jas_sum() used to be for optimizing Jastrow for periodic systems.
        call grad_hess_jas_sum(1.d0,0.d0,eoldw(iw,1),eoldw(iw,1),wt(iw)*fprod,wi_w(:,iw))
        call compute_averages_step !JT

        call systemflush(6)

!       write(6,'(''iw,xoldw(k,i,iw,1)'',i3,99d12.4)') iw,((xoldw(k,i,iw,1),k=1,3),i=1,nelec)
!       write(6,'(''iw,voldw(k,i,iw,1)'',i3,99d12.4)') iw,((voldw(k,i,iw,1),k=1,3),i=1,nelec)

  300 continue ! iw=1,nwalk

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
