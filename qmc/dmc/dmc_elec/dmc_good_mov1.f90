      subroutine dmc_good_mov1
! Written by Cyrus Umrigar and Claudia Filippi
! Uses the diffusion Monte Carlo algorithm described in:
! 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
!    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993)
! modified to do accept/reject after single-electron moves.
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
      use velratio_mod
      use pop_control_mod, only : ffn
      use determinants_mod
      use eloc_mod
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-10,huge=1.d+100,adrift0=0.1d0)

      dimension rvmino(3),rvminn(3),xstrech(3,nelec)
      dimension xnew(3),vnew(3,nelec)
      dimension xaxis(3),zaxis(3),xbac(3)
      dimension itryo(nelec),itryn(nelec),unacp(nelec)
      dimension dewto(nparm),dewtn(nparm),dexponent(nparm)

      data ncall,ipr_sav /0,0/
      save ipr_sav

!     gauss()=dcos(two*pi*rannyu(0))*sqrt(-two*dlog(rannyu(0)))
      e_sigma(x,x2,w)=sqrt(max((x2/w-(x/w)**2)*nconf,0.d0))

      if(icuspg.ge.1) stop 'exact cusp in G not implemented yet in 1-electron move algorithm'
      if(idiv_v.ge.1) stop 'div_v not implemented yet in 1-electron move algorithm'

      term=(sqrt(two*pi*tau))**3/pi

! Temporarily create wt_lambda_tau here, but it should be done at start of program
      wt_lambda_tau=wt_lambda**tau

! Undo products
      ipmod=mod(ipass,nfprod)
      ipmod2=mod(ipass+1,nfprod)
      if(idmc.gt.0) then
        if (l_population_control) then
         ginv=min(1.d0,tau)
        else
         ginv=0
        endif
        ffn=eigv*(wdsumo/nconf_global)**ginv
        ffn=ffn**wt_lambda
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

! Undo weights
      iwmod=mod(ipass,nwprod)

! Store (well behaved velocity/velocity)
      if(ncall.eq.0.and.irstar.eq.0) then
        do 7 iw=1,nwalk
          do 7 ifr=1,nforce
! Set nuclear coordinates (0 flag = no strech e-coord)
            if(nforce.gt.1) call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,ifr),ajacob,ifr,0)
            vav2sumo=zero
            v2sumo=zero
            do 6 i=1,nelec
! Find the nearest nucleus
              ren2mn=huge
              do 1 icent=1,ncent
                ren2=(xoldw(1,i,iw,ifr)-cent(1,icent))**2 &
     &              +(xoldw(2,i,iw,ifr)-cent(2,icent))**2 &
     &              +(xoldw(3,i,iw,ifr)-cent(3,icent))**2
                if(ren2.lt.ren2mn) then
                  ren2mn=ren2
                  iwnuc=icent
                endif
    1         continue
              rmino=zero
              do 2 k=1,ndim
                rvmino(k)=xoldw(k,i,iw,ifr)-cent(k,iwnuc)
    2           rmino=rmino+rvmino(k)**2
              rmino=sqrt(rmino)
              zeta=dsqrt(one/tau+znuc(iwctype(iwnuc))**2)

              voldr=zero
              do 3 k=1,ndim
    3           voldr=voldr+voldw(k,i,iw,ifr)*rvmino(k)
              voldr=voldr/rmino


! Tau secondary in drift is one (first time around)
              tratio=one

              hafzr2=(half*znuc(iwctype(iwnuc))*rmino)**2
              v2old=0
              do 5 k=1,ndim
    5           v2old=v2old+voldw(k,i,iw,ifr)**2
              volda=sqrt(v2old)
              adrift=(half*(one+eps+voldr/volda)) &
     &        +adrift0*hafzr2/(one+hafzr2)

              vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/ &
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
! Loop over primary walker
        current_walker = iw !JT
        call object_modified_by_index (current_walker_index) !JT

! Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
        if(nforce.gt.1) call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)

        if(ibasis.eq.3) then             !complex basis set
          call cwalkstrdet(iw)
         else
          call walkstrdet(iw)
        endif
        call walkstrjas(iw)

! Sample Green function for forward move
        r2sume=zero
        risume=zero
        dfus2ac=zero
        dfus2un=zero
        dr2ac=zero
        dr2un=zero
        drifdif=zero
        iaccept=0
        psum=0
        do 200 i=1,nelec
! Find the nearest nucleus & vector from that nucleus to electron
! & component of velocity in that direction
          ren2mn=huge
          do 10 icent=1,ncent
            ren2=0
            do 9 k=1,ndim
   9          ren2=ren2+(xoldw(k,i,iw,1)-cent(k,icent))**2
            if(ren2.lt.ren2mn) then
              ren2mn=ren2
              iwnuc=icent
            endif
   10     continue
          rmino=zero
          voldr=zero
          v2old=zero
          do 20 k=1,ndim
            rvmino(k)=xoldw(k,i,iw,1)-cent(k,iwnuc)
            rmino=rmino+rvmino(k)**2
            voldr=voldr+voldw(k,i,iw,1)*rvmino(k)
   20       v2old=v2old+voldw(k,i,iw,1)**2
          rmino=sqrt(rmino)
          voldr=voldr/rmino
          volda=sqrt(v2old)
          zeta=dsqrt(one/tau+znuc(iwctype(iwnuc))**2)

! Place zaxis along direction from nearest nucleus to electron and
! x-axis along direction of angular component of velocity.
! Calculate the velocity in the phi direction
          voldp=zero
          do 40 k=1,ndim
            zaxis(k)=rvmino(k)/rmino
            xaxis(k)=voldw(k,i,iw,1)-voldr*zaxis(k)
   40       voldp=voldp+xaxis(k)**2
          voldp=sqrt(voldp)
          if(voldp.lt.eps) then
            xaxis(1)=eps*(one-zaxis(1)**2)
            xaxis(2)=eps*(-zaxis(1)*zaxis(2))
            xaxis(3)=eps*(-zaxis(1)*zaxis(3))
            voldp=eps*dsqrt(one+eps-zaxis(1)**2)
          endif
          do 50 k=1,ndim
   50       xaxis(k)=xaxis(k)/voldp

! Use more accurate formula for the drift
          hafzr2=(half*znuc(iwctype(iwnuc))*rmino)**2
          adrift=(half*(1+eps+voldr/volda))+adrift0*hafzr2/(1+hafzr2)

! Tau primary -> tratio=one
          vavvt=(dsqrt(one+two*adrift*v2old*tau)-one)/(adrift*v2old)

          driftr=vavvt*voldr ! driftr can be < 0 because voldr can be < 0
          rtry=rmino+driftr  ! if rtry<0 then electron has drifted past nucleus

! Prob. of sampling exponential rather than gaussian is
! half*derfc(rtry/dsqrt(two*tau)) = half*(two-derfc(-rtry/dsqrt(two*tau)))
! We use both expressions because under AIX the former is faster if rtry>0
! and the latter is faster if rtry<0.
! Note that if adrift is always set to 1 then it may be better to use
! vavvt rather than tau since the max drift distance is dsqrt(2*tau/adrift),
! so both the position of the gaussian and its width are prop to dsqrt(tau)
! if tau is used in derfc, and so qgaus does not tend to 1 for large tau.
! However, we are using a variable adrift that can be very small and then
! using tau is the better choice.

          if(rtry.gt.zero) then
            qgaus=half*derfc(rtry/dsqrt(two*tau))

! Calculate drifted x and y coordinates in local coordinate system centered
! on nearest nucleus
            xprime=vavvt*voldp*rtry/(half*(rmino+rtry))
            zprime=rtry

! Convert back to original coordinate system
            do 60 k=1,ndim
   60         xnew(k)=cent(k,iwnuc)+xaxis(k)*xprime &
     &                             +zaxis(k)*zprime
           else
            qgaus=half*(two-derfc(-rtry/dsqrt(two*tau)))
            rtry=zero
            do 70 k=1,ndim
   70         xnew(k)=cent(k,iwnuc)
          endif
          pgaus=one-qgaus

          if(ipr.ge.1) write(6,'(''xnewdr'',2i4,9f8.5)') iw,i,(xnew(k),k=1,ndim)

! Do the diffusion.  Actually, this diffusion contains part of the drift too
! if idiv_v.ge.1

! Sample gaussian with prob pgaus, exponential with prob. qgaus
          dr2=zero
          dfus2a=zero
          if(rannyu(0).lt.pgaus) then
            dfus2b=zero
            do 80 k=1,ndim
              drift=xnew(k)-xoldw(k,i,iw,1)
              dfus=gauss()*rttau
              dx=drift+dfus
              dr2=dr2+dx**2
              dfus2a=dfus2a+dfus**2
              xnew(k)=xnew(k)+dfus
   80         dfus2b=dfus2b+(xnew(k)-cent(k,iwnuc))**2
            dfusb=sqrt(dfus2b)
           else
            dfusb=(-half/zeta)*dlog(rannyu(0)*rannyu(0)*rannyu(0))
            costht=two*(rannyu(0)-half)
            sintht=sqrt(one-costht*costht)
            phi=two*pi*rannyu(0)
            do 90 k=1,ndim
              drift=xnew(k)-xoldw(k,i,iw,1)
              if(k.eq.1) then
                dfus=dfusb*sintht*cos(phi)
               elseif(k.eq.2) then
                dfus=dfusb*sintht*sin(phi)
               else
                dfus=dfusb*costht
              endif
              dx=drift+dfus
              dr2=dr2+dx**2
              dfus2a=dfus2a+(cent(k,iwnuc)+dfus-xnew(k))**2
   90         xnew(k)=cent(k,iwnuc)+dfus
          endif
          dfus2o=dfus2a
! Radial probability density from atom centered part is 4*zeta^3*dfusb**2*exp(-two*zeta*dfusb), so
! Volume probability density from atom centered part is (zeta^3/pi)*exp(-two*zeta*dfusb).
! Pull out the gaussian since it is the same for the forward and backward moves.
! So, fnormo is what multiplies the gaussian.
          fnormo = pgaus + qgaus*term*(zeta**3)*dexp(-two*zeta*dfusb+half*dfus2a/tau)

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

! Find the nearest nucleus & vector from that nucleus to electron
! & component of velocity in that direction
          ren2mn=huge
          do 110 icent=1,ncent
            ren2=0
            do 100 k=1,ndim
  100           ren2=ren2+(xnew(k)-cent(k,icent))**2
            if(ren2.lt.ren2mn) then
              ren2mn=ren2
              iwnuc=icent
            endif
  110     continue
          rminn=zero
          vnewr=zero
          v2new=zero
          do 120 k=1,ndim
            rvminn(k)=xnew(k)-cent(k,iwnuc)
            rminn=rminn+rvminn(k)**2
            vnewr=vnewr+vnew(k,i)*rvminn(k)
  120       v2new=v2new+vnew(k,i)**2
          rminn=sqrt(rminn)
          vnewr=vnewr/rminn
          vnewa=sqrt(v2new)
          zeta=dsqrt(one/tau+znuc(iwctype(iwnuc))**2)

! Place zaxis along direction from nearest nucleus to electron and
! x-axis along direction of angular component of velocity.
! Calculate the velocity in the phi direction
          vnewp=zero
          do 140 k=1,ndim
            zaxis(k)=rvminn(k)/rminn
            xaxis(k)=vnew(k,i)-vnewr*zaxis(k)
  140       vnewp=vnewp+xaxis(k)**2
          vnewp=sqrt(vnewp)
          if(vnewp.lt.eps) then
            xaxis(1)=eps*(one-zaxis(1)**2)
            xaxis(2)=eps*(-zaxis(1)*zaxis(2))
            xaxis(3)=eps*(-zaxis(1)*zaxis(3))
            vnewp=eps*dsqrt(one+eps-zaxis(1)**2)
          endif
          do 150 k=1,ndim
  150       xaxis(k)=xaxis(k)/vnewp

! Use more accurate formula for the drift
          hafzr2=(half*znuc(iwctype(iwnuc))*rminn)**2
          adrift=(half*(1+eps+vnewr/vnewa))+adrift0*hafzr2/(1+hafzr2)

          vavvt=(dsqrt(one+two*adrift*v2new*tau)-one)/(adrift*v2new)

          driftr=vavvt*vnewr
          rtry=rminn+driftr
          dfus2a=zero
          dfus2b=zero
          if(rtry.gt.zero) then
            qgaus=half*derfc(rtry/dsqrt(two*tau))

! Calculate drifted x and y coordinates in local coordinate system centered
! on nearest nucleus
            xprime=vavvt*vnewp*rtry/(half*(rminn+rtry))
            zprime=rtry

! Convert back to original coordinate system
            do 153 k=1,ndim
              xbac(k)=cent(k,iwnuc)+xaxis(k)*xprime+zaxis(k)*zprime
              dfus2b=dfus2b+(cent(k,iwnuc)-xoldw(k,i,iw,1))**2
              dfus=xbac(k)-xoldw(k,i,iw,1)
  153         dfus2a=dfus2a+dfus**2
           else
            qgaus=half*(two-derfc(-rtry/dsqrt(two*tau)))
            rtry=zero
            do 155 k=1,ndim
              xbac(k)=cent(k,iwnuc)
              dfus=xbac(k)-xoldw(k,i,iw,1)
  155         dfus2b=dfus2b+dfus**2
              dfus2a=dfus2b
          endif
          dfus2n=dfus2a
          pgaus=one-qgaus
          dfusb=sqrt(dfus2b)

          fnormn = pgaus + qgaus*term*(zeta**3)*dexp(-two*zeta*dfusb+half*dfus2a/tau)

          if(ipr.ge.1) then
            write(6,'(''xoldw'',9f10.6)')(xoldw(k,i,iw,1),k=1,ndim), &
     &      (xnew(k),k=1,ndim), (xbac(k),k=1,ndim)
            write(6,'(''dfus2o'',9f10.6)')dfus2o,dfus2n, &
     &      psidow(iw,1),psidn,psijow(iw,1),psijn,fnormo,fnormn
          endif

          p=(psidn/psidow(iw,1))**2*exp(2*(psijn-psijow(iw,1)))* &
     &    exp((dfus2o-dfus2n)/(two*tau))*fnormn/fnormo

          if(ipr.ge.1) write(6,'(''p'',11f10.6)') &
     &    p,(psidn/psidow(iw,1))**2*exp(2*(psijn-psijow(iw,1))), &
     &    exp((dfus2o-dfus2n)/(two*tau)),psidn,psidow(iw,1), &
     &    psijn,psijow(iw,1),dfus2o,dfus2n,fnormo,fnormn


! The following is one reasonable way to cure persistent configurations
! Not needed if itau_eff <=0 and in practice we have never needed it even
! otherwise
          if(iage(iw).gt.50) p=p*1.1d0**(iage(iw)-50)

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
          itryn(i)=int(min(delri*rmino+1,dfloat(NRAD))+eps)

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
  170           voldw(k,l,iw,1)=vnew(k,l)
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

! Calculate moments of r and save rejection probability for primary walk
          r2sume=r2sume+(q*r2o+p*r2n)
          risume=risume+(q/dsqrt(r2o)+p/dsqrt(r2n))
          unacp(i)=q

  200   continue


! Effective tau for branching
        tauprim=tau*dfus2ac/dfus2un
!       tauprim=tau*dr2ac/dr2un
!       write(6,'(''dfus2ac,dfus2un,dr2ac,dr2un,dfus2ac/dfus2un,dr2ac/dr2un,dfus2ac/dfus2un/(dr2ac/dr2un)'',9f8.5)')
!    &  dfus2ac,dfus2un,dr2ac,dr2un,dfus2ac/dfus2un,dr2ac/dr2un,dfus2ac/dfus2un/(dr2ac/dr2un)

        do 280 ifr=1,nforce

          if(iaccept.eq.1) then   ! If any of the 1-electron moves were accepted then
            if(ifr.eq.1) then
! Primary configuration
              drifdifr=one
              if(nforce.gt.1) call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)
              call hpsi(xoldw(1,1,iw,1),psidn,psijn,voldw(1,1,iw,1),div_vow(1,iw),d2n,pen,pein,enew,denergy,1)

              psi_det = psidn                             !JT
              psi_jas = exp(psijn)
              call object_modified_by_index (voldw_index) !JT
              call object_modified_by_index (psi_det_index) !JT
              call object_modified_by_index (psi_jas_index) !JT
              call object_modified_by_index (div_vow_index) !JT

              if(ibasis.eq.3) then                  !complex calculations
                call cwalksav_det(iw)
               else
                call walksav_det(iw)
              endif
              call walksav_jas(iw)
             else
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
              call hpsi(xoldw(1,1,iw,ifr),psidn,psijn,voldw(1,1,iw,ifr),div_vow(1,iw),d2n,pen,pein,enew,denergy,ifr)
            endif

            vav2sumn=zero
            v2sumn=zero
            do 260 i=1,nelec
! Find the nearest nucleus & vector from that nucleus to electron
! & component of velocity in that direction
              ren2mn=huge
              do 230 icent=1,ncent
                ren2=0
                do 225 k=1,ndim
  225             ren2=ren2+(xoldw(k,i,iw,ifr)-cent(k,icent))**2
                if(ren2.lt.ren2mn) then
                  ren2mn=ren2
                  iwnuc=icent
                endif
  230         continue
              rminn=zero
              do 240 k=1,ndim
                rvminn(k)=xoldw(k,i,iw,ifr)-cent(k,iwnuc)
  240           rminn=rminn+rvminn(k)**2
              rminn=sqrt(rminn)
              zeta=dsqrt(one/tau+znuc(iwctype(iwnuc))**2)

              voldr=zero
              do 250 k=1,ndim
  250           voldr=voldr+voldw(k,i,iw,ifr)*rvminn(k)
              voldr=voldr/rminn

! Tau secondary in drift
              tratio=one
              if(ifr.gt.1.and.itausec.eq.1) tratio=drifdifr

! Use more accurate formula for the drift
              hafzr2=(half*znuc(iwctype(iwnuc))*rminn)**2
              v2old=0
              do 255 k=1,ndim
  255           v2old=v2old+voldw(k,i,iw,ifr)**2
              volda=sqrt(v2old)
              adrift=(half*(one+eps+voldr/volda)) &
     &        +adrift0*hafzr2/(one+hafzr2)

              vavvt=(dsqrt(one+two*adrift*v2old*tau*tratio)-one)/ &
     &            (adrift*v2old)
              vavvn=vavvt/(tau*tratio)

              vav2sumn=vav2sumn+vavvn**2*v2old
  260         v2sumn=v2sumn+v2old
            fration=dsqrt(vav2sumn/v2sumn)
           else
            drifdifr=one
            fration=fratio(iw,ifr)
            enew=eoldw(iw,ifr)
          endif


          taunow=tauprim*drifdifr

          if(ipr.ge.1)write(6,'(''wt'',9f10.5)') wt(iw),etrial,eest

!For many years I have been placing an upper limit on the reweighting factor of 1+10*sigma*tau
!(but no lower limit).  This is not in our 1993 paper because I started doing that only after
!we put in the capability to do pseudopotentials, which was after 1993.
!
!Alternatively, one could use for S_bar
!S_bar = min(E_cut,max(-E_cut,E_est-E_L))
!where E_cut = 10*sigma
!which limits both the maximum and the minimum reweighting.
!
!Both of these create a bias in the tau \to 0 limit, but the bias is negigibly small.
!If you object to having any bias, then you could use something like
!E_cut = 3*sigma/sqrt(tau)
!I am assuming that the largest tau one would want to use is about 0.1, so this would
!give a bound of about 10*sigma at tau=.1 and a larger bound at smaller tau values.
!
!All of these have the problem that at the very beginning of the run, one does
!not have a good estimate of sigma.  If that bothers you, you could use
!E_cut = 0.2 sqrt(N_elec/tau)
!which is what Alfe does.  The downside of this is that it assumes that the adhoc value
!of 0.2 is reasonable for all systems.
!
!Another possibility is to multiply (eest-e) by a function of fratio that is 0 at 0, 1 at 1 and deviates from 1 at 1 as some power, e.g.
!UNR93:        f(x)=x                     deviates from 1 linearly.
!new_ene_int:  f(x)=x**2*(3-2*x)          deviates from 1 quadratically.
!new_ene_int2: f(x)=x**3*(10-15*x+6*x**2) deviates from 1 cubically.
!no_ene_int:   f(x)=1                     does not deviate from 1     
!These will give the same energy at tau=0, but progressively lower energies at small finite values of tau

! Warning: Change UNR93 reweighting factor because it gives large time-step error at small tau for pseudo systems as pointed out by Alfe
!         ewto=eest-(eest-eoldw(iw,ifr))*fratio(iw,ifr)                                               ! UNR93
!         ewtn=eest-(eest-enew)*fration                                                               ! UNR93
!         ewto=eoldw(iw,ifr)                                                                          ! no_ene_int
!         ewtn=enew                                                                                   ! no_ene_int
          ewto=eest-(eest-eoldw(iw,ifr))*fratio(iw,ifr)**2*(3-2*fratio(iw,ifr))                       ! new_ene_int
          ewtn=eest-(eest-enew)*fration**2*(3-2*fration)                                              ! new_ene_int
!         ewto=eest-(eest-eoldw(iw,ifr))*fratio(iw,ifr)**3*(10-15*fratio(iw,ifr)+6*fratio(iw,ifr)**2) ! new_ene_int2
!         ewtn=eest-(eest-enew)*fration**3*(10-15*fration+6*fration**2)                               ! new_ene_int2
!         ecut=0.2*sqrt(nelec/tau)                                                                    ! Alfe
!         ewto=etrial+min(ecut,max(-ecut,eoldw(iw,ifr)-eest))                                         ! Alfe
!         ewtn=etrial+min(ecut,max(-ecut,enew-eest))                                                  ! Alfe
! Note in the 2 lines above we use etrial instead of eest for the first term because we do population control with a fixed etrial
! so that we keep track of the fluctuating factor and undo the population control

          do 262 iparm=1,nparm
            dewto(iparm)=denergy_old_dmc(iparm,iw)*fratio(iw,ifr)
  262       dewtn(iparm)=denergy(iparm)*fration

!Warning: tmp
!         write(6,'(''dfus2ac,dfus2un,dr2ac,dr2un,dfus2ac/dfus2un,dr2ac/dr2un,tau,tauprim,taunnow'',99f9.4)') dfus2ac,dfus2un,dr2ac,dr2un,dfus2ac/dfus2un,dr2ac/dr2un,tau,tauprim,taunow

          psum=psum/nelec
          qsum=1-psum

          if(idmc.gt.0) then
! Warning: tmp If instead of reweighting with .5(eold+enew) we use eold we get strong +ve bias at small tau.
            expon=(etrial-half*(ewto+ewtn))*taunow
!           expon=(etrial-ewto)*taunow
!           expon=(etrial-ewtn)*taunow
!           expon=(etrial-((half*psum+qsum)*ewto+half*psum*ewtn))*taunow
! Warning we are temporarily ignoring the term that comes from the derivative of (V_av/V) because
! it should be small compared to the term that we keep.
            do 264 iparm=1,nparm
  264         dexponent(iparm)=-half*(dewto(iparm)+dewtn(iparm))*taunow
            if(icut_br.le.0) then
              dwt=dexp(expon)
             elseif(icut_br.eq.1) then
              if(expon.le.0.d0) then
                dwt=dexp(expon)
               else
! Warning: tmp
!               dwt=1+expon+0.5d0*expon**2
                dwt=1+expon/(1+expon)
              endif
             else
              dwt=0.5d0+1/(1+exp(-4*expon))
            endif
          endif

! Warning: These lines were added to reduce the probability of population explosions.
! These occur mostly for nonlocal psps., and these are cured by our slightly modified version of Casula et al'
! size-consistent tmoves version 1.  So, these lines are no longer needed.
! We truncate wts that come from energies that are too low by more than 10*energy_sigma.
! This gives a DMC energy that is too high even in the tau->0 limit, but by a really negligible amount.
! For mpi1 runs a different energy_sigma is calculated on each processor because I did not want to add new MPI calls.
! For mpi2/3 runs a mpi_allreduce is done in the acues1 routine so that it is summed over the processors.
! So, multiply energy_sigma by sqrt(float(nproc)).
! It is more stable to use the energy_sigma with the population control bias than the one with the bias removed.
!         if(iblk.ge.2. or. (iblk.ge.1 .and. nstep.ge.2)) then
          if(ipass-nstep*2*nblkeq .gt. 5) then
            energy_sigma=e_sigma(ecum1,ecm21,wcum1)
            if(mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') energy_sigma=energy_sigma*sqrt(float(nproc))
            if(dwt.gt.1+10*energy_sigma*tau) then
              ipr_sav=ipr_sav+1
              if(ipr_sav.le.3) then
                write(6,'(''Warning: dwt>1+10*energy_sigma*tau: nwalk,energy_sigma,dwt,ewto,ewtn,fratio(iw,ifr),fration='',i5,9d12.4 &
     &          )') nwalk,energy_sigma,dwt,ewto,ewtn,fratio(iw,ifr),fration
                if(ipr_sav.eq.1) write(6,'(''This should add a totally negligible positive bias to the energy'')')
               elseif(ipr_sav.eq.4) then
                write(6,'(''Warning: Additional warning msgs. of dwt>1+10*energy_sigma*tau suppressed'')')
              endif
              dwt=1+10*energy_sigma*tau
            endif
          endif

! ffi has already been raised to wt_lambda.  Do the same for dwt.  We do this even for the current move so that wt_lambda can serve to limit size of move.
          dwt=dwt**wt_lambda

! Exercise population control if dmc or vmc with weights
          if(idmc.gt.0.or.iacc_rej.eq.0) dwt=dwt*ffi

! Set weights and product of weights over last nwprod steps
          if(ifr.eq.1) then

! Raise product of previous generation wts to power wt_lambda_tau to keep product under control if branching is turned off
            wt(iw)=(wt(iw)**wt_lambda_tau)*dwt
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

          if(ipr.ge.1)write(6,'(''eoldw,enew,wt'',9f10.5)') &
     &    eoldw(iw,ifr),enew,wtnow

          wtg=wtnow*fprod
          current_walker_weight = wt(iw) * fprod !JT
          call object_modified_by_index (current_walker_weight_index) !JT

          if(ifr.eq.1) then

            r2sum=r2sum+wtg*r2sume
            risum=risum+wtg*risume
            do 270 i=1,nelec
              if(i.le.nup) then
                rprobup(itryo(i))=rprobup(itryo(i))+wtg*unacp(i)
                rprobup(itryn(i))=rprobup(itryn(i))+wtg*(one-unacp(i))
               else
                rprobdn(itryo(i))=rprobdn(itryo(i))+wtg*unacp(i)
                rprobdn(itryn(i))=rprobdn(itryn(i))+wtg*(one-unacp(i))
              endif
              rprob(itryo(i))=rprob(itryo(i))+wtg*unacp(i)
  270         rprob(itryn(i))=rprob(itryn(i))+wtg*(one-unacp(i))
!             r2sum=r2sum+wtg*(unacp(i)*r2o+(one-unacp(i)*r2n)
! 270         risum=risum+wtg*(unacp(i)/dsqrt(r2o)+(one-unacp(i)/dsqrt(r2n))

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
!JT         call object_provide('ovlp_ovlp_fn')
!JT         ovlp_ovlp_fn_sum = ovlp_ovlp_fn_sum + ovlp_ovlp_fn

!           local energy for current walker
            eloc_dmc = eoldw(iw,1)
            call object_modified_by_index (eloc_dmc_index)

           else
            ro=ajacold(iw,ifr)*psidow(iw,ifr)**2* &
     &         exp(2*psijow(iw,ifr)-psi2savo)

            wsum1(ifr)=wsum1(ifr)+wtnow*ro
            esum1(ifr)=esum1(ifr)+wtnow*eoldw(iw,ifr)*ro
            pesum(ifr)=pesum(ifr)+wtg*peow(iw,ifr)*ro
            peisum(ifr)=peisum(ifr)+wtg*peiow(iw,ifr)*ro
            tpbsum(ifr)=tpbsum(ifr)+wtg*(eoldw(iw,ifr)-peow(iw,ifr))*ro
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*d2ow(iw,ifr)*ro
          endif

  280   continue
! Call to grad_hess_jas_sum() used to be for optimizing Jastrow for periodic systems.
        call grad_hess_jas_sum(1.d0,0.d0,eoldw(iw,1),eoldw(iw,1),wt(iw)*fprod,wi_w(:,iw))

        call compute_averages_step !JT
  300 continue

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
