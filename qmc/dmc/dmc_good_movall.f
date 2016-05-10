      subroutine dmc_good_movall
! Written by Cyrus Umrigar
! Uses the diffusion Monte Carlo algorithm described in:
! 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
!    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993).
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
! 2 1 1 1 0 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
! 2 1 0 1 1 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
! Another reasonable choice is:
! 2 1 0 1 1 1 1 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      use constants_mod
      use atom_mod
      use dets_mod
      use contrl_mod
      use average_mod
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
      use pop_control_mod, only : ffn
      implicit real*8(a-h,o-z)

      parameter (eps=1.d-10,huge=1.d+100,adrift0=0.1d0)

      common /tmp/ eacc,enacc,macc,mnacc

      dimension rmino(nelec),rminn(nelec),rvmin(3)
      dimension xnew(3,nelec,nforce),vnew(3,nelec,nforce),psidn(nforce),
     &psijn(nforce),enew(nforce),pen(nforce),pein(nforce),d2n(nforce),div_vn(nelec)
      dimension xaxis(3),zaxis(3),xbac(3),ajacnew(nforce)


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
        current_walker = iw !JT CU
        call object_modified_by_index (current_walker_index) !JT CU

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
! Find the nearest nucleus & vector from that nucleus to electron
! & component of velocity in that direction
            ren2mn=huge
            do 10 icent=1,ncent
              ren2=(xoldw(1,i,iw,ifr)-cent(1,icent))**2
     &            +(xoldw(2,i,iw,ifr)-cent(2,icent))**2
     &            +(xoldw(3,i,iw,ifr)-cent(3,icent))**2
              if(ren2.lt.ren2mn) then
                ren2mn=ren2
                iwnuc=icent
              endif
   10       continue
            rmino(i)=zero
            voldr=zero
            v2old=zero
            do 20 k=1,ndim
              rvmin(k)=xoldw(k,i,iw,ifr)-cent(k,iwnuc)
              rmino(i)=rmino(i)+rvmin(k)**2
              voldr=voldr+voldw(k,i,iw,ifr)*rvmin(k)
   20         v2old=v2old+voldw(k,i,iw,ifr)**2
            rmino(i)=sqrt(rmino(i))
            voldr=voldr/rmino(i)
            volda=sqrt(v2old)

! Place zaxis along direction from nearest nucleus to electron and
! x-axis along direction of angular component of velocity.
! Calculate the velocity in the phi direction
            voldp=zero
            do 40 k=1,ndim
              zaxis(k)=rvmin(k)/rmino(i)
              xaxis(k)=voldw(k,i,iw,ifr)-voldr*zaxis(k)
   40         voldp=voldp+xaxis(k)**2
            voldp=sqrt(voldp)
            if(voldp.lt.eps) then
              xaxis(1)=eps*(one-zaxis(1)**2)
              xaxis(2)=eps*(-zaxis(1)*zaxis(2))
              xaxis(3)=eps*(-zaxis(1)*zaxis(3))
              voldp=eps*dsqrt(one+eps-zaxis(1)**2)
            endif
            do 50 k=1,ndim
   50         xaxis(k)=xaxis(k)/voldp

! Use more accurate formula for the drift
            hafzr2=(half*znuc(iwctype(iwnuc))*rmino(i))**2
            adrift=(half*(1+eps+voldr/volda))+adrift0*hafzr2/(1+hafzr2)

! Tau secondary in drift in order to compute tau_s in second equil. block
            if(ifr.gt.1.and.itausec.eq.1) then
              if(itau_eff.ge.1) then
                tauu=tau*taueff(ifr)/taueff(1)
               else
                tauu=taueff(ifr)
              endif
             else
              tauu=tau
            endif

! Note: vavvo=vav2sumo/v2sumo appears in the branching
            vav=(dsqrt(1+2*adrift*v2old*tauu)-1)/(adrift*volda*tauu)
            vav2sumo=vav2sumo+vav**2
            v2sumo=v2sumo+v2old

            driftr=voldr*tauu*vav/volda
            rtry=rmino(i)+driftr
            rtrya=max(0.d0,rtry)
            rtrya2=rtrya**2

! Prob. of sampling exponential rather than gaussian is
! half*derfc(rtry/dsqrt(two*tau)) = half*(two-derfc(-rtry/dsqrt(two*tau)))
! Note that if adrift is always set to 1 then it may be better to use
! vavvt rather than tau since the max drift distance is dsqrt(2*tau/adrift),
! so both the position of the gaussian and its width are prop to dsqrt(tau)
! if tau is used in derfc, and so qgaus does not tend to 1 for large tau.
! However, we are using a variable adrift that can be very small and then
! using tau is the better choice.
            qgaus=half*derfc(rtry/dsqrt(two*tau))
            pgaus=1-qgaus

! Calculate drifted x and y coordinates in local coordinate system centered
! on nearest nucleus
            xprime=(vav/volda)*voldp*tauu*rtrya/(half*(rmino(i)+rtry))
            zprime=rtrya

! Convert back to original coordinate system
            do 60 k=1,ndim
              xbac(k)=cent(k,iwnuc)+xaxis(k)*xprime+zaxis(k)*zprime
              if(ifr.gt.1) then
                dfus=xnew(k,i,ifr)-xbac(k)
                dfus2o=dfus2o+dfus**2
              endif
   60       continue

            if(ipr.ge.1) then
              write(6,'(''vav,volda,driftr,rtry'',9f9.5)') vav,volda,driftr,rtry,xprime,zprime
              write(6,'(''xnewdr'',2i4,9f8.5)') iw,i,(xbac(k),k=1,ndim)
              write(6,'(''xaxis'',9f9.6)') (xaxis(k),k=1,ndim)
              write(6,'(''zaxis'',9f9.6)') (zaxis(k),k=1,ndim)
            endif

! Do the diffusion.  Actually, this diffusion contains part of the drift too
! if idiv_v.ge.1
! Use div_v_hom to set tau_hom
            div_v_hom=div_vow(i,iw)/3

            if(ifr.eq.1) then
              if(idiv_v.ge.1) then
                if(div_v_hom.lt.0.d0) then
                  tau_hom=3*tau - 1/div_v_hom - ((tau-1/div_v_hom)*sqrt(-div_v_hom) * (1-derfc(1/sqrt(-2*div_v_hom*tau)))
     &            + sqrt(2*tau/pi)*exp(1/(2*div_v_hom*tau)))**2
                 else
                  tau_hom=tau*(1+2*div_v_hom*tau)/(1+div_v_hom*tau)
                endif
               else
                tau_hom=tau
              endif
              rttau_hom=dsqrt(tau_hom)

! First set zeta to approx. value and then call zeta_cusp to set it to
! value that imposes the cusp condition on G.
! Warning I probably need to put zeta in the ifr.ne.1 part too
              zeta=dsqrt(one/tau+znuc(iwctype(iwnuc))**2)
              if(icuspg.ge.1) then
                zeta=znuc(iwctype(iwnuc)) * (1+derfc(-rtry/dsqrt(2*tau_hom))/(200*tau_hom))
!               zeta=sqrt(znuc(iwctype(iwnuc))**2 + 0.2*derfc(-rtry/dsqrt(2*tau_hom))/tau_hom)
                term2=pgaus*pi*znuc(iwctype(iwnuc))*exp(-rtrya2/tau_hom)/(qgaus*sqrt(2*pi*tau_hom)**3)
                if(qgaus.gt.eps) zeta=zeta_cusp(zeta,znuc(iwctype(iwnuc)),term2)
              endif

! Check cusp
!             qgaus2=2*znuc(iwctype(iwnuc))*exp(-(rtrya2/(2*tau_hom)))
!    &        /(2*pi*tau_hom)**1.5/
!    &        (2*zeta**4/pi+2*znuc(iwctype(iwnuc))*
!    &        (exp(-(rtrya2/(2*tau_hom)))/(2*pi*tau_hom)**1.5
!    &        -zeta**3/pi))
!             write(6,'(''o pgaus,qgaus='',4d9.2)') pgaus,qgaus,qgaus2

! Sample gaussian with prob pgaus, exponential with prob. qgaus
              dfus2a=0
              dfus2b=0
              if(rannyu(0).lt.pgaus) then
                do 80 k=1,ndim
                  dfus=gauss()*rttau_hom
                  xnew(k,i,ifr)=xbac(k)+dfus
                  dr2=dr2+(xnew(k,i,ifr)-xoldw(k,i,iw,ifr))**2
                  dfus2a=dfus2a+dfus**2
                  dfus2o=dfus2o+dfus**2
   80             dfus2b=dfus2b+(xnew(k,i,ifr)-cent(k,iwnuc))**2
                dfusb=sqrt(dfus2b)
               else
                dfusb=(-half/zeta)*dlog(rannyu(0)*rannyu(0)*rannyu(0))
                costht=two*(rannyu(0)-half)
                sintht=sqrt(one-costht*costht)
                phi=two*pi*rannyu(0)
                do 90 k=1,ndim
                  drift=xbac(k)-xoldw(k,i,iw,ifr)
                  if(k.eq.1) then
                    dfus=dfusb*sintht*cos(phi)
                   elseif(k.eq.2) then
                    dfus=dfusb*sintht*sin(phi)
                   else
                    dfus=dfusb*costht
                  endif
                  dx=drift+dfus
                  dr2=dr2+dx**2
                  xnew(k,i,ifr)=cent(k,iwnuc)+dfus
                  dfus2a=dfus2a+(xnew(k,i,ifr)-xbac(k))**2
   90             dfus2o=dfus2o+(xnew(k,i,ifr)-xbac(k))**2
              endif
            endif
            gaus_norm=1/sqrt(two*pi*tau_hom)**3
  100       fnormo=fnormo * (pgaus*gaus_norm*exp(-half*dfus2a/tau_hom) + qgaus*(zeta**3/pi)*dexp(-two*zeta*dfusb))
          vavvo=sqrt(vav2sumo/v2sumo)

! calculate psi etc. at new configuration
          call hpsi(xnew(1,1,ifr),psidn(ifr),psijn(ifr),vnew(1,1,ifr),div_vn,d2n(ifr),pen(ifr),pein(ifr),enew(ifr),denergy,ifr)

          psi_det = psidn(1)                          !JT CU
          psi_jas = exp(psijn(1))                     !JT CU
          call object_modified_by_index (voldw_index) !JT CU
          call object_modified_by_index (psi_det_index) !JT CU
          call object_modified_by_index (psi_jas_index) !JT CU
          call object_modified_by_index (div_vow_index) !JT CU

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
! Find the nearest nucleus & vector from that nucleus to electron
! & component of velocity in that direction
            ren2mn=huge
            do 110 icent=1,ncent
              ren2=(xnew(1,i,ifr)-cent(1,icent))**2
     &            +(xnew(2,i,ifr)-cent(2,icent))**2
     &            +(xnew(3,i,ifr)-cent(3,icent))**2
              if(ren2.lt.ren2mn) then
                ren2mn=ren2
                iwnuc=icent
              endif
  110       continue
            rminn(i)=zero
            vnewr=zero
            v2new=zero
            do 120 k=1,ndim
              rvmin(k)=xnew(k,i,ifr)-cent(k,iwnuc)
              rminn(i)=rminn(i)+rvmin(k)**2
              vnewr=vnewr+vnew(k,i,ifr)*rvmin(k)
  120         v2new=v2new+vnew(k,i,ifr)**2
            rminn(i)=sqrt(rminn(i))
            vnewr=vnewr/rminn(i)
            vnewa=sqrt(v2new)

! Place zaxis along direction from nearest nucleus to electron and
! x-axis along direction of angular component of velocity.
! Calculate the velocity in the phi direction
            vnewp=zero
            do 140 k=1,ndim
              zaxis(k)=rvmin(k)/rminn(i)
              xaxis(k)=vnew(k,i,ifr)-vnewr*zaxis(k)
  140         vnewp=vnewp+xaxis(k)**2
            vnewp=sqrt(vnewp)
            if(vnewp.lt.eps) then
              xaxis(1)=eps*(one-zaxis(1)**2)
              xaxis(2)=eps*(-zaxis(1)*zaxis(2))
              xaxis(3)=eps*(-zaxis(1)*zaxis(3))
              vnewp=eps*dsqrt(one+eps-zaxis(1)**2)
            endif
            do 150 k=1,ndim
  150         xaxis(k)=xaxis(k)/vnewp

! Use more accurate formula for the drift
            hafzr2=(half*znuc(iwctype(iwnuc))*rminn(i))**2
            adrift=(half*(1+eps+vnewr/vnewa))+adrift0*hafzr2/(1+hafzr2)

            vav=(dsqrt(1+2*adrift*v2new*tauu)-1)/(adrift*vnewa*tauu)
            vav2sumn=vav2sumn+vav**2
            v2sumn=v2sumn+v2new

            driftr=vnewr*tauu*vav/vnewa
            rtry=rminn(i)+driftr
            rtrya=max(0.d0,rtry)
            rtrya2=rtrya**2

! Prob. of sampling exponential rather than gaussian is
! half*derfc(rtry/dsqrt(two*tau)) = half*(two-derfc(-rtry/dsqrt(two*tau)))
            qgaus=half*derfc(rtry/dsqrt(two*tau))
            pgaus=1-qgaus

! Calculate drifted x and y coordinates in local coordinate system centered
! on nearest nucleus
            xprime=(vav/vnewa)*vnewp*tauu*rtrya/(half*(rminn(i)+rtry))
            zprime=rtrya

! Convert back to original coordinate system
            dfus2a=0
            dfus2b=0
            do 160 k=1,ndim
              xbac(k)=cent(k,iwnuc)+xaxis(k)*xprime+zaxis(k)*zprime
              dfus2a=dfus2a+(xoldw(k,i,iw,ifr)-xbac(k))**2
  160         dfus2b=dfus2b+(xoldw(k,i,iw,ifr)-cent(k,iwnuc))**2
            dfusb=sqrt(dfus2b)

! Use div_v_hom to set tau_hom
            div_v_hom=div_vn(i)/3

            if(idiv_v.ge.1) then
              if(div_v_hom.lt.0.d0) then
                tau_hom=3*tau - 1/div_v_hom - ((tau-1/div_v_hom)*sqrt(-div_v_hom) * (1-derfc(1/sqrt(-2*div_v_hom*tau)))
     &          + sqrt(2*tau/pi)*exp(1/(2*div_v_hom*tau)))**2
               else
                tau_hom=tau*(1+2*div_v_hom*tau)/(1+div_v_hom*tau)
              endif
             else
              tau_hom=tau
            endif

! First set zeta to approx. value and then call zeta_cusp to set it to
! value that imposes the cusp condition on G.
            zeta=dsqrt(one/tau+znuc(iwctype(iwnuc))**2)
            if(icuspg.ge.1) then
              zeta=znuc(iwctype(iwnuc)) * (1+derfc(-rtry/dsqrt(2*tau_hom))/(200*tau_hom))
!             zeta=sqrt(znuc(iwctype(iwnuc))**2 + 0.2*derfc(-rtry/dsqrt(2*tau_hom))/tau_hom)
              term2=pgaus*pi*znuc(iwctype(iwnuc))*exp(-rtrya2/tau_hom) / (qgaus*sqrt(2*pi*tau_hom)**3)
              if(qgaus.gt.1.d-10) zeta=zeta_cusp(zeta,znuc(iwctype(iwnuc)),term2)
            endif

! Check cusp
!           qgaus2=2*znuc(iwctype(iwnuc))*exp(-(rtrya2**2/(2*tau_hom)))/
!    &      (2*pi*tau_hom)**1.5/
!    &      (2*zeta**4/pi+2*znuc(iwctype(iwnuc))*
!    &      (exp(-(rtrya2**2/(2*tau_hom)))/(2*pi*tau_hom)**1.5
!    &      -zeta**3/pi))
!           write(6,'(''n pgaus,qgaus='',4d9.2)') pgaus,qgaus,qgaus2

            gaus_norm=1/sqrt(two*pi*tau_hom)**3
            fnormn=fnormn*(pgaus*gaus_norm*exp(-half*dfus2a/tau_hom) + qgaus*(zeta**3/pi)*dexp(-two*zeta*dfusb))

            if(ipr.ge.1) then
              write(6,'(''xoldw'',9f10.6)')(xoldw(k,i,iw,ifr),k=1,ndim), (xnew(k,i,ifr),k=1,ndim), (xbac(k),k=1,ndim)
              write(6,'(''dfus2o'',9f10.6)')dfus2o,dfus2n,psidow(iw,ifr), psidn(ifr),psijow(iw,ifr),psijn(ifr),fnormo,fnormn
            endif

  200     continue
          vavvn=sqrt(vav2sumn/v2sumn)

          p=(psidn(ifr)/psidow(iw,ifr))**2*exp(2*(psijn(ifr) - psijow(iw,ifr)))*fnormn/fnormo

          if(ipr.ge.1) write(6,'(''p'',9f10.6)')
     &    p, (psidn(ifr)/psidow(iw,ifr))**2*exp(2*(psijn(ifr)-psijow(iw,ifr))), exp((dfus2o-dfus2n)/(two*tau)),psidn(ifr),
     &    psidow(iw,ifr), psijn(ifr), psijow(iw,ifr), dfus2o, dfus2n

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
!           expon=taunow*min(max((etrial-half*(eoldw(iw,ifr)+enew)),-.2*sqrt(nelec/tau)),.2*sqrt(nelec/tau)) ! Alfe
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

          if(ipr.ge.1)write(6,'(''eoldw,enew,wt'',9f10.5)') eoldw(iw,ifr),enew(ifr),wtnow

          wtg=wtnow*fprod

          current_walker_weight = wt(iw) * fprod !JT CU
          call object_modified_by_index (current_walker_weight_index) !JT CU

          enowo=eoldw(iw,ifr)
          enown=enew(ifr)

          if(ifr.eq.1) then
            psi2savo=2*(psijow(iw,1)+dlog(dabs(psidow(iw,1))))
            psi2savn=2*(psijn(1)+dlog(dabs(psidn(1))))

            wsum1(ifr)=wsum1(ifr)+wtnow
            esum1(ifr)=esum1(ifr)+wtnow*(q*enowo+p*enown)
            pesum(ifr)=pesum(ifr)+  wtg*(q*peow(iw,ifr) +p*pen(ifr))
            peisum(ifr)=peisum(ifr)+  wtg*(q*peiow(iw,ifr) +p*pein(ifr))
            tpbsum(ifr)=tpbsum(ifr)+wtg*(q*(enowo-peow(iw,ifr)) + p*(enown-pen(ifr)))
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*(q*d2ow(iw,ifr) + p*d2n(ifr))

!           local energy for current walker
!           Warning: why is this here?
            eloc_dmc = eoldw(iw,1) ! JT CU
            call object_modified_by_index (eloc_dmc_index) ! JT CU

           else

            ro=ajacold(iw,ifr)*psidow(iw,ifr)**2*exp(2*psijow(iw,ifr)-psi2savo)
            rn=ajacnew(ifr)*psidn(ifr)**2*exp(2*psijn(ifr)-psi2savn)

            wsum1(ifr)=wsum1(ifr) + wtnow*(qsav*ro + psav*rn)
            esum1(ifr)=esum1(ifr) + wtnow*(qsav*enowo*ro + psav*enown*rn)
            pesum(ifr)=pesum(ifr) + wtg*(qsav*peow(iw,ifr)*ro + psav*pen(ifr)*rn)
            peisum(ifr)=peisum(ifr) + wtg*(qsav*peiow(iw,ifr)*ro + psav*pein(ifr)*rn)
            tpbsum(ifr)=tpbsum(ifr) + wtg*(qsav*(enowo-peow(iw,ifr))*ro + psav*(enown-pen(ifr))*rn)
            tjfsum(ifr)=tjfsum(ifr) - wtg*half*hb*(qsav*d2ow(iw,ifr)*ro + psav*d2n(ifr)*rn)
          endif

! Collect density only for primary walk
          if(ifr.eq.1) then
            do 250 i=1,nelec
              r2o=zero
              r2n=zero
              do 240 k=1,ndim
                r2n=r2n+xnew(k,i,ifr)**2
  240           r2o=r2o+xoldw(k,i,iw,ifr)**2
!             itryo=min(int(delri*rmino(i))+1,NRAD)
!             itryn=min(int(delri*rminn(i))+1,NRAD)
              itryo=int(min(delri*rmino(i)+1,dfloat(NRAD))+eps)
              itryn=int(min(delri*rminn(i)+1,dfloat(NRAD))+eps)
              if(i.le.nup) then
                rprobup(itryo)=rprobup(itryo)+wtg*q
                rprobup(itryn)=rprobup(itryn)+wtg*p
               else
                rprobdn(itryo)=rprobdn(itryo)+wtg*q
                rprobdn(itryn)=rprobdn(itryn)+wtg*p
              endif
              rprob(itryo)=rprob(itryo)+wtg*q
              rprob(itryn)=rprob(itryn)+wtg*p
              r2sum=r2sum+wtg*(q*r2o+p*r2n)
  250         risum=risum+wtg*(q/dsqrt(r2o)+p/dsqrt(r2n))
          endif

          if(iaccept.eq.1) then

!           if(dabs((enew(ifr)-etrial)/etrial).gt.0.2d+0) then
            if(dabs((enew(ifr)-etrial)/etrial).gt.5.0d+0) then
               write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass, enew(ifr)-etrial,psidn(ifr),psijn(ifr),
     &         ((xnew(k,jj,ifr),k=1,ndim),jj=1,nelec)
            endif

            if(wt(iw).gt.3) write(18,'(i6,i4,3f8.2,30f8.4)') ipass,iw, wt(iw), enew(ifr)-etrial, eoldw(iw,ifr)-etrial,
     &      ((xnew(k,jj,ifr),k=1,ndim),jj=1,nelec)

            psidow(iw,ifr)=psidn(ifr)
            psijow(iw,ifr)=psijn(ifr)
            eoldw(iw,ifr)=enew(ifr)
            peow(iw,ifr)=pen(ifr)
            peiow(iw,ifr)=pein(ifr)
            d2ow(iw,ifr)=d2n(ifr)
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

        call compute_averages_step !JT CU
  300 continue

      if(wsum1(1).gt.1.1*nconf_global) write(18,'(i6,9d12.4)') ipass, ffn, fprod, fprodd, wsum1(1), wgdsumo

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
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv, eest, egcum(1), egsum(1), wgcum(1), wgsum(1), fprod
       else
        accavn=acc_int/try_int
        if(icut_br.le.0) then
          eigv=accavn*exp((etrial-eest)*tau)+(one-accavn)
         else
          eigv=one-accavn+accavn*(0.5d0+1/(1+exp(-4*(etrial-eest)*tau)))
        endif
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest,accavn, egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wdsumn
      wgdsumo=wgdsumn
      wtgen(ipmod)=wdsumn

!     Warning: why are these here?
      call object_modified_by_index (eoldw_index) !JT CU
      call object_modified_by_index (wt_index) !JT CU
      call object_modified_by_index (fprod_index) !JT CU

      return
      end
