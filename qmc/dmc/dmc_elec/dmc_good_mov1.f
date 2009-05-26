      subroutine dmc_good_mov1
c Written by Cyrus Umrigar and Claudia Filippi
c Uses the diffusion Monte Carlo algorithm described in:
c 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
c    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993)
c modified to do accept/reject after single-electron moves.
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
      use control_mod
      use average_mod
      use atom_mod
      use dets_mod
      use basis1_mod
      use contrl_mod
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
      use div_v_dmc_mod
      use contrldmc_mod
      use stats_mod
      use age_mod
      implicit real*8(a-h,o-z)

!JT      parameter (zero=0.d0,one=1.d0,two=2.d0,half=.5d0)
      parameter (eps=1.d-10,huge=1.d+100,adrift0=0.1d0)

      common /velratio/ fratio(MWALK,MFORCE)
      common /branch_dmc_opt/ denergy_old_dmc(MPARM,MWALK),wi_w(MPARM,MWALK)
!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere

      dimension rvmino(3),rvminn(3),xstrech(3,MELEC)
      dimension xnew(3),vnew(3,MELEC)
      dimension xaxis(3),zaxis(3),xbac(3)
      dimension itryo(MELEC),itryn(MELEC),unacp(MELEC)
      dimension dewto(MPARM),dewtn(MPARM),dexponent(MPARM)

      data ncall,ipr_sav /0,0/
      save ipr_sav

c     gauss()=dcos(two*pi*rannyu(0))*sqrt(-two*dlog(rannyu(0)))


      if(icuspg.ge.1) stop 'exact cusp in G not implemented yet in 1-el
     &ectron move algorithm'
      if(idiv_v.ge.1) stop 'div_v not implemented yet in 1-electron move
     &algorithm'

      term=(sqrt(two*pi*tau))**3/pi

c Undo products
      ipmod=mod(ipass,nfprod)
      ipmod2=mod(ipass+1,nfprod)
      if(idmc.gt.0) then
        ginv=min(1.d0,tau)
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
c Set nuclear coordinates (0 flag = no strech e-coord)
            if(nforce.gt.1)
     &      call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,ifr),ajacob,ifr,0)
            vav2sumo=zero
            v2sumo=zero
            do 6 i=1,nelec
c Find the nearest nucleus
              ren2mn=huge
              do 1 icent=1,ncent
                ren2=(xoldw(1,i,iw,ifr)-cent(1,icent))**2
     &              +(xoldw(2,i,iw,ifr)-cent(2,icent))**2
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


c Tau secondary in drift is one (first time around)
              tratio=one

              hafzr2=(half*znuc(iwctype(iwnuc))*rmino)**2
              v2old=0
              do 5 k=1,ndim
    5           v2old=v2old+voldw(k,i,iw,ifr)**2
              volda=sqrt(v2old)
              adrift=(half*(one+eps+voldr/volda))
     &        +adrift0*hafzr2/(one+hafzr2)

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

c Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
        if(nforce.gt.1)
     &  call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)

        if(ibasis.eq.3) then             !complex basis set
          call cwalkstrdet(iw)
         else
          call walkstrdet(iw)
        endif
        call walkstrjas(iw)

c Sample Green function for forward move
        r2sume=zero
        risume=zero
        dfus2ac=zero
        dfus2un=zero
        dr2ac=zero
        dr2un=zero
        drifdif=zero
        iaccept=0
        do 200 i=1,nelec
c Find the nearest nucleus & vector from that nucleus to electron
c & component of velocity in that direction
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

c Place zaxis along direction from nearest nucleus to electron and
c x-axis along direction of angular component of velocity.
c Calculate the velocity in the phi direction
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

c Use more accurate formula for the drift
          hafzr2=(half*znuc(iwctype(iwnuc))*rmino)**2
          adrift=(half*(1+eps+voldr/volda))+adrift0*hafzr2/(1+hafzr2)

c Tau primary -> tratio=one
          vavvt=(dsqrt(one+two*adrift*v2old*tau)-one)/(adrift*v2old)

          driftr=vavvt*voldr
          rtry=rmino+driftr

c Prob. of sampling exponential rather than gaussian is
c half*derfc(rtry/dsqrt(two*tau)) = half*(two-derfc(-rtry/dsqrt(two*tau)))
c We use both expressions because under AIX the former is faster if rtry>0
c and the latter is faster if rtry<0.
c Note that if adrift is always set to 1 then it may be better to use
c vavvt rather than tau since the max drift distance is dsqrt(2*tau/adrift),
c so both the position of the gaussian and its width are prop to dsqrt(tau)
c if tau is used in derfc, and so qgaus does not tend to 1 for large tau.
c However, we are using a variable adrift that can be very small and then
c using tau is the better choice.

          if(rtry.gt.zero) then
            qgaus=half*derfc(rtry/dsqrt(two*tau))

c Calculate drifted x and y coordinates in local coordinate system centered
c on nearest nucleus
            xprime=vavvt*voldp*rtry/(half*(rmino+rtry))
            zprime=rtry

c Convert back to original coordinate system
            do 60 k=1,ndim
   60         xnew(k)=cent(k,iwnuc)+xaxis(k)*xprime
     &                             +zaxis(k)*zprime
           else
            qgaus=half*(two-derfc(-rtry/dsqrt(two*tau)))
            rtry=zero
            do 70 k=1,ndim
   70         xnew(k)=cent(k,iwnuc)
          endif
          pgaus=one-qgaus

          if(ipr.ge.1)
     &    write(6,'(''xnewdr'',2i4,9f8.5)') iw,i,(xnew(k),k=1,ndim)

c Do the diffusion.  Actually, this diffusion contains part of the drift too
c if idiv_v.ge.1

c Sample gaussian with prob pgaus, exponential with prob. qgaus
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
          fnormo=(pgaus+qgaus*term*(zeta**3)*
     &           dexp(-two*zeta*dfusb+half*dfus2a/tau))

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

c Find the nearest nucleus & vector from that nucleus to electron
c & component of velocity in that direction
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

c Place zaxis along direction from nearest nucleus to electron and
c x-axis along direction of angular component of velocity.
c Calculate the velocity in the phi direction
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

c Use more accurate formula for the drift
          hafzr2=(half*znuc(iwctype(iwnuc))*rminn)**2
          adrift=(half*(1+eps+vnewr/vnewa))+adrift0*hafzr2/(1+hafzr2)

          vavvt=(dsqrt(one+two*adrift*v2new*tau)-one)/(adrift*v2new)

          driftr=vavvt*vnewr
          rtry=rminn+driftr
          dfus2a=zero
          dfus2b=zero
          if(rtry.gt.zero) then
            qgaus=half*derfc(rtry/dsqrt(two*tau))

c Calculate drifted x and y coordinates in local coordinate system centered
c on nearest nucleus
            xprime=vavvt*vnewp*rtry/(half*(rminn+rtry))
            zprime=rtry

c Convert back to original coordinate system
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

          fnormn=pgaus+qgaus*term*(zeta**3)*
     &    dexp(-two*zeta*dfusb+half*dfus2a/tau)

          if(ipr.ge.1) then
            write(6,'(''xoldw'',9f10.6)')(xoldw(k,i,iw,1),k=1,ndim),
     &      (xnew(k),k=1,ndim), (xbac(k),k=1,ndim)
            write(6,'(''dfus2o'',9f10.6)')dfus2o,dfus2n,
     &      psidow(iw,1),psidn,psijow(iw,1),psijn,fnormo,fnormn
          endif

          p=(psidn/psidow(iw,1))**2*exp(2*(psijn-psijow(iw,1)))*
     &    exp((dfus2o-dfus2n)/(two*tau))*fnormn/fnormo

          if(ipr.ge.1) write(6,'(''p'',11f10.6)')
     &    p,(psidn/psidow(iw,1))**2*exp(2*(psijn-psijow(iw,1))),
     &    exp((dfus2o-dfus2n)/(two*tau)),psidn,psidow(iw,1),
     &    psijn,psijow(iw,1),dfus2o,dfus2n,fnormo,fnormn


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

c Calculate moments of r and save rejection probability for primary walk
          r2sume=r2sume+(q*r2o+p*r2n)
          risume=risume+(q/dsqrt(r2o)+p/dsqrt(r2n))
          unacp(i)=q

  200   continue


c Effective tau for branching
        tauprim=tau*dfus2ac/dfus2un
c       tauprim=tau*dr2ac/dr2un
c       write(6,'(''dfus2ac,dfus2un,dr2ac,dr2un,dfus2ac/dfus2un,dr2ac/dr2un,dfus2ac/dfus2un/(dr2ac/dr2un)'',9f8.5)')
c    &  dfus2ac,dfus2un,dr2ac,dr2un,dfus2ac/dfus2un,dr2ac/dr2un,dfus2ac/dfus2un/(dr2ac/dr2un)

        do 280 ifr=1,nforce

          if(iaccept.eq.1) then   ! If any of the 1-electron moves were accepted then
            if(ifr.eq.1) then
c Primary configuration
              drifdifr=one
              if(nforce.gt.1)
     &        call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)
              call hpsi(xoldw(1,1,iw,1),psidn,psijn,voldw(1,1,iw,1),div_vow(1,iw),d2n,pen,pein,enew,denergy,1)

              call object_modified_by_index (voldw_index) !JT
              call object_modified_by_index (div_vow_index) !JT

              if(ibasis.eq.3) then                  !complex calculations
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
            endif

            vav2sumn=zero
            v2sumn=zero
            do 260 i=1,nelec
c Find the nearest nucleus & vector from that nucleus to electron
c & component of velocity in that direction
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

c Tau secondary in drift
              tratio=one
              if(ifr.gt.1.and.itausec.eq.1) tratio=drifdifr

c Use more accurate formula for the drift
              hafzr2=(half*znuc(iwctype(iwnuc))*rminn)**2
              v2old=0
              do 255 k=1,ndim
  255           v2old=v2old+voldw(k,i,iw,ifr)**2
              volda=sqrt(v2old)
              adrift=(half*(one+eps+voldr/volda))
     &        +adrift0*hafzr2/(one+hafzr2)

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
          endif


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
c A better solution would be to employ a better way of treating nonlocal psps. in DMC.
c At minimum, if one uses this cutoff, instead of having the factor 5 below, one should replace it by a few times sigma.
c Otherwise, this gives a DMC energy that is too high even in the tau->0 limit, actually especially in this limit,
c if sigma is large.
c         if(dwt.gt.1+5*tau) then
c           if(ipr_sav.eq.0) then
c             ipr_sav=1
c             write(6,'(''Warning: dwt>1+5*tau, nwalk,dwt,ewto,ewtn,fratio(iw,ifr),fration='',i5,9d12.4)')
c    &        nwalk,dwt,ewto,ewtn,fratio(iw,ifr),fration
c           endif
c           dwt=1+5*tau
c         endif

c Exercise population control if dmc or vmc with weights
          if(idmc.gt.0.or.iacc_rej.eq.0) dwt=dwt*ffi

c Set weights and product of weights over last nwprod steps
          if(ifr.eq.1) then

c Temporarily hard-wire the damping of wts to .9
c           rlambda_tau=0.9d0
            rlambda_tau=1.d0
c           wt(iw)=wt(iw)*dwt
            wt(iw)=(wt(iw)**rlambda_tau)*dwt
            do 266 iparm=1,nparm
  266         wi_w(iparm,iw)=rlambda_tau*wi_w(iparm,iw)+dexponent(iparm)
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
c             r2sum=r2sum+wtg*(unacp(i)*r2o+(one-unacp(i)*r2n)
c 270         risum=risum+wtg*(unacp(i)/dsqrt(r2o)+(one-unacp(i)/dsqrt(r2n))

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
           else
            ro=ajacold(iw,ifr)*psidow(iw,ifr)**2*
     &         exp(2*psijow(iw,ifr)-psi2savo)

            wsum1(ifr)=wsum1(ifr)+wtnow*ro
            esum1(ifr)=esum1(ifr)+wtnow*eoldw(iw,ifr)*ro
            pesum(ifr)=pesum(ifr)+wtg*peow(iw,ifr)*ro
            peisum(ifr)=peisum(ifr)+wtg*peiow(iw,ifr)*ro
            tpbsum(ifr)=tpbsum(ifr)+wtg*(eoldw(iw,ifr)-peow(iw,ifr))*ro
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*d2ow(iw,ifr)*ro
          endif

  280   continue
c Call to grad_hess_jas_sum() used to be for optimizing Jastrow for periodic systems.
        call grad_hess_jas_sum(1.d0,0.d0,eoldw(iw,1),eoldw(iw,1),wt(iw)*fprod,wi_w(1,iw))
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
      if(ipr.gt.-2) write(11,'(i8,f9.6,f12.5,f11.6,i5)') ipass,ffn,
     &wsum1(1),esum1(1)/wsum1(1),nwalk


      return
      end
