      subroutine dmc_good_ps_movall
c Written by Cyrus Umrigar and Claudia Filippi starting from dmc_good
c Uses the diffusion Monte Carlo algorithm described in:
c 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
c    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993).
c with portions related to nuclear cusps removed.
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
c             >= 1 *   use smooth formulae to limit branching to (1/2,2)
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
      implicit real*8(a-h,o-z)
      include '../vmc/vmc.h'
      include 'dmc.h'
      include '../vmc/force.h'
      include '../fit/fit.h'

      common /dim/ ndim
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /force_dmc/ itausec,nwprod
      parameter (zero=0.d0,one=1.d0,two=2.d0,half=.5d0)
      parameter (adrift=0.5d0)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /config_dmc/ xoldw(3,MELEC,MWALK,MFORCE),voldw(3,MELEC,MWALK,MFORCE),
     &psidow(MWALK,MFORCE),psijow(MWALK,MFORCE),peow(MWALK,MFORCE),peiow(MWALK,MFORCE),d2ow(MWALK,MFORCE)
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /stats/ dfus2ac,dfus2un(MFORCE),dr2ac,dr2un,acc,acc_int,try_int,
     &nbrnch,nodecr
      common /delocc/ denergy(MPARM)
      common /estsum_dmc/ wsum,w_acc_sum,wfsum,wgsum(MFORCE),wg_acc_sum,wdsum,
     &wgdsum, wsum1(MFORCE),w_acc_sum1,wfsum1,wgsum1(MFORCE),wg_acc_sum1,
     &wdsum1, esum,efsum,egsum(MFORCE),esum1(MFORCE),efsum1,egsum1(MFORCE),
     &ei1sum,ei2sum,ei3sum, pesum(MFORCE),peisum(MFORCE),tpbsum(MFORCE),tjfsum(MFORCE),r2sum,
     &risum,tausum(MFORCE)
      common /estcum_dmc/ wcum,w_acc_cum,wfcum,wgcum(MFORCE),wg_acc_cum,wdcum,
     &wgdcum, wcum1,w_acc_cum1,wfcum1,wgcum1(MFORCE),wg_acc_cum1,
     &wdcum1, ecum,efcum,egcum(MFORCE),ecum1,efcum1,egcum1(MFORCE),
     &ei1cum,ei2cum,ei3cum, pecum(MFORCE),peicum(MFORCE),tpbcum(MFORCE),tjfcum(MFORCE),r2cum,
     &ricum,taucum(MFORCE)
      common /stepv/ try(NRAD),suc(NRAD),trunfb(NRAD),rprob(NRAD),
     &ekin(NRAD),ekin2(NRAD)
      common /denupdn/ rprobup(NRAD),rprobdn(NRAD)
      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_global,nconf_new,isite,idump,irstar
      common /contrldmc/ tau,rttau,taueff(MFORCE),tautot,nfprod,idmc,ipq
     &,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
      common /iterat/ ipass,iblk
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eoldw(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)
      common /div_v_dmc/ div_vow(MELEC,MWALK)
      common /tmp/ eacc,enacc,macc,mnacc
      common /pairden/ xx0probut(0:NAX,-NAX:NAX,-NAX:NAX),xx0probuu(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probud(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdt(0:NAX,-NAX:NAX,-NAX:NAX),
     &xx0probdu(0:NAX,-NAX:NAX,-NAX:NAX),xx0probdd(0:NAX,-NAX:NAX,-NAX:NAX),
     &den2d_t(-NAX:NAX,-NAX:NAX),den2d_d(-NAX:NAX,-NAX:NAX),den2d_u(-NAX:NAX,-NAX:NAX),
     &delxi,xmax,xfix(3),ifixe
      common /fourier/ fourierrk_u(0:NAX,0:NAK1),fourierrk_d(0:NAX,0:NAK1)
     &,fourierrk_t(0:NAX,0:NAK1),fourierkk_u(-NAK2:NAK2,-NAK2:NAK2),fourierkk_d(-NAK2:NAK2,-NAK2:NAK2)
     &,fourierkk_t(-NAK2:NAK2,-NAK2:NAK2),delk1,delk2,fmax1,fmax2,ifourier
      dimension xnc(3,MELEC),xoc(3,MELEC)
      dimension ixo(3),ixn(3)

      dimension xnew(3,MELEC,MFORCE),vnew(3,MELEC,MFORCE),psidn(MFORCE),
     &psijn(MFORCE),enew(MFORCE),pen(MFORCE),pein(MFORCE),d2n(MFORCE),div_vn(MELEC)
      dimension xbac(3),ajacnew(MFORCE)


      gauss()=dcos(two*pi*rannyu(0))*sqrt(-two*dlog(rannyu(0)))

c     term=(sqrt(two*pi*tau))**3/pi

      ajacnew(1)=one

c Undo products
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

c Undo weights
      iwmod=mod(ipass,nwprod)

      ioldest=0
      do 300 iw=1,nwalk
c Loop over primary+secondary paths to compute forces
        do 300 ifr=1,nforce

c Stretch nuclei and electrons if calculating forces
          if(nforce.gt.1) then
            if(ifr.eq.1.or.istrech.eq.0) then
c Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
              call strech(xnew(1,1,1),xnew(1,1,ifr),ajacob,ifr,0)
              ajacnew(ifr)=one
             else
              call strech(xnew(1,1,1),xnew(1,1,ifr),ajacob,ifr,1)
              ajacnew(ifr)=ajacob
            endif
          endif

c Sample Green function for forward move
          dr2=zero
          dfus2o=zero
          vav2sumo=zero
          v2sumo=zero
          fnormo=one
          do 100 i=1,nelec
c Tau secondary in drift in order to compute tau_s in second equil. block
c Note: vavvo=vav2sumo/v2sumo appears in the branching
            if(ifr.gt.1.and.itausec.eq.1) then
              if(itau_eff.ge.1) then
                tauu=tau*taueff(ifr)/taueff(1)
               else
                tauu=taueff(ifr)
              endif
             else
              tauu=tau
            endif

c Use more accurate formula for the drift
            v2old=0
            do 5 k=1,ndim
    5         v2old=v2old+voldw(k,i,iw,ifr)**2
            volda=sqrt(v2old)
            vav=(dsqrt(1+2*adrift*v2old*tauu)-1)/(adrift*volda*tauu)
            vav2sumo=vav2sumo+vav**2
            v2sumo=v2sumo+v2old

c Calculate drifted position
            do 60 k=1,ndim
              xbac(k)=xoldw(k,i,iw,ifr)+voldw(k,i,iw,ifr)*tauu*vav/volda
              if(ifr.gt.1) then
                dfus=xnew(k,i,ifr)-xbac(k)
                dfus2o=dfus2o+dfus**2
              endif
   60       continue

c Do the diffusion.  Actually, this diffusion contains part of the drift too
c if idiv_v.ge.1
c Use div_v_hom to set tau_hom
            div_v_hom=div_vow(i,iw)/3

            if(ifr.eq.1) then
              if(idiv_v.ge.1) then
                if(div_v_hom.lt.0.d0) then
                  tau_hom=3*tau-1/div_v_hom
     &            -((tau-1/div_v_hom)*sqrt(-div_v_hom)
     &            *(1-derfc(1/sqrt(-2*div_v_hom*tau)))
     &            +sqrt(2*tau/pi)*exp(1/(2*div_v_hom*tau)))**2
                 else
                  tau_hom=tau*(1+2*div_v_hom*tau)/(1+div_v_hom*tau)
                endif
               else
                tau_hom=tau
              endif
              rttau_hom=dsqrt(tau_hom)

c Sample gaussian
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

c calculate psi etc. at new configuration
          call hpsi(xnew(1,1,ifr),psidn(ifr),psijn(ifr),vnew(1,1,ifr),div_vn,d2n(ifr),pen(ifr),pein(ifr),enew(ifr),denergy,ifr)

c Check for node crossings
          if(psidn(ifr)*psidow(iw,ifr).le.zero.and.ifr.eq.1) then
            nodecr=nodecr+1
            if(icross.le.0) then
              p=zero
              goto 210
            endif
          endif

c Calculate Green function for the reverse move
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

c Calculate drifted position
            dfus2a=0
            do 160 k=1,ndim
              xbac(k)=xnew(k,i,ifr)+vnew(k,i,ifr)*tauu*vav/vnewa
  160         dfus2a=dfus2a+(xoldw(k,i,iw,ifr)-xbac(k))**2

c Use div_v_hom to set tau_hom
            div_v_hom=div_vn(i)/3

            if(idiv_v.ge.1) then
              if(div_v_hom.lt.0.d0) then
                tau_hom=3*tau-1/div_v_hom
     &          -((tau-1/div_v_hom)*sqrt(-div_v_hom)
     &          *(1-derfc(1/sqrt(-2*div_v_hom*tau)))
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
              write(6,'(''xold'',9f10.6)')(xoldw(k,i,iw,ifr),k=1,ndim),
     &        (xnew(k,i,ifr),k=1,ndim), (xbac(k),k=1,ndim)
              write(6,'(''iel,fnormn,fnormn/fnormo'',i3,f10.6,d12.4)')i,fnormn,fnormn/fnormo
            endif
  200     continue
          vavvn=sqrt(vav2sumn/v2sumn)

          p=(psidn(ifr)/psidow(iw,ifr))**2*exp(2*(psijn(ifr)-
     &    psijow(iw,ifr)))*fnormn/fnormo

          if(ipr.ge.1) write(6,'(''p'',9f10.6)')
     &    p,(psidn(ifr)/psidow(iw,ifr))**2*exp(2*(psijn(ifr)-
     &    psijow(iw,ifr))),exp((dfus2o-dfus2n)/(two*tau)),psidn(ifr),
     &    psidow(iw,ifr),psijn(ifr),psijow(iw,ifr),dfus2o,dfus2n

c The following is one reasonable way to cure persistent configurations
c Not needed if itau_eff <=0 and in practice we have never needed it even
c otherwise
          if(iage(iw).gt.50) p=p*1.1d0**(iage(iw)-50)

          pp=p
          p=dmin1(one,p)
  210     q=one-p

          dfus2un(ifr)=dfus2un(ifr)+dfus2o

c Set taunow for branching.  Note that for primary walk taueff(1) always
c is tau*dfus2ac/dfus2un so if itau_eff.le.0 then we reset taunow to tau
c On the other hand, secondary walk only has dfus2ac/dfus2un if
c itau_eff.ge.1 so there is no need to reset taunow.
          if(ifr.eq.1) then
            acc=acc+p
            try_int=try_int+1
            dfus2ac=dfus2ac+p*dfus2o
            dr2ac=dr2ac+p*dr2
            dr2un=dr2un+dr2
            tautot=tautot+tau*dfus2ac/dfus2un(1)

c If we are using weights rather than accept/reject
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

          ewto=eest-(eest-eoldw(iw,ifr))*vavvo
          ewtn=eest-(eest   -enew(ifr))*vavvn

          if(idmc.gt.0) then
            expon=(etrial-half*((one+qsav)*ewto+psav*ewtn))*taunow
            if(icut_br.le.0) then
              dwt=dexp(expon)
             else
              dwt=0.5d0+1/(1+exp(-4*expon))
c             if(expon.gt.0) then
c               dwt=(1+2*expon)/(1+expon)
c              else
c               dwt=(1-expon)/(1-2*expon)
c             endif
            endif
          endif

c If we are using weights rather than accept/reject
          if(iacc_rej.eq.0) dwt=dwt*pp

c Exercise population control if dmc or vmc with weights
          if(idmc.gt.0.or.iacc_rej.eq.0) dwt=dwt*ffi

c Set weights and product of weights over last nwprod steps
          pwt(iw,ifr)=pwt(iw,ifr)*dwt/wthist(iw,iwmod,ifr)
          wthist(iw,iwmod,ifr)=dwt
          if(ifr.eq.1) then
            wt(iw)=wt(iw)*dwt
            wtnow=wt(iw)
           else
            wtnow=wt(iw)*pwt(iw,ifr)/pwt(iw,1)
          endif

          if(ipr.ge.1)write(6,'(''eold,enew,wt'',9f10.5)')
     &    eoldw(iw,ifr),enew(ifr),wtnow

          wtg=wtnow*fprod

          enowo=eoldw(iw,ifr)
          enown=enew(ifr)

          if(ifr.eq.1) then
            psi2savo=2*(psijow(iw,1)+dlog(dabs(psidow(iw,1))))
            psi2savn=2*(psijn(1)+dlog(dabs(psidn(1))))

            wsum1(ifr)=wsum1(ifr)+wtnow
            esum1(ifr)=esum1(ifr)+wtnow*(q*enowo+p*enown)
            pesum(ifr)=pesum(ifr)+  wtg*(q*peow(iw,ifr) +p*pen(ifr))
            peisum(ifr)=peisum(ifr)+  wtg*(q*peiow(iw,ifr) +p*pein(ifr))
            tpbsum(ifr)=tpbsum(ifr)+wtg*(q*(enowo-peow(iw,ifr))+
     &                                   p*(enown-pen(ifr)))
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*(q*d2ow(iw,ifr)+
     &                                           p*d2n(ifr))

           else

            ro=ajacold(iw,ifr)*psidow(iw,ifr)**2*
     &         exp(2*psijow(iw,ifr)-psi2savo)
            rn=ajacnew(ifr)*psidn(ifr)**2*
     &         exp(2*psijn(ifr)-psi2savn)

            wsum1(ifr)=wsum1(ifr)+wtnow*(qsav*ro+psav*rn)
            esum1(ifr)=esum1(ifr)+wtnow*(qsav*enowo*ro+psav*enown*rn)
            pesum(ifr)=pesum(ifr)+  wtg*(qsav*peow(iw,ifr)*ro+
     &                                   psav*pen(ifr)*rn)
            peisum(ifr)=peisum(ifr)+  wtg*(qsav*peiow(iw,ifr)*ro+
     &                                   psav*pein(ifr)*rn)
            tpbsum(ifr)=tpbsum(ifr)+wtg*(qsav*(enowo-peow(iw,ifr))*ro+
     &                                   psav*(enown-pen(ifr))*rn)
            tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*(qsav*d2ow(iw,ifr)*ro+
     &                                           psav*d2n(ifr)*rn)
          endif

c Collect density only for primary walk
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
              itryo=min(int(delri*rmino)+1,NRAD)
              itryn=min(int(delri*rminn)+1,NRAD)
              if(i.le.nup) then
                rprobup(itryo)=rprobup(itryo)+wtgq
                rprobup(itryn)=rprobup(itryn)+wtgp
               else
                rprobdn(itryo)=rprobdn(itryo)+wtgq
                rprobdn(itryn)=rprobdn(itryn)+wtgp
              endif
              rprob(itryo)=rprob(itryo)+wtgq
              rprob(itryn)=rprob(itryn)+wtgp
              r2sum=r2sum+wtg*(q*r2o+p*r2n)
              risum=risum+wtg*(q/dsqrt(r2o)+p/dsqrt(r2n))

c calculate 2density related functions:
              if(ifixe.eq.-1 .or. ifixe.ne.-3) then
                do 247 idim=1,ndim
c note that ix can be negative or positive. nint is a better choice.
                  ixo(idim)=nint(delxi*xoldw(idim,i,iw,ifr))
  247             ixn(idim)=nint(delxi*xnew(idim,i,ifr))
                if(abs(ixo(1)).le.NAX .and. abs(ixo(2)).le.NAX) then
                  den2d_t(ixo(1),ixo(2))=den2d_t(ixo(1),ixo(2))+wtgq
                  if(i.le.nup) then
                    den2d_u(ixo(1),ixo(2))=den2d_u(ixo(1),ixo(2))+wtgq
                  else
                    den2d_d(ixo(1),ixo(2))=den2d_d(ixo(1),ixo(2))+wtgq
                  endif
                endif
                if(abs(ixn(1)).le.NAX .and. abs(ixn(2)).le.NAX) then
                  den2d_t(ixn(1),ixn(2))=den2d_t(ixn(1),ixn(2))+wtgp
                  if(i.le.nup) then
                    den2d_u(ixn(1),ixn(2))=den2d_u(ixn(1),ixn(2))+wtgp
                  else
                    den2d_d(ixn(1),ixn(2))=den2d_d(ixn(1),ixn(2))+wtgp
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
c             if(ifourier.eq.1) call fourier2d(wtgp,wtgq,xoc,xnc)
            endif

          endif

          if(iaccept.eq.1) then

c           if(dabs((enew(ifr)-etrial)/etrial).gt.0.2d+0) then
            if(dabs((enew(ifr)-etrial)/etrial).gt.5.0d+0) then
               write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass,
     &         enew(ifr)-etrial,psidn(ifr),psijn(ifr),
     &         ((xnew(k,jj,ifr),k=1,ndim),jj=1,nelec)
            endif

            if(wt(iw).gt.3) write(18,'(i6,i4,3f8.2,30f8.4)') ipass,iw,
     &      wt(iw),enew(ifr)-etrial,eoldw(iw,ifr)-etrial,
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

  300 continue

      if(wsum1(1).gt.1.1*nconf_global) write(18,'(i6,9d12.4)') ipass,ffn,fprod,
     &fprodd,wsum1(1),wgdsumo

      wdsumn=wsum1(1)
      wdsum1=wdsumo
c     wgdsum1=wgdsumo
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

c Do split-join
!JT      call splitj ! moved outside the routine
      if(ipr.gt.-2) write(11,'(i8,f9.6,f12.5,f11.6,i5)') ipass,ffn,
     &wsum1(1),esum1(1)/wsum1(1),nwalk

c Estimate eigenvalue of G from the energy
      eest=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
      if(itau_eff.ge.1) then
        eigv=dexp((etrial-eest)*taueff(1))
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest,
     &  egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
       else
        accavn=acc_int/try_int
        if(icut_br.le.0) then
          eigv=accavn*exp((etrial-eest)*tau)+(one-accavn)
         else
          eigv=one-accavn+accavn*(0.5d0+1/(1+exp(-4*(etrial-eest)*tau)))
        endif
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest,accavn,
     &  egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wdsumn
      wgdsumo=wgdsumn
      wtgen(ipmod)=wdsumn

      return
      end
