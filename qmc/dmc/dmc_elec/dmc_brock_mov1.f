      subroutine dmc_brock_mov1
c Written by Cyrus Umrigar
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  Does either VMC or simple DMC
c  if  idmc < 0  VMC
c  if  idmc >=0  DMC  (Note until Jan92 idmc=0 did VMC)
c           = 0  DMC without accept/reject, kill   , not limit E,F
c           = 1  DMC without accept/reject, kill   , limit E,F
c           = 2  DMC with    accept/reject, kill   , not limit E,F
c           = 3  DMC with    accept/reject, reject , not limit E,F
c           = 4  DMC with    accept/reject, reject , limit E,F
c           = 5  DMC with    accept/reject, cross  , limit E,F
c  If idmc =2 walkers that cross nodes while doing DMC are killed.
c  Hence detailed balance is violated.  This is done only for comparing
c  with Reynolds et. al. J. Chem. Phys. 1982.
c  However, this is not quite the Reynolds algorithm since they do the
c  accept/reject on each electron and calculate taueff differently.
c  He and Ceperley told me that they have also now switched to rejecting
c  node crossing moves.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      use all_tools_mod
      use control_mod
      use average_mod
      use atom_mod

      implicit real*8(a-h,o-z)

      common /dim/ ndim
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod
!JT      parameter (zero=0.d0,one=1.d0,two=2.d0,half=.5d0)
      parameter (eps=1.d-10)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /config_dmc/ xoldw(3,MELEC,MWALK,MFORCE),voldw(3,MELEC,MWALK,MFORCE),
     &psidow(MWALK,MFORCE),psijow(MWALK,MFORCE),peow(MWALK,MFORCE),peiow(MWALK,MFORCE),d2ow(MWALK,MFORCE)
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /stats/ dfus2ac,dfus2un,dr2ac,dr2un,acc,acc_int,try_int,
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
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_global,nconf_new,isite,idump,irstar
      common /contrl_per/ iperiodic,ibasis
      common /contrldmc/ tau,rttau,taueff(MFORCE),tautot,nfprod,idmc,ipq
     &,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
      common /iterat/ ipass,iblk
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eoldw(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /branch_dmc_opt/ denergy_old_dmc(MPARM,MWALK),wi_w(MPARM,MWALK)
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere
      common /basis/ zex(MBASIS,MWF),betaq
     &,n1s(MCTYPE)
     &,n2s(MCTYPE),n2p(-1:1,MCTYPE)
     &,n3s(MCTYPE),n3p(-1:1,MCTYPE),n3d(-2:2,MCTYPE)
     &,n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE),n4f(-3:3,MCTYPE)
     &,n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE),n5f(-3:3,MCTYPE)
     &,n5g(-4:4,MCTYPE)
     &,n6d(-2:2,MCTYPE),n6f(-3:3,MCTYPE),n6g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
     &,n7g(-4:4,MCTYPE),n7h(-5:5,MCTYPE),n7i(-6:6,MCTYPE)
     &,n8i(-6:6,MCTYPE),n8j(-7:7,MCTYPE)
     &,n9k(-8:8,MCTYPE)
     &,n10l(-9:9,MCTYPE)
     &,n11m(-10:10,MCTYPE)
     &,n12n(-11:11,MCTYPE)
     &,n13o(-12:12,MCTYPE)
     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)
      common /div_v_dmc/ div_vow(MELEC,MWALK)

      dimension xnew(3),vnew(3,MELEC),xbac(3),xstrech(3,MELEC)
c     dimension dewto(MPARM),dewtn(MPARM),dexponent(MPARM)


c     gauss()=dcos(two*pi*rannyu(0))*dsqrt(-two*dlog(rannyu(0)))

      stop 'dmc_brock not yet updated for new /contrldmc/ variables'

      ipmod=mod(ipass,nfprod)
      ipmod2=mod(ipass+1,nfprod)
      ffn=eigv*wdsumo/nconf_global
      ffi=one/ffn
      fprod=fprod*ffn/ff(ipmod)
      ff(ipmod)=ffn

      fprodd=fprod/ff(ipmod2)

c Undo weights
      iwmod=mod(ipass,nwprod)

      idmcab=iabs(idmc)
      ioldest=0
      do 200 iw=1,nwalk
c primary path

c Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
        if(nforce.gt.1) call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)

        if(ibasis.eq.3) then                  !complex basis set
          call cwalkstrdet(iw)
        else
          call walkstrdet(iw)
        endif
        call walkstrjas(iw)

        dfus2ac=zero
        dfus2un=zero
        drifdif=zero
        iaccept=0
        do 160 i=1,nelec
          v2old=0
          do 5 k=1,ndim
    5       v2old=v2old+voldw(k,i,iw,1)**2
c Note that it is better to limit the drift term for VMC also,
c so idmc=-1,-4,-5 are preferable
          if(idmcab.eq.1 .or. idmcab.eq.4. or. idmcab.eq.5 ) then
            term=min(one/(dsqrt(v2old)*tau),one)
           else
            term=one
          endif

          dr2=zero
          dfus2o=zero
          do 10 k=1,ndim
            drift=term*voldw(k,i,iw,1)*tau
            dfus=gauss()*rttau
            dx=drift+dfus
            dr2=dr2+dx**2
            dfus2o=dfus2o+dfus**2
   10       xnew(k)=xoldw(k,i,iw,1)+dx

          if(ipr.ge.1)
     &    write(6,'(''xnewdr'',2i4,9f8.5)') iw,i,(xnew(k),k=1,ndim)

c calculate psi and velocity at new configurations
          call hpsiedmc(i,iw,xnew,psidn,psijn,vnew)

c Check if walker has crossed a node.
c Note that the better thing to do is to reject the move (p=zero)
c rather than to kill the walker (wt(iw)=zero).  The latter is only
c done for the sake of consistency with Reynolds et. al. 1982.
          if((psidn*psidow(iw,1)).le.zero .and.ifr.eq.1) then
            nodecr=nodecr+1
            if(idmc.ge.0 .and. idmc.le. 2) then
              wt(iw)=zero
             elseif(idmc.eq.3 .or. idmc.eq.4 ) then
              p=zero
              goto 165
            endif
          endif

          if(idmc.eq.0 .or. idmc.eq.1) then
            p=one
           else
c Calculate gaussian for the reverse move
            dfus2n=zero
c Note that it is better to limit the drift term for VMC also,
c so idmc=-1,-4 are preferable
            if(idmcab.eq.1 .or. idmcab.eq.4 .or. idmcab.eq.5 ) then
              v2new=0
                do 151 k=1,ndim
  151             v2new=v2new+vnew(k,i)**2
              term=min(one/(dsqrt(v2new)*tau),one)
             else
              term=one
            endif
            do 153 k=1,ndim
              drift=term*vnew(k,i)*tau
              xbac(k)=xnew(k)+drift
              dfus=xbac(k)-xoldw(k,i,iw,1)
  153         dfus2n=dfus2n+dfus**2

            p=(psidn/psidow(iw,1))**2*exp(2*(psijn-psijow(iw,1)))*
     &      exp((dfus2o-dfus2n)/(two*tau))
            p=dmin1(one,p)
          endif
  165     q=one-p

          if(ipr.ge.1) write(6,'(15d15.5)') p,psidow(iw,1),psidn,
     &    psijow(iw,1),psijn,dfus2o,dfus2n,tau,vnew(1,i),vnew(2,i),
     &    vnew(3,i),(xnew(k),k=1,ndim)

          acc=acc+p
          try_int=try_int+1
          dfus2ac=dfus2ac+p*dfus2o
          dfus2un=dfus2un+dfus2o
          dr2ac=dr2ac+p*dr2
          dr2un=dr2un+dr2

          if(rannyu(0).lt.p) then
            iaccept=1
            acc_int=acc_int+1

            iage(iw)=0
            do 158 k=1,ndim
              drifdif=drifdif+(xoldw(k,i,iw,1)-xnew(k))**2
              xoldw(k,i,iw,1)=xnew(k)
              do 158 l=1,nelec
  158           voldw(k,l,iw,1)=vnew(k,l)
            psidow(iw,1)=psidn
            psijow(iw,1)=psijn
            call jassav(i)
            if(ibasis.eq.3) then                        ! complex calculations
                call cdetsav(i)
             else
                call detsav(i)
            endif

          endif

  160   continue


        tauprim=tau*dfus2ac/dfus2un

        do 190 ifr=1,nforce

          wtmult=one
          if(iaccept.eq.1) then   ! If any of the 1-electron moves were accepted then

            if(ifr.eq.1) then
c Primary configuration
              drifdifr=one
              if(nforce.gt.1)
     &        call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,1),ajacob,1,0)
              call hpsi(xoldw(1,1,iw,1),psidn,psijn,voldw(1,1,iw,1),div_vow(1,iw),d2n,pen,pein,enew,denergy,1)
              if(ibasis.eq.3) then                        ! complex calculations
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
                do 161 i=1,nelec
                  do 161 k=1,ndim
  161               xoldw(k,i,iw,ifr)=xoldw(k,i,iw,1)
                ajacold(iw,ifr)=one
               else
c Compute streched electronic positions for all nucleus displacement
                call strech(xoldw(1,1,iw,1),xstrech,ajacob,ifr,1)
                drifdifs=zero
                do 162 i=1,nelec
                  do 162 k=1,ndim
                    drifdifs=drifdifs+(xstrech(k,i)-xoldw(k,i,iw,ifr))**2
  162               xoldw(k,i,iw,ifr)=xstrech(k,i)
                ajacold(iw,ifr)=ajacob
                drifdifr=drifdifs/drifdif
              endif
              call hpsi(xoldw(1,1,iw,ifr),psidn,psijn,voldw(1,1,iw,ifr),div_vow(1,iw),d2n,pen,pein,enew,denergy,ifr)
            endif
           else
            drifdifr=one
            enew=eoldw(iw,ifr)
          endif
c Truncate energy as in Brock if iabs(idmc) = 1 or 4 or 5
          if((idmcab.eq.1 .or. idmcab.eq.4 .or. idmcab.eq.5 ) .and.
     &    dabs(enew-etrial).gt.two/rttau)
     &    enew=etrial+sign(two/rttau,enew-etrial)

          taunow=tauprim*drifdifr

          if(idmc.ge.0)
     &    wtmult=dexp((etrial-half*(eoldw(iw,ifr)+enew))*taunow)

          dwt=wtmult*ffi
c Set weights and product of weights over last nwprod steps
          if(ifr.eq.1) then

            wt(iw)=wt(iw)*dwt
            wtnow=wt(iw)
            pwt(iw,ifr)=pwt(iw,ifr)*dwt/wthist(iw,iwmod,ifr)
            wthist(iw,iwmod,ifr)=dwt

           elseif(ifr.gt.1) then

            pwt(iw,ifr)=pwt(iw,ifr)*dwt/wthist(iw,iwmod,ifr)
            wthist(iw,iwmod,ifr)=dwt
            wtnow=wt(iw)*pwt(iw,ifr)/pwt(iw,1)

          endif

          wtg=wtnow*fprod

          if(iaccept.eq.1) then   ! If any of the 1-electron moves were accepted then
            psidow(iw,ifr)=psidn
            psijow(iw,ifr)=psijn
            eoldw(iw,ifr)=enew
            peow(iw,ifr)=pen
            peiow(iw,ifr)=pein
            d2ow(iw,ifr)=d2n
           else
            if(ifr.eq.1) then
              iage(iw)=iage(iw)+1
              ioldest=max(ioldest,iage(iw))
              ioldestmx=max(ioldestmx,iage(iw))
            endif
          endif

          ro=one
          if(ifr.gt.1) ro=(psidow(iw,ifr)/psidow(iw,1))**2*
     &    exp(2*(psijow(iw,ifr)-psijow(iw,1)))*ajacold(iw,ifr)
          wsum1(ifr)=wsum1(ifr)+wtnow*ro
          esum1(ifr)=esum1(ifr)+wtnow*eoldw(iw,ifr)*ro
          pesum(ifr)=pesum(ifr)+wtg*peow(iw,ifr)
          peisum(ifr)=peisum(ifr)+wtg*peiow(iw,ifr)
          tpbsum(ifr)=tpbsum(ifr)+wtg*(eoldw(iw,ifr)-peow(iw,ifr))
          tjfsum(ifr)=tjfsum(ifr)-wtg*half*hb*d2ow(iw,ifr)

          if(ifr.eq.1) then
            do 180 i=1,nelec
              r2=zero
              do 170 k=1,ndim
  170           r2=r2+xoldw(k,i,iw,ifr)**2
              r=sqrt(r2)
              itry=min(int(delri*r)+1,NRAD)
              rprob(itry)=rprob(itry)+wtg
              r2sum=r2sum+wtg*r2
  180         risum=risum+wtg/r
          endif
          tausum(ifr)=tausum(ifr)+wtg*taunow

  190   continue

  200 continue

      wfsum1=wsum1(1)*ffn
      efsum1=esum1(1)*ffn
      do 205 ifr=1,nforce
        wgsum1(ifr)=wsum1(ifr)*fprod
 205    egsum1(ifr)=esum1(ifr)*fprod

      call object_modified_by_index (eoldw_index) !JT
      call object_modified_by_index (wt_index) !JT
      call object_modified_by_index (fprod_index) !JT

!JT      call splitj ! moved outside the routine
      if(ipr.gt.-2) write(11,'(i8,f9.6,f12.5,f11.6,i5)') ipass,ffn,
     &wsum1(1),esum1(1)/wsum1(1),nwalk

      nfpro=min(nfprod,ipass)
      eigv=(wgsum1(1)/wtgen(ipmod))**(one/nfpro)

      return
      end
