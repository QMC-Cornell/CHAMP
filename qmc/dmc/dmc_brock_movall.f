      subroutine dmc_brock_movall
c Written by Cyrus Umrigar
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  Does either VMC or simple DMC
c  if  idmc < 0  VMC
c  if  idmc >=0  DMC  (Note until Jan92 idmc=0 did VMC)
c           = 0  DMC without accept/reject, kill   , not limit E,V
c           = 1  DMC without accept/reject, kill   , limit E,V
c           = 2  DMC with    accept/reject, kill   , not limit E,V
c           = 3  DMC with    accept/reject, reject , not limit E,V
c           = 4  DMC with    accept/reject, reject , limit E,V
c           = 5  DMC with    accept/reject, cross  , limit E,V
c  If idmc =2 walkers that cross nodes while doing DMC are killed.
c  Hence detailed balance is violated.  This is done only for comparing
c  with Reynolds et. al. J. Chem. Phys. 1982.
c  However, this is not quite the Reynolds algorithm since they do the
c  accept/reject on each electron and calculate taueff differently.
c  He and Ceperley told me that they have also now switched to rejecting
c  node crossing moves.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      use constants_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use forcepar_mod
      use delocc_mod
      use force_dmc_mod
      use iterat_mod
      use jacobsave_mod
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

      dimension xnew(3,nelec,nforce),vnew(3,nelec,nforce),enew(nforce)
     &,psidn(nforce),psijn(nforce),d2n(nforce),pen(nforce),pein(nforce),xbac(3)
      dimension ajacnew(nforce),div_vn(nelec)


      gauss()=dcos(two*pi*rannyu(0))*dsqrt(-two*dlog(rannyu(0)))

      stop 'dmc_brock not yet updated for new /cntrld/ variables'

      ajacnew(1)=one

      ipmod=mod(ipass,nfprod)
      ipmod2=mod(ipass+1,nfprod)
      ffn=eigv*wdsumo/nconf_global
      ffi=one/ffn
      fprod=fprod*ffn/ff(ipmod)
      ff(ipmod)=ffn

      fprodd=fprod/ff(ipmod2)

c Undo weights
      iwmod=mod(ipass,nwprod)

      idmcabs=iabs(idmc)
      ioldest=0
      do 200 iw=1,nwalk
c Loop over primary+secondary paths to compute forces
        do 200 ifr=1,nforce
c Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
          if(nforce.gt.1) call strech(xnew(1,1,1),xnew(1,1,ifr),ajacob,ifr,0)
          dr2=zero
          dfus2o=zero
          do 60 i=1,nelec
            v2old=0
            do 5 k=1,ndim
    5         v2old=v2old+voldw(k,i,iw,ifr)**2
c Note that it is better to limit the drift term for VMC also,
c so idmc=-1,-4,-5 are preferable
            if(idmcabs.eq.1 .or. idmcabs.eq.4. or. idmcabs.eq.5 ) then
              term=min(one/(dsqrt(v2old)*tau),one)
             else
              term=one
            endif
            if(ifr.eq.1) then
              do 10 k=1,ndim
                drift=term*voldw(k,i,iw,ifr)*tau
                dfus=gauss()*rttau
                dx=drift+dfus
                dr2=dr2+dx**2
                dfus2o=dfus2o+dfus**2
   10           xnew(k,i,ifr)=xoldw(k,i,iw,ifr)+dx
             else
c Tau secondary in drift in order to compute tau_s in second equil. block
              tratio=one
              if(itausec.eq.1) tratio=taueff(ifr)/tau

              do 15 k=1,ndim
                drift=term*voldw(k,i,iw,ifr)*tau*tratio
                dx=xnew(k,i,ifr)-xoldw(k,i,iw,ifr)
                dfus=dx-drift
                dr2=dr2+dx**2
   15           dfus2o=dfus2o+dfus**2
             endif

            if(ipr.ge.1)
     &      write(6,'(''xnewdr'',2i4,9f8.5)') iw,i,(xnew(k,i,ifr),k=1,ndim)

   60     continue

c calculate psi etc. at new configuration
          call hpsi(xnew(1,1,ifr),psidn(ifr),psijn(ifr),vnew(1,1,ifr),div_vn,d2n(ifr),pen(ifr),pein(ifr),enew(ifr),denergy,ifr)

c Truncate energy as in Brock if iabs(idmc) = 1 or 4 or 5
          if((idmcabs.eq.1 .or. idmcabs.eq.4 .or. idmcabs.eq.5 ) .and.
     &    dabs(enew(ifr)-etrial).gt.two/rttau)
     &    enew(ifr)=etrial+sign(two/rttau,enew(ifr)-etrial)

c Check if walker has crossed a node.
c Note that the better thing to do is to reject the move (p=zero)
c rather than to kill the walker (wt(iw)=zero).  The latter is only
c done for the sake of consistency with Reynolds et. al. 1982.
          if((psidn(ifr)*psidow(iw,ifr)).le.zero .and.ifr.eq.1) then
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
            do 160 i=1,nelec
c Note that it is better to limit the drift term for VMC also,
c so idmc=-1,-4 are preferable
              if(idmcabs.eq.1 .or. idmcabs.eq.4 .or. idmcabs.eq.5 ) then
              v2new=0
                do 151 k=1,ndim
  151             v2new=v2new+vnew(k,i,ifr)**2
                term=min(one/(dsqrt(v2new)*tau),one)
               else
                term=one
              endif
              do 153 k=1,ndim
                drift=term*vnew(k,i,ifr)*tau
                xbac(k)=xnew(k,i,ifr)+drift
                dfus=xbac(k)-xoldw(k,i,iw,ifr)
  153           dfus2n=dfus2n+dfus**2
  160         continue

              p=(psidn(ifr)/psidow(iw,ifr))**2*exp(2*(psijn(ifr)-
     &        psijow(iw,ifr)))*exp((dfus2o-dfus2n)/(two*tau))

              p=dmin1(one,p)
          endif

          if(ipr.ge.1) write(6,'(9d13.5)') p,psidow(iw,ifr),psidn(ifr),
     &    psijow(iw,ifr),psijn(ifr),dfus2o,dfus2n,tau

c 165     q=one-p
  165     continue

          if(ifr.eq.1) then
            do 36 ifs=2,nforce
              if(istrech.eq.0) then
c No streched positions for electrons
                do 35 i=1,nelec
                  do 35 k=1,ndim
   35               xnew(k,i,ifs)=xnew(k,i,1)
                ajacnew(ifs)=one
               else
c Compute streched electronic positions for all nucleus displacement
                call strech(xnew(1,1,1),xnew(1,1,ifs),ajacob,ifs,1)
                ajacnew(ifs)=ajacob
              endif
   36       continue
          endif

          dfus2unf(ifr)=dfus2unf(ifr)+dfus2o
          if(ifr.eq.1) then
            acc=acc+p
            try_int=try_int+1
            dfus2ac=dfus2ac+p*dfus2o
            dr2ac=dr2ac+p*dr2
            dr2un=dr2un+dr2
            tautot=tautot+taueff(ifr)

            if(rannyu(0).lt.p) then
              iaccept=1
              iage(iw)=0
             else
              iaccept=0
              iage(iw)=iage(iw)+1
              ioldest=max(ioldest,iage(iw))
              ioldestmx=max(ioldestmx,iage(iw))
            endif

          endif

          wtmult=one
          if(iaccept.eq.1) then
            if(idmc.ge.0) then
              wtmult=dexp((etrial-half*(eoldw(iw,ifr)+enew(ifr)))
     &        *taueff(1))
              if(ifr.gt.1.and.itausec.eq.1)
     &        wtmult=dexp((etrial-half*(eoldw(iw,ifr)+enew(ifr)))
     &        *taueff(ifr))
            endif

            psidow(iw,ifr)=psidn(ifr)
            psijow(iw,ifr)=psijn(ifr)
            eoldw(iw,ifr)=enew(ifr)
            peow(iw,ifr)=pen(ifr)
            peiow(iw,ifr)=pein(ifr)
            d2ow(iw,ifr)=d2n(ifr)
            do 190 i=1,nelec
              do 190 k=1,ndim
                xoldw(k,i,iw,ifr)=xnew(k,i,ifr)
  190           voldw(k,i,iw,ifr)=vnew(k,i,ifr)

            ajacold(iw,ifr)=ajacnew(ifr)
           else
            if(idmc.ge.0) then
              wtmult=dexp((etrial-eoldw(iw,ifr))*taueff(1))
              if(ifr.gt.1.and.itausec.eq.1)
     &        wtmult=dexp((etrial-eoldw(iw,ifr))*taueff(ifr))
            endif
          endif

          dwt=wtmult*ffi

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

          ro=one
          if(ifr.gt.1) ro=(psidow(iw,ifr)/psidow(iw,1))**2*
     &    exp(2*(psijow(iw,ifr)-psijow(iw,1)))*ajacold(iw,ifr)

          wtg=wtnow*fprod

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

  200 continue

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

      do 205 ifr=1,nforce
        wgsum1(ifr)=wsum1(ifr)*fprod
        egsum1(ifr)=esum1(ifr)*fprod
        wgsum(ifr)=wgsum(ifr)+wgsum1(ifr)
 205    egsum(ifr)=egsum(ifr)+egsum1(ifr)

      wsum=wsum+wsum1(1)
      wfsum=wfsum+wfsum1
      wdsum=wdsum+wdsumo
      wgdsum=wgdsum+wgdsumo
      esum=esum+esum1(1)
      efsum=efsum+efsum1
      eisum=eisum+wfsum1/wdsumo

!JT      call splitj ! moved outside the routine
      if(ipr.gt.-2) write(11,'(i8,f9.6,f12.5,f11.6,i5)') ipass,ffn,
     &wsum1(1),esum1(1)/wsum1(1),nwalk

      nfpro=min(nfprod,ipass)
      eigv=(wgsum1(1)/wtgen(ipmod))**(one/nfpro)

      wdsumo=wdsumn
      wgdsumo=wgdsumn
      wtgen(ipmod)=wdsumn

      return
      end
