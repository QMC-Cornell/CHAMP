      subroutine splitj_movall
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)
      include '../vmc/vmc.h'
      include 'dmc.h'
      include '../vmc/force.h'

      parameter (zero=0.d0,two=2.d0,half=.5d0)

      common /dim/ ndim
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /force_dmc/ itausec,nwprod
      common /config_dmc/ xoldw(3,MELEC,MWALK,MFORCE),voldw(3,MELEC,MWALK,MFORCE),
     &psidow(MWALK,MFORCE),psijow(MWALK,MFORCE),peow(MWALK,MFORCE),peiow(MWALK,MFORCE),d2ow(MWALK,MFORCE)
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /stats/ dfus2ac,dfus2un(MFORCE),dr2ac,dr2un,acc,acc_int,try_int,
     &nbrnch,nodecr
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eoldw(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /div_v/ div_vo(MELEC,MWALK)

      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)

      dimension iwundr(MWALK)

      iunder=0
      ipair=0
      wtsm=zero
      do 10 iw=1,nwalk
        wtsm=wtsm+wt(iw)
        if(wt(iw).lt.half) then
          if(wt(iw).eq.zero) then
            nbrnch=nbrnch+1
            iunder=iunder+1
            iwundr(iunder)=iw
           else
            if(ipair.eq.0) then
              ipair=1
              iw2=iw
             else
              nbrnch=nbrnch+1
              ipair=0
              iunder=iunder+1
              wttot=wt(iw)+wt(iw2)
              if(rannyu(0).gt.(wt(iw)/wttot)) then
                wt(iw2)=wttot
                iwundr(iunder)=iw
               else
                wt(iw)=wttot
                iwundr(iunder)=iw2
              endif
            endif
          endif
        endif
   10 continue

      nwalk2=nwalk
      do 20 iw=1,nwalk
        if(wt(iw).ge.two) then
          nbrnch=nbrnch+1
          if(iunder.gt.0) then
            iw2=iwundr(iunder)
            iunder=iunder-1
           else
            nwalk2=nwalk2+1
            iw2=nwalk2
            if(nwalk2.gt.MWALK) stop 'MWALK exceeded in splitj'
          endif
          wt(iw)=wt(iw)*half
          wt(iw2)=wt(iw)
          iage(iw2)=iage(iw)
          do 15 ifr=1,nforce
            ajacold(iw2,ifr)=ajacold(iw,ifr)
            eoldw(iw2,ifr)=eoldw(iw,ifr)
            psidow(iw2,ifr)=psidow(iw,ifr)
            psijow(iw2,ifr)=psijow(iw,ifr)
            peow(iw2,ifr)=peow(iw,ifr)
            peiow(iw2,ifr)=peiow(iw,ifr)
            d2ow(iw2,ifr)=d2ow(iw,ifr)
            pwt(iw2,ifr)=pwt(iw,ifr)
            do 12 ip=0,nwprod-1
   12         wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
            do 15 i=1,nelec
              div_vo(i,iw2)=div_vo(i,iw)
              do 15 k=1,ndim
                voldw(k,i,iw2,ifr)=voldw(k,i,iw,ifr)
   15           xoldw(k,i,iw2,ifr)=xoldw(k,i,iw,ifr)
        endif
   20 continue

      do 30 j=iunder,1,-1
        iw2=iwundr(j)
        iw=nwalk2
        nwalk2=nwalk2-1
        wt(iw2)=wt(iw)
        iage(iw2)=iage(iw)
        do 30 ifr=1,nforce
          ajacold(iw2,ifr)=ajacold(iw,ifr)
          eoldw(iw2,ifr)=eoldw(iw,ifr)
          psidow(iw2,ifr)=psidow(iw,ifr)
          psijow(iw2,ifr)=psijow(iw,ifr)
          peow(iw2,ifr)=peow(iw,ifr)
          peiow(iw2,ifr)=peiow(iw,ifr)
          d2ow(iw2,ifr)=d2ow(iw,ifr)
          pwt(iw2,ifr)=pwt(iw,ifr)
          do 25 ip=0,nwprod-1
   25       wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
          do 30 i=1,nelec
            div_vo(i,iw2)=div_vo(i,iw)
            do 30 k=1,ndim
              voldw(k,i,iw2,ifr)=voldw(k,i,iw,ifr)
   30         xoldw(k,i,iw2,ifr)=xoldw(k,i,iw,ifr)
      nwalk=nwalk2

      wtsm2=zero
      do 40 iw=1,nwalk
        wtsm2=wtsm2+wt(iw)
c       if(wt(iw).lt.half) write(11,'(i4,9d12.5)') iw,wt(iw),eoldw(iw)
c       if(wt(iw).gt.two) write(11,'(i4,9d12.5)') iw,wt(iw),eoldw(iw)
   40 continue
c     if(dabs(wtsm-wtsm2).gt.1.d-10) write(11,'(2f12.6)') wtsm,wtsm2

      return
      end
