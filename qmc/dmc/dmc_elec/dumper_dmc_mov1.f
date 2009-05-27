      subroutine dumper_dmc_mov1
c Written by Cyrus Umrigar, modified by Claudia Filippi
c routine to pick up and dump everything needed to restart
c job where it left off
      use constants_mod
      use control_mod
      use atom_mod
      use coefs_mod
      use dets_mod
      use basis1_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use numbas_mod
      use forcepar_mod
      use pseudo_mod
      use contrl_per_mod
      use delocc_mod
      use jaspar_mod
      use force_dmc_mod
      use iterat_mod
      use qua_mod
      use jacobsave_mod
      use forcest_dmc_mod
      use denupdn_mod
      use stepv_mod
      use jaspar1_mod
      use config_dmc_mod
      use branch_mod
      use estsum_dmc_mod
      use estcum_dmc_mod
      use div_v_dmc_mod
      use contrldmc_mod
      use estcm2_mod
      use stats_mod
      use age_mod
      implicit real*8(a-h,o-z)

      parameter (small=1.d-6)

      common /velratio/ fratio(MWALK,MFORCE)

      dimension irn(4)
      dimension coefx(nbasis,norb),zexx(nbasis),centx(3,ncent),znucx(ncent)
     &,n1sx(MCTYPE),n2sx(MCTYPE),n2px(-1:1,MCTYPE)
     &,n3sx(MCTYPE),n3px(-1:1,MCTYPE),n3dx(-2:2,MCTYPE)
c    &,n4sx(MCTYPE),n4px(-1:1,MCTYPE),n4dx(-2:2,MCTYPE)
     &,n4sx(MCTYPE),n4px(-1:1,MCTYPE)
     &,n4fx(-3:3,MCTYPE),n5gx(-4:4,MCTYPE),n6hx(-5:5,MCTYPE)
     &,nsax(MCTYPE),npax(-1:1,MCTYPE),ndax(-2:2,MCTYPE)
     &,csf_coefx(MDET)

      if(nforce.gt.1) call strech(xoldw,xoldw,ajacob,1,0)

      call savern(irn)
      rewind 10
      write(10) irn
      write(10) hb
      write(10) tau,rttau,taueff(1),tautot,idmc,nfprod
      write(10) nelec,nconf_global,nforce
      write(10) nwalk
      write(10) (wtgen(i),ff(i),i=0,nfprod),(wt(i),i=1,nwalk)
     &,eigv,eest,wdsumo,wgdsumo,fprod
      if(nforce.gt.1) write(10) (taueff(i),i=2,nforce)
      if(nforce.gt.1) write(10) nwprod
     &,((pwt(i,j),i=1,nwalk),j=1,nforce)
     &,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
      write(10) (iage(i),i=1,nwalk),ioldest,ioldestmx
      write(10) (((xoldw(k,i,iw,1),k=1,ndim),i=1,nelec),iw=1,nwalk)
      write(10) ((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
      write(10) wcum,wfcum,(wgcum(i),i=1,nforce),wdcum,wgdcum,wcum1
     &,wfcum1,(wgcum1(i),i=1,nforce),wdcum1, ecum,efcum
     &,(egcum(i),i=1,nforce), ecum1,efcum1,(egcum1(i),i=1,nforce)
     &,ei1cum,ei2cum,ei3cum, (pecum(i),i=1,nforce)
     &,(tpbcum(i),i=1,nforce),(tjfcum(i),i=1,nforce),r2cum,ricum
     &,(taucum(i),i=1,nforce)
      write(10) ipass,iblk
      write(10) wcm2,wfcm2,(wgcm2(i),i=1,nforce),wdcm2,wgdcm2,wcm21
     &,wfcm21,(wgcm21(i),i=1,nforce),wdcm21, ecm2,efcm2
     &,(egcm2(i),i=1,nforce), ecm21,efcm21,(egcm21(i),i=1,nforce)
     &,ei1cm2,ei2cm2,ei3cm2, (pecm2(i),i=1,nforce)
     &,(tpbcm2(i),i=1,nforce),(tjfcm2(i),i=1,nforce),r2cm2,ricm2
      write(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
      write(10) (rprob(i),rprobup(i),rprobdn(i),i=1,NRAD)
      write(10) dfus2ac,dfus2un,dr2ac,dr2un,acc
     &,acc_int,try_int,nbrnch,nodecr
      if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)

      write(10) ((coef(ib,i,1),ib=1,nbasis),i=1,norb)
      write(10) nbasis
      write(10) (zex(ib,1),ib=1,nbasis)
      write(10) nctype,ncent,(iwctype(i),i=1,ncent)
      write(10) ((cent(k,i),k=1,ndim),i=1,ncent)
      write(10) pecent
      write(10) (znuc(i),i=1,nctype)
      write(10) (n1s(i),i=1,nctype)
      if(numr.le.0) write(10) (n2s(i),i=1,nctype)
      write(10) ((n2p(m,i),m=-1,1),i=1,nctype)
      if(numr.le.0) then
        write(10) (n3s(i),i=1,nctype)
        write(10) ((n3p(m,i),m=-1,1),i=1,nctype)
      endif
      write(10) ((n3d(m,i),m=-2,2),i=1,nctype)
      if(numr.le.0) then
        write(10) (n4s(i),i=1,nctype)
        write(10) ((n4p(m,i),m=-1,1),i=1,nctype)
c       write(10) ((n4d(m,i),m=-2,2),i=1,nctype)
        write(10) (nsa(i),i=1,nctype)
        write(10) ((npa(m,i),m=-1,1),i=1,nctype)
        write(10) ((nda(m,i),m=-2,2),i=1,nctype)
       else
        write(10) ((n4f(m,i),m=-3,3),i=1,nctype)
        write(10) ((n5g(m,i),m=-4,4),i=1,nctype)
        write(10) ((n6h(m,i),m=-5,5),i=1,nctype)
      endif
      write(10) (csf_coef(i,1),i=1,ncsf)
      write(10) ncsf,ndet,nup,ndn
      write(10) cjas1(1),cjas2(1)
      rewind 10
      write(6,'(1x,''successful dump to restart_dmc'')')
      return

      entry startr_dmc_mov1
      write(6,'(1x,''attempting restart from restart_dmc'')')
      rewind 10
      read(10) irn
      call setrn(irn)
      read(10) hbx
      read(10) taux,rttau,taueff(1),tautot,idmc,nfprod
      read(10) nelecx,nconf_global,nforce
      if(dabs(hbx-hb).gt.small) stop 'hb'
      if(dabs(taux-tau).gt.small) stop 'tau'
      if(nelecx.ne.nelec) stop 'nelec'
      read(10) nwalk
      read(10) (wtgen(i),ff(i),i=0,nfprod),(wt(i),i=1,nwalk)
     &,eigv,eest,wdsumo,wgdsumo,fprod
      if(nforce.gt.1) read(10) (taueff(i),i=2,nforce)
      if(nforce.gt.1) read(10) nwprod
     &,((pwt(i,j),i=1,nwalk),j=1,nforce)
     &,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
      read(10) (iage(i),i=1,nwalk),ioldest,ioldestmx
      read(10) (((xoldw(k,i,iw,1),k=1,ndim),i=1,nelec),iw=1,nwalk)
      read(10) ((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
      read(10) wcum,wfcum,(wgcum(i),i=1,nforce),wdcum,wgdcum,wcum1
     &,wfcum1,(wgcum1(i),i=1,nforce),wdcum1, ecum,efcum
     &,(egcum(i),i=1,nforce), ecum1,efcum1,(egcum1(i),i=1,nforce)
     &,ei1cum,ei2cum,ei3cum, (pecum(i),i=1,nforce)
     &,(tpbcum(i),i=1,nforce),(tjfcum(i),i=1,nforce),r2cum,ricum
     &,(taucum(i),i=1,nforce)
      read(10) ipass,iblk
      read(10) wcm2,wfcm2,(wgcm2(i),i=1,nforce),wdcm2,wgdcm2,wcm21
     &,wfcm21,(wgcm21(i),i=1,nforce),wdcm21, ecm2,efcm2
     &,(egcm2(i),i=1,nforce), ecm21,efcm21,(egcm21(i),i=1,nforce)
     &,ei1cm2,ei2cm2,ei3cm2, (pecm2(i),i=1,nforce)
     &,(tpbcm2(i),i=1,nforce),(tjfcm2(i),i=1,nforce),r2cm2,ricm2
      read(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
      read(10) (rprob(i),rprobup(i),rprobdn(i),i=1,NRAD)
      read(10) dfus2ac,dfus2un,dr2ac,dr2un,acc
     &,acc_int,try_int,nbrnch,nodecr
      if(nloc.gt.0) then
        read(10) nqx,(xq(i),yq(i),zq(i),wq(i),i=1,nqx)
        if(nqx.ne.nquad) stop 'nquad'
      endif

      read(10) ((coefx(ib,i),ib=1,nbasis),i=1,norb)
      read(10) nbasx
      do 10 j=1,norb
      do 10 i=1,nbasis
      if(dabs(coefx(i,j)-coef(i,j,1)).gt.small) stop 'coef'
   10 continue
      if(nbasx.ne.nbasis) stop 'nbasis'
      read(10) (zexx(ib),ib=1,nbasis)
      read(10) nctypex,ncentx,(iwctype(i),i=1,ncentx)
      read(10) ((centx(k,i),k=1,ndim),i=1,ncentx)
      read(10) pecent
      read(10) (znucx(i),i=1,nctypex)
      read(10) (n1sx(i),i=1,nctype)
      if(numr.le.0) read(10) (n2sx(i),i=1,nctype)
      read(10) ((n2px(m,i),m=-1,1),i=1,nctype)
      if(numr.le.0) then
        read(10) (n3sx(i),i=1,nctype)
        read(10) ((n3px(m,i),m=-1,1),i=1,nctype)
      endif
      read(10) ((n3dx(m,i),m=-2,2),i=1,nctype)
      if(numr.le.0) then
        read(10) (n4sx(i),i=1,nctype)
        read(10) ((n4px(m,i),m=-1,1),i=1,nctype)
c       read(10) ((n4dx(m,i),m=-2,2),i=1,nctype)
        read(10) (nsax(i),i=1,nctype)
        read(10) ((npax(m,i),m=-1,1),i=1,nctype)
        read(10) ((ndax(m,i),m=-2,2),i=1,nctype)
       else
        read(10) ((n4fx(m,i),m=-3,3),i=1,nctype)
        read(10) ((n5gx(m,i),m=-4,4),i=1,nctype)
        read(10) ((n6hx(m,i),m=-5,5),i=1,nctype)
      endif

      if(ncentx.ne.ncent) stop 'ncent'
      if(nctypex.ne.nctype) stop 'nctype'
      do 20 i=1,nbasis
      if(dabs(zexx(i)-zex(i,1)).gt.small) stop 'zex'
   20 continue
      do 30 i=1,ncent
      do 30 k=1,ndim
      if(dabs(cent(k,i)-centx(k,i)).gt.small) stop 'cent'
   30 continue
      do 49 i=1,nctype
        if(n1s(i).ne.n1sx(i)) stop 'n1s'
        if(numr.le.0) then
          if(n2s(i).ne.n2sx(i)) stop 'n2s'
          if(n3s(i).ne.n3sx(i)) stop 'n3s'
          if(n4s(i).ne.n4sx(i)) stop 'n4s'
          if(nsa(i).ne.nsax(i)) stop 'nsa'
        endif
        do 40 m=-1,1
          if(n2p(m,i).ne.n2px(m,i)) stop 'n2p'
          if(numr.le.0) then
            if(n3p(m,i).ne.n3px(m,i)) stop 'n3p'
            if(n4p(m,i).ne.n4px(m,i)) stop 'n4p'
            if(npa(m,i).ne.npax(m,i)) stop 'npa'
          endif
   40   continue
        do 45 m=-2,2
          if(n3d(m,i).ne.n3dx(m,i)) stop 'n3d'
          if(numr.le.0) then
c           if(n4d(m,i).ne.n4dx(m,i)) stop 'n4d'
            if(nda(m,i).ne.ndax(m,i)) stop 'nda'
          endif
   45   continue
        if(numr.gt.0) then
          do 46 m=-3,3
   46       if(n4f(m,i).ne.n4fx(m,i)) stop 'n4f'
          do 47 m=-4,4
   47       if(n5g(m,i).ne.n5gx(m,i)) stop 'n5g'
          do 48 m=-5,5
   48       if(n6h(m,i).ne.n6hx(m,i)) stop 'n6h'
        endif
   49 continue
      read(10) (csf_coefx(i),i=1,ncsf)
      read(10) ncsfx,ndetx,nupx,ndnx
      do 50 i=1,ncsf
      if(dabs(csf_coefx(i)-csf_coef(i,1)).gt.small) stop 'csf_coef'
   50 continue
      if(ncsfx.ne.ncsf) stop 'ncsf'
      if(ndetx.ne.ndet) stop 'ndet'
      if(nupx.ne.nup) stop 'nup'
      if(ndnx.ne.ndn) stop 'ndn'
      read(10) cjas1x,cjas2x
      if(dabs(cjas1x-cjas1(1)).gt.small) stop 'cjas1'
      if(dabs(cjas2x-cjas2(1)).gt.small) stop 'cjas2'
      write(6,'(1x,''succesful read from restart_dmc'')')
      write(6,'(t5,''egnow'',t15,''egave'',t21
     &,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55
     &,''(tpberr)'' ,t66,''tjfave'',t72,''(tjferr)'',t83,''npass'',t93
     &,''wgsum'',t103 ,''ioldest'')')

      do 70 iw=1,nwalk
        if(istrech.eq.0) then
          do 60 ifr=2,nforce
            do 60 ie=1,nelec
              do 60 k=1,ndim
   60           xoldw(k,ie,iw,ifr)=xoldw(k,ie,iw,1)
        endif
        do 70 ifr=1,nforce
          if(nforce.gt.1) then
            if(ifr.eq.1.or.istrech.eq.0) then
              call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,ifr),ajacob,ifr,0)
               else
              call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,ifr),ajacob,ifr,1)
            endif
           else
            ajacob=one
          endif
          ajacold(iw,ifr)=ajacob
          call hpsi(xoldw(1,1,iw,ifr),psidow(iw,ifr),psijow(iw,ifr),voldw(1,1,iw,ifr),div_vow(1,iw),d2ow(iw,ifr),
     &    peow(iw,ifr),peiow(iw,ifr),eoldw(iw,ifr),denergy,ifr)
          if(ifr.eq.1) then
            if(ibasis.eq.3) then                ! complex calculation
              call cwalksav_det(iw)
            else
              call walksav_det(iw)
            endif
            call walksav_jas(iw)
          endif
   70 continue

c zero out xsum variables
      wsum=0
      wfsum=0
      wdsum=0
      wgdsum=0
      esum=0
      efsum=0
      ei1sum=0
      ei2sum=0
      r2sum=0
      risum=0

      do 75 ifr=1,nforce
        egsum(ifr)=0
        wgsum(ifr)=0
        pesum(ifr)=0
        peisum(ifr)=0
        tpbsum(ifr)=0
   75   tjfsum(ifr)=0

      if(ipr.gt.-2) then
        open(11,file='walkalize',status='old')
        do 80 i=1,2000000000
   80     read(11,fmt=*,end=90)
      endif
   90 backspace 11
      backspace 11

      return
      end
