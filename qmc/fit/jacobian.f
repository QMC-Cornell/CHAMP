      subroutine jacobian(ndata2,nparm,nanalytic,parm,ajac)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use basic_tools_mod
      use fitdet_mod
      use atom_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use basis1_mod
      implicit real*8(a-h,o-z)
      character*16 mode


c epsder1f=sqrt(max(eps_diff,dbl_epsilon))
c     parameter (eps_diff=1.d-15,dbl_epsilon=2.2204460492503131d-16
c    &,epsder1f=3.16227766016838d-8)

c complex common:
c      complex*16 cvd_sav,cvk_sav

      common /pars/ a00,a20,a21,eps_fock,c0000,c1110,c2000,
     &   xm1,xm2,xm12,xms,xma,Zfock
      common /contr2/ ijas,icusp,icusp2,isc,inum_orb,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch
     &,ipos,idcds,idcdu,idcdt,id2cds,id2cdu,id2cdt,idbds,idbdu,idbdt
      common /contr3/ mode
      common /contrl_opt2/ igradhess,iadd_diag_opt
      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata
!JT      common /optim/ lo(MORB),npoint(MORB),
!JT     &iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE),
!JT     &iwjasf(15,MCTYPE),iwbase(MBASIS),iwbasi(MPARM),iworb(MPARM),
!JT     &iwcsf(MCSF),iebase(2,MBASIS),iebasi(2,MPARM),ieorb(2,MPARM),
!JT     &imnbas(MCENT),
!JT     &nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj,
!JT     &nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE),
!JT     &necn,nebase

!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
!MS Declare arrays upto o-orbitals (l=12) for Jellium sphere
!JT      common /basis/ zex(MBASIS,MWF),betaq
!JT     &,n1s(MCTYPE)
!JT     &,n2s(MCTYPE),n2p(-1:1,MCTYPE)
!JT     &,n3s(MCTYPE),n3p(-1:1,MCTYPE),n3d(-2:2,MCTYPE)
!JT     &,n4s(MCTYPE),n4p(-1:1,MCTYPE),n4d(-2:2,MCTYPE),n4f(-3:3,MCTYPE)
!JT     &,n5s(MCTYPE),n5p(-1:1,MCTYPE),n5d(-2:2,MCTYPE),n5f(-3:3,MCTYPE)
!JT     &,n5g(-4:4,MCTYPE)
!JT     &,n6d(-2:2,MCTYPE),n6f(-3:3,MCTYPE),n6g(-4:4,MCTYPE),n6h(-5:5,MCTYPE)
!JT     &,n7g(-4:4,MCTYPE),n7h(-5:5,MCTYPE),n7i(-6:6,MCTYPE)
!JT     &,n8i(-6:6,MCTYPE),n8j(-7:7,MCTYPE)
!JT     &,n9k(-8:8,MCTYPE)
!JT     &,n10l(-9:9,MCTYPE)
!JT     &,n11m(-10:10,MCTYPE)
!JT     &,n12n(-11:11,MCTYPE)
!JT     &,n13o(-12:12,MCTYPE)
!JT     &,nsa(MCTYPE),npa(-1:1,MCTYPE),nda(-2:2,MCTYPE)
      common /jaspar1/ cjas1(MWF),cjas2(MWF)
      common /jaspar2/ a1(MPARMJ,3,MWF),a2(MPARMJ,3,MWF)
      common /jaspar3/ a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(MPARMJ,MCTYPE,MWF)
     &,fck(15,MCTYPE,MWF),scalek(MWF),nord
      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
!JT      common /dets/ csf_coef(MCSF,MWF),cdet_in_csf(MDET_CSF,MCSF),ndet_in_csf(MCSF),iwdet_in_csf(MDET_CSF,MCSF),ncsf,ndet,nup,ndn
      common /bparm/ nspin2b,nocuspb
      common /delocc/ denergy(MPARM)
c      common /fitdet/ cvd_sav(3,MELEC,MDATA),vd_sav(3,MELEC,MDATA),psid_sav(MDATA)
c     &               ,d2d_sav(MDATA),div_vd_sav(MELEC,MDATA),cvk_sav(3,MELEC,MDATA),psik_sav(MDATA)
c     &               ,div_vk_sav(MELEC,MDATA),d2k_sav(MDATA),iconfg,isaved


      common /mpioffset/ ircounts(0:MPROC),idispls(0:MPROC)

c     dimension parm(nparm),ajac(ndata2,nparm),denergy(MPARM)
      dimension parm(nparm),ajac(ndata2,nparm)
     &,velocity(3,MELEC),div_v(MELEC)

      do 10 iparm=1,nparml
   10   coef(iwbasi(iparm),iworb(iparm),1)=parm(iparm)
      do 20 iparm=1,nparme
   20   zex(iwbase(iparm),1)=parm(nparml+iparm)
c     do 22 iparm=1,nparmd
c  22   cdet(iwdet(iparm),1)=parm(nparml+nparme+iparm)
      do 22 iparm=1,nparmcsf
   22   csf_coef(iwcsf(iparm),1)=parm(nparml+nparme+iparm)

c     if(nparms.eq.1) scalek(1)=parm(nparml+nparme+nparmd+1)
c     if(nparmg.eq.1) a21   =parm(nparml+nparme+nparmd+nparms+1)
      if(nparms.eq.1) scalek(1)=parm(nparml+nparme+nparmcsf+1)
      if(nparmg.eq.1) a21   =parm(nparml+nparme+nparmcsf+nparms+1)
      if(ijas.eq.1) then
        if(nparmj.ge.1) cjas2(1)=parm(nparm)
       elseif(ijas.eq.2) then
        ntmp=nparmj
        do 26 isp=nspin1,nspin2
          do 25 iparm=1,nparma(isp)
   25       a1(iwjasa(iparm,isp),isp,1)=parm(nparm-ntmp+iparm)
   26     ntmp=ntmp-nparma(isp)
        do 28 isp=nspin1,nspin2
          do 27 iparm=1,nparmb(isp)
   27       a2(iwjasb(iparm,isp),isp,1)=parm(nparm-ntmp+iparm)
   28     ntmp=ntmp-nparmb(isp)
       elseif(ijas.eq.3) then
        ntmp=nparmj
        do 29 iparm=1,nparma(1)
   29     a(iwjasa(iparm,1),1)=parm(nparm-ntmp+iparm)
        ntmp=ntmp-nparma(1)
        do 31 isp=nspin1,nspin2b
          do 30 iparm=1,nparmb(isp)
   30       b(iwjasb(iparm,isp),isp,1)=parm(nparm-ntmp+iparm)
   31     ntmp=ntmp-nparmb(isp)
        do 33 it=1,nctype
          do 32 iparm=1,nparmc(it)
   32       c(iwjasc(iparm,it),it,1)=parm(nparm-ntmp+iparm)
   33     ntmp=ntmp-nparmc(it)
        if(ifock.gt.0) then
          do 38 it=1,nctype
            do 37 iparm=1,nparmf(it)
   37         fck(iwjasf(iparm,it),it,1)=parm(nparm-ntmp+iparm)
            ntmp=ntmp-nparmf(it)
            if(ifock.gt.2) then
              call scale3(1,it)
              if(ifock.eq.4) call scale20(1,it)
            endif
   38     continue
        endif
       elseif(ijas.ge.4.and.ijas.le.6) then
        ntmp=nparmj
        do 40 it=1,nctype
          do 39 iparm=1,nparma(it)
   39       a4(iwjasa(iparm,it),it,1)=parm(nparm-ntmp+iparm)
   40     ntmp=ntmp-nparma(it)
        do 42 isp=nspin1,nspin2b
          do 41 iparm=1,nparmb(isp)
   41       b(iwjasb(iparm,isp),isp,1)=parm(nparm-ntmp+iparm)
   42     ntmp=ntmp-nparmb(isp)
        do 44 it=1,nctype
          do 43 iparm=1,nparmc(it)
   43       c(iwjasc(iparm,it),it,1)=parm(nparm-ntmp+iparm)
   44     ntmp=ntmp-nparmc(it)
      endif

      if(icusp2.ge.1.and.isc.le.7) then
        do 45 isp=nspin1,nspin2
   45     if(ijas.eq.2) call cuspexact2(a1(1,isp,1),a2(1,isp,1))
        if(ijas.eq.3) call cuspexact3(0)
        if(ijas.ge.4.and.ijas.le.6) call cuspexact4(0,1)
      endif

c Here we are calculating analytical derivs. wrt. wavefn. params so turn igradhess on before calling hpsi.
      igradhess=1

      nnumerical=nparm-nanalytic
c     do 125 i=1,ndata
      do 125 i=idispls(idtask)+1,idispls(idtask+1)
        iconfg=i
c       call deriv_hpsi(x(1,1,i),psid(i),psij(i),energy,denergy,1)
        call hpsi(x(1,1,i),psid(i),psij(i),velocity,div_v,d2psi,pe,pei,energy,denergy,1)
        do 125 iparm=1,nanalytic
          ajac(i,nnumerical+iparm)=denergy(iparm)
  125   continue

      if(index(mode,'mpi').ne.0) call jacobian_mpi(ndata2,nanalytic,nparm,ajac)

      return
      end
