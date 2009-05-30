      subroutine jacobian(ndata2,nparm,nanalytic,parm,ajac)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use basic_tools_mod
      use fitdet_mod
      use atom_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use basis1_mod
      use contr2_mod
      use contrl_opt2_mod
      use delocc_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use bparm_mod
      use contr3_mod
      use pars_mod
      use jaspar1_mod
      use jaspar2_mod
      use confg_mod
      use const_mod
      implicit real*8(a-h,o-z)

c epsder1f=sqrt(max(eps_diff,dbl_epsilon))
c     parameter (eps_diff=1.d-15,dbl_epsilon=2.2204460492503131d-16
c    &,epsder1f=3.16227766016838d-8)

      common /mpioffset/ ircounts(0:MPROC),idispls(0:MPROC)

c     dimension parm(nparm),ajac(ndata2,nparm),denergy(MPARM)
      dimension parm(nparm),ajac(ndata2,nparm),velocity(3,nelec),div_v(nelec)

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
