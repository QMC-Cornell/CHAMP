      function func(ndata2,nparm,parm,diff,iflag)
c Written by Cyrus Umrigar, modified by Claudia Filippi
c If iflag=0 then it computes diffs otherwise it just does sum of squares
      use constants_mod
      use basic_tools_mod
      use fitdet_mod
      use atom_mod
      use coefs_mod
      use dets_mod
      use optim_mod
      use basis1_mod
      use numbas_mod
      use basis2_mod
      use contr2_mod
      use contrl_opt2_mod
      use pseudo_mod
      use jaspar_mod
      use jaspar3_mod
      use jaspar4_mod
      use bparm_mod
      use contr3_mod
      use pars_mod
      use jaspar1_mod
      use jaspar2_mod
      use ncusp_mod
      implicit real*8(a-h,o-z)

c      complex*16 cvd_sav,cvk_sav

      common /confg/ x(3,MELEC,MDATA),eguess,psid(MDATA),psij(MDATA),
     &psio(MDATA),eold(MDATA),uwdiff(MDATA),wght(MDATA),wghtsm,cuspwt,
     &dvpdv(MDATA),ndata

c     common /wcsf/ frac(ICX,MDET),icsf(ICSFX)

      common /fcn_calls/icalls
c      common /fitdet/ cvd_sav(3,MELEC,MDATA),vd_sav(3,MELEC,MDATA),psid_sav(MDATA)
c     &               ,d2d_sav(MDATA),div_vd_sav(MELEC,MDATA),cvk_sav(3,MELEC,MDATA),psik_sav(MDATA)
c     &               ,div_vk_sav(MELEC,MDATA),d2k_sav(MDATA),iconfg,isaved


      common /mpioffset/ ircounts(0:MPROC),idispls(0:MPROC)

      dimension velocity(3,MELEC)
     &,div_v(MELEC)

      dimension parm(*),diff(*)

      if(iflag.eq.0) then

      icalls=icalls+1

      do 10 iparm=1,nparml
   10   coef(iwbasi(iparm),iworb(iparm),1)=parm(iparm)
      do 20 iparm=1,nparme
        zex(iwbase(iparm),1)=parm(nparml+iparm)
        ict=ictype_basis(iwbase(iparm))
        irb=iwrwf2(iwbase(iparm))
   20   zex2(irb,ict,1)=zex(iwbase(iparm),1)
c     do 22 iparm=1,nparmd
c  22   cdet(iwdet(iparm),1)=parm(nparml+nparme+iparm)
      do 22 iparm=1,nparmcsf
   22   csf_coef(iwcsf(iparm),1)=parm(nparml+nparme+iparm)

cc Calculate coefs to construct the piece of the orbital that comes
cc from basis fns that are related by symmetry.  This is needed to
cc impose cusp conditions when there is more than one atom.
c      if(nloc.eq.0.and.numr.le.0) call equiv_bas

      if(nparms.eq.1) then
c       scalek(1)=parm(nparml+nparme+nparmd+1)
        scalek(1)=parm(nparml+nparme+nparmcsf+1)
c If we are varying scalek reset dependent constants
c       if(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17) then
c         call set_scale_dist(0,1)
c       endif
      endif

c     if(nparmg.eq.1) a21   =parm(nparml+nparme+nparmd+nparms+1)
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
        write(6,'(''fck'',i3,f14.8)') iparm,parm(nparm-ntmp+iparm)
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

c Foll. is now in cuspco
c     if(ijas.eq.2) then
c       if(nspin1.eq.1 .and. nspin2.eq.1) then
c         if(nelec.eq.2) then
c           aa1=a1(2,1,1)
c          else
c           aa1=d1b4*(three*(nup+ndn)-2)*a1(2,1,1)
c         endif
c        elseif(nspin1.eq.2 .and. nspin2.eq.2) then
c         aa1=(nup-1)*a1(2,2,1)
c        elseif(nspin1.eq.1 .and. nspin2.eq.2) then
c         aa1=half*((nup+ndn)*a1(2,1,1)+(nup+ndn-2)*a1(2,2,1))
c        elseif(nspin1.eq.1 .and. nspin2.eq.3) then
c         aa1=nup*a1(2,1,1)+(ndown-1)*a1(2,3,1)
c       endif
c      elseif(ijas.eq.3) then
c       aa1=a(1,1)
c      elseif(ijas.ge.4.and.ijas.le.6) then
c       aa1=a4(1,1,1)
c      else
c       aa1=zero
c     endif

c If we are varying scalek, a or b parms reset dependent constants
c      if((ijas.eq.4.or.ijas.eq.5).and.(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17)) then
        call set_scale_dist(0,1)    ! always call this in the analytical scalek opt
c      endif

      if(icusp2.ge.1.and.isc.le.7) then
        do 45 isp=nspin1,nspin2
   45     if(ijas.eq.2) call cuspexact2(a1(1,isp,1),a2(1,isp,1))
        if(ijas.eq.3) call cuspexact3(0)
        if(ijas.ge.4.and.ijas.le.6) call cuspexact4(0,1)
      endif

      if(iabs(icusp2).ge.2) then
        if(ijas.eq.2) then
          do 46 isp=nspin1,nspin2
            ishft=ncuspc*(isp-nspin1)
            if(nspin2.eq.2. and. (nup-ndn).ne.0) a1(2,2,1)=a1(2,1,1)
            if(nspin2.eq.3 .and. isp.eq.nspin2)
     &      a1(2,3,1)=((ndn-nup)*a1(2,1,1)+(nup-1)*a1(2,2,1))/(ndn-1)
            call cuspcheck2(scalek(1),a1(1,isp,1),a2(1,isp,1),
     &      diff(ndata+ishft+1),isp,nspin1,ncuspc,0)
            do 46 i=1,ncuspc
   46         diff(ndata+ishft+i)=diff(ndata+ishft+i)*cuspwt
         elseif(ijas.eq.3) then
          call cuspcheck3(diff(ndata+1),0)
          do 47 i=1,ncuspc+nfockc
   47       diff(ndata+i)=diff(ndata+i)*cuspwt
        endif
      endif

c If icusp>=0 impose cusp condition for each s orbital.
c In any case calculate cusp-violation penalty. Should be 0 if icusp>=0.
      ishft=ncuspc*(nspin2-nspin1+1)+nfockc
      if((nloc.eq.0. .or. nloc.eq.6) .and. numr.le.0) call cuspco(diff(ndata+ishft+1),0)

      do 50 i=1,ncent*norbc
   50   diff(ndata+ishft+i)=diff(ndata+ishft+i)*cuspwt
      do 60 i=1,necn
   60   coef(iebasi(1,i),ieorb(1,i),1)=sign(one,dfloat(ieorb(2,i)))*
     &  coef(iebasi(2,i),iabs(ieorb(2,i)),1)
      do 70 i=1,nebase
        zex(iebase(1,i),1)=zex(iebase(2,i),1)
        ict=ictype_basis(iebase(1,i))
        irb=iwrwf2(iebase(1,i))
   70   zex2(irb,ict,1)=zex(iebase(2,i),1)

c     do 80 i=1,iabs(nedet)
c       cdet(iedet(1,i),1)=zero
c       do 80 j=2,icsf(i)
c  80     cdet(iedet(1,i),1)=cdet(iedet(1,i),1)+frac(j,i)*
c    &    sign(one,dfloat(iedet(j,i)))*cdet(iabs(iedet(j,i)),1)

      if(ipos+idcds+idcdu+idcdt+id2cds+id2cdu+id2cdt+idbds+idbdu+idbdt.
     &gt.0) then
        ishft=ndata+ncuspc*(nspin2-nspin1+1)+nfockc+ncent*norbc+1
        do 121 isp=nspin1,nspin2
          call checkjas2(scalek(1),isp,ncnstr,diff(ishft),0,0)
  121     ishft=ishft+ncnstr
      endif

c Data distribution for parallel run.  Modified by RGH.
      isize=ndata/nproc
      imod=mod(ndata,nproc)

      if(imod.eq.0) then
        do i=0,nproc-1
          ircounts(i)=isize
          idispls(i)=i*isize
        enddo
       else
        do i=0,imod-1
          ircounts(i)=isize+1
          idispls(i)=i*(isize+1)
        enddo
        do i=imod,nproc-1
          ircounts(i)=isize
          idispls(i)=i*isize+imod
        enddo
      endif
      idispls(nproc)=ndata
c RGH: End modification

c Here we are calculating numerical derivs. wrt. wavefn. params so turn igradhess off before calling hpsi.
      igradhess=0

c For a serial run this reduces to a "do 123 i=1,ndata"
      do 123 i=idispls(idtask)+1,idispls(idtask+1)
        iconfg=i
        call hpsi(x(1,1,i),psid(i),psij(i),velocity,div_v,d2psi,pe,pei,energy,denergy,1)
  123   diff(i)=energy-eguess
c 123   uwdiff(i)=diff(i)
      if(nparml+nparme+nparmcsf.eq.0) isaved=1

      if(index(mode,'mpi').ne.0) call func_qmc_mpi(diff)

      do 125 i=1,ndata
  125   uwdiff(i)=diff(i)

      if(mod(irewgt,100).ne.0) call rewght(diff)

      eavr=zero
      wavr=zero
      do 127 i=1,ndata
        wavr=wavr+wght(i)
  127   eavr=eavr+(uwdiff(i)+eguess)*wght(i)
      eavr=eavr/wavr

      eavri=one
      if(irewgt.ge.100) eavri=-one/eavr
      do 128 i=1,ndata
  128   diff(i)=diff(i)*eavri

c     err=zero
c     do 130 i=1,ndata2
c 130   err=err+diff(i)**2
c     func=err

      else

      err=zero
      do 130 i=1,ndata2
  130   err=err+diff(i)**2
      func=err

      endif

      return
      end
