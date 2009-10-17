      subroutine acuest_dmc_mov1_mpi1
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

# if defined (MPI)

      use all_tools_mod
      use control_mod
      use montecarlo_mod
      use variables_mod
      use mpi_mod
      use atom_mod
      use dets_mod
      use basis1_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use forcepar_mod
      use pseudo_mod
      use contrl_per_mod
      use delocc_mod
      use force_dmc_mod
      use iterat_mod
      use qua_mod
      use jacobsave_mod
      use mpiblk_mod
      use forcest_dmc_mod
      use denupdn_mod
      use stepv_mod
      use config_dmc_mod
      use branch_mod
      use estsum_dmc_mod
      use estcum_dmc_mod
      use div_v_dmc_mod
      use contrldmc_mod
      use estcm2_mod
      use stats_mod
      use age_mod
      use pairden_mod
      use pop_control_mod, only : ffn
      implicit real*8(a-h,o-z)

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4
      common /compferm/ emagv,nv,idot

      dimension egcollect(nforce),wgcollect(nforce),pecollect(nforce),peicollect(nforce),
     &tpbcollect(nforce),tjfcollect(nforce),eg2collect(nforce),wg2collect(nforce),
     &pe2collect(nforce),pei2collect(nforce),tpb2collect(nforce),tjf2collect(nforce),fsum(nforce),
     &f2sum(nforce),eg2sum(nforce),wg2sum(nforce),pe2sum(nforce),pei2sum(nforce),tpb2sum(nforce),
     &tjf2sum(nforce),taucollect(nforce)

c statement function for error calculation
      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

c wt   = weight of configurations
c xsum = sum of values of x from dmc
c xnow = average of values of x from dmc
c xcum = accumulated sums of xnow
c xcm2 = accumulated sums of xnow**2
c xave = current average value of x
c xerr = current error of x

      iblk=iblk+1
      iblk_proc=iblk_proc+nproc

      npass=iblk_proc*nstep

c     wnow=wsum/nstep
c     wfnow=wfsum/nstep
      enow=esum/wsum
      efnow=efsum/wfsum
      ei1now=wfsum/wdsum
      ei2now=wgsum(1)/wgdsum
      rinow=risum/wgsum(1)
      r2now=r2sum/wgsum(1)

      ei1cm2=ei1cm2+ei1now**2
      ei2cm2=ei2cm2+ei2now**2
      r2cm2=r2cm2+r2sum*r2now
      ricm2=ricm2+risum*rinow

      wdcum=wdcum+wdsum
      wgdcum=wgdcum+wgdsum
      ei1cum=ei1cum+ei1now
      ei2cum=ei2cum+ei2now
      r2cum=r2cum+r2sum
      ricum=ricum+risum

      w2sum=wsum**2
      wf2sum=wfsum**2
      e2sum=esum*enow
      ef2sum=efsum*efnow

      do 10 ifr=1,nforce
c       wgnow=wgsum(ifr)/nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum(ifr)/wgsum(ifr)
        peinow=peisum(ifr)/wgsum(ifr)
        tpbnow=tpbsum(ifr)/wgsum(ifr)
        tjfnow=tjfsum(ifr)/wgsum(ifr)

        wg2sum(ifr)=wgsum(ifr)**2
        eg2sum(ifr)=egsum(ifr)*egnow
        pe2sum(ifr)=pesum(ifr)*penow
        pei2sum(ifr)=peisum(ifr)*peinow
        tpb2sum(ifr)=tpbsum(ifr)*tpbnow
        tjf2sum(ifr)=tjfsum(ifr)*tjfnow

        if(ifr.gt.1) then
          fsum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))
          f2sum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
        endif
  10  continue

      call mpi_allreduce(wgsum,wgcollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(egsum,egcollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tausum,taucollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 12 ifr=1,nforce
c Warning temp fix
        wgsum(ifr)=wgcollect(ifr)
        egsum(ifr)=egcollect(ifr)
        egnow=egsum(ifr)/wgsum(ifr)
        wgcum(ifr)=wgcum(ifr)+wgcollect(ifr)
        egcum(ifr)=egcum(ifr)+egcollect(ifr)
        taucum(ifr)=taucum(ifr)+taucollect(ifr)
  12  continue

      walker_weights_sum_block = wgsum(1)
      walker_weights_sum = wgcum(1)
      call object_modified_by_index (walker_weights_sum_block_index) !JT
      call object_modified_by_index (walker_weights_sum_index) !JT

      call mpi_allreduce(pesum,pecollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(peisum,peicollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tpbsum,tpbcollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tjfsum,tjfcollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      call mpi_allreduce(wg2sum,wg2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(eg2sum,eg2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(pe2sum,pe2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(pei2sum,pei2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tpb2sum,tpb2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tjf2sum,tjf2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      call mpi_allreduce(fsum,fcollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(f2sum,f2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      call mpi_allreduce(esum,ecollect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wsum,wcollect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(efsum,efcollect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wfsum,wfcollect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      call mpi_allreduce(e2sum,e2collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(w2sum,w2collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ef2sum,ef2collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wf2sum,wf2collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

!JT   if(idtask.ne.0) goto 17 ! The slaves also have to calculate egerr so that they can stop when egerr < error_threshold

      wcm2=wcm2+w2collect
      wfcm2=wfcm2+wf2collect
      ecm2=ecm2+e2collect
      efcm2=efcm2+ef2collect

      wcum=wcum+wcollect
      wfcum=wfcum+wfcollect
      ecum=ecum+ecollect
      efcum=efcum+efcollect

      do 15 ifr=1,nforce
        wgcm2(ifr)=wgcm2(ifr)+wg2collect(ifr)
        egcm2(ifr)=egcm2(ifr)+eg2collect(ifr)
        pecm2(ifr)=pecm2(ifr)+pe2collect(ifr)
        peicm2(ifr)=peicm2(ifr)+pei2collect(ifr)
        tpbcm2(ifr)=tpbcm2(ifr)+tpb2collect(ifr)
        tjfcm2(ifr)=tjfcm2(ifr)+tjf2collect(ifr)

        pecum(ifr)=pecum(ifr)+pecollect(ifr)
        peicum(ifr)=peicum(ifr)+peicollect(ifr)
        tpbcum(ifr)=tpbcum(ifr)+tpbcollect(ifr)
        tjfcum(ifr)=tjfcum(ifr)+tjfcollect(ifr)

        call grad_hess_jas_cum(wgsum(ifr),egnow)

        if(iblk.eq.1) then
          egerr=0
          peerr=0
          peierr=0
          tpberr=0
          tjferr=0
         else
          egerr=errg(egcum(ifr),egcm2(ifr),ifr)
          peerr=errg(pecum(ifr),pecm2(ifr),ifr)
          peierr=errg(peicum(ifr),peicm2(ifr),ifr)
          tpberr=errg(tpbcum(ifr),tpbcm2(ifr),ifr)
          tjferr=errg(tjfcum(ifr),tjfcm2(ifr),ifr)
        endif

        egave=egcum(ifr)/wgcum(ifr)
        peave=pecum(ifr)/wgcum(ifr)
        peiave=peicum(ifr)/wgcum(ifr)
        tpbave=tpbcum(ifr)/wgcum(ifr)
        tjfave=tjfcum(ifr)/wgcum(ifr)

        if(ifr.eq.1) then                                !JT
         eloc_bav = egnow                                !JT
         eloc_av = egave                                 !JT
         call object_modified_by_index (eloc_bav_index)  !JT
         call object_modified_by_index (eloc_av_index)   !JT
        endif                                            !JT

        if(ifr.gt.1) then
          fgcum(ifr)=fgcum(ifr)+fsum(ifr)
          fgcm2(ifr)=fgcm2(ifr)+f2sum(ifr)
          fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
          if(iblk.eq.1) then
            fgerr=0
           else
            fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
          endif
          if(deltot(ifr).ne.0.d0) then
            fgave=fgave/abs(deltot(ifr))
            fgerr=fgerr/abs(deltot(ifr))
          endif
        endif

c write out header first time

        if(iblk.eq.1.and.ifr.eq.1) then
          if(ndim.eq.2) then
            write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66
     &      ,''tjfave'',t72,''(tjferr)'',t83,''emave'',t89,''(emerr)'',t100,''fgave'',t106,''(fgerr)'',
     &      t118,''npass'',t128,''wgsum'',t138,''ioldest'')')
           else
            write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66
     &      ,''tjfave'',t72,''(tjferr)'',t83,''fgave'',t89,''(fgerr)'',t101,''npass'',t111,''wgsum'',t121,''ioldest'')')
          endif
        endif

c write out current values of averages etc.

        if(ndim.eq.2) then
          iegerr=nint(10000000*egerr)
          ipeerr=nint(10000000*peerr)
          ipeierr=nint(10000000*peierr)
          itpber=nint(10000000*tpberr)
          itjfer=nint(10000000*tjferr)
          if(ifr.gt.1) ifgerr=nint(10000000*fgerr)
         else
          iegerr=nint(100000*egerr)
          ipeerr=nint(100000*peerr)
          itpber=nint(100000*tpberr)
          itjfer=nint(100000*tjferr)
          if(ifr.gt.1) ifgerr=nint(100000*fgerr)
        endif

c magnetic energy for quantum dots...
c right definition of the potential energy does not include magnetic energy.
        if(ndim.eq.2) then
c         emave=0.125*bext*bext*r2cum/wgcum(ifr)+emaglz+emagsz+emagv
          temp=0.25d0*bext*bext/(we*we)
          emave=(peave-peiave-emag)*temp+emag
          emerr=(peerr+peierr)*temp
          iemerr=nint(10000000*emerr)
          peave=peave-emave
c         ipeerr=ipeerr+iemerr
          ipeerr=nint(10000000*(peerr*(1-temp)+temp*peierr))
        endif

        if(ifr.eq.1) then
          if(ndim.eq.2) then
            write(6,'(f12.7,5(f12.7,''('',i7,'')''),17x,3i10)') egcollect(ifr)/wgcollect(ifr),
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,emave,iemerr,npass,nint(wgcollect(ifr)/nproc),ioldest
           else
            write(6,'(f10.5,4(f10.5,''('',i5,'')''),17x,3i10)') egcollect(ifr)/wgcollect(ifr),
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,npass,nint(wgcollect(ifr)/nproc),ioldest
          endif
         else
          if(ndim.eq.2) then
            write(6,'(f12.7,5(f12.7,''('',i7,'')''),17x,3i10)') egcollect(ifr)/wgcollect(ifr),
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,emave,iemerr,fgave,ifgerr,nint(wgcollect(ifr)/nproc)
           else
            write(6,'(f10.5,5(f10.5,''('',i5,'')''),10x,i10)') egcollect(ifr)/wgcollect(ifr),
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,fgave,ifgerr,nint(wgcollect(ifr)/nproc)
          endif
        endif
   15 continue

!      moved up
!      eloc_av = egave                                 !JT
!      call object_modified_by_index (eloc_av_index)   !JT


c I have changed the dwt limit in the dmc routines so it does not depend on etrial so  there is no need for these warning msgs.
c The dwt limit is there to prevent population explosions with nonlocal psps. but there are
c better solutions than dwt limits.
c     if(wgsum(1).gt.1.5d0*nstep*nconf*nproc .or. (iblk.gt.2*nblkeq+5 .and. etrial .gt. egave+200*egerr))
c    &write(6,'(''Warning: etrial too high? It should be reasonably close to DMC energy because of dwt in dmc'')')
c     if(wgsum(1).lt.0.7d0*nstep*nconf*nproc .or. (iblk.gt.2*nblkeq+5 .and. etrial .lt. egave-200*egerr))
c    &write(6,'(''Warning: etrial too low?  It should be reasonably close to DMC energy because of dwt in dmc'')')

      call systemflush(6)

c zero out xsum variables for metrop

   17 wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r2sum=zero
      risum=zero

      do 20 ifr=1,nforce
        wgsum(ifr)=zero
        egsum(ifr)=zero
        pesum(ifr)=zero
        peisum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
   20   tausum(ifr)=zero

      return

      entry acues1_dmc_mov1_mpi1
c statistical fluctuations without blocking

      if(ipr.gt.-2) then
         write(11,'(i8,f11.8,f15.8,f13.8,i5)') ipass,ffn,wsum1(1),esum1(1)/wsum1(1),nwalk
      end if


      wdsum1=wdsumo
      wgdsum1=wgdsumo

      wcum1=wcum1+wsum1(1)
      wfcum1=wfcum1+wfsum1
      ecum1=ecum1+esum1(1)
      efcum1=efcum1+efsum1
      ei3cum=ei3cum+wfsum1/wdsum1

      wcm21=wcm21+wsum1(1)**2
      wfcm21=wfcm21+wfsum1**2
      ecm21=ecm21+esum1(1)**2/wsum1(1)
      efcm21=efcm21+efsum1**2/wfsum1
      ei3cm2=ei3cm2+(wfsum1/wdsum1)**2
      do 30 ifr=1,nforce
        wgcum1(ifr)=wgcum1(ifr)+wgsum1(ifr)
        egcum1(ifr)=egcum1(ifr)+egsum1(ifr)
        wgcm21(ifr)=wgcm21(ifr)+wgsum1(ifr)**2
        if(wgsum1(ifr).ne.0.d0) then
          egcm21(ifr)=egcm21(ifr)+egsum1(ifr)**2/wgsum1(ifr)
         else
          egcm21(ifr)=0
        endif
   30 continue

c collect block averages
      wsum=wsum+wsum1(1)
      wfsum=wfsum+wfsum1
      wdsum=wdsum+wdsumo
      wgdsum=wgdsum+wgdsum1
      esum=esum+esum1(1)
      efsum=efsum+efsum1
c     eisum=eisum+wfsum1/wdsum1
      do 35 ifr=1,nforce
        wgsum(ifr)=wgsum(ifr)+wgsum1(ifr)
   35   egsum(ifr)=egsum(ifr)+egsum1(ifr)

c Estimate eigenvalue of G from the energy
      ipmod=mod(ipass,nfprod)
      if(iabs(idmc).eq.1) then
        nfpro=min(nfprod,ipass)
        eigv=(wgsum1(1)/wtgen(ipmod))**(one/nfpro)
       else
        eest=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
        eigv=dexp((etrial-eest)*(taucum(1)+tausum(1))/(wgcum(1)+wgsum(1)))
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest,egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wsum1(1)
      wgdsumo=wsum1(1)*fprod/ff(mod(ipass+1,nfprod))
      wtgen(ipmod)=wsum1(1)

c zero out step averages
      wfsum1=zero
      wdsum1=zero
      efsum1=zero
      do 40 ifr=1,nforce
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        esum1(ifr)=zero
   40   egsum1(ifr)=zero

      return

      entry zeres0_dmc_mov1_mpi1
c Initialize various quantities at beginning of run
c the initial values of energy psi etc. are calculated here

      ipass=0

c set quadrature points

c     if(nloc.gt.0) call gesqua(nquad,xq,yq,zq,wq)
      if(nloc.gt.0) call rotqua

      eigv=one
      eest=etrial
      nwalk=nconf
      wdsumo=nconf_global
      wgdsumo=nconf_global
      fprod=one

      call object_modified_by_index (nwalk_index)

      do 70 i=0,nfprod
        wtgen(i)=nconf_global
  70    ff(i)=one

! JT: it seems that the code remains stuck around here runs when compiled with ifort 10.1 with optimization option -O3.
! It works with optimization option -O2. It also works when a write statement is added in the loop as done below!
!     write(6,'(/,''These lines are printed out just because otherwise it gets stuck here with ifort 10.1 -O3'')')
      do 80 iw=1,nconf
!       write(6,'(''iw='',i4)') iw
        wt(iw)=one
        if(istrech.eq.0) then
          do 71 ifr=2,nforce
            do 71 ie=1,nelec
              do 71 k=1,ndim
   71           xoldw(k,ie,iw,ifr)=xoldw(k,ie,iw,1)
        endif
        do 72 ifr=1,nforce
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
          pwt(iw,ifr)=one
          do 72 ip=0,nwprod-1
   72       wthist(iw,ip,ifr)=one
   80 continue

      entry zerest_dmc_mov1_mpi1
c entry point to zero out all averages etc. after equilibration runs

      iblk=0
      iblk_proc=0

c zero out estimators

      wcum1=zero
      wfcum1=zero
      wcum=zero
      wfcum=zero
      wdcum=zero
      wgdcum=zero
      ecum1=zero
      efcum1=zero
      ecum=zero
      efcum=zero
      ei1cum=zero
      ei2cum=zero
      ei3cum=zero
      r2cum=zero
      ricum=zero

      wcm21=zero
      wfcm21=zero
      wcm2=zero
      wfcm2=zero
      wdcm2=zero
      wgdcm2=zero
      ecm21=zero
      efcm21=zero
      ecm2=zero
      efcm2=zero
      ei1cm2=zero
      ei2cm2=zero
      ei3cm2=zero
      r2cm2=zero
      ricm2=zero

      wfsum1=zero
      wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      efsum1=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      ei3sum=zero
      r2sum=zero
      risum=zero

      call grad_hess_jas_init

      call alloc ('fgcum', fgcum, nforce)
      call alloc ('fgcm2', fgcm2, nforce)
      call alloc ('wgcm2', wgcm2, nforce)
      call alloc ('wgcm21', wgcm21, nforce)
      call alloc ('egcm2', egcm2, nforce)
      call alloc ('egcm21', egcm21, nforce)
      call alloc ('pecm2', pecm2, nforce)
      call alloc ('tpbcm2', tpbcm2, nforce)
      call alloc ('tjfcm2', tjfcm2, nforce)
      call alloc ('peicm2', peicm2, nforce)
      call alloc ('wgcum', wgcum, nforce)
      call alloc ('wgcum1', wgcum1, nforce)
      call alloc ('egcum', egcum, nforce)
      call alloc ('egcum1', egcum1, nforce)
      call alloc ('pecum', pecum, nforce)
      call alloc ('peicum', peicum, nforce)
      call alloc ('tpbcum', tpbcum, nforce)
      call alloc ('tjfcum', tjfcum, nforce)
      call alloc ('taucum', taucum, nforce)
      call alloc ('wgsum', wgsum, nforce)
      call alloc ('wsum1', wsum1, nforce)
      call alloc ('wgsum1', wgsum1, nforce)
      call alloc ('egsum', egsum, nforce)
      call alloc ('egsum1', egsum1, nforce)
      call alloc ('pesum', pesum, nforce)
      call alloc ('peisum', peisum, nforce)
      call alloc ('tpbsum', tpbsum, nforce)
      call alloc ('tjfsum', tjfsum, nforce)
      call alloc ('tausum', tausum, nforce)
      call alloc ('esum1', esum1, nforce)
c Do it for MFORCE rather than nforce because in optimization at the start nforce=1 but later nforce=3
!JT: this should not be necessary and it is annoying for dynamic allocation!
      do 85 ifr=1,nforce
!JT     do 85 ifr=1,MFORCE
        tausum(ifr)=zero
        taucum(ifr)=zero
        wgcum1(ifr)=zero
        wgcum(ifr)=zero
        egcum1(ifr)=zero
        egcum(ifr)=zero
        wgcm21(ifr)=zero
        wgcm2(ifr)=zero
        egcm21(ifr)=zero
        egcm2(ifr)=zero
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        wgsum(ifr)=zero
        esum1(ifr)=zero
        egsum1(ifr)=zero
        egsum(ifr)=zero
        pecum(ifr)=zero
        peicum(ifr)=zero
        tpbcum(ifr)=zero
        tjfcum(ifr)=zero
        pecm2(ifr)=zero
        peicm2(ifr)=zero
        tpbcm2(ifr)=zero
        tjfcm2(ifr)=zero
        pesum(ifr)=zero
        peisum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
        fgcum(ifr)=zero
   85   fgcm2(ifr)=zero

      nbrnch=0

      try_int=0
      acc=0
      acc_int=0
      nodecr=0

c Zero out estimators for charge density of atom.
      do 90 i=1,NRAD
        rprobup(i)=zero
        rprobdn(i)=zero
   90   rprob(i)=zero

c Zero out estimators for pair densities:
      if (ifixe.ne.0) then
      allocate (den2d_t(-NAX:NAX,-NAX:NAX))
      allocate (den2d_u(-NAX:NAX,-NAX:NAX))
      allocate (den2d_d(-NAX:NAX,-NAX:NAX))
      allocate (xx0probdt(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probdu(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probdd(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probut(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probuu(0:NAX,-NAX:NAX,-NAX:NAX))
      allocate (xx0probud(0:NAX,-NAX:NAX,-NAX:NAX))
      do 100 i2=-NAX,NAX
        do 100 i3=-NAX,NAX
          den2d_t(i2,i3)=0
          den2d_u(i2,i3)=0
          den2d_d(i2,i3)=0
          do 100 i1=0,NAX
            xx0probdt(i1,i2,i3)=0
            xx0probdu(i1,i2,i3)=0
            xx0probdd(i1,i2,i3)=0
            xx0probut(i1,i2,i3)=0
            xx0probuu(i1,i2,i3)=0
  100       xx0probud(i1,i2,i3)=0
      endif

      call grad_hess_jas_save
# endif

      return
      end
