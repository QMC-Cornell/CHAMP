      subroutine acuest_dmc_mov1_mpi2
! MPI version created by Claudia Filippi starting from serial version
! routine to accumulate estimators for energy etc.

# if defined (MPI)

      use all_tools_mod
      use control_mod
      use montecarlo_mod
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
      use fourier_mod
      use pop_control_mod, only : ffn
      use zigzag_mod
      implicit real*8(a-h,o-z)



      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot

      dimension pecollect(nforce),peicollect(nforce), &
     &tpbcollect(nforce),tjfcollect(nforce),taucollect(nforce), &
     &collect(2*nforce+5),collect_t(2*nforce+5), &
     &zzsum_collect(nzzvars),zznow(nzzvars)

! statement function for error calculation
      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

! wt   = weight of configurations
! xsum = sum of values of x from dmc
! xnow = average of values of x from dmc
! xcum = accumulated sums of xnow
! xcm2 = accumulated sums of xnow**2
! xave = current average value of x
! xerr = current error of x

      iblk=iblk+1

      npass=iblk*nstep

      call mpi_reduce(pesum,pecollect,nforce,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(peisum,peicollect,nforce,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpbsum,tpbcollect,nforce,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tjfsum,tjfcollect,nforce,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
!      call mpi_reduce(tausum,taucollect,nforce,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tausum,taucollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ioldest,ioldest_collect,1,mpi_integer,mpi_max,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ioldestmx,ioldestmx_collect,1,mpi_integer,mpi_max,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(r1sum,r1sum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(r2sum,r2sum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(r3sum,r3sum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(r4sum,r4sum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(risum,risum_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      if(izigzag.gt.0) then
       call mpi_allreduce(zzsum,zzsum_collect,nzzvars,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      endif

      ioldest=ioldest_collect
      ioldestmx=ioldestmx_collect

      call mpi_barrier(MPI_COMM_WORLD,ierr)

!     if(idtask.ne.0) goto 17 ! The slaves also have to calculate egerr so that they can stop when egerr < error_threshold

!     wnow=wsum/nstep
!     wfnow=wfsum/nstep
      r1sum = r1sum_collect
      r2sum = r2sum_collect
      r3sum = r3sum_collect
      r4sum = r4sum_collect
      risum = risum_collect
      if(izigzag.gt.0) then
       zzsum(:) = zzsum_collect(:)
      endif

      enow=esum/wsum
      efnow=efsum/wfsum
      ei1now=wfsum/wdsum
      ei2now=wgsum(1)/wgdsum
      rinow=risum/wgsum(1)
      r1now=r1sum/wgsum(1)
      r2now=r2sum/wgsum(1)
      r3now=r3sum/wgsum(1)
      r4now=r4sum/wgsum(1)
      if(izigzag.gt.0) then
       zznow(:)=zzsum(:)/wgsum(1)
      endif

      ei1cm2=ei1cm2+ei1now**2
      ei2cm2=ei2cm2+ei2now**2
      r1cm2=r1cm2+r1sum*r1now
      r2cm2=r2cm2+r2sum*r2now
      r3cm2=r3cm2+r3sum*r3now
      r4cm2=r4cm2+r4sum*r4now
      ricm2=ricm2+risum*rinow
      if(izigzag.gt.0) then
       zzcm2(:)=zzcm2(:)+zzsum(:)*zznow(:)
      endif

      wdcum=wdcum+wdsum
      wgdcum=wgdcum+wgdsum
      ei1cum=ei1cum+ei1now
      ei2cum=ei2cum+ei2now
      r1cum=r1cum+r1sum
      r2cum=r2cum+r2sum
      r3cum=r3cum+r3sum
      r4cum=r4cum+r4sum
      ricum=ricum+risum
      if(izigzag.gt.0) then
       zzcum(:)=zzcum(:)+zzsum(:)
      endif

      wcm2=wcm2+wsum**2
      wfcm2=wfcm2+wfsum**2
      ecm2=ecm2+esum*enow
      efcm2=efcm2+efsum*efnow

      wcum=wcum+wsum
      wfcum=wfcum+wfsum
      ecum=ecum+esum
      efcum=efcum+efsum

      do 15 ifr=1,nforce

        pesum(ifr)=pecollect(ifr)
        peisum(ifr)=peicollect(ifr)
        tpbsum(ifr)=tpbcollect(ifr)
        tjfsum(ifr)=tjfcollect(ifr)
        tausum(ifr)=taucollect(ifr)

!       wgnow=wgsum(ifr)/nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum(ifr)/wgsum(ifr)
        peinow=peisum(ifr)/wgsum(ifr)
        tpbnow=tpbsum(ifr)/wgsum(ifr)
        tjfnow=tjfsum(ifr)/wgsum(ifr)

        wgcm2(ifr)=wgcm2(ifr)+wgsum(ifr)**2
        egcm2(ifr)=egcm2(ifr)+egsum(ifr)*egnow
        pecm2(ifr)=pecm2(ifr)+pesum(ifr)*penow
        peicm2(ifr)=peicm2(ifr)+peisum(ifr)*peinow
        tpbcm2(ifr)=tpbcm2(ifr)+tpbsum(ifr)*tpbnow
        tjfcm2(ifr)=tjfcm2(ifr)+tjfsum(ifr)*tjfnow

        wgcum(ifr)=wgcum(ifr)+wgsum(ifr)
        egcum(ifr)=egcum(ifr)+egsum(ifr)
        pecum(ifr)=pecum(ifr)+pesum(ifr)
        peicum(ifr)=peicum(ifr)+peisum(ifr)
        tpbcum(ifr)=tpbcum(ifr)+tpbsum(ifr)
        tjfcum(ifr)=tjfcum(ifr)+tjfsum(ifr)
        taucum(ifr)=taucum(ifr)+tausum(ifr)

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

        eest=egave
        eigv=dexp((etrial-eest)*taucum(1)/wgcum(1)) !TA

        if(ifr.gt.1) then
          fgcum(ifr)=fgcum(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))
          fgcm2(ifr)=fgcm2(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
          fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
          if(iblk.eq.1) then
!           ifgerr=0
            fgerr=0
           else
            fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
          endif
          if(deltot(ifr).ne.0.d0) then
            fgave=fgave/abs(deltot(ifr))
            fgerr=fgerr/abs(deltot(ifr))
          endif
        endif

! write out header first time

        if(iblk.eq.1.and.ifr.eq.1) then
          if(ndim.eq.2) then
            write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66 &
     &      ,''tjfave'',t72,''(tjferr)'',t83,''emave'',t89,''(emerr)'',t100,''fgave'',t106,''(fgerr)'', &
     &      t118,''npass'',t128,''wgsum'',t138,''ioldest'')')
           else
            write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66 &
     &      ,''tjfave'',t72,''(tjferr)'',t83,''fgave'',t89,''(fgerr)'',t101,''npass'',t111,''wgsum'',t121,''ioldest'')')
          endif
        endif

! write out current values of averages etc.

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

! magnetic energy for quantum dots...
! right definition of the potential energy does not include magnetic energy.
        if(ndim.eq.2) then
!         emave=0.125*bext*bext*r2cum/wgcum(ifr)+emaglz+emagsz+emagv
          temp=0.25d0*bext*bext/(we*we)
          emave=(peave-peiave-emag)*temp+emag
          emerr=(peerr+peierr)*temp
          iemerr=nint(10000000*emerr)
          peave=peave-emave
!         ipeerr=ipeerr+iemerr
          ipeerr=nint(10000000*(peerr*(1-temp)+temp*peierr))
        endif

        if(ifr.eq.1) then
          if(ndim.eq.2) then
            write(6,'(f12.7,5(f12.7,''('',i7,'')''),17x,3i10)') egsum(ifr)/wgsum(ifr), &
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,emave,iemerr,npass,nint(wgsum(ifr)/nproc),ioldest
           else
            write(6,'(f10.5,4(f10.5,''('',i5,'')''),17x,3i10)') egsum(ifr)/wgsum(ifr), &
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,npass,nint(wgsum(ifr)/nproc),ioldest
          endif
         else
          if(ndim.eq.2) then
            write(6,'(f12.7,5(f12.7,''('',i7,'')''),17x,3i10)') egsum(ifr)/wgsum(ifr), &
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,emave,iemerr,fgave,ifgerr,nint(wgsum(ifr)/nproc)
           else
            write(6,'(f10.5,5(f10.5,''('',i5,'')''),10x,i10)') egsum(ifr)/wgsum(ifr), &
     &      egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,fgave,ifgerr,nint(wgsum(ifr)/nproc)
          endif
        endif
   15 continue

! The reason why having etrial far from e_DMC creates a bias is that in dwt we use taunow, but in
! removing the choice of etrial from dwt we use taueff.  So, the bias is roughly a function of (eest-etrial)*(taunow-taueff).
! However, here we use a simpler condition.
      if (.not.l_reset_etrial .and. iblk>1 .and. iblk<5) then
        if(abs(eest-etrial).gt.3*egerr) then
          write(6,'(''Warning: abs(eest-etrial).gt.3*egerr, eest,etrial,egerr='',3es12.4,'' Fix: use better etrial'')') eest,etrial,egerr
        endif
      elseif(.not.l_reset_etrial .and. iblk>1 .and. iblk==5) then
        if(abs(eest-etrial).gt.5*egerr) then
          write(6,'(''Warning: abs(eest-etrial).gt.5*egerr, eest,etrial,egerr='',3es12.4,'' Fix: use better etrial'')') eest,etrial,egerr
          call die ('acuest_dmc_mov1_mpi2', 'Warning: abs(eest-etrial).gt.5*egerr. Fix: use better etrial')
        endif
      endif

! I have changed the dwt limit in the dmc routines so it does not depend on etrial so  there is no need for these warning msgs.
! The dwt limit is there to prevent population explosions with nonlocal psps. but there are
! better solutions than dwt limits.
!     if(wgsum(1).gt.1.5d0*nstep*nconf*nproc .or. (iblk.gt.2*nblkeq+5 .and. etrial .gt. egave+200*egerr))
!    &write(6,'(''Warning: etrial too high? It should be reasonably close to DMC energy because of dwt in dmc'')')
!     if(wgsum(1).lt.0.7d0*nstep*nconf*nproc .or. (iblk.gt.2*nblkeq+5 .and. etrial .lt. egave-200*egerr))
!    &write(6,'(''Warning: etrial too low?  It should be reasonably close to DMC energy because of dwt in dmc'')')

      call systemflush(6)

! zero out xsum variables for metrop

   17 wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r1sum=zero
      r2sum=zero
      r3sum=zero
      r4sum=zero
      risum=zero
      if(izigzag.gt.0) then
       zzsum(:)=zero
      endif

      do 20 ifr=1,nforce
        wgsum(ifr)=zero
        egsum(ifr)=zero
        pesum(ifr)=zero
        peisum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
   20   tausum(ifr)=zero

      return

      entry acues1_dmc_mov1_mpi2

! Warning: next 7 lines were attempt to save communication time
! if running in vmc mode.  But they do not work.
!     eigv=1
!     eest=etrial
!     wdsumo=nwalk
!     wgdsumo=nwalk
!     ipmod=mod(ipass,nfprod)
!     wtgen(ipmod)=nwalk
!     if(idmc.lt.0) goto 38

      call systemflush(11)
! statistical fluctuations without blocking
      wdsum1=wdsumo
      wgdsum1=wgdsumo

! Put all variable in one array, so we can do a single MPI call
      collect(1)=esum1(1)
      collect(2)=wsum1(1)
      collect(3)=efsum1
      collect(4)=wfsum1
      collect(5)=tausum(1)
      do 22 ifr=1,nforce
        collect(5+ifr)=wgsum1(ifr)
   22   collect(5+nforce+ifr)=egsum1(ifr)

! Need mpi_allreduce so that all processors can put upper bounds on weights for psp. calculations
!     call mpi_reduce(collect,collect_t,2*nforce+5,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(collect,collect_t,2*nforce+5,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

!     if(idtask.ne.0) goto 36 ! The slaves also have to calculate either wcum1,ecum1,ecm21 or wgcum1,egcum1,egcm21 to put bounds on branching in psp. calcul.

! Warning: The tausum is being done in an ugly way, and should be cleaned up
      esum1(1)=collect_t(1)
      wsum1(1)=collect_t(2)
      efsum1=collect_t(3)
      wfsum1=collect_t(4)
!     tausum(1)=collect_t(5)

      if(ipr.gt.-2 .and. idtask .eq. 0) then
         write(11,'(i8,f11.8,f15.8,f13.8)') ipass,ffn,wsum1(1),esum1(1)/wsum1(1)
      end if

      do 24 ifr=1,nforce
        wgsum1(ifr)=collect_t(5+ifr)
   24   egsum1(ifr)=collect_t(5+nforce+ifr)

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
          egcm21allprocs(ifr)=egcm21allprocs(ifr)+egsum1(ifr)**2/wgsum1(ifr)
         else
          egcm21(ifr)=0
          egcm21allprocs(ifr)=0
        endif
   30 continue
      call object_modified('wgcum1') !worry about speed
      call object_modified('wgcm21')

! collect block averages
      wsum=wsum+wsum1(1)
      wfsum=wfsum+wfsum1
      wdsum=wdsum+wdsum1
      wgdsum=wgdsum+wgdsum1
      esum=esum+esum1(1)
      efsum=efsum+efsum1
!     eisum=eisum+wfsum1/wdsum1
      do 35 ifr=1,nforce
        wgsum(ifr)=wgsum(ifr)+wgsum1(ifr)
   35   egsum(ifr)=egsum(ifr)+egsum1(ifr)

! Estimate eigenvalue of G from the energy
      ipmod=mod(ipass,nfprod)
      if(iabs(idmc).eq.1) then
        nfpro=min(nfprod,ipass)
        eigv=(wgsum1(1)/wtgen(ipmod))**(one/nfpro)
       else
        eest=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
!       eigv=dexp((etrial-eest)*(taucum(1)+taublock)/
        eigv=dexp((etrial-eest)*(taucum(1)+collect_t(5))/(wgcum(1)+wgsum(1)))
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest,egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wsum1(1)
      wgdsumo=wsum1(1)*fprod/ff(mod(ipass+1,nfprod))
      wtgen(ipmod)=wsum1(1)

   36 call mpi_bcast(eigv,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(eest,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(wdsumo,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      call redistribute

      call mpi_barrier(MPI_COMM_WORLD,ierr)

! zero out step averages
   38 wfsum1=zero
      wdsum1=zero
      efsum1=zero
!     tausum(1)=zero
      do 40 ifr=1,nforce
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        esum1(ifr)=zero
   40   egsum1(ifr)=zero

      return

      entry zeres0_dmc_mov1_mpi2
! Initialize various quantities at beginning of run
! the initial values of energy psi etc. are calculated here

      ipass=0

! set quadrature points

      if(nloc.gt.0) call rotqua

      write(6,'(''nconf,nconf_global='',9i8)') nconf,nconf_global
      nwalk=nconf
      wdsumo=nconf_global
      wgdsumo=nconf_global

      eest=0d0
      do 80 iw=1,nconf
        current_walker = iw !TA
        call object_modified_by_index (current_walker_index) !TA
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
          call hpsi(xoldw(1,1,iw,ifr),psidow(iw,ifr),psijow(iw,ifr),voldw(1,1,iw,ifr),div_vow(1,iw),d2ow(iw,ifr), &
     &    peow(iw,ifr),peiow(iw,ifr),eoldw(iw,ifr),denergy,ifr)
          if(ifr.eq.1) then
            eest=eest+eoldw(iw,ifr)   !TA
            if(ibasis.eq.3) then      ! complex calculation
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
      eest=eest/nconf !TA
#if defined(MPI)
      call MPI_Allreduce(MPI_IN_PLACE,eest,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR) !TA
      eest=eest/nproc
#endif

      eigv=dexp((etrial-eest)*tau) !TA - I do this so that the weights do not depend on etrial on the first step
      fprod=1d0 !TA
      do i=0,nfprod-1
        wtgen(i)=nconf_global
        ff(i)=eigv
        fprod=fprod*ff(i)
      enddo

      entry zerest_dmc_mov1_mpi2
! entry point to zero out all averages etc. after equilibration runs

      if (l_reset_etrial .and. ipass.ne.0) then
        dff=dexp((eest-etrial)*taucum(1)/wgcum(1)) !TA - reset etrial after equilibration runs
        do i=0,nfprod-1
          ff(i)=ff(i)*dff
          fprod=fprod*dff
        enddo
        write(6,'(/,''etrial changed from'',f11.6,'' to'',f11.6,/)') etrial, eest
        etrial=eest
!        eigv=1d0
      endif
      eigv=dexp((etrial-eest)*tau) !TA

      iblk=0

! zero out estimators

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
      r1cum=zero
      r2cum=zero
      r3cum=zero
      r4cum=zero
      ricum=zero
      dr2un=zero
      dr2ac=zero
      if(izigzag.gt.0) then
       zzcum(:)=zero
      endif

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
      r1cm2=zero
      r2cm2=zero
      r3cm2=zero
      r4cm2=zero
      ricm2=zero
      if(izigzag.gt.0) then
       zzcm2(:)=zero
      endif

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
      r1sum=zero
      r2sum=zero
      r3sum=zero
      r4sum=zero
      risum=zero
      if(izigzag.gt.0) then
       zzsum(:)=zero
      endif

      call grad_hess_jas_init

      call alloc ('fgcum', fgcum, nforce)
      call alloc ('fgcm2', fgcm2, nforce)
      call alloc ('wgcm2', wgcm2, nforce)
      call alloc ('wgcm21', wgcm21, nforce)
      call alloc ('egcm2', egcm2, nforce)
      call alloc ('egcm21', egcm21, nforce)
      call alloc ('egcm21allprocs', egcm21allprocs, nforce)
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
! Do it for MFORCE rather than nforce because in optimization at the start nforce=1 but later nforce=3
!JT: this should not be necessary and it is annoying for dynamic allocation!
      do 85 ifr=1,nforce
!JT   do 85 ifr=1,MFORCE
        tausum(ifr)=zero
        taucum(ifr)=zero
        wgcum1(ifr)=zero
        wgcum(ifr)=zero
        egcum1(ifr)=zero
        egcum(ifr)=zero
        wgcm21(ifr)=zero
        wgcm2(ifr)=zero
        egcm21(ifr)=zero
        egcm21allprocs(ifr)=zero
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

! Zero out estimators for charge density of atom.
      do 90 i=1,NRAD
        rprobup(i)=zero
        rprobdn(i)=zero
   90   rprob(i)=zero

! Zero out estimators for pair densities:
      if (ifixe.ne.0) then
      call alloc_range ('den2d_t', den2d_t, -NAX, NAX, -NAX, NAX)
      call alloc_range ('den2d_u', den2d_u, -NAX, NAX, -NAX, NAX)
      call alloc_range ('den2d_d', den2d_d, -NAX, NAX, -NAX, NAX)
      call alloc_range ('pot_ee2d_t', pot_ee2d_t, -NAX, NAX, -NAX, NAX)
      call alloc_range ('pot_ee2d_u', pot_ee2d_u, -NAX, NAX, -NAX, NAX)
      call alloc_range ('pot_ee2d_d', pot_ee2d_d, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probdt', xx0probdt, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probdu', xx0probdu, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probdd', xx0probdd, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probut', xx0probut, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probuu', xx0probuu, 0, NAX, -NAX, NAX, -NAX, NAX)
      call alloc_range ('xx0probud', xx0probud, 0, NAX, -NAX, NAX, -NAX, NAX)
      do 100 i2=-NAX,NAX
        do 100 i3=-NAX,NAX
          den2d_t(i2,i3)=0
          den2d_u(i2,i3)=0
          den2d_d(i2,i3)=0
          pot_ee2d_t(i2,i3)=0
          pot_ee2d_u(i2,i3)=0
          pot_ee2d_d(i2,i3)=0
          do 100 i1=0,NAX
            xx0probdt(i1,i2,i3)=0
            xx0probdu(i1,i2,i3)=0
            xx0probdd(i1,i2,i3)=0
            xx0probut(i1,i2,i3)=0
            xx0probuu(i1,i2,i3)=0
  100       xx0probud(i1,i2,i3)=0
      endif
      if (ifourier.ne.0) then
      call alloc_range ('fourierrk_t', fourierrk_t, -NAX, NAX, 0, NAK1)
      call alloc_range ('fourierrk_u', fourierrk_u, -NAX, NAX, 0, NAK1)
      call alloc_range ('fourierrk_d', fourierrk_d, -NAX, NAX, 0, NAK1)
      call alloc_range ('fourierkk_t', fourierkk_t, -NAK2, NAK2, -NAK2, NAK2)
      call alloc_range ('fourierkk_u', fourierkk_u, -NAK2, NAK2, -NAK2, NAK2)
      call alloc_range ('fourierkk_d', fourierkk_d, -NAK2, NAK2, -NAK2, NAK2)
      do 110 i1=-NAX,NAX
        do 110 i2=0,NAK1
          fourierrk_t(i1,i2)=0
          fourierrk_u(i1,i2)=0
  110     fourierrk_d(i1,i2)=0
      do 120 i1=-NAK2,NAK2
        do 120 i2=-NAK2,NAK2
          fourierkk_t(i1,i2)=0
          fourierkk_u(i1,i2)=0
  120     fourierkk_d(i1,i2)=0
      endif
      if (izigzag.ne.0) then
        call alloc_range ('zzpairden_t', zzpairden_t, -NAX, NAX, -NAX, NAX)
        call alloc_range ('zzpairdenij_t', zzpairdenij_t, -NAX, NAX, 0, (nelec-1))
        call alloc_range ('zzcorr', zzcorr, 0, NAX)
        call alloc_range ('zzcorr1', zzcorr1, 0, NAX)
        call alloc_range ('zzcorr2', zzcorr2, 0, NAX)
        call alloc_range ('zzcorrij', zzcorrij, 0, (nelec-1))
        call alloc_range ('zzcorrij1', zzcorrij1, 0, (nelec-1))
        call alloc_range ('zzcorrij2', zzcorrij2, 0, (nelec-1))
        call alloc_range ('yycorr', yycorr, 0, NAX)
        call alloc_range ('yycorr1', yycorr1, 0, NAX)
        call alloc_range ('yycorr2', yycorr2, 0, NAX)
        call alloc_range ('yycorrij', yycorrij, 0, (nelec-1))
        call alloc_range ('yycorrij1', yycorrij1, 0, (nelec-1))
        call alloc_range ('yycorrij2', yycorrij2, 0, (nelec-1))
        call alloc_range ('znncorr', znncorr, 0, NAX)
        call alloc_range ('zn2ncorr', zn2ncorr, 0, NAX)
        zzpairden_t(:,:) = 0
        zzpairdenij_t(:,:) = 0
        zzcorr(:) = 0
        zzcorr1(:) = 0
        zzcorr2(:) = 0
        zzcorrij(:) = 0
        zzcorrij1(:) = 0
        zzcorrij2(:) = 0
        yycorr(:) = 0
        yycorr1(:) = 0
        yycorr2(:) = 0
        yycorrij(:) = 0
        yycorrij1(:) = 0
        yycorrij2(:) = 0
        znncorr(:) = 0
        zn2ncorr(:) = 0
      endif




      call grad_hess_jas_save
# endif

      return
      end
