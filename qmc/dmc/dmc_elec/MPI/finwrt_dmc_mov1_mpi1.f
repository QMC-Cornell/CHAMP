      subroutine finwrt_dmc_mov1_mpi1
c MPI version created by Claudia Filippi starting from serial version
c routine to print out final results

# if defined (MPI)

      use main_menu_mod
      use control_mod
      use mpi_mod
      use atom_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use gradhess_mod
      use forcepar_mod
      use contrl_per_mod
      use contr3_mod
      use iterat_mod
      use mpiblk_mod
      use forcest_dmc_mod
      use denupdn_mod
      use stepv_mod
      use config_dmc_mod
      use branch_mod
      use estsum_dmc_mod
      use estcum_dmc_mod
      use contrldmc_mod
      use estcm2_mod
      use stats_mod
      use age_mod
      use pairden_mod
      use opt_ovlp_fn_mod
      implicit real*8(a-h,o-z)

c     common /force_dmc/ itausec,nwprod

c /config_dmc/ included to print out xoldw and voldw for old walkers
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4
      common /compferm/ emagv,nv,idot

      dimension eg1collect(nforce),eg21collect(nforce),wg1collect(nforce),wg21collect(nforce),rprobcollect(NRAD)
      dimension xx0probt(0:NAX,-NAX:NAX,-NAX:NAX),den2dt(-NAX:NAX,-NAX:NAX),pot_ee2dt(-NAX:NAX,-NAX:NAX)

      character*80 fmt
!JT      character*80 title,fmt
!JT      character*24 date

c statement functions for error calculation
      rn_eff(w,w2)=w**2/w2

      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errorn(x,x2,rn)=dsqrt(max((x2/rn-(x/rn)**2)/(rn-1),0.d0))
!JT      errori(x,x2,w,w2,rn)=dsqrt(max((x2/rn-(x/rn)**2)/(rn_eff(w,w2)-1),
      errc(x,x2)=error(x,x2,wcum,wcm2)
      errf(x,x2)=error(x,x2,wfcum,wfcm2)
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))
      errc1(x,x2)=error(x,x2,wcum1,wcm21)
      errf1(x,x2)=error(x,x2,wfcum1,wfcm21)
      errg1(x,x2,i)=error(x,x2,wgcum1(i),wgcm21(i))
      errw(x,x2)=errorn(x,x2,dfloat(iblk_proc))/nstep
      errw1(x,x2)=errorn(x,x2,pass_proc)
!JT      erric(x,x2)=errori(x,x2,wcum,wcm2,dfloat(iblk_proc))
!JT      erric1(x,x2)=errori(x,x2,wcum1,wcm21,pass_proc)
!JT      errig(x,x2)=errori(x,x2,wgcum(1),wgcm2(1),dfloat(iblk_proc))

      passes=dfloat(iblk)*dfloat(nstep)
      eval=nconf*passes
      pass_proc=dfloat(iblk_proc)*dfloat(nstep)
      eval_proc=nconf_global*pass_proc
c Either the next 3 lines or the 3 lines following them could be used.
c They should give nearly (but not exactly) the same result.
c Strictly the 1st 3 are for step-by-step quantities and the last 3 for blk-by-blk
c     eval_proc_eff=dfloat(nconf_global)*rn_eff(wcum1,wcm21)
c     evalf_proc_eff=dfloat(nconf_global)*rn_eff(wfcum1,wfcm21)
c     evalg_proc_eff=dfloat(nconf_global)*rn_eff(wgcum1(1),wgcm21(1))
      eval_proc_eff=dfloat(nconf_global)*nstep*rn_eff(wcum,wcm2)
      evalf_proc_eff=dfloat(nconf_global)*nstep*rn_eff(wfcum,wfcm2)
      evalg_proc_eff=dfloat(nconf_global)*nstep*rn_eff(wgcum(1),wgcm2(1))
c     rtpass1=dsqrt(passes-1)
      rtpass_proc1=dsqrt(pass_proc-1)
c     rteval=dsqrt(eval)
c     rteval_proc=dsqrt(eval_proc)
      rteval_proc_eff1=dsqrt(eval_proc_eff-1)
      rtevalf_proc_eff1=dsqrt(evalf_proc_eff-1)
      rtevalg_proc_eff1=dsqrt(evalg_proc_eff-1)

      write(6,'(/,''Final write:'')')

      write(6,'(''passes, eval, pass_proc, eval_proc, eval_proc_eff, evalf_proc_eff, evalg_proc_eff'',19f11.0)')
     & passes,eval,pass_proc,eval_proc,eval_proc_eff,evalf_proc_eff,evalg_proc_eff

c For the next 15 calls we switched from using mpi_reduce to mpi_allreduce
      call mpi_allreduce(ecum1,e1collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ecm21,e21collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wcum1,w1collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wcm21,w21collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      ecum1=e1collect
      wcum1=w1collect
      ecm21=e21collect
      wcm21=w21collect

      call mpi_allreduce(efcum1,ef1collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(efcm21,ef21collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wfcum1,wf1collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wfcm21,wf21collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

       efcum1=ef1collect
       wfcum1=wf1collect
       efcm21=ef21collect
       wfcm21=wf21collect

      call mpi_allreduce(egcum1,eg1collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(egcm21,eg21collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wgcum1,wg1collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(wgcm21,wg21collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do 1 ifr=1,nforce
        egcum1(ifr)=eg1collect(ifr)
        wgcum1(ifr)=wg1collect(ifr)
        egcm21(ifr)=eg21collect(ifr)
    1   wgcm21(ifr)=wg21collect(ifr)

c Collect radial charge density for atoms
      if(iperiodic.eq.0) then
        call mpi_reduce(rprob,rprobcollect,NRAD,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 2 i=1,NRAD
    2     rprob(i)=rprobcollect(i)
        call mpi_reduce(rprobup,rprobcollect,NRAD,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 3 i=1,NRAD
    3     rprobup(i)=rprobcollect(i)
        call mpi_reduce(rprobdn,rprobcollect,NRAD,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 4 i=1,NRAD
    4     rprobdn(i)=rprobcollect(i)
      endif

      if(ifixe.ne.0) then

        naxt=(NAX+1)*(2*NAX+1)*(2*NAX+1)
        call mpi_reduce(xx0probut,xx0probt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 5 i1=0,NAX
          do 5 i2=-NAX,NAX
            do 5 i3=-NAX,NAX
    5         xx0probut(i1,i2,i3)=xx0probt(i1,i2,i3)

        call mpi_reduce(xx0probud,xx0probt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 6 i1=0,NAX
          do 6 i2=-NAX,NAX
            do 6 i3=-NAX,NAX
    6         xx0probud(i1,i2,i3)=xx0probt(i1,i2,i3)

        call mpi_reduce(xx0probuu,xx0probt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 7 i1=0,NAX
          do 7 i2=-NAX,NAX
            do 7 i3=-NAX,NAX
    7         xx0probuu(i1,i2,i3)=xx0probt(i1,i2,i3)

        call mpi_reduce(xx0probdu,xx0probt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 8 i1=0,NAX
          do 8 i2=-NAX,NAX
            do 8 i3=-NAX,NAX
    8         xx0probdu(i1,i2,i3)=xx0probt(i1,i2,i3)

        call mpi_reduce(xx0probdd,xx0probt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 9 i1=0,NAX
          do 9 i2=-NAX,NAX
            do 9 i3=-NAX,NAX
    9         xx0probdd(i1,i2,i3)=xx0probt(i1,i2,i3)

        call mpi_reduce(xx0probdt,xx0probt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 10 i1=0,NAX
          do 10 i2=-NAX,NAX
            do 10 i3=-NAX,NAX
   10         xx0probdt(i1,i2,i3)=xx0probt(i1,i2,i3)

        naxt=(2*NAX+1)*(2*NAX+1)
        call mpi_reduce(den2d_t,den2dt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 11 i1=-NAX,NAX
          do 11 i2=-NAX,NAX
   11       den2d_t(i1,i2)=den2dt(i1,i2)

        call mpi_reduce(den2d_d,den2dt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 12 i1=-NAX,NAX
          do 12 i2=-NAX,NAX
   12       den2d_d(i1,i2)=den2dt(i1,i2)

        call mpi_reduce(den2d_u,den2dt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 13 i1=-NAX,NAX
          do 13 i2=-NAX,NAX
   13       den2d_u(i1,i2)=den2dt(i1,i2)

        call mpi_reduce(pot_ee2d_t,pot_ee2dt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 14 i1=-NAX,NAX
          do 14 i2=-NAX,NAX
   14       pot_ee2d_t(i1,i2)=pot_ee2dt(i1,i2)

        call mpi_reduce(pot_ee2d_d,pot_ee2dt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 15 i1=-NAX,NAX
          do 15 i2=-NAX,NAX
   15       pot_ee2d_d(i1,i2)=pot_ee2dt(i1,i2)

        call mpi_reduce(pot_ee2d_u,pot_ee2dt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do 16 i1=-NAX,NAX
          do 16 i2=-NAX,NAX
   16       pot_ee2d_u(i1,i2)=pot_ee2dt(i1,i2)

      endif

      call mpi_allreduce(nodecr,nodecr_collect,1,mpi_integer,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(try_int,try_int_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(acc,acc_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(acc_int,acc_int_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      nodecr=nodecr_collect
      try_int=try_int_collect
      acc=acc_collect
      acc_int=acc_int_collect

      call grad_hess_jas_mpi

c     write(11,'(4i5,f11.5,f7.4,f10.7,'' nstep,nblk,nblkeq,nconf_global,etrial,tau,taueff'')')
c    &nstep,nblk,nblkeq,nconf_global,etrial,tau,taucum(1)/wgcum(1)

      if(ipr.gt.-2) write(11,'(3i5,f11.5,f7.4,f10.7,'' nstep,nblk,nconf_global,etrial,tau,taueff'')')
     &  nstep,iblk,nconf_global,etrial,tau,taucum(1)/wgcum(1)

c Warning tmp
c     if(idtask.ne.0) return

      if(print_radial_probability .and. iperiodic.eq.0 .and. ncent.eq.1 .and. ipr.ge.-4) then
        if(ndim.eq.3) write(6,'(/,'' r  4*pi*r^2*rho 4*pi*r^2*rhoup 4*pi*r^2*rhodn'')')
        if(ndim.eq.2) write(6,'(/,'' r  2*pi*r*rho 2*pi*r*rhoup 2*pi*r*rhodn'')')
        delr=one/delri
        term=one/(wgcum(1)*delr)
        do 18 i=1,NRAD
   18     write(6,'(f5.3,3f10.6)') delr*(i-half),rprob(i)*term,rprobup(i)*term,rprobdn(i)*term
      endif

      if(idmc.ge.0) then
        write(6,'(/,''ages of walkers are:'')')
        write(6,'(10i4)') (iage(i),i=1,nwalk)
        do 19 i=1,nwalk
          if(iage(i).gt.50) then
            write(6,'(i4,i6,f10.4,99f8.4)') i,iage(i),eoldw(i,1),((xoldw(k,j,i,1),k=1,ndim),j=1,nelec)
            write(6,'(99f8.4)') ((voldw(k,j,i,1),k=1,ndim),j=1,nelec)
          endif
   19   continue

        write(6,'(''age of oldest walker (this generation, any gen)='',3i9)') ioldest,ioldestmx
      endif

      write(6,'(''average of the squares of the accepted step-size='',f10.5)') dr2ac/try_int

      write(6,'(''taueff'',20f7.4)') (taucum(ifr)/wgcum(ifr),ifr=1,nforce)

      accav=acc/try_int
      accavn=acc_int/try_int

      wave=wcum/pass_proc
      wfave=wfcum/pass_proc
      eave=ecum/wcum
      efave=efcum/wfcum
c     ei1ave=wfcum/wdcum
c     ei2ave=wgcum(1)/wgdcum
c     ei3ave=ei3cum/passes

c     r2ave=r2cum/(wgcum(1)*nelec)
c     riave=ricum/(wgcum(1)*nelec)
c     if(itau_eff.ge.1) then
c       e1ave=etrial-dlog(ei1ave)/(taucum(1)/wgcum(1))
c       e2ave=etrial-dlog(ei2ave)/(taucum(1)/wgcum(1))
c       e3ave=etrial-dlog(ei3ave)/(taucum(1)/wgcum(1))
c      else
c       if(icut_br.le.0) then
c         e1ave=etrial-dlog((ei1ave-one+accavn)/accavn)/tau
c         e2ave=etrial-dlog((ei2ave-one+accavn)/accavn)/tau
c         e3ave=etrial-dlog((ei3ave-one+accavn)/accavn)/tau
c        else
c         e1ave=etrial+dlog((accavn+2-2*ei1ave)/(accavn-2+2*ei1ave))/(4*tau)
c         e2ave=etrial+dlog((accavn+2-2*ei2ave)/(accavn-2+2*ei2ave))/(4*tau)
c         e3ave=etrial+dlog((accavn+2-2*ei3ave)/(accavn-2+2*ei3ave))/(4*tau)
c       endif
c     endif

      werr=errw(wcum,wcm2)
      wferr=errw(wfcum,wfcm2)
      werr1=errw1(wcum1,wcm21)
      wferr1=errw1(wfcum1,wfcm21)
      eerr=errc(ecum,ecm2)
      eferr=errf(efcum,efcm2)
      eerr1=errc1(ecum1,ecm21)
      eferr1=errf1(efcum1,efcm21)
c     ei1err=erri(ei1cum,ei1cm2)
c     ei2err=erri(ei2cum,ei2cm2)
c     ei3err=erri1(ei3cum,ei3cm2)
c     r2err=errg(r2cum,r2cm2,1)/nelec
c     rierr=errg(ricum,ricm2,1)/nelec
c     if(itau_eff.ge.1) then
c       e1err=dlog((ei1ave+ei1err)/(ei1ave-ei1err))/(2*taucum(1)/wgcum(1))
c       e2err=dlog((ei2ave+ei2err)/(ei2ave-ei2err))/(2*taucum(1)/wgcum(1))
c       e3err=dlog((ei3ave+ei3err)/(ei3ave-ei3err))/(2*taucum(1)/wgcum(1))
c      else
c       e1err=dlog((ei1ave+ei1err)/(ei1ave-ei1err))/(2*tau)
c       e2err=dlog((ei2ave+ei2err)/(ei2ave-ei2err))/(2*tau)
c       e3err=dlog((ei3ave+ei3err)/(ei3ave-ei3err))/(2*tau)
c     endif

c     write(6,'(''dmc_mov1_mpi '',2a10)') title
      write(fmt,'(''(/,a16,2x,a'',i3,'')'')') len_trim(title)
      write(6,fmt) mode,title
      write(6,'(''No/frac. of node crossings,acceptance='',i9,3f10.6)') nodecr,dfloat(nodecr)/try_int,accav,accavn
      if(idmc.lt.0.and.accav.lt.0.3d0) write(6,'(''Warning: Low acceptance, reduce time-step tau'')')
      if(idmc.ge.0.and.accav.lt.0.7d0) write(6,'(''Warning: Low acceptance, reduce time-step tau'')')

      if(idmc.ge.0) then
c       write(6,'(''Actual, expected # of branches for 0, inf corr time='',2f9.0)') nbrnch,nconf_global*passes*
c    &  (dlog(one+eerr1*rteval*taucum(1)/wgcum(1))/dlog(two))**2
c    &  ,nconf_global*passes*(dlog(one+eerr1*rteval*taucum(1)/wgcum(1))/dlog(two))
        write(6,'(''No. of walkers at end of run='',i5)') nwalk

        write(6,'(''nwalk_eff/nwalk         ='',2f6.3)') rn_eff(wcum1,wcm21)/pass_proc,rn_eff(wcum,wcm2)/iblk_proc
        write(6,'(''nwalk_eff/nwalk with f  ='',2f6.3)') rn_eff(wfcum1,wfcm21)/pass_proc,rn_eff(wfcum,wfcm2)/iblk_proc
        write(6,'(''nwalk_eff/nwalk with fs ='',2f6.3)')
     &  rn_eff(wgcum1(1),wgcm21(1))/pass_proc,rn_eff(wgcum(1),wgcm2(1))/iblk_proc
      endif

      write(6,'(''nconf*nproc*passes'',t20,''nconf nproc     passes nstep  nblk nblkeq
     &     tau   taueff'',/,f18.0,2i6,f11.0,3i6,2f9.5)')
     &eval*nproc,nconf,nproc,passes,nstep,iblk,nblkeq,tau,taucum(1)/wgcum(1)
c     write(6,'(''nconf_global*passes'',t19,''passes  nconf_global nstep  nblk nblkeq
c    & nproc  tau    taueff'',/,2f12.0,2i6,i7,2i5,2f9.5)')
c    &eval,passes,nconf_global,nstep,iblk,nblkeq,nproc,tau,taucum(1)/wgcum(1)
      write(6,'(''physical variable         average     rms error   sig
     &ma*T_cor  sigma   T_cor'')')
      if(idmc.ge.0) then
        write(6,'(''weights ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &  wave,werr,werr*rtpass_proc1,werr1*rtpass_proc1,(werr/werr1)**2
        write(6,'(''wts with f ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &  wfave,wferr,wferr*rtpass_proc1,wferr1*rtpass_proc1,(wferr/wferr1)**2
        do 20 ifr=1,nforce
          wgave=wgcum(ifr)/pass_proc
          wgerr=errw(wgcum(ifr),wgcm2(ifr))
          wgerr1=errw1(wgcum1(ifr),wgcm21(ifr))
          write(6,'(''wts with fs ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &    wgave,wgerr,wgerr*rtpass_proc1,wgerr1*rtpass_proc1,(wgerr/wgerr1)**2
  20    continue
        call object_provide('ovlp_trial_fn')
        write(6,'(a,f7.5)') 'approx to normalized overlap of FN and trial wave functions=',ovlp_trial_fn
        call object_provide('ovlp_trial_fn_over_ovlp_trial')
        write(6,'(a,f7.5)') 'unnormalized overlap of FN and trial wave functions', ovlp_trial_fn_over_ovlp_trial

c Mixed energy estimators
        write(6,'(''total energy (   0) ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &  eave,eerr,eerr*rteval_proc_eff1,eerr1*rteval_proc_eff1,(eerr/eerr1)**2
        write(6,'(''total energy (   1) ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &  efave,eferr,eferr*rtevalf_proc_eff1,eferr1*rtevalf_proc_eff1,(eferr/eferr1)**2
      endif
      call alloc ('eloc_tc', eloc_tc, nforce) !JT
      do 30 ifr=1,nforce
        egave=egcum(ifr)/wgcum(ifr)
        egerr=errg(egcum(ifr),egcm2(ifr),ifr)
        egerr1=errg1(egcum1(ifr),egcm21(ifr),ifr)
        eloc_tc (ifr) = (egerr/egerr1)**2 !JT
c save energy, energy_sigma and energy_err for optimization
        energy(ifr)=egave
        energy_sigma(ifr)=egerr1*rtevalg_proc_eff1
        energy_err(ifr)=egerr
        if(nforce.eq.1) then
          write(6,'(''total energy ('',i4,'') ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &    nfprod,egave,egerr,egerr*rtevalg_proc_eff1,egerr1*rtevalg_proc_eff1,(egerr/egerr1)**2
         else
          write(6,'(''total energy ('',i4,'')'',i1,''='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &    nfprod,ifr,egave,egerr,egerr*rtevalg_proc_eff1,egerr1*rtevalg_proc_eff1,(egerr/egerr1)**2
        endif
  30  continue
      call object_modified ('eloc_tc') !JT
c     write(6,'(''total energy (   0) ='',t22,f14.7,'' +-'',f11.7,f9.5)') e1ave,e1err,e1err*rteval
c     write(6,'(''total energy ('',i4,'') ='',t22,f14.7,'' +-'',f11.7,f9.5)') nfprod-1,e2ave,e2err,e2err*rteval
c     write(6,'(''total energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') e3ave,e3err,e3err*rteval
      do 40 ifr=1,nforce
        peave=pecum(ifr)/wgcum(ifr)
        peiave=peicum(ifr)/wgcum(ifr)
        tpbave=tpbcum(ifr)/wgcum(ifr)
        tjfave=tjfcum(ifr)/wgcum(ifr)

        peerr=errg(pecum(ifr),pecm2(ifr),ifr)
        peierr=errg(peicum(ifr),peicm2(ifr),ifr)
        tpberr=errg(tpbcum(ifr),tpbcm2(ifr),ifr)
        tjferr=errg(tjfcum(ifr),tjfcm2(ifr),ifr)

        if(ndim.eq.2) then
          temp=0.25d0*bext*bext/(we*we)
          tmave=(peave-peiave-emag)*temp
          tmerr=(peerr+peierr)*temp
          peave=peave-tmave-emag
c         peerr=peerr+tmerr     ! is this correct?
          peerr=peerr*(1-temp)+peierr*temp
        endif

        write(6,'(''potential energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') peave,peerr,peerr*rtevalg_proc_eff1
        write(6,'(''interaction energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') peiave,peierr,peierr*rtevalg_proc_eff1
        write(6,'(''jf kinetic energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') tjfave,tjferr,tjferr*rtevalg_proc_eff1
        write(6,'(''pb kinetic energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') tpbave,tpberr,tpberr*rtevalg_proc_eff1

        if(ndim.eq.2) then
          write(6,'(''radial mag. energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') tmave,tmerr,tmerr*rtevalg_proc_eff1
          write(6,'(''orbital mag. energy ='',t22,f14.7)') emaglz+emagv
          write(6,'(''Zeeman energy ='',t22,f14.7)') emagsz
        endif
  40  continue
      do 50 ifr=2,nforce
        fgave=egcum(ifr)/wgcum(ifr)-egcum(1)/wgcum(1)
        fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
c       fgave=fgave/deltot(ifr)
c       fgerr=fgerr/abs(deltot(ifr))
c save energy difference and error in energy difference for optimization
        force(ifr)=fgave
        force_err(ifr)=fgerr
        write(6,'(''total energy diff'',i2,t22,f14.7,'' +-'',f11.7,f9.5)') ifr,fgave,fgerr,fgerr*rtevalg_proc_eff1
  50  continue

c These are not being collected at the moment.
c     if(iperiodic.eq.0 .and. ncent.eq.1) then
c       write(6,'(''<r2>_av ='',t22,f14.7,'' +-'',f11.7,f9.5)') r2ave,r2err,r2err*rtevalg_proc_eff1
c       write(6,'(''<ri>_av ='',t22,f14.7,'' +-'',f11.7,f9.5)') riave,rierr,rierr*rtevalg_eff1
c     endif

      if(ifixe.ne.0) call den2dwrt(wgcum(1))

      call routines_write_final
      call reinit_routines_write_block
      call reinit_routines_write_final
# endif

      return
      end
