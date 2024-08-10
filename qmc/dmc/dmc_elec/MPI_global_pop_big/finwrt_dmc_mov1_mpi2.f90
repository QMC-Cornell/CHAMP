      subroutine finwrt_dmc_mov1_mpi2
! MPI version created by Claudia Filippi starting from serial version
! routine to print out final results

# if defined (MPI)

!     use all_tools_mod
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
      use fourier_mod
      use opt_ovlp_fn_mod
      use zigzag_mod
      implicit real*8(a-h,o-z)

! /config_dmc/ included to print out xoldw and voldw for old walkers
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot

      dimension rprobcollect(NRAD)
      dimension xx0probt(0:NAX,-NAX:NAX,-NAX:NAX),den2dt(-NAX:NAX,-NAX:NAX),pot_ee2dt(-NAX:NAX,-NAX:NAX)
      dimension fouriert(-NAX:NAX,0:NAK1), fourierkkt(-NAK2:NAK2,-NAK2:NAK2)
      dimension zzpairtot(-NAX:NAX,-NAX:NAX),zzdenijtot(-NAX:NAX,0:(nelec-1))
      dimension zzcorrtot(0:NAX), zzcorrijtot(0:(nelec-1))
      dimension zzave(nzzvars), zzerr(nzzvars)

      character*80 fmt
!JT      character*80 title,fmt
!JT      character*24 date

! statement functions for error calculation
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
      errw(x,x2)=errorn(x,x2,dfloat(iblk))/nstep
      errw1(x,x2)=errorn(x,x2,passes)
!JT      erric(x,x2)=errori(x,x2,wcum,wcm2,dfloat(iblk)) ! unused
!JT      erric1(x,x2)=errori(x,x2,wcum1,wcm21,passes)  ! unused
!JT      errig(x,x2)=errori(x,x2,wgcum(1),wgcm2(1),dfloat(iblk))  ! unused

      passes=dfloat(iblk)*dfloat(nstep)
      eval=nconf_global*passes
! Either the next 3 lines or the 3 lines following them could be used.
! They should give nearly (but not exactly) the same result.
! Strictly the 1st 3 are for step-by-step quantities and the last 3 for blk-by-blk
!     eval_eff=dfloat(nconf_global)*rn_eff(wcum1,wcm21)
!     evalf_eff=dfloat(nconf_global)*rn_eff(wfcum1,wfcm21)
!     evalg_eff=dfloat(nconf_global)*rn_eff(wgcum1(1),wgcm21(1))
      eval_eff=dfloat(nconf_global)*nstep*rn_eff(wcum,wcm2)
      evalf_eff=dfloat(nconf_global)*nstep*rn_eff(wfcum,wfcm2)
      write(6,'(''evalf_eff,nconf_global,nstep,wfcum,wfcm2'',d12.4,2i5,9d12.4)') evalf_eff,nconf_global,nstep,wfcum,wfcm2
      evalg_eff=nconf_global*nstep*rn_eff(wgcum(1),wgcm2(1))
      rtpass1=dsqrt(passes-1)
      rteval=dsqrt(eval)
      rteval_eff1=dsqrt(eval_eff-1)
      rtevalf_eff1=dsqrt(evalf_eff-1)
      rtevalg_eff1=dsqrt(evalg_eff-1)

      write(6,'(/,''Final write:'')')

!     write(6,'(''wfcum,wfcm2='',9d12.4)') wfcum,wfcm2

      write(6,'(''passes, eval, eval_eff, evalf_eff, evalg_eff'',19f11.0)') passes,eval,eval_eff,evalf_eff,evalg_eff

! Collect radial charge density for atoms
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

      if(ifourier .ne. 0) then
        naxt = (2*NAX + 1) * (NAK1 + 1)
        call mpi_reduce(fourierrk_t,fouriert,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do i1=-NAX,NAX
          do i2=0,NAK1
            fourierrk_t(i1,i2)=fouriert(i1,i2)
          enddo
        enddo

        call mpi_reduce(fourierrk_u,fouriert,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do i1=-NAX,NAX
          do i2=0,NAK1
            fourierrk_u(i1,i2)=fouriert(i1,i2)
          enddo
        enddo

        call mpi_reduce(fourierrk_d,fouriert,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do i1=-NAX,NAX
          do i2=0,NAK1
            fourierrk_d(i1,i2)=fouriert(i1,i2)
          enddo
        enddo

        naxt = (2*NAK2 + 1) * (2*NAK2 + 1)

        call mpi_reduce(fourierkk_t,fourierkkt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do i1=-NAK2,NAK2
          do i2=-NAK2,NAK2
            fourierkk_t(i1,i2)=fourierkkt(i1,i2)
          enddo
        enddo

        call mpi_reduce(fourierkk_u,fourierkkt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do i1=-NAK2,NAK2
          do i2=-NAK2,NAK2
            fourierkk_u(i1,i2)=fourierkkt(i1,i2)
          enddo
        enddo

        call mpi_reduce(fourierkk_d,fourierkkt,naxt,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        do i1=-NAK2,NAK2
          do i2=-NAK2,NAK2
            fourierkk_d(i1,i2)=fourierkkt(i1,i2)
          enddo
        enddo

      endif

      if(izigzag.gt.0) then
        call mpi_allreduce(zzcorr, zzcorrtot, NAX+1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        zzcorr(:) = zzcorrtot(:)
        call mpi_allreduce(zzcorr1, zzcorrtot, NAX+1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        zzcorr1(:) = zzcorrtot(:)
        call mpi_allreduce(zzcorr2, zzcorrtot, NAX+1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        zzcorr2(:) = zzcorrtot(:)
        call mpi_allreduce(yycorr, zzcorrtot, NAX+1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        yycorr(:) = zzcorrtot(:)
        call mpi_allreduce(yycorr1, zzcorrtot, NAX+1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        yycorr1(:) = zzcorrtot(:)
        call mpi_allreduce(yycorr2, zzcorrtot, NAX+1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        yycorr2(:) = zzcorrtot(:)
        call mpi_allreduce(znncorr, zzcorrtot, NAX+1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        znncorr(:) = zzcorrtot(:)
        call mpi_allreduce(zn2ncorr, zzcorrtot, NAX+1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        zn2ncorr(:) = zzcorrtot(:)
        call mpi_allreduce(zzcorrij, zzcorrijtot, nelec,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        zzcorrij(:) = zzcorrijtot(:)
        call mpi_allreduce(zzcorrij1, zzcorrijtot, nelec,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        zzcorrij1(:) = zzcorrijtot(:)
        call mpi_allreduce(zzcorrij2, zzcorrijtot, nelec,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        zzcorrij2(:) = zzcorrijtot(:)
        call mpi_allreduce(yycorrij, zzcorrijtot, nelec,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        yycorrij(:) = zzcorrijtot(:)
        call mpi_allreduce(yycorrij1, zzcorrijtot, nelec,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        yycorrij1(:) = zzcorrijtot(:)
        call mpi_allreduce(yycorrij2, zzcorrijtot, nelec,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        yycorrij2(:) = zzcorrijtot(:)

        if(izigzag.eq.2) then
          naxt = (2*NAX + 1) * (2*NAX + 1)
          call mpi_allreduce(zzpairden_t, zzpairtot, naxt,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
          zzpairden_t(:,:) = zzpairtot(:,:)

          naxt = (2*NAX + 1) * nelec
          call mpi_allreduce(zzpairdenij_t, zzdenijtot, naxt,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
          zzpairdenij_t(:,:) = zzdenijtot(:,:)
        endif
      endif

      call mpi_allreduce(nodecr,nodecr_collect,1,mpi_integer,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(try_int,try_int_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(acc,acc_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(acc_int,acc_int_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(dr2ac,dr2ac_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(dr2un,dr2un_collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      nodecr=nodecr_collect
      try_int=try_int_collect
      acc=acc_collect
      acc_int=acc_int_collect
      dr2ac=dr2ac_collect
      dr2un=dr2un_collect

      call grad_hess_jas_mpi

!     write(11,'(4i5,f11.5,f7.4,f10.7,'' nstep,nblk,nblkeq,nconf_global,etrial,tau,taueff'')')
!    &nstep,nblk,nblkeq,nconf_global,etrial,tau,taucum(1)/wgcum(1)

      if(ipr.gt.-2 .and. idtask .eq. 0) write(11,'(3i5,f11.5,f7.4,f10.7,'' nstep,nblk,nconf_global,etrial,tau,taueff'')') &
     &  nstep,iblk,nconf_global,etrial,tau,taucum(1)/wgcum(1)

! Warning tmp
!     if(idtask.ne.0) return

      if(print_radial_probability .and. iperiodic.eq.0 .and. ncent.eq.1 .and. ipr.ge.-4) then
        if(ndim.eq.3) write(6,'(/,'' r  4*pi*r^2*rho 4*pi*r^2*rhoup 4*pi*r^2*rhodn'')')
        if(ndim.eq.2) write(6,'(/,'' r  2*pi*r*rho 2*pi*r*rhoup 2*pi*r*rhodn'')')
        delr=one/delri
        term=one/(wgcum(1)*delr)
        do 18 i=1,NRAD
   18     write(6,'(f5.3,3f10.6)') delr*(i-half),rprob(i)*term,rprobup(i)*term,rprobdn(i)*term
      endif
      write(6,'(''5'')')

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

      write(6,'(a,f10.5)') 'average of the squares drift-dif moves for accepted steps = ',dr2ac/try_int
      write(6,'(a,f10.5)') 'average of the squares drift-dif moves for all      steps = ',dr2un/try_int

      write(6,'(''taueff'',20f7.4)') (taucum(ifr)/wgcum(ifr),ifr=1,nforce)

      accav=acc/try_int
      accavn=acc_int/try_int

      wave=wcum/passes
      wfave=wfcum/passes
      eave=ecum/wcum
      efave=efcum/wfcum
!     ei1ave=wfcum/wdcum
!     ei2ave=wgcum(1)/wgdcum
!     ei3ave=ei3cum/passes

      r1ave=r1cum/(wgcum(1)*nelec)
      r2ave=r2cum/(wgcum(1)*nelec)
      r3ave=r3cum/(wgcum(1)*nelec)
      r4ave=r4cum/(wgcum(1)*nelec)
      riave=ricum/(wgcum(1)*nelec)
      if(izigzag.gt.0) then
       zzave(:)=zzcum(:)/wgcum(1)
      endif
!     if(itau_eff.ge.1) then
!       e1ave=etrial-dlog(ei1ave)/(taucum(1)/wgcum(1))
!       e2ave=etrial-dlog(ei2ave)/(taucum(1)/wgcum(1))
!       e3ave=etrial-dlog(ei3ave)/(taucum(1)/wgcum(1))
!      else
!       if(icut_br.le.0) then
!         e1ave=etrial-dlog((ei1ave-one+accavn)/accavn)/tau
!         e2ave=etrial-dlog((ei2ave-one+accavn)/accavn)/tau
!         e3ave=etrial-dlog((ei3ave-one+accavn)/accavn)/tau
!        else
!         e1ave=etrial+dlog((accavn+2-2*ei1ave)/(accavn-2+2*ei1ave))/(4*tau)
!         e2ave=etrial+dlog((accavn+2-2*ei2ave)/(accavn-2+2*ei2ave))/(4*tau)
!         e3ave=etrial+dlog((accavn+2-2*ei3ave)/(accavn-2+2*ei3ave))/(4*tau)
!       endif
!     endif

      werr=errw(wcum,wcm2)
      wferr=errw(wfcum,wfcm2)
      werr1=errw1(wcum1,wcm21)
      wferr1=errw1(wfcum1,wfcm21)
      eerr=errc(ecum,ecm2)
      eferr=errf(efcum,efcm2)
      eerr1=errc1(ecum1,ecm21)
      eferr1=errf1(efcum1,efcm21)
!     ei1err=erri(ei1cum,ei1cm2)
!     ei2err=erri(ei2cum,ei2cm2)
!     ei3err=erri1(ei3cum,ei3cm2)
      r1err=errg(r1cum,r1cm2,1)/nelec
      r2err=errg(r2cum,r2cm2,1)/nelec
      r3err=errg(r3cum,r3cm2,1)/nelec
      r4err=errg(r4cum,r4cm2,1)/nelec
      rierr=errg(ricum,ricm2,1)/nelec
      if(izigzag.gt.0) then
       do iz=1,nzzvars
        zzerr(iz)=errg(zzcum(iz),zzcm2(iz),1)
       enddo
      endif
!     if(itau_eff.ge.1) then
!       e1err=dlog((ei1ave+ei1err)/(ei1ave-ei1err))/(2*taucum(1)/wgcum(1))
!       e2err=dlog((ei2ave+ei2err)/(ei2ave-ei2err))/(2*taucum(1)/wgcum(1))
!       e3err=dlog((ei3ave+ei3err)/(ei3ave-ei3err))/(2*taucum(1)/wgcum(1))
!      else
!       e1err=dlog((ei1ave+ei1err)/(ei1ave-ei1err))/(2*tau)
!       e2err=dlog((ei2ave+ei2err)/(ei2ave-ei2err))/(2*tau)
!       e3err=dlog((ei3ave+ei3err)/(ei3ave-ei3err))/(2*tau)
!     endif

!     write(6,'(''dmc_mov1_mpi_globalpop '',2a10)') title,date
      write(fmt,'(''(/,a16,2x,a'',i3,'')'')') len_trim(title)
      write(6,fmt) mode,title
      write(6,'(''No/frac. of node crossings,acceptance='',i9,3f10.6)')nodecr,dfloat(nodecr)/try_int,accav,accavn
      if(idmc.lt.0.and.accav.lt.0.3d0) write(6,'(''Warning: Low acceptance, reduce time-step tau'')')
      if(idmc.ge.0.and.accav.lt.0.7d0) write(6,'(''Warning: Low acceptance, reduce time-step tau'')')

      if(idmc.ge.0) then
!       write(6,'(''Actual, expected # of branches for 0, inf corr time='',i6,2f9.0)') nbrnch,nconf_global*passes*
!    &  (dlog(one+eerr1*rteval*taucum(1)/wgcum(1))/dlog(two))**2,nconf_global*passes*
!    &  (dlog(one+eerr1*rteval*taucum(1)/wgcum(1))/dlog(two))

        write(6,'(''No. of walkers at end of run='',i5)') nwalk

        write(6,'(''nwalk_eff/nwalk         ='',2f6.3)') rn_eff(wcum1,wcm21)/passes,rn_eff(wcum,wcm2)/iblk
        write(6,'(''nwalk_eff/nwalk with f  ='',2f6.3)') rn_eff(wfcum1,wfcm21)/passes,rn_eff(wfcum,wfcm2)/iblk
        write(6,'(''nwalk_eff/nwalk with fs ='',2f6.3)') rn_eff(wgcum1(1),wgcm21(1))/passes,rn_eff(wgcum(1),wgcm2(1))/iblk
      endif

      write(6,'(''nconf*nproc*passes'',t20,''nconf nproc     passes nstep  nblk nblkeq &
          &tau   taueff'',/,f18.0,2i6,f11.0,3i6,2f9.5)') &
     &eval*nproc,nconf,nproc,passes,nstep,iblk,nblkeq,tau,taucum(1)/wgcum(1)
!     write(6,'(''nconf_global*passes'',t19,''passes  nconf_global nstep  nblk nblkeq
!    & nproc  tau    taueff'',/,2f12.0,2i6,i7,2i5,2f9.5)')
!    &eval,passes,nconf_global,nstep,iblk,nblkeq,nproc,tau,taucum(1)/wgcum(1)
      write(6,'(''physical variable         average     rms error   sigma*T_cor  sigma   T_cor'')')
      if(idmc.ge.0) then
        write(6,'(''weights ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') wave,werr,werr*rtpass1,werr1*rtpass1,(werr/werr1)**2
        write(6,'(''wts with f ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') wfave,wferr,wferr*rtpass1,wferr1*rtpass1,(wferr/wferr1)**2
        do 20 ifr=1,nforce
          wgave=wgcum(ifr)/passes
          wgerr=errw(wgcum(ifr),wgcm2(ifr))
          wgerr1=errw1(wgcum1(ifr),wgcm21(ifr))
          write(6,'(''wts with fs ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
     &    wgave,wgerr,wgerr*rtpass1,wgerr1*rtpass1,(wgerr/wgerr1)**2
  20    continue
        call object_provide('ovlp_trial_fn')
        write(6,'(a,f10.8)') 'approx. normalized overlap of FN and trial wave functions= ',ovlp_trial_fn
        call object_provide('ovlp_trial_fn_over_ovlp_trial')
        write(6,'(a,f10.8)') 'unnormalized overlap of FN and trial wave functions= ', ovlp_trial_fn_over_ovlp_trial

! Mixed energy estimators
        write(6,'(''total energy (   0) ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
     &  eave,eerr,eerr*rteval_eff1,eerr1*rteval_eff1,(eerr/eerr1)**2
        write(6,'(''total energy (   1) ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
     &  efave,eferr,eferr*rtevalf_eff1,eferr1*rtevalf_eff1,(eferr/eferr1)**2
      endif
      call alloc ('eloc_tc', eloc_tc, nforce) !JT
      do 30 ifr=1,nforce
        egave=egcum(ifr)/wgcum(ifr)
        egerr=errg(egcum(ifr),egcm2(ifr),ifr)
        egerr1=errg1(egcum1(ifr),egcm21(ifr),ifr)
        eloc_tc (ifr) = (egerr/egerr1)**2 !JT
! save energy, energy_sigma and energy_err for optimization
        energy(ifr)=egave
        energy_sigma(ifr)=egerr1*rtevalg_eff1
        energy_err(ifr)=egerr
        if(nforce.eq.1) then
          write(6,'(''total energy ('',i4,'') ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
     &    nfprod,egave,egerr,egerr*rtevalg_eff1,egerr1*rtevalg_eff1,(egerr/egerr1)**2
         else
          write(6,'(''total energy ('',i4,'')'',i1,''='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
     &    nfprod,ifr,egave,egerr,egerr*rtevalg_eff1,egerr1*rtevalg_eff1,(egerr/egerr1)**2
        endif
  30  continue
      call object_modified ('eloc_tc') !JT
!     write(6,'(''total energy (   0) ='',t22,f14.7,'' +-'',f11.7,f10.5)') e1ave,e1err,e1err*rteval
!     write(6,'(''total energy ('',i4,'') ='',t22,f14.7,'' +-'',f11.7,f10.5)') nfprod-1,e2ave,e2err,e2err*rteval
!     write(6,'(''total energy ='',t22,f14.7,'' +-'',f11.7,f10.5)') e3ave,e3err,e3err*rteval
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
!         peerr=peerr+tmerr     ! is this correct?
          peerr=peerr*(1-temp)+peierr*temp
        endif

        write(6,'(''potential energy ='',t22,f14.7,'' +-'',f11.7,f10.5)') peave,peerr,peerr*rtevalg_eff1
        write(6,'(''interaction energy ='',t22,f14.7,'' +-'',f11.7,f10.5)') peiave,peierr,peierr*rtevalg_eff1
        write(6,'(''jf kinetic energy ='',t22,f14.7,'' +-'',f11.7,f10.5)') tjfave,tjferr,tjferr*rtevalg_eff1
        write(6,'(''pb kinetic energy ='',t22,f14.7,'' +-'',f11.7,f10.5)') tpbave,tpberr,tpberr*rtevalg_eff1

        if(ndim.eq.2) then
          write(6,'(''radial mag. energy ='',t22,f14.7,'' +-'',f11.7,f10.5)') tmave,tmerr,tmerr*rtevalg_proc_eff1
          write(6,'(''orbital mag. energy ='',t22,f14.7)') emaglz+emagv
          write(6,'(''Zeeman energy ='',t22,f14.7)') emagsz
        endif
  40  continue
      do 50 ifr=2,nforce
        fgave=egcum(ifr)/wgcum(ifr)-egcum(1)/wgcum(1)
        fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
!       fgave=fgave/deltot(ifr)
!       fgerr=fgerr/abs(deltot(ifr))
! save energy difference and error in energy difference for optimization
        force(ifr)=fgave
        force_err(ifr)=fgerr
        write(6,'(''total energy diff'',i2,t22,f14.7,'' +-'',f11.7,f10.5)') ifr,fgave,fgerr,fgerr*rtevalg_eff1
  50  continue

      if(iperiodic.eq.0 .and. ncent.eq.1) then
        write(6,'(''<r>_av ='',t22,f14.7,'' +-'',f11.7,f10.5)') r1ave,r1err,r1err*rtevalg_eff1
        write(6,'(''<r2>_av ='',t22,f14.7,'' +-'',f11.7,f10.5)') r2ave,r2err,r2err*rtevalg_eff1
        write(6,'(''<r3>_av ='',t22,f14.6,'' +-'',f11.6,f10.4)') r3ave,r3err,r3err*rtevalg_eff1
        write(6,'(''<r4>_av ='',t22,f14.5,'' +-'',f11.5,f10.3)') r4ave,r4err,r4err*rtevalg_eff1
        write(6,'(''<ri>_av ='',t22,f14.7,'' +-'',f11.7,f10.5)') riave,rierr,rierr*rtevalg_eff1
      endif

      if(izigzag.ge.1) then
        call print_zigzag_vars(zzave,zzerr,rtevalg_eff1)
!       write(6,'(''<ZigZag Amp> ='',t17,f12.7,'' +-'',f11.7,f10.5)') zzave(3),zzerr(3),zzerr(3)*rtevalg_eff1
!       write(6,'(''<|ZigZag Amp|> ='',t17,f12.7,'' +-'',f11.7,f10.5)') zzave(1),zzerr(1),zzerr(1)*rtevalg_eff1
!       write(6,'(''<ZigZag Amp^2> ='',t17,f12.7,'' +-'',f11.7,f10.5)') zzave(2),zzerr(2),zzerr(2)*rtevalg_eff1
!       write(6,'(''<ZigZag Amp (red)>='',t22,f12.7,'' +-'',f11.7,f10.5)') zzave(6),zzerr(6),zzerr(6)*rtevalg_eff1
!       write(6,'(''<|ZigZag Amp| (red)>='',t22,f12.7,'' +-'',f11.7,f10.5)') zzave(4),zzerr(4),zzerr(4)*rtevalg_eff1
!       write(6,'(''<ZigZag Amp^2 (red)>='',t22,f12.7,'' +-'',f11.7,f10.5)') zzave(5),zzerr(5),zzerr(5)*rtevalg_eff1
!       write(6,'(''<ZigZag rand Amp>='',t22,f12.7,'' +-'',f11.7,f10.5)') zzave(9),zzerr(9),zzerr(9)*rtevalg_eff1
!       write(6,'(''<|ZigZag rand Amp|>='',t22,f12.7,'' +-'',f11.7,f10.5)') zzave(7),zzerr(7),zzerr(7)*rtevalg_eff1
!       write(6,'(''<ZigZag rand Amp^2>='',t22,f12.7,'' +-'',f11.7,f10.5)') zzave(8),zzerr(8),zzerr(8)*rtevalg_eff1
      endif

      if(ifixe.ne.0 .or. ifourier.ne.0 .or. izigzag.ne.0) call den2dwrt(wgcum(1),r1ave)

      call routines_write_final
      call reinit_routines_write_block
      call reinit_routines_write_final
# endif

      return
      end
