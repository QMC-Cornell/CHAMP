      subroutine finwrt_dmc_mov1 !JT
c Written by Cyrus Umrigar, modified by Claudia Filippi
c routine to print out final results

      use all_tools_mod
      use constants_mod
      use control_mod
      use main_menu_mod
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

c     common /force_dmc/ itausec,nwprod

c /config_dmc/ included to print out xoldw and voldw for old walkers
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot

      character*80 fmt
  
      dimension zzave(nzzvars), zzerr(nzzvars)
c statement functions for error calculation
      rn_eff(w,w2)=w**2/w2

      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errorn(x,x2,rn)=dsqrt(max((x2/rn-(x/rn)**2)/(rn-1),0.d0))
      errori(x,x2,w,w2,rn)=dsqrt(max((x2/rn-(x/rn)**2)/(rn_eff(w,w2)-1),0.d0))
      errc(x,x2)=error(x,x2,wcum,wcm2)
      errf(x,x2)=error(x,x2,wfcum,wfcm2)
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))
      errc1(x,x2)=error(x,x2,wcum1,wcm21)
      errf1(x,x2)=error(x,x2,wfcum1,wfcm21)
      errg1(x,x2,i)=error(x,x2,wgcum1(i),wgcm21(i))
      errw(x,x2)=errorn(x,x2,dfloat(iblk))/nstep
      errw1(x,x2)=errorn(x,x2,passes)
      erric(x,x2)=errori(x,x2,wcum,wcm2,dfloat(iblk))
      erric1(x,x2)=errori(x,x2,wcum1,wcm21,passes)
      errig(x,x2)=errori(x,x2,wgcum(1),wgcm2(1),dfloat(iblk))

      passes=dfloat(iblk)*dfloat(nstep)
      eval=nconf_global*passes
c Either the next 3 lines or the 3 lines following them could be used.
c They should give nearly (but not exactly) the same result.
c Strictly the 1st 3 are for step-by-step quantities and the last 3 for blk-by-blk
c     eval_eff=dfloat(nconf_global)*rn_eff(wcum1,wcm21)
c     evalf_eff=dfloat(nconf_global)*rn_eff(wfcum1,wfcm21)
c     evalg_eff=dfloat(nconf_global)*rn_eff(wgcum1(1),wgcm21(1))
      eval_eff=dfloat(nconf_global)*nstep*rn_eff(wcum,wcm2)
      evalf_eff=dfloat(nconf_global)*nstep*rn_eff(wfcum,wfcm2)
      evalg_eff=dfloat(nconf_global)*nstep*rn_eff(wgcum(1),wgcm2(1))
      rtpass1=dsqrt(passes-1)
      rteval=dsqrt(eval)
      rteval_eff1=dsqrt(eval_eff-1)
      rtevalf_eff1=dsqrt(evalf_eff-1)
      rtevalg_eff1=dsqrt(evalg_eff-1)

      write(6,'(/,''Final write:'')')

c Write out radial charge density for atoms
      if(print_radial_probability .and. iperiodic.eq.0 .and. ncent.eq.1 .and. ipr.ge.-4) then
        if(ndim.eq.3) write(6,'(/,'' r  4*pi*r^2*rho 4*pi*r^2*rhoup 4*pi*r^2*rhodn'')')
        if(ndim.eq.2) write(6,'(/,'' r  2*pi*r*rho 2*pi*r*rhoup 2*pi*r*rhodn'')')
        delr=one/delri
        term=one/(wgcum(1)*delr)
        do 5 i=1,NRAD
    5     write(6,'(f12.8,3f10.6)') delr*(i-half),rprob(i)*term,rprobup(i)*term,rprobdn(i)*term
      endif

      if(idmc.ge.0) then
        write(6,'(/,''ages of walkers are:'')')
        write(6,'(10i4)') (iage(i),i=1,nwalk)
        do 10 i=1,nwalk
          if(iage(i).gt.50) then
            write(6,'(i4,i6,f10.4,99f8.4)') i,iage(i),eoldw(i,1),((xoldw(k,j,i,1),k=1,ndim),j=1,nelec)
            write(6,'(99f8.4)') ((voldw(k,j,i,1),k=1,ndim),j=1,nelec)
          endif
   10   continue

        write(6,'(''age of oldest walker (this generation, any gen)='',3i9)') ioldest,ioldestmx
      endif

      write(6,'(''average of the squares of the accepted step-size='',f10.5)') dr2ac/try_int

      write(6,'(''taueff'',20f7.4)') (taucum(ifr)/wgcum(ifr),ifr=1,nforce)

      accav=acc/try_int
      accavn=acc_int/try_int

      wave=wcum/passes
      wfave=wfcum/passes
      eave=ecum/wcum
      efave=efcum/wfcum
      ei1ave=wfcum/wdcum
      ei2ave=wgcum(1)/wgdcum
      ei3ave=ei3cum/passes

      r1ave=r1cum/(wgcum(1)*nelec)
      r2ave=r2cum/(wgcum(1)*nelec)
      riave=ricum/(wgcum(1)*nelec)
      if(izigzag.gt.0) then
       zzave(:)=zzcum(:)/wgcum(1)
      endif
      e1ave=etrial-dlog(ei1ave)/(taucum(1)/wgcum(1))
      e2ave=etrial-dlog(ei2ave)/(taucum(1)/wgcum(1))
      e3ave=etrial-dlog(ei3ave)/(taucum(1)/wgcum(1))

      if(iblk.eq.1) then
        werr=0
        wferr=0
        werr1=0
        wferr1=0
        eerr=0
        eferr=0
        eerr1=0
        eferr1=0
        ei1err=0
        ei2err=0
        ei3err=0
        r1err=0
        r2err=0
        rierr=0 
        if(izigzag.gt.0) then
         zzerr(:)=0.d0
        endif
       else
        werr=errw(wcum,wcm2)
        wferr=errw(wfcum,wfcm2)
        werr1=errw1(wcum1,wcm21)
        wferr1=errw1(wfcum1,wfcm21)
        eerr=errc(ecum,ecm2)
        eferr=errf(efcum,efcm2)
        eerr1=errc1(ecum1,ecm21)
        eferr1=errf1(efcum1,efcm21)
        ei1err=erric(ei1cum,ei1cm2)
        ei2err=errig(ei2cum,ei2cm2)
        ei3err=erric1(ei3cum,ei3cm2)
        r1err=errg(r1cum,r1cm2,1)/nelec
        r2err=errg(r2cum,r2cm2,1)/nelec
        rierr=errg(ricum,ricm2,1)/nelec
        if(izigzag.gt.0) then
         do iz=1,nzzvars
          zzerr(iz)=errg(zzcum(iz),zzcm2(iz),1)
         enddo
        endif
      endif

      e1err=dlog((ei1ave+ei1err)/(ei1ave-ei1err))/(2*taucum(1)/wgcum(1))
      e2err=dlog((ei2ave+ei2err)/(ei2ave-ei2err))/(2*taucum(1)/wgcum(1))
      e3err=dlog((ei3ave+ei3err)/(ei3ave-ei3err))/(2*taucum(1)/wgcum(1))

c     write(6,'(''dmc_mov1 '',2a10)') title
      write(fmt,'(''(/,a16,2x,a'',i3,'')'')') len_trim(title)
      write(6,fmt) mode,title

      write(6,'(''No/frac. of node crossings,acceptance='',i9,3f10.6)') nodecr,dfloat(nodecr)/try_int,accav,accavn
      if(idmc.lt.0.and.accav.lt.0.3) write(6,'(''Warning: Low acceptance, reduce time-step tau'')')
      if(idmc.ge.0.and.accav.lt.0.7) write(6,'(''Warning: Low acceptance, reduce time-step tau'')')
      if(idmc.ge.0) then
        write(6,'(''Actual, expected # of branches for 0, inf corr time='',i6,2f9.0)') nbrnch
     &  ,nconf_global*passes*(dlog(one+eerr1*rteval*taucum(1)/wgcum(1))/dlog(two))**2
     &  ,nconf_global*passes*(dlog(one+eerr1*rteval*taucum(1)/wgcum(1))/dlog(two))
        write(6,'(''No. of walkers at end of run='',i5)') nwalk

        write(6,'(''nwalk_eff/nwalk         ='',2f6.3)') rn_eff(wcum1,wcm21)/passes,rn_eff(wcum,wcm2)/iblk
        write(6,'(''nwalk_eff/nwalk with f  ='',2f6.3)') rn_eff(wfcum1,wfcm21)/passes,rn_eff(wfcum,wfcm2)/iblk
        write(6,'(''nwalk_eff/nwalk with fs ='',2f6.3)') rn_eff(wgcum1(1),wgcm21(1))/passes,rn_eff(wgcum(1),wgcm2(1))/iblk
      endif

      write(6,'(''nconf*passes     passes   nconf nstep  nblk nblkeq  tau    taueff'',/,
     & 2f12.0,2i6,i7,i5,2f9.5)') eval,passes,nconf_global,nstep,iblk,nblkeq,tau,taucum(1)/wgcum(1)
      write(6,'(''physical variable         average     rms error   sigma*T_cor  sigma   T_cor'')')
      if(idmc.ge.0) then
        write(6,'(''weights ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)') wave,werr,werr*rtpass1,werr1*rtpass1,(werr/werr1)**2
        write(6,'(''wts with f ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)') wfave,wferr,wferr*rtpass1,wferr1*rtpass1,(wferr/wferr1)**2
        do 20 ifr=1,nforce
          wgave=wgcum(ifr)/passes
          wgerr=errw(wgcum(ifr),wgcm2(ifr))
          wgerr1=errw1(wgcum1(ifr),wgcm21(ifr))
          write(6,'(''wts with fs ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &    wgave,wgerr,wgerr*rtpass1,wgerr1*rtpass1,(wgerr/wgerr1)**2
  20    continue
        call object_provide('ovlp_trial_fn')
        write(6,'(''Appr overlap psi_t psi_fn='',f10.8)') ovlp_trial_fn
        write(6,'(a,f7.5)') 'approx to normalized overlap of FN and trial wave functions=',ovlp_trial_fn
        call object_provide('ovlp_trial_fn_over_ovlp_trial')
        write(6,'(a,f7.5)') 'unnormalized overlap of FN and trial wave functions', ovlp_trial_fn_over_ovlp_trial

c Mixed energy estimators
        write(6,'(''total energy (   0) ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &  eave,eerr,eerr*rteval_eff1,eerr1*rteval_eff1,(eerr/eerr1)**2
        write(6,'(''total energy (   1) ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &  efave,eferr,eferr*rtevalf_eff1,eferr1*rtevalf_eff1,(eferr/eferr1)**2
      endif
      call alloc ('eloc_tc', eloc_tc, nforce) !JT
      do 30 ifr=1,nforce
        egave=egcum(ifr)/wgcum(ifr)
        egerr=errg(egcum(ifr),egcm2(ifr),ifr)
        egerr1=errg1(egcum1(ifr),egcm21(ifr),ifr)
        eloc_tc (ifr) = (egerr/egerr1)**2 !JT
c save energy, energy_sigma and energy_err for optimization
        energy(ifr)=egave
        energy_sigma(ifr)=egerr1*rtevalg_eff1
        energy_err(ifr)=egerr
        if(nforce.eq.1) then
          write(6,'(''total energy ('',i4,'') ='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &    nfprod,egave,egerr,egerr*rtevalg_eff1,egerr1*rtevalg_eff1,(egerr/egerr1)**2
         else
          write(6,'(''total energy ('',i4,'')'',i1,''='',t22,f14.7,'' +-'',f11.7,2f9.5,f8.2)')
     &    nfprod,ifr,egave,egerr,egerr*rtevalg_eff1,egerr1*rtevalg_eff1,(egerr/egerr1)**2
        endif
  30  continue
      call object_modified ('eloc_tc') !JT
c Growth energy estimators
      if(idmc.ge.0) then
        write(6,'(''total energy (   0) ='',t22,f14.7,'' +-'',f11.7,f9.5)') e1ave,e1err,e1err*rteval_eff1
        write(6,'(''total energy ('',i4,'') ='',t22,f14.7,'' +-'',f11.7,f9.5)') nfprod-1,e2ave,e2err,e2err*rtevalg_eff1
        write(6,'(''total energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') e3ave,e3err,e3err*rteval_eff1
      endif

!JT      ovlp_ovlp_fn_av = ovlp_ovlp_fn_cum / (passes*nwalk)
!JT      call object_modified('ovlp_ovlp_fn_av')

      do 40 ifr=1,nforce
        peave=pecum(ifr)/wgcum(ifr)
        peiave=peicum(ifr)/wgcum(ifr)
        tpbave=tpbcum(ifr)/wgcum(ifr)
        tjfave=tjfcum(ifr)/wgcum(ifr)

        peerr=errg(pecum(ifr),pecm2(ifr),ifr)
        peierr=errg(peicum(ifr),peicm2(ifr),ifr)
        tpberr=errg(tpbcum(ifr),tpbcm2(ifr),ifr)
        tjferr=errg(tjfcum(ifr),tjfcm2(ifr),ifr)

c separate "magnetic energy" for quantum dots:
        if(ndim.eq.2) then
          temp=0.25d0*bext*bext/(we*we)
          tmave=(peave-peiave-emag)*temp
          tmerr=(peerr+peierr)*temp
          peave=peave-tmave-emag
c         peerr=peerr+tmerr     ! is this correct?
          peerr=peerr*(1-temp)+peierr*temp
        endif

        write(6,'(''potential energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') peave,peerr,peerr*rtevalg_eff1
        write(6,'(''interaction energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') peiave,peierr,peierr*rtevalg_eff1
        write(6,'(''jf kinetic energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') tjfave,tjferr,tjferr*rtevalg_eff1
        write(6,'(''pb kinetic energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') tpbave,tpberr,tpberr*rtevalg_eff1

        if(ndim.eq.2) then
          write(6,'(''radial mag. energy ='',t22,f14.7,'' +-'',f11.7,f9.5)') tmave,tmerr,tmerr*rtevalg_eff1
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
        write(6,'(''total energy diff'',i2,t22,f14.7,'' +-'',f11.7,f9.5)') ifr,fgave,fgerr,fgerr*rtevalg_eff1
  50  continue

      if(iperiodic.eq.0 .and. ncent.eq.1) then
        write(6,'(''<r>_av ='',t22,f14.7,'' +-'',f11.7,f9.5)') r1ave,r1err,r1err*rtevalg_eff1
        write(6,'(''<r2>_av ='',t22,f14.7,'' +-'',f11.7,f9.5)') r2ave,r2err,r2err*rtevalg_eff1
        write(6,'(''<ri>_av ='',t22,f14.7,'' +-'',f11.7,f9.5)') riave,rierr,rierr*rtevalg_eff1
      endif

      if(izigzag.ge.1) then
        write(6,'(''<ZigZag Amp> ='',t17,f12.7,'' +-'',f11.7,f9.5)') zzave(3),zzerr(3),zzerr(3)*rtevalg_eff1
        write(6,'(''<|ZigZag Amp|> ='',t17,f12.7,'' +-'',f11.7,f9.5)') zzave(1),zzerr(1),zzerr(1)*rtevalg_eff1
        write(6,'(''<ZigZag Amp^2> ='',t17,f12.7,'' +-'',f11.7,f9.5)') zzave(2),zzerr(2),zzerr(2)*rtevalg_eff1
        write(6,'(''<ZigZag Amp (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(6),zzerr(6),zzerr(6)*rtevalg_eff1
        write(6,'(''<|ZigZag Amp| (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(4),zzerr(4),zzerr(4)*rtevalg_eff1
        write(6,'(''<ZigZag Amp^2 (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(5),zzerr(5),zzerr(5)*rtevalg_eff1
        write(6,'(''<ZigZag rand Amp>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(9),zzerr(9),zzerr(9)*rtevalg_eff1
        write(6,'(''<|ZigZag rand Amp|>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(7),zzerr(7),zzerr(7)*rtevalg_eff1
        write(6,'(''<ZigZag rand Amp^2>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(8),zzerr(8),zzerr(8)*rtevalg_eff1
      endif

      if(ipr.gt.-2) write(11,'(3i5,f11.5,f7.4,f10.7,'' nstep,nblk,nconf,etrial,tau,taueff'')')
     &  nstep,iblk,nconf_global,etrial,tau,taucum(1)/wgcum(1)

      if(ifixe.ne.0 .or. ifourier.ne.0 .or. izigzag.ne.0) call den2dwrt(wgcum(1))

      call routines_write_final
      call reinit_routines_write_block
      call reinit_routines_write_final

      return
      end
