module accumulation_mod
! Not used at present

  use all_tools_mod
  use control_mod
  use montecarlo_mod
  use average_mod
  use print_mod

! Declaration of global variables and default values

  contains

! ==============================================================================
  subroutine zeres0_dmc_clean
! ------------------------------------------------------------------------------
! Description   : Initialize various quantities at beginning of run
! Description   : the initial values of energy psi etc. are calculated here
!
! Created       :
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, iw, ifr, ie, k, ip

! begin
  ipass=0

! set quadrature points
  if(nloc.gt.0) call gesqua (nquad,xq,yq,zq,wq)

  eigv=one
  eest=etrial
  nwalk=nconf
  wdsumo=nconf
  wgdsumo=nconf
  fprod=one

  call object_modified_by_index (nwalk_index)

  do i=0,nfprod
     wtgen(i)=nconf
     ff(i)=one
  enddo

      do 80 iw=1,nconf
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
          call hpsi(xoldw(1,1,iw,ifr),psidow(iw,ifr),psijow(iw,ifr),voldw(1,1,iw,ifr),div_vow(1,iw),d2ow(iw,ifr),peow(iw,ifr),eoldw(iw,ifr),denergy,ifr)
          pot_ee_oldw(:,iw,ifr) = pot_ee(:)
          if(ifr.eq.1) then
            call walksav_det(iw)
            call walksav_jas(iw)
          endif
          pwt(iw,ifr)=one
          do 72 ip=0,nwprod-1
   72       wthist(iw,ip,ifr)=one
   80 continue

  call zerest_dmc_clean

  end subroutine zeres0_dmc_clean

! ==============================================================================
  subroutine zerest_dmc_clean
! ------------------------------------------------------------------------------
! Description   : zero out all averages etc. after equilibration runs
!
! Created       :
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer ifr, i

! begin
  iblk=0
  iblk_proc=0

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
      zzcum(:)=zero

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
      zzcm2(:)=zero

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
      zzsum(:)=zero

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
      call alloc ('wgcum', wgcum, nforce)
      call alloc ('wgcum1', wgcum1, nforce)
      call alloc ('egcum', egcum, nforce)
      call alloc ('egcum1', egcum1, nforce)
      call alloc ('pecum', pecum, nforce)
      call alloc ('tpbcum', tpbcum, nforce)
      call alloc ('tjfcum', tjfcum, nforce)
      call alloc ('taucum', taucum, nforce)
      call alloc ('wgsum', wgsum, nforce)
      call alloc ('wsum1', wsum1, nforce)
      call alloc ('wgsum1', wgsum1, nforce)
      call alloc ('egsum', egsum, nforce)
      call alloc ('egsum1', egsum1, nforce)
      call alloc ('pesum', pesum, nforce)
      call alloc ('tpbsum', tpbsum, nforce)
      call alloc ('tjfsum', tjfsum, nforce)
      call alloc ('tausum', tausum, nforce)

      do 85 ifr=1,nforce
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
        tpbcum(ifr)=zero
        tjfcum(ifr)=zero
        pecm2(ifr)=zero
        tpbcm2(ifr)=zero
        tjfcm2(ifr)=zero
        pesum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
        fgcum(ifr)=zero
   85   fgcm2(ifr)=zero

      nbrnch=0

      try_int=0
      acc=0
      acc_int=0
      dr2ac=0
      dr2un=0
      nodecr=0

! Zero out estimators for charge density of atom.
      do 90 i=1,NRAD
        rprobup(i)=zero
        rprobdn(i)=zero
   90   rprob(i)=zero

      call grad_hess_jas_save

      return
  end subroutine zerest_dmc_clean

! ==============================================================================
  subroutine acues1_dmc_clean
! ------------------------------------------------------------------------------
! Description   :
!
! Created       :
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer ifr, ipmod
  real(dp) wgdsum1, nfpro

! begin

! statistical fluctuations without blocking
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

! collect block averages
      wsum=wsum+wsum1(1)
      wfsum=wfsum+wfsum1
      wdsum=wdsum+wdsumo
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

# if defined (MPI)
        eest=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
        eigv=dexp((etrial-eest)*(taucum(1)+tausum(1))/(wgcum(1)+wgsum(1)))
# else
        eest=egcum1(1)/wgcum1(1)
        eigv=dexp((etrial-eest)*(taucum(1)+tausum(1))/wgcum1(1))
# endif

        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest,egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wsum1(1)
      wgdsumo=wsum1(1)*fprod/ff(mod(ipass+1,nfprod))
      wtgen(ipmod)=wsum1(1)

! zero out step averages
      wfsum1=zero
      wdsum1=zero
      efsum1=zero
      do 40 ifr=1,nforce
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        esum1(ifr)=zero
   40   egsum1(ifr)=zero

      return
  end subroutine acues1_dmc_clean

! ==============================================================================
  subroutine acuest_dmc_clean
! ------------------------------------------------------------------------------
! Description   : routine to accumulate estimators for energy etc.
!
! Created       :
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  real(dp) :: rn_eff, w, w2
  real(dp) :: error, x, x2
  real(dp) :: errg
  integer i, npass, ifr
  real(dp) :: enow, efnow, ei1now, ei2now
  real(dp) :: rinow, r2now, r1now, r3now, r4now
  real(dp) :: zznow(nzzvars)
  real(dp) :: egnow, penow, tpbnow, tjfnow
  real(dp) :: peerr, tpberr, tjferr, fgerr
  real(dp) :: egave, peave, tpbave, tjfave, fgave
  integer  :: iegerr, ipeerr, itpber, itjfer, ifgerr
# if defined (MPI)
  integer ierr
  real(dp) :: w2sum, wf2sum, e2sum, ef2sum
  real(dp) :: ecollect, efcollect, e2collect, ef2collect
  real(dp) :: wcollect, wfcollect, w2collect, wf2collect
  real(dp) :: fcollect, f2collect
  real(dp) :: egcollect(nforce),wgcollect(nforce),pecollect(nforce)
  real(dp) :: tpbcollect(nforce),tjfcollect(nforce),eg2collect(nforce),wg2collect(nforce)
  real(dp) :: pe2collect(nforce),tpb2collect(nforce),tjf2collect(nforce),fsum(nforce)
  real(dp) :: f2sum(nforce),eg2sum(nforce),wg2sum(nforce),pe2sum(nforce),tpb2sum(nforce)
  real(dp) :: tjf2sum(nforce),taucollect(nforce)
# endif

! begin

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
      iblk_proc=iblk_proc+nproc
      npass=iblk_proc*nstep

!     wnow=wsum/nstep
!     wfnow=wfsum/nstep
      enow=esum/wsum
      efnow=efsum/wfsum
      ei1now=wfsum/wdsum
      ei2now=wgsum(1)/wgdsum
      rinow=risum/wgsum(1)
      r1now=r1sum/wgsum(1)
      r2now=r2sum/wgsum(1)
      r3now=r3sum/wgsum(1)
      r4now=r4sum/wgsum(1)
      zznow(:)=zzsum(:)/wgsum(1)

# if !defined (MPI)
      wcm2=wcm2+wsum**2
      wfcm2=wfcm2+wfsum**2
      ecm2=ecm2+esum*enow
      efcm2=efcm2+efsum*efnow
# endif
      ei1cm2=ei1cm2+ei1now**2
      ei2cm2=ei2cm2+ei2now**2
      r1cm2=r1cm2+r1sum*r1now
      r2cm2=r2cm2+r2sum*r2now
      r3cm2=r3cm2+r3sum*r3now
      r4cm2=r4cm2+r4sum*r4now
      ricm2=ricm2+risum*rinow
      zzcm2(:)=zzcm2(:)+zzsum(:)*zznow(:)

# if !defined (MPI)
      wcum=wcum+wsum
      wfcum=wfcum+wfsum
# endif
      wdcum=wdcum+wdsum
      wgdcum=wgdcum+wgdsum
# if !defined (MPI)
      ecum=ecum+esum
      efcum=efcum+efsum
# endif
      ei1cum=ei1cum+ei1now
      ei2cum=ei2cum+ei2now
      r1cum=r1cum+r1sum
      r2cum=r2cum+r2sum
      r3cum=r3cum+r3sum
      r4cum=r4cum+r4sum
      ricum=ricum+risum
      zzcum(:)=zzcum(:)+zzsum(:)

# if defined (MPI)
      w2sum=wsum**2
      wf2sum=wfsum**2
      e2sum=esum*enow
      ef2sum=efsum*efnow
# endif

      do 10 ifr=1,nforce

!       wgnow=wgsum(ifr)/nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum(ifr)/wgsum(ifr)
        tpbnow=tpbsum(ifr)/wgsum(ifr)
        tjfnow=tjfsum(ifr)/wgsum(ifr)

# if !defined (MPI)
        wgcm2(ifr)=wgcm2(ifr)+wgsum(ifr)**2
        egcm2(ifr)=egcm2(ifr)+egsum(ifr)*egnow
        pecm2(ifr)=pecm2(ifr)+pesum(ifr)*penow
        tpbcm2(ifr)=tpbcm2(ifr)+tpbsum(ifr)*tpbnow
        tjfcm2(ifr)=tjfcm2(ifr)+tjfsum(ifr)*tjfnow

        wgcum(ifr)=wgcum(ifr)+wgsum(ifr)
        egcum(ifr)=egcum(ifr)+egsum(ifr)
        pecum(ifr)=pecum(ifr)+pesum(ifr)
        tpbcum(ifr)=tpbcum(ifr)+tpbsum(ifr)
        tjfcum(ifr)=tjfcum(ifr)+tjfsum(ifr)

# else
        wg2sum(ifr)=wgsum(ifr)**2
        eg2sum(ifr)=egsum(ifr)*egnow
        pe2sum(ifr)=pesum(ifr)*penow
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
! Warning temp fix
        wgsum(ifr)=wgcollect(ifr)
        egsum(ifr)=egcollect(ifr)
        egnow=egsum(ifr)/wgsum(ifr)
        wgcum(ifr)=wgcum(ifr)+wgcollect(ifr)
        egcum(ifr)=egcum(ifr)+egcollect(ifr)
        taucum(ifr)=taucum(ifr)+taucollect(ifr)
  12  continue


!      call compute_averages_walk_block !JT
!      call compute_errors  !JT

      call mpi_allreduce(pesum,pecollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tpbsum,tpbcollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tjfsum,tjfcollect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      call mpi_allreduce(wg2sum,wg2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(eg2sum,eg2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(pe2sum,pe2collect,nforce,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
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

!JT      if(idtask.ne.0) goto 17

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
        tpbcm2(ifr)=tpbcm2(ifr)+tpb2collect(ifr)
        tjfcm2(ifr)=tjfcm2(ifr)+tjf2collect(ifr)

        pecum(ifr)=pecum(ifr)+pecollect(ifr)
        tpbcum(ifr)=tpbcum(ifr)+tpbcollect(ifr)
        tjfcum(ifr)=tjfcum(ifr)+tjfcollect(ifr)
# endif

        call grad_hess_jas_cum(wgsum(ifr),egnow)

        if(iblk.eq.1) then
          egerr=0
          peerr=0
          tpberr=0
          tjferr=0
         else
          egerr=errg(egcum(ifr),egcm2(ifr),ifr)
          peerr=errg(pecum(ifr),pecm2(ifr),ifr)
          tpberr=errg(tpbcum(ifr),tpbcm2(ifr),ifr)
          tjferr=errg(tjfcum(ifr),tjfcm2(ifr),ifr)
        endif

        egave=egcum(ifr)/wgcum(ifr)
        peave=pecum(ifr)/wgcum(ifr)
        tpbave=tpbcum(ifr)/wgcum(ifr)
        tjfave=tjfcum(ifr)/wgcum(ifr)

        if(ifr.gt.1) then
# if !defined (MPI)
          fgcum(ifr)=fgcum(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))
          fgcm2(ifr)=fgcm2(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
# else
          fgcum(ifr)=fgcum(ifr)+fsum(ifr)
          fgcm2(ifr)=fgcm2(ifr)+f2sum(ifr)
# endif
          fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
          fgave=fgave/deltot(ifr)
          if(iblk.eq.1) then
            fgerr=0
            ifgerr=0
           else
            fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
            fgerr=fgerr/abs(deltot(ifr))
            ifgerr=nint(100000* fgerr)
          endif
        endif

# if !defined (MPI)
        taucum(ifr)=taucum(ifr)+tausum(ifr)
# endif

! write out header first time
  if (iblk==1 .and. ifr ==1) then
    if (nforce == 1) then
      write(6,'(a)') '    passes    energy      averaged energy      weights  oldest walker'
    else
    write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66,''tjfave'',t72,''(tjferr)'',t83,''fgave'',t89,''(fgerr)'',t101,''npass'',t111,''wgsum'',t121,''ioldest'')')
    endif
  endif

! write out current values of averages etc.
  iegerr=nint(100000* egerr)
  ipeerr=nint(100000* peerr)
  itpber=nint(100000*tpberr)
  itjfer=nint(100000*tjferr)

# if !defined (MPI)
  if (ifr == 1) then
    if (nforce == 1) then
     write(6,'(i10,1x,f10.5,1x,f10.5,a,f8.5,2i10)') npass,egnow,egave,' +-',egerr,nint(wgsum(ifr)/nstep_total),ioldest
    else
     write(6,'(f10.5,4(f10.5,''('',i5,'')''),17x,3i10)') egnow,egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,npass,nint(wgsum(ifr)),ioldest
    endif
  else
    write(6,'(f10.5,5(f10.5,''('',i5,'')''),10x,i10)') egnow,egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,fgave,ifgerr,nint(wgsum(ifr))
  endif

   10 continue

# else

  if(ifr.eq.1) then
    if (nforce == 1) then
     write(6,'(i10,1x,f10.5,1x,f10.5,a,f8.5,3i10)') npass,egcollect(ifr)/wgcollect(ifr),egave,' +-',egerr,nint(wgcollect(ifr)/nproc),ioldest
    else
          write(6,'(f10.5,4(f10.5,''('',i5,'')''),17x,3i10)') egcollect(ifr)/wgcollect(ifr),egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,npass,nint(wgcollect(ifr)/nproc),ioldest
    endif
   else
          write(6,'(f10.5,5(f10.5,''('',i5,'')''),10x,i10)') egcollect(ifr)/wgcollect(ifr),egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,fgave,ifgerr,nint(wgcollect(ifr)/nproc)
   endif

   15 continue
# endif


      eloc_av = egave                                 !JT
      call object_modified_by_index (eloc_av_index)   !JT

      call systemflush(6)

! zero out xsum variables for metrop

      wsum=zero
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
      zzsum(:)=zero

      do 20 ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
   20   tausum(ifr)=zero

      call systemflush(6)

      return
      end subroutine acuest_dmc_clean

! ==============================================================================
  subroutine finwrt_dmc_clean
! ------------------------------------------------------------------------------
! Description   : routine to print out final results
!
! Created       :
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j, k, ifr, iz
  character*80 fmt
  real (dp) :: rn_eff, error, errorn, errc, errf, errg, errc1
  real (dp) :: errf1, errg1, errw, errw1
  real (dp) :: x,x2,w,w2,rn
  real (dp) :: eval, eval_proc_eff,evalf_proc_eff, evalg_proc_eff
  real (dp) :: rtpass1
  real (dp) :: rteval, rteval_proc_eff1, rtevalf_proc_eff1, rtevalg_proc_eff1
!  real (dp) :: delr, delri, term
  real (dp) :: accav, accavn, wave, wfave, eave, efave
  real (dp) :: werr, wferr, werr1, wferr1
  real (dp) :: eerr, eferr, eerr1, eferr1
  real (dp) :: wgave, wgerr, wgerr1
  real (dp) :: egave, egerr1
  real (dp) :: peave, tpbave, tjfave
  real (dp) :: peerr, tpberr, tjferr
  real (dp) :: fgave, fgerr
  real(dp)  :: passes, pass_proc, eval_proc
  real(dp)  :: rtpass_proc1, rteval_proc
  real(dp)  :: tcor_g

# if !defined (MPI)
  real (dp) :: erric, erric1, errig, errori
  real (dp) :: e1ave, e2ave, e3ave, ei1ave, ei2ave, ei3ave
  real (dp) :: e1err, e2err, e3err, ei1err, ei2err, ei3err
  real (dp) :: r2ave, r2err, riave, rierr, r1ave, r1err, r3ave, r3err, r4ave, r4err
  real (dp) :: zzave(nzzvars), zzerr(nzzvars)
# endif

# if defined (MPI)
  integer  :: ierr
  real(dp) :: e1collect, e21collect, ef21collect, ef1collect
  real(dp) :: w1collect, w21collect, wf21collect, wf1collect
  integer  :: nodecr_collect
  real(dp) :: try_int_collect, acc_collect, acc_int_collect
  real(dp) :: eg1collect(nforce),eg21collect(nforce),wg1collect(nforce),wg21collect(nforce),rprobcollect(NRAD)
# endif

! statement functions for error calculation
  rn_eff(w,w2)=w**2/w2

  error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
  errorn(x,x2,rn)=dsqrt(max((x2/rn-(x/rn)**2)/(rn-1),0.d0))
  errc(x,x2)=error(x,x2,wcum,wcm2)
  errf(x,x2)=error(x,x2,wfcum,wfcm2)
  errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))
  errc1(x,x2)=error(x,x2,wcum1,wcm21)
  errf1(x,x2)=error(x,x2,wfcum1,wfcm21)
  errg1(x,x2,i)=error(x,x2,wgcum1(i),wgcm21(i))
  errw(x,x2)=errorn(x,x2,dfloat(iblk_proc))/nstep
  errw1(x,x2)=errorn(x,x2,pass_proc)
# if !defined (MPI)
  errori(x,x2,w,w2,rn)=dsqrt(max((x2/rn-(x/rn)**2)/(rn_eff(w,w2)-1),0.d0))
  erric(x,x2)=errori(x,x2,wcum,wcm2,dfloat(iblk_proc))
  erric1(x,x2)=errori(x,x2,wcum1,wcm21,pass_proc)
  errig(x,x2)=errori(x,x2,wgcum(1),wgcm2(1),dfloat(iblk_proc))
# endif

  passes=dfloat(iblk)*dfloat(nstep)
  eval=nconf*passes
  pass_proc=dfloat(iblk_proc)*dfloat(nstep)
  eval_proc=nconf*pass_proc

! Either the next 3 lines or the 3 lines following them could be used.
! They should give nearly (but not exactly) the same result.
! Strictly the 1st 3 are for step-by-step quantities and the last 3 for blk-by-blk
!     eval_proc_eff=nconf*rn_eff(wcum1,wcm21)
!     evalf_proc_eff=nconf*rn_eff(wfcum1,wfcm21)
!     evalg_proc_eff=nconf*rn_eff(wgcum1(1),wgcm21(1))
      eval_proc_eff=nconf*nstep*rn_eff(wcum,wcm2)
      evalf_proc_eff=nconf*nstep*rn_eff(wfcum,wfcm2)
      evalg_proc_eff=nconf*nstep*rn_eff(wgcum(1),wgcm2(1))
      rtpass1=dsqrt(passes-1)
      rtpass_proc1=dsqrt(pass_proc-1)
      rteval=dsqrt(eval)
      rteval_proc=dsqrt(eval_proc)
      rteval_proc_eff1=dsqrt(eval_proc_eff-1)
      rtevalf_proc_eff1=dsqrt(evalf_proc_eff-1)
      rtevalg_proc_eff1=dsqrt(evalg_proc_eff-1)


# if defined (MPI)

!JT      write(6,'(''passes,eval,pass_proc,eval_proc,eval_proc_eff,evalf_proc_eff,evalg_proc_eff'',19f9.0)') passes,eval,pass_proc,eval_proc,eval_proc_eff,evalf_proc_eff,evalg_proc_eff

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

! Collect radial charge density for atoms
      if(iperiodic.eq.0) then
        call mpi_allreduce(rprob,rprobcollect,NRAD,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        do 2 i=1,NRAD
    2     rprob(i)=rprobcollect(i)
        call mpi_allreduce(rprobup,rprobcollect,NRAD,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        do 3 i=1,NRAD
    3     rprobup(i)=rprobcollect(i)
        call mpi_allreduce(rprobdn,rprobcollect,NRAD,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        do 4 i=1,NRAD
    4     rprobdn(i)=rprobcollect(i)
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

!     write(11,'(4i5,f11.5,f7.4,f10.7,
!    &'' nstep,nblk,nblkeq,nconf,etrial,tau,taueff'')')
!    &nstep,nblk,nblkeq,nconf,etrial,tau,taucum(1)/wgcum(1)

!JT      if(ipr.gt.-2) write(11,'(3i5,f11.5,f7.4,f10.7,'' nstep,nblk,nconf,etrial,tau,taueff'')') nstep,iblk,nconf,etrial,tau,taucum(1)/wgcum(1)

!JT      if(idtask.ne.0) return
# endif

! Write out radial charge density for atoms
!JT      if (print_radial_probability) then  !JT
!JT      if(iperiodic.eq.0 .and. ncent.eq.1) then
!JT        if(ndim.eq.3) write(6,*)'  r  4*pi*r^2*rho'
!JT        if(ndim.eq.2) write(6,*)'  r  2*pi*r*rho'
!JT        delr=one/delri
!JT        term=one/(wgcum(1)*delr)
!JT        do 5 i=1,NRAD
!JT    5     write(6,'(f5.3,3f10.6)') delr*(i-half),rprob(i)*term,rprobup(i)*term,rprobdn(i)*term
!JT      endif
!JT      endif


! Write final results
  write(6,'(a)') 'Final results:'

! Age of walkers:
  if(idmc >= 0) then
       write(6,'(a,2i4)') 'Age of oldest walker in the last generation, in any generation:',ioldest,ioldestmx
!JT        write(6,'(/,''ages of walkers are:'')')
!JT        write(6,'(10i4)') (iage(i),i=1,nwalk)
        do i=1,nwalk
          if(iage(i).gt.50) then
            write(6,'(a)') 'warning: some walkers are old:'
            write(6,'(i4,i6,f10.4,99f8.4)') i,iage(i),eoldw(i,1),((xoldw(k,j,i,1),k=1,ndim),j=1,nelec)
            write(6,'(99f8.4)') ((voldw(k,j,i,1),k=1,ndim),j=1,nelec)
          endif
        enddo
  endif

! Diffusion analysis
      write(6,'(a,f10.5)') 'Average of the squares drift-dif moves for accepted steps = ',dr2ac/try_int
      write(6,'(a,f10.5)') 'Average of the squares drift-dif moves for all      steps = ',dr2un/try_int

      write(6,'(a,20f7.4)') 'Effective time step of diffusion = ',(taucum(ifr)/wgcum(ifr),ifr=1,nforce)

      accav=acc/try_int
      accavn=acc_int/try_int

      wave=wcum/pass_proc
      wfave=wfcum/pass_proc
      eave=ecum/wcum
      efave=efcum/wfcum

# if !defined (MPI)
      ei1ave=wfcum/wdcum
      ei2ave=wgcum(1)/wgdcum
      ei3ave=ei3cum/passes

      r1ave=r1cum/(wgcum(1)*nelec)
      r2ave=r2cum/(wgcum(1)*nelec)
      r3ave=r3cum/(wgcum(1)*nelec)
      r4ave=r4cum/(wgcum(1)*nelec)
      riave=ricum/(wgcum(1)*nelec)
      zzave(:)=zzcum(:)/wgcum(1)
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
        r3err=0
        r4err=0
        rierr=0
        zzerr(:)=0.d0
       else
# endif
        werr=errw(wcum,wcm2)
        wferr=errw(wfcum,wfcm2)
        werr1=errw1(wcum1,wcm21)
        wferr1=errw1(wfcum1,wfcm21)
        eerr=errc(ecum,ecm2)
        eferr=errf(efcum,efcm2)
        eerr1=errc1(ecum1,ecm21)
        eferr1=errf1(efcum1,efcm21)
# if !defined (MPI)
        ei1err=erric(ei1cum,ei1cm2)
        ei2err=errig(ei2cum,ei2cm2)
        ei3err=erric1(ei3cum,ei3cm2)
        r1err=errg(r1cum,r1cm2,1)/nelec
        r2err=errg(r2cum,r2cm2,1)/nelec
        r3err=errg(r3cum,r3cm2,1)/nelec
        r4err=errg(r4cum,r4cm2,1)/nelec
        rierr=errg(ricum,ricm2,1)/nelec
        do iz = 1,nzzvars
          zzerr(iz)=errg(zzcum(iz),zzcm2(iz),1)
        enddo
      endif

      e1err=dlog((ei1ave+ei1err)/(ei1ave-ei1err))/(2*taucum(1)/wgcum(1))
      e2err=dlog((ei2ave+ei2err)/(ei2ave-ei2err))/(2*taucum(1)/wgcum(1))
      e3err=dlog((ei3ave+ei3err)/(ei3ave-ei3err))/(2*taucum(1)/wgcum(1))
# endif

      write(6,'(a,i9,f10.6)') 'Number and fraction of node crossings:',nodecr,dfloat(nodecr)/try_int
      write(6,'(a,2f10.6)') 'Acceptance:', accav,accavn

      if(idmc.ge.0) then

# if !defined (MPI)
        write(6,'(a,i6,2f9.0)') 'Actual number of branches and expected numbers for 0 and inf corr time:',nbrnch,nconf*passes*(dlog(one+eerr1*rteval*taucum(1)/wgcum(1))/dlog(two))**2,nconf*passes*(dlog(one+eerr1*rteval*taucum(1)/wgcum(1))/dlog(two))
# endif
        write(6,'(a,i5)') 'Number of walkers at end of run =',nwalk

        write(6,'(a)') 'Effective/actual number of walkers ratio (with population-control correction):'
        write(6,'(''nwalk_eff/nwalk (   0) ='',2f6.3)') rn_eff(wcum1,wcm21)/pass_proc,rn_eff(wcum,wcm2)/iblk_proc
        write(6,'(''nwalk_eff/nwalk (   1) ='',2f6.3)') rn_eff(wfcum1,wfcm21)/pass_proc,rn_eff(wfcum,wfcm2)/iblk_proc
        write(6,'(''nwalk_eff/nwalk ('',i4,'') ='',2f6.3)') nfprod, rn_eff(wgcum1(1),wgcm21(1))/pass_proc,rn_eff(wgcum(1),wgcm2(1))/iblk_proc

      endif

      write(fmt,'(''(/,a16,2x,a'',i3,'')'')') len_trim(title)
      write(6,fmt) mode,title
      write(6,'(''nconf*passes'',t19,''passes  nconf nstep  nblk nblkeq nproc  tau    taueff'',/,2f12.0,2i6,i7,2i5,2f9.5)') eval,passes,nconf,nstep,iblk,nblkeq,nproc,tau,taucum(1)/wgcum(1)

      write(6,'(a)') 'quantity (pop. corr.)     average       rms error sigma*T_cor  sigma   T_cor'

      if(idmc.ge.0) then

        write(6,'(''weights      (   0) ='',t22,f14.7,'' +-'',f11.7,1x,2f9.5,f8.2)') wave,werr,werr*rtpass_proc1,werr1*rtpass_proc1,(werr/werr1)**2
        write(6,'(''weights      (   1) ='',t22,f14.7,'' +-'',f11.7,1x,2f9.5,f8.2)') wfave,wferr,wferr*rtpass_proc1,wferr1*rtpass_proc1,(wferr/wferr1)**2

        do 20 ifr=1,nforce
          wgave=wgcum(ifr)/pass_proc
          wgerr=errw(wgcum(ifr),wgcm2(ifr))
          wgerr1=errw1(wgcum1(ifr),wgcm21(ifr))
          write(6,'(''weights      ('',i4,'') ='',t22,f14.7,'' +-'',f11.7,1x,2f9.5,f8.2)') nfprod,wgave,wgerr,wgerr*rtpass_proc1,wgerr1*rtpass_proc1,(wgerr/wgerr1)**2
  20    continue

        write(6,'(''total energy (   0) ='',t24,f12.7,'' +-'',f11.7,1x,2f9.5,f8.2)') eave,eerr,eerr*rteval_proc_eff1,eerr1*rteval_proc_eff1,(eerr/eerr1)**2
        write(6,'(''total energy (   1) ='',t24,f12.7,'' +-'',f11.7,1x,2f9.5,f8.2)') efave,eferr,eferr*rtevalf_proc_eff1,eferr1*rtevalf_proc_eff1,(eferr/eferr1)**2
      endif

      do ifr=1,nforce
        egave=egcum(ifr)/wgcum(ifr)
        egerr=errg(egcum(ifr),egcm2(ifr),ifr)
        egerr1=errg1(egcum1(ifr),egcm21(ifr),ifr)
        tcor_g = (egerr/egerr1)**2
        write(6,'(''total energy ('',i4,'') ='',t24,f12.7,'' +-'',f11.7,1x,2f9.5,f8.2)') nfprod,egave,egerr,egerr*rtevalg_proc_eff1,egerr1*rtevalg_proc_eff1,tcor_g
      enddo

# if !defined (MPI)
      if(idmc.ge.0) then
        write(6,'(''total energy (   0) ='',t24,f12.7,'' +-'',f11.7,1x,f9.5, ''  (growth estimator)'')') e1ave,e1err,e1err*rteval_proc_eff1
        write(6,'(''total energy ('',i4,'') ='',t24,f12.7,'' +-'',f11.7,1x,f9.5, ''  (growth estimator)'')') nfprod-1,e2ave,e2err,e2err*rtevalg_proc_eff1
        write(6,'(''total energy        ='',t24,f12.7,'' +-'',f11.7,1x,f9.5, ''  (growth estimator)'')') e3ave,e3err,e3err*rteval_proc_eff1
      endif
# endif

      do 40 ifr=1,nforce
        peave=pecum(ifr)/wgcum(ifr)
        tpbave=tpbcum(ifr)/wgcum(ifr)
        tjfave=tjfcum(ifr)/wgcum(ifr)

        peerr=errg(pecum(ifr),pecm2(ifr),ifr)
        tpberr=errg(tpbcum(ifr),tpbcm2(ifr),ifr)
        tjferr=errg(tjfcum(ifr),tjfcm2(ifr),ifr)

        write(6,'(a,t24,f12.7,a,f11.7,1x,f9.5)') 'potential  energy   =',peave,' +-',peerr,peerr*rtevalg_proc_eff1
        write(6,'(a,t24,f12.7,a,f11.7,1x,f9.5,a,f9.5)') 'JF kinetic energy   =',tjfave,' +-',tjferr,tjferr*rtevalg_proc_eff1,'  Virial ratio =',peave/tjfave
        write(6,'(a,t24,f12.7,a,f11.7,1x,f9.5,a,f9.5)') 'PB kinetic energy   =',tpbave,' +-',tpberr,tpberr*rtevalg_proc_eff1,'  Virial ratio =',peave/tpbave
  40  continue

      do 50 ifr=2,nforce
        fgave=egcum(ifr)/wgcum(ifr)-egcum(1)/wgcum(1)
        fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
!       fgave=fgave/deltot(ifr)
!       fgerr=fgerr/abs(deltot(ifr))
        write(6,'(''total energy diff'',i2,t24,f12.7,'' +-'',f11.7,1x,f9.5)') ifr,fgave,fgerr,fgerr*rtevalg_proc_eff1
  50  continue

# if !defined (MPI)
!     These are not being collected at the moment (they are for ndim=2)
      if(iperiodic.eq.0 .and. ncent.eq.1) then
        write(6,'(''<r>_av ='',t24,f12.7,'' +-'',f11.7,1x,f9.5)') r1ave,r1err,r1err*rtevalg_proc_eff1
        write(6,'(''<r2>_av ='',t24,f12.7,'' +-'',f11.7,1x,f9.5)') r2ave,r2err,r2err*rtevalg_proc_eff1
        write(6,'(''<r3>_av ='',t24,f12.4,'' +-'',f11.5,1x,f9.2)') r3ave,r3err,r3err*rtevalg_proc_eff1
        write(6,'(''<r4>_av ='',t24,f12.2,'' +-'',f11.3,1x,f9.1)') r4ave,r4err,r4err*rtevalg_proc_eff1
        write(6,'(''<ri>_av ='',t24,f12.7,'' +-'',f11.7,1x,f9.5)') riave,rierr,rierr*rtevalg_proc_eff1
      endif

      if(izigzag.ge.1) then
        call print_zigzag_vars(zzave,zzerr,rtevalg_proc_eff1)
!       write(6,'(''<ZigZag Amp> ='',t17,f12.7,'' +-'',f11.7,f9.5)') zzave(3),zzerr(3),zzerr(3)*rtevalg_proc_eff1
!       write(6,'(''<|ZigZag Amp|> ='',t17,f12.7,'' +-'',f11.7,f9.5)') zzave(1),zzerr(1),zzerr(1)*rtevalg_proc_eff1
!       write(6,'(''<ZigZag Amp^2> ='',t17,f12.7,'' +-'',f11.7,f9.5)') zzave(2),zzerr(2),zzerr(2)*rtevalg_proc_eff1
!       write(6,'(''<ZigZag Amp (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(6),zzerr(6),zzerr(6)*rtevalg_proc_eff1
!       write(6,'(''<|ZigZag Amp| (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(4),zzerr(4),zzerr(4)*rtevalg_proc_eff1
!       write(6,'(''<ZigZag Amp^2 (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(5),zzerr(5),zzerr(5)*rtevalg_proc_eff1
!       write(6,'(''<ZigZag rand Amp>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(9),zzerr(9),zzerr(9)*rtevalg_proc_eff1
!       write(6,'(''<|ZigZag rand Amp|>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(7),zzerr(7),zzerr(7)*rtevalg_proc_eff1
!       write(6,'(''<ZigZag rand Amp^2>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzave(8),zzerr(8),zzerr(8)*rtevalg_proc_eff1
      endif

      if(ipr.gt.-2) write(11,'(3i5,f11.5,f7.4,f10.7,'' nstep,nblk,nconf,etrial,tau,taueff'')')nstep,iblk,nconf,etrial,tau,taucum(1)/wgcum(1)
# endif

! Print warnings
  if ((idmc < 0 .and. accav < 0.3) .or. (idmc >= 0 .and. accav < 0.7)) then
    write(6,'(a,f10.6,a)') 'Warning: low acceptance =', accav,'. Reduce time step tau.'
  endif

  if (nstep < 10.d0*tcor_g) then
    write(6,'(a,i6,a,f8.2)') 'Warning: number of step per block =',nstep,' is less than 10 times the autocorrelation time =',tcor_g
    write(6,'(a)') 'Warning: the errors are underestimated. Increase nstep.'
  endif

  end subroutine finwrt_dmc_clean

end module accumulation_mod
