function func(ndata2,nparm,parm,diff,iflag)
! Written by Cyrus Umrigar, modified by Claudia Filippi
! iflag=ioffset was introduced by Peter to allow quench_anneal to work in parallel.
! Note that in the parallel version of CHAMP we are doing the parallelism in CHAMP
! and so we are using just the serial version of quench_anneal.
! In the parallel version of quench_anneal func is called with iflag>=0 to compute diffs offset by iflag=ioffset,
! and with -1 to do the sum of squares.
! In the serial version of quench_anneal func is called with iflag=0 to compute diffs,
! and with -1 to do the sum of squares.
  use constants_mod
  use mpi_mod
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
  use confg_mod
  use const_mod
  use mpioffset_mod
  implicit real*8(a-h,o-z)

  common /fcn_calls/icalls

  logical called_by_qa
  common /quenchsim_pr/ ipr_com,called_by_qa

  dimension velocity(3,nelec),div_v(nelec)

  dimension parm(*),diff(*)

  data func_sav/1.d99/
  save func_sav, eavri

  if(iflag.eq.0) then

    icalls=icalls+1

    do iparm=1,nparml
      coef(iwbasi(iparm),iworb(iparm),1)=parm(iparm)
    enddo
    do iparm=1,nparme
      zex(iwbase(iparm),1)=parm(nparml+iparm)
      ict=ictype_basis(iwbase(iparm))
      irb=iwrwf2(iwbase(iparm))
      zex2(irb,ict,1)=zex(iwbase(iparm),1)
    enddo
    do iparm=1,nparmcsf
      csf_coef(iwcsf(iparm),1)=parm(nparml+nparme+iparm)
    enddo

!c Calculate coefs to construct the piece of the orbital that comes
!c from basis fns that are related by symmetry.  This is needed to
!c impose cusp conditions when there is more than one atom.
!    if(nloc.eq.0.and.numr.le.0) call equiv_bas

    if(nparms.eq.1) then
!     scalek(1)=parm(nparml+nparme+nparmd+1)
      scalek(1)=parm(nparml+nparme+nparmcsf+1)
! If we are varying scalek reset dependent constants
!     if(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17) then
!       call set_scale_dist(0,1)
!     endif
    endif

!   if(nparmg.eq.1) a21   =parm(nparml+nparme+nparmd+nparms+1)
    if(nparmg.eq.1) a21   =parm(nparml+nparme+nparmcsf+nparms+1)
    if(ijas.eq.1) then
      if(nparmj.ge.1) cjas2(1)=parm(nparm)
    elseif(ijas.eq.2) then
      ntmp=nparmj
      do isp=nspin1,nspin2
        do iparm=1,nparma(isp)
          a1(iwjasa(iparm,isp),isp,1)=parm(nparm-ntmp+iparm)
        enddo
        ntmp=ntmp-nparma(isp)
      enddo
      do isp=nspin1,nspin2
        do iparm=1,nparmb(isp)
          a2(iwjasb(iparm,isp),isp,1)=parm(nparm-ntmp+iparm)
        enddo
        ntmp=ntmp-nparmb(isp)
      enddo
    elseif(ijas.eq.3) then
      ntmp=nparmj
      do iparm=1,nparma(1)
        a(iwjasa(iparm,1),1)=parm(nparm-ntmp+iparm)
      enddo
      ntmp=ntmp-nparma(1)
      do isp=nspin1,nspin2b
        do iparm=1,nparmb(isp)
          b(iwjasb(iparm,isp),isp,1)=parm(nparm-ntmp+iparm)
        enddo
        ntmp=ntmp-nparmb(isp)
      enddo
      do it=1,nctype
        do iparm=1,nparmc(it)
          c(iwjasc(iparm,it),it,1)=parm(nparm-ntmp+iparm)
        enddo
        ntmp=ntmp-nparmc(it)
      enddo
      if(ifock.gt.0) then
        do it=1,nctype
          do iparm=1,nparmf(it)
            fck(iwjasf(iparm,it),it,1)=parm(nparm-ntmp+iparm)
          enddo
          ntmp=ntmp-nparmf(it)
          if(ifock.gt.2) then
            call scale3(1,it)
            if(ifock.eq.4) call scale20(1,it)
          endif
        enddo ! nctype
      endif
    elseif(ijas.ge.4.and.ijas.le.6) then
      ntmp=nparmj
      do it=1,nctype
        do iparm=1,nparma(it)
          a4(iwjasa(iparm,it),it,1)=parm(nparm-ntmp+iparm)
        enddo
        ntmp=ntmp-nparma(it)
      enddo
      do isp=nspin1,nspin2b
        do iparm=1,nparmb(isp)
          b(iwjasb(iparm,isp),isp,1)=parm(nparm-ntmp+iparm)
        enddo
        ntmp=ntmp-nparmb(isp)
      enddo
      do it=1,nctype
        do iparm=1,nparmc(it)
          c(iwjasc(iparm,it),it,1)=parm(nparm-ntmp+iparm)
        enddo
        ntmp=ntmp-nparmc(it)
      enddo
    endif

! Foll. is now in cuspco
!   if(ijas.eq.2) then
!     if(nspin1.eq.1 .and. nspin2.eq.1) then
!       if(nelec.eq.2) then
!         aa1=a1(2,1,1)
!        else
!         aa1=d1b4*(three*(nup+ndn)-2)*a1(2,1,1)
!       endif
!      elseif(nspin1.eq.2 .and. nspin2.eq.2) then
!       aa1=(nup-1)*a1(2,2,1)
!      elseif(nspin1.eq.1 .and. nspin2.eq.2) then
!       aa1=half*((nup+ndn)*a1(2,1,1)+(nup+ndn-2)*a1(2,2,1))
!      elseif(nspin1.eq.1 .and. nspin2.eq.3) then
!       aa1=nup*a1(2,1,1)+(ndown-1)*a1(2,3,1)
!     endif
!    elseif(ijas.eq.3) then
!     aa1=a(1,1)
!    elseif(ijas.ge.4.and.ijas.le.6) then
!     aa1=a4(1,1,1)
!    else
!     aa1=zero
!   endif

! If we are varying scalek, a or b parms reset dependent constants
!    if((ijas.eq.4.or.ijas.eq.5).and.(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17)) then
      call set_scale_dist(-1,1)    ! always call this in the analytical scalek opt
!    endif

    if(icusp2.ge.1.and.isc.le.7) then
      do isp=nspin1,nspin2
        if(ijas.eq.2) call cuspexact2(a1(1,isp,1),a2(1,isp,1))
      enddo
      if(ijas.eq.3) call cuspexact3(0)
      if(ijas.ge.4.and.ijas.le.6) call cuspexact4(0,1)
    endif

    if(ijas.eq.2 .and. iabs(icusp2).ge.2) then
      do isp=nspin1,nspin2
        ishft=ncuspc*(isp-nspin1)
        if(nspin2.eq.2 .and. (nup-ndn).ne.0) a1(2,2,1)=a1(2,1,1)
        if(nspin2.eq.3 .and. isp.eq.nspin2) a1(2,3,1)=((ndn-nup)*a1(2,1,1)+(nup-1)*a1(2,2,1))/(ndn-1)
        call cuspcheck2(scalek(1),a1(1,isp,1),a2(1,isp,1),diff(ndata+ishft+1),isp,nspin1,ncuspc,0)
        do i=1,ncuspc
          diff(ndata+ishft+i)=diff(ndata+ishft+i)*cuspwt
        enddo
      enddo
    elseif(ijas.eq.3) then
      call cuspcheck3(diff(ndata+1),0)
      do i=1,ncuspc+nfockc
        diff(ndata+i)=diff(ndata+i)*cuspwt
      enddo
    endif

! If icusp>=0 impose cusp condition for each s orbital.
! In any case calculate cusp-violation penalty. Should be 0 if icusp>=0.
    ishft=ncuspc*(nspin2-nspin1+1)+nfockc
!   if((nloc.eq.0. .or. nloc.eq.6) .and. numr.le.0) call cuspco(diff(ndata+ishft+1),0)
    if(icusp.ge.0) call cuspco(diff(ndata+ishft+1),0)

    do i=1,ncent*norbc
      diff(ndata+ishft+i)=diff(ndata+ishft+i)*cuspwt
    enddo
    do i=1,necn
      coef(iebasi(1,i),ieorb(1,i),1)=sign(one,dfloat(ieorb(2,i))) * coef(iebasi(2,i),iabs(ieorb(2,i)),1)
     enddo
    do i=1,nebase
      zex(iebase(1,i),1)=zex(iebase(2,i),1)
      ict=ictype_basis(iebase(1,i))
      irb=iwrwf2(iebase(1,i))
      zex2(irb,ict,1)=zex(iebase(2,i),1)
    enddo

!   do 80 i=1,iabs(nedet)
!     cdet(iedet(1,i),1)=zero
!     do 80 j=2,icsf(i)
!  80   cdet(iedet(1,i),1)=cdet(iedet(1,i),1)+frac(j,i)*
!  &    sign(one,dfloat(iedet(j,i)))*cdet(iabs(iedet(j,i)),1)

    if(ipos+idcds+idcdu+idcdt+id2cds+id2cdu+id2cdt+idbds+idbdu+idbdt.gt.0 .and. ijas.eq.2) then
      ishft=ndata+ncuspc*(nspin2-nspin1+1)+nfockc+ncent*norbc+1
      icalcul_diff=1
      do isp=nspin1,nspin2
        call checkjas2(scalek(1),isp,ncnstr,diff(ishft),0,0,icalcul_diff)
        ishft=ishft+ncnstr
      enddo
    endif

! Data distribution for parallel run.
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

! Here we are calculating numerical derivs. wrt. wavefn. params so turn igradhess off before calling hpsi.
    igradhess=0

! For a serial run this reduces to a "do 123 i=1,ndata"
    do i=idispls(idtask)+1,idispls(idtask+1)
      iconfg=i
      call hpsi(x(1,1,i),psid(i),psij(i),velocity,div_v,d2psi,pe,pei,energy,denergy,1)
      diff(i)=energy-eguess
    enddo
    if(nparml+nparme+nparmcsf.eq.0) isaved=1

    if(index(mode,'mpi').ne.0) call func_qmc_mpi(diff)

    do i=1,ndata
      uwdiff(i)=diff(i)
    enddo

    if(mod(irewgt,100).ne.0) call rewght(diff)

    eavr=zero
    wavr=zero
    do i=1,ndata
      wavr=wavr+wght(i)
      eavr=eavr+(uwdiff(i)+eguess)*wght(i)
    enddo
    eavr=eavr/wavr

    eavri=one
    if(irewgt.ge.100) eavri=-one/eavr
    do i=1,ndata
      diff(i)=diff(i)*eavri
    enddo

!   err=zero
!   do 130 i=1,ndata2
! 130   err=err+diff(i)**2
!   func=err
    func=0  ! needed for ifort compiler but not for gfortran

  else ! iflag /= 0

    err=zero
    do i=1,ndata2
      err=err+diff(i)**2
    enddo
    func=err

! We are optimizing fluctuations from eguess rather than from the average.
! However, if the foll. lines are not commented out, we reset eguess to the average when the fluctuations improve.
!   if(called_by_qa .and. (func.lt.func_sav)) then
!     uwdiff_av=sum(uwdiff(1:ndata))/ndata
!     write(6,'(''eguess being reset from'',f14.6,'' to'',f14.6,'' change='',es12.4)') eguess, eguess+uwdiff_av, uwdiff_av
!     eguess=eguess+uwdiff_av
!     uwdiff(1:ndata)=uwdiff(1:ndata)-uwdiff_av
!     diff(1:ndata)=diff(1:ndata)-uwdiff_av*eavri
!     func=sum(diff(1:ndata2)**2)
!     func_sav=func
!   endif

  endif

  return
end function func
