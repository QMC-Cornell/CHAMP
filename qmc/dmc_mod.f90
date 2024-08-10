module dmc_mod

  use all_tools_mod
  use control_mod
  use montecarlo_mod
  use average_mod
  use print_mod
  use restart_mod
  use allocations_mod
  use walkers_mod
  use projector, only : interface_projector

! Declaration of global variables and default values


  contains

! ==============================================================================
  subroutine dmc_init
! ------------------------------------------------------------------------------
! Description   : Initialize various quantities for Diffusion Monte Carlo and wave function optimization
! Description   : Used both for all-electon and 1-electron move versions.
! Description   : Taken from C. Umrigar
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  use branch_dmc_opt_mod
  implicit none

! local

! If add_diag(1) <= 0 turn OFF add_diag optimization and use fixed, add_diag
! equal to abs(add_diag) specified in input.  Usually you do NOT want to do this
! because optimization of add_diag adds only about 20% to the time since the
! number of blocks used for the correlated sampling is one-tenth of those used
! for calculating the gradient vector and hessian, overlap and hamiltonian matrices,
! and because optimization of add_diag makes the method totally stable.
  if(add_diag(1).gt.0) then
    iadd_diag_opt=1
  else
    add_diag(1)=abs(add_diag(1))
    iadd_diag_opt=0
  endif

  if(nforce.gt.1) then
! force parameters
    if (.not. use_parser) then
      call readforce
    endif
! parameters for secondary geometry wave function
    call wf_secondary
  else
    nwprod=1
    call object_modified ('nwprod')
    nwftype=1
    iwftype(1)=1
  endif
  call common_allocations

! allocations for DMC
  call object_provide ('nforce')
  call object_provide ('nelec')
  call object_provide ('ndetup')
  call object_provide ('ndetdn')
  call object_provide ('ndetupdn')
  call object_provide ('nup_square')
  call object_provide ('ndn_square')
  call object_provide ('nparm')

!  call alloc ('iage', iage, MWALK)
  call alloc ('xoldw', xoldw, 3, nelec, MWALK, nforce)
!  call alloc ('voldw', voldw, 3, nelec, MWALK, nforce)
!  call alloc ('voldw', voldw, 3, nelec, MWALK, nforce)
!  call alloc ('psidow', psidow, MWALK, nforce)
!  call alloc ('psijow', psijow, MWALK, nforce)
!  call alloc ('peow', peow, MWALK, nforce)
!  call alloc ('peiow', peiow, MWALK, nforce)
!  call alloc ('d2ow', d2ow, MWALK, nforce)
!
!  call alloc ('eoldw', eoldw, MWALK, nforce)
!  call alloc ('pwt', pwt, MWALK, nforce)
!  call alloc ('wt', wt, MWALK)
!
!  call alloc ('div_vow', div_vow, nelec, MWALK)
!  
!  call alloc ('slmuiw', slmuiw, nup_square, ndetup, MWALK)
!  call alloc ('slmdiw', slmdiw, ndn_square, ndetdn, MWALK)
!  call alloc ('fpuw', fpuw, 3, nup_square, ndetup, MWALK)
!  call alloc ('fpdw', fpdw, 3, ndn_square, ndetdn, MWALK)
!  call alloc ('detuw', detuw, ndetup, MWALK)
!  call alloc ('detdw', detdw, ndetdn, MWALK)
!  call alloc ('ddeti_detiw', ddeti_detiw, 3, nelec, ndetupdn, MWALK)
!  
!  call alloc ('fsow', fsow, nelec, nelec, MWALK)
!  call alloc ('fijow', fijow, 3, nelec, nelec, MWALK)
!  call alloc ('fsumow', fsumow, MWALK)
!  call alloc ('fjow', fjow, 3, nelec, MWALK)
!
!  call alloc ('denergy_old_dmc', denergy_old_dmc, nparm, MWALK)
!  call alloc ('wi_w', wi_w, nparm, MWALK)
!
!  call object_modified ('nfprod')
!  allocate(wtgen(0:nfprod))
!  allocate(ff(0:nfprod))
!
!  call object_modified ('nwprod')
!  allocate(wthist(MWALK,0:nwprod,nforce))
!
!  call alloc ('ajacold', ajacold, MWALK, nforce)
!  call alloc ('fratio', fratio, MWALK, nforce)
 
! get initial walkers from files or generate them
  if (l_mode_mpi) then
   call open_files_mpi
  else
   call open_files
  endif

  if(idmc.eq.3) call interface_projector(nctype, ncent, iwctype, cent, znuc, nelec, nup, tau, etrial)

  end subroutine dmc_init

! ==============================================================================
  subroutine dmc_alloc
! Created       : T. Anderson, 6 Feb 2022
! ------------------------------------------------------------------------------
  use all_modules_mod
  use branch_dmc_opt_mod
  use gamma_mod, only: noccup, noccdn
!  use contrldmc_mod, only: nfprod
!  use rundmc_mod, only: rundmc
!  use branch_mod, only: MWALK, nwalk
  implicit none

  call object_provide ('nforce')
  call object_provide ('nelec')
  call object_provide ('ndetup')
  call object_provide ('ndetdn')
  call object_provide ('ndetupdn')
  call object_provide ('nup_square')
  call object_provide ('ndn_square')
  call object_provide ('nparm')
  call object_provide ('noccup')
  call object_provide ('noccdn')

  call alloc ('iage', iage, MWALK)
  call alloc ('xoldw', xoldw, 3, nelec, MWALK, nforce)
  call alloc ('voldw', voldw, 3, nelec, MWALK, nforce)
  call alloc ('voldw', voldw, 3, nelec, MWALK, nforce)
  call alloc ('psidow', psidow, MWALK, nforce)
  call alloc ('psijow', psijow, MWALK, nforce)
  call alloc ('peow', peow, MWALK, nforce)
  call alloc ('peiow', peiow, MWALK, nforce)
  call alloc ('d2ow', d2ow, MWALK, nforce)

!  call alloc ('quadr', quadr,       nquad*ncent, nelec) !TA
!  call alloc ('quadx', quadx, ndim, nquad*ncent, nelec) !TA
  nullify(quadx,quadr) !TA

  call alloc ('eoldw', eoldw, MWALK, nforce)
  call alloc ('pwt', pwt, MWALK, nforce)
  call alloc ('wt', wt, MWALK)

  call alloc ('div_vow', div_vow, nelec, MWALK)

  call alloc ('slmuiw', slmuiw, nup_square, ndetup, MWALK)
  call alloc ('slmdiw', slmdiw, ndn_square, ndetdn, MWALK)
  call alloc ('fpuw', fpuw, 3, nup_square, ndetup, MWALK)
  call alloc ('fpdw', fpdw, 3, ndn_square, ndetdn, MWALK)
  call alloc ('detuw', detuw, ndetup, MWALK)
  call alloc ('detdw', detdw, ndetdn, MWALK)
  call alloc ('ddeti_detiw', ddeti_detiw, 3, nelec, ndetupdn, MWALK)
  
  call alloc ('fsow', fsow, nelec, nelec, MWALK)
  call alloc ('fijow', fijow, 3, nelec, nelec, MWALK)
  call alloc ('fsumow', fsumow, MWALK)
  call alloc ('fjow', fjow, 3, nelec, MWALK)
  call alloc ('lapjijow', lapjijow, nelec, nelec, MWALK)
  call alloc ('lapjow', lapjow, nelec, MWALK)

  call alloc ('denergy_old_dmc', denergy_old_dmc, nparm, MWALK)
  call alloc ('wi_w', wi_w, nparm, MWALK)

  call object_provide ('nfprod')
  call alloc_range ('wtgen', wtgen, 0, nfprod)
  call alloc_range ('ff', ff, 0, nfprod-1)

  call object_provide ('nwprod')
  call alloc_range ('wthist', wthist, 1, MWALK, 0, nwprod, 1, nforce)

  call alloc ('ajacold', ajacold, MWALK, nforce)
  call alloc ('fratio', fratio, MWALK, nforce)

  call alloc ('pot_ee', pot_ee, nelec) ! not sure why this is needed, since we
    ! already did this in dmc_init (common_allocations), but we get an error if we
    ! don't.  ACM 8/27/11
  call alloc ('pot_ee_old', pot_ee_old, nelec)
  call alloc ('pot_ee_new', pot_ee_new, nelec)
  call alloc ('pot_ee_oldw', pot_ee_oldw, nelec, MWALK, nforce)

  !Allocate walker arrays used for fast derivatives TA
  call alloc('orbw', orbw, nelec, orb_tot_nb, MWALK)
  call alloc('dorbw', dorbw, ndim, nelec, orb_tot_nb, MWALK)
  call alloc('ddorbw', ddorbw, nelec, orb_tot_nb, MWALK)
  call alloc('aiupw', aiupw, nup, nup, MWALK)
  call alloc('aidnw', aidnw, ndn, ndn, MWALK)
  call alloc('deta_upw', deta_upw, MWALK)
  call alloc('deta_dnw', deta_dnw, MWALK)
  call alloc('tupw', tupw, nup,orb_tot_nb,MWALK)
  call alloc('tdnw', tdnw, ndn,orb_tot_nb,MWALK)
  call alloc('yupw', yupw, noccup,nup,MWALK)
  call alloc('ydnw', ydnw, noccdn,ndn,MWALK)
  call alloc('chiw', chiw, MWALK)
  call alloc('invupw', invupw, nup*nup,ndetup,MWALK)
  call alloc('invdnw', invdnw, ndn*ndn,ndetdn,MWALK)
  call alloc('detupw', detupw, ndetup,MWALK)
  call alloc('detdnw', detdnw, ndetdn,MWALK)
  end subroutine dmc_alloc

! ==============================================================================
  subroutine dmc
! ------------------------------------------------------------------------------
! Description   : routine cleaned from original dmc routine
! written by Cyrus Umrigar with major contributions by Claudia Filippi.
! Uses the diffusion Monte Carlo algorithm described in:
! A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
!    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993).
!
! common block variables:
!
! /const/
!   nelec       = number of electrons
!   pi          = 3.14159...
!   hb          = hbar**2/(2m)
!   delta       = side of box in which metropolis steps are made
!   deltai      = 1/delta
!   fbias       = force bias parameter
! /contrl/
!   nstep       = number of metropolis steps/block
!   nblk        = number of blocks of nstep steps after the equilibrium blocks
!   nblkeq      = number of equilibrium blocks
!   nconf       = initial and target number of dmc configurations per processor
!   nconf_global= nconf*nproc in dmc_mov1_mpi2 and dmc_mov1_mpi3 modes
!   nconf_new   = number of new configurations saved per processor
!   idump       = 1 dump out stuff for a restart
!   irstar      = 1 pick up stuff for a restart
! /contrldmc/ tau,rttau,taueff(MFORCE),tautot,nfprod,idmc,ipq
!       ,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
!   tau         = time-step
!   nfprod      = number of f's used for removing finite popul. bias
! Control variables are:
! idmc         < 0     VMC
!              > 0     DMC
! abs(idmc)    = 1 **  simple kernel using dmc.brock.f
!              = 2     good kernel using dmc_good or dmc_good_inhom
! ipq         <= 0 *   do not use expected averages
!             >= 1     use expected averages (mostly in all-electron move algorithm)
! itau_eff    <=-1 *   always use tau in branching (not implemented)
!              = 0     use 0 / tau for acc /acc_int moves in branching
!             >= 1     use tau_eff (calcul in equilibration runs) for all moves
! iacc_rej    <=-1 **  accept all moves (except possibly node crossings)
!              = 0 **  use weights rather than accept/reject
!             >= 1     use accept/reject
! icross      <=-1 **  kill walkers that cross nodes (not implemented)
!              = 0     reject walkers that cross nodes
!             >= 1     allow walkers to cross nodes
!                      (OK since crossing prob. goes as tau^(3/2))
! icuspg      <= 0     approximate cusp in Green function
!             >= 1     impose correct cusp condition on Green function
! idiv_v      <= 0     do not use div_v to modify diffusion in good kernel
!              = 1     use homog. div_v to modify diffusion in good kernel
!             >= 2     use inhomog. div_v to modify diffusion in good kernel
! icut_br     <= 0     do not limit branching
!             >= 1 *   use smooth formulae to limit branching to (1/2,2)
!                      (bad because it makes energies depend on E_trial)
! icut_e      <= 0     do not limit energy
!             >= 1 *   use smooth formulae to limit energy (not implemented)
!
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  use branch_dmc_opt_mod
  use gamma_mod, only: noccup, noccdn
  use walker_type_mod
  use stats_type_mod
  use contrldmc_mod, only: nfprod
  use rundmc_mod, only: rundmc
  use branch_mod, only: MWALK, nwalk
  use fragments_mod, only: nfrag, etrialfrag
  use stats_index_mod
  implicit none

! local
  integer i, istep, jj, k, ifrag, IERROR
  integer ngfmc
  character (len=27) fmt
  integer iwalk
  real(dp) passes, efin
  integer iw
  real(dp) pexp, undo_pop_ctrl, dff
  type(stats_t) :: stats
  type(walker_t), target :: walkers(MWALK)
  integer :: nest_tot
  real(dp) :: eestfrag(nfrag+1)

! initial printing
  write(6,*)
  write(6,'(a)') '*********** START DMC CALCULATION  ***********'

! request average of local energy
  call object_average_request ('eloc_av')

! request variance on averaged local energy
  call object_variance_request ('eloc_av_var')

! sigma
  call object_average_request ('eloc_sq_av')
  call object_error_request ('error_sigma')

  call print_list_of_averages_and_errors

! allocations (must be here because nforce can changed)
  call dmc_alloc

! initialize sums and averages
  if(irstar /= 1) then
     call zeres0_dmc
  endif

! forces implemented only for certain dmc control options
  if(nforce > 1) write(6,'(''Possible Warning: force implemented for certain dmc control options'')')

! If nconf_new > 0 then we want to write nconf_new configurations from each processor for a future
! optimization or dmc calculation. So figure out how often we need to write a
! configuration to produce nconf_new configurations. If nconf_new = 0
! then set up so no configurations are written.
!  if(nconf_new.eq.0) then
!     ngfmc=2*nstep*nblk
!  else
!     ngfmc=max(1,(nstep*nblk)*nconf/nconf_new)
!  endif

  eest=0
  if (l_fragments) eestfrag=0
  do iw=1,MWALK
    if (iw.LE.nwalk) then
      call walkers(iw)%init(xoldw(:,:,iw,1),wt(iw))
      eest=eest+walkers(iw)%psi%eloc
      if (l_fragments) eestfrag=eestfrag+walkers(iw)%psi%enefrag
    else
      call walkers(iw)%init
    endif
  enddo
  eest=eest/nwalk
  if (l_fragments) eestfrag=eestfrag/nwalk

#if defined(MPI)
  call MPI_Allreduce(MPI_IN_PLACE,eest,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
  eest=eest/nproc
#endif

  ff=dexp((etrial-eest)*tau)
  fprod=product(ff)
  pexp=exp(-1d0/nfprod)
  undo_pop_ctrl=dexp((eest-etrial)*tau/(1-pexp))

  nest_tot=nest_static
  if (l_fragments) then
    nest_tot=nest_tot+2*(nfrag+1)
    offset_enefrag_=nest_static
    offset_taufrag_=nest_static+nfrag+1
  endif
  call stats%init(nest_tot)

  call print_cpu_time_in_seconds ('Beginning of equilibration')
  call systemflush(6)

  !equilibration run 1
  call rundmc(nblkeq, nstep, &
              tau, etrial, eest, eestfrag, undo_pop_ctrl, &
              walkers, nwalk, stats, .TRUE.)

  call print_cpu_time_in_seconds ('End        of equilibration')

  if (l_reset_etrial) then 
      if (l_fragments) then
        write(6,'('''')')
        do ifrag=1,nfrag
          ff=ff*dexp((eestfrag(ifrag)-etrialfrag(ifrag))*stats%avg(offset_taufrag_+ifrag))
          write(6,'(''etrial of fragment '',I2,'' changed from'',f11.6,'' to'',f11.6)') ifrag, etrialfrag(ifrag), eestfrag(ifrag)
          etrialfrag=eestfrag
        enddo
        write(6,'('''')')
      else
        ff=ff*dexp((eest-etrial)*stats%avg(taueff_))
        write(6,'(/,''etrial changed from'',f11.6,'' to'',f11.6,/)') etrial, eest
        etrial=eest
      endif
      !ff=ff*dexp((eest-etrial)*stats%avg(taueff_))
      !write(6,'(/,''etrial changed from'',f11.6,'' to'',f11.6,/)') etrial, eest
      !etrial=eest
  endif
  call stats%clear

  !equilibration run 2
  call rundmc(nblkeq, nstep, &
              tau, etrial, eest, eestfrag, undo_pop_ctrl, &
              walkers, nwalk, stats, .TRUE.)

  call print_cpu_time_in_seconds ('End        of equilibration')

  if (l_reset_etrial) then 
      if (l_fragments) then
        do ifrag=1,nfrag
          ff=ff*dexp((eestfrag(ifrag)-etrialfrag(ifrag))*stats%avg(offset_taufrag_+ifrag))
          write(6,'(/,''etrial of fragment '',I2,'' changed from'',f11.6,'' to'',f11.6,/)') ifrag, etrialfrag(ifrag), eestfrag(ifrag)
          etrialfrag=eestfrag
        enddo
      else
        ff=ff*dexp((eest-etrial)*stats%avg(taueff_))
        write(6,'(/,''etrial changed from'',f11.6,'' to'',f11.6,/)') etrial, eest
        etrial=eest
      endif
      !ff=ff*dexp((eest-etrial)*stats%avg(taueff_))
      !write(6,'(/,''etrial changed from'',f11.6,'' to'',f11.6,/)') etrial, eest
      !etrial=eest
  endif
  call stats%clear

  call print_cpu_time_in_seconds ('Beginning of  accumulation')

  !main dmc run
  call rundmc(nblk, nstep, &
              tau, etrial, eest, eestfrag, undo_pop_ctrl, &
              walkers, nwalk, stats, .FALSE.)

  call print_cpu_time_in_seconds ('End       of  accumulation')

  write (6,*)
  write(6,'(''passes,wgcum='',9es12.4)') dfloat(nblk*nstep),stats%wgt_tsav(3)
  call final_write(nstep, nblk, etrial, walkers, nwalk, stats)
  end subroutine dmc

! ==============================================================================
  subroutine final_write(nstep, nblock, etrial, walkers, nwalk, stats)

! Created       : T. Anderson, Jul 2023
! ------------------------------------------------------------------------------
      use walker_type_mod, only: walker_t
      use stats_type_mod, only: stats_t
      use contrldmc_mod, only: idmc
      use age_mod, only: ioldest, ioldestmx
      use contrl_per_mod, only: iperiodic
      use contr3_mod, only: mode
      use contrl_mod, only: nblkeq, nconf, nconf_global
      use atom_mod, only: ncent
      use dim_mod, only: ndim
      use const_mod, only: nelec
      use contrldmc_mod, only: nfprod, tau
      use stats_mod, only: nodecr, try_int, dr2ac, dr2un, acc, acc_int
      use stats_index_mod
      implicit none
      integer           :: nstep
      integer           :: nblock
      real(dp)          :: etrial
      type(walker_t)    :: walkers(:)
      integer           :: nwalk
      type(stats_t)     :: stats

      character*80  :: fmt
      integer       :: i, j, k
      integer       :: ipr
      real(dp)      :: passes
      real(dp)      :: eval
      integer       :: nblock_proc
      real(dp)      :: pass_proc
      real(dp)      :: eval_proc
      real(dp)      :: eval0_proc_eff
      real(dp)      :: evalf_proc_eff
      real(dp)      :: evalg_proc_eff
      real(dp)      :: rtpass_proc1
      real(dp)      :: rteval0_proc_eff1
      real(dp)      :: rtevalf_proc_eff1
      real(dp)      :: rtevalg_proc_eff1
      real(dp)      :: avg(stats%nest), sig(stats%nest)
      real(dp)      :: err(stats%nest), err1(stats%nest)
      real(dp)      :: sgtc(stats%nest), tcorr(stats%nest)
      real(dp)      :: wgcum

      call stats % calculate(avg, err, err1)

      passes=dfloat(nblock)*dfloat(nstep)
      eval=nconf*passes
      if (l_mode_dmc_mov1_mpi2) then
          nblock_proc=nblock !in mpi2 mode, all processes contribute to a single step
      else
          nblock_proc=nblock*nproc
      endif
      pass_proc=dfloat(nblock_proc)*dfloat(nstep)
      eval_proc=nconf_global*pass_proc
      wgcum =stats % wgt_tsav(elocalg_) !db
      eval0_proc_eff=dfloat(nconf_global)*nstep*stats % nblck_eff(elocal0_) !db
      evalf_proc_eff=dfloat(nconf_global)*nstep*stats % nblck_eff(elocalf_) !db
      evalg_proc_eff=dfloat(nconf_global)*nstep*stats % nblck_eff(elocalg_) !db
      rtpass_proc1=dsqrt(pass_proc-1)
      rteval0_proc_eff1=dsqrt(eval0_proc_eff-1)
      rtevalf_proc_eff1=dsqrt(evalf_proc_eff-1)
      rtevalg_proc_eff1=dsqrt(evalg_proc_eff-1)

      sig=err1*rtevalg_proc_eff1
      sig(elocal0_)=err1(elocal0_)*rteval0_proc_eff1
      sig(elocalf_)=err1(elocalf_)*rtevalf_proc_eff1
      sgtc=err*rtevalg_proc_eff1
      sgtc(elocal0_)=err(elocal0_)*rteval0_proc_eff1
      sgtc(elocalf_)=err(elocalf_)*rtevalf_proc_eff1
      err1=max(err1,tiny(0d0))
      tcorr=(err/err1)**2

      write(6,'(/,''Final write:'')')

      write(6,'(''passes, eval, pass_proc, eval_proc, eval_proc_eff, evalf_proc_eff, evalg_proc_eff'',19f11.0)') &
     & passes,eval,pass_proc,eval_proc,eval0_proc_eff,evalf_proc_eff,evalg_proc_eff

      if(ipr.gt.-2) write(11,'(3i5,f11.5,f7.4,f10.7,'' nstep,nblk,nconf_global,etrial,tau,taueff'')') &
    & nstep,nblock,nconf_global,etrial,tau,avg(taueff_)

!TODO: Implement radial probability histogram
!      if(print_radial_probability .and. iperiodic.eq.0 .and. ncent.eq.1 .and. ipr.ge.-4) then
!        if(ndim.eq.3) write(6,'(/,'' r  4*pi*r^2*rho 4*pi*r^2*rhoup 4*pi*r^2*rhodn'')')
!        if(ndim.eq.2) write(6,'(/,'' r  2*pi*r*rho 2*pi*r*rhoup 2*pi*r*rhodn'')')
!        delr=one/delri
!        term=one/(wgcum(1)*delr)
!        !TODO: calculate rprob, etc.
!        do 18 i=1,NRAD
!   18     write(6,'(f5.3,3f10.6)') delr*(i-half),rprob(i)*term,rprobup(i)*term,rprobdn(i)*term
!      endif

      if(idmc.ge.0) then
        write(6,'(/,''ages of walkers are:'')')
        write(6,'(10i4)') (walkers(i)%age,i=1,nwalk)
        do 19 i=1,nwalk
          if(walkers(i)%age.gt.50) then
            write(6,'(i4,i6,f10.4,99f8.4)') i,walkers(i)%age,walkers%psi%eloc,((walkers(i)%x(k,j),k=1,ndim),j=1,nelec)
            write(6,'(99f8.4)') ((walkers(i)%psi%grad(k,j),k=1,ndim),j=1,nelec)
          endif
   19   continue

        write(6,'(''age of oldest walker (this generation, any gen)='',3i9)') ioldest,ioldestmx
      endif

#if defined(MPI)
      call MPI_Allreduce(MPI_IN_PLACE,nodecr, 1,MPI_INTEGER,         MPI_SUM,MPI_COMM_WORLD,I)
      call MPI_Allreduce(MPI_IN_PLACE,try_int,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,I)
      call MPI_Allreduce(MPI_IN_PLACE,dr2ac,  1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,I)
      call MPI_Allreduce(MPI_IN_PLACE,dr2un,  1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,I)
      call MPI_Allreduce(MPI_IN_PLACE,acc,    1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,I)
      call MPI_Allreduce(MPI_IN_PLACE,acc_int,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,I)
#endif

      write(6,'(a,f10.5)') 'average of the squares drift-dif moves for accepted steps = ',dr2ac/try_int
      write(6,'(a,f10.5)') 'average of the squares drift-dif moves for all      steps = ',dr2un/try_int

      write(6,'(''taueff'',20f7.4)') avg(taueff_)

      write(fmt,'(''(/,a16,2x,a'',i3,'')'')') len_trim(title)
      write(6,fmt) mode,title
      write(6,'(''No/frac. of node crossings,acceptance='',i9,3f10.6)') nodecr,dfloat(nodecr)/try_int,acc/try_int,acc_int/try_int
      if(idmc.lt.0.and.(acc/try_int).lt.0.3d0) write(6,'(''Warning: Low acceptance, reduce time-step tau'')')
      if(idmc.ge.0.and.(acc/try_int).lt.0.7d0) write(6,'(''Warning: Low acceptance, reduce time-step tau'')')

      if(idmc.ge.0) then
        write(6,'(''No. of walkers at end of run='',i5)') nwalk
        write(6,'(''nwalk_eff/nwalk         ='',2f6.3)') stats%nstep_eff(elocal0_)/pass_proc,stats%nblck_eff(elocal0_)/nblock_proc
        write(6,'(''nwalk_eff/nwalk with f  ='',2f6.3)') stats%nstep_eff(elocalf_)/pass_proc,stats%nblck_eff(elocalf_)/nblock_proc
        write(6,'(''nwalk_eff/nwalk with fs ='',2f6.3)') stats%nstep_eff(elocalg_)/pass_proc,stats%nblck_eff(elocalg_)/nblock_proc
      endif

      write(6,'(''nconf*nproc*passes'',t20,''nconf nproc     passes nstep  nblk nblkeq    tau   taueff'' &
     &,/,f18.0,2i6,f11.0,3i6,2f9.5)') eval*nproc,nconf,nproc,passes,nstep,nblock,nblkeq,tau,avg(taueff_)
      write(6,'(''physical variable         average     rms error   sigma*T_cor  sigma   T_cor'')')
      if(idmc.ge.0) then
        write(6,'(''weights ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
     &  avg(ws0_),err(ws0_),err(ws0_)*rtpass_proc1,err1(ws0_)*rtpass_proc1,(err(ws0_)/err1(ws0_))**2
        write(6,'(''wts with f ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
     &  avg(wsf_),err(wsf_),err(wsf_)*rtpass_proc1,err1(wsf_)*rtpass_proc1,(err(wsf_)/err1(wsf_))**2
        write(6,'(''wts with fs ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
     &  avg(wsg_),err(wsg_),err(wsg_)*rtpass_proc1,err1(wsg_)*rtpass_proc1,(err(wsg_)/err1(wsg_))**2
        write(6,'(a,2f10.8)') 'approx. normalized overlap of FN and trial wave functions= ',wgcum/dsqrt(try_int*stats%est_tsav(wtsq_)/nelec)
        write(6,'(a,f10.8)') 'unnormalized overlap of FN and trial wave functions= ', wgcum/(try_int/nelec)

! Mixed energy estimators
        write(6,'(''total energy (   0) ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
      & avg(elocal0_), err(elocal0_), sgtc(elocal0_), sig(elocal0_), tcorr(elocal0_)
        write(6,'(''total energy (   1) ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
      & avg(elocalf_), err(elocalf_), sgtc(elocalf_), sig(elocalf_), tcorr(elocalf_)
      endif

      write(6,'(''total energy ('',i4,'') ='',t22,f14.7,'' +-'',f11.7,2f10.5    ,f8.2)') &
    & nfprod, avg(elocalg_), err(elocalg_), sgtc(elocalg_), sig(elocalg_), tcorr(elocalg_)

      write(6,'(''potential energy ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)')   & 
    & avg(epot_), err(epot_), sgtc(epot_)
      write(6,'(''interaction energy ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
    & avg(epot_ee_), err(epot_ee_), sgtc(epot_ee_)
      write(6,'(''jf kinetic energy ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)')  & 
    & avg(ekinjf_), err(ekinjf_), sgtc(ekinjf_)
      write(6,'(''pb kinetic energy ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)')  & 
    & avg(ekinpb_), err(ekinpb_), sgtc(ekinpb_)

      if(iperiodic.eq.0 .and. ncent.eq.1) then
        write(6,'(''<r>_av ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') &
    &   avg(R1_), err(R1_), sgtc(R1_)
        write(6,'(''<r2>_av ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') & 
    &   avg(R2_), err(R2_), sgtc(R2_)
        write(6,'(''<r3>_av ='',t22,f14.6,'' +-'',f11.6,2f10.4,f8.2)') & 
    &   avg(R3_), err(R3_), sgtc(R3_)
        write(6,'(''<r4>_av ='',t22,f14.5,'' +-'',f11.5,2f10.3,f8.2)') & 
    &   avg(R4_), err(R4_), sgtc(R4_)
        write(6,'(''<ri>_av ='',t22,f14.7,'' +-'',f11.7,2f10.5,f8.2)') & 
    &   avg(RI_), err(RI_), sgtc(RI_)
      endif

    end subroutine final_write

! ==============================================================================
  subroutine dmc_old
! ------------------------------------------------------------------------------
! Description   : routine cleaned from original dmc routine
! written by Cyrus Umrigar with major contributions by Claudia Filippi.
! Uses the diffusion Monte Carlo algorithm described in:
! A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
!    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993).
!
! common block variables:
!
! /const/
!   nelec       = number of electrons
!   pi          = 3.14159...
!   hb          = hbar**2/(2m)
!   delta       = side of box in which metropolis steps are made
!   deltai      = 1/delta
!   fbias       = force bias parameter
! /contrl/
!   nstep       = number of metropolis steps/block
!   nblk        = number of blocks of nstep steps after the equilibrium blocks
!   nblkeq      = number of equilibrium blocks
!   nconf       = initial and target number of dmc configurations per processor
!   nconf_global= nconf*nproc in dmc_mov1_mpi2 and dmc_mov1_mpi3 modes
!   nconf_new   = number of new configurations saved per processor
!   idump       = 1 dump out stuff for a restart
!   irstar      = 1 pick up stuff for a restart
! /contrldmc/ tau,rttau,taueff(MFORCE),tautot,nfprod,idmc,ipq
!       ,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
!   tau         = time-step
!   nfprod      = number of f's used for removing finite popul. bias
! Control variables are:
! idmc         < 0     VMC
!              > 0     DMC
! abs(idmc)    = 1 **  simple kernel using dmc.brock.f
!              = 2     good kernel using dmc_good or dmc_good_inhom
! ipq         <= 0 *   do not use expected averages
!             >= 1     use expected averages (mostly in all-electron move algorithm)
! itau_eff    <=-1 *   always use tau in branching (not implemented)
!              = 0     use 0 / tau for acc /acc_int moves in branching
!             >= 1     use tau_eff (calcul in equilibration runs) for all moves
! iacc_rej    <=-1 **  accept all moves (except possibly node crossings)
!              = 0 **  use weights rather than accept/reject
!             >= 1     use accept/reject
! icross      <=-1 **  kill walkers that cross nodes (not implemented)
!              = 0     reject walkers that cross nodes
!             >= 1     allow walkers to cross nodes
!                      (OK since crossing prob. goes as tau^(3/2))
! icuspg      <= 0     approximate cusp in Green function
!             >= 1     impose correct cusp condition on Green function
! idiv_v      <= 0     do not use div_v to modify diffusion in good kernel
!              = 1     use homog. div_v to modify diffusion in good kernel
!             >= 2     use inhomog. div_v to modify diffusion in good kernel
! icut_br     <= 0     do not limit branching
!             >= 1 *   use smooth formulae to limit branching to (1/2,2)
!                      (bad because it makes energies depend on E_trial)
! icut_e      <= 0     do not limit energy
!             >= 1 *   use smooth formulae to limit energy (not implemented)
!
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  use branch_dmc_opt_mod
  use gamma_mod, only: noccup, noccdn
  implicit none

! local
  integer i, istep, jj, k
  integer ngfmc
  character (len=27) fmt
  integer iwalk
  real(dp) passes, efin

! initial printing
  write(6,*)
  write(6,'(a)') '*********** START DMC CALCULATION  ***********'

! request average of local energy
  call object_average_request ('eloc_av')

! request variance on averaged local energy
  call object_variance_request ('eloc_av_var')

! sigma
  call object_average_request ('eloc_sq_av')
  call object_error_request ('error_sigma')

  call print_list_of_averages_and_errors


! allocations (must be here because nforce can changed)
  call object_provide ('nforce')
  call object_provide ('nelec')
  call object_provide ('ndetup')
  call object_provide ('ndetdn')
  call object_provide ('ndetupdn')
  call object_provide ('nup_square')
  call object_provide ('ndn_square')
  call object_provide ('nparm')
  call object_provide ('noccup')
  call object_provide ('noccdn')

  call alloc ('iage', iage, MWALK)
  call alloc ('xoldw', xoldw, 3, nelec, MWALK, nforce)
  call alloc ('voldw', voldw, 3, nelec, MWALK, nforce)
  call alloc ('voldw', voldw, 3, nelec, MWALK, nforce)
  call alloc ('psidow', psidow, MWALK, nforce)
  call alloc ('psijow', psijow, MWALK, nforce)
  call alloc ('peow', peow, MWALK, nforce)
  call alloc ('peiow', peiow, MWALK, nforce)
  call alloc ('d2ow', d2ow, MWALK, nforce)

!  call alloc ('quadr', quadr,       nquad*ncent, nelec, MWALK) !TA
!  call alloc ('quadx', quadx, ndim, nquad*ncent, nelec, MWALK) !TA

  call alloc ('eoldw', eoldw, MWALK, nforce)
  call alloc ('pwt', pwt, MWALK, nforce)
  call alloc ('wt', wt, MWALK)

  call alloc ('div_vow', div_vow, nelec, MWALK)

  call alloc ('slmuiw', slmuiw, nup_square, ndetup, MWALK)
  call alloc ('slmdiw', slmdiw, ndn_square, ndetdn, MWALK)
  call alloc ('fpuw', fpuw, 3, nup_square, ndetup, MWALK)
  call alloc ('fpdw', fpdw, 3, ndn_square, ndetdn, MWALK)
  call alloc ('detuw', detuw, ndetup, MWALK)
  call alloc ('detdw', detdw, ndetdn, MWALK)
  call alloc ('ddeti_detiw', ddeti_detiw, 3, nelec, ndetupdn, MWALK)
  
  call alloc ('fsow', fsow, nelec, nelec, MWALK)
  call alloc ('fijow', fijow, 3, nelec, nelec, MWALK)
  call alloc ('fsumow', fsumow, MWALK)
  call alloc ('fjow', fjow, 3, nelec, MWALK)

  call alloc ('denergy_old_dmc', denergy_old_dmc, nparm, MWALK)
  call alloc ('wi_w', wi_w, nparm, MWALK)

  call object_provide ('nfprod')
  call alloc_range ('wtgen', wtgen, 0, nfprod)
  call alloc_range ('ff', ff, 0, nfprod)

  call object_provide ('nwprod')
  call alloc_range ('wthist', wthist, 1, MWALK, 0, nwprod, 1, nforce)

  call alloc ('ajacold', ajacold, MWALK, nforce)
  call alloc ('fratio', fratio, MWALK, nforce)

  call alloc ('pot_ee', pot_ee, nelec) ! not sure why this is needed, since we
    ! already did this in dmc_init (common_allocations), but we get an error if we
    ! don't.  ACM 8/27/11
  call alloc ('pot_ee_old', pot_ee_old, nelec)
  call alloc ('pot_ee_new', pot_ee_new, nelec)
  call alloc ('pot_ee_oldw', pot_ee_oldw, nelec, MWALK, nforce)

  !Allocate walker arrays used for fast derivatives TA
  call alloc('orbw', orbw, nelec, orb_tot_nb, MWALK)
  call alloc('dorbw', dorbw, ndim, nelec, orb_tot_nb, MWALK)
  call alloc('aiupw', aiupw, nup, nup, MWALK)
  call alloc('aidnw', aidnw, ndn, ndn, MWALK)
  call alloc('deta_upw', deta_upw, MWALK)
  call alloc('deta_dnw', deta_dnw, MWALK)
  call alloc('tupw', tupw, nup,orb_tot_nb,MWALK)
  call alloc('tdnw', tdnw, ndn,orb_tot_nb,MWALK)
  call alloc('yupw', yupw, noccup,nup,MWALK)
  call alloc('ydnw', ydnw, noccdn,ndn,MWALK)
  call alloc('chiw', chiw, MWALK)
  call alloc('invupw', invupw, nup*nup,ndetup,MWALK)
  call alloc('invdnw', invdnw, ndn*ndn,ndetdn,MWALK)
  call alloc('detupw', detupw, ndetup,MWALK)
  call alloc('detdnw', detdnw, ndetdn,MWALK)

! initialize sums and averages
  if(irstar /= 1) then
     call zeres0_dmc
  endif

! forces implemented only for certain dmc control options
  if(nforce > 1) write(6,'(''Possible Warning: force implemented for certain dmc control options'')')

! If nconf_new > 0 then we want to write nconf_new configurations from each processor for a future
! optimization or dmc calculation. So figure out how often we need to write a
! configuration to produce nconf_new configurations. If nconf_new = 0
! then set up so no configurations are written.
  if(nconf_new.eq.0) then
     ngfmc=2*nstep*nblk
  else
     ngfmc=max(1,(nstep*nblk)*nconf/nconf_new)
  endif

  call print_cpu_time_in_seconds ('Beginning of equilibration')
  call systemflush(6)

! Equilibration phase
  l_equilibration = .true.
  do i = 1, 2*nblkeq

     if ((i == nblkeq+1) .and. irstar /= 1) then
        call print_cpu_time_in_seconds ('End        of equilibration')
        call zerest_dmc
     endif

     do istep=1,nstep

        ipass=ipass+1
        call dmc_algorithm

        call stochastic_reconfiguration2

!       averages at each step
        call acues1_dmc

!       walkers reconfiguration
        if(l_branching) call splitj

      enddo ! istep

!   averages at each block
    call acuest_dmc

  enddo ! end of equilibration
  l_equilibration = .false.

  call print_cpu_time_in_seconds ('End       of equilibration')
  write (6,*)

! reinitialization of averages
  call zerest_dmc

! End of equilibration

  call print_cpu_time_in_seconds ('Beginning of  accumulation')

! Actual DMC calculation
! loop over blocks
  do i=1,100000000
    block_iterations_nb = block_iterations_nb + 1
    do istep=1,nstep
      step_iterations_nb = step_iterations_nb + 1
      ipass=ipass+1
      call dmc_algorithm

      call object_modified_by_index (xoldw_index)

      call stochastic_reconfiguration2

!     compute averages
      call acues1_dmc
      call compute_averages_walk_step

!     walkers reconfiguration
      if(l_branching) call splitj

! Write out configuration for optimization/dmc/gfmc here
! We would like to:
! Reduce each electron to central simulation cell before writing position.
! Warning: This may result in the sign of the wavefunction being wrong depending
! on the k-pt and the number of cells moved.  Test done in si_cub_.5.0.0 shows this
! gives wrong energies and moves, possibly because the velocity is no longer consistent
! with the move, though I would have thought the velocity would be OK since both the
! wavfn and its gradients would change sign simultaneously.  May be I need to move
! in only even multiples of the sim. lattice.
          if(nconf_new.gt.0 .and. (mod(ipass,ngfmc).eq.1 .or. ngfmc.eq.1)) then
            if(ndim*nelec.lt.10) then
              write(fmt,'(a1,i1,a22)')'(',ndim*nelec,'f14.8,i3,es12.4,f13.5)'
             elseif(ndim*nelec.lt.100) then
              write(fmt,'(a1,i2,a22)')'(',ndim*nelec,'f14.8,i3,es12.4,f13.5)'
             elseif(ndim*nelec.lt.1000) then
              write(fmt,'(a1,i3,a22)')'(',ndim*nelec,'f14.8,i3,es12.4,f13.5)'
             elseif(ndim*nelec.lt.10000) then
              write(fmt,'(a1,i4,a22)')'(',ndim*nelec,'f14.8,i3,es12.4,f13.5)'
             else
              stop 'fmt not defined for ndim*nelec>10000'
            endif

            do 352 iwalk=1,nwalk
              if(iperiodic.ne.0) then
                do 351 jj=1,nelec
! 351             call reduce_sim_cell(xoldw(1,jj,iwalk,1),rlatt_sim,rlatt_sim_inv)
  351             call reduce_sim_cell(xoldw(1,jj,iwalk,1))
              endif
  352         write(7,fmt) ((xoldw(k,jj,iwalk,1),k=1,ndim),jj=1,nelec),int(sign(1.d0,psidow(iwalk,1))),log(dabs(psidow(iwalk,1)))+psijow(iwalk,1),eoldw(iwalk,1)
          endif

!       write out configurations in Scemama's format
        if (l_write_walkers) then
         if(mod(ipass,write_walkers_step) == 0) then
          do iwalk=1,nwalk
           do jj = 1, nelec
            if (l_write_walkers_modified_format) then
             write(file_walkers_out_unit,'(3(F11.7,X))') (xoldw(k,jj,iwalk,1),k=1,ndim)
            else
             write(file_walkers_out_unit,'(3(F7.3,X))') (xoldw(k,jj,iwalk,1),k=1,ndim)
            endif
           enddo
           write(file_walkers_out_unit,*) dexp(psijow(iwalk,1))*psidow(iwalk,1)
           if (l_write_walkers_modified_format) then
            write(file_walkers_out_unit,*) eoldw(iwalk,1)
           endif
          enddo
         endif
        endif

        enddo ! istep

        call acuest_dmc
        call compute_averages_block ! new 
        call compute_averages_walk_block
        call compute_covariances
        call compute_variances
        call compute_errors

!       write at each block
        call objects_print_at_each_block
        call routines_write_block

!        l_end_of_block = .false.

!     exit loop if nblk and threshold on statistical error reached
! In dmc_mov1_mpi2/3 modes egerr is always 0 since slaves never calculate it.
!     if (i .ge. nblk .and. egerr .gt. 0 .and. egerr .le. error_threshold) then
!     egerr .ne. egerr is true if egerr is NaN
      if (i .ge. nblk .and. ((egerr .le. error_threshold) .or. (egerr .ne. egerr))) then
        exit
      endif

      enddo ! i

      call print_cpu_time_in_seconds ('End       of  accumulation')
      write (6,*)

!     reinitilization at the end of MC iterations
      block_nb = block_iterations_nb  !save final block number
      step_iterations_nb  = 0
      block_iterations_nb = 0
      call reinit_averages_and_errors
      call reinit_routines_write_block

      passes=dble(iblk)*dble(nstep)
!     efin=ecum1/passes
      efin=egcum(1)/wgcum(1)

!     if(igradhess.ge.1) call grad_hess_jas_fin(passes,efin)
      if(igradhess.ge.1) call grad_hess_jas_fin(wgcum,efin)
      write(6,'(''passes,wgcum='',9es12.4)') passes,(wgcum(i),i=1,nforce)

!     Final writing
      call finwrt_dmc

      if(idump.eq.1) then
        if(idmc.lt.0) then
          if(index(mode,'mov1').ne.0) then
            open(10,file='restart_dmcvmc_mov1',status='unknown',form='unformatted')
           else
            open(10,file='restart_dmcvmc',status='unknown',form='unformatted')
          endif
         else
          if(index(mode,'mov1').ne.0) then
            open(10,file='restart_dmc_mov1',status='unknown',form='unformatted')
           else
            open(10,file='restart_dmc',status='unknown',form='unformatted')
          endif
        endif
        rewind 10
        call dumper_dmc
        close (unit=10)
      endif
      close (unit=9)

  end subroutine dmc_old

! ==============================================================================
  subroutine dmc_algorithm
! ------------------------------------------------------------------------------
! Description   : Choice of DMC algoritm
!
! Created       : J. Toulouse, 27 Mar 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character (len=max_string_len_rout), save :: lhere = 'dmc_algorithm'

! begin
  if(nloc.gt.0) call rotqua

  if (iabs(idmc) == 1) then
     call dmc_brock
  elseif(iabs(idmc) == 2) then
     if (abs(nloc) > 0) then
        call dmc_good_ps
     elseif(idiv_v <= 1) then
        call dmc_good
     else
        call dmc_good_inhom
     endif
  elseif(iabs(idmc) == 3) then
     if (abs(nloc) > 0) then
        stop 'pseudopotential version of antisym. pair-proj. not implemented'
     elseif(idiv_v <= 0) then
        call dmc_good_ap
     else
        stop 'inhom version of antisym. pair-proj. not implemented'
     endif
  else
     call die (lhere, 'iabs(idmc) must be 1 or 2 or 3.')
  endif

  end subroutine dmc_algorithm

! ==============================================================================
  subroutine dmc_brock
! ------------------------------------------------------------------------------
! Description   : Choice of dmc_brock routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
  if (l_mode_dmc_mov1) then
   call dmc_brock_mov1
  else
   call dmc_brock_movall
  endif

  end subroutine dmc_brock

! ==============================================================================
  subroutine dmc_good_ps
! ------------------------------------------------------------------------------
! Description   : Choice of dmc_good_ps routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
  if (l_mode_dmc_mov1) then
   call dmc_good_ps_mov1
  else
   call dmc_good_ps_movall
  endif

  end subroutine dmc_good_ps

! ==============================================================================
  subroutine dmc_good
! ------------------------------------------------------------------------------
! Description   : Choice of dmc_good routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
  if (l_mode_dmc_mov1) then
    call dmc_good_mov1
  else
    call dmc_good_movall
  endif

  end subroutine dmc_good

! ==============================================================================
  subroutine dmc_good_inhom
! ------------------------------------------------------------------------------
! Description   : Choice of dmc_good_inhom routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
  if (l_mode_dmc_mov1) then
   call dmc_good_inhom_mov1
  else
   call dmc_good_inhom_movall
  endif

  end subroutine dmc_good_inhom

! ==============================================================================
  subroutine splitj
! ------------------------------------------------------------------------------
! Description   : Choice of dmc splitj routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
    if (l_mode_dmc_mov1) then
     call splitj_mov1
    else
     call splitj_movall
    endif

  end subroutine splitj

! ==============================================================================
  subroutine acuest_dmc
! ------------------------------------------------------------------------------
! Description   : Choice of dmc acuest routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
# if defined (MPI)
    if (l_mode_dmc_mov1_mpi1) then
     call acuest_dmc_mov1_mpi1
    else
     call acuest_dmc_mov1_mpi2
    endif
# else
    if (l_mode_dmc_mov1) then
     call acuest_dmc_mov1
    else
     call acuest_dmc_movall
    endif
# endif

  end subroutine acuest_dmc

! ==============================================================================
  subroutine acues1_dmc
! ------------------------------------------------------------------------------
! Description   : Choice of dmc acues1 routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
# if defined (MPI)
    if (l_mode_dmc_mov1_mpi1) then
     call acues1_dmc_mov1_mpi1
    else
     call acues1_dmc_mov1_mpi2
    endif
# else
    if (l_mode_dmc_mov1) then
     call acues1_dmc_mov1
    else
     call acues1_dmc_movall
    endif
# endif

  end subroutine acues1_dmc

! ==============================================================================
  subroutine zeres0_dmc
! ------------------------------------------------------------------------------
! Description   : Choice of dmc zeres0 routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
# if defined (MPI)
    if (l_mode_dmc_mov1_mpi1) then
     call zeres0_dmc_mov1_mpi1
    else
     call zeres0_dmc_mov1_mpi2
    endif
# else
    if (l_mode_dmc_mov1) then
     call zeres0_dmc_mov1
    else
     call zeres0_dmc_movall
    endif
# endif

  end subroutine zeres0_dmc

! ==============================================================================
  subroutine zerest_dmc
! ------------------------------------------------------------------------------
! Description   : Choice of dmc zerest routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
# if defined (MPI)
    if (l_mode_dmc_mov1_mpi1) then
     call zerest_dmc_mov1_mpi1
    else
     call zerest_dmc_mov1_mpi2
    endif
# else
    if (l_mode_dmc_mov1) then
     call zerest_dmc_mov1
    else
     call zerest_dmc_movall
    endif
# endif

  end subroutine zerest_dmc

! ==============================================================================
  subroutine finwrt_dmc
! ------------------------------------------------------------------------------
! Description   : Choice of dmc finwrt routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! begin
# if defined (MPI)
    if (l_mode_dmc_mov1_mpi1) then
     call finwrt_dmc_mov1_mpi1
    else
     call finwrt_dmc_mov1_mpi2
    endif
# else
    if (l_mode_dmc_mov1) then
     call finwrt_dmc_mov1
    else
     call finwrt_dmc_movall
    endif
# endif

  end subroutine finwrt_dmc

! ==============================================================================
  subroutine send_jas(i)
! ------------------------------------------------------------------------------
! Description   : Choice of send_jas routine
!
! Created       : J. Toulouse, 09 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input
# if defined (MPI)
  character (len=max_string_len_rout), save :: lhere = 'send_jas'
# endif
  integer i

! local

! begin
# if defined (MPI)
    if (l_mode_dmc_mov1_mpi2) then
     call send_jas_mpi2(i)
    elseif (l_mode_dmc_mov1_mpi3) then
     call send_jas_mpi3(i)
    else
     call die (lhere,'unknown mode >'+trim(mode)+'<.')
    endif
# endif

  end subroutine send_jas

! ==============================================================================
  subroutine recv_jas(i)
! ------------------------------------------------------------------------------
! Description   : Choice of recv_jas routine
!
! Created       : J. Toulouse, 09 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input
# if defined (MPI)
  character (len=max_string_len_rout), save :: lhere = 'recv_jas'
# endif
  integer i

! local

! begin
# if defined (MPI)
    if (l_mode_dmc_mov1_mpi2) then
     call recv_jas_mpi2(i)
    elseif (l_mode_dmc_mov1_mpi3) then
     call recv_jas_mpi3(i)
    else
     call die (lhere,'unknown mode >'+trim(mode)+'<.')
    endif
# endif

  end subroutine recv_jas

end module dmc_mod
