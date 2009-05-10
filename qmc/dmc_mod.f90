module dmc_mod

  use all_tools_mod
  use control_mod
  use montecarlo_mod
  use average_mod
  use print_mod
  use restart_mod
  use allocations_mod

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
  include 'modules.h'
  implicit none
  include 'commons.h'

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
    call readforce
! parameters for secondary geometry wave function
    call wf_secondary
  else
    nwprod=1
    nwftype=1
    iwftype(1)=1
  endif

! allocations for DMC
  call common_allocations
  
  end subroutine dmc_init

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
  implicit none
  include 'commons.h'

! local
  integer i, istep, jj, k
  integer ngfmc
  character (len=27) fmt
  integer iwalk
  real(dp) passes, efin

! initial printing
  write(6,*)
  write(6,'(a)') '*********** START DMC CALCULATION  ***********'

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
  do i = 1, 2*nblkeq

     if ((i == nblkeq+1) .and. irstar /= 1) then
        call print_cpu_time_in_seconds ('End        of equilibration')
        call zerest_dmc
     endif

     do istep=1,nstep

        ipass=ipass+1
        call dmc_algorithm

!       averages at each step
        call acues1_dmc

!       walkers reconfiguration
        call splitj

      enddo ! istep

!   averages at each block
    call acuest_dmc

  enddo ! end of equilibration

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

!     compute averages
      call acues1_dmc
      call compute_averages_walk_step

!     walkers reconfiguration
      call splitj

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
              write(fmt,'(a1,i1,a21)')'(',ndim*nelec,'f14.8,i3,d12.4,f12.5)'
             elseif(ndim*nelec.lt.100) then
              write(fmt,'(a1,i2,a21)')'(',ndim*nelec,'f14.8,i3,d12.4,f12.5)'
             elseif(ndim*nelec.lt.1000) then
              write(fmt,'(a1,i3,a21)')'(',ndim*nelec,'f14.8,i3,d12.4,f12.5)'
             elseif(ndim*nelec.lt.10000) then
              write(fmt,'(a1,i4,a21)')'(',ndim*nelec,'f14.8,i3,d12.4,f12.5)'
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

         enddo ! istep

        call acuest_dmc
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
      write(6,'(''passes,wgcum='',9d12.4)') passes,(wgcum(i),i=1,nforce)

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

  end subroutine dmc

! ==============================================================================
  subroutine dmc_algorithm
! ------------------------------------------------------------------------------
! Description   : Choice of DMC algoritm
!
! Created       : J. Toulouse, 27 Mar 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

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
  else
     call die (lhere, 'iabs(idmc) must be 1 or 2.')
  endif

  end subroutine dmc_algorithm

! ==============================================================================
  subroutine dmc_brock
! ------------------------------------------------------------------------------
! Description   : Choice of dmc_brock routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
