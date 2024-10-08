module vmc_mod

  use all_tools_mod
  use allocations_mod

! Declaration of global variables and default values
  logical :: l_walkers_available = .false.

  contains

! ==============================================================================
  subroutine vmc_init
! ------------------------------------------------------------------------------
! Description   : initialization for vmc run
!
! Created       : J. Toulouse, 10 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
!  use psi_type_mod !TA
  use psi_utils_mod, only: set_oldvar_pointers, set_newvar_pointers
  use psi_dummy_mod
  implicit none

  call object_provide ('nelec')
  call object_provide ('nforce')

! local
  call my_second(0,'begin ')

! correlated sampling
  if (nforce > 1) then
!   force parameters
    if (.not. use_parser) then
      call readforce
    endif
!   parameters for secondary geometry wave function
    call wf_secondary
  else
    nwftype=1
    iwftype(1)=1
  endif

! Initialize starting MC configuration if warkers not already available from a previous run
  if (.not. l_walkers_available) then
   call mc_configs_read
   l_walkers_available = .true.
  endif

! If we are moving one electron at a time, then we need to initialize
! xnew, since only the first electron gets initialized in metrop
  call alloc ('xnew', xnew, 3, nelec)
  if (irstar /= 1) then
     xnew (1:ndim,1:nelec) = xold (1:ndim,1:nelec)
  endif

! default: no gradient nor Hessian
  igradhess=0

! various allocations for VMC
  call alloc ('vold', vold, 3, nelec)
  call alloc ('vnew', vnew, 3, nelec)
  call alloc ('psi2o', psi2o, nwf)
  call alloc ('psi2n', psi2n, nwf)
  call alloc ('eold', eold, nwf)
  call alloc ('enew', enew, nwf)
  call alloc ('rmino', rmino, nelec)
  call alloc ('rminn', rminn, nelec)
  call alloc ('rvmino', rvmino, 3, nelec)
  call alloc ('rvminn', rvminn, 3, nelec)
  call alloc ('rminon', rminon, nelec)
  call alloc ('rminno', rminno, nelec)
  call alloc ('rvminon', rvminon, 3, nelec)
  call alloc ('rvminno', rvminno, 3, nelec)
  call alloc ('delttn', delttn, nelec)
  call alloc ('nearesto', nearesto, nelec)
  call alloc ('nearestn', nearestn, nelec)
  call alloc ('pot_ee_old', pot_ee_old, nelec)
  call alloc ('pot_ee_new', pot_ee_new, nelec)
  call alloc ('pot_ee', pot_ee, nelec)
  call common_allocations

  call set_newvar_pointers(1,psi_dummy2)
  call set_oldvar_pointers(psi_dummy)

! set current walker and current walker weight to 1 for VMC
  current_walker = 1
  call object_modified ('current_walker')
! Warning: tmp moved to vmc/vmc.f
! current_walker_weight = 1
  call object_modified ('current_walker_weight')

  end subroutine vmc_init

! ==============================================================================
  subroutine vmc_release
! ------------------------------------------------------------------------------
! Description   : release arrays allocated for VMC
!
! Created       : J. Toulouse, 16 Oct 2009
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  call release ('xnew', xnew)
  call release ('vold', vold)
  call release ('vnew', vnew)
  call release ('psi2o', psi2o)
  call release ('psi2n', psi2n)
  call release ('eold', eold)
  call release ('enew', enew)
  call release ('rmino', rmino)
  call release ('rminn', rminn)
  call release ('rvmino', rvmino)
  call release ('rvminn', rvminn)
  call release ('rminon', rminon)
  call release ('rminno', rminno)
  call release ('rvminon', rvminon)
  call release ('rvminno', rvminno)
  call release ('delttn', delttn)
  call release ('nearesto', nearesto)
  call release ('nearestn', nearestn)
  call release ('pot_ee_old', pot_ee_old)
  call release ('pot_ee_new', pot_ee_new)
  call release ('pot_ee', pot_ee)

  end subroutine vmc_release

! ==============================================================================
  subroutine vmc_run
! ------------------------------------------------------------------------------
! Description   : launch a vmc run
!
! Created       : J. Toulouse, 12 Feb 2007
! ------------------------------------------------------------------------------
  implicit none

! local

  call vmc_init
  call vmc

  end subroutine vmc_run

! ==============================================================================
  subroutine metrop(l)
! ------------------------------------------------------------------------------
! Description   : Choice of metrop routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input
  integer l

! local

! begin
  if (l_mode_vmc_mov1) then
   call metrop_mov1(l)
  else
   call metrop_movall(l)
  endif

  end subroutine metrop

! ==============================================================================
  subroutine metrop_polar(l)
! ------------------------------------------------------------------------------
! Description   : Choice of metrop_polar routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input
  integer l

! local

! begin
  if (l_mode_vmc_mov1) then
   call metrop_polar_mov1(l)
  else
   call metrop_polar_movall(l)
  endif

  end subroutine metrop_polar

! ==============================================================================
  subroutine metrop_slat(l)
! ------------------------------------------------------------------------------
! Description   : Choice of metrop_slat routine
!
! Created       : J. Toulouse, 08 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input
  integer l

! local

! begin
  if (l_mode_vmc_mov1) then
   call metrop_slat_mov1(l)
  else
   call metrop_slat_movall(l)
  endif

  end subroutine metrop_slat

end module vmc_mod
