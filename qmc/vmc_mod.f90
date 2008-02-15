module vmc_mod

  use all_tools_mod

! Declaration of global variables and default values

  contains

! ==============================================================================
  subroutine vmc_init
! ------------------------------------------------------------------------------
! Description   : initialization for vmc run
!
! Created       : J. Toulouse, 10 Jan 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character (len=max_string_len_rout), save :: lhere = 'vmc_init'

  call my_second(0,'begin ')

! EFP
  if (nefp >0 ) call readbas ! to remove

! correlated sampling
  if (nforce > 1) then
!   force parameters
    call readforce
!   parameters for secondary geometry wave function
    call wf_secondary
  else
    nwftype=1
    iwftype(1)=1
  endif

! Initialize starting MC configuration
  call mc_configs_read

! If we are moving one electron at a time, then we need to initialize
! xnew, since only the first electron gets initialized in metrop
  if (irstar /= 1) then
     xnew (1:ndim,1:nelec) = xold (1:ndim,1:nelec)
  endif

! default: no gradient nor Hessian
  igradhess=0

  end subroutine vmc_init

! ==============================================================================
  subroutine vmc_run
! ------------------------------------------------------------------------------
! Description   : launch a vmc run
!
! Created       : J. Toulouse, 12 Feb 2007
! ------------------------------------------------------------------------------
  implicit none

! local
  character (len=max_string_len_rout), save :: lhere = 'vmc_run'

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
  implicit none
  include 'commons.h'

! input
  integer l

! local
  character (len=max_string_len_rout), save :: lhere = 'metrop'

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
  implicit none
  include 'commons.h'

! input
  integer l

! local
  character (len=max_string_len_rout), save :: lhere = 'metrop_polar'

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
  implicit none
  include 'commons.h'

! input
  integer l

! local
  character (len=max_string_len_rout), save :: lhere = 'metrop_slat'

! begin
  if (l_mode_vmc_mov1) then
   call metrop_slat_mov1(l)
  else
   call metrop_slat_movall(l)
  endif

  end subroutine metrop_slat

end module vmc_mod
