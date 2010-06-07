module opt_ovlp_fn_mod

  use all_tools_mod
  use jastrow_mod
  use csfs_mod
!  use opt_common_mod
!  use deriv_mod

! Declaration of global variables and default values
  logical                         :: l_weighted_overlap = .true.
  real(dp), allocatable           :: csf_over_psit_j (:)
  real(dp), allocatable           :: csf_over_psit_j_av (:)
  real(dp), allocatable           :: delta_ovlp_fn (:)

  contains

!===========================================================================
  subroutine opt_ovlp_fn_menu
!---------------------------------------------------------------------------
! Description : menu for overlap with FN wavefn. optimization
!
! Created     : Cyrus Umrigar and Frank Petruzielo 7 Jun 2010
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'opt_ovlp_fn_menu'

! begin

! loop over menu lines
  do
  call get_next_word(word)

  select case (trim(word))
  case ('help')
   write(6,'(a)') ' HELP for overlap fixed-node optimization menu'
   write(6,'(a)') '    weighted_overlap = [logical] : have zero variance estimator by calculating the difference between the FN and the trial wavefn. (default=true)'
   write(6,'(a)') ' end'

  case ('weighted_overlap')
   call get_next_value(l_weighted_overlap)

  case default
   call die(lhere, 'unknown word >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  end subroutine opt_ovlp_fn_menu

! ==============================================================================
  subroutine csf_over_psit_j_bld
! ------------------------------------------------------------------------------
! Description   : csf divided by Psi_T times Jastrow^2
!
! Created       : Cyrus Umrigar and Frank Petruzielo, 7 Jun 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer csf_i

! header
  if(header_exe) then

   call object_create('csf_over_psit_j')
   call object_average_define ('csf_over_psit_j', 'csf_over_psit_j_av')

   call object_needed('ncsf')
   call object_needed('csf_over_psid')
   call object_needed('psi_jas')

   return

  endif

! begin

! allocations
  call object_alloc('csf_over_psit_j', csf_over_psit_j, ncsf)
  call object_alloc('csf_over_psit_j_av', csf_over_psit_j_av, ncsf)

  do csf_i = 1, ncsf
    csf_over_psit_j(csf_i) = csf_over_psid (csf_i)/psi_jas**2
  enddo ! csf_i

  end subroutine csf_over_psit_j_bld

! ==============================================================================
  subroutine delta_ovlp_fn_bld
! ------------------------------------------------------------------------------
! Description   : variation of parameters for overlap FN method
!
! Created       : Cyrus Umrigar and Frank Petruzielo, 7 Jun 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_bld'
  integer iparm


! header
  if(header_exe) then

   call object_create('delta_ovlp_fn')

   call object_needed('csf_over_psit_j_av')
   call object_needed ('param_nb')

   return

  endif

! begin

! allocation
  call object_alloc('delta_ovlp_fn', delta_ovlp_fn, param_nb)

  call object_provide ('csf_coef')

  do iparm = 1, param_nb
    delta_ovlp_fn(iparm) = csf_over_psit_j_av(iparm) - csf_coef(iparm, 1)
    write(6,*) "temp: :",csf_coef(iparm, 1), delta_ovlp_fn(iparm)
  enddo

  end subroutine delta_ovlp_fn_bld

end module opt_ovlp_fn_mod
