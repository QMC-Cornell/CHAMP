module opt_ptb_mod

  use all_tools_mod
  use montecarlo_mod
  use deriv_mod
  use opt_lin_mod

! Declaration of global variables and default values
  real(dp), allocatable          :: delta_ptb (:)
  real(dp), allocatable          :: e_ptb (:)
  real(dp), allocatable          :: delta_e_ptb (:)

  contains

!===========================================================================
  subroutine opt_ptb_menu
!---------------------------------------------------------------------------
! Description : menu for perturbative optimization method
!
! Created     : J. Toulouse, 24 Apr 2006
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'opt_ptb_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)')  ' HELP for perturbative menu:'
   write(6,'(a)')  ' perturbative'
   write(6,'(a)')  '   use_orbital_eigenvalues = [logical] : use orbital eigenvalues for energy denominators? (default=false)'
   write(6,'(a)')  '   diagonal_overlap = [logical] : approximate overlap matrix of derivatives by its diagonal? (default=false)'
   write(6,'(a)')  '   lambda = [real] : scaling factor for correction of MCSCF orbital eigenvalues (default=0.3)'
!   write(6,'(a)')  ':   delta_e_ptb 0.1 3. 5. end : read energy denominators'
   write(6,'(a)')  ': end'

  case ('use_orbital_eigenvalues')
   call get_next_value (l_opt_orb_eig)

  case ('diagonal_overlap')
   call get_next_value (l_diagonal_overlap)

  case ('lambda')
   call get_next_value (lambda)

!  case ('delta_e_ptb')
!   call delta_e_ptb_rd

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines

  end subroutine opt_ptb_menu

! ==============================================================================
  subroutine e_ptb_bld
! ------------------------------------------------------------------------------
! Description   : energies for perturbative method
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer i, j

! header
  if (header_exe) then

   call object_create ('e_ptb')

   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('dpsi_av')
   call object_needed ('dpsi_eloc_av')
   call object_needed ('dpsi_deloc_covar')
   call object_needed ('dpsi_sq_eloc_av')
   call object_needed ('dpsi_dpsi_covar')
   call object_needed ('eloc_av')

   return

  endif

! begin
  call object_alloc ('e_ptb', e_ptb, param_nb)

  do i = 1, param_nb
    e_ptb (i) = (dpsi_deloc_covar (i, i) + dpsi_sq_eloc_av (i)                  &
                       - dpsi_av (i) * dpsi_eloc_av (i) - dpsi_av (i) * dpsi_eloc_av (i) &
                       + dpsi_av (i) * dpsi_av (i) * eloc_av ) / dpsi_dpsi_covar (i,i)
  enddo

  end subroutine e_ptb_bld

!===========================================================================
  subroutine delta_e_ptb_rd
!---------------------------------------------------------------------------
! Description : read delta_e_ptb
!
! Created     : J. Toulouse, 24 May 2006
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_e_ptb_rd'
  integer delta_e_ptb_nb

! begin
  call get_next_value_list ('delta_e_ptb', delta_e_ptb, delta_e_ptb_nb)

  call object_provide ('param_nb')
  if (delta_e_ptb_nb /= param_nb) then
   call die (lhere, 'delta_e_ptb_nb='+delta_e_ptb_nb+' /= param_nb='+param_nb)
  endif

  call object_modified ('delta_e_ptb')

  call object_write (lhere, 'delta_e_ptb')

  end subroutine delta_e_ptb_rd

! ==============================================================================
  subroutine delta_e_ptb_bld
! ------------------------------------------------------------------------------
! Description   : energy denominators for perturbative method
! Description   : calculated only once
!
! Created       : J. Toulouse, 04 Feb 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_e_ptb_bld'

! header
  if (header_exe) then

   call object_create ('delta_e_ptb')

   call object_needed ('param_nb')

   return

  endif

! begin

! allocations
  call object_alloc ('delta_e_ptb', delta_e_ptb, param_nb)

  if (l_opt_orb_eig) then
    if (param_nb /= param_orb_nb) then
     call die (here, 'param_nb='+param_nb+' /=  param_orb_nb='+param_orb_nb)
    endif
    call object_provide ('delta_eps')
    delta_e_ptb (:) = delta_eps (:)
  else
   call object_provide ('e_ptb')
   call object_provide ('eloc_av')
   delta_e_ptb (:) = e_ptb (:) - eloc_av
  endif

  call object_write (lhere, 'delta_e_ptb')

 end subroutine delta_e_ptb_bld

! ==============================================================================
  subroutine delta_ptb_bld
! ------------------------------------------------------------------------------
! Description   : find variation of parameters for perturbative method
!
! Created       : J. Toulouse, 04 Feb 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_ptb_bld'
  integer parm_i, parm_j

! header
  if (header_exe) then

   call object_create ('delta_ptb')

   call object_needed ('param_nb')
   call object_needed ('gradient')
   call object_needed ('delta_e_ptb')
   call object_needed ('diag_stab')

   return

  endif

! begin
! allocations
  call object_alloc ('delta_ptb', delta_ptb, param_nb)

  delta_ptb (:) = 0.d0

! with only diagonal of overlap
  if (l_diagonal_overlap) then

   if (param_nb /= param_orb_nb) then
     call die (here, 'param_nb='+param_nb+' /=  param_orb_nb='+param_orb_nb)
   endif
   call object_provide_in_node (lhere, 'dpsi_sq_covar')
   do parm_i = 1, param_nb
     delta_ptb (parm_i) = - (1.d0/(dpsi_sq_covar (parm_i) * (delta_e_ptb (parm_i) + diag_stab))) * gradient (parm_i)/2.d0
   enddo

! with full overlap
  else

   call object_provide_in_node (lhere, 'dpsi_dpsi_covar_inv')
   do parm_i = 1, param_nb
    do parm_j = 1, param_nb
     delta_ptb (parm_i) = delta_ptb (parm_i) - (dpsi_dpsi_covar_inv (parm_i, parm_j)/(delta_e_ptb (parm_i) + diag_stab)) * gradient (parm_j)/2.d0
    enddo
   enddo

  endif

  end subroutine delta_ptb_bld

end module opt_ptb_mod
