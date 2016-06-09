 module forces_mod

  use all_tools_mod
  use grid_mod
  use electrons_mod
  use nuclei_mod
  use psi_mod
  use montecarlo_mod
  use eloc_mod
  use opt_lin_mod
  use forces_pulay_mod

! Declaration of global variables and default values
# if defined (PATHSCALE)
   character(len=max_string_len) :: forces_list (max_string_array_len) ! for pathscale compiler
# else
   character(len=max_string_len), allocatable  :: forces_list (:)
# endif
  logical                                     :: l_forces_bare = .false.
  logical                                     :: l_forces_zv = .false.
  logical                                     :: l_forces_zv_deriv = .false.
  logical                                     :: l_forces_zv_linear = .false.
  logical                                     :: l_forces_zv_deriv_linear = .false.
  logical                                     :: l_forces_zvzb = .false.
  logical                                     :: l_forces_pulay = .false.
  logical                                     :: l_forces_zv_pulay = .false.
  logical                                     :: l_forces_zv_deriv_pulay = .true.
!  logical                                     :: l_forces_zv_deriv_pulay = .false.
  real(dp), allocatable                       :: forces_nn (:)
  real(dp), allocatable                       :: forces_bare (:)
  real(dp), allocatable                       :: forces_bare_av (:)
  real(dp), allocatable                       :: forces_bare_av_err (:)
  real(dp), allocatable                       :: forces_zv (:)
  real(dp), allocatable                       :: forces_zv_av (:)
  real(dp), allocatable                       :: forces_zv_av_var (:)
  real(dp), allocatable                       :: forces_zv_av_err (:)
  real(dp), allocatable                       :: forces_zv_sq (:)
  real(dp), allocatable                       :: forces_zv_sq_av (:)
  real(dp), allocatable                       :: forces_zv_var (:)
  real(dp), allocatable                       :: forces_q (:)
  real(dp), allocatable                       :: forces_q_av (:)
  real(dp), allocatable                       :: forces_q_av_var (:)
  real(dp), allocatable                       :: forces_q_av_err (:)
  real(dp), allocatable                       :: forces_q_eloc (:)
  real(dp), allocatable                       :: forces_q_eloc_av (:)
  real(dp), allocatable                       :: forces_q_eloc_av_var (:)
  real(dp), allocatable                       :: forces_q_eloc_av_err (:)
  real(dp), allocatable                       :: forces_zvzb (:)
  real(dp), allocatable                       :: forces_zvzb_av (:)
  real(dp), allocatable                       :: forces_zvzb_av_var (:)
  real(dp), allocatable                       :: forces_zvzb_av_err (:)
  real(dp), allocatable                       :: forces_zv_av_eloc_av_covar (:)
  real(dp), allocatable                       :: forces_zv_av_q_eloc_av_covar (:)
  real(dp), allocatable                       :: forces_zv_av_q_av_covar (:)
  real(dp), allocatable                       :: forces_q_eloc_av_eloc_av_covar (:)
  real(dp), allocatable                       :: forces_q_eloc_av_q_av_covar (:)
  real(dp), allocatable                       :: forces_q_av_eloc_av_covar (:)
  real(dp), allocatable                       :: forces_zv_linear_av (:)
  real(dp), allocatable                       :: forces_zv_linear_av_var (:)
  real(dp), allocatable                       :: forces_zv_linear_av_err (:)
  real(dp), allocatable                       :: forces_zv_linear_var (:)
  real(dp), allocatable                       :: forces_zv_linear_coef (:,:)
  real(dp), allocatable                       :: forces_zv_deriv_linear_av (:)
  real(dp), allocatable                       :: forces_zv_deriv_linear_av_var (:)
  real(dp), allocatable                       :: forces_zv_deriv_linear_av_err (:)
  real(dp), allocatable                       :: forces_zv_deriv_linear_var (:)
  real(dp), allocatable                       :: forces_zv_deriv_linear_coef (:,:)
  real(dp), allocatable                       :: forces_zv_deriv (:)
  real(dp), allocatable                       :: forces_zv_deriv_av (:)
  real(dp), allocatable                       :: forces_zv_deriv_var (:)
  real(dp), allocatable                       :: forces_zv_deriv_av_var (:)
  real(dp), allocatable                       :: forces_zv_deriv_av_err (:)
  real(dp), allocatable                       :: forces_zv_deriv_sq (:)
  real(dp), allocatable                       :: forces_zv_deriv_sq_av (:)
  real(dp), allocatable                       :: forces_zv_deloc (:,:)
  real(dp), allocatable                       :: forces_zv_deloc_av (:,:)
  real(dp), allocatable                       :: forces_zv_deloc_covar (:,:)
  real(dp), allocatable                       :: forces_zv_av_deloc_av_covar (:,:)
  real(dp), allocatable                       :: forces_zv_deriv_deloc (:,:)
  real(dp), allocatable                       :: forces_zv_deriv_deloc_av (:,:)
  real(dp), allocatable                       :: forces_zv_deriv_deloc_covar (:,:)
  real(dp), allocatable                       :: forces_zv_deriv_av_deloc_av_covar (:,:)
  real(dp), allocatable                       :: forces_zv_pulay_av (:)
  real(dp), allocatable                       :: forces_zv_pulay_av_var (:)
  real(dp), allocatable                       :: forces_zv_pulay_av_err (:)
  real(dp), allocatable                       :: forces_zv_av_forces_pulay_av_covar (:)
  real(dp), allocatable                       :: forces_zv_deriv_pulay_av (:)
  real(dp), allocatable                       :: forces_zv_deriv_pulay_av_var (:)
  real(dp), allocatable                       :: forces_zv_deriv_pulay_av_err (:)
  real(dp), allocatable                       :: forces_zv_deriv_av_forces_pulay_av_covar (:)

  logical                                     :: l_eloc_av_fixed = .false.
  logical                                     :: l_forces_q_av_fixed = .false.
  real(dp)                                    :: eloc_av_fixed = 0.d0
  real(dp), allocatable                       :: forces_q_av_fixed (:)

  contains

!===========================================================================
  subroutine forces_menu
!---------------------------------------------------------------------------
! Description : menu for calculation of forces
!
! Created     : J. Toulouse, 26 Jul 2007
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'forces_menu'
  character(len=max_string_len) estimator

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of forces menu ---------------------------------------------------------------------------------'

! not for a pseudopotential
  if (nloc > 0) then
    call die (here, 'calculation of forces not implemented for pseudopotentials!')
  endif

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))
   case ('help')
    write(6,'(a)') 'HELP for menu forces'
    write(6,'(a)') 'forces'
    write(6,'(a)') ' components 1x 1y 1z 2x 2y 2z end'
    write(6,'(a)') ' estimator = [string] estimator to use:'
    write(6,'(a)') '           = bare : bare Hellmann-Feynamn estimator'
    write(6,'(a)') '           = zv : simplest zero-variance estimator'
    write(6,'(a)') '           = zv_deriv : zero-variance estimator using wave function derivatives wrt nuclear coordinates'
    write(6,'(a)') '           = zv_linear : zero-variance estimator using wave function derivatives wrt parameters'
    write(6,'(a)') '           = zv_deriv_linear : zero-variance estimator using wave function derivatives wrt nuclear coordinates and parameters'
    write(6,'(a)') '           = zvzb : simplest zero-variance zero-bias estimator'
    write(6,'(a)') '           = zv_pulay : simplest zero-variance + Pulay estimator'
    write(6,'(a)') '           = zv_deriv_pulay : zero-variance using wave function derivatives wrt nuclear coordinates + Pulay estimator (default)'
    write(6,'(a)') '           = pulay : Pulay contribution to the force'
    write(6,'(a)') ' eloc_av_fixed  = [real] fixed value of average of local energy to use in ZB term'
    write(6,'(a)') ' forces_q_av_fixed  list of reals end: fixed value of average of Q to use in ZB term'
    write(6,'(a)') 'end'

   case ('components')
    call forces_list_rd

   case ('estimator')
    call get_next_value (estimator)
    select case (estimator)
     case ('bare')
       l_forces_bare = .true.
     case ('zv')
       l_forces_zv = .true.
     case ('zv_deriv')
       l_forces_zv_deriv = .true.
     case ('zv_linear')
       l_forces_zv_linear = .true.
       l_forces_zv = .true.
     case ('zv_deriv_linear')
       l_forces_zv_deriv_linear = .true.
       l_forces_zv_deriv = .true.
     case ('zvzb')
       l_forces_zvzb = .true.
       l_forces_zv   = .true.
     case ('zv_pulay')
       l_forces_zv = .true.
       l_forces_zv_pulay = .true.
       l_forces_pulay = .true.
     case ('zv_deriv_pulay')
       l_forces_zv_deriv = .true.
       l_forces_zv_deriv_pulay = .true.
       l_forces_pulay = .true.
     case ('pulay')
       l_forces_pulay = .true.
     case default; call die (lhere, 'unknown keyword >'+trim(estimator)+'<.')
    end select

   case ('pulay')
    call get_next_value (l_forces_pulay)

   case ('eloc_av_fixed')
    call get_next_value (eloc_av_fixed)
    call object_modified ('eloc_av_fixed')
    l_eloc_av_fixed = .true.

   case ('forces_q_av_fixed')
    call get_next_value_list_object ('forces_q_av_fixed', forces_q_av_fixed, forces_nb)
    l_forces_q_av_fixed = .true.

   case ('end')
    exit

   case default
    call die (lhere, 'unknown keyword >'+trim(word)+'<.')

  end select

  enddo ! end loop over menu lines

! dependency between estimators
  if (l_forces_zv_linear) then
   l_forces_zv = .true.
  endif
  if (l_forces_zv_deriv_linear) then
   l_forces_zv_deriv = .true.
  endif
  if (l_forces_zvzb) then
   l_forces_zv = .true.
  endif
  if (l_forces_zv_pulay) then
   l_forces_zv = .true.
   l_forces_pulay = .true.
  endif
  if (l_forces_zv_deriv_pulay) then
   l_forces_zv_deriv = .true.
   l_forces_pulay = .true.
  endif

  if (l_forces_bare) then
   write(6,'(a)') ' Forces will be calculated with bare (Hellmann-Feynman) estimator.'
  endif
  if (l_forces_zv) then
   write(6,'(a)') ' Forces will be calculated with simplest zero-variance estimator.'
  endif
  if (l_forces_zv_deriv) then
   write(6,'(a)') ' Forces will be calculated with zero-variance estimator using wave function derivatives wrt nuclear coordinates.'
  endif
  if (l_forces_zv_linear) then
   write(6,'(a)') ' Forces will be calculated with zero-variance estimator using wave function derivatives wrt parameters.'
  endif
  if (l_forces_zv_deriv_linear) then
   write(6,'(a)') ' Forces will be calculated with zero-variance estimator using wave function derivatives wrt nuclear coordinates and parameters.'
  endif
  if (l_forces_zvzb) then
   write(6,'(a)') ' Forces will be calculated with simplest zero-variance zero-bias estimator.'
  endif
  if (l_forces_pulay) then
   write(6,'(a)') ' Pulay contribution to forces will be calculated.'
  endif
  if (l_forces_zv_pulay) then
   write(6,'(a)') ' Forces will be calculated with zero-variance + Pulay estimator.'
  endif
  if (l_forces_zv_deriv_pulay) then
   write(6,'(a)') ' Forces will be calculated with zero-variance using wave function derivatives wrt nuclear coordinates + Pulay estimator.'
  endif

! request averages and statistical errors
  if (l_forces_bare) then
   call object_average_request ('forces_bare_av')
   call object_error_request ('forces_bare_av_err')
  endif

  if (l_forces_zv) then
   call object_average_request ('forces_zv_av')
   call object_error_request ('forces_zv_av_err')
   call object_average_request ('forces_zv_sq_av')
  endif

  if (l_forces_zv_deriv) then
   call object_average_request ('forces_zv_deriv_av')
   call object_error_request ('forces_zv_deriv_av_err')
   call object_average_request ('forces_zv_deriv_sq_av')
  endif

  if (l_forces_zv_linear) then
   call object_average_request ('deloc_av')
   call object_average_request ('deloc_deloc_av')
   call object_average_request ('forces_zv_deloc_av')
   call object_covariance_request ('deloc_av_deloc_av_covar')
   call object_covariance_request ('forces_zv_av_deloc_av_covar')
   call object_error_request ('forces_zv_linear_av_err')
  endif

  if (l_forces_zv_deriv_linear) then
   call object_average_request ('deloc_av')
   call object_average_request ('deloc_deloc_av')
   call object_average_request ('forces_zv_deriv_deloc_av')
   call object_covariance_request ('deloc_av_deloc_av_covar')
   call object_covariance_request ('forces_zv_deriv_av_deloc_av_covar')
   call object_error_request ('forces_zv_deriv_linear_av_err')
  endif

  if (l_forces_zvzb) then
    call object_average_request ('forces_q_av')
    call object_average_request ('forces_q_eloc_av')
    call object_variance_request ('forces_q_eloc_av_var')
    call object_variance_request ('forces_q_av_var')
    call object_error_request ('forces_zvzb_av_err')
    call object_covariance_request ('forces_zv_av_eloc_av_covar')
    call object_covariance_request ('forces_zv_av_q_eloc_av_covar')
    call object_covariance_request ('forces_zv_av_q_av_covar')
    call object_covariance_request ('forces_q_eloc_av_eloc_av_covar')
    call object_covariance_request ('forces_q_eloc_av_q_av_covar')
    call object_covariance_request ('forces_q_av_eloc_av_covar')
  endif

  if (l_forces_pulay) then
   call object_average_request ('dpsi_rn_av')
   call object_average_request ('dpsi_rn_eloc_av')
   call object_error_request ('forces_pulay_av_err')
   call object_variance_request ('dpsi_rn_av_var')
   call object_variance_request ('dpsi_rn_eloc_av_var')
   call object_covariance_request ('dpsi_rn_eloc_av_dpsi_rn_av_covar')
   call object_covariance_request ('dpsi_rn_eloc_av_eloc_av_covar')
   call object_covariance_request ('dpsi_rn_av_eloc_av_covar')
  endif

  if (l_forces_zv_pulay) then
   call object_error_request ('forces_zv_pulay_av_err')
   call object_covariance_request ('forces_zv_av_forces_pulay_av_covar')
  endif

  if (l_forces_zv_deriv_pulay) then
   call object_error_request ('forces_zv_deriv_pulay_av_err')
   call object_covariance_request ('forces_zv_deriv_av_forces_pulay_av_covar')
  endif

  call routine_write_final_request ('forces_wrt')

  write(6,'(a)') 'End of forces menu ---------------------------------------------------------------------------------------'

  end subroutine forces_menu

! ==============================================================================
  subroutine forces_list_rd
! ------------------------------------------------------------------------------
! Description   : read values for list of force components to calculate
!
! Created       : J. Toulouse, 28 Mar 2009
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'forces_list_rd'
  integer force_i

! begin
  call get_next_value_list_object ('forces_list', forces_list, forces_nb)
  call object_modified ('forces_nb')

! determine forces' centers and directions
  call object_alloc ('forces_cent', forces_cent, forces_nb)
  call object_alloc ('forces_direct', forces_direct, forces_nb)
  do force_i = 1, forces_nb
      forces_cent (force_i) = string_to_integer (forces_list (force_i)(1:1))
      select case (forces_list (force_i)(2:2))
      case ('x')
       forces_direct (force_i) = 1
      case ('y')
       forces_direct (force_i) = 2
      case ('z')
       forces_direct (force_i) = 3
      case default
       call die (lhere, 'Second character of the forces_list (force_i)='+forces_list (force_i)(2:2)+' must be x, y or z.')
      end select
  enddo
  call object_modified ('forces_cent')
  call object_modified ('forces_direct')

  end subroutine forces_list_rd

! ==============================================================================
  subroutine forces_list_bld
! ------------------------------------------------------------------------------
! Description   : default values for list of force components to calculate
!
! Created       : J. Toulouse, 28 Mar 2009
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, cent_i

! begin

! header
  if (header_exe) then

   call object_create ('forces_nb')
   call object_create ('forces_list')
   call object_create ('forces_cent')
   call object_create ('forces_direct')

   call object_needed ('ncent')

   return

  endif

! allocations
  forces_nb = 3 * ncent
  call object_alloc ('forces_list', forces_list, forces_nb)
  call object_alloc ('forces_cent', forces_cent, forces_nb)
  call object_alloc ('forces_direct', forces_direct, forces_nb)

  force_i = 0
  do cent_i = 1, ncent

!  x component
   force_i = force_i + 1
   forces_list (force_i) = string(cent_i) + 'x'
   forces_cent (force_i) = cent_i
   forces_direct (force_i) = 1

!  y component
   force_i = force_i + 1
   forces_list (force_i) = string(cent_i) + 'y'
   forces_cent (force_i) = cent_i
   forces_direct (force_i) = 2

!  z component
   force_i = force_i + 1
   forces_list (force_i) = string(cent_i) + 'z'
   forces_cent (force_i) = cent_i
   forces_direct (force_i) = 3

  enddo ! cent_i

  call require (here, 'force_i = forces_nb', force_i == forces_nb)

  end subroutine forces_list_bld

! ==============================================================================
  subroutine forces_nn_bld
! ------------------------------------------------------------------------------
! Description   : nuclei-nuclei contribution to forces
!
! Created       : J. Toulouse, 26 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, dim_i, cent_i, cent_j

! begin

! header
  if (header_exe) then

   call object_create ('forces_nn')

   call object_needed ('forces_nb')
   call object_needed ('forces_cent')
   call object_needed ('forces_direct')
   call object_needed ('ncent')
   call object_needed ('cent')
   call object_needed ('dist_nn')
   call object_needed ('iwctype')
   call object_needed ('znuc')

   return

  endif

! allocations
  call object_alloc ('forces_nn', forces_nn, forces_nb)

  do force_i = 1, forces_nb

   cent_i = forces_cent (force_i)
   dim_i = forces_direct (force_i)

   forces_nn (force_i) = 0.d0

!  loop over other nuclei
   do cent_j = 1, ncent
    if (cent_j == cent_i) cycle
    forces_nn (force_i) = forces_nn (force_i) + znuc(iwctype(cent_i)) * znuc(iwctype(cent_j)) * (cent (dim_i, cent_i) - cent (dim_i, cent_j))/ dist_nn (cent_i, cent_j)**3
   enddo ! cent_j

  enddo ! force_i

  end subroutine forces_nn_bld

! ==============================================================================
  subroutine forces_bare_bld
! ------------------------------------------------------------------------------
! Description   : Bare Hellmann-Feynman forces
!
! Created       : J. Toulouse, 26 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, elec_i, cent_i, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('forces_bare')
   call object_average_define ('forces_bare', 'forces_bare_av')
   call object_error_define ('forces_bare_av', 'forces_bare_av_err')

   call object_needed ('forces_nb')
   call object_needed ('forces_nn')
   call object_needed ('znuc')
   call object_needed ('iwctype')
   call object_needed ('vec_en_xyz_wlk')
   call object_needed ('dist_en_wlk')

   return

  endif

! allocations
  call object_alloc ('forces_bare', forces_bare, forces_nb)
  call object_alloc ('forces_bare_av', forces_bare_av, forces_nb)
  call object_alloc ('forces_bare_av_err', forces_bare_av_err, forces_nb)


  do force_i = 1, forces_nb
   cent_i = forces_cent (force_i)
   dim_i = forces_direct (force_i)

   forces_bare (force_i) = forces_nn (force_i)

!  loop over electrons
   do elec_i = 1, nelec
    forces_bare (force_i) = forces_bare (force_i) + znuc(iwctype(cent_i)) * vec_en_xyz_wlk (dim_i, elec_i, cent_i, 1) / dist_en_wlk (elec_i, cent_i, 1)**3
   enddo ! elec_i

  enddo ! force_i

  end subroutine forces_bare_bld

! ==============================================================================
  subroutine forces_zv_bld
! ------------------------------------------------------------------------------
! Description   : Simple ZV estimator for forces
! Description   : Assaraf and Caffarel, J. Chem. Phys. 113, 4028 (2000)
!
! Created       : J. Toulouse, 26 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, elec_i, cent_i, dim_i, dim_k
  real(dp) dotproduct

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv')
   call object_average_define ('forces_zv', 'forces_zv_av')
   call object_variance_define ('forces_zv_av', 'forces_zv_av_var')
   call object_error_define ('forces_zv_av', 'forces_zv_av_err')
   call object_covariance_define ('forces_zv_av','forces_pulay_av','forces_zv_av_forces_pulay_av_covar')

   call object_needed ('forces_nb')
   call object_needed ('ndim')
   call object_needed ('forces_nn')
   call object_needed ('znuc')
   call object_needed ('iwctype')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('dist_en_wlk')
   call object_needed ('vec_en_xyz_wlk')

   return

  endif

! allocations
  call object_alloc ('forces_zv', forces_zv, forces_nb)
  call object_alloc ('forces_zv_av', forces_zv_av, forces_nb)
  call object_alloc ('forces_zv_av_var', forces_zv_av_var, forces_nb)
  call object_alloc ('forces_zv_av_err', forces_zv_av_err, forces_nb)
  call object_alloc ('forces_zv_av_eloc_av_covar', forces_zv_av_eloc_av_covar, forces_nb)
  call object_alloc ('forces_zv_av_q_eloc_av_covar', forces_zv_av_q_eloc_av_covar, forces_nb)
  call object_alloc ('forces_zv_av_q_av_covar', forces_zv_av_q_av_covar, forces_nb)
  call object_alloc ('forces_zv_av_forces_pulay_av_covar', forces_zv_av_forces_pulay_av_covar, forces_nb)

  do force_i = 1, forces_nb
   cent_i = forces_cent (force_i)
   dim_i = forces_direct (force_i)

   forces_zv (force_i) = forces_nn (force_i)

!  loop over electrons
   do elec_i = 1, nelec

    dotproduct = 0.d0
    do dim_k = 1, ndim
     dotproduct = dotproduct + vec_en_xyz_wlk (dim_k, elec_i, cent_i, 1) * grd_psi_over_psi_wlk (dim_k, elec_i, 1)
    enddo ! dim_k

    forces_zv (force_i) = forces_zv (force_i) + znuc(iwctype(cent_i)) * (grd_psi_over_psi_wlk (dim_i, elec_i, 1) / dist_en_wlk (elec_i, cent_i, 1) - vec_en_xyz_wlk (dim_i, elec_i, cent_i, 1) * dotproduct / dist_en_wlk (elec_i, cent_i, 1)**3)
   enddo ! elec_i

  enddo ! force_i

  end subroutine forces_zv_bld

! ==============================================================================
  subroutine forces_zv_sq_bld
! ------------------------------------------------------------------------------
! Description   : square of forces_zv
!
! Created       : J. Toulouse, 24 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_sq')
   call object_average_define ('forces_zv_sq', 'forces_zv_sq_av')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv')

   return

  endif

! allocations
  call object_alloc ('forces_zv_sq', forces_zv_sq, forces_nb)
  call object_alloc ('forces_zv_sq_av', forces_zv_sq_av, forces_nb)

  forces_zv_sq (:) = forces_zv (:)**2

  end subroutine forces_zv_sq_bld

! ==============================================================================
  subroutine forces_zv_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of forces_zv
!
! Created       : J. Toulouse, 24 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_var')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv_sq_av')
   call object_needed ('forces_zv_av')

   return

  endif

! allocations
  call object_alloc ('forces_zv_var', forces_zv_var, forces_nb)

  forces_zv_var (:) = forces_zv_sq_av (:) - forces_zv_av (:)**2

  end subroutine forces_zv_var_bld

! ==============================================================================
  subroutine forces_zv_deriv_bld
! ------------------------------------------------------------------------------
! Description   : ZV estimator for forces using derivaive wrt nuclear coordinate
!
! Created       : J. Toulouse, 04 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv')
   call object_variance_define ('forces_zv_deriv_av', 'forces_zv_deriv_av_var')
   call object_average_define ('forces_zv_deriv', 'forces_zv_deriv_av')
   call object_error_define ('forces_zv_deriv_av', 'forces_zv_deriv_av_err')
   call object_covariance_define ('forces_zv_deriv_av','forces_pulay_av','forces_zv_deriv_av_forces_pulay_av_covar')

   call object_needed ('forces_nb')
!   call object_needed ('forces_bare')
   call object_needed ('deloc_rn_num')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv', forces_zv_deriv, forces_nb)
  call object_alloc ('forces_zv_deriv_av', forces_zv_deriv_av, forces_nb)
  call object_alloc ('forces_zv_deriv_av_var', forces_zv_deriv_av_var, forces_nb)
  call object_alloc ('forces_zv_deriv_av_err', forces_zv_deriv_av_err, forces_nb)
  call object_alloc ('forces_zv_deriv_av_forces_pulay_av_covar', forces_zv_deriv_av_forces_pulay_av_covar, forces_nb)

!  forces_zv_deriv (:) = forces_bare (:) - deloc_rn_num (:)
! deloc_rn_num includes forces_bare !
  forces_zv_deriv (:) = - deloc_rn_num (:)

  end subroutine forces_zv_deriv_bld

! ==============================================================================
  subroutine forces_zv_deriv_sq_bld
! ------------------------------------------------------------------------------
! Description   : square of forces_zv_deriv
!
! Created       : J. Toulouse, 05 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_sq')
   call object_average_define ('forces_zv_deriv_sq', 'forces_zv_deriv_sq_av')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv_deriv')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_sq', forces_zv_deriv_sq, forces_nb)
  call object_alloc ('forces_zv_deriv_sq_av', forces_zv_deriv_sq_av, forces_nb)

  forces_zv_deriv_sq (:) = forces_zv_deriv (:)**2

  end subroutine forces_zv_deriv_sq_bld

! ==============================================================================
  subroutine forces_zv_deriv_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of forces_zv_deriv
!
! Created       : J. Toulouse, 05 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_var')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv_deriv_sq_av')
   call object_needed ('forces_zv_deriv_av')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_var', forces_zv_deriv_var, forces_nb)

  forces_zv_deriv_var (:) = forces_zv_deriv_sq_av (:) - forces_zv_deriv_av (:)**2

  end subroutine forces_zv_deriv_var_bld

! ==============================================================================
  subroutine forces_zv_deriv_pulay_av_bld
! ------------------------------------------------------------------------------
! Description   : forces_zv_deriv_av + forces_pulay_av
!
! Created       : J. Toulouse, 06 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_pulay_av')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv_deriv_av')
   call object_needed ('forces_pulay_av')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_pulay_av', forces_zv_deriv_pulay_av, forces_nb)

  forces_zv_deriv_pulay_av (:) = forces_zv_deriv_av (:) + forces_pulay_av (:)

  end subroutine forces_zv_deriv_pulay_av_bld

! ==============================================================================
  subroutine forces_zv_deriv_pulay_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of forces_zv_deriv_pulay_av
!
! Created       : J. Toulouse, 06 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_pulay_av_var')
   call object_error_define_from_variance ('forces_zv_deriv_pulay_av_var', 'forces_zv_deriv_pulay_av_err')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv_deriv_av_var')
   call object_needed ('forces_pulay_av_var')
   call object_needed ('forces_zv_deriv_av_forces_pulay_av_covar')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_pulay_av_var', forces_zv_deriv_pulay_av_var, forces_nb)
  call object_alloc ('forces_zv_deriv_pulay_av_err', forces_zv_deriv_pulay_av_err, forces_nb)

  forces_zv_deriv_pulay_av_var (:) = forces_zv_deriv_av_var (:) + forces_pulay_av_var (:) + 2.d0 * forces_zv_deriv_av_forces_pulay_av_covar (:)

  end subroutine forces_zv_deriv_pulay_av_var_bld

! ==============================================================================
  subroutine forces_q_bld
! ------------------------------------------------------------------------------
! Description   : Simple Q term for forces
! Description   : Assaraf and Caffarel, J. Chem. Phys. (2003)
!
! Created       : J. Toulouse, 09 Aug 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, elec_i, cent_i, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('forces_q')
   call object_average_define ('forces_q', 'forces_q_av')
   call object_variance_define ('forces_q_av', 'forces_q_av_var')
   call object_error_define ('forces_q_av', 'forces_q_av_err')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv')
   call object_needed ('znuc')
   call object_needed ('iwctype')
   call object_needed ('dist_en_wlk')
   call object_needed ('vec_en_xyz_wlk')

   return

  endif

! allocations
  call object_alloc ('forces_q', forces_q, forces_nb)
  call object_alloc ('forces_q_av', forces_q_av, forces_nb)
  call object_alloc ('forces_q_av_var', forces_q_av_var, forces_nb)
  call object_alloc ('forces_q_av_err', forces_q_av_err, forces_nb)
  call object_alloc ('forces_q_av_eloc_av_covar', forces_q_av_eloc_av_covar, forces_nb)

  do force_i = 1, forces_nb
   cent_i = forces_cent (force_i)
   dim_i = forces_direct (force_i)

   forces_q (force_i) = 0.d0

!  loop over electrons
   do elec_i = 1, nelec
    forces_q (force_i) = forces_q (force_i) - znuc(iwctype(cent_i)) * vec_en_xyz_wlk (dim_i, elec_i, cent_i, 1) / dist_en_wlk (elec_i, cent_i, 1)
   enddo ! elec_i

  enddo ! force_i

  end subroutine forces_q_bld

! ==============================================================================
  subroutine forces_q_eloc_bld
! ------------------------------------------------------------------------------
! Description   : Simple Q E_L term for forces
! Description   : Assaraf and Caffarel, J. Chem. Phys. (2003)
!
! Created       : J. Toulouse, 19 Sep 2007
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_q_eloc')
   call object_average_define ('forces_q_eloc', 'forces_q_eloc_av')
   call object_variance_define ('forces_q_eloc_av', 'forces_q_eloc_av_var')
   call object_error_define ('forces_q_eloc_av', 'forces_q_eloc_av_err')

   call object_needed ('forces_nb')
   call object_needed ('forces_q')
   call object_needed ('eloc')

   return

  endif

! allocations
  call object_alloc ('forces_q_eloc', forces_q_eloc, forces_nb)
  call object_alloc ('forces_q_eloc_av', forces_q_eloc_av, forces_nb)
  call object_alloc ('forces_q_eloc_av_var', forces_q_eloc_av_var, forces_nb)
  call object_alloc ('forces_q_eloc_av_err', forces_q_eloc_av_err, forces_nb)
  call object_alloc ('forces_q_eloc_av_eloc_av_covar', forces_q_eloc_av_eloc_av_covar, forces_nb)
  call object_alloc ('forces_q_eloc_av_q_av_covar', forces_q_eloc_av_q_av_covar, forces_nb)

  forces_q_eloc (:) = forces_q (:) * eloc

  end subroutine forces_q_eloc_bld

! ==============================================================================
  subroutine forces_zvzb_av_bld
! ------------------------------------------------------------------------------
! Description   : Simple averaged ZVZB estimator for forces
! Description   : Assaraf and Caffarel, J. Chem. Phys. (2003)
!
! Created       : J. Toulouse, 19 Sep 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zvzb_av')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv_av')
   call object_needed ('forces_q_eloc_av')
   call object_needed ('forces_q_av')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('forces_zvzb_av', forces_zvzb_av, forces_nb)

  forces_zvzb_av (:) = forces_zv_av (:) + 2.d0 * (forces_q_eloc_av (:) - eloc_av * forces_q_av (:))

  end subroutine forces_zvzb_av_bld

! ==============================================================================
  subroutine forces_zvzb_av_var_bld
! ------------------------------------------------------------------------------
! Description   : Variance of averaged simple ZVZB estimator for forces
!
! Created       : J. Toulouse, 20 Sep 2007
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zvzb_av_var')
   call object_error_define_from_variance ('forces_zvzb_av_var', 'forces_zvzb_av_err')
   call object_covariance_define ('forces_zv_av', 'eloc_av', 'forces_zv_av_eloc_av_covar')
   call object_covariance_define ('forces_zv_av', 'forces_q_eloc_av', 'forces_zv_av_q_eloc_av_covar')
   call object_covariance_define ('forces_zv_av', 'forces_q_av', 'forces_zv_av_q_av_covar')
   call object_covariance_define ('forces_q_eloc_av', 'forces_q_av', 'forces_q_eloc_av_q_av_covar')
   call object_covariance_define ('forces_q_eloc_av', 'eloc_av', 'forces_q_eloc_av_eloc_av_covar')
   call object_covariance_define ('forces_q_av', 'eloc_av', 'forces_q_av_eloc_av_covar')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv_av_var')
   call object_needed ('forces_q_eloc_av_var')
   call object_needed ('forces_q_av')
   call object_needed ('forces_q_av_var')
   call object_needed ('eloc_av')
   call object_needed ('eloc_av_var')
   call object_needed ('forces_zv_av_q_eloc_av_covar')
   call object_needed ('forces_zv_av_q_av_covar')
   call object_needed ('forces_zv_av_eloc_av_covar')
   call object_needed ('forces_q_eloc_av_q_av_covar')
   call object_needed ('forces_q_eloc_av_eloc_av_covar')
   call object_needed ('forces_q_av_eloc_av_covar')

   return

  endif

! allocations
  call object_alloc ('forces_zvzb_av_var', forces_zvzb_av_var, forces_nb)
  call object_alloc ('forces_zvzb_av_err', forces_zvzb_av_err, forces_nb)

  forces_zvzb_av_var (:) = forces_zv_av_var (:) + 4.d0 * forces_q_eloc_av_var (:)                                  &
                          + 4.d0 * (eloc_av**2) * forces_q_av_var (:) + 4.d0 * (forces_q_av (:)**2) * eloc_av_var  &
                          + 4.d0 * forces_zv_av_q_eloc_av_covar (:) - 4.d0 * eloc_av * forces_zv_av_q_av_covar (:) &
                          - 4.d0 * forces_q_av (:) * forces_zv_av_eloc_av_covar (:)                                &
                          - 8.d0 * eloc_av * forces_q_eloc_av_q_av_covar (:)                                       &
                          - 8.d0 * forces_q_av (:) * forces_q_eloc_av_eloc_av_covar (:)                            &
                          + 8.d0 * forces_q_av (:) * eloc_av * forces_q_av_eloc_av_covar (:)

  end subroutine forces_zvzb_av_var_bld

! ==============================================================================
  subroutine forces_zvzb_old1_bld
! ------------------------------------------------------------------------------
! Description   : Simple ZVZB estimator for forces
! Description   : Assaraf and Caffarel, J. Chem. Phys. (2003)
!
! Created       : J. Toulouse, 27 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, cent_i, dim_i
  real(dp) eloc_c

! begin

! header
  if (header_exe) then

   call object_create ('forces_zvzb')
   call object_average_define ('forces_zvzb', 'forces_zvzb_av')
   call object_error_define ('forces_zvzb_av', 'forces_zvzb_av_err')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv')
   call object_needed ('forces_q')
   call object_needed ('znuc')
   call object_needed ('iwctype')
   call object_needed ('dist_en_wlk')
   call object_needed ('vec_en_xyz_wlk')
   call object_needed ('eloc')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('forces_zvzb', forces_zvzb, forces_nb)
  call object_alloc ('forces_zvzb_av', forces_zvzb_av, forces_nb)
  call object_alloc ('forces_zvzb_av_err', forces_zvzb_av_err, forces_nb)

  if (l_eloc_av_fixed) then
   eloc_c = eloc - eloc_av_fixed
  else
   eloc_c = eloc - eloc_av
  endif

  do force_i = 1, forces_nb
   cent_i = forces_cent (force_i)
   dim_i = forces_direct (force_i)

   forces_zvzb (force_i) = forces_zv (force_i)

   if (l_forces_q_av_fixed) then
    forces_zvzb (force_i) = forces_zvzb (force_i) + 2.d0*eloc_c*(forces_q (force_i) - forces_q_av_fixed (force_i))
   else
    forces_zvzb (force_i) = forces_zvzb (force_i) + 2.d0*eloc_c*(forces_q (force_i))
   endif

  enddo ! force_i

  end subroutine forces_zvzb_old1_bld

! ==============================================================================
  subroutine forces_zvzb_old2_bld
! ------------------------------------------------------------------------------
! Description   : Simple ZVZB estimator for forces
! Description   : Assaraf and Caffarel, J. Chem. Phys. (2003)
!
! Created       : J. Toulouse, 27 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, elec_i, cent_i, dim_i
  real(dp) eloc_c

! begin

! header
  if (header_exe) then

   call object_create ('forces_zvzb')
   call object_average_define ('forces_zvzb', 'forces_zvzb_av')
   call object_error_define ('forces_zvzb_av', 'forces_zvzb_av_err')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv')
   call object_needed ('znuc')
   call object_needed ('iwctype')
   call object_needed ('dist_en_wlk')
   call object_needed ('vec_en_xyz_wlk')
   call object_needed ('eloc')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('forces_zvzb', forces_zvzb, forces_nb)
  call object_alloc ('forces_zvzb_av', forces_zvzb_av, forces_nb)
  call object_alloc ('forces_zvzb_av_err', forces_zvzb_av_err, forces_nb)

  if (l_eloc_av_fixed) then
   eloc_c = eloc - eloc_av_fixed
  else
   eloc_c = eloc - eloc_av
  endif

  do force_i = 1, forces_nb
   cent_i = forces_cent (force_i)
   dim_i = forces_direct (force_i)

   forces_zvzb (force_i) = forces_zv (force_i)

!  loop over electrons
   do elec_i = 1, nelec
    forces_zvzb (force_i) = forces_zvzb (force_i) - 2.d0*eloc_c*znuc(iwctype(cent_i)) * vec_en_xyz_wlk (dim_i, elec_i, cent_i, 1) / dist_en_wlk (elec_i, cent_i, 1)
   enddo ! elec_i

  enddo ! force_i

  end subroutine forces_zvzb_old2_bld

! ==============================================================================
  subroutine forces_zv_pulay_av_bld
! ------------------------------------------------------------------------------
! Description   : simple averaged ZV + Pulay estimator for forces
!
! Created       : J. Toulouse, 31 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_pulay_av')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv_av')
   call object_needed ('forces_pulay_av')

   return

  endif

! allocations
  call object_alloc ('forces_zv_pulay_av', forces_zv_pulay_av, forces_nb)

  forces_zv_pulay_av (:) = forces_zv_av (:) + forces_pulay_av (:)

  end subroutine forces_zv_pulay_av_bld

! ==============================================================================
  subroutine forces_zv_pulay_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of forces_zv_pulay_av
!
! Created       : J. Toulouse, 31 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_pulay_av_var')
   call object_error_define_from_variance ('forces_zv_pulay_av_var', 'forces_zv_pulay_av_err')

   call object_needed ('forces_nb')
   call object_needed ('forces_zv_av_var')
   call object_needed ('forces_pulay_av_var')
   call object_needed ('forces_zv_av_forces_pulay_av_covar')

   return

  endif

! allocations
  call object_alloc ('forces_zv_pulay_av_var', forces_zv_pulay_av_var, forces_nb)
  call object_alloc ('forces_zv_pulay_av_err', forces_zv_pulay_av_err, forces_nb)

  forces_zv_pulay_av_var (:) = forces_zv_av_var (:) + forces_pulay_av_var (:) + 2.d0 * forces_zv_av_forces_pulay_av_covar (:)

  end subroutine forces_zv_pulay_av_var_bld

! ==============================================================================
  subroutine forces_zv_deloc_bld
! ------------------------------------------------------------------------------
! Description   : forces_zv * deloc
!
! Created       : J. Toulouse, 23 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i, force_i

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deloc')
   call object_average_define ('forces_zv_deloc', 'forces_zv_deloc_av')
   call object_covariance_define ('forces_zv_av', 'deloc_av', 'forces_zv_av_deloc_av_covar')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv')
   call object_needed ('deloc')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deloc', forces_zv_deloc, forces_nb, param_nb)
  call object_alloc ('forces_zv_deloc_av', forces_zv_deloc_av, forces_nb, param_nb)
  call object_alloc ('forces_zv_av_deloc_av_covar', forces_zv_av_deloc_av_covar, forces_nb, param_nb)

  do force_i = 1, forces_nb
   do param_i = 1, param_nb
    forces_zv_deloc (force_i, param_i) = forces_zv (force_i) * deloc (param_i)
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_deloc_bld

! ==============================================================================
  subroutine forces_zv_deloc_covar_bld
! ------------------------------------------------------------------------------
! Description   : covariance : < forces_zv * deloc > - < forces_zv > * < deloc >
!
! Created       : J. Toulouse, 23 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i, force_i

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deloc_covar')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_deloc_av')
   call object_needed ('forces_zv_av')
   call object_needed ('deloc_av')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deloc_covar', forces_zv_deloc_covar, forces_nb, param_nb)

  do force_i = 1, forces_nb
   do param_i = 1, param_nb
    forces_zv_deloc_covar (force_i, param_i) = forces_zv_deloc_av (force_i, param_i) - forces_zv_av (force_i) * deloc_av (param_i)
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_deloc_covar_bld

! ==============================================================================
  subroutine forces_zv_linear_coef_bld
! ------------------------------------------------------------------------------
! Description   : coefficient minimizing the variance of zero-variance estimator of total dipole moment
!
! Created       : J. Toulouse, 23 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer force_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_linear_coef')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_deloc_covar')
   call object_needed ('deloc_deloc_covar_inv')

   return

  endif

! allocations
  call object_alloc ('forces_zv_linear_coef', forces_zv_linear_coef, forces_nb, param_nb)

  do force_i = 1, forces_nb
   do param_i = 1, param_nb
    forces_zv_linear_coef (force_i, param_i) = 0.d0
    do param_j = 1, param_nb
     forces_zv_linear_coef (force_i, param_i) = forces_zv_linear_coef (force_i, param_i) - deloc_deloc_covar_inv (param_i, param_j) * forces_zv_deloc_covar (force_i, param_j)
    enddo ! param_j
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_linear_coef_bld

! ==============================================================================
  subroutine forces_zv_linear_av_bld
! ------------------------------------------------------------------------------
! Description   : average of ZV estimator for forces using wave function derivatives wrt parameters
!
! Created       : J. Toulouse, 12 Feb 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer force_i, param_i

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_linear_av')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_linear_coef')
   call object_needed ('forces_zv_av')
   call object_needed ('deloc_av')

   return

  endif

! allocations
  call object_alloc ('forces_zv_linear_av', forces_zv_linear_av, forces_nb)

  do force_i = 1, forces_nb
   forces_zv_linear_av (force_i) = forces_zv_av (force_i)
   do param_i = 1, param_nb
    forces_zv_linear_av (force_i) = forces_zv_linear_av (force_i) + forces_zv_linear_coef (force_i, param_i) * deloc_av (param_i)
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_linear_av_bld

! ==============================================================================
  subroutine forces_zv_linear_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of zero-variance estimator of force
!
! Created       : J. Toulouse, 24 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer force_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_linear_var')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_var')
   call object_needed ('forces_zv_linear_coef')
   call object_needed ('forces_zv_deloc_covar')
   call object_needed ('deloc_deloc_covar')

   return

  endif

! allocations
  call object_alloc ('forces_zv_linear_var', forces_zv_linear_var, forces_nb)

  do force_i = 1, forces_nb
   forces_zv_linear_var (force_i) = forces_zv_var (force_i)
   do param_i = 1, param_nb
    forces_zv_linear_var (force_i) = forces_zv_linear_var (force_i) + 2.d0 * forces_zv_linear_coef (force_i, param_i) * forces_zv_deloc_covar (force_i, param_i)
     do param_j = 1, param_nb
      forces_zv_linear_var (force_i) = forces_zv_linear_var (force_i) + forces_zv_linear_coef (force_i, param_i) * forces_zv_linear_coef (force_i, param_j) * deloc_deloc_covar (param_i, param_j)
     enddo ! param_j
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_linear_var_bld

! ==============================================================================
  subroutine forces_zv_linear_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of average of zero-variance estimator of force
!
! Created       : J. Toulouse, 24 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer force_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_linear_av_var')
   call object_error_define_from_variance ('forces_zv_linear_av_var', 'forces_zv_linear_av_err')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_av_var')
   call object_needed ('forces_zv_linear_coef')
   call object_needed ('forces_zv_av_deloc_av_covar')
   call object_needed ('deloc_av_deloc_av_covar')

   return

  endif

! allocations
  call object_alloc ('forces_zv_linear_av_var', forces_zv_linear_av_var, forces_nb)
  call object_alloc ('forces_zv_linear_av_err', forces_zv_linear_av_err, forces_nb)

  do force_i = 1, forces_nb
   forces_zv_linear_av_var (force_i) = forces_zv_av_var (force_i)
   do param_i = 1, param_nb
    forces_zv_linear_av_var (force_i) = forces_zv_linear_av_var (force_i) + 2.d0 * forces_zv_linear_coef (force_i, param_i) * forces_zv_av_deloc_av_covar (force_i, param_i)
     do param_j = 1, param_nb
      forces_zv_linear_av_var (force_i) = forces_zv_linear_av_var (force_i) + forces_zv_linear_coef (force_i, param_i) * forces_zv_linear_coef (force_i, param_j) * deloc_av_deloc_av_covar (param_i, param_j)
     enddo ! param_j
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_linear_av_var_bld

! ==============================================================================
  subroutine forces_zv_deriv_deloc_bld
! ------------------------------------------------------------------------------
! Description   : forces_zv_deriv * deloc
!
! Created       : J. Toulouse, 10 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i, force_i

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_deloc')
   call object_average_define ('forces_zv_deriv_deloc', 'forces_zv_deriv_deloc_av')
   call object_covariance_define ('forces_zv_deriv_av', 'deloc_av', 'forces_zv_deriv_av_deloc_av_covar')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_deriv')
   call object_needed ('deloc')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_deloc', forces_zv_deriv_deloc, forces_nb, param_nb)
  call object_alloc ('forces_zv_deriv_deloc_av', forces_zv_deriv_deloc_av, forces_nb, param_nb)
  call object_alloc ('forces_zv_deriv_av_deloc_av_covar', forces_zv_deriv_av_deloc_av_covar, forces_nb, param_nb)

  do force_i = 1, forces_nb
   do param_i = 1, param_nb
    forces_zv_deriv_deloc (force_i, param_i) = forces_zv_deriv (force_i) * deloc (param_i)
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_deriv_deloc_bld

! ==============================================================================
  subroutine forces_zv_deriv_deloc_covar_bld
! ------------------------------------------------------------------------------
! Description   : covariance : < forces_zv_deriv * deloc > - < forces_zv_deriv > * < deloc >
!
! Created       : J. Toulouse, 10 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i, force_i

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_deloc_covar')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_deriv_deloc_av')
   call object_needed ('forces_zv_deriv_av')
   call object_needed ('deloc_av')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_deloc_covar', forces_zv_deriv_deloc_covar, forces_nb, param_nb)

  do force_i = 1, forces_nb
   do param_i = 1, param_nb
    forces_zv_deriv_deloc_covar (force_i, param_i) = forces_zv_deriv_deloc_av (force_i, param_i) - forces_zv_deriv_av (force_i) * deloc_av (param_i)
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_deriv_deloc_covar_bld

! ==============================================================================
  subroutine forces_zv_deriv_linear_coef_bld
! ------------------------------------------------------------------------------
! Description   : coefficient minimizing the variance of zero-variance estimator
!
! Created       : J. Toulouse, 10 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer force_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_linear_coef')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_deriv_deloc_covar')
   call object_needed ('deloc_deloc_covar_inv')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_linear_coef', forces_zv_deriv_linear_coef, forces_nb, param_nb)

  do force_i = 1, forces_nb
   do param_i = 1, param_nb
    forces_zv_deriv_linear_coef (force_i, param_i) = 0.d0
    do param_j = 1, param_nb
     forces_zv_deriv_linear_coef (force_i, param_i) = forces_zv_deriv_linear_coef (force_i, param_i) - deloc_deloc_covar_inv (param_i, param_j) * forces_zv_deriv_deloc_covar (force_i, param_j)
    enddo ! param_j
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_deriv_linear_coef_bld

! ==============================================================================
  subroutine forces_zv_deriv_linear_av_bld
! ------------------------------------------------------------------------------
! Description   : average of ZV estimator for forces using wave function derivatives
! Description   : wrt to nuclear coordinates and parameters
!
! Created       : J. Toulouse, 10 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer force_i, param_i

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_linear_av')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_deriv_linear_coef')
   call object_needed ('forces_zv_deriv_av')
   call object_needed ('deloc_av')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_linear_av', forces_zv_deriv_linear_av, forces_nb)

  do force_i = 1, forces_nb
   forces_zv_deriv_linear_av (force_i) = forces_zv_deriv_av (force_i)
   do param_i = 1, param_nb
    forces_zv_deriv_linear_av (force_i) = forces_zv_deriv_linear_av (force_i) + forces_zv_deriv_linear_coef (force_i, param_i) * deloc_av (param_i)
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_deriv_linear_av_bld

! ==============================================================================
  subroutine forces_zv_deriv_linear_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of zero-variance estimator of force
!
! Created       : J. Toulouse, 10 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer force_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_linear_var')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_deriv_var')
   call object_needed ('forces_zv_deriv_linear_coef')
   call object_needed ('forces_zv_deriv_deloc_covar')
   call object_needed ('deloc_deloc_covar')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_linear_var', forces_zv_deriv_linear_var, forces_nb)

  do force_i = 1, forces_nb
   forces_zv_deriv_linear_var (force_i) = forces_zv_deriv_var (force_i)
   do param_i = 1, param_nb
    forces_zv_deriv_linear_var (force_i) = forces_zv_deriv_linear_var (force_i) + 2.d0 * forces_zv_deriv_linear_coef (force_i, param_i) * forces_zv_deriv_deloc_covar (force_i, param_i)
     do param_j = 1, param_nb
      forces_zv_deriv_linear_var (force_i) = forces_zv_deriv_linear_var (force_i) + forces_zv_deriv_linear_coef (force_i, param_i) * forces_zv_deriv_linear_coef (force_i, param_j) * deloc_deloc_covar (param_i, param_j)
     enddo ! param_j
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_deriv_linear_var_bld

! ==============================================================================
  subroutine forces_zv_deriv_linear_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of average of zero-variance estimator of force
!
! Created       : J. Toulouse, 10 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer force_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('forces_zv_deriv_linear_av_var')
   call object_error_define_from_variance ('forces_zv_deriv_linear_av_var', 'forces_zv_deriv_linear_av_err')

   call object_needed ('forces_nb')
   call object_needed ('param_nb')
   call object_needed ('forces_zv_deriv_av_var')
   call object_needed ('forces_zv_deriv_linear_coef')
   call object_needed ('forces_zv_deriv_av_deloc_av_covar')
   call object_needed ('deloc_av_deloc_av_covar')

   return

  endif

! allocations
  call object_alloc ('forces_zv_deriv_linear_av_var', forces_zv_deriv_linear_av_var, forces_nb)
  call object_alloc ('forces_zv_deriv_linear_av_err', forces_zv_deriv_linear_av_err, forces_nb)

  do force_i = 1, forces_nb
   forces_zv_deriv_linear_av_var (force_i) = forces_zv_deriv_av_var (force_i)
   do param_i = 1, param_nb
    forces_zv_deriv_linear_av_var (force_i) = forces_zv_deriv_linear_av_var (force_i) + 2.d0 * forces_zv_deriv_linear_coef (force_i, param_i) * forces_zv_deriv_av_deloc_av_covar (force_i, param_i)
     do param_j = 1, param_nb
      forces_zv_deriv_linear_av_var (force_i) = forces_zv_deriv_linear_av_var (force_i) + forces_zv_deriv_linear_coef (force_i, param_i) * forces_zv_deriv_linear_coef (force_i, param_j) * deloc_av_deloc_av_covar (param_i, param_j)
     enddo ! param_j
   enddo ! param_i
  enddo ! force_i

  end subroutine forces_zv_deriv_linear_av_var_bld
!===========================================================================
  subroutine forces_wrt
!---------------------------------------------------------------------------
! Description : write forces
!
! Created     : J. Toulouse, 23 Jul 2008
!---------------------------------------------------------------------------
  implicit none

! local
  integer force_i

! begin
  write(6,*)
  write(6,'(a)') 'Forces results:'

  write(6,'(a)') 'Nuclei-nuclei contribution to force:'
  call object_provide ('forces_nb')
  call object_provide ('forces_list')
  call object_provide ('forces_nn')
  do force_i = 1, forces_nb
   write(6,'(a,a4,a,f15.8)') 'component # ',forces_list (force_i),' : ', forces_nn (force_i)
  enddo

  if (l_forces_bare) then
   write(6,*)
   write(6,'(a)') 'Total force with bare (Hellmann-Feynman) estimator:'
   call object_provide ('forces_nb')
   call object_provide ('forces_list')
   call object_provide ('forces_bare_av')
   call object_provide ('forces_bare_av_err')
   do force_i = 1, forces_nb
    write(6,'(a,a4,a,f15.8,a,f15.8)') 'component # ',forces_list (force_i),' : ', forces_bare_av (force_i), ' +-', forces_bare_av_err (force_i)
   enddo
  endif

  if (l_forces_zv) then
   write(6,*)
   write(6,'(a)') 'Total force with simplest zero-variance estimator:'
   call object_provide ('forces_nb')
   call object_provide ('forces_list')
   call object_provide ('forces_zv_av')
   call object_provide ('forces_zv_av_err')
   call object_provide ('forces_zv_var')
   do force_i = 1, forces_nb
    write(6,'(a,a4,a,f15.8,a,f15.8,a,f15.8,a)') 'component # ',forces_list (force_i),' : ', forces_zv_av (force_i), ' +-', forces_zv_av_err (force_i),' (variance =',forces_zv_var (force_i),')'
   enddo
  endif

  if (l_forces_zv_deriv) then
   write(6,*)
   write(6,'(a)') 'Total force with zero-variance estimator using wave function derivatives wrt nuclear coordinates:'
   call object_provide ('forces_nb')
   call object_provide ('forces_list')
   call object_provide ('forces_zv_deriv_av')
   call object_provide ('forces_zv_deriv_av_err')
   call object_provide ('forces_zv_deriv_var')
   do force_i = 1, forces_nb
    write(6,'(a,a4,a,f15.8,a,f15.8,a,f15.8,a)') 'component # ',forces_list (force_i),' : ', forces_zv_deriv_av (force_i), ' +-', forces_zv_deriv_av_err (force_i),' (variance =',forces_zv_deriv_var (force_i),')'
   enddo
  endif

  if (l_forces_zv_linear) then
   write(6,*)
   write(6,'(a)') 'Total force with zero-variance estimator using wave function derivatives wrt parameters:'
   call object_provide ('forces_nb')
   call object_provide ('forces_list')
   call object_provide ('forces_zv_linear_av')
   call object_provide ('forces_zv_linear_av_err')
   call object_provide ('forces_zv_linear_var')
   do force_i = 1, forces_nb
    write(6,'(a,a4,a,f15.8,a,f15.8,a,f15.8,a)') 'component # ',forces_list (force_i),' : ', forces_zv_linear_av (force_i), ' +-', forces_zv_linear_av_err (force_i),' (variance =',forces_zv_linear_var (force_i),')'
   enddo
  endif

  if (l_forces_zv_deriv_linear) then
   write(6,*)
   write(6,'(a)') 'Total force with zero-variance estimator using wave function derivatives wrt nuclear coordinates and parameters:'
   call object_provide ('forces_nb')
   call object_provide ('forces_list')
   call object_provide ('forces_zv_deriv_linear_av')
   call object_provide ('forces_zv_deriv_linear_av_err')
   call object_provide ('forces_zv_deriv_linear_var')
   do force_i = 1, forces_nb
    write(6,'(a,a4,a,f15.8,a,f15.8,a,f15.8,a)') 'component # ',forces_list (force_i),' : ', forces_zv_deriv_linear_av (force_i), ' +-', forces_zv_deriv_linear_av_err (force_i),' (variance =',forces_zv_deriv_linear_var (force_i),')'
   enddo
  endif

  if (l_forces_pulay) then
   write(6,*)
   write(6,'(a)') 'Pulay contribution to force:'
   call object_provide ('forces_nb')
   call object_provide ('forces_list')
   call object_provide ('forces_pulay_av')
   call object_provide ('forces_pulay_av_err')
   do force_i = 1, forces_nb
    write(6,'(a,a4,a,f15.8,a,f15.8)') 'component # ',forces_list (force_i),' : ', forces_pulay_av (force_i), ' +-', forces_pulay_av_err (force_i)
   enddo
  endif

  if (l_forces_zvzb) then
   write(6,*)
   write(6,'(a)') 'Total force with simplest zero-variance zero-bias estimator:'
   call object_provide ('forces_nb')
   call object_provide ('forces_list')
   call object_provide ('forces_zvzb_av')
   call object_provide ('forces_zvzb_av_err')
   do force_i = 1, forces_nb
    write(6,'(a,a4,a,f15.8,a,f15.8)') 'component # ',forces_list (force_i),' : ', forces_zvzb_av (force_i), ' +-', forces_zvzb_av_err (force_i)
   enddo
  endif

  if (l_forces_zv_pulay) then
   write(6,*)
   write(6,'(a)') 'Total force with simplest zero-variance estimator + Pulay contribution:'
   call object_provide ('forces_nb')
   call object_provide ('forces_list')
   call object_provide ('forces_zv_pulay_av')
   call object_provide ('forces_zv_pulay_av_err')
   do force_i = 1, forces_nb
    write(6,'(a,a4,a,f15.8,a,f15.8)') 'component # ',forces_list (force_i),' : ', forces_zv_pulay_av (force_i), ' +-', forces_zv_pulay_av_err (force_i)
   enddo
  endif

  if (l_forces_zv_deriv_pulay) then
   write(6,*)
   write(6,'(a)') 'Total force with zero-variance estimator using wave function derivatives wrt nuclear coordinates + Pulay contribution:'
   call object_provide ('forces_nb')
   call object_provide ('forces_list')
   call object_provide ('forces_zv_deriv_pulay_av')
   call object_provide ('forces_zv_deriv_pulay_av_err')
   do force_i = 1, forces_nb
    write(6,'(a,a4,a,f15.8,a,f15.8)') 'component # ',forces_list (force_i),' : ', forces_zv_deriv_pulay_av (force_i), ' +-', forces_zv_deriv_pulay_av_err (force_i)
   enddo
  endif

  end subroutine forces_wrt

end module forces_mod
