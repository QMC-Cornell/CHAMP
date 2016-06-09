module deriv_mod

  use all_tools_mod
  use deriv_jas_mod
  use deriv_csf_mod
  use deriv_orb_mod
  use deriv_exp_mod
  use pjas_mod, only: param_pjas_nb
  use periodic_jastrow_mod
  use deriv_geo_mod
  use eloc_mod

! Declaration of global variables and default values
  integer                        :: param_aug_nb
  integer                        :: param_pairs_nb
  integer, allocatable           :: param_pairs (:,:)
  character(max_string_len), allocatable :: param_type (:)
  logical, allocatable           ::  is_param_type_csf (:)
  logical, allocatable           ::  is_param_type_jas (:)
  logical, allocatable           ::  is_param_type_pjas (:)
  logical, allocatable           ::  is_param_type_exp (:)
  logical, allocatable           ::  is_param_type_orb (:)
  logical, allocatable           ::  is_param_type_geo (:)

  real(dp), allocatable          :: dpsi (:)
  real(dp), allocatable          :: dpsi_bav (:)
  real(dp), allocatable          :: dpsi_av (:)
  real(dp), allocatable          :: dpsi_av_var (:)
  real(dp), allocatable          :: dpsi_uwav (:)
  real(dp), allocatable          :: dpsi_sq (:)
  real(dp), allocatable          :: dpsi_sq_av (:)
  real(dp), allocatable          :: dpsi_var (:)
  real(dp), allocatable          :: dpsi_sq_covar (:)
  real(dp), allocatable          :: dpsi_dpsi (:)
  real(dp), allocatable          :: dpsi_dpsi_av (:)
  real(dp), allocatable          :: dpsi_dpsi_uwav (:)
  real(dp), allocatable          :: dpsi_dpsi_covar (:,:)
  real(dp), allocatable          :: dpsi_dpsi_uwcovar (:,:)
  real(dp), allocatable          :: dpsi_dpsi_covar_inv (:,:)
  real(dp), allocatable          :: dpsi_dpsi_uwcovar_inv (:,:)
  real(dp), allocatable          :: dpsi_sq_eloc   (:)
  real(dp), allocatable          :: dpsi_sq_eloc_av(:)
  real(dp), allocatable          :: dpsi_dpsi_eloc (:)
  real(dp), allocatable          :: dpsi_dpsi_eloc_av (:)
  real(dp), allocatable          :: dpsi_dpsi_c_eloc_av (:,:)
  real(dp), allocatable          :: dpsi_dpsi_eloc_sq (:)
  real(dp), allocatable          :: dpsi_dpsi_eloc_sq_av (:)
  real(dp), allocatable          :: dpsi_eloc (:)
  real(dp), allocatable          :: dpsi_eloc_bav (:)
  real(dp), allocatable          :: dpsi_eloc_av (:)
  real(dp), allocatable          :: dpsi_eloc_av_var (:)
  real(dp), allocatable          :: dpsi_eloc_covar (:)
  real(dp), allocatable          :: dpsi_eloc_covar_var (:)
  real(dp), allocatable          :: dpsi_eloc_covar_err (:)
  real(dp), allocatable          :: dpsi_eloc_covar_err_2 (:)
  real(dp), allocatable          :: dpsi_eloc_covar_deloc_av (:)
  real(dp), allocatable          :: dpsi_eloc_covar_deloc_av_var (:)
  real(dp), allocatable          :: dpsi_eloc_covar_deloc_av_err (:)
  real(dp), allocatable          :: dpsi_eloc_covar_deloc_av_err_2 (:)
  real(dp), allocatable          :: dpsi_eloc_blk_covar (:)
  real(dp), allocatable          :: dpsi_eloc_sq (:)
  real(dp), allocatable          :: dpsi_eloc_sq_av (:)
  real(dp), allocatable          :: dpsi_eloc_sq_covar (:)
  real(dp), allocatable          :: dpsi_deloc (:,:)
  real(dp), allocatable          :: dpsi_deloc_av (:,:)
  real(dp), allocatable          :: dpsi_deloc_covar (:,:)
  real(dp), allocatable          :: dpsi_deloc_eloc (:,:)
  real(dp), allocatable          :: dpsi_deloc_eloc_av (:,:)
  real(dp), allocatable          :: dpsi_deloc_eloc_tc_av (:,:)
  real(dp), allocatable          :: dpsi_dpsi_eloc_eloc_qc_av (:,:)
  real(dp), allocatable          :: dpsi_eloc_av_dpsi_eloc_av_covar (:,:)
  real(dp), allocatable          :: dpsi_eloc_av_dpsi_av_covar (:,:)
  real(dp), allocatable          :: dpsi_eloc_av_eloc_av_covar (:)
  real(dp), allocatable          :: dpsi_av_dpsi_av_covar (:,:)
  real(dp), allocatable          :: dpsi_av_eloc_av_covar (:)
  real(dp), allocatable          :: deloc (:)
  real(dp), allocatable          :: deloc_av (:)
  real(dp), allocatable          :: deloc_av_var (:)
  real(dp), allocatable          :: deloc_av_err (:)
  real(dp), allocatable          :: deloc_bav (:)
  real(dp), allocatable          :: deloc_sq (:)
  real(dp), allocatable          :: deloc_sq_av (:)
  real(dp), allocatable          :: deloc_var (:)
  real(dp)                       :: deloc_av_abs_max
  real(dp), allocatable          :: deloc_eloc (:)
  real(dp), allocatable          :: deloc_eloc_av (:)
  real(dp), allocatable          :: deloc_eloc_covar (:)
  real(dp), allocatable          :: deloc_deloc (:)
  real(dp), allocatable          :: deloc_deloc_av (:)
  real(dp), allocatable          :: deloc_deloc_bav (:)
  real(dp), allocatable          :: deloc_av_deloc_av_covar (:,:)
  real(dp), allocatable          :: deloc_deloc_covar (:,:)
  real(dp), allocatable          :: deloc_deloc_blk_covar (:,:)
  real(dp), allocatable          :: deloc_deloc_covar_inv (:,:)
  real(dp), allocatable          :: deloc_deloc_blk_covar_inv (:,:)
  real(dp), allocatable          :: deloc_av_eloc_av_covar (:)
  real(dp), allocatable          :: deloc_av_dpsi_eloc_av_covar (:,:)
  real(dp), allocatable          :: deloc_av_dpsi_av_covar (:,:)
  real(dp), allocatable          :: d2psi (:)
  real(dp), allocatable          :: d2psi_av (:)
  real(dp), allocatable          :: d2psi_eloc (:)
  real(dp), allocatable          :: d2psi_eloc_av (:)

  real(dp), allocatable          :: gradient_energy (:)
  real(dp), allocatable          :: gradient_energy_blk (:)
  real(dp), allocatable          :: gradient_variance (:)
  real(dp), allocatable          :: gradient (:)
  real(dp)                       :: gradient_norm
  real(dp)                       :: gradient_norm_err

  character(len=max_string_len)  :: hessian_variance_type = 'levenberg_marquardt'
  real(dp), allocatable          :: hessian_variance_lm (:,:)
  real(dp), allocatable          :: hessian_variance_lmcov (:,:)
  real(dp), allocatable          :: hessian_variance_lin (:,:)
  real(dp), allocatable          :: hessian_variance (:,:)

  contains

! ==============================================================================
  subroutine param_nb_bld
! ------------------------------------------------------------------------------
! Description   : number of parameters to optimize
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j, pair
  integer shift

! header
  if (header_exe) then

   call object_create ('param_nb')
   call object_create ('param_aug_nb')
   call object_create ('param_pairs_nb')
   call object_create ('param_pairs')
   call object_create ('param_type')
   call object_create ('is_param_type_csf')
   call object_create ('is_param_type_jas')
   call object_create ('is_param_type_pjas')
   call object_create ('is_param_type_exp')
   call object_create ('is_param_type_orb')
   call object_create ('is_param_type_geo', is_param_type_geo_index)

   return

  endif

  call object_associate ('param_nb', param_nb)
  call object_associate ('param_pairs_nb', param_pairs_nb)

! begin
  param_nb = 0

  if (l_opt_csf) then
   call object_provide_by_index (param_nb_bld_index, nparmcsf_index)
   param_nb = param_nb + nparmcsf
  endif

  if (l_opt_jas) then
   call object_provide_by_index (param_nb_bld_index, nparmj_index)
   param_nb = param_nb + nparmj
  endif

  if (l_opt_pjas) then
   call object_provide_by_index (param_nb_bld_index, param_pjas_nb_index)
   param_nb = param_nb + param_pjas_nb
  endif

  if (l_opt_exp) then
   call object_provide_by_index (param_nb_bld_index, param_exp_nb_index)
   param_nb = param_nb + param_exp_nb
  endif

  if (l_opt_orb) then
   call object_provide_by_index (param_nb_bld_index, param_orb_nb_index)
   param_nb = param_nb + param_orb_nb
  endif


  if (l_opt_geo) then
   call object_provide_by_index (param_nb_bld_index, param_geo_nb_index)
   param_nb = param_nb + param_geo_nb
  endif

  call require (here, 'param_nb > 0', param_nb > 0)

  param_pairs_nb = param_nb * (param_nb + 1 ) / 2

  call object_alloc ('param_pairs', param_pairs, param_nb, param_nb)
  pair = 0
  do i = 1, param_nb
   do j = i, param_nb
      pair = pair + 1
      param_pairs (i,j) = pair
      param_pairs (j,i) = pair
   enddo
  enddo

! Type of each optimized parameter
  call object_alloc ('param_type', param_type, param_nb)
  call object_alloc ('is_param_type_csf', is_param_type_csf, param_nb)
  call object_alloc ('is_param_type_jas', is_param_type_jas, param_nb)
  call object_alloc ('is_param_type_pjas', is_param_type_pjas, param_nb)
  call object_alloc ('is_param_type_exp', is_param_type_exp, param_nb)
  call object_alloc ('is_param_type_orb', is_param_type_orb, param_nb)
  call object_alloc ('is_param_type_geo', is_param_type_geo, param_nb)
  is_param_type_csf (:) = .false.
  is_param_type_jas (:) = .false.
  is_param_type_pjas (:) = .false.
  is_param_type_exp (:) = .false.
  is_param_type_orb (:) = .false.
  is_param_type_geo (:) = .false.

  shift = 0
  param_type (shift+1:shift+nparmcsf) = 'CSF'
  is_param_type_csf (shift+1:shift+nparmcsf) = .true.
  shift = shift + nparmcsf

  param_type (shift+1:shift+nparmj) = 'Jastrow'
  is_param_type_jas (shift+1:shift+nparmj) = .true.
  shift = shift + nparmj

  param_type (shift+1:shift+param_pjas_nb) = 'periodic Jastrow'
  is_param_type_pjas (shift+1:shift+param_pjas_nb) = .true.
  shift = shift + param_pjas_nb

  param_type (shift+1:shift+param_exp_nb) = 'exponent'
  is_param_type_exp (shift+1:shift+param_exp_nb) = .true.
  shift = shift + param_exp_nb

  param_type (shift+1:shift+param_orb_nb) = 'orbitals'
  is_param_type_orb (shift+1:shift+param_orb_nb) = .true.
  shift = shift + param_orb_nb

  param_type (shift+1:shift+param_geo_nb) = 'geometry'
  is_param_type_geo (shift+1:shift+param_geo_nb) = .true.
  shift = shift + param_geo_nb

  param_aug_nb = param_nb + 1

  end subroutine param_nb_bld

! ==============================================================================
  subroutine dpsi_bld
! ------------------------------------------------------------------------------
! Description   :  Logarithmic derivatives of Psi with respect to
! Description   :  all parameters to optimize
! Description   :  d ln Psi / d c = (1/Psi) * d Psi / d c
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer shift, param_i

! header
  if (header_exe) then

   call object_create ('dpsi')
   call object_block_average_define ('dpsi', 'dpsi_bav')
   call object_average_define ('dpsi', 'dpsi_av')
   call object_variance_define ('dpsi_av', 'dpsi_av_var')
   call object_average_unweighted_define ('dpsi', 'dpsi_uwav')
   call object_covariance_define ('dpsi_av', 'dpsi_av', 'dpsi_av_dpsi_av_covar')
   call object_covariance_define ('dpsi_av', 'eloc_av', 'dpsi_av_eloc_av_covar')

   call object_needed ('param_nb')
   call object_needed ('is_param_type_jas') ! for bounds
   call object_needed ('is_param_type_pjas') ! for bounds

   return

  endif

! begin
  call object_alloc ('dpsi', dpsi, param_nb)
  call object_alloc ('dpsi_bav', dpsi_bav, param_nb)
  call object_alloc ('dpsi_av', dpsi_av, param_nb)
  call object_alloc ('dpsi_av_var', dpsi_av_var, param_nb)
  call object_alloc ('dpsi_uwav', dpsi_uwav, param_nb)
  call object_alloc ('dpsi_av_dpsi_av_covar', dpsi_av_dpsi_av_covar, param_nb, param_nb)
  call object_alloc ('dpsi_av_eloc_av_covar', dpsi_av_eloc_av_covar, param_nb)

  shift = 0

! CSFs contribution
  if (l_opt_csf) then
   call object_provide_by_index (dpsi_bld_index, dpsi_csf_index)
   dpsi (shift+1:shift+nparmcsf) = dpsi_csf (:)
   shift = shift + nparmcsf
  endif

! Jastrow contribution
  if (l_opt_jas) then
   call object_provide_by_index (dpsi_bld_index, dpsi_jas_index)
   dpsi (shift+1:shift+nparmj) = dpsi_jas (:)
   shift = shift + nparmj
  endif

! Periodic Jastrow contribution
  if (l_opt_pjas) then
     call object_provide_by_index (dpsi_bld_index, dpsi_pjas_index)
     dpsi (shift+1:shift+param_pjas_nb) = dpsi_pjas (:)
     shift = shift + param_pjas_nb
  endif

! exponents contribution
  if (l_opt_exp) then
   if (l_optimize_log_exp) then
    call object_provide_by_index (dpsi_bld_index, dpsi_lnexp_index)
    dpsi (shift+1:shift+param_exp_nb) = dpsi_lnexp (:)
   else
    call object_provide_by_index (dpsi_bld_index, dpsi_exp_index)
    dpsi (shift+1:shift+param_exp_nb) = dpsi_exp (:)
   endif
   shift = shift + param_exp_nb
  endif

! orbitals contribution
  if (l_opt_orb) then
   call object_provide_by_index (dpsi_bld_index, dpsi_orb_index)
   dpsi (shift+1:shift+param_orb_nb) = dpsi_orb (:)
   shift = shift + param_orb_nb
  endif

! geometry contribution
  if (l_opt_geo) then
   call object_provide_by_index (dpsi_bld_index, dpsi_geo_index)
   dpsi (shift+1:shift+param_geo_nb) = dpsi_geo (:)
   shift = shift + param_geo_nb
  endif

! applying bounds to dpsi for non-Jastrow parameters
!  if (l_deriv_bound) then
!   do param_i = 1, param_nb
!    if ((.not. is_param_type_jas (param_i) .and. .not. is_param_type_pjas (param_i)) .and. dabs(dpsi (param_i)) > deriv_bound_value) then
!      dpsi (param_i) = sign (deriv_bound_value, dpsi (param_i))
!    endif
!   enddo ! param_i
!  endif

  if (l_deriv_bound .and. block_iterations_nb > 1) then
   call object_provide_by_index (dpsi_av_index)
   call object_provide_by_index (dpsi_var_index)
   do param_i = 1, param_nb
    if ((.not. is_param_type_jas (param_i) .and. .not. is_param_type_pjas (param_i))  &
         .and. dabs(dpsi (param_i) - dpsi_av (param_i)) > deriv_bound_value * dsqrt(dpsi_var (param_i))) then
      if(dabs(dpsi (param_i) - dpsi_av (param_i)) > 10*dsqrt(dpsi_var (param_i))) write (6,'(a,i8,a,i4,4(a,f12.6))') 'step=', step_iterations_nb, ' param_i =', param_i, ' dpsi =',dpsi (param_i),' dpsi_av =',dpsi_av (param_i), ' dpsi_var =',dpsi_var (param_i), ' ratio=', dabs(dpsi (param_i) - dpsi_av (param_i)) / dsqrt(dpsi_var (param_i))
      dpsi (param_i) = dpsi_av (param_i) + sign (deriv_bound_value * dsqrt(dpsi_var (param_i)), dpsi (param_i) - dpsi_av (param_i))
!     write (6,'(a,i3,a,f12.6,a)') 'param_i =', param_i, 'dpsi =',dpsi (param_i),' after bound applied'
    endif
   enddo ! param_i
  endif

!  call require (lhere, 'shift == param_nb', shift == param_nb) !fp

  end subroutine dpsi_bld

! ==============================================================================
  subroutine d2psi_bld
! ------------------------------------------------------------------------------
! Description   :  2nd derivatives of Psi with respect to
! Description   :  all parameters to optimize
! Description   :  (1/Psi) * d^2 Psi / d c^2
!
! Created       : J. Toulouse, 19 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, iparmj_pair

! header
  if (header_exe) then

   call object_create ('d2psi')

   call object_needed ('param_pairs_nb')

   return

  endif

! begin
  call object_alloc ('d2psi', d2psi, param_pairs_nb)
  call object_alloc ('d2psi_av', d2psi_av, param_pairs_nb)

  d2psi = 0.d0

! zero second-order derivatives
  if (.not. l_deriv2nd) then
   return
  endif

  i = 0

! CSFs contribution
  if (l_opt_csf) then
    i = i + nparmcsf * (nparmcsf + 1) / 2
  endif

! Jastrow contribution
  if (l_opt_jas) then
   call object_provide_by_index (d2psi_bld_index, jas_pairs_nb_index)
   call object_provide_by_index (d2psi_bld_index, d2psi_jas_index)
   do iparmj_pair = 1, jas_pairs_nb
    i = i + 1
    d2psi (i) = d2psi_jas (iparmj_pair)
   enddo
  endif

! orbitals contribution
!  if (l_opt_exp) then
!    call die ('d2psi not yet implemented for exponents')
!  endif

! orbitals contribution
!  if (l_opt_orb) then
!    call die ('d2psi not yet implemented for orbitals')
!  endif

! geometry contribution
!  if (l_opt_geo) then
!    call die ('d2psi not yet implemented for geometry parameters')
!  endif

  end subroutine d2psi_bld

! ==============================================================================
  subroutine deloc_bld
! ------------------------------------------------------------------------------
! Description   :  derivative of local energy with respect to
! Description   :  all parameters to optimize
! Description   :  d eloc / d c
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer shift, param_i

! header
  if (header_exe) then

   call object_create ('deloc')
   call object_block_average_define ('deloc', 'deloc_bav')
   call object_average_define ('deloc', 'deloc_av')
   call object_variance_define ('deloc_av', 'deloc_av_var')
   call object_error_define ('deloc_av', 'deloc_av_err')
   call object_covariance_define ('deloc_av', 'eloc_av', 'deloc_av_eloc_av_covar')
   call object_covariance_define ('deloc_av', 'dpsi_eloc_av', 'deloc_av_dpsi_eloc_av_covar')
   call object_covariance_define ('deloc_av', 'dpsi_av', 'deloc_av_dpsi_av_covar')

   call object_needed ('param_nb')
   call object_needed ('is_param_type_jas') !for bounds
   call object_needed ('is_param_type_pjas') ! for bounds
   call object_needed ('dpsi') ! for bounds

   return

  endif

! begin
  call object_alloc ('deloc', deloc, param_nb)
  call object_alloc ('deloc_av', deloc_av, param_nb)
  call object_alloc ('deloc_av_var', deloc_av_var, param_nb) !!! temporary for testing DMC gradients
  call object_alloc ('deloc_av_err', deloc_av_err, param_nb) !!! temporary for testing DMC gradients
  call object_alloc ('deloc_bav', deloc_bav, param_nb)
  call object_alloc ('deloc_av_eloc_av_covar', deloc_av_eloc_av_covar, param_nb)
  call object_alloc ('deloc_av_dpsi_eloc_av_covar', deloc_av_dpsi_eloc_av_covar, param_nb, param_nb)
  call object_alloc ('deloc_av_dpsi_av_covar', deloc_av_dpsi_av_covar, param_nb, param_nb)

  shift = 0

! CSFs contribution
  if (l_opt_csf) then
   call object_provide_by_index (deloc_bld_index, deloc_csf_index)
   deloc (shift+1:shift+nparmcsf) = deloc_csf (:)
   shift = shift + nparmcsf
  endif

! Jastrow contribution
  if (l_opt_jas) then
   call object_provide_by_index (deloc_bld_index, deloc_jas_index)
   deloc (shift+1:shift+nparmj) = deloc_jas (:)
   shift = shift + nparmj
  endif

! Periodic Jastrow contribution
  if (l_opt_pjas) then
   call object_provide_by_index (deloc_bld_index, deloc_pjas_index)
   deloc (shift+1:shift+param_pjas_nb) = deloc_pjas (:)
   shift = shift + param_pjas_nb
  endif

! Exponent contribution
  if (l_opt_exp) then
   if (l_optimize_log_exp) then
    call object_provide_by_index (deloc_bld_index, deloc_lnexp_index)
    deloc (shift+1:shift+param_exp_nb) = deloc_lnexp (:)
   else
    call object_provide_by_index (deloc_bld_index, deloc_exp_index)
    deloc (shift+1:shift+param_exp_nb) = deloc_exp (:)
   endif
   shift = shift + param_exp_nb
  endif

! Orbitals contribution
  if (l_opt_orb .and. .not. l_opt_orb_eig) then
   call object_provide_by_index (deloc_bld_index, deloc_orb_index)
   deloc (shift+1:shift+param_orb_nb) = deloc_orb (:)
   shift = shift + param_orb_nb
  endif

! geometry contribution
  if (l_opt_geo) then
   call object_provide_by_index (deloc_bld_index, deloc_geo_index)
   deloc (shift+1:shift+param_geo_nb) = deloc_geo (:)
   shift = shift + param_geo_nb
  endif

!  call require (lhere, 'shift == param_nb', shift == param_nb) !fp

! applying bounds to deloc for non-Jastrow parameters
!  if (l_deriv_bound) then
!   do param_i = 1, param_nb
!    if ((.not. is_param_type_jas(param_i) .and. .not. is_param_type_pjas(param_i)) .and. dabs(dpsi(param_i)) > deriv_bound_value) then
!      deloc(param_i) = (deriv_bound_value/dabs(dpsi(param_i))) * deloc(param_i)
!    endif
!   enddo ! param_i
!  endif

  if (l_deriv_bound .and. block_iterations_nb > 1) then
   call object_provide_by_index (deloc_av_index)
   call object_provide_by_index (deloc_var_index)
   do param_i = 1, param_nb
    if ((.not. is_param_type_jas (param_i) .and. .not. is_param_type_pjas (param_i))  &
         .and. dabs(deloc (param_i) - deloc_av (param_i)) > deriv_bound_value * dsqrt(deloc_var (param_i))) then
      if(dabs(deloc (param_i) - deloc_av (param_i)) > 10*dsqrt(deloc_var (param_i))) write (6,'(a,i8,a,i3,4(a,f12.6))') 'step=', step_iterations_nb, ' param_i =', param_i, ' deloc =',deloc (param_i),' deloc_av =',deloc_av (param_i), ' deloc_var =',deloc_var (param_i), ' ratio=', dabs(deloc (param_i) - deloc_av (param_i)) / dsqrt(deloc_var (param_i))
      deloc (param_i) = deloc_av (param_i) + sign (deriv_bound_value * dsqrt(deloc_var (param_i)), deloc (param_i) - deloc_av (param_i))
!     write (6,'(a,i3,a,f12.6,a)') 'param_i =', param_i, 'deloc =',deloc (param_i),' after bound applied'
    endif
   enddo ! param_i
  endif

  end subroutine deloc_bld

! ==============================================================================
  subroutine deloc_av_abs_max_bld
! ------------------------------------------------------------------------------
! Description   : absolute maximal value of deloc_av
! Description   : must be zero in VMC with statistical noise because of Hermiticity of Hamiltonian
!
! Created       : J. Toulouse, 23 Jan 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i

! header
  if (header_exe) then

   call object_create ('deloc_av_abs_max')

   call object_needed ('deloc_av')
   call object_needed ('param_nb')

   return

  endif

! begin
  deloc_av_abs_max = 0.d0
  do param_i = 1, param_nb
   deloc_av_abs_max = max(abs(deloc_av (param_i)), deloc_av_abs_max)
  enddo

  end subroutine deloc_av_abs_max_bld

! ==============================================================================
  subroutine dpsi_sq_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_sq = dpsi (i) * dpsi (i)
!
! Created       : J. Toulouse, 17 Oct 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_sq')

   call object_needed ('param_nb')
   call object_needed ('dpsi')

   return

  endif

! begin
  call object_alloc ('dpsi_sq', dpsi_sq, param_nb)
  call object_alloc ('dpsi_sq_av', dpsi_sq_av, param_nb)

  dpsi_sq (:) = dpsi (:) * dpsi (:)

  end subroutine dpsi_sq_bld

! ==============================================================================
  subroutine dpsi_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of dpsi
!
! Created       : J. Toulouse, 11 Aug 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_var', dpsi_var_index)

   call object_needed ('param_nb')
   call object_needed ('dpsi_sq_av')
   call object_needed ('dpsi_av')

   return

  endif

! begin
  call object_alloc ('dpsi_var', dpsi_var, param_nb)

  dpsi_var (:) = dpsi_sq_av (:) - dpsi_av (:)**2

  end subroutine dpsi_var_bld

! ==============================================================================
  subroutine dpsi_sq_covar_bld
! ------------------------------------------------------------------------------
! Description   : the following covariance
! Description   : dpsi_sq_covar = < dpsi * dpsi >  - < dpsi > * < dpsi >
!
! Created       : J. Toulouse, 17 Oct 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_sq_covar')

   call object_needed ('param_nb')
   call object_needed ('dpsi_sq_av')
   call object_needed ('dpsi_av')

   return

  endif

! begin
  call object_alloc ('dpsi_sq_covar', dpsi_sq_covar, param_nb)

  dpsi_sq_covar (:) = dpsi_sq_av (:) - dpsi_av (:) * dpsi_av (:)

  end subroutine dpsi_sq_covar_bld

! ==============================================================================
  subroutine dpsi_dpsi_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_dpsi = dpsi (i) * dpsi (j)
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j

! header
  if (header_exe) then

   call object_create ('dpsi_dpsi')
   call object_average_define ('dpsi_dpsi', 'dpsi_dpsi_av')
   call object_average_unweighted_define ('dpsi_dpsi', 'dpsi_dpsi_uwav')

   call object_needed ('dpsi')
   call object_needed ('param_nb')
   call object_needed ('param_pairs_nb')
   call object_needed ('param_pairs')

   return

  endif

! begin

  call object_alloc ('dpsi_dpsi', dpsi_dpsi, param_pairs_nb)
  call object_alloc ('dpsi_dpsi_av', dpsi_dpsi_av, param_pairs_nb)
  call object_alloc ('dpsi_dpsi_uwav', dpsi_dpsi_uwav, param_pairs_nb)

  do i = 1, param_nb
   do j = i, param_nb
    dpsi_dpsi (param_pairs(i,j))= dpsi (i) *  dpsi (j)
   enddo
  enddo

  end subroutine dpsi_dpsi_bld

! ==============================================================================
  subroutine dpsi_dpsi_covar_bld
! ------------------------------------------------------------------------------
! Description   : the following covariance
! Description   : dpsi_dpsi_covar = < dpsi * dpsi > - <dpsi > * <dpsi >
!
! Created       : J. Toulouse, 27 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  integer i, j

! header
  if (header_exe) then

   call object_create ('dpsi_dpsi_covar')

   call object_needed ('param_pairs')
   call object_needed ('dpsi_av')
   call object_needed ('dpsi_dpsi_av')

   return

  endif

! begin
  call object_alloc ('dpsi_dpsi_covar', dpsi_dpsi_covar, param_nb, param_nb)

  do i = 1, param_nb
   do j = 1, param_nb
    dpsi_dpsi_covar (i,j)= dpsi_dpsi_av (param_pairs(i,j)) - dpsi_av (i) * dpsi_av (j)
   enddo
  enddo

  end subroutine dpsi_dpsi_covar_bld

! ==============================================================================
  subroutine dpsi_dpsi_covar_inv_bld
! ------------------------------------------------------------------------------
! Description   : inverse of dpsi_dpsi_covar
!
! Created       : J. Toulouse, 04 Feb 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  real(dp) threshold

! header
  if (header_exe) then

   call object_create ('dpsi_dpsi_covar_inv')

   call object_needed ('param_nb')
   call object_needed ('dpsi_dpsi_covar')

   return

  endif

! begin
  call object_alloc ('dpsi_dpsi_covar_inv', dpsi_dpsi_covar_inv, param_nb, param_nb)

  threshold = 1.d-10
  call inverse_by_svd (dpsi_dpsi_covar, dpsi_dpsi_covar_inv, param_nb, threshold)

  end subroutine dpsi_dpsi_covar_inv_bld

! ==============================================================================
  subroutine dpsi_dpsi_uwcovar_bld
! ------------------------------------------------------------------------------
! Description   : unweighted version of dpsi_dpsi_covar for DMC optimization
!
! Created       : J. Toulouse, 07 Jul 2010
! ------------------------------------------------------------------------------
  implicit none

! local
  integer i, j

! header
  if (header_exe) then

   call object_create ('dpsi_dpsi_uwcovar')

   call object_needed ('param_pairs')
   call object_needed ('dpsi_uwav')
   call object_needed ('dpsi_dpsi_uwav')

   return

  endif

! begin
  call object_alloc ('dpsi_dpsi_uwcovar', dpsi_dpsi_uwcovar, param_nb, param_nb)

  do i = 1, param_nb
   do j = 1, param_nb
    dpsi_dpsi_uwcovar (i,j)= dpsi_dpsi_uwav (param_pairs(i,j)) - dpsi_uwav (i) * dpsi_uwav (j)
   enddo
  enddo

  end subroutine dpsi_dpsi_uwcovar_bld

! ==============================================================================
  subroutine dpsi_dpsi_uwcovar_inv_bld
! ------------------------------------------------------------------------------
! Description   : inverse of dpsi_dpsi_uwcovar
!
! Created       : J. Toulouse, 07 Jul 2010
! ------------------------------------------------------------------------------
  implicit none

! local
  real(dp) threshold

! header
  if (header_exe) then

   call object_create ('dpsi_dpsi_uwcovar_inv')

   call object_needed ('param_nb')
   call object_needed ('dpsi_dpsi_uwcovar')

   return

  endif

! begin
  call object_alloc ('dpsi_dpsi_uwcovar_inv', dpsi_dpsi_uwcovar_inv, param_nb, param_nb)

  threshold = 1.d-10
  call inverse_by_svd (dpsi_dpsi_uwcovar, dpsi_dpsi_uwcovar_inv, param_nb, threshold)

  end subroutine dpsi_dpsi_uwcovar_inv_bld

! ==============================================================================
  subroutine dpsi_sq_eloc_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_sq_eloc (i)= dpsi_dpsi (i,i) * eloc
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  integer i

! header
  if (header_exe) then

   call object_create ('dpsi_sq_eloc')

   call object_needed ('param_nb')
   call object_needed ('dpsi_dpsi')
   call object_needed ('eloc')

   return

  endif

! begin
  call object_alloc ('dpsi_sq_eloc', dpsi_sq_eloc, param_nb)
  call object_alloc ('dpsi_sq_eloc_av', dpsi_sq_eloc_av, param_nb)

  do i= 1, param_nb
   dpsi_sq_eloc (i) = dpsi_dpsi (param_pairs(i,i)) * eloc
  enddo

  end subroutine dpsi_sq_eloc_bld

! ==============================================================================
  subroutine dpsi_dpsi_eloc_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_dpsi_eloc = dpsi_dpsi * eloc
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_dpsi_eloc')

   call object_needed ('param_pairs_nb')
   call object_needed ('dpsi_dpsi')
   call object_needed ('eloc')

   return

  endif

! begin

  call object_alloc ('dpsi_dpsi_eloc', dpsi_dpsi_eloc, param_pairs_nb)
  call object_alloc ('dpsi_dpsi_eloc_av', dpsi_dpsi_eloc_av, param_pairs_nb)

  dpsi_dpsi_eloc (:) = dpsi_dpsi (:) * eloc

  end subroutine dpsi_dpsi_eloc_bld

! ==============================================================================
  subroutine dpsi_dpsi_eloc_sq_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_dpsi_eloc_sq = dpsi_dpsi * eloc_sq
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_dpsi_eloc_sq')
   call object_average_define ('dpsi_dpsi_eloc_sq', 'dpsi_dpsi_eloc_sq_av')

   call object_needed ('param_pairs_nb')
   call object_needed ('dpsi_dpsi')
   call object_needed ('eloc_sq')

   return

  endif

! begin
  call object_alloc ('dpsi_dpsi_eloc_sq', dpsi_dpsi_eloc_sq, param_pairs_nb)
  call object_alloc ('dpsi_dpsi_eloc_sq_av', dpsi_dpsi_eloc_sq_av, param_pairs_nb)

  dpsi_dpsi_eloc_sq (:) = dpsi_dpsi (:) * eloc_sq

  end subroutine dpsi_dpsi_eloc_sq_bld

! ==============================================================================
  subroutine dpsi_eloc_bld
! ------------------------------------------------------------------------------
! Description   :  dpsi_eloc = dpsi  * eloc
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_eloc')

   call object_needed ('param_nb')
   call object_needed ('dpsi')
   call object_needed ('eloc')
   call object_block_average_define ('dpsi_eloc', 'dpsi_eloc_bav')
   call object_average_define ('dpsi_eloc', 'dpsi_eloc_av')
   call object_variance_define ('dpsi_eloc_av', 'dpsi_eloc_av_var')
   call object_covariance_define ('dpsi_eloc_av', 'dpsi_eloc_av', 'dpsi_eloc_av_dpsi_eloc_av_covar')
   call object_covariance_define ('dpsi_eloc_av', 'dpsi_av', 'dpsi_eloc_av_dpsi_av_covar')
   call object_covariance_define ('dpsi_eloc_av', 'eloc_av', 'dpsi_eloc_av_eloc_av_covar')

   return

  endif

! begin

  call object_alloc ('dpsi_eloc', dpsi_eloc, param_nb)
  call object_alloc ('dpsi_eloc_bav', dpsi_eloc_bav, param_nb)
  call object_alloc ('dpsi_eloc_av', dpsi_eloc_av, param_nb)
  call object_alloc ('dpsi_eloc_av_var', dpsi_eloc_av_var, param_nb)
  call object_alloc ('dpsi_eloc_av_dpsi_eloc_av_covar', dpsi_eloc_av_dpsi_eloc_av_covar, param_nb, param_nb)
  call object_alloc ('dpsi_eloc_av_dpsi_av_covar', dpsi_eloc_av_dpsi_av_covar, param_nb, param_nb)
  call object_alloc ('dpsi_eloc_av_eloc_av_covar', dpsi_eloc_av_eloc_av_covar, param_nb)

  dpsi_eloc = dpsi * eloc

  end subroutine dpsi_eloc_bld

! ==============================================================================
  subroutine dpsi_eloc_sq_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_eloc_sq = dpsi * eloc^2
!
! Created       : J. Toulouse, 24 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_eloc_sq')
   call object_average_define ('dpsi_eloc_sq', 'dpsi_eloc_sq_av')

   call object_needed ('param_nb')
   call object_needed ('dpsi')
   call object_needed ('eloc_sq')

   return

  endif

! begin
  call object_alloc ('dpsi_eloc_sq', dpsi_eloc_sq, param_nb)
  call object_alloc ('dpsi_eloc_sq_av', dpsi_eloc_sq_av, param_nb)

  dpsi_eloc_sq (:) = dpsi (:) * eloc_sq

  end subroutine dpsi_eloc_sq_bld

! ==============================================================================
  subroutine dpsi_eloc_sq_covar_bld
! ------------------------------------------------------------------------------
! Description   : the following covariance
! Description   : dpsi_eloc_sq_covar = < dpsi * eloc^2 > - < dpsi > * < eloc^2 >
!
! Created       : J. Toulouse, 27 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_eloc_sq_covar')

   call object_needed ('param_nb')
   call object_needed ('dpsi_eloc_sq_av')
   call object_needed ('dpsi_av')
   call object_needed ('eloc_sq_av')

   return

  endif

! begin
  call object_alloc ('dpsi_eloc_sq_covar', dpsi_eloc_sq_covar, param_nb)
  dpsi_eloc_sq_covar (:) = dpsi_eloc_sq_av (:) - dpsi_av (:) * eloc_sq_av

  end subroutine dpsi_eloc_sq_covar_bld

! ==============================================================================
  subroutine deloc_sq_bld
! ------------------------------------------------------------------------------
! Description   : deloc_sq = deloc (i) * deloc (i)
!
! Created       : J. Toulouse, 11 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! header
  if (header_exe) then

   call object_create ('deloc_sq')
   call object_average_define ('deloc_sq', 'deloc_sq_av')

   call object_needed ('param_nb')
   call object_needed ('deloc')

   return

  endif

! begin
  call object_alloc ('deloc_sq', deloc_sq, param_nb)
  call object_alloc ('deloc_sq_av', deloc_sq_av, param_nb)

  deloc_sq (:) = deloc (:) * deloc (:)

  end subroutine deloc_sq_bld

! ==============================================================================
  subroutine deloc_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of deloc
!
! Created       : J. Toulouse, 11 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! header
  if (header_exe) then

   call object_create ('deloc_var', deloc_var_index)

   call object_needed ('param_nb')
   call object_needed ('deloc_sq_av')
   call object_needed ('deloc_av')

   return

  endif

! begin
  call object_alloc ('deloc_var', deloc_var, param_nb)

  deloc_var (:) = deloc_sq_av (:) - deloc_av (:)**2

  end subroutine deloc_var_bld

! ==============================================================================
  subroutine deloc_eloc_bld
! ------------------------------------------------------------------------------
! Description   : deloc_eloc = deloc * eloc
!
! Created       : J. Toulouse, 24 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('deloc_eloc')
   call object_average_define ('deloc_eloc', 'deloc_eloc_av')

   call object_needed ('param_nb')
   call object_needed ('deloc')
   call object_needed ('eloc')

   return

  endif

! begin
  call object_alloc ('deloc_eloc', deloc_eloc, param_nb)
  call object_alloc ('deloc_eloc_av', deloc_eloc_av, param_nb)

  deloc_eloc = deloc * eloc

  end subroutine deloc_eloc_bld

! ==============================================================================
  subroutine deloc_eloc_covar_bld
! ------------------------------------------------------------------------------
! Description   : the following covariance
! Description   : deloc_eloc_covar = < deloc * eloc > - < deloc > * < eloc >
!
! Created       : J. Toulouse, 24 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('deloc_eloc_covar')

   call object_needed ('param_nb')
   call object_needed ('deloc_eloc_av')
   call object_needed ('deloc_av')
   call object_needed ('eloc_av')

   return

  endif

! begin
  call object_alloc ('deloc_eloc_covar', deloc_eloc_covar, param_nb)
  deloc_eloc_covar (:) = deloc_eloc_av (:) - deloc_av (:) * eloc_av

  end subroutine deloc_eloc_covar_bld

! ==============================================================================
  subroutine deloc_deloc_bld
! ------------------------------------------------------------------------------
! Description   : deloc_deloc = deloc * deloc
!
! Created       : J. Toulouse, 25 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('deloc_deloc')
   call object_block_average_define ('deloc_deloc', 'deloc_deloc_bav')
   call object_average_define ('deloc_deloc', 'deloc_deloc_av')
   call object_covariance_define ('deloc_av', 'deloc_av', 'deloc_av_deloc_av_covar')

   call object_needed ('param_nb')
   call object_needed ('param_pairs_nb')
   call object_needed ('param_pairs')
   call object_needed ('deloc')

   return

  endif

! begin
  call object_alloc ('deloc_deloc', deloc_deloc, param_pairs_nb)
  call object_alloc ('deloc_deloc_av', deloc_deloc_av, param_pairs_nb)
  call object_alloc ('deloc_deloc_bav', deloc_deloc_bav, param_pairs_nb)
  call object_alloc ('deloc_av_deloc_av_covar', deloc_av_deloc_av_covar, param_nb, param_nb)

  do param_i = 1, param_nb
   do param_j = param_i, param_nb
     deloc_deloc (param_pairs (param_i, param_j)) = deloc (param_i) * deloc (param_j)
   enddo
  enddo

  end subroutine deloc_deloc_bld

! ==============================================================================
  subroutine deloc_deloc_covar_bld
! ------------------------------------------------------------------------------
! Description   : the following covariance
! Description   : deloc_deloc_covar =  < deloc * deloc > - < deloc > * < deloc >
!
! Created       : J. Toulouse, 27 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('deloc_deloc_covar')


   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('deloc_deloc_av')
   call object_needed ('deloc_av')

   return

  endif

! begin
  call object_alloc ('deloc_deloc_covar', deloc_deloc_covar, param_nb, param_nb)

  do param_i = 1, param_nb
   do param_j = 1, param_nb
     deloc_deloc_covar (param_i, param_j) = deloc_deloc_av (param_pairs (param_i, param_j)) - deloc_av (param_i) * deloc_av (param_j)
   enddo
  enddo

  end subroutine deloc_deloc_covar_bld

! ==============================================================================
  subroutine deloc_deloc_blk_covar_bld
! ------------------------------------------------------------------------------
! Description   : deloc_deloc_covar calculated on current block
!
! Created       : J. Toulouse, 20 Apr 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('deloc_deloc_blk_covar')


   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('deloc_deloc_bav')
   call object_needed ('deloc_bav')

   return

  endif

! begin
  call object_alloc ('deloc_deloc_blk_covar', deloc_deloc_blk_covar, param_nb, param_nb)

  do param_i = 1, param_nb
   do param_j = 1, param_nb
     deloc_deloc_blk_covar (param_i, param_j) = deloc_deloc_bav (param_pairs (param_i, param_j)) - deloc_bav (param_i) * deloc_bav (param_j)
   enddo
  enddo

  end subroutine deloc_deloc_blk_covar_bld

! ==============================================================================
  subroutine deloc_deloc_covar_inv_bld
! ------------------------------------------------------------------------------
! Description   : inverse of deloc_deloc_covar
!
! Created       : J. Toulouse, 12 Feb 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  real(dp) threshold

! header
  if (header_exe) then

   call object_create ('deloc_deloc_covar_inv')

   call object_needed ('param_nb')
   call object_needed ('deloc_deloc_covar')

   return

  endif

! begin
  call object_alloc ('deloc_deloc_covar_inv', deloc_deloc_covar_inv, param_nb, param_nb)

  threshold = 1.d-10
!!  threshold = 1.d-8 !!!!!!
  call inverse_by_svd (deloc_deloc_covar, deloc_deloc_covar_inv, param_nb, threshold)

  end subroutine deloc_deloc_covar_inv_bld

! ==============================================================================
  subroutine deloc_deloc_blk_covar_inv_bld
! ------------------------------------------------------------------------------
! Description   : inverse of deloc_deloc_blk_covar
!
! Created       : J. Toulouse, 21 Apr 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  real(dp) threshold

! header
  if (header_exe) then

   call object_create ('deloc_deloc_blk_covar_inv')

   call object_needed ('param_nb')
   call object_needed ('deloc_deloc_blk_covar')

   return

  endif

! begin
  call object_alloc ('deloc_deloc_blk_covar_inv', deloc_deloc_blk_covar_inv, param_nb, param_nb)

  threshold = 1.d-10
  call inverse_by_svd (deloc_deloc_blk_covar, deloc_deloc_blk_covar_inv, param_nb, threshold)

  end subroutine deloc_deloc_blk_covar_inv_bld

! ==============================================================================
  subroutine dpsi_deloc_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_deloc = dpsi * deloc
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j

! header
  if (header_exe) then

   call object_create ('dpsi_deloc')

   call object_needed ('param_nb')
   call object_needed ('dpsi')
   call object_needed ('deloc')

   return

  endif

! begin
  call object_alloc ('dpsi_deloc', dpsi_deloc, param_nb, param_nb)
  call object_alloc ('dpsi_deloc_av', dpsi_deloc_av, param_nb, param_nb)

  do i = 1, param_nb
   do j = 1, param_nb
!     if (.not. (l_opt_orb_eig .and. j > nparmcsf+nparmj)) then
      dpsi_deloc (i,j) = dpsi (i) * deloc (j)
!     endif
   enddo
  enddo

  end subroutine dpsi_deloc_bld

! ==============================================================================
  subroutine dpsi_deloc_covar_bld
! ------------------------------------------------------------------------------
! Description   : the following covariance
! Description   : dpsi_deloc_covar = < dpsi * deloc > - < dpsi > * < deloc >
!
! Created       : J. Toulouse, 27 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('dpsi_deloc_covar')

   call object_needed ('param_nb')
   call object_needed ('dpsi_deloc_av')
   call object_needed ('dpsi_av')
   call object_needed ('deloc_av')

   return

  endif

! begin
  call object_alloc ('dpsi_deloc_covar', dpsi_deloc_covar, param_nb, param_nb)

  do param_i = 1, param_nb
   do param_j = 1, param_nb
     dpsi_deloc_covar (param_i, param_j) = dpsi_deloc_av (param_i, param_j) - dpsi_av (param_i) * deloc_av (param_j)
   enddo
  enddo

  end subroutine dpsi_deloc_covar_bld

! ==============================================================================
  subroutine d2psi_eloc_bld
! ------------------------------------------------------------------------------
! Description   : d2psi_eloc (i,j)= d2psi (i,j) * eloc
!
! Created       : J. Toulouse, 19 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none


! header
  if (header_exe) then

   call object_create ('d2psi_eloc')

   call object_needed ('param_pairs_nb')
   call object_needed ('d2psi')
   call object_needed ('eloc')

   return

  endif

! begin
  call object_alloc ('d2psi_eloc', d2psi_eloc, param_pairs_nb)
  call object_alloc ('d2psi_eloc_av', d2psi_eloc_av, param_pairs_nb)

  d2psi_eloc = d2psi * eloc

  end subroutine d2psi_eloc_bld

! ==============================================================================
  subroutine dpsi_eloc_covar_bld
! ------------------------------------------------------------------------------
! Description   : the following covariance
! Description   : <dpsi_eloc> - <dpsi> <eloc>
!
! Created       : J. Toulouse, 27 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_eloc_covar')
   call object_error_define ('dpsi_eloc_covar', 'dpsi_eloc_covar_err') ! temporary for testing DMC gradient

   call object_needed ('param_nb')
   call object_needed ('dpsi_eloc_av')
   call object_needed ('dpsi_av')
   call object_needed ('eloc_av')

   return

  endif

! begin
  call object_alloc ('dpsi_eloc_covar', dpsi_eloc_covar, param_nb)
  call object_alloc ('dpsi_eloc_covar_err', dpsi_eloc_covar_err, param_nb)

  dpsi_eloc_covar =  dpsi_eloc_av - dpsi_av * eloc_av

  end subroutine dpsi_eloc_covar_bld

! ==============================================================================
  subroutine dpsi_eloc_covar_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of dpsi_eloc_covar
! Description   : this is to calculate the error more rigourously than from block averages
! Description   : but it almost do not change the error anyway (really?)
!
! Created       : J. Toulouse, 24 Apr 2015
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i

! header
  if (header_exe) then

   call object_create ('dpsi_eloc_covar_var')
   call object_error_define_from_variance ('dpsi_eloc_covar_var', 'dpsi_eloc_covar_err_2')

   call object_needed ('param_nb')
   call object_needed ('dpsi_av')
   call object_needed ('eloc_av')
   call object_needed ('eloc_av_var')
   call object_needed ('dpsi_eloc_av_var')
   call object_needed ('dpsi_av_var')
   call object_needed ('dpsi_eloc_av_dpsi_av_covar')
   call object_needed ('dpsi_eloc_av_eloc_av_covar')
   call object_needed ('dpsi_av_eloc_av_covar')

   return

  endif

! begin
  call object_alloc ('dpsi_eloc_covar_var', dpsi_eloc_covar_var, param_nb)
  call object_alloc ('dpsi_eloc_covar_err_2', dpsi_eloc_covar_err_2, param_nb)

  do param_i = 1, param_nb
   dpsi_eloc_covar_var (param_i) = dpsi_eloc_av_var (param_i) + eloc_av**2 * dpsi_av_var (param_i) + dpsi_av (param_i)**2 * eloc_av_var  &
                                 - 2 * eloc_av * dpsi_eloc_av_dpsi_av_covar (param_i, param_i)    &
                                 - 2 * dpsi_av (param_i) * dpsi_eloc_av_eloc_av_covar (param_i)   &
                                 + 2 * dpsi_av (param_i) * eloc_av * dpsi_av_eloc_av_covar (param_i)
  enddo

  end subroutine dpsi_eloc_covar_var_bld


! ==============================================================================
  subroutine dpsi_eloc_blk_covar_bld
! ------------------------------------------------------------------------------
! Description   : covariance <dpsi_eloc> - <dpsi> <eloc> for current block
! Description   : <dpsi_eloc> - <dpsi> <eloc>
!
! Created       : J. Toulouse, 08 May 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_eloc_blk_covar')

   call object_needed ('param_nb')
   call object_needed ('dpsi_eloc_bav')
   call object_needed ('dpsi_bav')
   call object_needed ('eloc_bav')

   return

  endif

! begin
  call object_alloc ('dpsi_eloc_blk_covar', dpsi_eloc_blk_covar, param_nb)

  dpsi_eloc_blk_covar =  dpsi_eloc_bav - dpsi_bav * eloc_bav

  end subroutine dpsi_eloc_blk_covar_bld

! ==============================================================================
  subroutine dpsi_eloc_covar_deloc_av_bld
! ------------------------------------------------------------------------------
! Description   : 
!
! Created       : J. Toulouse, 27 Aug 2014
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_eloc_covar_deloc_av')
   call object_error_define ('dpsi_eloc_covar_deloc_av', 'dpsi_eloc_covar_deloc_av_err') ! temporary for testing DMC gradient

   call object_needed ('param_nb')
   call object_needed ('dpsi_eloc_covar')
   call object_needed ('deloc_av')

   return

  endif

! begin
  call object_alloc ('dpsi_eloc_covar_deloc_av', dpsi_eloc_covar_deloc_av, param_nb)
  call object_alloc ('dpsi_eloc_covar_deloc_av_err', dpsi_eloc_covar_deloc_av_err, param_nb)

  dpsi_eloc_covar_deloc_av =  dpsi_eloc_covar + deloc_av

  end subroutine dpsi_eloc_covar_deloc_av_bld


! ==============================================================================
  subroutine dpsi_eloc_covar_deloc_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of dpsi_eloc_covar_deloc_av
! Description   : this is to calculate the error more rigourously than from block averages
! Description   : but it almost do not change the error anyway
!
! Created       : J. Toulouse, 24 Apr 2015
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i

! header
  if (header_exe) then

   call object_create ('dpsi_eloc_covar_deloc_av_var')
   call object_error_define_from_variance ('dpsi_eloc_covar_deloc_av_var', 'dpsi_eloc_covar_deloc_av_err_2')

   call object_needed ('param_nb')
   call object_needed ('deloc_av_var')
   call object_needed ('dpsi_av')
   call object_needed ('eloc_av')
   call object_needed ('eloc_av_var')
   call object_needed ('dpsi_eloc_av_var')
   call object_needed ('dpsi_av_var')
   call object_needed ('dpsi_eloc_av_dpsi_av_covar')
   call object_needed ('dpsi_eloc_av_eloc_av_covar')
   call object_needed ('dpsi_av_eloc_av_covar')
   call object_needed ('deloc_av_dpsi_eloc_av_covar')
   call object_needed ('deloc_av_dpsi_av_covar')
   call object_needed ('deloc_av_eloc_av_covar')

   return

  endif

! begin
  call object_alloc ('dpsi_eloc_covar_deloc_av_var', dpsi_eloc_covar_deloc_av_var, param_nb)
  call object_alloc ('dpsi_eloc_covar_deloc_av_err_2', dpsi_eloc_covar_deloc_av_err_2, param_nb)

  do param_i = 1, param_nb
   dpsi_eloc_covar_deloc_av_var (param_i) = dpsi_eloc_av_var (param_i) + eloc_av**2 * dpsi_av_var (param_i)  &
                                 + dpsi_av (param_i)**2 * eloc_av_var  + deloc_av_var (param_i)      &
                                 - 2 * eloc_av * dpsi_eloc_av_dpsi_av_covar (param_i, param_i)       &
                                 - 2 * dpsi_av (param_i) * dpsi_eloc_av_eloc_av_covar (param_i)      &
                                 + 2 * dpsi_av (param_i) * eloc_av * dpsi_av_eloc_av_covar (param_i) &
                                 + 2 * deloc_av_dpsi_eloc_av_covar (param_i, param_i)                &
                                 - 2 * eloc_av * deloc_av_dpsi_av_covar (param_i, param_i)           &
                                 - 2 * dpsi_av (param_i) * deloc_av_eloc_av_covar (param_i) 
  enddo

  end subroutine dpsi_eloc_covar_deloc_av_var_bld

! ==============================================================================
  subroutine dpsi_dpsi_c_eloc_av_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_dpsi_c_eloc_av = < ( dpsi - < dpsi > ) * ( dpsi - < dpsi > ) * eloc >
! Description   : This is the Hamitonian for energy mimization with linear method
!
! Created       : J. Toulouse, 27 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i, param_j
! header
  if (header_exe) then

   call object_create ('dpsi_dpsi_c_eloc_av')

   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('dpsi_dpsi_eloc_av')
   call object_needed ('dpsi_eloc_av')
   call object_needed ('dpsi_av')
   call object_needed ('eloc_av')

   return

  endif

! begin
  call object_alloc ('dpsi_dpsi_c_eloc_av', dpsi_dpsi_c_eloc_av, param_nb, param_nb)

  do param_i = 1, param_nb
    do param_j = 1, param_nb
      dpsi_dpsi_c_eloc_av (param_i, param_j) = dpsi_dpsi_eloc_av (param_pairs (param_i, param_j))  &
      - dpsi_av (param_i) * dpsi_eloc_av (param_j) - dpsi_av (param_j) * dpsi_eloc_av (param_i)    &
      + dpsi_av (param_i) * dpsi_av (param_j) * eloc_av
    enddo
  enddo

  end subroutine dpsi_dpsi_c_eloc_av_bld

! ==============================================================================
  subroutine dpsi_deloc_eloc_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_deloc_eloc = dpsi * deloc * eloc
!
! Created       : J. Toulouse, 27 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('dpsi_deloc_eloc')
   call object_average_define ('dpsi_deloc_eloc', 'dpsi_deloc_eloc_av')

   call object_needed ('param_nb')
   call object_needed ('dpsi')
   call object_needed ('deloc')
   call object_needed ('eloc')

   return

  endif

! begin
  call object_alloc ('dpsi_deloc_eloc', dpsi_deloc_eloc, param_nb, param_nb)
  call object_alloc ('dpsi_deloc_eloc_av', dpsi_deloc_eloc_av, param_nb, param_nb)

  do param_i = 1, param_nb
    do param_j = 1, param_nb
      dpsi_deloc_eloc (param_i, param_j) =  dpsi (param_i) * deloc (param_j) * eloc
    enddo
  enddo

  end subroutine dpsi_deloc_eloc_bld

! ==============================================================================
  subroutine dpsi_deloc_eloc_tc_av_bld
! ------------------------------------------------------------------------------
! Description   : the following tricovariance
! Description   : dpsi_deloc_eloc_tc_av =  < (< dpsi - <dpsi>) * (deloc - <deloc>) * (eloc - <eloc>) >
!
! Created       : J. Toulouse, 27 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('dpsi_deloc_eloc_tc_av')

   call object_needed ('param_nb')
   call object_needed ('dpsi_deloc_eloc_av')
   call object_needed ('dpsi_av')
   call object_needed ('deloc_eloc_av')
   call object_needed ('dpsi_eloc_covar')
   call object_needed ('deloc_av')
   call object_needed ('dpsi_deloc_covar')
   call object_needed ('eloc_av')

   return

  endif

! begin
  call object_alloc ('dpsi_deloc_eloc_tc_av', dpsi_deloc_eloc_tc_av, param_nb, param_nb)

  do param_i = 1, param_nb
    do param_j = 1, param_nb
      dpsi_deloc_eloc_tc_av (param_i, param_j) =  dpsi_deloc_eloc_av (param_i, param_j) - dpsi_av (param_i) * deloc_eloc_av (param_j) &
         - dpsi_eloc_covar (param_i) * deloc_av (param_j) - dpsi_deloc_covar (param_i, param_j) * eloc_av
    enddo
  enddo

  end subroutine dpsi_deloc_eloc_tc_av_bld

! ==============================================================================
  subroutine dpsi_dpsi_eloc_eloc_qc_av_bld
! ------------------------------------------------------------------------------
! Description   : the following quadricovariance
! Description   : dpsi_dpsi_eloc_eloc_qc_av = <(<dpsi-<dpsi>)*(<dpsi-<dpsi>)*(eloc-<eloc>)*(eloc-<eloc>)>
!
! Created       : J. Toulouse, 27 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('dpsi_dpsi_eloc_eloc_qc_av')

   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('dpsi_dpsi_eloc_sq_av')
   call object_needed ('dpsi_av')
   call object_needed ('dpsi_eloc_sq_av')
   call object_needed ('eloc_sq_av')
   call object_needed ('dpsi_dpsi_c_eloc_av')
   call object_needed ('eloc_av')
   call object_needed ('dpsi_dpsi_covar')

   return

  endif

! begin
  call object_alloc ('dpsi_dpsi_eloc_eloc_qc_av', dpsi_dpsi_eloc_eloc_qc_av, param_nb, param_nb)

  do param_i = 1, param_nb
    do param_j = 1, param_nb
      dpsi_dpsi_eloc_eloc_qc_av (param_i, param_j) = dpsi_dpsi_eloc_sq_av (param_pairs (param_i, param_j)) &
      - dpsi_av (param_i) * dpsi_eloc_sq_av (param_j) - dpsi_av (param_j) * dpsi_eloc_sq_av (param_i)  &
      + dpsi_av (param_i) * dpsi_av (param_j) * eloc_sq_av                                       &
      - 2.d0 * dpsi_dpsi_c_eloc_av (param_i, param_j) * eloc_av  + dpsi_dpsi_covar (param_i, param_j) * (eloc_av**2)
    enddo
  enddo

  end subroutine dpsi_dpsi_eloc_eloc_qc_av_bld

! ==============================================================================
  subroutine gradient_energy_bld
! ------------------------------------------------------------------------------
! Description   : gradient_energy =  2 [ <dpsi_eloc> - <dpsi_lin> <eloc> ]
!
! Created       : J. Toulouse, 15 Jan 2006
! Modified      : J. Toulouse, 06 Aug 2008, add geometry parameters case
! Modified      : J. Toulouse, 29 Aug 2014, add DMC case
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i

! header
  if (header_exe) then

   call object_create ('gradient_energy')

   call object_needed ('param_nb')
   call object_needed ('dpsi_eloc_covar')
   call object_needed ('is_param_type_geo')

   return

  endif

! begin
  call object_alloc ('gradient_energy', gradient_energy, param_nb)

  if (l_opt_geo .or. l_mode_dmc) then
    call object_provide_by_index (gradient_energy_bld_index, deloc_av_index)
  endif
  
  do param_i = 1, param_nb
   gradient_energy (param_i) = 2.d0 * dpsi_eloc_covar (param_i)

!  add derivative of local energy for geometry parameters
   if (l_opt_geo .and. is_param_type_geo (param_i)) then
     gradient_energy (param_i) = gradient_energy (param_i) + deloc_av (param_i) 
   endif

!  add twice derivative of local energy for approximate DMC gradient
   if (l_mode_dmc) then
     gradient_energy (param_i) = gradient_energy (param_i) + 2*deloc_av (param_i) 
   endif

  enddo ! param_i

  end subroutine gradient_energy_bld

! ==============================================================================
  subroutine gradient_energy_blk_bld
! ------------------------------------------------------------------------------
! Description   : gradient_energy on one block
!
! Created       : J. Toulouse, 08 May 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('gradient_energy_blk')

   call object_needed ('param_nb')
   call object_needed ('dpsi_eloc_blk_covar')

   return

  endif

! begin
  call object_alloc ('gradient_energy_blk', gradient_energy_blk, param_nb)
  gradient_energy_blk =  2.d0 * dpsi_eloc_blk_covar

  end subroutine gradient_energy_blk_bld

! ==============================================================================
  subroutine gradient_variance_bld
! ------------------------------------------------------------------------------
! Description   : gradient =  2 [ <dpsi_eloc> - <dpsi_lin> <eloc> ]
!
! Created       : J. Toulouse, 25 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('gradient_variance', gradient_variance_index)

   call object_needed ('param_nb')
   call object_needed ('deloc_eloc_covar')
   call object_needed ('dpsi_eloc_sq_covar')
   call object_needed ('eloc_av')
   call object_needed ('gradient_energy')

   return

  endif

! begin
  call object_alloc ('gradient_variance', gradient_variance, param_nb)
  gradient_variance =  2.d0 * (deloc_eloc_covar + dpsi_eloc_sq_covar - eloc_av * gradient_energy)

  end subroutine gradient_variance_bld

! ==============================================================================
  subroutine gradient_bld
! ------------------------------------------------------------------------------
! Description   : gradient of linear combinaison of energy and variance
!
! Created       : J. Toulouse, 25 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('gradient')

   call object_needed ('param_nb')
   call object_needed ('gradient_energy')
   call object_needed ('p_var')

   return

  endif

! begin
  call object_alloc ('gradient', gradient, param_nb)

  if (p_var /= 0.d0) then
    call object_provide_by_index (gradient_bld_index, gradient_variance_index)
    gradient =  (1.d0 - p_var) * gradient_energy + p_var * gradient_variance
  else
    gradient =  gradient_energy
  endif

  end subroutine gradient_bld

! ==============================================================================
  subroutine gradient_norm_bld
! ------------------------------------------------------------------------------
! Description   : norm of gradient
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  integer i

! header
  if (header_exe) then

   call object_create ('gradient_norm')

   call object_needed ('param_nb')
   call object_needed ('gradient')

   return

  endif

! begin
  call object_associate ('gradient_norm', gradient_norm)
  call object_associate ('gradient_norm_err', gradient_norm_err)

  gradient_norm = 0.d0
  do i = 1, param_nb
   gradient_norm =  gradient_norm + gradient (i)**2
  enddo
  gradient_norm = dsqrt(gradient_norm)

  end subroutine gradient_norm_bld

! ==============================================================================
  subroutine hessian_variance_lm_bld
! ------------------------------------------------------------------------------
! Description   : Levenberg-Marquardt approximation to Hessian of variance
!
! Created       : J. Toulouse, 25 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('hessian_variance_lm', hessian_variance_lm_index)

   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('deloc_deloc_av')
   call object_needed ('deloc_av')
   call object_needed ('gradient_energy')

   return

  endif

! begin
  call object_alloc ('hessian_variance_lm', hessian_variance_lm, param_nb, param_nb)

  do param_i = 1, param_nb
   do param_j = 1, param_nb
    hessian_variance_lm (param_i, param_j) =  2.d0 * (deloc_deloc_av (param_pairs (param_i, param_j)) - deloc_av (param_i) * gradient_energy (param_j) &
         - deloc_av (param_j) * gradient_energy (param_i) + gradient_energy (param_i) * gradient_energy (param_j) )
   enddo ! param_j
  enddo ! param_i

  end subroutine hessian_variance_lm_bld

! ==============================================================================
  subroutine hessian_variance_lmcov_bld
! ------------------------------------------------------------------------------
! Description   : Levenberg-Marquardt approximation to Hessian of variance
! Description   : expressed in terms of covariances
!
! Created       : J. Toulouse, 25 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('hessian_variance_lmcov', hessian_variance_lmcov_index)

   call object_needed ('param_nb')
   call object_needed ('deloc_deloc_covar')
   call object_needed ('gradient_energy')

   return

  endif

! begin
  call object_alloc ('hessian_variance_lmcov', hessian_variance_lmcov, param_nb, param_nb)

  do param_i = 1, param_nb
   do param_j = 1, param_nb
    hessian_variance_lmcov (param_i, param_j) =  2.d0 * (deloc_deloc_covar (param_i, param_j) + gradient_energy (param_i) * gradient_energy (param_j) )
   enddo ! param_j
  enddo ! param_i

  end subroutine hessian_variance_lmcov_bld

! ==============================================================================
  subroutine hessian_variance_lin_bld
! ------------------------------------------------------------------------------
! Description   : Hessian of variance for normalized linear wave function
! Description   : does not work as well as the hessian_variance_lm
! Description   : probably because more noisy
!
! Created       : J. Toulouse, 27 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i, param_j

! header
  if (header_exe) then

   call object_create ('hessian_variance_lin', hessian_variance_lin_index)

   call object_needed ('param_nb')
   call object_needed ('deloc_deloc_covar')
   call object_needed ('dpsi_deloc_eloc_tc_av')
   call object_needed ('dpsi_dpsi_eloc_eloc_qc_av')
   call object_needed ('dpsi_eloc_covar')
   call object_needed ('eloc_var')
   call object_needed ('dpsi_dpsi_covar')

   return

  endif

! begin
  call object_alloc ('hessian_variance_lin', hessian_variance_lin, param_nb, param_nb)

  do param_i = 1, param_nb
   do param_j = 1, param_nb
    hessian_variance_lin (param_i, param_j) = 2.d0 * (deloc_deloc_covar (param_i, param_j)  &
         + dpsi_deloc_eloc_tc_av (param_i, param_j) + dpsi_deloc_eloc_tc_av (param_j, param_i) &
         + dpsi_dpsi_eloc_eloc_qc_av (param_i, param_j) - 4.d0 * dpsi_eloc_covar (param_i) * dpsi_eloc_covar (param_j) &
         - eloc_var * dpsi_dpsi_covar (param_i, param_j) )
   enddo ! param_j
  enddo ! param_i

  end subroutine hessian_variance_lin_bld

! ==============================================================================
  subroutine hessian_variance_bld
! ------------------------------------------------------------------------------
! Description   : Hessian of variance of local energy
!
! Created       : J. Toulouse, 27 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('hessian_variance', hessian_variance_index)

   return

  endif

! begin
  call object_alloc ('hessian_variance', hessian_variance, param_nb, param_nb)

  select case(trim(hessian_variance_type))
   case ('levenberg_marquardt')
    call object_provide_by_index (hessian_variance_bld_index, hessian_variance_lm_index)
    hessian_variance (:,:) = hessian_variance_lm (:,:)
   case ('levenberg_marquardt_cov')
    call object_provide_by_index (hessian_variance_bld_index, hessian_variance_lmcov_index)
    hessian_variance (:,:) = hessian_variance_lmcov (:,:)
   case ('linear')
    call object_provide_by_index (hessian_variance_bld_index, hessian_variance_lin_index)
    hessian_variance (:,:) = hessian_variance_lin (:,:)
   case default
    call die (here, 'unknown variance hessian type >'+trim(hessian_variance_type)+'<.')
  end select

  end subroutine hessian_variance_bld

end module deriv_mod

