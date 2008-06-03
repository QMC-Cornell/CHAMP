 module dipole_moment_mod

  use all_tools_mod
  use electrons_mod
  use nuclei_mod
  use deriv_mod
  use opt_lin_mod

! Declaration of global variables and default values
  logical                                     :: l_dipole_moment = .false.
  logical                                     :: l_dipole_moment_debye = .true.
  logical                                     :: l_dipole_moment_hf = .true.
  logical                                     :: l_dipole_moment_zv = .false.
  logical                                     :: l_dipole_moment_zvzb = .false.
  logical                                     :: l_dipole_moment_lin = .false.
  real(dp)                                    :: dipole_moment_units = 1.d0
  real(dp), parameter                         :: atomic_unit_to_debye = 1.d0/0.393456d0
  real(dp), allocatable                       :: dipole_moment_origin (:)
  real(dp), allocatable                       :: dipole_moment_nucl (:)
  real(dp), allocatable                       :: dipole_moment (:,:)
  real(dp), allocatable                       :: dipole_moment_av (:)
  real(dp), allocatable                       :: dipole_moment_bav (:)
  real(dp), allocatable                       :: dipole_moment_av_var (:)
  real(dp), allocatable                       :: dipole_moment_av_err (:)
  real(dp), allocatable                       :: dipole_moment_dpsi (:,:)
  real(dp), allocatable                       :: dipole_moment_dpsi_av (:,:)
  real(dp), allocatable                       :: dipole_moment_dpsi_dpsi (:,:)
  real(dp), allocatable                       :: dipole_moment_dpsi_dpsi_av (:,:)
  real(dp), allocatable                       :: dipole_moment_deloc (:,:)
  real(dp), allocatable                       :: dipole_moment_deloc_av (:,:)
  real(dp), allocatable                       :: dipole_moment_deloc_bav (:,:)
  real(dp), allocatable                       :: dipole_moment_deloc_covar (:,:)
  real(dp), allocatable                       :: dipole_moment_deloc_blk_covar (:,:)
  real(dp), allocatable                       :: dipole_moment_av_deloc_av_covar (:,:)
  real(dp), allocatable                       :: dipole_moment_av_dpsi_eloc_av_covar (:,:)
  real(dp), allocatable                       :: dipole_moment_av_dpsi_av_covar (:,:)
  real(dp), allocatable                       :: dipole_moment_av_eloc_av_covar (:)
  real(dp), allocatable                       :: dipole_moment_zv_coef (:,:)
  real(dp), allocatable                       :: dipole_moment_zv_coef_blk (:,:)
  real(dp), allocatable                       :: dipole_moment_zv_av (:)
  real(dp), allocatable                       :: dipole_moment_zv_bav (:)
  real(dp), allocatable                       :: dipole_moment_zv_bav_av (:)
  real(dp), allocatable                       :: dipole_moment_zv_bav_av_err (:)
  real(dp), allocatable                       :: dipole_moment_zv_bav_sq (:)
  real(dp), allocatable                       :: dipole_moment_zv_bav_sq_av (:)
  real(dp), allocatable                       :: dipole_moment_zv_av_var (:)
  real(dp), allocatable                       :: dipole_moment_zv_av_err (:)
  real(dp), allocatable                       :: dipole_moment_zvzb_av (:)
  real(dp), allocatable                       :: dipole_moment_delta_zb_av (:)
  real(dp), allocatable                       :: dipole_moment_delta_zb_av_var (:)
  real(dp), allocatable                       :: dipole_moment_delta_zb_av_err (:)
  real(dp), allocatable                       :: dipole_moment_zvzb_av_var (:)
  real(dp), allocatable                       :: dipole_moment_zvzb_av_err (:)
  real(dp), allocatable                       :: dipole_moment_zvzb_bav (:)
  real(dp), allocatable                       :: dipole_moment_zvzb_bav_av (:)
  real(dp), allocatable                       :: dipole_moment_zvzb_bav_av_err (:)
  real(dp), allocatable                       :: dipole_moment_lin_av (:)

  contains

!===========================================================================
  subroutine dipole_moment_menu
!---------------------------------------------------------------------------
! Description : menu for calculation of dipole moment
!
! Created     : J. Toulouse, 02 Feb 2008
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'dipole_moment_menu'
  integer dipole_moment_origin_dim
  character(len=max_string_len) estimator, units
  logical l_origin_given

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of dipole_moment menu --------------------------------------------------------------------------'

! initialization
  l_dipole_moment = .true.
  l_origin_given = .false.

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))
   case ('help')
    write(6,'(a)') 'HELP for menu dipole_moment'
    write(6,'(a)') 'dipole_moment'
    write(6,'(a)') ' units = {debye,ua}: units for printing dipole moment (default=debye)'
    write(6,'(a)') ' origin real real real end : origin with respect to which dipole moment is calculated (default = center of mass)'
    write(6,'(a)') ' estimator = {hf,zv,zvzb,linear}: estimator to use (default=hf)'
    write(6,'(a)') 'end'

   case ('units')
    call get_next_value (units)
    select case (units)
     case ('debye');   l_dipole_moment_debye = .true.
     case ('ua');      l_dipole_moment_debye = .false.
     case default; call die (lhere, 'unknown keyword >'+trim(word)+'<.')
    end select

   case ('origin')
    call get_next_value_list ('dipole_moment_origin', dipole_moment_origin, dipole_moment_origin_dim)
    call require ('dipole_moment_origin_dim == ndim', dipole_moment_origin_dim == ndim)
    l_origin_given = .true.

   case ('estimator')
    call get_next_value (estimator)
    select case (estimator)
     case ('hf');   l_dipole_moment_hf = .true.
     case ('zv');   l_dipole_moment_zv = .true.
     case ('zvzb'); l_dipole_moment_zv = .true.; l_dipole_moment_zvzb = .true.
     case ('linear'); l_dipole_moment_lin = .true.; l_dipole_moment_lin = .true.
     case default; call die (lhere, 'unknown keyword >'+trim(word)+'<.')
    end select

   case ('end')
    exit

   case default
    call die (lhere, 'unknown keyword >'+trim(word)+'<.')

  end select

  enddo ! end loop over menu lines

  if (l_dipole_moment_debye) then
   dipole_moment_units = atomic_unit_to_debye
  else
   dipole_moment_units = 1.d0
  endif

  if (l_origin_given) then
   write(6,'(a,3f)') ' Dipole moment will be calculated with respect to the origin: ', dipole_moment_origin (:)
  else
   write(6,'(a)') ' Dipole moment will be calculated with respect to the center of mass.'
  endif
  if (l_dipole_moment_hf) then
   write(6,'(a)') ' Dipole moment will be calculated with Hellmann-Feynman estimator.'
  endif
  if (l_dipole_moment_zv) then
   write(6,'(a)') ' Dipole moment will be calculated with zero-variance estimator.'
  endif
  if (l_dipole_moment_zvzb) then
   write(6,'(a)') ' Dipole moment will be calculated with zero-variance zero-bias estimator.'
  endif
  if (l_dipole_moment_lin) then
   write(6,'(a)') ' Dipole moment will be calculated with linear wave function.'
  endif


! request averages and statistical errors
  call object_average_request ('dipole_moment_av')
  call object_error_request ('dipole_moment_av_err')

  if (l_dipole_moment_zv .or. l_dipole_moment_zvzb) then
   call object_average_request ('deloc_av')
   call object_average_request ('deloc_deloc_av')
   call object_average_request ('dipole_moment_deloc_av')
   call object_covariance_request ('dipole_moment_av_deloc_av_covar')
   call object_covariance_request ('deloc_av_deloc_av_covar')
!   call object_average_request ('dipole_moment_zv_bav_sq_av')

  endif

  if (l_dipole_moment_zv) then
   call object_error_request ('dipole_moment_zv_av_err')
   call object_average_request ('dipole_moment_zv_bav_av') !!!!
   call object_error_request ('dipole_moment_zv_bav_av_err') !!!
  endif

  if (l_dipole_moment_zvzb) then
   call object_average_request ('dpsi_av')
   call object_average_request ('dpsi_eloc_av')

   call object_covariance_request ('dipole_moment_av_dpsi_eloc_av_covar')
   call object_covariance_request ('dipole_moment_av_dpsi_av_covar')
   call object_covariance_request ('dipole_moment_av_eloc_av_covar')
   call object_covariance_request ('deloc_av_dpsi_eloc_av_covar')
   call object_covariance_request ('deloc_av_dpsi_av_covar')
   call object_covariance_request ('deloc_av_eloc_av_covar')
   call object_covariance_request ('dpsi_eloc_av_dpsi_eloc_av_covar')
   call object_covariance_request ('dpsi_eloc_av_dpsi_av_covar')
   call object_covariance_request ('dpsi_eloc_av_eloc_av_covar')
   call object_covariance_request ('dpsi_av_dpsi_av_covar')
   call object_covariance_request ('dpsi_av_eloc_av_covar')

   call object_error_request ('dipole_moment_zvzb_av_err')
   call object_error_request ('dipole_moment_delta_zb_av_err') !

   call object_average_request ('dipole_moment_zvzb_bav_av') !!!!
   call object_error_request ('dipole_moment_zvzb_bav_av_err') !!!
  endif

  if (l_dipole_moment_lin) then
   call object_average_request ('dpsi_av')
   call object_average_request ('dpsi_eloc_av')
   call object_average_request ('dpsi_dpsi_av')
   call object_average_request ('deloc_av')
   call object_average_request ('dpsi_deloc_av')
   call object_average_request ('dpsi_dpsi_eloc_av')
   call object_average_request ('dipole_moment_dpsi_av')
   call object_average_request ('dipole_moment_dpsi_dpsi_av')
  endif

  write(6,'(a)') 'End of dipole_moment menu --------------------------------------------------------------------------------'

  end subroutine dipole_moment_menu

! ==============================================================================
  subroutine dipole_moment_origin_bld
! ------------------------------------------------------------------------------
! Description   : default origin for dipole moment
!
! Created       : J. Toulouse, 02 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_origin')

   call object_needed ('ndim')
   call object_needed ('mass_nucl_center')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_origin', dipole_moment_origin, ndim)
  dipole_moment_origin (:) = mass_nucl_center (:)

  end subroutine dipole_moment_origin_bld

! ==============================================================================
  subroutine dipole_moment_nucl_bld
! ------------------------------------------------------------------------------
! Description   : nuclear contribution to dipole moment
!
! Created       : J. Toulouse, 02 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, cent_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_nucl')

   call object_needed ('ndim')
   call object_needed ('znuc')
   call object_needed ('iwctype')
   call object_needed ('ncent')
   call object_needed ('cent')
   call object_needed ('dipole_moment_origin')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_nucl', dipole_moment_nucl, ndim)

  do dim_i = 1, ndim
   dipole_moment_nucl (dim_i) = 0.d0
   do cent_i = 1, ncent
    dipole_moment_nucl (dim_i) = dipole_moment_nucl (dim_i) + znuc(iwctype(cent_i)) * (cent (dim_i, cent_i) - dipole_moment_origin (dim_i))
   enddo ! cent_i
  enddo ! dim_i

  end subroutine dipole_moment_nucl_bld

! ==============================================================================
  subroutine dipole_moment_bld
! ------------------------------------------------------------------------------
! Description   : total dipole moment
!
! Created       : J. Toulouse, 02 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer elec_i, dim_i, walk_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment')
   call object_block_average_define ('dipole_moment', 'dipole_moment_bav')
   call object_average_walk_define ('dipole_moment', 'dipole_moment_av')
   call object_variance_define ('dipole_moment_av', 'dipole_moment_av_var')
   call object_error_define ('dipole_moment_av', 'dipole_moment_av_err')

   call object_needed ('nwalk')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('dipole_moment_nucl')
   call object_needed ('coord_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('dipole_moment', dipole_moment, ndim, nwalk)
  call object_alloc ('dipole_moment_av', dipole_moment_av, ndim)
  call object_alloc ('dipole_moment_bav', dipole_moment_bav, ndim)
  call object_alloc ('dipole_moment_av_var', dipole_moment_av_var, ndim)
  call object_alloc ('dipole_moment_av_err', dipole_moment_av_err, ndim)

  do walk_i = 1, nwalk
   do dim_i = 1, ndim
    dipole_moment (dim_i, walk_i) = dipole_moment_nucl (dim_i)
    do elec_i = 1, nelec
     dipole_moment (dim_i, walk_i) = dipole_moment (dim_i, walk_i) - (coord_elec_wlk (dim_i, elec_i, walk_i) - dipole_moment_origin (dim_i))
    enddo ! elec_i
   enddo ! dim_i
  enddo ! walk_i

  end subroutine dipole_moment_bld

! ==============================================================================
  subroutine dipole_moment_deloc_bld
! ------------------------------------------------------------------------------
! Description   : dipole moment * deloc
!
! Created       : J. Toulouse, 12 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer param_i, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_deloc')
   call object_block_average_define ('dipole_moment_deloc', 'dipole_moment_deloc_bav')
   call object_average_define ('dipole_moment_deloc', 'dipole_moment_deloc_av')
   call object_covariance_define ('dipole_moment_av', 'deloc_av', 'dipole_moment_av_deloc_av_covar')
   call object_covariance_define ('dipole_moment_av', 'dpsi_eloc_av', 'dipole_moment_av_dpsi_eloc_av_covar')
   call object_covariance_define ('dipole_moment_av', 'dpsi_av', 'dipole_moment_av_dpsi_av_covar')
   call object_covariance_define ('dipole_moment_av', 'eloc_av', 'dipole_moment_av_eloc_av_covar')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment')
   call object_needed ('deloc')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_deloc', dipole_moment_deloc, ndim, param_nb)
  call object_alloc ('dipole_moment_deloc_av', dipole_moment_deloc_av, ndim, param_nb)
  call object_alloc ('dipole_moment_deloc_bav', dipole_moment_deloc_bav, ndim, param_nb)
  call object_alloc ('dipole_moment_av_deloc_av_covar', dipole_moment_av_deloc_av_covar, ndim, param_nb)
  call object_alloc ('dipole_moment_av_dpsi_eloc_av_covar', dipole_moment_av_dpsi_eloc_av_covar, ndim, param_nb)
  call object_alloc ('dipole_moment_av_dpsi_av_covar', dipole_moment_av_dpsi_av_covar, ndim, param_nb)
  call object_alloc ('dipole_moment_av_eloc_av_covar', dipole_moment_av_eloc_av_covar, ndim)

  do dim_i = 1, ndim
   do param_i = 1, param_nb
    dipole_moment_deloc (dim_i, param_i) = dipole_moment (dim_i, 1) * deloc (param_i)
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_deloc_bld

! ==============================================================================
  subroutine dipole_moment_dpsi_bld
! ------------------------------------------------------------------------------
! Description   : dipole moment * dpsi
!
! Created       : J. Toulouse, 18 May 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer param_i, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_dpsi')
   call object_average_define ('dipole_moment_dpsi', 'dipole_moment_dpsi_av')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment')
   call object_needed ('dpsi')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_dpsi', dipole_moment_dpsi, ndim, param_nb)
  call object_alloc ('dipole_moment_dpsi_av', dipole_moment_dpsi_av, ndim, param_nb)

  do dim_i = 1, ndim
   do param_i = 1, param_nb
    dipole_moment_dpsi (dim_i, param_i) = dipole_moment (dim_i, 1) * dpsi (param_i)
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_dpsi_bld

! ==============================================================================
  subroutine dipole_moment_dpsi_dpsi_bld
! ------------------------------------------------------------------------------
! Description   : dipole moment * dpsi_dpsi
!
! Created       : J. Toulouse, 18 May 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer pair_i, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_dpsi_dpsi')
   call object_average_define ('dipole_moment_dpsi_dpsi', 'dipole_moment_dpsi_dpsi_av')

   call object_needed ('ndim')
   call object_needed ('param_pairs_nb')
   call object_needed ('dipole_moment')
   call object_needed ('dpsi_dpsi')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_dpsi_dpsi', dipole_moment_dpsi_dpsi, ndim, param_pairs_nb)
  call object_alloc ('dipole_moment_dpsi_dpsi_av', dipole_moment_dpsi_dpsi_av, ndim, param_pairs_nb)

  do dim_i = 1, ndim
   do pair_i = 1, param_pairs_nb
    dipole_moment_dpsi_dpsi (dim_i, pair_i) = dipole_moment (dim_i, 1) * dpsi_dpsi (pair_i)
   enddo ! pair_i
  enddo ! dim_i

  end subroutine dipole_moment_dpsi_dpsi_bld

! ==============================================================================
  subroutine dipole_moment_deloc_covar_bld
! ------------------------------------------------------------------------------
! Description   : covariance : < dipole moment * deloc > - < dipole moment > * < deloc >
!
! Created       : J. Toulouse, 12 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer param_i, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_deloc_covar')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_deloc_av')
   call object_needed ('dipole_moment_av')
   call object_needed ('deloc_av')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_deloc_covar', dipole_moment_deloc_covar, ndim, param_nb)

  do dim_i = 1, ndim
   do param_i = 1, param_nb
    dipole_moment_deloc_covar (dim_i, param_i) = dipole_moment_deloc_av (dim_i, param_i) - dipole_moment_av (dim_i) * deloc_av (param_i)
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_deloc_covar_bld

! ==============================================================================
  subroutine dipole_moment_deloc_blk_covar_bld
! ------------------------------------------------------------------------------
! Description   : dipole_moment_deloc_covar calculated on current block
!
! Created       : J. Toulouse, 20 Apr 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer param_i, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_deloc_blk_covar')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_deloc_bav')
   call object_needed ('dipole_moment_bav')
   call object_needed ('deloc_bav')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_deloc_blk_covar', dipole_moment_deloc_blk_covar, ndim, param_nb)

  do dim_i = 1, ndim
   do param_i = 1, param_nb
    dipole_moment_deloc_blk_covar (dim_i, param_i) = dipole_moment_deloc_bav (dim_i, param_i) - dipole_moment_bav (dim_i) * deloc_bav (param_i)
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_deloc_blk_covar_bld

! ==============================================================================
  subroutine dipole_moment_zv_coef_bld
! ------------------------------------------------------------------------------
! Description   : coefficient minimizing the variance of zero-variance estimator of total dipole moment
!
! Created       : J. Toulouse, 12 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_zv_coef')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_deloc_covar')
   call object_needed ('deloc_deloc_covar_inv')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zv_coef', dipole_moment_zv_coef, ndim, param_nb)

  do dim_i = 1, ndim
   do param_i = 1, param_nb
    dipole_moment_zv_coef (dim_i, param_i) = 0.d0
    do param_j = 1, param_nb
     dipole_moment_zv_coef (dim_i, param_i) = dipole_moment_zv_coef (dim_i, param_i) - deloc_deloc_covar_inv (param_i, param_j) * dipole_moment_deloc_covar (dim_i, param_j)
    enddo ! param_j
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_zv_coef_bld

! ==============================================================================
  subroutine dipole_moment_zv_coef_blk_bld
! ------------------------------------------------------------------------------
! Description   : coefficient minimizing the variance of zero-variance estimator of total dipole moment
! Description   : calculated on current block average
!
! Created       : J. Toulouse, 20 Apr 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_zv_coef_blk')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_deloc_blk_covar')
   call object_needed ('deloc_deloc_blk_covar_inv')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zv_coef_blk', dipole_moment_zv_coef_blk, ndim, param_nb)

  do dim_i = 1, ndim
   do param_i = 1, param_nb
    dipole_moment_zv_coef_blk (dim_i, param_i) = 0.d0
    do param_j = 1, param_nb
     dipole_moment_zv_coef_blk (dim_i, param_i) = dipole_moment_zv_coef_blk (dim_i, param_i) - deloc_deloc_blk_covar_inv (param_i, param_j) * dipole_moment_deloc_blk_covar (dim_i, param_j)
    enddo ! param_j
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_zv_coef_blk_bld

! ==============================================================================
  subroutine dipole_moment_zv_av_bld
! ------------------------------------------------------------------------------
! Description   : average of zero-variance estimator total dipole moment
!
! Created       : J. Toulouse, 12 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_zv_av')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_zv_coef')
   call object_needed ('dipole_moment_av')
   call object_needed ('deloc_av')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zv_av', dipole_moment_zv_av, ndim)

  do dim_i = 1, ndim
   dipole_moment_zv_av (dim_i) = dipole_moment_av (dim_i)
   do param_i = 1, param_nb
    dipole_moment_zv_av (dim_i) = dipole_moment_zv_av (dim_i) + dipole_moment_zv_coef (dim_i, param_i) * deloc_av (param_i)
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_zv_av_bld

! ==============================================================================
  subroutine dipole_moment_zv_bav_bld
! ------------------------------------------------------------------------------
! Description   : block average of zero-variance estimator total dipole moment
!
! Created       : J. Toulouse, 20 Apr 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_zv_bav')
   call object_average_define_from_block_average ('dipole_moment_zv_bav',  'dipole_moment_zv_bav_av')
   call object_error_define ('dipole_moment_zv_bav_av',  'dipole_moment_zv_bav_av_err')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_zv_coef_blk')
   call object_needed ('dipole_moment_bav')
   call object_needed ('deloc_bav')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zv_bav', dipole_moment_zv_bav, ndim)
  call object_alloc ('dipole_moment_zv_bav_av', dipole_moment_zv_bav_av, ndim)
  call object_alloc ('dipole_moment_zv_bav_av_err', dipole_moment_zv_bav_av_err, ndim)

  do dim_i = 1, ndim
   dipole_moment_zv_bav (dim_i) = dipole_moment_bav (dim_i)
   do param_i = 1, param_nb
    dipole_moment_zv_bav (dim_i) = dipole_moment_zv_bav (dim_i) + dipole_moment_zv_coef_blk (dim_i, param_i) * deloc_bav (param_i)
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_zv_bav_bld

! ==============================================================================
  subroutine dipole_moment_zv_bav_sq_bld
! ------------------------------------------------------------------------------
! Description   : square of dipole_moment_zv_bav
!
! Created       : J. Toulouse, 20 Apr 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! begin

! header
  if (header_exe) then

   call object_average_define_from_block_average ('dipole_moment_zv_bav_sq',  'dipole_moment_zv_bav_sq_av')
   call object_create ('dipole_moment_zv_bav_sq')

   call object_needed ('ndim')
   call object_needed ('dipole_moment_zv_bav')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zv_bav_sq', dipole_moment_zv_bav_sq, ndim)
  call object_alloc ('dipole_moment_zv_bav_sq_av', dipole_moment_zv_bav_sq_av, ndim)

  dipole_moment_zv_bav_sq (:) = dipole_moment_zv_bav (:)**2

  end subroutine dipole_moment_zv_bav_sq_bld

! ==============================================================================
  subroutine dipole_moment_zv_av_var_2_bld
! ------------------------------------------------------------------------------
! Description   : variance of average of zero-variance estimator total dipole moment
! Description   : to be checked
!
! Created       : J. Toulouse, 20 Apr 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_zv_av_var')
   call object_error_define_from_variance ('dipole_moment_zv_av_var', 'dipole_moment_zv_av_err')

   call object_needed ('ndim')
   call object_needed ('dipole_moment_zv_bav_sq_av')
   call object_needed ('dipole_moment_zv_av')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zv_av_var', dipole_moment_zv_av_var, ndim)
  call object_alloc ('dipole_moment_zv_av_err', dipole_moment_zv_av_err, ndim)

  dipole_moment_zv_av_var (:) = dipole_moment_zv_bav_sq_av (:)/block_iterations_nb - dipole_moment_zv_av (:)**2

  end subroutine dipole_moment_zv_av_var_2_bld

! ==============================================================================
  subroutine dipole_moment_zv_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of average of zero-variance estimator total dipole moment
!
! Created       : J. Toulouse, 14 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_zv_av_var')
   call object_error_define_from_variance ('dipole_moment_zv_av_var', 'dipole_moment_zv_av_err')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_av_var')
   call object_needed ('dipole_moment_zv_coef')
   call object_needed ('dipole_moment_av_deloc_av_covar')
   call object_needed ('deloc_av_deloc_av_covar')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zv_av_var', dipole_moment_zv_av_var, ndim)
  call object_alloc ('dipole_moment_zv_av_err', dipole_moment_zv_av_err, ndim)

  do dim_i = 1, ndim
   dipole_moment_zv_av_var (dim_i) = dipole_moment_av_var (dim_i)
   do param_i = 1, param_nb
    dipole_moment_zv_av_var (dim_i) = dipole_moment_zv_av_var (dim_i) + 2.d0 * dipole_moment_zv_coef (dim_i, param_i) * dipole_moment_av_deloc_av_covar (dim_i, param_i)
     do param_j = 1, param_nb
      dipole_moment_zv_av_var (dim_i) = dipole_moment_zv_av_var (dim_i) + dipole_moment_zv_coef (dim_i, param_i) * dipole_moment_zv_coef (dim_i, param_j) * deloc_av_deloc_av_covar (param_i, param_j)
     enddo ! param_j
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_zv_av_var_bld

! ==============================================================================
  subroutine dipole_moment_delta_zb_av_bld
! ------------------------------------------------------------------------------
! Description   : average of zero-bias term for total dipole moment
!
! Created       : J. Toulouse, 16 May 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_delta_zb_av')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_zv_coef')
   call object_needed ('gradient_energy')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_delta_zb_av', dipole_moment_delta_zb_av, ndim)

  do dim_i = 1, ndim
   dipole_moment_delta_zb_av (dim_i) = 0.d0
   do param_i = 1, param_nb
    dipole_moment_delta_zb_av (dim_i) = dipole_moment_delta_zb_av (dim_i) + dipole_moment_zv_coef (dim_i, param_i) * gradient_energy (param_i)
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_delta_zb_av_bld

! ==============================================================================
  subroutine dipole_moment_delta_zb_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of average of zero-bias term for total dipole moment
!
! Created       : J. Toulouse, 18 May 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_delta_zb_av_var')
   call object_error_define_from_variance ('dipole_moment_delta_zb_av_var', 'dipole_moment_delta_zb_av_err')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_zv_coef')
   call object_needed ('dpsi_eloc_av_dpsi_eloc_av_covar')
   call object_needed ('eloc_av')
   call object_needed ('dpsi_eloc_av_dpsi_av_covar')
   call object_needed ('dpsi_av')
   call object_needed ('dpsi_eloc_av_eloc_av_covar')
   call object_needed ('dpsi_av_eloc_av_covar')
   call object_needed ('eloc_av_var')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_delta_zb_av_var', dipole_moment_delta_zb_av_var, ndim)
  call object_alloc ('dipole_moment_delta_zb_av_err', dipole_moment_delta_zb_av_err, ndim)

  do dim_i = 1, ndim
   dipole_moment_delta_zb_av_var (dim_i) = 0.d0
   do param_i = 1, param_nb
    do param_j = 1, param_nb
     dipole_moment_delta_zb_av_var (dim_i) = dipole_moment_delta_zb_av_var (dim_i) + 4.d0 * dipole_moment_zv_coef (dim_i, param_i) * dipole_moment_zv_coef (dim_i, param_j) *  &
     (dpsi_eloc_av_dpsi_eloc_av_covar (param_i, param_j) - 2.d0 * eloc_av * dpsi_eloc_av_dpsi_av_covar (param_i, param_j)         &
      - 2.d0 * dpsi_av (param_j) * dpsi_eloc_av_eloc_av_covar (param_i) + (eloc_av**2) * dpsi_av_dpsi_av_covar (param_i, param_j) &
      + 2.d0 * eloc_av * dpsi_av (param_j) * dpsi_av_eloc_av_covar (param_i) + dpsi_av (param_i) * dpsi_av (param_j) * eloc_av_var)

    enddo ! param_j
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_delta_zb_av_var_bld

! ==============================================================================
  subroutine dipole_moment_zvzb_av_bld
! ------------------------------------------------------------------------------
! Description   : average of zero-variance zero-bias estimator total dipole moment
!
! Created       : J. Toulouse, 12 Feb 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_zvzb_av')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_zv_av')
   call object_needed ('dipole_moment_zv_coef')
   call object_needed ('gradient_energy')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zvzb_av', dipole_moment_zvzb_av, ndim)

  do dim_i = 1, ndim
   dipole_moment_zvzb_av (dim_i) = dipole_moment_zv_av (dim_i)
   do param_i = 1, param_nb
    dipole_moment_zvzb_av (dim_i) = dipole_moment_zvzb_av (dim_i) + dipole_moment_zv_coef (dim_i, param_i) * gradient_energy (param_i)
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_zvzb_av_bld

! ==============================================================================
  subroutine dipole_moment_zvzb_bav_bld
! ------------------------------------------------------------------------------
! Description   : block average of zero-variance zero-bias estimator total dipole moment
!
! Created       : J. Toulouse, 08 May 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_zvzb_bav')
   call object_average_define_from_block_average ('dipole_moment_zvzb_bav',  'dipole_moment_zvzb_bav_av')
   call object_error_define ('dipole_moment_zvzb_bav_av',  'dipole_moment_zvzb_bav_av_err')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_zv_bav')
   call object_needed ('dipole_moment_zv_coef_blk')
   call object_needed ('gradient_energy_blk')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zvzb_bav', dipole_moment_zvzb_bav, ndim)
  call object_alloc ('dipole_moment_zvzb_bav_av', dipole_moment_zvzb_bav_av, ndim)
  call object_alloc ('dipole_moment_zvzb_bav_av_err', dipole_moment_zvzb_bav_av_err, ndim)

  do dim_i = 1, ndim
   dipole_moment_zvzb_bav (dim_i) = dipole_moment_zv_bav (dim_i)
   do param_i = 1, param_nb
    dipole_moment_zvzb_bav (dim_i) = dipole_moment_zvzb_bav (dim_i) + dipole_moment_zv_coef_blk (dim_i, param_i) * gradient_energy_blk (param_i)
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_zvzb_bav_bld

! ==============================================================================
  subroutine dipole_moment_zvzb_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of average of zero-variance zero-bias estimator total dipole moment
!
! Created       : J. Toulouse, 22 Apr 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_zvzb_av_var')
   call object_error_define_from_variance ('dipole_moment_zvzb_av_var', 'dipole_moment_zvzb_av_err')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_zv_av_var')
   call object_needed ('dipole_moment_zv_coef')
   call object_needed ('dipole_moment_av_dpsi_eloc_av_covar')
   call object_needed ('eloc_av')
   call object_needed ('dipole_moment_av_dpsi_av_covar')
   call object_needed ('dipole_moment_av_eloc_av_covar')
   call object_needed ('deloc_av_dpsi_eloc_av_covar')
   call object_needed ('deloc_av_dpsi_av_covar')
   call object_needed ('dpsi_av')
   call object_needed ('deloc_av_eloc_av_covar')
   call object_needed ('dpsi_eloc_av_dpsi_eloc_av_covar')
   call object_needed ('dpsi_eloc_av_dpsi_av_covar')
   call object_needed ('dpsi_eloc_av_eloc_av_covar')
   call object_needed ('dpsi_av_dpsi_av_covar')
   call object_needed ('dpsi_av_eloc_av_covar')
   call object_needed ('eloc_av_var')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_zvzb_av_var', dipole_moment_zvzb_av_var, ndim)
  call object_alloc ('dipole_moment_zvzb_av_err', dipole_moment_zvzb_av_err, ndim)

  do dim_i = 1, ndim
   dipole_moment_zvzb_av_var (dim_i) = dipole_moment_zv_av_var (dim_i)
   do param_i = 1, param_nb
    dipole_moment_zvzb_av_var (dim_i) = dipole_moment_zvzb_av_var (dim_i) + 4.d0 * dipole_moment_zv_coef (dim_i, param_i) *  &
         (dipole_moment_av_dpsi_eloc_av_covar (dim_i, param_i) - eloc_av * dipole_moment_av_dpsi_av_covar (dim_i, param_i)   &
          - dpsi_av (param_i) * dipole_moment_av_eloc_av_covar (dim_i))

     do param_j = 1, param_nb
      dipole_moment_zvzb_av_var (dim_i) = dipole_moment_zvzb_av_var (dim_i) + 4.d0 * dipole_moment_zv_coef (dim_i, param_i) * dipole_moment_zv_coef (dim_i, param_j) *  &
          (deloc_av_dpsi_eloc_av_covar (param_i, param_j) - eloc_av * deloc_av_dpsi_av_covar (param_i, param_j) - dpsi_av (param_j) * deloc_av_eloc_av_covar (param_i)   &
           + dpsi_eloc_av_dpsi_eloc_av_covar (param_i, param_j) - 2.d0 * eloc_av * dpsi_eloc_av_dpsi_av_covar (param_i, param_j) - 2.d0 * dpsi_av (param_j) * dpsi_eloc_av_eloc_av_covar (param_i)  &
           + (eloc_av**2) * dpsi_av_dpsi_av_covar (param_i, param_j) + 2.d0 * eloc_av * dpsi_av (param_j) * dpsi_av_eloc_av_covar (param_i) + dpsi_av (param_i) * dpsi_av (param_j) * eloc_av_var )

     enddo ! param_j
   enddo ! param_i
  enddo ! dim_i

  end subroutine dipole_moment_zvzb_av_var_bld

! ==============================================================================
  subroutine dipole_moment_lin_av_bld
! ------------------------------------------------------------------------------
! Description   : average of total dipole moment over linear wave function
!
! Created       : J. Toulouse, 18 May 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, param_i, param_j, pair

! begin

! header
  if (header_exe) then

   call object_create ('dipole_moment_lin_av')

   call object_needed ('ndim')
   call object_needed ('param_nb')
   call object_needed ('dipole_moment_av')
   call object_needed ('delta_lin')
   call object_needed ('dipole_moment_dpsi_av')
   call object_needed ('dipole_moment_dpsi_dpsi_av')
   call object_needed ('param_pairs')
   call object_needed ('psi_lin_norm_sq')

   return

  endif

! allocations
  call object_alloc ('dipole_moment_lin_av', dipole_moment_lin_av, ndim)

  do dim_i = 1, ndim
   dipole_moment_lin_av (dim_i) = dipole_moment_av (dim_i) 
   do param_i = 1, param_nb
    dipole_moment_lin_av (dim_i) = dipole_moment_lin_av (dim_i) + 2.d0 * delta_lin (param_i) * dipole_moment_dpsi_av (dim_i, param_i)
    do param_j = 1, param_nb
      pair = param_pairs (param_i, param_j)
      dipole_moment_lin_av (dim_i) = dipole_moment_lin_av (dim_i) + delta_lin (param_i) * delta_lin (param_j) * dipole_moment_dpsi_dpsi_av (dim_i, pair)
    enddo ! param_j
   enddo ! param_i
   dipole_moment_lin_av (dim_i) = dipole_moment_lin_av (dim_i) / psi_lin_norm_sq
  enddo ! dim_i


  end subroutine dipole_moment_lin_av_bld

!===========================================================================
  subroutine print_dipole_moment
!---------------------------------------------------------------------------
! Description : print dipole moment
!
! Created     : J. Toulouse, 12 Feb 2008
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'print_dipole_moment'
  integer dim_i

! begin
  write(6,*)
  if (l_dipole_moment_debye) then
   write(6,'(a)') 'Dipole moment results in Debye:'
  else
   write(6,'(a)') 'Dipole moment results in atomic units:'
  endif

  call object_provide ('ndim')
  call object_provide ('dipole_moment_nucl')
  write(6,'(a)') 'Nuclear dipole moment:'
  do dim_i = 1, ndim
   write(6,'(a,i1,a,f)') 'component # ',dim_i,' : ', dipole_moment_units * dipole_moment_nucl (dim_i)
  enddo

  if (l_dipole_moment_hf) then
   call object_provide ('dipole_moment_av')
   call object_provide ('dipole_moment_av_err')
   write(6,'(a)') 'Total dipole moment using Hellmann-Feynman estimator:'
   do dim_i = 1, ndim
    write(6,'(a,i1,a,f,a,f)') 'component # ',dim_i,' : ', dipole_moment_units * dipole_moment_av (dim_i), ' +-', dipole_moment_units * dipole_moment_av_err (dim_i)
   enddo
  endif

  if (l_dipole_moment_zv) then
   call object_provide ('dipole_moment_zv_av')
   call object_provide ('dipole_moment_zv_av_err')
   call object_provide ('dipole_moment_zv_bav_av')
   call object_provide ('dipole_moment_zv_bav_av_err')
   write(6,'(a)') 'Total dipole moment using zero-variance estimator (average calculated with global averages and error with covariances):'
   do dim_i = 1, ndim
    write(6,'(a,i1,a,f,a,f)') 'component # ',dim_i,' : ', dipole_moment_units * dipole_moment_zv_av (dim_i), ' +-', dipole_moment_units * dipole_moment_zv_av_err (dim_i)
   enddo
   write(6,'(a)') 'Total dipole moment using zero-variance estimator (average and error calculated with block averages):'
   do dim_i = 1, ndim
    write(6,'(a,i1,a,f,a,f)') 'component # ',dim_i,' : ', dipole_moment_units * dipole_moment_zv_bav_av (dim_i), ' +-', dipole_moment_units * dipole_moment_zv_bav_av_err (dim_i)
   enddo
  endif

  if (l_dipole_moment_zvzb) then
   call object_provide ('dipole_moment_zvzb_av')
   call object_provide ('dipole_moment_zvzb_av_err')
   call object_provide ('dipole_moment_zvzb_bav_av')
   call object_provide ('dipole_moment_zvzb_bav_av_err')
   call object_provide ('dipole_moment_delta_zb_av')
   call object_provide ('dipole_moment_delta_zb_av_err')
   write(6,'(a)') 'Total dipole moment using zero-variance zero-bias estimator (average calculated with global averages and error with covariances):'
   do dim_i = 1, ndim
    write(6,'(a,i1,a,f,a,f)') 'component # ',dim_i,' : ', dipole_moment_units * dipole_moment_zvzb_av (dim_i), ' +-', dipole_moment_units * dipole_moment_zvzb_av_err (dim_i)
   enddo
   write(6,'(a)') 'Total dipole moment using zero-variance zero-bias estimator (average and error calculated with block averages):'
   do dim_i = 1, ndim
    write(6,'(a,i1,a,f,a,f)') 'component # ',dim_i,' : ', dipole_moment_units * dipole_moment_zvzb_bav_av (dim_i), ' +-', dipole_moment_units * dipole_moment_zvzb_bav_av_err (dim_i)
   enddo
   write(6,'(a)') 'Zero-bias contribution (average calculated with global averages and error with covariances):'
   do dim_i = 1, ndim
    write(6,'(a,i1,a,f,a,f)') 'component # ',dim_i,' : ', dipole_moment_units * dipole_moment_delta_zb_av (dim_i), ' +-', dipole_moment_units * dipole_moment_delta_zb_av_err (dim_i)
   enddo
  endif

  if (l_dipole_moment_lin) then
   call object_provide ('dipole_moment_lin_av')
   write(6,'(a)') 'Total dipole moment using linear wave function:'
   do dim_i = 1, ndim
    write(6,'(a,i1,a,f,a,f)') 'component # ',dim_i,' : ', dipole_moment_units * dipole_moment_lin_av (dim_i)
   enddo
  endif

  end subroutine print_dipole_moment

end module dipole_moment_mod
