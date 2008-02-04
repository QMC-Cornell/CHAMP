 module forces_mod

  use all_tools_mod
  use grid_mod
  use electrons_mod
  use nuclei_mod
  use psi_mod
  use montecarlo_mod

! Declaration of global variables and default values
  integer                                     :: forces_nb = 0
  integer, allocatable                        :: forces_cent (:)
  integer, allocatable                        :: forces_direct (:)
  real(dp), allocatable                       :: forces_nn (:)
  real(dp), allocatable                       :: forces_bare (:)
  real(dp), allocatable                       :: forces_bare_av (:)
  real(dp), allocatable                       :: forces_bare_av_err (:)
  real(dp), allocatable                       :: forces_zv (:)
  real(dp), allocatable                       :: forces_zv_av (:)
  real(dp), allocatable                       :: forces_zv_av_var (:)
  real(dp), allocatable                       :: forces_zv_av_err (:)
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
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'forces_menu'
# if defined (PATHSCALE)
   character(len=max_string_len) :: forces_list (max_string_array_len) ! for pathscale compiler
# else
   character(len=max_string_len), allocatable  :: forces_list (:)
# endif
  integer force_i

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))
   case ('help')
    write(6,'(a)') 'HELP for menu forces'
    write(6,'(a)') ': forces'
    write(6,'(a)') ':  components 1x 1y 1z 2x 2y 2z end'
    write(6,'(a)') ':  eloc_av_fixed  = [real] fixed value of average of local energy to use in ZB term'
    write(6,'(a)') ':  forces_q_av_fixed  list of reals end: fixed value of average of Q to use in ZB term'
    write(6,'(a)') ': end'

   case ('components')
# if defined (PATHSCALE)
    call get_next_value_list_string ('forces_list', forces_list, forces_nb) ! for pathscale compiler
# else
    call get_next_value_list ('forces_list', forces_list, forces_nb)
# endif
    call object_modified ('forces_nb')

   case ('eloc_av_fixed')
    call get_next_value (eloc_av_fixed)
    call object_modified ('eloc_av_fixed')
    l_eloc_av_fixed = .true.

   case ('forces_q_av_fixed')
    call get_next_value_list ('forces_q_av_fixed', forces_q_av_fixed, forces_nb)
    call object_modified ('forces_q_av_fixed')
    l_forces_q_av_fixed = .true.

   case ('end')
    exit

   case default
    call die (lhere, 'unknown keyword >'+trim(word)+'<.')

  end select

  enddo ! end loop over menu lines

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

  call object_average_request ('forces_bare_av')
  call object_error_request ('forces_bare_av_err')

  call object_average_request ('forces_zv_av')
  call object_error_request ('forces_zv_av_err')

  call object_average_request ('forces_q_av')
  call object_average_request ('forces_q_eloc_av')
  call object_variance_request ('forces_q_eloc_av_var')

  call object_error_request ('forces_zvzb_av_err')

  call object_covariance_request ('forces_zv_av_eloc_av_covar')
  call object_covariance_request ('forces_zv_av_q_eloc_av_covar')
  call object_covariance_request ('forces_zv_av_q_av_covar')
  call object_covariance_request ('forces_q_eloc_av_eloc_av_covar')
  call object_covariance_request ('forces_q_eloc_av_q_av_covar')
  call object_covariance_request ('forces_q_av_eloc_av_covar')

  end subroutine forces_menu

! ==============================================================================
  subroutine forces_nn_bld
! ------------------------------------------------------------------------------
! Description   : nuclei-nuclei contribution to forces
!
! Created       : J. Toulouse, 26 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

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
! Description   : Bare (Hellmann-Feynman) forces
!
! Created       : J. Toulouse, 26 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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
  subroutine forces_q_bld
! ------------------------------------------------------------------------------
! Description   : Simple Q term for forces
! Description   : Assaraf and Caffarel, J. Chem. Phys. (2003)
!
! Created       : J. Toulouse, 09 Aug 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer force_i, elec_i, cent_i, dim_i
  real(dp) eloc_c

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
  implicit none
  include 'commons.h'

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
                          + 4.d0 * (eloc_av**2) * forces_q_av_var (:) + 2.d0 * (forces_q_av (:)**2) * eloc_av_var  &
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
  implicit none
  include 'commons.h'

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
  implicit none
  include 'commons.h'

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

end module forces_mod
