module montecarlo_mod

  use all_tools_mod

! Declaration of global variables and default values

  integer             :: block_nb

! local energy
  real(dp)    :: eloc
  real(dp)    :: eloc_bav
  real(dp)    :: eloc_av
  real(dp)    :: eloc_av_var
  real(dp)    :: eloc_av_err

! local energy over walkers
  real(dp), allocatable   :: eloc_wlk (:)

! local energy for tests
  real(dp), allocatable   :: eloc_test (:)          ! for tests
  real(dp), allocatable   :: eloc_test_bav (:)         ! for tests
  real(dp), allocatable   :: eloc_test_av (:)         ! for tests
  real(dp), allocatable   :: eloc_test_av_err (:)     ! for tests
  real(dp), allocatable   :: eloc_wlk_test (:,:)          ! for tests
  real(dp), allocatable   :: eloc_wlk_test_av (:)         ! for tests
  real(dp), allocatable   :: eloc_wlk_test_av_err (:)     ! for tests
  real(dp), allocatable   :: eloc_wlk_test2 (:,:,:)          ! for tests
  real(dp), allocatable   :: eloc_wlk_test2_av (:,:)         ! for tests
  real(dp), allocatable   :: eloc_wlk_test2_av_err (:,:)     ! for tests

! error on average of local energy
  real(dp)    :: eerr
  real(dp)    :: egerr

! for sigma
  real(dp)    :: eloc_sq
  real(dp)    :: eloc_sq_av
  real(dp)    :: eloc_var

! standard deviation of local energy
  real(dp)    :: sigma

! error on sigma
  real(dp)    :: error_sigma

! autocorrelation time on energy
  real(dp)    :: eloc_tc (MFORCE)

!  real(dp)                  :: e2_eloc_av

  contains

! ==============================================================================
  subroutine eloc_wlk_bld
! ------------------------------------------------------------------------------
! Description   : local energy over walkers
!
! Created       : J. Toulouse, 15 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('eloc_wlk')

   call object_needed ('nwalk')

   return

  endif

! begin
  call object_alloc ('eloc_wlk', eloc_wlk, nwalk)

! allocations
  if (index(mode, 'vmc') /= 0) then
   call object_provide_by_index (eloc_wlk_bld_index, eold_index)
   eloc_wlk (1) = eold (1)

  elseif (index(mode, 'dmc') /= 0) then
   call object_provide_by_index (eloc_wlk_bld_index, eoldw_index)
   eloc_wlk (:) = eoldw (1:nwalk, 1)

  else
   call die (here, 'mode='+trim(mode)+' should contain either vmc or dmc.')

  endif

  end subroutine eloc_wlk_bld

! ==============================================================================
  subroutine eloc_test_bld
! ------------------------------------------------------------------------------
! Description   : local energy for tests
!
! Created       : J. Toulouse, 08 May 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('eloc_test')
   call object_block_average_define ('eloc_test', 'eloc_test_bav')
   call object_average_define_from_block_average ('eloc_test_bav', 'eloc_test_av')
!   call object_average_define ('eloc_test', 'eloc_test_av')
   call object_error_define ('eloc_test_av', 'eloc_test_av_err')

   call object_needed ('eold')

   return

  endif

! begin

! allocations
  call object_alloc ('eloc_test', eloc_test, 3)
  call object_alloc ('eloc_test_bav', eloc_test_bav, 3)
  call object_alloc ('eloc_test_av', eloc_test_av, 3)
  call object_alloc ('eloc_test_av_err', eloc_test_av_err, 3)

  eloc_test (1) = eold (1)

  eloc_test (2) = 2.d0 * eloc_test (1)
  eloc_test (3) = - eloc_test (1)

  end subroutine eloc_test_bld

! ==============================================================================
  subroutine eloc_wlk_test_bld
! ------------------------------------------------------------------------------
! Description   : local energy over walkers for tests
!
! Created       : J. Toulouse, 15 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer walk_i

! header
  if (header_exe) then

   call object_create ('eloc_wlk_test')
   call object_average_walk_define ('eloc_wlk_test', 'eloc_wlk_test_av')
   call object_error_define ('eloc_wlk_test_av', 'eloc_wlk_test_av_err')

   call object_needed ('nwalk')

   return

  endif

! begin

! allocations
  call object_alloc ('eloc_wlk_test', eloc_wlk_test, 3, nwalk)
  call object_alloc ('eloc_wlk_test_av', eloc_wlk_test_av, 3)
  call object_alloc ('eloc_wlk_test_av_err', eloc_wlk_test_av_err, 3)

  if (index(mode, 'vmc') /= 0) then
   call object_provide_by_index (eloc_wlk_test_bld_index, eold_index)
   eloc_wlk_test (1,1) = eold (1)

  elseif (index(mode, 'dmc') /= 0) then
   call object_provide_by_index (eloc_wlk_test_bld_index, eoldw_index)
   do walk_i = 1, nwalk
    eloc_wlk_test (1,walk_i) = eoldw (walk_i, 1)
   enddo

  else
   call die (here, 'mode='+trim(mode)+' should contain either vmc or dmc.')
  endif

   eloc_wlk_test (2,1) = 2.d0 * eloc_wlk_test (1,1)
   eloc_wlk_test (3,1) = - eloc_wlk_test (1,1)


  end subroutine eloc_wlk_test_bld

! ==============================================================================
  subroutine eloc_wlk_test2_bld
! ------------------------------------------------------------------------------
! Description   : local energy over walkers for tests
!
! Created       : J. Toulouse, 15 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer walk_i

! header
  if (header_exe) then

   call object_create ('eloc_wlk_test2')
   call object_average_walk_define ('eloc_wlk_test2', 'eloc_wlk_test2_av')
   call object_error_define ('eloc_wlk_test2_av', 'eloc_wlk_test2_av_err')

   call object_needed ('nwalk')

   return

  endif

! begin

! allocations
  call object_alloc ('eloc_wlk_test2', eloc_wlk_test2, 2, 3, nwalk)
  call object_alloc ('eloc_wlk_test2_av', eloc_wlk_test2_av, 2, 3)
  call object_alloc ('eloc_wlk_test2_av_err', eloc_wlk_test2_av_err, 2, 3)

  if (index(mode, 'vmc') /= 0) then
   call object_provide_by_index (eloc_wlk_test2_bld_index, eold_index)
   eloc_wlk_test2 (1,1,1) = eold (1)

  elseif (index(mode, 'dmc') /= 0) then
   call object_provide_by_index (eloc_wlk_test2_bld_index, eoldw_index)
   do walk_i = 1, nwalk
    eloc_wlk_test2 (1,1,walk_i) = eoldw (walk_i, 1)
   enddo


  else
   write(6,'(4a)') trim(here),': mode=',trim(mode),' should contain either vmc or dmc.'
   call die (here)

  endif

   eloc_wlk_test2 (1,2,:) = -eloc_wlk_test2 (1,1,:)
   eloc_wlk_test2 (1,3,:) = 2.d0*eloc_wlk_test2 (1,1,:)
   eloc_wlk_test2 (2,1,:) = -2.d0*eloc_wlk_test2 (1,1,:)
   eloc_wlk_test2 (2,2,:) = 1.d0
   eloc_wlk_test2 (2,3,:) = -6.d0

  end subroutine eloc_wlk_test2_bld

! ==============================================================================
  subroutine eloc_sq_bld
! ------------------------------------------------------------------------------
! Description   : energy^2
!
! Created       : J. Toulouse, 02 Nov 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('eloc_sq')

   call object_needed ('eloc')

   return

  endif

! begin

! allocations
  call object_associate ('eloc_sq', eloc_sq)
  call object_associate ('eloc_sq_av', eloc_sq_av)

  eloc_sq = eloc**2

  end subroutine eloc_sq_bld

! ==============================================================================
!  subroutine eloc_av_bld
! ------------------------------------------------------------------------------
! Description   : Average of local energy
!
! Created       : J. Toulouse, 24 Sep 2007
! ------------------------------------------------------------------------------
!  implicit none
!  include 'commons.h'
!
!! header
!  if (header_exe) then
!
!   call object_create ('eloc_av')
!   call object_variance_define ('eloc_av', 'eloc_var')
!
!   call object_needed ('eloc')
!
!   return
!
!  endif
!
!! begin
!
!
!  end subroutine eloc_av_bld

! ==============================================================================
  subroutine eloc_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of local energy
!
! Created       : J. Toulouse, 02 Nov 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('eloc_var')

   call object_needed ('eloc_sq_av')
   call object_needed ('eloc_av')

   return

  endif

! begin

! allocations
  call object_associate ('eloc_var', eloc_var)

  eloc_var = eloc_sq_av - eloc_av**2

  end subroutine eloc_var_bld

! ==============================================================================
  subroutine sigma_bld
! ------------------------------------------------------------------------------
! Description   : standard deviation of local energy
!
! Created       : J. Toulouse, 02 Nov 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('sigma')

   call object_needed ('eloc_var')

   return

  endif

! begin

! allocations
  call object_associate ('sigma', sigma)

  sigma = dsqrt(eloc_var)

  end subroutine sigma_bld

! ==============================================================================
  subroutine walker_weights_bld
! ------------------------------------------------------------------------------
! Description   : weights of walkers
!
! Created       : J. Toulouse, 14 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer walk_i

! header
  if (header_exe) then

   call object_create ('walker_weights', walker_weights_index)

   call object_needed ('nwalk')

   return

  endif

! begin

! allocations
  call object_alloc ('walker_weights', walker_weights, nwalk)

  if (index(mode, 'vmc') /= 0) then
   walker_weights (:) = 1.d0

  elseif (index(mode, 'dmc') /= 0) then
   call object_provide_by_index (walker_weights_bld_index, wt_index)
   call object_provide_by_index (walker_weights_bld_index, fprod_index) ! fprod is the finite population correction
   do walk_i = 1, nwalk
    walker_weights (walk_i) = wt (walk_i) * fprod
   enddo

  else
   write(6,'(4a)') trim(here),': mode=',trim(mode),' should contain either vmc or dmc.'
   call die (here)

  endif

  end subroutine walker_weights_bld

! ==============================================================================
!  subroutine walker_weights_sum_bld
! ------------------------------------------------------------------------------
! Description   : sum of weights of walkers over step and blocks
!
! Created       : J. Toulouse, 14 Nov 2006
! ------------------------------------------------------------------------------
!  implicit none
!  include 'commons.h'
!
!! local
!  integer walk_i
!
!! header
!  if (header_exe) then
!
!   call object_create ('walker_weights_sum', walker_weights_sum_index)
!
!   return
!
!  endif
!
!! begin
!
!! allocations
!  if (index(mode, 'vmc') /= 0) then
!   call object_provide_by_index (walker_weights_sum_index, block_iterations_nb_index)
!   walker_weights_sum  = nstep_total*block_iterations_nb
!
!  elseif (index(mode, 'dmc') /= 0) then
!   call object_provide_by_index (walker_weights_sum_bld_index, wgcum_index)
!   walker_weights_sum  = wgcum(1)
!
!  else
!   write(6,'(4a)') trim(here),': mode=',trim(mode),' should contain either vmc or dmc.'
!   call die (here)
!
!  endif
!
!  end subroutine walker_weights_sum_bld

end module montecarlo_mod
