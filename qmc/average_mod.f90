module average_mod

  use basic_tools_mod
  use constants_mod
  use variables_mod
  use objects_mod
  use indexes_mod
  use parser_tools_mod
  use strings_tools_mod

! current MC interation number
  integer                      :: step_iterations_nb = 0
  integer                      :: block_iterations_nb = 0
  real(dp)                     :: total_iterations_nb = 0.d0
  real(dp)                     :: total_iterations_block_nb = 0.d0
  integer                      :: current_walker 
  real(dp)                     :: current_walker_weight

  integer                      :: average_routines_nb = 0
  integer, allocatable         :: average_routines_index (:)

  integer                      :: block_averages_defined_nb = 0
  integer, allocatable         :: block_averages_defined_object_index (:)
  integer, allocatable         :: block_averages_defined_object_bav_index (:)
  integer                      :: block_averages_nb = 0
  integer, allocatable         :: block_averages_object_index (:)
  integer, allocatable         :: block_averages_object_bav_index (:)

  integer                      :: averages_defined_nb = 0
  integer, allocatable         :: averages_defined_object_index (:)
  integer, allocatable         :: averages_defined_object_bav_index (:)
  integer, allocatable         :: averages_defined_object_av_index (:)
  integer                      :: averages_nb = 0
  integer, allocatable         :: averages_object_index (:)
  integer, allocatable         :: averages_object_bav_index (:)
  integer, allocatable         :: averages_object_av_index (:)
  integer                      :: averages_walk_nb = 0
  integer, allocatable         :: averages_walk_object_index (:)
  integer, allocatable         :: averages_walk_object_bav_index (:)
  integer, allocatable         :: averages_walk_object_av_index (:)

  integer                      :: variances_defined_nb = 0
  integer, allocatable         :: variances_defined_object_av_index (:)
  integer, allocatable         :: variances_defined_object_var_index (:)
  integer                      :: variances_nb = 0
  integer, allocatable         :: variances_object_av_index (:)
  integer, allocatable         :: variances_object_var_index (:)

  integer                      :: covariances_defined_nb = 0
  integer, allocatable         :: covariances_defined_object_av1_index (:)
  integer, allocatable         :: covariances_defined_object_av2_index (:)
  integer, allocatable         :: covariances_defined_object_covar_index (:)
  integer                      :: covariances_nb = 0
  integer, allocatable         :: covariances_object_av1_index (:)
  integer, allocatable         :: covariances_object_av2_index (:)
  integer, allocatable         :: covariances_object_covar_index (:)

  integer                      :: errors_defined_nb = 0
  integer, allocatable         :: errors_defined_object_av_index (:)
  integer, allocatable         :: errors_defined_object_var_index (:)
  integer, allocatable         :: errors_defined_object_err_index (:)
  integer                      :: errors_nb = 0
  integer, allocatable         :: errors_object_av_index (:)
  integer, allocatable         :: errors_object_var_index (:)
  integer, allocatable         :: errors_object_err_index (:)

  integer                      :: block_averages_nb_save
  integer, allocatable         :: block_averages_object_index_save (:)
  integer, allocatable         :: block_averages_object_bav_index_save (:)
  integer                      :: averages_nb_save
  integer, allocatable         :: averages_object_index_save (:)
  integer, allocatable         :: averages_object_bav_index_save (:)
  integer, allocatable         :: averages_object_av_index_save (:)
  integer                      :: variances_nb_save
  integer, allocatable         :: variances_object_av_index_save (:)
  integer, allocatable         :: variances_object_var_index_save (:)
  integer                      :: covariances_nb_save
  integer, allocatable         :: covariances_object_av1_index_save (:)
  integer, allocatable         :: covariances_object_av2_index_save (:)
  integer, allocatable         :: covariances_object_covar_index_save (:)
  integer                      :: errors_nb_save
  integer, allocatable         :: errors_object_av_index_save (:)
  integer, allocatable         :: errors_object_var_index_save (:)
  integer, allocatable         :: errors_object_err_index_save (:)
  integer                      :: average_routines_nb_save
  integer, allocatable         :: average_routines_index_save (:)

  contains

!===========================================================================
  subroutine average_menu
!---------------------------------------------------------------------------
! Description : menu for calculating averages and statistical errors
!
! Created     : J. Toulouse, 15 Nov 2005
!---------------------------------------------------------------------------
  implicit none

! local
  character (len=max_string_len_rout), save :: lhere = 'average_menu'
  integer block_averages_list_nb, averages_list_nb, errors_list_nb, obj_i
# if defined (PATHSCALE)
   character (len=max_string_len) :: block_averages_list (max_string_array_len) ! for pathscale compiler
   character (len=max_string_len) :: averages_list (max_string_array_len) ! for pathscale compiler
   character (len=max_string_len) :: errors_list (max_string_array_len) ! for pathscale compiler
# else
   character (len=max_string_len), allocatable :: block_averages_list (:)
   character (len=max_string_len), allocatable :: averages_list (:)
   character (len=max_string_len), allocatable :: errors_list (:)
# endif

! begin
  block_averages_list_nb = 0
  averages_list_nb = 0
  errors_list_nb = 0

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for average menu:'
   write(6,'(a)') 'average'
   write(6,'(a)') '  block_averages ... end: list of block averages to compute'
   write(6,'(a)') '  averages ... end: list of global averages to compute'
   write(6,'(a)') '  errors ... end: list of statistical errors to compute'
   write(6,'(a)') 'end'
   write(6,*)

  case ('block_averages')
# if defined (PATHSCALE)
   call get_next_value_list_string ('block_averages_list', block_averages_list, block_averages_list_nb) ! for pathscale compiler
# else
   call get_next_value_list ('block_averages_list', block_averages_list, block_averages_list_nb)
# endif

  case ('averages')
# if defined (PATHSCALE)
   call get_next_value_list_string ('averages_list', averages_list, averages_list_nb) ! for pathscale compiler
# else
   call get_next_value_list ('averages_list', averages_list, averages_list_nb)
# endif

  case ('errors')
# if defined (PATHSCALE)
   call get_next_value_list_string ('errors_list', errors_list, errors_list_nb) ! for pathscale compiler
# else
   call get_next_value_list ('errors_list', errors_list, errors_list_nb)
# endif

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

! request block averages to compute
  do obj_i = 1, block_averages_list_nb
   call object_block_average_request (block_averages_list (obj_i))
  enddo

! request averages to compute
  do obj_i = 1, averages_list_nb
   call object_average_request (averages_list (obj_i))
  enddo

! request statistical errors to compute
  do obj_i = 1, errors_list_nb
   call object_error_request (errors_list (obj_i))
  enddo

  end subroutine average_menu

! ===================================================================================
  subroutine routine_average (routine_name)
! -----------------------------------------------------------------------------------
! Description   : define routine computing the average of an object
! Description   : store index of routine
!
! Created       : J. Toulouse, 20 Dec 2005
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! local
  character(len=max_string_len_rout), save :: lhere = 'routine_average'
  integer routine_ind, rtn_i

! begin

! index of routine
  routine_ind = routine_index (routine_name)
  if (routine_ind == 0) then
    call die (lhere, 'routine >'+trim(routine_name)+'< is not catalogued.')
  endif

! if routine already defined as average routine, do nothing
  do rtn_i = 1, average_routines_nb
   if (routine_ind == average_routines_index (rtn_i)) then
       return
   endif
  enddo

  average_routines_nb = average_routines_nb + 1
  call alloc ('average_routines_index', average_routines_index, average_routines_nb)
  average_routines_index (average_routines_nb) = routine_ind

  write(6,'(4a)') trim(lhere),': ', trim(routine_name),' defined as average routine'

 end subroutine routine_average

! ===================================================================================
  subroutine object_block_average_define (object_name, object_bav_name, l_unweighted)
! -----------------------------------------------------------------------------------
! Description   : define block average of object to be computed in MC iterations
! Description   : store indexes of couple (object, block average)
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name, object_bav_name
  logical, optional :: l_unweighted

! local
  character(len=max_string_len_rout), save :: lhere = 'object_block_average_define'
  integer object_ind, object_bav_ind, obj_i

! begin

! get index of object, catalogue it if necessary
  call object_add_once_and_index (object_name, object_ind)

!!!comment out because of unweighted averages
! test if block average has already been defined through a global average
!!!  do obj_i = 1, averages_defined_nb
!!!   if (object_ind == averages_defined_object_index (obj_i)) then
!!!     write(6,'(6a)') trim(lhere), ': block average >'+object_bav_name+'< of object >'+object_name+'< has already been automatically defined through the corresponding global average.'
!!!     write(6,'(2a)') trim(lhere), ': To define the block average, object_block_average_define has to be placed before object_average_define.'
!!!     call die (lhere)
!!!   endif
!!!  enddo

! get index of block average, catalogue it if necessary
  call object_add_once_and_index (object_bav_name, object_bav_ind)

! test if object and its block average are different
  if (object_ind == object_bav_ind) then
    call die (lhere, 'object >'+trim(objects(object_bav_ind)%name)+'< is defined as its own block average!')
  endif

! test if block average not already defined
  do obj_i = 1, block_averages_defined_nb
   if (object_bav_ind == block_averages_defined_object_bav_index (obj_i)) then
    call die (lhere, 'block average object >'+trim(objects(object_bav_ind)%name)+'< defined more than once.')
   endif
  enddo

  block_averages_defined_nb = block_averages_defined_nb + 1
  call alloc ('block_averages_defined_object_index', block_averages_defined_object_index, block_averages_defined_nb)
  call alloc ('block_averages_defined_object_bav_index', block_averages_defined_object_bav_index, block_averages_defined_nb)
  block_averages_defined_object_index (block_averages_defined_nb) = object_ind
  block_averages_defined_object_bav_index (block_averages_defined_nb) = object_bav_ind

! for unweighted averages
  if (present(l_unweighted)) then
   objects(object_bav_ind)%unweighted=l_unweighted
  endif

 end subroutine object_block_average_define

! ===================================================================================
  subroutine object_block_average_request (object_bav_name)
! -----------------------------------------------------------------------------------
! Description   : request calculation of block average object
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_bav_name

! local
  integer object_bav_ind

! begin

! index of block average
  object_bav_ind = object_index_or_die (object_bav_name)

  call object_block_average_request_by_index (object_bav_ind)

 end subroutine object_block_average_request

! ===================================================================================
  subroutine object_block_average_request_by_index (object_bav_ind)
! -----------------------------------------------------------------------------------
! Description   : request calculation of block average object by its index
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_bav_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_block_average_request_by_index'
  integer object_ind, obj_i
  logical object_found

! begin

! test if object is a defined block average object
  object_found = .false.
  do obj_i = 1, block_averages_defined_nb
   if (object_bav_ind == block_averages_defined_object_bav_index (obj_i)) then
    object_found = .true.
    object_ind = block_averages_defined_object_index (obj_i)
    exit
   endif
  enddo
  if (.not. object_found) then
    call die (lhere, ': object >'+trim(objects(object_bav_ind)%name)+'< has never been defined as a block average.')
  endif

! do nothing if object already recorded for block average
  do obj_i = 1, block_averages_nb
   if (object_bav_ind == block_averages_object_bav_index (obj_i)) then
     return
   endif
  enddo

! add block average to the list of block averages
  block_averages_nb = block_averages_nb + 1
  call alloc ('block_averages_object_index', block_averages_object_index, block_averages_nb)
  call alloc ('block_averages_object_bav_index', block_averages_object_bav_index, block_averages_nb)
  block_averages_object_index (block_averages_nb) = object_ind
  block_averages_object_bav_index (block_averages_nb) = object_bav_ind

! invalidate block average
  call object_invalidate_by_index (object_bav_ind)

 end subroutine object_block_average_request_by_index

! ===================================================================================
  subroutine object_average_define (object_name, object_av_name, l_unweighted)
! -----------------------------------------------------------------------------------
! Description   : define couple (object, average of object)
!
! Created       : J. Toulouse, 15 Jan 2006
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name, object_av_name
  logical, optional :: l_unweighted

! local
  character(len=max_string_len_rout), save :: lhere = 'object_define_average'
  integer object_ind, object_bav_ind, object_av_ind, obj_i
  logical block_average_found

! begin

! index of object, catalogue it if necessary
  call object_add_once_and_index (object_name, object_ind)

! index of average, catalogue it if necessary
  call object_add_once_and_index (object_av_name, object_av_ind)

! test if object and its average are different
  if (object_ind == object_av_ind) then
    call die (lhere, 'object >'+trim(objects(object_av_ind)%name)+'< is defined as its own object average!')
  endif

! test if average not already defined
  do obj_i = 1, averages_defined_nb
   if (object_av_ind == averages_defined_object_av_index (obj_i)) then
    call die (lhere, 'object average >'+trim(objects(object_av_ind)%name)+'< defined more than once.')
   endif
  enddo

! for unweighted averages
  if (present(l_unweighted)) then
   objects(object_av_ind)%unweighted=l_unweighted
  endif

! block average !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test if a corresponding block average is defined (must check also whether they are both weighted or unweighted)
  block_average_found = .false.
  do obj_i = 1, block_averages_defined_nb
   if (object_ind == block_averages_defined_object_index (obj_i) .and. &
      (objects(object_av_ind)%unweighted .eqv. objects(block_averages_defined_object_bav_index(obj_i))%unweighted)) then
     block_average_found = .true.
     object_bav_ind = block_averages_defined_object_bav_index (obj_i)
     exit
   endif
  enddo

! if corresponding block average not found, add it as a new object and define as a block average
  if (.not. block_average_found) then
   call object_add_and_index (object_bav_ind)  
   call object_block_average_define (object_name, objects(object_bav_ind)%name, l_unweighted)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  averages_defined_nb = averages_defined_nb + 1
  call alloc ('averages_defined_object_index', averages_defined_object_index, averages_defined_nb)
  call alloc ('averages_defined_object_bav_index', averages_defined_object_bav_index, averages_defined_nb) ! block average
  call alloc ('averages_defined_object_av_index', averages_defined_object_av_index, averages_defined_nb)
  averages_defined_object_index (averages_defined_nb) = object_ind
  averages_defined_object_bav_index (averages_defined_nb) = object_bav_ind  ! block average
  averages_defined_object_av_index (averages_defined_nb) = object_av_ind

!  write(6,'(5a)') trim(lhere),': ', trim(objects(object_av_ind)%name),' is average of ',trim(objects(object_ind)%name)

  end subroutine object_average_define

! ===================================================================================
  subroutine object_average_define_from_block_average (object_bav_name, object_av_name)
! -----------------------------------------------------------------------------------
! Description   : define global average associated to a block average
!
! Created       : J. Toulouse, 20 Apr 2008
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_bav_name, object_av_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_average_define_from_block_average'
  integer  object_bav_ind, object_av_ind, obj_i


! begin

! index of block average, catalogue it if necessary
  call object_add_once_and_index (object_bav_name, object_bav_ind)

! index of average, catalogue it if necessary
  call object_add_once_and_index (object_av_name, object_av_ind)

! test if block average and average are different
  if (object_bav_ind == object_av_ind) then
    call die (lhere, 'block average >'+trim(objects(object_bav_ind)%name)+'< is identical to its associated global average!')
  endif

! test if average not already defined
  do obj_i = 1, averages_defined_nb
   if (object_av_ind == averages_defined_object_av_index (obj_i)) then
    call die (lhere, 'average object >'+trim(objects(object_av_ind)%name)+'< defined more than once.')
   endif
  enddo

  averages_defined_nb = averages_defined_nb + 1
  call alloc ('averages_defined_object_index', averages_defined_object_index, averages_defined_nb)
  call alloc ('averages_defined_object_bav_index', averages_defined_object_bav_index, averages_defined_nb)
  call alloc ('averages_defined_object_av_index', averages_defined_object_av_index, averages_defined_nb)
  averages_defined_object_index (averages_defined_nb) = 0 ! no object
  averages_defined_object_bav_index (averages_defined_nb) = object_bav_ind
  averages_defined_object_av_index (averages_defined_nb) = object_av_ind

  end subroutine object_average_define_from_block_average

! ===================================================================================
  subroutine object_average_unweighted_define (object_name, object_av_name)
! -----------------------------------------------------------------------------------
! Description   : define couple (object, average of object) for unweigthed averages
! Description   : in DMC. 
!
! Created       : J. Toulouse, 07 Jul 2010
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name, object_av_name

! local
  integer object_av_ind

  call object_average_define (object_name, object_av_name, .true.)

  end subroutine object_average_unweighted_define

! ===================================================================================
  subroutine object_average_walk_define (object_name, object_av_name)
! -----------------------------------------------------------------------------------
! Description   : define couple (object, average of object) when object is given for
! Description   : a list of walkers
!
! Created       : J. Toulouse, 11 Nov 2006
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_name, object_av_name

! local
  integer object_ind

  call object_average_define (object_name, object_av_name)

  object_ind = object_index (object_name)
  objects(object_ind)%walkers = .true.

  end subroutine object_average_walk_define

! ===================================================================================
  subroutine object_average_request (object_av_name)
! -----------------------------------------------------------------------------------
! Description   : request for calculation of average object 'object_av_name' over MC iterations
!
! Created       : J. Toulouse, 15 Jan 2006
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_av_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_average_request'
  integer object_ind, object_bav_ind, object_av_ind, obj_i
  logical object_found

! begin

! index of average
  object_av_ind = object_index_or_die (object_av_name)

! test if object is a defined average object
  object_found = .false.
  do obj_i = 1, averages_defined_nb
   if (object_av_ind == averages_defined_object_av_index (obj_i)) then
    object_found = .true.
    object_ind = averages_defined_object_index (obj_i)
    object_bav_ind = averages_defined_object_bav_index (obj_i) ! block average
    exit
   endif
  enddo
  if (.not. object_found) then
    call die (lhere, 'object >'+trim(objects(object_av_ind)%name)+'< has never been defined as an average.')
  endif

! do nothing if object already recorded for average
  do obj_i = 1, averages_nb
   if (object_av_ind == averages_object_av_index (obj_i) ) then
     return
   endif
  enddo
  do obj_i = 1, averages_walk_nb
   if (object_av_ind == averages_walk_object_av_index (obj_i)) then
     return
   endif
  enddo

! add average to the list of averages
  if (.not. objects(object_ind)%walkers) then
   averages_nb = averages_nb + 1
   call alloc ('averages_object_index', averages_object_index, averages_nb)
   call alloc ('averages_object_bav_index', averages_object_bav_index, averages_nb) ! block average
   call alloc ('averages_object_av_index', averages_object_av_index, averages_nb)
   averages_object_index (averages_nb) = object_ind
   averages_object_bav_index (averages_nb) = object_bav_ind ! block average
   averages_object_av_index (averages_nb) = object_av_ind
  else
   averages_walk_nb = averages_walk_nb + 1
   call alloc ('averages_walk_object_index', averages_walk_object_index, averages_walk_nb)
   call alloc ('averages_walk_object_bav_index', averages_walk_object_bav_index, averages_walk_nb)
   call alloc ('averages_walk_object_av_index', averages_walk_object_av_index, averages_walk_nb)
   averages_walk_object_index (averages_walk_nb) = object_ind
   averages_walk_object_bav_index (averages_walk_nb) = object_bav_ind ! block average
   averages_walk_object_av_index (averages_walk_nb) = object_av_ind
  endif

! invalidate average
  call object_invalidate_by_index (object_av_ind)

! request corresponding block average, except if block average is an independent object or walk average
  if (object_ind /= 0 .and. .not. objects(object_ind)%walkers) then
   call object_block_average_request_by_index (object_bav_ind)
  endif

 end subroutine object_average_request

! ===================================================================================
  subroutine object_variance_define (object_av_name, object_var_name)
! -----------------------------------------------------------------------------------
! Description   : define variance of averaged object to be computed in MC iterations
! Description   : store indexes of couple (average, variance)
!
! Created       : J. Toulouse, 04 Sep 2007
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_av_name, object_var_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_variance_define'
  integer object_av_ind, object_var_ind
  integer obj_i

! begin

! get index of average, catalogue average if necessary
  call object_add_once_and_index (object_av_name, object_av_ind)

! test if variance has already been defined through an statistical error
  do obj_i = 1, errors_defined_nb
   if (object_av_ind == errors_defined_object_av_index (obj_i)) then
     write(6,'(6a)') lhere, ': variance >'+object_var_name+'< of average >'+object_av_name+'< has already been automatically defined though the corresponding statistical error.'
     write(6,'(2a)') lhere, ': To define the variance, object_variance_define has to be placed before object_error_define.'
     call die (lhere)
   endif
  enddo

! get index of variance, catalogue average if necessary
  call object_add_once_and_index (object_var_name, object_var_ind)

! test if average and its variance are different
  if (object_av_ind == object_var_ind) then
    call die (lhere, 'object >'+trim(objects(object_var_ind)%name)+'< is defined as its own variance object!')
  endif

! test if variance not already defined
  do obj_i = 1, variances_defined_nb
   if (object_var_ind == variances_defined_object_var_index (obj_i) ) then
    call die (lhere, 'variance object >'+trim(objects(object_var_ind)%name)+'< defined more than once.')
   endif
  enddo

  variances_defined_nb = variances_defined_nb + 1
  call alloc ('variances_defined_object_av_index', variances_defined_object_av_index, variances_defined_nb)
  call alloc ('variances_defined_object_var_index', variances_defined_object_var_index, variances_defined_nb)
  variances_defined_object_av_index (variances_defined_nb) = object_av_ind
  variances_defined_object_var_index (variances_defined_nb) = object_var_ind

 end subroutine object_variance_define

! ===================================================================================
  subroutine object_variance_request (object_var_name)
! -----------------------------------------------------------------------------------
! Description   : request for calculation of variance object 'object_var_name' over MC iterations
!
! Created       : J. Toulouse, 04 Sep 2007
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_var_name

! local
  integer object_var_ind

! begin

! index of error
  object_var_ind = object_index_or_die (object_var_name)

  call object_variance_request_by_index (object_var_ind)

 end subroutine object_variance_request

! ===================================================================================
  subroutine object_variance_request_by_index (object_var_ind)
! -----------------------------------------------------------------------------------
! Description   : request calculation of variance object by its index
!
! Created       : J. Toulouse, 06 Sep 2007
! -----------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: object_var_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_variance_request_by_index'
  integer object_av_ind, obj_i
  logical object_found

! begin

! test if object is a defined variance object
  object_found = .false.
  do obj_i = 1, variances_defined_nb
   if (object_var_ind == variances_defined_object_var_index (obj_i)) then
    object_found = .true.
    object_av_ind = variances_defined_object_av_index (obj_i)
    exit
   endif
  enddo
  if (.not. object_found) then
    call die (lhere, ': object >'+trim(objects(object_var_ind)%name)+'< has never been defined as a variance.')
  endif

! do nothing if object already recorded for variance
  do obj_i = 1, variances_nb
   if (object_var_ind == variances_object_var_index (obj_i)) then
     return
   endif
  enddo

! add variance to the list of variances
  variances_nb = variances_nb + 1
  call alloc ('variances_object_av_index', variances_object_av_index, variances_nb)
  call alloc ('variances_object_var_index', variances_object_var_index, variances_nb)
  variances_object_av_index (variances_nb) = object_av_ind
  variances_object_var_index (variances_nb) = object_var_ind

! invalidate variance
  call object_invalidate_by_index (object_var_ind)

 end subroutine object_variance_request_by_index

! ===================================================================================
  subroutine object_covariance_define (object_av1_name, object_av2_name, object_covar_name)
! -----------------------------------------------------------------------------------
! Description   : define covariance of two averaged objects to be computed in MC iterations
!
! Created       : J. Toulouse, 07 Sep 2007
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_av1_name, object_av2_name, object_covar_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_covariance_define'
  integer object_av1_ind, object_av2_ind, object_covar_ind
  integer obj_i

! begin

! get index of average, catalogue average if necessary
  call object_add_once_and_index (object_av1_name, object_av1_ind)
  call object_add_once_and_index (object_av2_name, object_av2_ind)

! get index of variance, catalogue average if necessary
  call object_add_once_and_index (object_covar_name, object_covar_ind)

! test if covariance not already defined
  do obj_i = 1, covariances_defined_nb
   if (object_covar_ind == covariances_defined_object_covar_index (obj_i)) then
    call die (lhere, 'covariance object >'+trim(objects(object_covar_ind)%name)+'< defined more than once.')
   endif
  enddo

  covariances_defined_nb = covariances_defined_nb + 1
  call alloc ('covariances_defined_object_av1_index', covariances_defined_object_av1_index, covariances_defined_nb)
  call alloc ('covariances_defined_object_av2_index', covariances_defined_object_av2_index, covariances_defined_nb)
  call alloc ('covariances_defined_object_covar_index', covariances_defined_object_covar_index, covariances_defined_nb)
  covariances_defined_object_av1_index (covariances_defined_nb) = object_av1_ind
  covariances_defined_object_av2_index (covariances_defined_nb) = object_av2_ind
  covariances_defined_object_covar_index (covariances_defined_nb) = object_covar_ind

 end subroutine object_covariance_define

! ===================================================================================
  subroutine object_covariance_request (object_covar_name)
! -----------------------------------------------------------------------------------
! Description   : request for calculation of covariance object over MC iterations
!
! Created       : J. Toulouse, 17 Jul 2007
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_covar_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_covariance_request'
  integer object_av1_ind, object_av2_ind, object_covar_ind
  integer obj_i
  logical object_found

! begin

! index of covariance
  object_covar_ind = object_index_or_die (object_covar_name)

! test if object is a defined covariance
  object_found = .false.
  do obj_i = 1, errors_defined_nb
   if (object_covar_ind == covariances_defined_object_covar_index (obj_i)) then
    object_found = .true.
    object_av1_ind = covariances_defined_object_av1_index (obj_i)
    object_av2_ind = covariances_defined_object_av2_index (obj_i)
    exit
   endif
  enddo
  if (.not. object_found) then
    call die (lhere, ': object >'+trim(objects(object_covar_ind)%name)+'< has never been defined as a covariance.')
  endif

! do nothing if object already recorded for covariance
  do obj_i = 1, covariances_nb
   if (object_covar_ind == covariances_object_covar_index (obj_i) ) then
     return
   endif
  enddo

! add covariance to the list of covariances
  covariances_nb = covariances_nb + 1
  call alloc ('covariances_object_av1_index', covariances_object_av1_index, covariances_nb)
  call alloc ('covariances_object_av2_index', covariances_object_av2_index, covariances_nb)
  call alloc ('covariances_object_covar_index', covariances_object_covar_index, covariances_nb)
  covariances_object_av1_index (covariances_nb) = object_av1_ind
  covariances_object_av2_index (covariances_nb) = object_av2_ind
  covariances_object_covar_index (covariances_nb) = object_covar_ind

! invalidate covariance
  call object_invalidate_by_index (object_covar_ind)

 end subroutine object_covariance_request

! ===================================================================================
  subroutine object_error_define (object_av_name, object_err_name)
! -----------------------------------------------------------------------------------
! Description   : define statistical error of averaged object to be computed in MC iterations
! Description   : store indexes of couple (average, error)
!
! Created       : J. Toulouse, 20 Oct 2005
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_av_name, object_err_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_error_define'
  integer object_av_ind, object_var_ind, object_err_ind
  integer obj_i
  logical variance_found

! begin

! index of average, catalogue object if necessary
  call object_add_once_and_index (object_av_name, object_av_ind)

! index of error, catalogue object if necessary
  call object_add_once_and_index (object_err_name, object_err_ind)

! test if average and its error are different
  if (object_av_ind == object_err_ind) then
    call die (lhere, 'object >'+trim(objects(object_av_ind)%name)+'< is defined as its own object error!')
  endif

! test if error not already defined
  do obj_i = 1, errors_defined_nb
   if (object_err_ind == errors_defined_object_err_index (obj_i) ) then
    call die (lhere, 'statistical error object >'+trim(objects(object_err_ind)%name)+'< defined more than once.')
   endif
  enddo

! test if a corresponding variance is defined
  variance_found = .false.
  do obj_i = 1, variances_defined_nb
   if (object_av_ind == variances_defined_object_av_index (obj_i)) then
     variance_found = .true.
     object_var_ind = variances_defined_object_var_index (obj_i)
     exit
   endif
  enddo

! if corresponding variance not found, add it as a new object and define as a variance
  if (.not. variance_found) then
   call object_add_and_index (object_var_ind)
   call object_variance_define (object_av_name, objects(object_var_ind)%name)
  endif

  errors_defined_nb = errors_defined_nb + 1
  call alloc ('errors_defined_object_av_index', errors_defined_object_av_index, errors_defined_nb)
  call alloc ('errors_defined_object_var_index', errors_defined_object_var_index, errors_defined_nb)
  call alloc ('errors_defined_object_err_index', errors_defined_object_err_index, errors_defined_nb)
  errors_defined_object_av_index (errors_defined_nb) = object_av_ind
  errors_defined_object_var_index (errors_defined_nb) = object_var_ind
  errors_defined_object_err_index (errors_defined_nb) = object_err_ind

  end subroutine object_error_define

! ===================================================================================
  subroutine object_error_define_from_variance (object_var_name, object_err_name)
! -----------------------------------------------------------------------------------
! Description   : define statistical error associated to a variance
!
! Created       : J. Toulouse, 20 Sep 2007
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_var_name, object_err_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_error_define_from_variance'
  integer  object_var_ind, object_err_ind, obj_i

! begin

! index of variance, catalogue object if necessary
  call object_add_once_and_index (object_var_name, object_var_ind)

! index of error, catalogue object if necessary
  call object_add_once_and_index (object_err_name, object_err_ind)

! test if variance and error are different
  if (object_var_ind == object_err_ind) then
    call die (lhere, 'variance >'+trim(objects(object_var_ind)%name)+'< is identical to its associated statistical error!')
  endif

! test if error not already defined
  do obj_i = 1, errors_defined_nb
   if (object_err_ind == errors_defined_object_err_index (obj_i) ) then
    call die (lhere, 'statistical error object >'+trim(objects(object_err_ind)%name)+'< defined more than once.')
   endif
  enddo

  errors_defined_nb = errors_defined_nb + 1
  call alloc ('errors_defined_object_av_index', errors_defined_object_av_index, errors_defined_nb)
  call alloc ('errors_defined_object_var_index', errors_defined_object_var_index, errors_defined_nb)
  call alloc ('errors_defined_object_err_index', errors_defined_object_err_index, errors_defined_nb)
  errors_defined_object_av_index (errors_defined_nb) = 0 ! no average object
  errors_defined_object_var_index (errors_defined_nb) = object_var_ind
  errors_defined_object_err_index (errors_defined_nb) = object_err_ind

 end subroutine object_error_define_from_variance

! ===================================================================================
  subroutine object_error_request (object_err_name)
! -----------------------------------------------------------------------------------
! Description   : request for calculation of statistical error object 'object_err_name' over MC iterations
!
! Created       : J. Toulouse, 15 Jan 2006
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: object_err_name

! local
  character(len=max_string_len_rout), save :: lhere = 'object_error_request'
  integer object_av_ind, object_var_ind, object_err_ind, obj_i
  logical object_found

! begin

! index of error
  object_err_ind = object_index_or_die (object_err_name)

! test if object is a defined statistical error
  object_found = .false.
  do obj_i = 1, errors_defined_nb
   if (object_err_ind == errors_defined_object_err_index (obj_i)) then
    object_found = .true.
    object_av_ind = errors_defined_object_av_index (obj_i)
    object_var_ind = errors_defined_object_var_index (obj_i)
    exit
   endif
  enddo
  if (.not. object_found) then
    call die (lhere, ': object >'+trim(objects(object_err_ind)%name)+'< has never been defined as a statistical error.')
  endif

! do nothing if object already recorded for error
  do obj_i = 1, errors_nb
   if (object_err_ind == errors_object_err_index (obj_i) ) then
     return
   endif
  enddo


! add error to the list of erros
  errors_nb = errors_nb + 1
  call alloc ('errors_object_av_index', errors_object_av_index, errors_nb)
  call alloc ('errors_object_var_index', errors_object_var_index, errors_nb)
  call alloc ('errors_object_err_index', errors_object_err_index, errors_nb)
  errors_object_av_index (errors_nb) = object_av_ind
  errors_object_var_index (errors_nb) = object_var_ind
  errors_object_err_index (errors_nb) = object_err_ind

! invalidate error
  call object_invalidate_by_index (object_err_ind)

! request corresponding variance, except if variance is an independent object
  if (object_av_ind /= 0) then
   call object_variance_request_by_index (object_var_ind)
  endif

 end subroutine object_error_request

!! ===================================================================================
!  subroutine object_average_by_index_double_0 (object_ind, object_av_ind)
!! -----------------------------------------------------------------------------------
!! Description   : calculate average of object over MC iterations
!! Description   : no longer used
!!
!! Created       : J. Toulouse, 18 Oct 2005
!! Modified      : J. Toulouse, 05 Feb 2006: use fortran pointers, warning association not checked!
!! Modified      : J. Toulouse, 03 Oct 2006: MPI version
!! -----------------------------------------------------------------------------------
!  include 'modules.h'
!  implicit none
!
!! input
!  integer object_ind, object_av_ind
!
!! local
!  character(len=max_string_len_rout), save :: lhere = 'object_average_by_index_double_0'
!  character(len=max_string_len_type)   :: object_type, object_av_type
!
!# if defined (MPI)
!  real(dp) collect
!  integer ierr
!# endif
!
!! begin
!
!! only for first iteration
!  if (step_iterations_nb == 1) then
!
!!   test association
!    call object_associated_or_die_by_index (object_ind)
!    call object_associated_or_die_by_index (object_av_ind)
!
!!   test on type
!    object_type = objects(object_ind)%type
!    object_av_type = objects(object_av_ind)%type
!
!    if (object_type /= object_av_type) then
!     write(6,'(5a)') trim(lhere),': type of object ',trim(objects(object_ind)%name),' is ', object_type
!     write(6,'(5a)') trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is ', object_av_type
!     write(6,'(2a)') trim(lhere),': they should be identical'
!     call die (lhere)
!    endif
!
!  endif ! first iteration
!
!! invalidate by default the average
!  call object_invalidate_by_index (object_av_ind)
!
!! reinitialization of sum for each block
!  if (mod(step_iterations_nb - 1, nstep) == 0) then
!     objects(object_ind)%sum_double_0 = 0.d0
!  endif
!
!! at each step:
!
!! for one-electron move version
!!  if (index (trim(mode), 'mov1') /= 0) then
!  objects(object_ind)%sum_double_0 = objects(object_ind)%sum_double_0 + objects(object_ind)%pointer_double_0
!
!! for all-electron move version
!!   else
!!   objects(object_ind)%sum_double_0 = objects(object_ind)%sum_double_0 + prob_acc*objects(object_ind)%pointer_double_0 + prob_rej*objects(object_ind)%previous_double_0
!!   write(6,*) 'objects(object_ind)%sum_double_0=',objects(object_ind)%sum_double_0
!!   objects(object_ind)%previous_double_0 = objects(object_ind)%pointer_double_0
!
!!  endif
!
!!  at the end of each block
!   if (mod(step_iterations_nb , nstep) == 0) then
!
!!    initialization for first block
!     if (block_iterations_nb == 1 ) then
!       objects(object_ind)%sum_blk_double_0 = 0.d0
!     endif
!
!# if defined (MPI)
!!    sum values from all processes
!     call mpi_allreduce(objects(object_ind)%sum_double_0,collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
!     if (ierr /= 0) then
!        call die (lhere, 'error in mpi_allreduce')
!     endif
!     objects(object_ind)%sum_double_0 = collect
!# endif
!
!     objects(object_ind)%sum_blk_double_0 = objects(object_ind)%sum_blk_double_0 + objects(object_ind)%sum_double_0/nstep_total
!
!!    calculate average
!     objects(object_av_ind)%pointer_double_0  = objects(object_ind)%sum_blk_double_0 / block_iterations_nb
!     call object_modified_by_index (object_av_ind)
!
!   endif
!
! end subroutine object_average_by_index_double_0
!
!! ===================================================================================
!  subroutine object_average_by_index_double_1 (object_ind, object_av_ind)
!! -----------------------------------------------------------------------------------
!! Description   : calculate average of object over MC iterations
!!
!! Created       : J. Toulouse, 18 Oct 2005
!! Modified      : J. Toulouse, 05 Feb 2006: use fortran pointers, warning association not checked!
!! -----------------------------------------------------------------------------------
!  include 'modules.h'
!  implicit none
!
!! input
!  integer object_ind, object_av_ind
!
!! local
!  character(len=max_string_len_rout), save :: lhere = 'object_average_by_index_double_1'
!  character(len=max_string_len_type)   :: object_type, object_av_type
!  integer dim1, dim_av1
!
!# if defined (MPI)
!  real(dp), allocatable :: collect(:)
!  integer ierr
!# endif
!
!
!! begin
!
!! only for first iteration
!  if (step_iterations_nb == 1) then
!
!!   test association
!    call object_associated_or_die_by_index (object_ind)
!    call object_associated_or_die_by_index (object_av_ind)
!
!!   test on type
!    object_type = objects(object_ind)%type
!    object_av_type = objects(object_av_ind)%type
!
!    if (object_type /= object_av_type) then
!     write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is ', object_type
!     write(6,*) trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is ', object_av_type
!     write(6,*) trim(lhere),': they should be identical'
!     call die (lhere)
!    endif
!
!!   test on dimensions
!    dim1 = objects(object_ind)%dimensions(1)
!    dim_av1 = objects(object_av_ind)%dimensions(1)
!
!    if (dim1 /= dim_av1) then
!     write(6,'(4a,i8)') trim(lhere),': dimension of object', trim(objects(object_ind)%name),' is ', dim1
!     write(6,'(4a,i8)') trim(lhere),': dimension of object ',trim(objects(object_av_ind)%name),' is ', dim_av1
!     write(6,'(2a)') trim(lhere),': they should be identical'
!     call die (lhere)
!    endif
!
!!  initialization
!   call alloc ('objects(object_ind)%sum_double_1', objects(object_ind)%sum_double_1, dim1)
!
!! for all-electron move version
!!   if (index (trim(mode), 'mov1') == 0 ) then
!!       call alloc ('objects(object_ind)%previous_double_1', objects(object_ind)%previous_double_1, dim1)
!!   endif
!
!  endif ! first iteration
!
!! invalidate by default the average
!  call object_invalidate_by_index (object_av_ind)
!
!! reinitialization of sum for each block
!  if (mod(step_iterations_nb - 1, nstep) == 0) then
!     objects(object_ind)%sum_double_1 = 0.d0
!  endif
!
!! at each step
!
!! for one-electron move version
!!  if (index (trim(mode), 'mov1') /= 0) then
!  objects(object_ind)%sum_double_1 = objects(object_ind)%sum_double_1 + objects(object_ind)%pointer_double_1
!
!! for all-electron move version
!!   else
!!   objects(object_ind)%sum_double_1 = objects(object_ind)%sum_double_1 + prob_acc*objects(object_ind)%pointer_double_1 + prob_rej*objects(object_ind)%previous_double_1
!!   objects(object_ind)%previous_double_1 = objects(object_ind)%pointer_double_1
!!
!!  endif
!
!!  at the end of each block
!   if (mod(step_iterations_nb , nstep) == 0) then
!
!!    initialization for first block
!     if (block_iterations_nb == 1 ) then
!       call alloc ('objects(object_ind)%sum_blk_double_1', objects(object_ind)%sum_blk_double_1, objects(object_ind)%dimensions(1))
!       objects(object_ind)%sum_blk_double_1 = 0.d0
!     endif
!
!# if defined (MPI)
!!    sum values from all processes
!     call alloc ('collect', collect, objects(object_ind)%dimensions(1))
!     call mpi_allreduce(objects(object_ind)%sum_double_1,collect,objects(object_ind)%dimensions(1),mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
!     if (ierr /= 0) then
!        call die (lhere, 'error in mpi_allreduce')
!     endif
!     objects(object_ind)%sum_double_1 (:) = collect (:)
!# endif
!
!     objects(object_ind)%sum_blk_double_1 = objects(object_ind)%sum_blk_double_1 + objects(object_ind)%sum_double_1/nstep_total
!
!!    calculate average
!     objects(object_av_ind)%pointer_double_1 = objects(object_ind)%sum_blk_double_1 / block_iterations_nb
!     call object_modified_by_index (object_av_ind)
!
!   endif
!
! end subroutine object_average_by_index_double_1
!
!! ===================================================================================
!  subroutine object_average_by_index_double_2 (object_ind, object_av_ind)
!! -----------------------------------------------------------------------------------
!! Description   : calculate average of object over MC iterations
!!
!! Created       : J. Toulouse, 18 Oct 2005
!! Modified      : J. Toulouse, 05 Feb 2006: use fortran pointers, warning association not checked!
!! -----------------------------------------------------------------------------------
!  include 'modules.h'
!  implicit none
!
!! input
!  integer object_ind, object_av_ind
!
!! local
!  character(len=max_string_len_rout), save :: lhere = 'object_average_by_index_double_2'
!  character(len=max_string_len_type)   :: object_type, object_av_type
!  integer dim1, dim_av1, dim2, dim_av2
!
!
!# if defined (MPI)
!  real(dp), allocatable :: collect(:,:)
!  integer ierr
!# endif
!
!
!! begin
!
!! only for first iteration
!  if (step_iterations_nb == 1) then
!
!!   test association
!    call object_associated_or_die_by_index (object_ind)
!    call object_associated_or_die_by_index (object_av_ind)
!
!!   test on type
!    object_type = objects(object_ind)%type
!    object_av_type = objects(object_av_ind)%type
!
!    if ( object_type /= object_av_type) then
!     write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is ', object_type
!     write(6,*) trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is ', object_av_type
!     write(6,*) trim(lhere),': they should be identical'
!     call die (lhere)
!    endif
!
!!   test on dimensions
!    dim1 = objects(object_ind)%dimensions(1)
!    dim_av1 = objects(object_av_ind)%dimensions(1)
!    dim2 = objects(object_ind)%dimensions(2)
!    dim_av2 = objects(object_av_ind)%dimensions(2)
!
!    if ( dim1 /= dim_av1 .or. dim2 /= dim_av2 ) then
!     write(6,*) trim(lhere),': dimensions of object', trim(objects(object_ind)%name),' are ', dim1, dim2
!     write(6,*) trim(lhere),': dimensions of object ',trim(objects(object_av_ind)%name),' are ', dim_av1, dim_av2
!     write(6,*) trim(lhere),': they should be identical'
!     call die (lhere)
!    endif
!
!!  initialization
!   call alloc ('objects(object_ind)%sum_double_2', objects(object_ind)%sum_double_2, dim1, dim2)
!
!! for all-electron move version
!!   if (index (trim(mode), 'mov1') == 0 ) then
!!       call alloc ('objects(object_ind)%previous_double_2', objects(object_ind)%previous_double_2, dim1, dim2)
!!   endif
!
!  endif ! first iteration
!
!! invalidate by default the average
!  call object_invalidate_by_index (object_av_ind)
!
!! reinitialization of sum for each block
!  if (mod(step_iterations_nb - 1, nstep) == 0) then
!     objects(object_ind)%sum_double_2 = 0.d0
!  endif
!
!! at each step
!
!! for one-electron move version
!!  if (index (trim(mode), 'mov1') /= 0) then
!  objects(object_ind)%sum_double_2 = objects(object_ind)%sum_double_2 + objects(object_ind)%pointer_double_2
!
!! for all-electron move version
!!   else
!!   objects(object_ind)%sum_double_2 = objects(object_ind)%sum_double_2 + prob_acc*objects(object_ind)%pointer_double_2 + prob_rej*objects(object_ind)%previous_double_2
!!   objects(object_ind)%previous_double_2 = objects(object_ind)%pointer_double_2
!!
!!  endif
!
!!  at the end of each block
!   if (mod(step_iterations_nb , nstep) == 0) then
!
!!    initialization for first block
!     if (block_iterations_nb == 1 ) then
!       call alloc ('objects(object_ind)%sum_blk_double_2', objects(object_ind)%sum_blk_double_2,  objects(object_ind)%dimensions(1), objects(object_ind)%dimensions(2))
!       objects(object_ind)%sum_blk_double_2 = 0.d0
!     endif
!
!# if defined (MPI)
!!    sum values from all processes
!     call alloc ('collect', collect,  objects(object_ind)%dimensions(1),  objects(object_ind)%dimensions(2))
!     call mpi_allreduce(objects(object_ind)%sum_double_2,collect,objects(object_ind)%dimensions(1)*objects(object_ind)%dimensions(2),mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
!     if (ierr /= 0) then
!        call die (lhere, 'error in mpi_allreduce')
!     endif
!     objects(object_ind)%sum_double_2 (:,:) = collect (:,:)
!# endif
!
!     objects(object_ind)%sum_blk_double_2 = objects(object_ind)%sum_blk_double_2 + objects(object_ind)%sum_double_2/nstep_total
!
!!    calculate average
!     objects(object_av_ind)%pointer_double_2 = objects(object_ind)%sum_blk_double_2 / block_iterations_nb
!     call object_modified_by_index (object_av_ind)
!
!   endif
!
! end subroutine object_average_by_index_double_2

! ===================================================================================
  subroutine object_block_average_by_index_double_0 (object_ind, object_bav_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate block average of object
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_ind, object_bav_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_block_average_by_index_double_0'
  character(len=max_string_len_type)   :: object_type, object_bav_type

# if defined (MPI)
  integer ierr
  real(dp) collect
# endif

! begin

! only for first iteration
  if (step_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_ind)
    call object_associated_or_die_by_index (object_bav_ind)

!   test on type
    object_type = objects(object_ind)%type
    object_bav_type = objects(object_bav_ind)%type

    if (object_type /= object_bav_type) then
     write(6,'(5a)') trim(lhere),': type of object ',trim(objects(object_ind)%name),' is ', object_type
     write(6,'(5a)') trim(lhere),': type of object ',trim(objects(object_bav_ind)%name),' is ', object_bav_type
     write(6,'(2a)') trim(lhere),': they should be identical'
     call die (lhere)
    endif

  endif ! first iteration

! invalidate by default the block average
  call object_invalidate_by_index (object_bav_ind)

! reinitialization of sum for each block
  if (mod(step_iterations_nb - 1, nstep) == 0) then
     objects(object_ind)%sum_double_0 = 0.d0
  endif

! at each step:

! for one-electron move version
!  if (index (trim(mode), 'mov1') /= 0) then
  objects(object_ind)%sum_double_0 = objects(object_ind)%sum_double_0 + objects(object_ind)%pointer_double_0

! for all-electron move version
!   else
!   objects(object_ind)%sum_double_0 = objects(object_ind)%sum_double_0 + prob_acc*objects(object_ind)%pointer_double_0 + prob_rej*objects(object_ind)%previous_double_0
!   write(6,*) 'objects(object_ind)%sum_double_0=',objects(object_ind)%sum_double_0
!   objects(object_ind)%previous_double_0 = objects(object_ind)%pointer_double_0

!  endif

!  at the end of each block
   if (mod(step_iterations_nb, nstep) == 0) then

# if defined (MPI)
!    sum values from all processes
     call mpi_allreduce(objects(object_ind)%sum_double_0,collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
     if (ierr /= 0) then
        call die (lhere, 'error in mpi_allreduce')
     endif
     objects(object_ind)%sum_double_0 = collect
# endif

!    calculate block average
     objects(object_bav_ind)%pointer_double_0  = objects(object_ind)%sum_double_0 / nstep_total
     call object_modified_by_index (object_bav_ind)

   endif

 end subroutine object_block_average_by_index_double_0

! ===================================================================================
  subroutine object_block_average_by_index_double_1 (object_ind, object_bav_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate block average of object
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_ind, object_bav_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_block_average_by_index_double_1'
  character(len=max_string_len_type)   :: object_type, object_bav_type
  integer dim1, dim_bav1


# if defined (MPI)
  integer ierr
  real(dp), allocatable :: collect(:)
# endif



! begin

! only for first iteration
  if (step_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_ind)
    call object_associated_or_die_by_index (object_bav_ind)

!   test on type
    object_type = objects(object_ind)%type
    object_bav_type = objects(object_bav_ind)%type

    if (object_type /= object_bav_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is ', object_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_bav_ind)%name),' is ', object_bav_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   test on dimensions
    dim1 = objects(object_ind)%dimensions(1)
    dim_bav1 = objects(object_bav_ind)%dimensions(1)

    if (dim1 /= dim_bav1) then
     write(6,'(4a,i8)') trim(lhere),': dimension of object', trim(objects(object_ind)%name),' is ', dim1
     write(6,'(4a,i8)') trim(lhere),': dimension of object ',trim(objects(object_bav_ind)%name),' is ', dim_bav1
     write(6,'(2a)') trim(lhere),': they should be identical'
     call die (lhere)
    endif

!  initialization
   call alloc ('objects(object_ind)%sum_double_1', objects(object_ind)%sum_double_1, dim1)

! for all-electron move version
!   if (index (trim(mode), 'mov1') == 0 ) then
!       call alloc ('objects(object_ind)%previous_double_1', objects(object_ind)%previous_double_1, dim1)
!   endif

  endif ! first iteration

! invalidate by default the block average
  call object_invalidate_by_index (object_bav_ind)

! reinitialization of sum for each block
  if (mod(step_iterations_nb - 1, nstep) == 0) then
     objects(object_ind)%sum_double_1 = 0.d0
  endif

! at each step

! for one-electron move version
!  if (index (trim(mode), 'mov1') /= 0) then
  objects(object_ind)%sum_double_1 = objects(object_ind)%sum_double_1 + objects(object_ind)%pointer_double_1

! for all-electron move version
!   else
!   objects(object_ind)%sum_double_1 = objects(object_ind)%sum_double_1 + prob_acc*objects(object_ind)%pointer_double_1 + prob_rej*objects(object_ind)%previous_double_1
!   objects(object_ind)%previous_double_1 = objects(object_ind)%pointer_double_1
!
!  endif

!  at the end of each block
   if (mod(step_iterations_nb, nstep) == 0) then

# if defined (MPI)
!    sum values from all processes
     call alloc ('collect', collect, objects(object_ind)%dimensions(1))
     call mpi_allreduce(objects(object_ind)%sum_double_1,collect,objects(object_ind)%dimensions(1),mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
     if (ierr /= 0) then
        call die (lhere, 'error in mpi_allreduce')
     endif
     objects(object_ind)%sum_double_1 (:) = collect (:)
# endif

!    calculate block average
     objects(object_bav_ind)%pointer_double_1 = objects(object_ind)%sum_double_1 / nstep_total
     call object_modified_by_index (object_bav_ind)

   endif

 end subroutine object_block_average_by_index_double_1

! ===================================================================================
  subroutine object_block_average_by_index_double_2 (object_ind, object_bav_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate block average of object
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_ind, object_bav_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_block_average_by_index_double_2'
  character(len=max_string_len_type)   :: object_type, object_bav_type
  integer dim1, dim_bav1, dim2, dim_bav2

# if defined (MPI)
  integer ierr
  real(dp), allocatable :: collect(:,:)
# endif

! begin

! only for first iteration
  if (step_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_ind)
    call object_associated_or_die_by_index (object_bav_ind)

!   test on type
    object_type = objects(object_ind)%type
    object_bav_type = objects(object_bav_ind)%type

    if (object_type /= object_bav_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is ', object_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_bav_ind)%name),' is ', object_bav_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   test on dimensions
    dim1 = objects(object_ind)%dimensions(1)
    dim_bav1 = objects(object_bav_ind)%dimensions(1)
    dim2 = objects(object_ind)%dimensions(2)
    dim_bav2 = objects(object_bav_ind)%dimensions(2)

    if (dim1 /= dim_bav1 .or. dim2 /= dim_bav2) then
     write(6,*) trim(lhere),': dimensions of object', trim(objects(object_ind)%name),' are ', dim1, dim2
     write(6,*) trim(lhere),': dimensions of object ',trim(objects(object_bav_ind)%name),' are ', dim_bav1, dim_bav2
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!  initialization
   call alloc ('objects(object_ind)%sum_double_2', objects(object_ind)%sum_double_2, dim1, dim2)

! for all-electron move version
!   if (index (trim(mode), 'mov1') == 0 ) then
!       call alloc ('objects(object_ind)%previous_double_2', objects(object_ind)%previous_double_2, dim1, dim2)
!   endif

  endif ! first iteration

! invalidate by default the average
  call object_invalidate_by_index (object_bav_ind)

! reinitialization of sum for each block
  if (mod(step_iterations_nb - 1, nstep) == 0) then
     objects(object_ind)%sum_double_2 = 0.d0
  endif

! at each step

! for one-electron move version
!  if (index (trim(mode), 'mov1') /= 0) then
  objects(object_ind)%sum_double_2 = objects(object_ind)%sum_double_2 + objects(object_ind)%pointer_double_2

! for all-electron move version
!   else
!   objects(object_ind)%sum_double_2 = objects(object_ind)%sum_double_2 + prob_acc*objects(object_ind)%pointer_double_2 + prob_rej*objects(object_ind)%previous_double_2
!   objects(object_ind)%previous_double_2 = objects(object_ind)%pointer_double_2
!
!  endif

!  at the end of each block
   if (mod(step_iterations_nb, nstep) == 0) then

# if defined (MPI)
!    sum values from all processes
     call alloc ('collect', collect,  objects(object_ind)%dimensions(1),  objects(object_ind)%dimensions(2))
     call mpi_allreduce(objects(object_ind)%sum_double_2,collect,objects(object_ind)%dimensions(1)*objects(object_ind)%dimensions(2),mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
     if (ierr /= 0) then
        call die (lhere, 'error in mpi_allreduce')
     endif
     objects(object_ind)%sum_double_2 (:,:) = collect (:,:)
# endif

!    calculate average
     objects(object_bav_ind)%pointer_double_2 = objects(object_ind)%sum_double_2 / nstep_total
     call object_modified_by_index (object_bav_ind)

   endif

 end subroutine object_block_average_by_index_double_2

! ===================================================================================
  subroutine object_global_average_by_index_double_0 (object_bav_ind, object_av_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate global average of object over all blocks
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_bav_ind, object_av_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_global_average_by_index_double_0'
  character(len=max_string_len_type)   :: object_bav_type, object_av_type

! begin

! only for first block
  if (block_iterations_nb == 1 ) then

!   test association
    call object_associated_or_die_by_index (object_bav_ind)
    call object_associated_or_die_by_index (object_av_ind)

!   test on type
    object_bav_type = objects(object_bav_ind)%type
    object_av_type = objects(object_av_ind)%type

    if (object_bav_type /= object_av_type) then
     write(6,'(5a)') trim(lhere),': type of object ',trim(objects(object_bav_ind)%name),' is ', object_bav_type
     write(6,'(5a)') trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is ', object_av_type
     write(6,'(2a)') trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   initialization
    objects(object_bav_ind)%sum_blk_double_0 = 0.d0

  endif ! first block

! sum over blocks
  objects(object_bav_ind)%sum_blk_double_0 = objects(object_bav_ind)%sum_blk_double_0 + objects(object_bav_ind)%pointer_double_0

! calculate global average
  objects(object_av_ind)%pointer_double_0  = objects(object_bav_ind)%sum_blk_double_0 / block_iterations_nb
  call object_modified_by_index (object_av_ind)

 end subroutine object_global_average_by_index_double_0

! ===================================================================================
  subroutine object_global_average_by_index_double_1 (object_bav_ind, object_av_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate global average of object over all blocks
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_bav_ind, object_av_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_global_average_by_index_double_1'
  character(len=max_string_len_type)   :: object_bav_type, object_av_type
  integer dim_bav1, dim_av1

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_bav_ind)
    call object_associated_or_die_by_index (object_av_ind)

!   test on type
    object_bav_type = objects(object_bav_ind)%type
    object_av_type = objects(object_av_ind)%type

    if (object_bav_type /= object_av_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_bav_ind)%name),' is ', object_bav_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is ', object_av_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   test on dimensions
    dim_bav1 = objects(object_bav_ind)%dimensions(1)
    dim_av1 = objects(object_av_ind)%dimensions(1)

    if (dim_bav1 /= dim_av1) then
     write(6,'(4a,i8)') trim(lhere),': dimension of object', trim(objects(object_bav_ind)%name),' is ', dim_bav1
     write(6,'(4a,i8)') trim(lhere),': dimension of object ',trim(objects(object_av_ind)%name),' is ', dim_av1
     write(6,'(2a)') trim(lhere),': they should be identical'
     call die (lhere)
    endif

!  initialization
   call alloc ('objects(object_bav_ind)%sum_blk_double_1', objects(object_bav_ind)%sum_blk_double_1, objects(object_bav_ind)%dimensions(1))
   objects(object_bav_ind)%sum_blk_double_1 = 0.d0

  endif ! first block

! sum over blocks
  objects(object_bav_ind)%sum_blk_double_1 = objects(object_bav_ind)%sum_blk_double_1 + objects(object_bav_ind)%pointer_double_1

! calculate global average
  objects(object_av_ind)%pointer_double_1 = objects(object_bav_ind)%sum_blk_double_1 / block_iterations_nb
  call object_modified_by_index (object_av_ind)

 end subroutine object_global_average_by_index_double_1

! ===================================================================================
  subroutine object_global_average_by_index_double_2 (object_bav_ind, object_av_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate global average of object over all blocks
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) ::  object_bav_ind, object_av_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_global_average_by_index_double_2'
  character(len=max_string_len_type)   :: object_bav_type, object_av_type
  integer dim_bav1, dim_av1, dim_bav2, dim_av2

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_bav_ind)
    call object_associated_or_die_by_index (object_av_ind)

!   test on type
    object_bav_type = objects(object_bav_ind)%type
    object_av_type = objects(object_av_ind)%type

    if (object_bav_type /= object_av_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_bav_ind)%name),' is ', object_bav_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is ', object_av_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   test on dimensions
    dim_bav1 = objects(object_bav_ind)%dimensions(1)
    dim_av1 = objects(object_av_ind)%dimensions(1)
    dim_bav2 = objects(object_bav_ind)%dimensions(2)
    dim_av2 = objects(object_av_ind)%dimensions(2)

    if (dim_bav1 /= dim_av1 .or. dim_bav2 /= dim_av2) then
     write(6,*) trim(lhere),': dimensions of object', trim(objects(object_bav_ind)%name),' are ', dim_bav1, dim_bav2
     write(6,*) trim(lhere),': dimensions of object ',trim(objects(object_av_ind)%name),' are ', dim_av1, dim_av2
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   initialization
    call alloc ('objects(object_bav_ind)%sum_blk_double_2', objects(object_bav_ind)%sum_blk_double_2,  objects(object_bav_ind)%dimensions(1), objects(object_bav_ind)%dimensions(2))
    objects(object_bav_ind)%sum_blk_double_2 = 0.d0

  endif ! first block

! sum over blocks
  objects(object_bav_ind)%sum_blk_double_2 = objects(object_bav_ind)%sum_blk_double_2 + objects(object_bav_ind)%pointer_double_2

! calculate global average
  objects(object_av_ind)%pointer_double_2 = objects(object_bav_ind)%sum_blk_double_2 / block_iterations_nb
  call object_modified_by_index (object_av_ind)

 end subroutine object_global_average_by_index_double_2

! ===================================================================================
  subroutine object_variance_by_index_double_0 (object_av_ind, object_var_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate variance of an average object over MC iterations
!
! Created       : J. Toulouse, 05 Sep 2007
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_av_ind, object_var_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_variance_by_index_double_0'
  character(len=max_string_len_type) object_av_type, object_var_type

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_av_ind)
    call object_associated_or_die_by_index (object_var_ind)

!   test on type
    object_av_type = objects(object_av_ind)%type
    object_var_type = objects(object_var_ind)%type

    if (object_av_type /= object_var_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is ', object_av_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_var_ind)%name),' is ', object_var_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   intermediate objects
    objects(object_var_ind)%sum_blk_square_double_0 = 0.d0
    objects(object_var_ind)%previous_av1_double_0 = 0.d0

   endif ! first block

!  calculate sum of square of averages over blocks
   objects(object_var_ind)%sum_blk_square_double_0 = objects(object_var_ind)%sum_blk_square_double_0   &
       + (objects(object_av_ind)%pointer_double_0 * block_iterations_nb - objects(object_var_ind)%previous_av1_double_0 * (block_iterations_nb - 1))**2

!  calculate variance
   if (block_iterations_nb == 1) then
    objects(object_var_ind)%pointer_double_0 = 0.d0
   else
    objects(object_var_ind)%pointer_double_0 = (objects(object_var_ind)%sum_blk_square_double_0/block_iterations_nb - objects(object_av_ind)%pointer_double_0**2)/(block_iterations_nb-1)
   endif
   call object_modified_by_index (object_var_ind)

!  save current average value for next iteration
   objects(object_var_ind)%previous_av1_double_0 = objects(object_av_ind)%pointer_double_0

 end subroutine object_variance_by_index_double_0

! ===================================================================================
  subroutine object_variance_by_index_double_1 (object_av_ind, object_var_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate variance of an average object over MC iterations
!
! Created       : J. Toulouse, 07 Sep 2007
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_av_ind, object_var_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_variance_by_index_double_1'
  character(len=max_string_len_type) object_av_type, object_var_type
  integer dim_av1, dim_var1

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_av_ind)
    call object_associated_or_die_by_index (object_var_ind)

!   test on type
    object_av_type = objects(object_av_ind)%type
    object_var_type = objects(object_var_ind)%type

    if (object_av_type /= object_var_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is ', object_av_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_var_ind)%name),' is ', object_var_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   test on dimensions
    dim_av1 = objects(object_av_ind)%dimensions(1)
    dim_var1 = objects(object_var_ind)%dimensions(1)

    if (dim_av1 /= dim_var1) then
     write(6,*) trim(lhere),': dimension of object', trim(objects(object_av_ind)%name),' is ', dim_av1
     write(6,*) trim(lhere),': dimension of object ',trim(objects(object_var_ind)%name),' is ', dim_var1
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   intermediate objects
    call alloc ('objects(object_var_ind)%sum_blk_square_double_1', objects(object_var_ind)%sum_blk_square_double_1, dim_var1)
    call alloc('objects(object_var_ind)%previous_av1_double_1', objects(object_var_ind)%previous_av1_double_1, dim_av1)
    objects(object_var_ind)%sum_blk_square_double_1 = 0.d0
    objects(object_var_ind)%previous_av1_double_1 = 0.d0

   endif ! first block


!  calculate sum of square of averages over blocks
   objects(object_var_ind)%sum_blk_square_double_1 = objects(object_var_ind)%sum_blk_square_double_1   &
       + (objects(object_av_ind)%pointer_double_1 * block_iterations_nb - objects(object_var_ind)%previous_av1_double_1 * (block_iterations_nb - 1 ))**2

!  calculate variance
   if (block_iterations_nb == 1) then
    objects(object_var_ind)%pointer_double_1 = 0.d0
   else
    objects(object_var_ind)%pointer_double_1 = (objects(object_var_ind)%sum_blk_square_double_1/block_iterations_nb - objects(object_av_ind)%pointer_double_1**2)/(block_iterations_nb-1)
   endif
   call object_modified_by_index (object_var_ind)

!  save current average value for next iteration
   objects(object_var_ind)%previous_av1_double_1 = objects(object_av_ind)%pointer_double_1

 end subroutine object_variance_by_index_double_1

! ===================================================================================
  subroutine object_variance_by_index_double_2 (object_av_ind, object_var_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate variance of an average object over MC iterations
!
! Created       : J. Toulouse, 07 Sep 2007
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_av_ind, object_var_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_variance_by_index_double_2'
  character(len=max_string_len_type) object_av_type, object_var_type
  integer dim_av1, dim_var1, dim_av2, dim_var2

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_av_ind)
    call object_associated_or_die_by_index (object_var_ind)

!   test on type
    object_av_type = objects(object_av_ind)%type
    object_var_type = objects(object_var_ind)%type

    if (object_av_type /= object_var_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is ', object_av_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_var_ind)%name),' is ', object_var_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   test on dimensions
    dim_av1 = objects(object_av_ind)%dimensions(1)
    dim_var1 = objects(object_var_ind)%dimensions(1)
    dim_av2 = objects(object_av_ind)%dimensions(2)
    dim_var2 = objects(object_var_ind)%dimensions(2)

    if (dim_av1 /= dim_var1 .or. dim_av2 /= dim_var2) then
     write(6,*) trim(lhere),': dimensions of object', trim(objects(object_av_ind)%name),' are ', dim_av1, dim_av2
     write(6,*) trim(lhere),': dimensions of object ',trim(objects(object_var_ind)%name),' are ', dim_var1, dim_var2
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   intermediate objects
    call alloc ('objects(object_var_ind)%sum_blk_square_double_2', objects(object_var_ind)%sum_blk_square_double_2, dim_var1, dim_var2)
    call alloc('objects(object_var_ind)%previous_av1_double_2', objects(object_var_ind)%previous_av1_double_2, dim_av1, dim_av2)
    objects(object_var_ind)%sum_blk_square_double_2 = 0.d0
    objects(object_var_ind)%previous_av1_double_2 = 0.d0

   endif ! first block

!  calculate sum of square of averages over blocks
   objects(object_var_ind)%sum_blk_square_double_2 = objects(object_var_ind)%sum_blk_square_double_2   &
       + ( objects(object_av_ind)%pointer_double_2 * block_iterations_nb - objects(object_var_ind)%previous_av1_double_2 * (block_iterations_nb - 1 ) )**2

!  calculate variance
   if (block_iterations_nb == 1) then
    objects(object_var_ind)%pointer_double_2 = 0.d0
   else
    objects(object_var_ind)%pointer_double_2 = (objects(object_var_ind)%sum_blk_square_double_2/block_iterations_nb - objects(object_av_ind)%pointer_double_2**2)/(block_iterations_nb-1)
   endif
   call object_modified_by_index (object_var_ind)

!  save current average value for next iteration
   objects(object_var_ind)%previous_av1_double_2 = objects(object_av_ind)%pointer_double_2

 end subroutine object_variance_by_index_double_2

! ===================================================================================
  subroutine object_covariance_by_index_double_0_double_0 (object_av1_ind, object_av2_ind, object_covar_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate covariance of two averaged objects
!
! Created       : J. Toulouse, 17 Sep 2007
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_av1_ind, object_av2_ind, object_covar_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_covariance_by_index_double_0_double_0'
  character(len=max_string_len_type) object_covar_type

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_av1_ind)
    call object_associated_or_die_by_index (object_av2_ind)
    call object_associated_or_die_by_index (object_covar_ind)

!   test on type
    object_covar_type = objects(object_covar_ind)%type
    if (object_covar_type /= 'double_0') then
     write(6,'(6a)') trim(lhere),': type of object >',trim(objects(object_covar_ind)%name),'< is ', trim(object_covar_type),' /= double_0'
     call die (lhere)
    endif

!   intermediate objects
    objects(object_covar_ind)%sum_blk_square_double_0 = 0.d0
    objects(object_covar_ind)%previous_av1_double_0 = 0.d0
    objects(object_covar_ind)%previous_av2_double_0 = 0.d0

   endif ! first block

!  calculate sum of product of averages over blocks
   objects(object_covar_ind)%sum_blk_square_double_0 = objects(object_covar_ind)%sum_blk_square_double_0   &
       + (objects(object_av1_ind)%pointer_double_0 * block_iterations_nb - objects(object_covar_ind)%previous_av1_double_0 * (block_iterations_nb - 1)) &
       * (objects(object_av2_ind)%pointer_double_0 * block_iterations_nb - objects(object_covar_ind)%previous_av2_double_0 * (block_iterations_nb - 1))

!  calculate covariance
   if (block_iterations_nb == 1) then
    objects(object_covar_ind)%pointer_double_0 = 0.d0
   else
    objects(object_covar_ind)%pointer_double_0 = (objects(object_covar_ind)%sum_blk_square_double_0/block_iterations_nb - objects(object_av1_ind)%pointer_double_0*objects(object_av2_ind)%pointer_double_0)/(block_iterations_nb-1)
   endif
   call object_modified_by_index (object_covar_ind)

!  save current average value for next iteration
   objects(object_covar_ind)%previous_av1_double_0 = objects(object_av1_ind)%pointer_double_0
   objects(object_covar_ind)%previous_av2_double_0 = objects(object_av2_ind)%pointer_double_0

 end subroutine object_covariance_by_index_double_0_double_0

! ===================================================================================
  subroutine object_covariance_by_index_double_1_double_0 (object_av1_ind, object_av2_ind, object_covar_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate covariance of two averaged objects
!
! Created       : J. Toulouse, 19 Sep 2007
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_av1_ind, object_av2_ind, object_covar_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_covariance_by_index_double_1_double_0'
  character(len=max_string_len_type) object_covar_type
  integer dim_av11, dim_covar1

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_av1_ind)
    call object_associated_or_die_by_index (object_av2_ind)
    call object_associated_or_die_by_index (object_covar_ind)

!   test on type
    object_covar_type = objects(object_covar_ind)%type

    if (object_covar_type /= 'double_1') then
     write(6,'(6a)') trim(lhere),': type of object >',trim(objects(object_covar_ind)%name),'< is ', object_covar_type,' /= double_1'
     call die (lhere)
    endif

!   test on dimensions
    dim_av11 = objects(object_av1_ind)%dimensions(1)
    dim_covar1 = objects(object_covar_ind)%dimensions(1)

    if (dim_av11 /= dim_covar1) then
     write(6,*) trim(lhere),': dimension of object >',trim(objects(object_av1_ind)%name),'< is ', dim_av11
     write(6,*) trim(lhere),': dimension of object >',trim(objects(object_covar_ind)%name),'< is ', dim_covar1
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   intermediate objects
    call alloc ('objects(object_covar_ind)%sum_blk_square_double_1', objects(object_covar_ind)%sum_blk_square_double_1, dim_covar1)
    call alloc ('objects(object_covar_ind)%previous_av1_double_1', objects(object_covar_ind)%previous_av1_double_1, dim_av11)
    objects(object_covar_ind)%sum_blk_square_double_1 = 0.d0
    objects(object_covar_ind)%previous_av1_double_1 = 0.d0
    objects(object_covar_ind)%previous_av2_double_0 = 0.d0

   endif ! first block

!  calculate sum of product of averages over blocks
   objects(object_covar_ind)%sum_blk_square_double_1 = objects(object_covar_ind)%sum_blk_square_double_1   &
       + (objects(object_av1_ind)%pointer_double_1 * block_iterations_nb - objects(object_covar_ind)%previous_av1_double_1 * (block_iterations_nb - 1 )) &
       * (objects(object_av2_ind)%pointer_double_0 * block_iterations_nb - objects(object_covar_ind)%previous_av2_double_0 * (block_iterations_nb - 1 ))

!  calculate covariance
   if (block_iterations_nb == 1) then
    objects(object_covar_ind)%pointer_double_1 = 0.d0
   else
    objects(object_covar_ind)%pointer_double_1 = (objects(object_covar_ind)%sum_blk_square_double_1/block_iterations_nb - objects(object_av1_ind)%pointer_double_1*objects(object_av2_ind)%pointer_double_0)/(block_iterations_nb-1)
   endif
   call object_modified_by_index (object_covar_ind)

!  save current average value for next iteration
   objects(object_covar_ind)%previous_av1_double_1 = objects(object_av1_ind)%pointer_double_1
   objects(object_covar_ind)%previous_av2_double_0 = objects(object_av2_ind)%pointer_double_0

 end subroutine object_covariance_by_index_double_1_double_0

! ===================================================================================
  subroutine object_covariance_by_index_double_0_double_1 (object_av1_ind, object_av2_ind, object_covar_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate covariance of two averaged objects
!
! Created       : J. Toulouse, 19 Sep 2007
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_av1_ind, object_av2_ind, object_covar_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_covariance_by_index_double_0_double_1'
  character(len=max_string_len_type) object_covar_type
  integer dim_av21, dim_covar1

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_av1_ind)
    call object_associated_or_die_by_index (object_av2_ind)
    call object_associated_or_die_by_index (object_covar_ind)

!   test on type
    object_covar_type = objects(object_covar_ind)%type

    if (object_covar_type /= 'double_1') then
     write(6,'(6a)') trim(lhere),': type of object >',trim(objects(object_covar_ind)%name),'< is ', object_covar_type,' /= double_1'
     call die (lhere)
    endif

!   test on dimensions
    dim_av21 = objects(object_av2_ind)%dimensions(1)
    dim_covar1 = objects(object_covar_ind)%dimensions(1)

    if (dim_av21 /= dim_covar1) then
     write(6,*) trim(lhere),': dimension of object >',trim(objects(object_av2_ind)%name),'< is ', dim_av21
     write(6,*) trim(lhere),': dimension of object >',trim(objects(object_covar_ind)%name),'< is ', dim_covar1
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   intermediate objects
    call alloc ('objects(object_covar_ind)%sum_blk_square_double_1', objects(object_covar_ind)%sum_blk_square_double_1, dim_covar1)
    call alloc ('objects(object_covar_ind)%previous_av2_double_1', objects(object_covar_ind)%previous_av2_double_1, dim_av21)
    objects(object_covar_ind)%sum_blk_square_double_1 = 0.d0
    objects(object_covar_ind)%previous_av1_double_0 = 0.d0
    objects(object_covar_ind)%previous_av2_double_1 = 0.d0

   endif ! first block

!  calculate sum of product of averages over blocks
   objects(object_covar_ind)%sum_blk_square_double_1 = objects(object_covar_ind)%sum_blk_square_double_1   &
       + (objects(object_av1_ind)%pointer_double_0 * block_iterations_nb - objects(object_covar_ind)%previous_av1_double_0 * (block_iterations_nb - 1 )) &
       * (objects(object_av2_ind)%pointer_double_1 * block_iterations_nb - objects(object_covar_ind)%previous_av2_double_1 * (block_iterations_nb - 1 ))

!  calculate covariance
   if (block_iterations_nb == 1) then
    objects(object_covar_ind)%pointer_double_1 = 0.d0
   else
    objects(object_covar_ind)%pointer_double_1 = (objects(object_covar_ind)%sum_blk_square_double_1/block_iterations_nb - objects(object_av1_ind)%pointer_double_0*objects(object_av2_ind)%pointer_double_1)/(block_iterations_nb-1)
   endif
   call object_modified_by_index (object_covar_ind)

!  save current average value for next iteration
   objects(object_covar_ind)%previous_av1_double_0 = objects(object_av1_ind)%pointer_double_0
   objects(object_covar_ind)%previous_av2_double_1 = objects(object_av2_ind)%pointer_double_1

 end subroutine object_covariance_by_index_double_0_double_1

! ===================================================================================
  subroutine object_covariance_by_index_double_1_double_1_double_1 (object_av1_ind, object_av2_ind, object_covar_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate covariance of two averaged objects of corresponding dimmensions
!
! Created       : J. Toulouse, 19 Sep 2007
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_av1_ind, object_av2_ind, object_covar_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_covariance_by_index_double_1_double_1_double_1'
  character(len=max_string_len_type) object_covar_type
  integer dim_av11, dim_av21, dim_covar1

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_av1_ind)
    call object_associated_or_die_by_index (object_av2_ind)
    call object_associated_or_die_by_index (object_covar_ind)

!   test on type
    object_covar_type = objects(object_covar_ind)%type

    if (object_covar_type /= 'double_1') then
     write(6,'(6a)') trim(lhere),': type of object >',trim(objects(object_covar_ind)%name),'< is ', object_covar_type,' /= double_1'
     call die (lhere)
    endif

!   test on dimensions
    dim_av11 = objects(object_av1_ind)%dimensions(1)
    dim_av21 = objects(object_av2_ind)%dimensions(1)
    dim_covar1 = objects(object_covar_ind)%dimensions(1)

    if (dim_av11 /= dim_av21 .or. dim_av21 /= dim_covar1) then
     write(6,*) trim(lhere),': dimension of object >',trim(objects(object_av1_ind)%name),'< is ', dim_av11
     write(6,*) trim(lhere),': dimension of object >',trim(objects(object_av2_ind)%name),'< is ', dim_av21
     write(6,*) trim(lhere),': dimensions of object >',trim(objects(object_covar_ind)%name),'< is ', dim_covar1
     write(6,*) trim(lhere),': they should match'
     call die (lhere)
    endif

!   intermediate objects
    call alloc ('objects(object_covar_ind)%sum_blk_square_double_1', objects(object_covar_ind)%sum_blk_square_double_1, dim_covar1)
    call alloc ('objects(object_covar_ind)%previous_av1_double_1', objects(object_covar_ind)%previous_av1_double_1, dim_av11)
    call alloc ('objects(object_covar_ind)%previous_av2_double_1', objects(object_covar_ind)%previous_av2_double_1, dim_av21)
    objects(object_covar_ind)%sum_blk_square_double_1 = 0.d0
    objects(object_covar_ind)%previous_av1_double_1 = 0.d0
    objects(object_covar_ind)%previous_av2_double_1 = 0.d0

   endif ! first block

!  calculate sum of product of averages over blocks
   objects(object_covar_ind)%sum_blk_square_double_1 (:) = objects(object_covar_ind)%sum_blk_square_double_1 (:) &
       + (objects(object_av1_ind)%pointer_double_1 (:) * block_iterations_nb - objects(object_covar_ind)%previous_av1_double_1 (:) * (block_iterations_nb - 1 )) &
       * (objects(object_av2_ind)%pointer_double_1 (:) * block_iterations_nb - objects(object_covar_ind)%previous_av2_double_1 (:) * (block_iterations_nb - 1 ))

!  calculate covariance
   if (block_iterations_nb == 1) then
    objects(object_covar_ind)%pointer_double_1 (:) = 0.d0
   else
    objects(object_covar_ind)%pointer_double_1 (:) = (objects(object_covar_ind)%sum_blk_square_double_1 (:)/block_iterations_nb - objects(object_av1_ind)%pointer_double_1 (:) * objects(object_av2_ind)%pointer_double_1 (:))/(block_iterations_nb-1)
   endif
   call object_modified_by_index (object_covar_ind)

!  save current average value for next iteration
   objects(object_covar_ind)%previous_av1_double_1 (:) = objects(object_av1_ind)%pointer_double_1 (:)
   objects(object_covar_ind)%previous_av2_double_1 (:) = objects(object_av2_ind)%pointer_double_1 (:)

 end subroutine object_covariance_by_index_double_1_double_1_double_1

! ===================================================================================
  subroutine object_covariance_by_index_double_1_double_1_double_2 (object_av1_ind, object_av2_ind, object_covar_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate covariance of two averaged objects of corresponding dimmensions
!
! Created       : J. Toulouse, 15 Feb 2008
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_av1_ind, object_av2_ind, object_covar_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_covariance_by_index_double_1_double_1_double_2'
  character(len=max_string_len_type) object_covar_type
  integer dim_av11, dim_av21, dim_covar1, dim_covar2
  integer i,j

! begin
  dim_covar1 = objects(object_covar_ind)%dimensions(1)
  dim_covar2 = objects(object_covar_ind)%dimensions(2)

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_av1_ind)
    call object_associated_or_die_by_index (object_av2_ind)
    call object_associated_or_die_by_index (object_covar_ind)

!   test on type
    object_covar_type = objects(object_covar_ind)%type

    if (object_covar_type /= 'double_2') then
     write(6,'(6a)') trim(lhere),': type of object >',trim(objects(object_covar_ind)%name),'< is ', object_covar_type,' /= double_2'
     call die (lhere)
    endif

!   test on dimensions
    dim_av11 = objects(object_av1_ind)%dimensions(1)
    dim_av21 = objects(object_av2_ind)%dimensions(1)

    if (dim_av11 /= dim_covar1 .or. dim_av21 /= dim_covar2) then
     write(6,*) trim(lhere),': dimension of object >',trim(objects(object_av1_ind)%name),'< is ', dim_av11
     write(6,*) trim(lhere),': dimension of object >',trim(objects(object_av2_ind)%name),'< is ', dim_av21
     write(6,*) trim(lhere),': dimensions of object >',trim(objects(object_covar_ind)%name),'< are ', dim_covar1, dim_covar2
     write(6,*) trim(lhere),': they should match'
     call die (lhere)
    endif

!   intermediate objects
    call alloc ('objects(object_covar_ind)%sum_blk_square_double_2', objects(object_covar_ind)%sum_blk_square_double_2, dim_covar1, dim_covar2)
    call alloc ('objects(object_covar_ind)%previous_av1_double_1', objects(object_covar_ind)%previous_av1_double_1, dim_av11)
    call alloc ('objects(object_covar_ind)%previous_av2_double_1', objects(object_covar_ind)%previous_av2_double_1, dim_av21)
    objects(object_covar_ind)%sum_blk_square_double_2 = 0.d0
    objects(object_covar_ind)%previous_av1_double_1 = 0.d0
    objects(object_covar_ind)%previous_av2_double_1 = 0.d0

   endif ! first block

!  calculate sum of product of averages over blocks
   do i = 1, dim_covar1
    do j = 1, dim_covar2
     objects(object_covar_ind)%sum_blk_square_double_2 (i,j) = objects(object_covar_ind)%sum_blk_square_double_2 (i,j) &
       + (objects(object_av1_ind)%pointer_double_1 (i) * block_iterations_nb - objects(object_covar_ind)%previous_av1_double_1 (i) * (block_iterations_nb - 1 )) &
       * (objects(object_av2_ind)%pointer_double_1 (j) * block_iterations_nb - objects(object_covar_ind)%previous_av2_double_1 (j) * (block_iterations_nb - 1 ))
    enddo
   enddo

!  calculate covariance
   if (block_iterations_nb == 1) then
    objects(object_covar_ind)%pointer_double_2 (:,:) = 0.d0
   else
    do i = 1, dim_covar1
     do j = 1, dim_covar2
      objects(object_covar_ind)%pointer_double_2 (i,j) = (objects(object_covar_ind)%sum_blk_square_double_2 (i,j)/block_iterations_nb - objects(object_av1_ind)%pointer_double_1 (i) * objects(object_av2_ind)%pointer_double_1 (j))/(block_iterations_nb-1)
     enddo
    enddo
   endif
   call object_modified_by_index (object_covar_ind)

!  save current average value for next iteration
   objects(object_covar_ind)%previous_av1_double_1 (:) = objects(object_av1_ind)%pointer_double_1 (:)
   objects(object_covar_ind)%previous_av2_double_1 (:) = objects(object_av2_ind)%pointer_double_1 (:)

 end subroutine object_covariance_by_index_double_1_double_1_double_2

! ===================================================================================
  subroutine object_error_by_index_double_0 (object_var_ind, object_err_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate statistical error corresponding to a variance
!
! Created       : J. Toulouse, 19 Oct 2005
! Modified      : J. Toulouse, 05 Feb 2006: use fortran pointers, warning association not checked!
! Modified      : J. Toulouse, 06 Sep 2007: statistical error calculated from variance
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_var_ind, object_err_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_error_by_index_double_0'
  character(len=max_string_len_type) object_var_type, object_err_type

! begin

! only for first block
  if (block_iterations_nb == 1) then
!
!   test association
    call object_associated_or_die_by_index (object_var_ind)
    call object_associated_or_die_by_index (object_err_ind)

!   test on type
    object_var_type = objects(object_var_ind)%type
    object_err_type = objects(object_err_ind)%type

    if (object_var_type /= object_err_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_var_ind)%name),' is ', object_var_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_err_ind)%name),' is ', object_err_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

   endif ! first block

!  calculate error
   objects(object_err_ind)%pointer_double_0 = dsqrt(objects(object_var_ind)%pointer_double_0)
   call object_modified_by_index (object_err_ind)

 end subroutine object_error_by_index_double_0

! ===================================================================================
  subroutine object_error_by_index_double_1 (object_var_ind, object_err_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate statistical error corresponding to a variance
!
! Created       : J. Toulouse, 19 Oct 2005
! Modified      : J. Toulouse, 05 Feb 2006: use fortran pointers, warning association not checked!
! Modified      : J. Toulouse, 07 Sep 2007: statistical error calculated from variance
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_var_ind, object_err_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_error_by_index_double_1'
  character(len=max_string_len_type) object_var_type, object_err_type
  integer dim_var1, dim_err1

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_var_ind)
    call object_associated_or_die_by_index (object_err_ind)

!   test on type
    object_var_type = objects(object_var_ind)%type
    object_err_type = objects(object_err_ind)%type

    if (object_var_type /= object_err_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_var_ind)%name),' is ', object_var_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_err_ind)%name),' is ', object_err_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   test on dimensions
    dim_var1 = objects(object_var_ind)%dimensions(1)
    dim_err1 = objects(object_err_ind)%dimensions(1)

    if (dim_var1 /= dim_err1) then
     write(6,*) trim(lhere),': dimension of object', trim(objects(object_var_ind)%name),' is ', dim_var1
     write(6,*) trim(lhere),': dimension of object ',trim(objects(object_err_ind)%name),' is ', dim_err1
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

   endif ! first block

!  calculate error
   objects(object_err_ind)%pointer_double_1 = dsqrt(objects(object_var_ind)%pointer_double_1)
   call object_modified_by_index (object_err_ind)

 end subroutine object_error_by_index_double_1

! ===================================================================================
  subroutine object_error_by_index_double_2 (object_var_ind, object_err_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate statistical error corresponding to a variance
!
! Created       : J. Toulouse, 19 Oct 2005
! Modified      : J. Toulouse, 05 Feb 2006: use fortran pointers, warning association not checked!
! Modified      : J. Toulouse, 07 Sep 2007: statistical error calculated from variance
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_var_ind, object_err_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_error_by_index_double_2'
  character(len=max_string_len_type) object_var_type, object_err_type
  integer dim_var1, dim_err1, dim_var2, dim_err2

! begin

! only for first block
  if (block_iterations_nb == 1) then

!   test association
    call object_associated_or_die_by_index (object_var_ind)
    call object_associated_or_die_by_index (object_err_ind)

!   test on type
    object_var_type = objects(object_var_ind)%type
    object_err_type = objects(object_err_ind)%type

    if (object_var_type /= object_err_type) then
     write(6,*) trim(lhere),': type of object ',trim(objects(object_var_ind)%name),' is ', object_var_type
     write(6,*) trim(lhere),': type of object ',trim(objects(object_err_ind)%name),' is ', object_err_type
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

!   test on dimensions
    dim_var1 = objects(object_var_ind)%dimensions(1)
    dim_err1 = objects(object_err_ind)%dimensions(1)
    dim_var2 = objects(object_var_ind)%dimensions(2)
    dim_err2 = objects(object_err_ind)%dimensions(2)

    if (dim_var1 /= dim_err1 .or. dim_var2 /= dim_err2) then
     write(6,*) trim(lhere),': dimensions of object', trim(objects(object_var_ind)%name),' are ', dim_var1, dim_var2
     write(6,*) trim(lhere),': dimensions of object ',trim(objects(object_err_ind)%name),' are ', dim_err1, dim_err2
     write(6,*) trim(lhere),': they should be identical'
     call die (lhere)
    endif

   endif ! first block

!  calculate error
   objects(object_err_ind)%pointer_double_2 = dsqrt(objects(object_var_ind)%pointer_double_2)
   call object_modified_by_index (object_err_ind)

 end subroutine object_error_by_index_double_2

! ===================================================================================
  subroutine object_average_walk_step_by_index_double_0 (object_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate average of object over MC iterations
!
! Created       : J. Toulouse, 12 Nov 2006
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_ind

! local
  integer walk_i

! begin

! only for first iteration
  if (step_iterations_nb == 1) then

!   initialization of sum
    objects(object_ind)%sum_double_0 = 0.d0

  endif ! first iteration

! at each step
  do walk_i = 1, nwalk
   objects(object_ind)%sum_double_0 = objects(object_ind)%sum_double_0 + objects(object_ind)%pointer_double_1(walk_i) * walker_weights(walk_i)
  enddo

 end subroutine object_average_walk_step_by_index_double_0

! ===================================================================================
  subroutine object_average_walk_step_by_index_double_1 (object_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate average of object over MC iterations
!
! Created       : J. Toulouse, 15 Nov 2006
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_average_walk_step_by_index_double_1'
  integer walk_i

! begin

! only for first iteration
  if (step_iterations_nb == 1) then

!   initialization
    call require (lhere, 'objects(object_ind)%dimensions(1) > 0', objects(object_ind)%dimensions(1) > 0) !fp
    call alloc ('objects(object_ind)%sum_double_1', objects(object_ind)%sum_double_1, objects(object_ind)%dimensions(1))
    objects(object_ind)%sum_double_1 (:) = 0.d0

  endif ! first iteration

! at each step
  do walk_i = 1, nwalk
   objects(object_ind)%sum_double_1 (:) = objects(object_ind)%sum_double_1 (:) + objects(object_ind)%pointer_double_2 (:,walk_i) * walker_weights(walk_i)
  enddo

 end subroutine object_average_walk_step_by_index_double_1

! ===================================================================================
  subroutine object_average_walk_step_by_index_double_2 (object_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate average of object over MC iterations
!
! Created       : J. Toulouse, 15 Nov 2006
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_average_walk_step_by_index_double_2'
  integer walk_i

! begin

! only for first iteration
  if (step_iterations_nb == 1) then

!   initialization
    call require (lhere, 'objects(object_ind)%dimensions(1) > 0', objects(object_ind)%dimensions(1) > 0) !fp
    call require (lhere, 'objects(object_ind)%dimensions(2) > 0', objects(object_ind)%dimensions(2) > 0) !fp
    call alloc ('objects(object_ind)%sum_double_2', objects(object_ind)%sum_double_2, objects(object_ind)%dimensions(1), objects(object_ind)%dimensions(2))
    objects(object_ind)%sum_double_2 (:,:) = 0.d0

  endif ! first iteration

! at each step
  do walk_i = 1, nwalk
   objects(object_ind)%sum_double_2 (:,:) = objects(object_ind)%sum_double_2 (:,:) + objects(object_ind)%pointer_double_3 (:,:,walk_i) * walker_weights(walk_i)
  enddo

 end subroutine object_average_walk_step_by_index_double_2

! ===================================================================================
  subroutine object_average_walk_block_by_index_double_0 (object_ind, object_bav_ind, object_av_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate average of object over MC iterations
!
! Created       : J. Toulouse, 12 Nov 2006
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_ind, object_bav_ind, object_av_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_average_walk_block_by_index_double_0'
# if defined (MPI)
  integer ierr
  real(dp) collect
# endif

! begin

! for first block
  if (block_iterations_nb == 1 ) then

!   test association
    call object_associated_or_die_by_index (object_ind)
    call object_associated_or_die_by_index (object_bav_ind)
    call object_associated_or_die_by_index (object_av_ind)

!   test dimensions
    if (size(objects(object_ind)%dimensions) - 1 /= size(objects(object_av_ind)%dimensions) .or. size(objects(object_bav_ind)%dimensions) /= size(objects(object_av_ind)%dimensions)) then
      write(6,'(4a,i3)') trim(lhere),': object ',trim(objects(object_ind)%name),' is of dimension ',size(objects(object_ind)%dimensions)
      write(6,'(4a,i3)') trim(lhere),': object ',trim(objects(object_bav_ind)%name),' is of dimension ',size(objects(object_bav_ind)%dimensions)
      write(6,'(4a,i3)') trim(lhere),': object ',trim(objects(object_av_ind)%name),' is of dimension ',size(objects(object_av_ind)%dimensions)
      call die (lhere)
    endif

!   initialization for first block
    objects(object_ind)%sum_blk_double_0 = 0.d0
  endif

# if defined (MPI)
! sum values from all processes
  call mpi_allreduce(objects(object_ind)%sum_double_0,collect,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
  if (ierr /= 0) then
     call die (lhere, 'error in mpi_allreduce')
  endif
  objects(object_ind)%sum_double_0 = collect
# endif

! calculate block average
  objects(object_bav_ind)%pointer_double_0 = objects(object_ind)%sum_double_0 / walker_weights_sum_block
  call object_modified_by_index (object_bav_ind)

! calculate average
  objects(object_ind)%sum_blk_double_0 = objects(object_ind)%sum_blk_double_0 + objects(object_ind)%sum_double_0
  objects(object_av_ind)%pointer_double_0 = objects(object_ind)%sum_blk_double_0 / walker_weights_sum
  call object_modified_by_index (object_av_ind)

! reinitialization of sum for each block
  objects(object_ind)%sum_double_0 = 0.d0

 end subroutine object_average_walk_block_by_index_double_0

! ===================================================================================
  subroutine object_average_walk_block_by_index_double_1 (object_ind, object_bav_ind, object_av_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate average of object over MC iterations
!
! Created       : J. Toulouse, 15 Nov 2006
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_ind, object_bav_ind, object_av_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_average_walk_block_by_index_double_1'

# if defined (MPI)
  integer ierr
  real(dp), allocatable :: collect (:)
# endif

! begin

! for first block
  if (block_iterations_nb == 1 ) then

!   test association
    call object_associated_or_die_by_index (object_ind)
    call object_associated_or_die_by_index (object_bav_ind)
    call object_associated_or_die_by_index (object_av_ind)

!   test dimensions
    if (size(objects(object_ind)%dimensions) - 1 /= size(objects(object_av_ind)%dimensions) .or. size(objects(object_bav_ind)%dimensions) /= size(objects(object_av_ind)%dimensions)) then
      write(6,'(4a,i3)') trim(lhere),': object ',trim(objects(object_ind)%name),' is of dimension ',size(objects(object_ind)%dimensions)
      write(6,'(4a,i3)') trim(lhere),': object ',trim(objects(object_bav_ind)%name),' is of dimension ',size(objects(object_bav_ind)%dimensions)
      write(6,'(4a,i3)') trim(lhere),': object ',trim(objects(object_av_ind)%name),' is of dimension ',size(objects(object_av_ind)%dimensions)
      call die (lhere)
    endif

!   initialization for first block
    call alloc ('objects(object_ind)%sum_blk_double_1', objects(object_ind)%sum_blk_double_1, objects(object_ind)%dimensions(1))
    objects(object_ind)%sum_blk_double_1 (:) = 0.d0
  endif

# if defined (MPI)
!   sum values from all processes
    call alloc ('collect', collect, objects(object_ind)%dimensions(1))
    call mpi_allreduce(objects(object_ind)%sum_double_1,collect,objects(object_ind)%dimensions(1),mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
    if (ierr /= 0) then
       call die (lhere, 'error in mpi_allreduce')
    endif
    objects(object_ind)%sum_double_1 (:) = collect (:)
# endif

! calculate block average
  objects(object_bav_ind)%pointer_double_1 = objects(object_ind)%sum_double_1 / walker_weights_sum_block
  call object_modified_by_index (object_bav_ind)

! calculate global average
  objects(object_ind)%sum_blk_double_1 (:) = objects(object_ind)%sum_blk_double_1 (:) + objects(object_ind)%sum_double_1 (:)
  objects(object_av_ind)%pointer_double_1 (:) = objects(object_ind)%sum_blk_double_1 (:) / walker_weights_sum
  call object_modified_by_index (object_av_ind)

! reinitialization of sum for each block
  objects(object_ind)%sum_double_1 (:) = 0.d0

 end subroutine object_average_walk_block_by_index_double_1

! ===================================================================================
  subroutine object_average_walk_block_by_index_double_2 (object_ind, object_bav_ind, object_av_ind)
! -----------------------------------------------------------------------------------
! Description   : calculate average of object over MC iterations
!
! Created       : J. Toulouse, 15 Nov 2006
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! input
  integer, intent(in) :: object_ind, object_bav_ind, object_av_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'object_average_walk_block_by_index_double_2'
# if defined (MPI)
  integer ierr
  real(dp), allocatable :: collect (:,:)
# endif


! begin

!  for first block
   if (block_iterations_nb == 1 ) then

!   test association
    call object_associated_or_die_by_index (object_ind)
    call object_associated_or_die_by_index (object_av_ind)

!   test dimensions
    if (size(objects(object_ind)%dimensions) - 1 /= size(objects(object_av_ind)%dimensions) .or. size(objects(object_bav_ind)%dimensions) /= size(objects(object_av_ind)%dimensions)) then
      write(6,'(4a,i3)') trim(lhere),': object ',trim(objects(object_ind)%name),' is of dimension ',size(objects(object_ind)%dimensions)
      write(6,'(4a,i3)') trim(lhere),': object ',trim(objects(object_bav_ind)%name),' is of dimension ',size(objects(object_bav_ind)%dimensions)
      write(6,'(4a,i3)') trim(lhere),': object ',trim(objects(object_av_ind)%name),' is of dimension ',size(objects(object_av_ind)%dimensions)
      call die (lhere)
    endif

!  initialization for first block
     call alloc ('objects(object_ind)%sum_blk_double_2', objects(object_ind)%sum_blk_double_2, objects(object_ind)%dimensions(1),objects(object_ind)%dimensions(2))
     objects(object_ind)%sum_blk_double_2 (:,:) = 0.d0
   endif

# if defined (MPI)
!   sum values from all processes
    call alloc ('collect', collect, objects(object_ind)%dimensions(1), objects(object_ind)%dimensions(2))
    call mpi_allreduce(objects(object_ind)%sum_double_2,collect,objects(object_ind)%dimensions(1)*objects(object_ind)%dimensions(2),mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
     if (ierr /= 0) then
        call die (lhere, 'error in mpi_allreduce')
     endif
     objects(object_ind)%sum_double_2 (:,:) = collect (:,:)
# endif

! calculate block average
  objects(object_bav_ind)%pointer_double_2 = objects(object_ind)%sum_double_2 / walker_weights_sum_block
  call object_modified_by_index (object_bav_ind)

!  calculate average
   objects(object_ind)%sum_blk_double_2 (:,:) = objects(object_ind)%sum_blk_double_2 (:,:) + objects(object_ind)%sum_double_2 (:,:)
   objects(object_av_ind)%pointer_double_2 (:,:) = objects(object_ind)%sum_blk_double_2 (:,:) / walker_weights_sum
   call object_modified_by_index (object_av_ind)

! reinitialization of sum for each block
  objects(object_ind)%sum_double_2 (:,:) = 0.d0

 end subroutine object_average_walk_block_by_index_double_2

! ===================================================================================
  subroutine print_list_of_averages_and_errors
! -----------------------------------------------------------------------------------
! Description   : print list of averages, variances, covariances and statistcal errors to be calculated
!
! Created       : J. Toulouse, 25 Mar 2007
! Modified      : J. Toulouse, 21 Dec 2007 : add variances and covariances
! -----------------------------------------------------------------------------------
  implicit none

! local
  integer ind

! begin
  if (averages_nb > 0 .or. averages_walk_nb > 0) then
   write(6,'(a)') 'The following averages will be calculated:'
   do ind = 1, averages_nb
    if (averages_object_index (ind) /= 0) then
     if (objects(averages_object_av_index(ind))%unweighted) then
      write(6,'(4a)') '- ', objects(averages_object_av_index (ind))%name,' = unweighted average of ', trim(objects(averages_object_index (ind))%name)
     else
      write(6,'(4a)') '- ', objects(averages_object_av_index (ind))%name,' = average of ', trim(objects(averages_object_index (ind))%name)
     endif
    else
     write(6,'(4a)') '- ', objects(averages_object_av_index (ind))%name,' = average associated with block average ', trim(objects(averages_object_bav_index (ind))%name)
    endif
   enddo
   do ind = 1, averages_walk_nb
    write(6,'(4a)') '- ', objects(averages_walk_object_av_index (ind))%name,' = average of ', trim(objects(averages_walk_object_index (ind))%name)
   enddo
  endif

!  if (block_averages_nb > 0) then
!   write(6,*)
!   write(6,'(a)') 'The following block averages will be calculated:'
!   do ind = 1, block_averages_nb
!     if (objects(block_averages_object_bav_index(ind))%unweighted) then
!      write(6,'(4a)') '- ', objects(block_averages_object_bav_index (ind))%name,' = unweighted block average of ', trim(objects(block_averages_object_index (ind))%name)
!     else
!      write(6,'(4a)') '- ', objects(block_averages_object_bav_index (ind))%name,' = block average of ', trim(objects(block_averages_object_index (ind))%name)
!     endif
!   enddo
!  endif

  if (variances_nb > 0) then
   write(6,*)
   write(6,'(a)') 'The following variances will be calculated:'
   do ind = 1, variances_nb
    write(6,'(4a)') '- ', objects(variances_object_var_index (ind))%name,' = variance of average ', trim(objects(variances_object_av_index (ind))%name)
   enddo
  endif

  if (covariances_nb > 0) then
   write(6,*)
   write(6,'(a)') 'The following covariances will be calculated:'
   do ind = 1, covariances_nb
    write(6,'(6a)') '- ', objects(covariances_object_covar_index (ind))%name,' = covariance of averages ', trim(objects(covariances_object_av1_index (ind))%name), ' and ',trim(objects(covariances_object_av2_index (ind))%name)
   enddo
  endif

  if (errors_nb > 0) then
   write(6,*)
   write(6,'(a)') 'The following statistical errors will be calculated:'
   do ind = 1, errors_nb
    if (errors_object_av_index (ind) /= 0) then
     write(6,'(4a)') '- ', objects(errors_object_err_index (ind))%name,' = statistical error of average ', trim(objects(errors_object_av_index (ind))%name)
    else
     write(6,'(4a)') '- ', objects(errors_object_err_index (ind))%name,' = statistical error associated with variance ', trim(objects(errors_object_var_index (ind))%name)
    endif
   enddo
  endif
  write(6,*)

 end subroutine print_list_of_averages_and_errors

! ===================================================================================
  subroutine compute_block_averages
! -----------------------------------------------------------------------------------
! Description   : compute block averages in MC iterations
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'compute_block_averages'
  integer ind
  integer object_ind, object_bav_ind
  character(len=max_string_len_type) object_type

! begin

! block averages defined by objects
  do ind = 1, block_averages_nb

    object_ind = block_averages_object_index (ind)
    object_bav_ind = block_averages_object_bav_index (ind)

    call object_provide_by_index (object_ind)

    if (objects(object_ind)%walkers) then
      call die (lhere, 'object has been marked as an array over walkers')
    endif

    object_type = objects(object_ind)%type

    if (trim(object_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

    if (trim(object_type) == 'double_0') then
     call object_associate_by_index_double_0 (object_bav_ind)
     call object_block_average_by_index_double_0 (object_ind, object_bav_ind)

    elseif (trim(object_type) == 'double_1') then
     call object_associate_by_index_double_1 (object_bav_ind, objects(object_ind)%dimensions(1))
     call object_block_average_by_index_double_1 (object_ind, object_bav_ind)

    elseif (trim(object_type) == 'double_2') then
     call object_associate_by_index_double_2 (object_bav_ind, objects(object_ind)%dimensions(1), objects(object_ind)%dimensions(2))
     call object_block_average_by_index_double_2 (object_ind, object_bav_ind)

    else
     call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
    endif

  enddo ! ind

 end subroutine compute_block_averages

! ===================================================================================
  subroutine compute_global_averages
! -----------------------------------------------------------------------------------
! Description   : compute global averages over all blocks
!
! Created       : J. Toulouse, 19 Apr 2008
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'compute_global_averages'
  integer ind, object_bav_ind, object_av_ind
  character(len=max_string_len_type) object_type

! begin

! averages defined by objects
  do ind = 1, averages_nb

    object_bav_ind = averages_object_bav_index (ind)
    object_av_ind = averages_object_av_index (ind)

    call object_provide_by_index (object_bav_ind)

    object_type = objects(object_bav_ind)%type

    if (trim(object_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_bav_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

     if (trim(object_type) == 'double_0') then
      call object_global_average_by_index_double_0 (object_bav_ind, object_av_ind)

     elseif (trim(object_type) == 'double_1') then
      call object_global_average_by_index_double_1 (object_bav_ind, object_av_ind)

     elseif (trim(object_type) == 'double_2') then
      call object_global_average_by_index_double_2 (object_bav_ind, object_av_ind)

     else
      call die (lhere, 'object type >'+trim(object_type)+'< not handled.')

     endif

  enddo ! ind

 end subroutine compute_global_averages

!! ===================================================================================
!  subroutine compute_averages
!! -----------------------------------------------------------------------------------
!! Description   : compute averages in MC iterations
!! Description   : no longer used
!!
!! Created       : J. Toulouse, 20 Oct 2005
!! -----------------------------------------------------------------------------------
!  implicit none
!
!! local
!  character(len=max_string_len_rout), save :: lhere = 'compute_averages'
!  integer ind, rtn_i
!  integer object_ind, object_av_ind
!  character(len=max_string_len_type) object_type
!
!! begin
!
!! averages defined by objects
!  do ind = 1, averages_nb
!
!    object_ind = averages_object_index (ind)
!    object_av_ind = averages_object_av_index (ind)
!
!    call object_provide_by_index (object_ind)
!
!    if (objects(object_ind)%walkers) then
!      write(6,'(2a)') trim(lhere),': object has been marked as an array over walkers'
!      call die (lhere)
!    endif
!
!    object_type = objects(object_ind)%type
!
!    if (trim(object_type) == '') then
!      write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is unknown'
!      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
!      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
!      call die (lhere)
!    endif
!
!     if (trim(object_type) == 'double_0') then
!      call object_average_by_index_double_0 (object_ind, object_av_ind)
!
!     elseif (trim(object_type) == 'double_1') then
!      call object_average_by_index_double_1 (object_ind, object_av_ind)
!
!     elseif (trim(object_type) == 'double_2') then
!      call object_average_by_index_double_2 (object_ind, object_av_ind)
!
!     else
!      write(6,*) trim(lhere),': object type ',trim(object_type),' not handled'
!      call die (lhere)
!
!     endif
!
!  enddo ! ind
!
!! averages computing via defined routines
!  do rtn_i = 1, average_routines_nb
!    call exe_by_address_0 (routines(average_routines_index(rtn_i))%address)
!  enddo !rtn_i
!
! end subroutine compute_averages

! ===================================================================================
  subroutine compute_variances
! -----------------------------------------------------------------------------------
! Description   : compute variances in MC iterations
!
! Created       : J. Toulouse, 06 Sep 2007
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'compute_variances'
  integer ind
  integer object_av_ind, object_var_ind
  character(len=max_string_len_type) object_av_type

! begin
  do ind = 1, variances_nb

    object_av_ind = variances_object_av_index (ind)
    object_var_ind = variances_object_var_index (ind)
!    write(6,'(4a)') 'computing variance >',trim(objects(object_var_ind)%name),'< of average >',trim(objects(object_av_ind)%name),'<'

    call object_provide_by_index (object_av_ind)

    object_av_type = objects(object_av_ind)%type

    if (trim(object_av_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_av_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

    if (trim(object_av_type) == 'double_0') then

      call object_associate_by_index_double_0 (object_var_ind)
      call object_variance_by_index_double_0 (object_av_ind, object_var_ind)

    elseif (trim(object_av_type) == 'double_1') then

      call object_associate_by_index_double_1 (object_var_ind, objects(object_av_ind)%dimensions(1))
      call object_variance_by_index_double_1 (object_av_ind, object_var_ind)

    elseif (trim(object_av_type) == 'double_2') then

      call object_associate_by_index_double_2 (object_var_ind, objects(object_av_ind)%dimensions(1), objects(object_av_ind)%dimensions(2))
      call object_variance_by_index_double_2 (object_av_ind, object_var_ind)

    else
     call die (lhere, 'object type >'+trim(object_av_type)+'< not handled.')
    endif

  enddo ! ind

 end subroutine compute_variances

! ===================================================================================
  subroutine compute_covariances
! -----------------------------------------------------------------------------------
! Description   : compute covariances
!
! Created       : J. Toulouse, 19 Sep 2007
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'compute_covariances'
  integer ind
  integer object_av1_ind, object_av2_ind, object_covar_ind
  character(len=max_string_len_type) object_av1_type, object_av2_type, object_covar_type

! begin
  do ind = 1, covariances_nb

    object_av1_ind = covariances_object_av1_index (ind)
    object_av2_ind = covariances_object_av2_index (ind)
    object_covar_ind = covariances_object_covar_index (ind)
!    write(6,'(6a)') 'computing covariance >',trim(objects(object_covar_ind)%name),'< of averages >',trim(objects(object_av1_ind)%name),'< and >',trim(objects(object_av2_ind)%name),'<'

    call object_provide_by_index (object_av1_ind)
    call object_provide_by_index (object_av2_ind)

    object_av1_type = objects(object_av1_ind)%type
    object_av2_type = objects(object_av2_ind)%type
    object_covar_type = objects(object_covar_ind)%type

    if (trim(object_av1_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_av1_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif
    if (trim(object_av2_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_av2_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif
    if (trim(object_covar_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_covar_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

    if (trim(object_av1_type) == 'double_0' .and. trim(object_av2_type) == 'double_0') then
      call object_covariance_by_index_double_0_double_0 (object_av1_ind, object_av2_ind, object_covar_ind)

    elseif (trim(object_av1_type) == 'double_1' .and. trim(object_av2_type) == 'double_0') then
      call object_covariance_by_index_double_1_double_0 (object_av1_ind, object_av2_ind, object_covar_ind)

    elseif (trim(object_av1_type) == 'double_0' .and. trim(object_av2_type) == 'double_1') then
      call object_covariance_by_index_double_0_double_1 (object_av1_ind, object_av2_ind, object_covar_ind)

    elseif (trim(object_av1_type) == 'double_1' .and. trim(object_av2_type) == 'double_1' .and. trim(object_covar_type) == 'double_1') then
      call object_covariance_by_index_double_1_double_1_double_1 (object_av1_ind, object_av2_ind, object_covar_ind)

    elseif (trim(object_av1_type) == 'double_1' .and. trim(object_av2_type) == 'double_1' .and. trim(object_covar_type) == 'double_2') then
      call object_covariance_by_index_double_1_double_1_double_2 (object_av1_ind, object_av2_ind, object_covar_ind)

    else
     call die (lhere, 'objects types >'+trim(object_av1_type)+'<, >'+trim(object_av2_type)+'< and >'+trim(object_covar_type)+'< not handled.')
    endif

  enddo ! ind

 end subroutine compute_covariances

! ===================================================================================
  subroutine compute_errors
! -----------------------------------------------------------------------------------
! Description   : compute statistical errors in MC iterations
!
! Created       : J. Toulouse, 20 Oct 2005
! -----------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'compute_errors'
  integer ind, object_var_ind, object_err_ind
  character(len=max_string_len_type) object_var_type

! begin
  do ind = 1, errors_nb

    object_var_ind = errors_object_var_index (ind)
    object_err_ind = errors_object_err_index (ind)

    call object_provide_by_index (object_var_ind)

    object_var_type = objects(object_var_ind)%type

    if (trim(object_var_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_var_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

    if (trim(object_var_type) == 'double_0') then

      call object_error_by_index_double_0 (object_var_ind, object_err_ind)

    elseif (trim(object_var_type) == 'double_1') then

      call object_error_by_index_double_1 (object_var_ind, object_err_ind)

    elseif (trim(object_var_type) == 'double_2') then

      call object_error_by_index_double_2 (object_var_ind, object_err_ind)

    else

     call die (lhere, 'object type >'+trim(object_var_type)+'< not handled.')

    endif

  enddo ! ind

 end subroutine compute_errors

! ===================================================================================
  subroutine compute_averages_walk_step
! -----------------------------------------------------------------------------------
! Description   : compute sums at each steps of DMC iterations for averages
!
! Created       : J. Toulouse, 12 Nov 2006
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'compute_averages_walk_step'
  integer ind
  integer object_ind
  character(len=max_string_len_type) object_type

! begin
# if defined (DEBUG)
   call routine_enter (lhere)
# endif

  if (averages_walk_nb == 0) then
   return
  endif

! weights of walkers
  call object_provide_by_index (walker_weights_index)

! averages defined by objects
  do ind = 1, averages_walk_nb

    object_ind = averages_walk_object_index (ind)

    call object_provide_by_index (object_ind)

!    if (.not. objects(object_ind)%walkers) then
!      write(6,'(2a)') trim(lhere),': object as not be marked as a array over walkers'
!      call die (lhere)
!    endif

    object_type = objects(object_ind)%type
    if (trim(object_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

     if (trim(object_type) == 'double_1') then
      call object_average_walk_step_by_index_double_0 (object_ind)

     elseif (trim(object_type) == 'double_2') then
      call object_average_walk_step_by_index_double_1 (object_ind)

     elseif (trim(object_type) == 'double_3') then
      call object_average_walk_step_by_index_double_2 (object_ind)

     else
      call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
     endif

  enddo ! ind

# if defined (DEBUG)
   call routine_exit (lhere)
# endif

 end subroutine compute_averages_walk_step

! ===================================================================================
  subroutine compute_averages_walk_block
! -----------------------------------------------------------------------------------
! Description   : compute averages in DMC iterations for each block
!
! Created       : J. Toulouse, 12 Nov 2006
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'compute_averages_walk_block'
  integer object_ind, object_bav_ind, object_av_ind, ind
  character(len=max_string_len_type) object_type

! begin
# if defined (DEBUG)
   call routine_enter (lhere)
# endif

  if (averages_walk_nb == 0) then
   return
  endif

! sum of weights of walkers over current block
  call object_provide_by_index (walker_weights_sum_block_index)

! sum of weights of walkers over entire run
  call object_provide_by_index (walker_weights_sum_index)

! averages defined by objects
  do ind = 1, averages_walk_nb

    object_ind = averages_walk_object_index (ind)
    object_bav_ind = averages_walk_object_bav_index (ind) ! block average
    object_av_ind = averages_walk_object_av_index (ind)

!    if (.not. objects(object_ind)%walkers) then
!      write(6,'(2a)') trim(lhere),': object as not be marked as a array over walkers'
!      call die (lhere)
!    endif

    object_type = objects(object_ind)%type
    if (trim(object_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

     if (trim(object_type) == 'double_1') then
      call object_associate_by_index_double_0 (object_bav_ind)
      call object_average_walk_block_by_index_double_0 (object_ind, object_bav_ind, object_av_ind)

     elseif (trim(object_type) == 'double_2') then
      call object_associate_by_index_double_1 (object_bav_ind, objects(object_av_ind)%dimensions(1))
      call object_average_walk_block_by_index_double_1 (object_ind, object_bav_ind, object_av_ind)

     elseif (trim(object_type) == 'double_3') then
      call object_associate_by_index_double_2 (object_bav_ind, objects(object_av_ind)%dimensions(1), objects(object_av_ind)%dimensions(2))
      call object_average_walk_block_by_index_double_2 (object_ind, object_bav_ind, object_av_ind)

     else
      call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
     endif

  enddo ! ind

# if defined (DEBUG)
   call routine_exit (lhere)
# endif

 end subroutine compute_averages_walk_block

! ===================================================================================
  subroutine compute_averages_step
! -----------------------------------------------------------------------------------
! Description   : calculate weighted sums of objects over VMC or DMC iterations
! Description   : for DMC, the sum of walkers is done outside

! Description   : this routine should eventually replace all the other averages_step routines
!
! Created       : J. Toulouse, 27 Mar 2010
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'compute_averages_step'
  integer ind, object_ind, object_bav_ind
  character(len=max_string_len_type) object_type

  if (l_equilibration) return

! initialization for first walker and first step
  if (current_walker == 1 .and. step_iterations_nb == 1) then
   total_iterations_block_nb = 0.d0
   total_iterations_nb = 0.d0
  endif

  total_iterations_block_nb = total_iterations_block_nb + 1
  call object_provide_by_index (current_walker_index)
  call object_provide_by_index (current_walker_weight_index)

! loop over block averages
  do ind = 1, block_averages_nb

    object_ind = block_averages_object_index (ind)
    object_bav_ind = block_averages_object_bav_index (ind)

    call object_provide_by_index (object_ind)

    if (objects(object_ind)%walkers) then
      call die (lhere, 'object has been marked as an array over walkers')
    endif

    object_type = objects(object_ind)%type
    if (trim(object_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

!   initialization for first walker and first step
    if (current_walker == 1 .and. step_iterations_nb == 1) then
      select case (trim(object_type))
      case ('double_0')
        objects(object_bav_ind)%sum_double_0 = 0.d0
      case ('double_1')
        call require (lhere, 'objects(object_ind)%dimensions(1) > 0', objects(object_ind)%dimensions(1) > 0)
        call alloc ('objects(object_bav_ind)%sum_double_1', objects(object_bav_ind)%sum_double_1, objects(object_ind)%dimensions(1))
        objects(object_bav_ind)%sum_double_1 (:) = 0.d0
      case ('double_2')
        call require (lhere, 'objects(object_ind)%dimensions(1) > 0', objects(object_ind)%dimensions(1) > 0)
        call require (lhere, 'objects(object_ind)%dimensions(2) > 0', objects(object_ind)%dimensions(2) > 0)
        call alloc ('objects(object_bav_ind)%sum_double_2', objects(object_bav_ind)%sum_double_2, objects(object_ind)%dimensions(1), objects(object_ind)%dimensions(2))
        objects(object_bav_ind)%sum_double_2 (:,:) = 0.d0
      case default
        call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
      end select
    endif ! initialization

!   sum over walkers and steps
     select case (trim(object_type))
     case ('double_0')
      if (objects(object_bav_ind)%unweighted) then
       objects(object_bav_ind)%sum_double_0 = objects(object_bav_ind)%sum_double_0 + objects(object_ind)%pointer_double_0 
      else
       objects(object_bav_ind)%sum_double_0 = objects(object_bav_ind)%sum_double_0 + objects(object_ind)%pointer_double_0 * current_walker_weight
      endif
     case ('double_1')
      if (objects(object_bav_ind)%unweighted) then
       objects(object_bav_ind)%sum_double_1 = objects(object_bav_ind)%sum_double_1 + objects(object_ind)%pointer_double_1 
      else
       objects(object_bav_ind)%sum_double_1 = objects(object_bav_ind)%sum_double_1 + objects(object_ind)%pointer_double_1 * current_walker_weight
      endif
     case ('double_2')
      if (objects(object_bav_ind)%unweighted) then
       objects(object_bav_ind)%sum_double_2 = objects(object_bav_ind)%sum_double_2 + objects(object_ind)%pointer_double_2 
      else
       objects(object_bav_ind)%sum_double_2 = objects(object_bav_ind)%sum_double_2 + objects(object_ind)%pointer_double_2 * current_walker_weight
      endif
     case default
       call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
     end select

  enddo ! ind

 end subroutine compute_averages_step

! ===================================================================================
  subroutine compute_averages_block
! -----------------------------------------------------------------------------------
! Description   : compute block and global averages in VMC or DMC
! Description   : this routine must be called at each block
! Description   : the block average of a local object X_ik where k refers blocks, i to steps is
! Description   : X_k = sum_i (X_ik w_ik ) / sum_i w_ik
! Description   : where w_ik are weights (= 1 for VMC)
! Description   : the global average is
! Description   : X = sum_k sum_i ( X_ik w_ik ) / sum_k sum_i w_ik
! Description   :   = sum_k (X_k * sum_i w_ik ) / sum_k sum_i w_ik
! Description   :  
! Description   : In DMC, for block averages only defined from other bock averages (in this case object_ind=0),
! Description   : maybe it would make more sense to calculate them simply as
! Description   : X = sum_k X_k / Nblock ???
!
! Description   : this routine should eventually replace all the other averages_block routines
!
! Created       : J. Toulouse, 27 Mar 2010
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'compute_averages_block'
  integer object_ind, object_bav_ind, object_av_ind, ind
  character(len=max_string_len_type)   :: object_type, object_bav_type, object_av_type
# if defined (MPI)
  integer ierr
  real(dp) collect_double_0
  real(dp), allocatable :: collect_double_1 (:)
  real(dp), allocatable :: collect_double_2 (:,:)
# endif

# if defined (MPI)
!   sum values from all processes
    call mpi_allreduce(total_iterations_block_nb,collect_double_0,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
    if (ierr /= 0) call die (lhere, 'error in mpi_allreduce')
    total_iterations_block_nb = collect_double_0
# endif
  total_iterations_nb = total_iterations_nb + total_iterations_block_nb
  call object_modified_by_index (total_iterations_nb_index)

!JT  if (block_averages_nb == 0) return

! sum of weights of walkers over current block
  call object_provide_by_index (walker_weights_sum_block_index)

! sum of weights of walkers over entire run
  call object_provide_by_index (walker_weights_sum_index)

! loop over block averages ------------------------------------------------------------------------------------------
  do ind = 1, block_averages_nb

    object_ind = block_averages_object_index (ind)
    object_bav_ind = block_averages_object_bav_index (ind)

    object_type = objects(object_ind)%type
    if (trim(object_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

!   for first block: associate block average if necessary and test types and dimensions
    if (block_iterations_nb == 1) then

      select case (trim(object_type))
      case ('double_0') 
        call object_associate_by_index_double_0 (object_bav_ind)
        call require (lhere, 'type of '+trim(objects(object_ind)%name)+' = type of '+trim(objects(object_bav_ind)%name), objects(object_ind)%type == objects(object_bav_ind)%type)
      case ('double_1') 
        call object_associate_by_index_double_1 (object_bav_ind, objects(object_ind)%dimensions(1))
        call require (lhere, 'type of '+trim(objects(object_ind)%name)+' = type of '+trim(objects(object_bav_ind)%name), objects(object_ind)%type == objects(object_bav_ind)%type)
        call require (lhere, 'dim # 1 of '+trim(objects(object_ind)%name)+' = dim # 1 of '+trim(objects(object_bav_ind)%name), objects(object_ind)%dimensions(1) == objects(object_bav_ind)%dimensions(1))
      case ('double_2') 
        call object_associate_by_index_double_2 (object_bav_ind, objects(object_ind)%dimensions(1), objects(object_ind)%dimensions(2))
        call require (lhere, 'type of '+trim(objects(object_ind)%name)+' = type of '+trim(objects(object_bav_ind)%name), objects(object_ind)%type == objects(object_bav_ind)%type)
        call require (lhere, 'dim # 1 of '+trim(objects(object_ind)%name)+' = dim # 1 of '+trim(objects(object_bav_ind)%name), objects(object_ind)%dimensions(1) == objects(object_bav_ind)%dimensions(1))
        call require (lhere, 'dim # 2 of '+trim(objects(object_ind)%name)+' = dim # 2 of '+trim(objects(object_bav_ind)%name), objects(object_ind)%dimensions(2) == objects(object_bav_ind)%dimensions(2))
      case default
        call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
      end select

    endif ! first block

    select case (trim(object_type))
    case ('double_0') 
#   if defined (MPI)
!   sum values from all processes
    call mpi_allreduce(objects(object_bav_ind)%sum_double_0,collect_double_0,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
    if (ierr /= 0) call die (lhere, 'error in mpi_allreduce')
    objects(object_bav_ind)%sum_double_0 = collect_double_0
#   endif
!   calculate block average
    if (objects(object_bav_ind)%unweighted) then
     objects(object_bav_ind)%pointer_double_0 = objects(object_bav_ind)%sum_double_0 / total_iterations_block_nb
    else
     objects(object_bav_ind)%pointer_double_0 = objects(object_bav_ind)%sum_double_0 / walker_weights_sum_block
    endif
!   reinitialization of sum for each block
    objects(object_bav_ind)%sum_double_0 = 0.d0

    case ('double_1') 
#   if defined (MPI)
!   sum values from all processes
    call alloc ('collect_double_1', collect_double_1, objects(object_ind)%dimensions(1))
    call mpi_allreduce(objects(object_bav_ind)%sum_double_1,collect_double_1,objects(object_ind)%dimensions(1),mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
    if (ierr /= 0) call die (lhere, 'error in mpi_allreduce')
    objects(object_bav_ind)%sum_double_1 = collect_double_1
#   endif
!   calculate block average
    if (objects(object_bav_ind)%unweighted) then
     objects(object_bav_ind)%pointer_double_1 = objects(object_bav_ind)%sum_double_1 / total_iterations_block_nb
    else
     objects(object_bav_ind)%pointer_double_1 = objects(object_bav_ind)%sum_double_1 / walker_weights_sum_block
    endif
!   reinitialization of sum for each block
    objects(object_bav_ind)%sum_double_1 = 0.d0

    case ('double_2') 
#   if defined (MPI)
!   sum values from all processes
    call alloc ('collect_double_2', collect_double_2, objects(object_ind)%dimensions(1), objects(object_ind)%dimensions(2))
    call mpi_allreduce(objects(object_bav_ind)%sum_double_2,collect_double_2,objects(object_ind)%dimensions(1)*objects(object_ind)%dimensions(2),mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
    if (ierr /= 0) call die (lhere, 'error in mpi_allreduce')
    objects(object_bav_ind)%sum_double_2 = collect_double_2
#   endif
!   calculate block average
    if (objects(object_bav_ind)%unweighted) then
     objects(object_bav_ind)%pointer_double_2 = objects(object_bav_ind)%sum_double_2 / total_iterations_block_nb
    else
     objects(object_bav_ind)%pointer_double_2 = objects(object_bav_ind)%sum_double_2 / walker_weights_sum_block
    endif
!   reinitialization of sum for each block
    objects(object_bav_ind)%sum_double_2 = 0.d0
    case default
      call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
    end select

    call object_modified_by_index (object_bav_ind)

  enddo ! end loop over block averages ---------------------------------------------------------------------------------------

! loop over global averages --------------------------------------------------------------------------------------------------
! a global average is calculated as average of a block average
  do ind = 1, averages_nb

    object_bav_ind = averages_object_bav_index (ind)
    object_av_ind = averages_object_av_index (ind)

    call object_provide_by_index (object_bav_ind)

    object_type = objects(object_bav_ind)%type
    if (trim(object_type) == '') then
      write(6,*) trim(lhere),': type of object ',trim(objects(object_bav_ind)%name),' is unknown'
      write(6,*) trim(lhere),': the object name has probably not be associated with a address'
      write(6,*) trim(lhere),': by using either object_alloc or object_associate'
      call die (lhere)
    endif

!   for first block: associate block average if necessary and test types and dimensions
    if (block_iterations_nb == 1) then

!     test association
      call object_associated_or_die_by_index (object_bav_ind)
      call object_associated_or_die_by_index (object_av_ind)

!     test type
      call require (lhere, 'type of '+trim(objects(object_bav_ind)%name)+' = type of '+trim(objects(object_av_ind)%name), objects(object_bav_ind)%type == objects(object_av_ind)%type)

!     test dimensions and initialization of sum over blocks
      select case (trim(object_type))
      case ('double_0') 
        objects(object_av_ind)%sum_double_0 = 0.d0
      case ('double_1') 
        call require (lhere, 'dim # 1 of '+trim(objects(object_bav_ind)%name)+' = dim # 1 of '+trim(objects(object_av_ind)%name), objects(object_bav_ind)%dimensions(1) == objects(object_av_ind)%dimensions(1))
        call alloc ('objects(object_av_ind)%sum_double_1', objects(object_av_ind)%sum_double_1, objects(object_bav_ind)%dimensions(1))
        objects(object_av_ind)%sum_double_1 = 0.d0
      case ('double_2') 
        call require (lhere, 'dim # 1 of '+trim(objects(object_bav_ind)%name)+' = dim # 1 of '+trim(objects(object_av_ind)%name), objects(object_bav_ind)%dimensions(1) == objects(object_av_ind)%dimensions(1))
        call require (lhere, 'dim # 2 of '+trim(objects(object_bav_ind)%name)+' = dim # 2 of '+trim(objects(object_av_ind)%name), objects(object_bav_ind)%dimensions(2) == objects(object_av_ind)%dimensions(2))
        call alloc ('objects(object_av_ind)%sum_double_2', objects(object_av_ind)%sum_double_2, objects(object_bav_ind)%dimensions(1), objects(object_bav_ind)%dimensions(2))
        objects(object_av_ind)%sum_double_2 = 0.d0
      case default
        call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
      end select

    endif ! first block

    select case (trim(object_type))
    case ('double_0') 
     if (objects(object_bav_ind)%unweighted) then
!     sum over blocks
      objects(object_av_ind)%sum_double_0 = objects(object_av_ind)%sum_double_0 + objects(object_bav_ind)%pointer_double_0 * total_iterations_block_nb
!     calculate global average
      objects(object_av_ind)%pointer_double_0 = objects(object_av_ind)%sum_double_0 / total_iterations_nb
     else
!     sum over blocks
      objects(object_av_ind)%sum_double_0 = objects(object_av_ind)%sum_double_0 + objects(object_bav_ind)%pointer_double_0 * walker_weights_sum_block
!     calculate global average
      objects(object_av_ind)%pointer_double_0 = objects(object_av_ind)%sum_double_0 / walker_weights_sum
     endif

    case ('double_1') 
     if (objects(object_bav_ind)%unweighted) then
!     sum over blocks
      objects(object_av_ind)%sum_double_1 = objects(object_av_ind)%sum_double_1 + objects(object_bav_ind)%pointer_double_1 * total_iterations_block_nb
!     calculate global average
      objects(object_av_ind)%pointer_double_1 = objects(object_av_ind)%sum_double_1 / total_iterations_nb
     else
!     sum over blocks
      objects(object_av_ind)%sum_double_1 = objects(object_av_ind)%sum_double_1 + objects(object_bav_ind)%pointer_double_1 * walker_weights_sum_block
!     calculate global average
      objects(object_av_ind)%pointer_double_1 = objects(object_av_ind)%sum_double_1 / walker_weights_sum
     endif

    case ('double_2') 
     if (objects(object_bav_ind)%unweighted) then
!     sum over blocks
      objects(object_av_ind)%sum_double_2 = objects(object_av_ind)%sum_double_2 + objects(object_bav_ind)%pointer_double_2 * total_iterations_block_nb
!     calculate global average
      objects(object_av_ind)%pointer_double_2 = objects(object_av_ind)%sum_double_2 / total_iterations_nb
     else
!     sum over blocks
      objects(object_av_ind)%sum_double_2 = objects(object_av_ind)%sum_double_2 + objects(object_bav_ind)%pointer_double_2 * walker_weights_sum_block
!     calculate global average
      objects(object_av_ind)%pointer_double_2 = objects(object_av_ind)%sum_double_2 / walker_weights_sum
     endif

    case default
      call die (lhere, 'object type >'+trim(object_type)+'< not handled.')
    end select

    call object_modified_by_index (object_av_ind)

  enddo ! end loop over global averages --------------------------------------------------------------------------------------

! reinitilizations at the end of each block
  total_iterations_block_nb = 0.d0

 end subroutine compute_averages_block

! ===================================================================================
  subroutine reinit_averages_and_errors
! -----------------------------------------------------------------------------------
! Description   : reinitialize array of averages and errors
!
! Created       : J. Toulouse, 04 Nov 2005
! -----------------------------------------------------------------------------------
  implicit none

! local

! begin
  block_averages_nb = 0
  call release ('block_averages_object_index', block_averages_object_index)
  call release ('block_averages_object_bav_index', block_averages_object_bav_index)

  averages_nb = 0
  call release ('averages_object_index', averages_object_index)
  call release ('averages_object_bav_index', averages_object_bav_index)
  call release ('averages_object_av_index', averages_object_av_index)

  variances_nb = 0
  call release ('variances_object_av_index', variances_object_av_index)
  call release ('variances_object_var_index', variances_object_var_index)

  covariances_nb = 0
  call release ('covariances_object_av1_index', covariances_object_av1_index)
  call release ('covariances_object_av2_index', covariances_object_av2_index)
  call release ('covariances_object_covar_index', covariances_object_covar_index)

  errors_nb = 0
  call release ('errors_object_av_index', errors_object_av_index)
  call release ('errors_object_var_index', errors_object_var_index)
  call release ('errors_object_err_index', errors_object_err_index)

! reinitialize average routines
  average_routines_nb = 0
  call release ('average_routines_index', average_routines_index)

 end subroutine reinit_averages_and_errors

! ===================================================================================
  subroutine save_averages_and_errors
! -----------------------------------------------------------------------------------
! Description   : save arrays of averages and errors
!
! Created       : J. Toulouse, 27 Mar 2010
! -----------------------------------------------------------------------------------
  implicit none

! local

! begin
  block_averages_nb_save = block_averages_nb
  call copy (block_averages_object_index, block_averages_object_index_save)
  call copy (block_averages_object_bav_index, block_averages_object_bav_index_save)

  averages_nb_save = averages_nb
  call copy (averages_object_index, averages_object_index_save)
  call copy (averages_object_bav_index, averages_object_bav_index_save)
  call copy (averages_object_av_index, averages_object_av_index_save)

  variances_nb_save = variances_nb
  call copy (variances_object_av_index, variances_object_av_index_save)
  call copy (variances_object_var_index, variances_object_var_index_save)

  covariances_nb_save = covariances_nb
  call copy (covariances_object_av1_index, covariances_object_av1_index_save)
  call copy (covariances_object_av2_index, covariances_object_av2_index_save)
  call copy (covariances_object_covar_index, covariances_object_covar_index_save)

  errors_nb_save = errors_nb
  call copy (errors_object_av_index, errors_object_av_index_save)
  call copy (errors_object_var_index, errors_object_var_index_save)
  call copy (errors_object_err_index, errors_object_err_index_save)

  average_routines_nb_save = average_routines_nb
  call copy (average_routines_index, average_routines_index_save)

 end subroutine save_averages_and_errors

! ===================================================================================
  subroutine restore_averages_and_errors
! -----------------------------------------------------------------------------------
! Description   : restore arrays of averages and errors
!
! Created       : J. Toulouse, 27 Mar 2010
! -----------------------------------------------------------------------------------
  implicit none

! local

! begin
  block_averages_nb = block_averages_nb_save
  call move (block_averages_object_index_save, block_averages_object_index)
  call move (block_averages_object_bav_index_save, block_averages_object_bav_index)

  averages_nb = averages_nb_save
  call move (averages_object_index_save, averages_object_index)
  call move (averages_object_bav_index_save, averages_object_bav_index)
  call move (averages_object_av_index_save, averages_object_av_index)

  variances_nb = variances_nb_save
  call move (variances_object_av_index_save, variances_object_av_index)
  call move (variances_object_var_index_save, variances_object_var_index)

  covariances_nb = covariances_nb_save
  call move (covariances_object_av1_index_save, covariances_object_av1_index)
  call move (covariances_object_av2_index_save, covariances_object_av2_index)
  call move (covariances_object_covar_index_save, covariances_object_covar_index)

  errors_nb = errors_nb_save
  call move (errors_object_av_index_save, errors_object_av_index)
  call move (errors_object_var_index_save, errors_object_var_index)
  call move (errors_object_err_index_save, errors_object_err_index)

  average_routines_nb = average_routines_nb_save
  call move (average_routines_index_save, average_routines_index)

 end subroutine restore_averages_and_errors

! ===================================================================================
  subroutine invalidate_averages_and_errors
! -----------------------------------------------------------------------------------
! Description   : invalidate averages and errors, etc...
! Description   : at the beginning of a new qmc run
! Description   : not used because in an optimization run with correlated sampling stabilization
! Description   : we need to keep objects depending on averages such as ovlp_lin valid.
!
! Created       : J. Toulouse, 08 Aug 2008
! -----------------------------------------------------------------------------------
  implicit none

! local
  integer block_averages_defined_i, averages_defined_i, variances_defined_i
  integer covariances_defined_i, errors_defined_i

! begin
  do block_averages_defined_i = 1, block_averages_defined_nb
    call object_invalidate_by_index (block_averages_defined_object_bav_index (block_averages_defined_i))
  enddo

  do averages_defined_i = 1, averages_defined_nb
    call object_invalidate_by_index (averages_defined_object_av_index (averages_defined_i))
  enddo

  do variances_defined_i = 1, variances_defined_nb
    call object_invalidate_by_index (variances_defined_object_var_index (variances_defined_i))
  enddo

  do covariances_defined_i = 1, covariances_defined_nb
    call object_invalidate_by_index (covariances_defined_object_covar_index (covariances_defined_i))
  enddo

  do errors_defined_i = 1, errors_defined_nb
    call object_invalidate_by_index (errors_defined_object_err_index (errors_defined_i))
  enddo

 end subroutine invalidate_averages_and_errors

end module average_mod
