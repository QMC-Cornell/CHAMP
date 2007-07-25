module routines_mod

  use basic_tools_mod
  use constants_mod
  use nodes_mod

  type routine
   character(len=max_string_len) :: name
   logical                       :: debug = .false.
   logical                       :: inside = .false.
   integer                       :: calls_nb = 0
   real(dp)                      :: cpu_start
   real(dp)                      :: cpu_end
   real(dp)                      :: cpu_duration = 0.d0
   integer,allocatable           :: routines_up_index (:)
   integer                       :: recursion_level = 0
   integer                       :: skip = 0
   integer*8                     :: address
  end type routine

  integer                         :: routines_nb = 0
  type (routine), save            :: routines (max_routines_nb)
  integer                         :: routine_previous_index = 0

  integer, allocatable            :: routines_indexes_sort (:)
  integer, allocatable            :: routines_calls_nb_sort (:)
  real(dp), allocatable           :: routines_cpu_duration_sort (:)

  contains

! ==============================================================================
  subroutine catalog_one_routine (routine_name, routine, routine_ind)
! ------------------------------------------------------------------------------
! Description   : catalog a new routine (not a node!)
!
! Created       : J. Toulouse, 14 Dec 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)  :: routine_name
  external                      :: routine

! input/output
  integer, intent(inout), optional :: routine_ind

! local
  character(len=max_string_len_rout), save :: lhere = 'catalog_one_routine'
  logical routine_found
  integer rtn_i
  integer rout_ind

! begin

! make sure that routine_ind = 0
  if (present(routine_ind)) then
  if (routine_ind /= 0) then
   write(6,*) trim(lhere), ': routine_name=',trim(routine_name)
   write(6,*) trim(lhere), ': calling with routine_ind=',routine_ind, '/=0'
   call die (lhere)
  endif
  endif

! index of routine
  rout_ind = routine_index (routine_name)

! if routine already catalogued, dies
  if (rout_ind /= 0) then
    call die(here, 'routine >'+trim(routine_name)+'< catalogued more than once.')
  endif

! if current routine not found, catalog routine
  if (routines_nb >= max_routines_nb) then
   write(6,*) trim(lhere),': maximun number of routines max_routines_nb=',max_routines_nb,' reached'
   call die (lhere)
  endif

  routines_nb = routines_nb + 1
  rout_ind = routines_nb

  routines(rout_ind)%name = routine_name
  routines(rout_ind)%address = address(routine)

! save routine index
  if (present(routine_ind)) then
   routine_ind = rout_ind
  endif

 end subroutine catalog_one_routine

! ==============================================================================
  function routine_index (routine_name) result(result)
! ------------------------------------------------------------------------------
! Description   : return the index of a routine
! Description   : 0 if routine not catalogued
!
! Created       : J. Toulouse, 15 Dec 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! output
  integer                      :: result

! local
  integer rtn_i
  integer length, half
  character(len=max_string_len_obj) current_routine_name

! begin
  length = len_trim(routine_name)
  half   = max(1,length/2)

! check if routine is catalogued
   do rtn_i = 1, routines_nb

    current_routine_name = routines(rtn_i)%name

!    if (length /= len_trim(current_routine_name)) cycle  ! too slow
    if (routine_name(1:1) /= current_routine_name(1:1) ) cycle
    if (routine_name(length:length) /= current_routine_name(length:length) ) cycle
    if (routine_name(half:half) /= current_routine_name(half:half) ) cycle

    if (routine_name == current_routine_name ) then
      result = rtn_i
      return
    endif

   enddo

! if object not found
  result = 0

  return

 end function routine_index

! ==============================================================================
  subroutine routine_add (routine_name)
! ------------------------------------------------------------------------------
! Description   : catalogue a new routine
!
! Created       : J. Toulouse, 18 Dec 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! local
  character(len=max_string_len_rout), save :: lhere = 'routine_add'

! begin
  if (routines_nb == max_routines_nb) then
   write(6,*) trim(lhere),': maximun number of routines max_routines_nb=',max_routines_nb,' reached'
   call die (lhere)
  endif

  routines_nb = routines_nb + 1
  routines(routines_nb)%name = routine_name

 end subroutine routine_add

! ==============================================================================
  function routine_catalog_and_index (routine_name) result(routine_ind)
! ------------------------------------------------------------------------------
! Description   : catalog routine if not catalogued and return the index
!
! Created       : J. Toulouse, 13 Apr 2007
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! output
  integer                      :: routine_ind

! local
  routine_ind = routine_index (routine_name)

! if routine not found, catalog it
  if (routine_ind == 0) then
   call routine_add (routine_name)
   routine_ind = routines_nb
  endif

  return
  end function routine_catalog_and_index

!===============================================================================
  subroutine routine_enter_simple (routine_name)
! ------------------------------------------------------------------------------
! Description   : tool for managing routines
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout) :: here
  character(len=*) :: routine_name
  integer routine_ind

! begin
  here ='routine_enter_simple'

! index of current routine
  routine_ind = routine_index (routine_name)

! if routine not found, catalog it
   if (routine_ind == 0) then
    call routine_add (routine_name)
    routine_ind = routines_nb
   endif

  routines(routine_ind)%calls_nb = routines(routine_ind)%calls_nb + 1

! if we have already entered in this routine
   if (routines(routine_ind)%inside) then
      routines(routine_ind)%skip = routines(routine_ind)%skip + 1
   else
    routines(routine_ind)%inside = .true.
    call cpu_time (routines(routine_ind)%cpu_start)
   endif


 end subroutine routine_enter_simple

!===============================================================================
  subroutine routine_exit_simple (routine_name)
! ------------------------------------------------------------------------------
! Description   : tool for managing routines
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout) :: here
  character(len=*) :: routine_name
  integer routine_ind

! begin
  here ='routine_exit_simple'

! index of routine
  routine_ind = routine_index (routine_name)

! if routine not found, die
   if (routine_ind == 0) then
      write(6,*) trim(here),': exiting routine ', trim(routine_name),' that not catalogued'
      call die(here)
   endif

   if (routines(routine_ind)%skip > 0) then
    routines(routine_ind)%skip = routines(routine_ind)%skip - 1
   else
    routines(routine_ind)%inside = .false.
    call cpu_time (routines(routine_ind)%cpu_end)
    routines(routine_ind)%cpu_duration = routines(routine_ind)%cpu_duration + (routines(routine_ind)%cpu_end - routines(routine_ind)%cpu_start)
   endif

 end subroutine routine_exit_simple

!===============================================================================
  subroutine routine_enter (routine_name)
! ------------------------------------------------------------------------------
! Description   : tool for managing routines
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! local
  character(len=max_string_len_rout) :: lhere = 'routine_enter'
  integer routine_ind

! begin
  if (l_trace) then
   write (6,'(2a)') 'trace: entering routine ', trim(routine_name)
  endif

! if we have not exited previous routine yet, record temporary cpu time of previous routine
   if (routine_previous_index /= 0) then
    if (routines(routine_previous_index)%inside) then
     call cpu_time (routines(routine_previous_index)%cpu_end)
     routines(routine_previous_index)%cpu_duration = routines(routine_previous_index)%cpu_duration + (routines(routine_previous_index)%cpu_end - routines(routine_previous_index)%cpu_start)
    endif
   endif

! index of current routine
  routine_ind = routine_catalog_and_index (routine_name)

! set debug
  debug = routines(routine_ind)%debug

! initialize recursion level
  if (routines(routine_ind)%recursion_level == 0) then
    routines(routine_ind)%recursion_level = 1
  endif

! if we have already entered in this routine but not exited, increase the recursivity level
  if (routines(routine_ind)%inside) then
    routines(routine_ind)%recursion_level = routines(routine_ind)%recursion_level + 1
  endif
  call alloc ('routines(routine_ind)%routines_up_index', routines(routine_ind)%routines_up_index, routines(routine_ind)%recursion_level)

  if (routine_previous_index /= 0) then

!   record previous routine as one-level up routine
    if (routines(routine_previous_index)%inside) then
      routines(routine_ind)%routines_up_index (routines(routine_ind)%recursion_level) = routine_previous_index

!   previous routine was at the same level
    else
      routines(routine_ind)%routines_up_index (routines(routine_ind)%recursion_level) = 0
    endif
  endif

  routines(routine_ind)%inside = .true.
  routines(routine_ind)%calls_nb = routines(routine_ind)%calls_nb + 1
  routine_previous_index = routine_ind

  call cpu_time (routines(routine_ind)%cpu_start)

  end subroutine routine_enter

!===============================================================================
  subroutine routine_exit (routine_name)
! ------------------------------------------------------------------------------
! Description   : tool for managing routines
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! local
  character(len=max_string_len_rout), save :: lhere ='routine_exit'
  integer routine_ind

! begin
  if (l_trace) then
   write (6,'(2a)') 'trace: exiting routine ', trim(routine_name)
  endif

! index of routine
  routine_ind = routine_index (routine_name)

! if routine not found, die
   if (routine_ind == 0) then
      call die (lhere, 'exiting routine >'+trim(routine_name)+'< that is not catalogued.')
   endif

! if we have not entered this routine previously
   if (.not. routines(routine_ind)%inside) then
      call die (lhere, 'exiting routine >'+trim(routine_name)+'< but entering has not been notified.')
   endif

!  calculate cpu time
   call cpu_time (routines(routine_ind)%cpu_end)
   routines(routine_ind)%cpu_duration = routines(routine_ind)%cpu_duration + (routines(routine_ind)%cpu_end - routines(routine_ind)%cpu_start)

!  we are no longer inside the routine
   routines(routine_ind)%inside = .false.

! if we come back to the calling routine one level above, restart timing for this routine
  routine_previous_index = routines(routine_ind)%routines_up_index (routines(routine_ind)%recursion_level)
  if (routine_previous_index /= 0) then
     routines(routine_previous_index)%inside = .true.
!    set debug
     debug = routines(routine_previous_index)%debug

!    decrease level of recursion
     routines(routine_ind)%recursion_level = routines(routine_ind)%recursion_level - 1
     call alloc ('routines(routine_ind)%routines_up_index', routines(routine_ind)%routines_up_index, routines(routine_ind)%recursion_level)

!    start cpu time
     call cpu_time (routines(routine_previous_index)%cpu_start)
   endif

 end subroutine routine_exit

!===============================================================================
  subroutine routines_statistics
! ------------------------------------------------------------------------------
! Description   : print statisitics for routines (that are not routines)
!
! Created       : J. Toulouse, 18 Dec 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout) :: here
  character(len=30) :: routine_name
  integer rn_i, rn_j, rn_srt_i
  real(dp) temp
  real(dp) cpu_duration_total

! begin
  here ='routines_statistics'

!  call die (here, 'routine not tested')

  cpu_duration_total = 0.d0

! sort by decreasing cpu times
  call alloc ('routines_indexes_sort', routines_indexes_sort, routines_nb)
  call alloc ('routines_cpu_duration_sort', routines_cpu_duration_sort, routines_nb)
  call alloc ('routines_calls_nb_sort', routines_calls_nb_sort, routines_nb)

  do rn_i = 1, routines_nb
   routines_indexes_sort (rn_i) = rn_i
   routines_cpu_duration_sort (rn_i) = routines(rn_i)%cpu_duration
   routines_calls_nb_sort (rn_i) = routines(rn_i)%calls_nb
  enddo

  do rn_i = 1, routines_nb
   do rn_j = rn_i+1, routines_nb

   if (routines_cpu_duration_sort (rn_j) > routines_cpu_duration_sort (rn_i)) then
    temp = routines_cpu_duration_sort (rn_i)
    routines_cpu_duration_sort (rn_i) = routines_cpu_duration_sort (rn_j)
    routines_cpu_duration_sort (rn_j) = temp
    temp = routines_indexes_sort (rn_i)
    routines_indexes_sort (rn_i) = routines_indexes_sort (rn_j)
    routines_indexes_sort (rn_j) = temp
    temp = routines_calls_nb_sort (rn_i)
    routines_calls_nb_sort (rn_i) = routines_calls_nb_sort (rn_j)
    routines_calls_nb_sort (rn_j) = temp
   endif

   if (routines_cpu_duration_sort (rn_j) == routines_cpu_duration_sort (rn_i)) then
   if (routines_calls_nb_sort (rn_j) > routines_calls_nb_sort (rn_i)) then
    temp = routines_cpu_duration_sort (rn_i)
    routines_cpu_duration_sort (rn_i) = routines_cpu_duration_sort (rn_j)
    routines_cpu_duration_sort (rn_j) = temp
    temp = routines_indexes_sort (rn_i)
    routines_indexes_sort (rn_i) = routines_indexes_sort (rn_j)
    routines_indexes_sort (rn_j) = temp
    temp = routines_calls_nb_sort (rn_i)
    routines_calls_nb_sort (rn_i) = routines_calls_nb_sort (rn_j)
    routines_calls_nb_sort (rn_j) = temp
   endif
   endif

  enddo
  enddo

    write(6,'(a)') '----------------------------------------------------------------------'
    write(6,'(a)') '                     ROUTINE STATISTICS'
    write(6,'(a)') ''
    write(6,'(a)') '    routine                        calls           cpu time'
  do rn_i = 1, routines_nb
   rn_srt_i = routines_indexes_sort (rn_i)
   if (routines(rn_srt_i)%calls_nb >= 1) then
   routine_name = trim(routines(rn_srt_i)%name)
   write(6,'(a,i10,f)') routine_name, routines(rn_srt_i)%calls_nb, routines(rn_srt_i)%cpu_duration
    cpu_duration_total = cpu_duration_total +  routines(rn_srt_i)%cpu_duration
   endif
  enddo
   write(6,'(a)') '------------------------------------------------------------------'
   write(6,'(a,f)') 'TOTAL                                   ',cpu_duration_total
   write(6,'(a)') ''

 end subroutine routines_statistics

end module routines_mod
