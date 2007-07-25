module write_mod

  use basic_tools_mod
  use constants_mod
  use variables_mod
  use objects_mod

  integer                      :: write_routines_nb = 0
  integer, allocatable         :: write_routines_index (:)

! Interfaces

  contains

! ===================================================================================
  subroutine routine_write (routine_name)
! -----------------------------------------------------------------------------------
! Description   : define writing routine to be called at each block
! Description   : store index of routine
!
! Created       : J. Toulouse, 04 Mar 2006
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! local
  character(len=max_string_len_rout), save :: lhere = 'routine_write'
  integer routine_ind
  integer rtn_i


! begin
! index of routine
  routine_ind = routine_index (routine_name)
  if (routine_ind == 0) then
     write(6,'(4a)') trim(lhere),': routine ', trim(routine_name), ' is not catalogued.'
    call die (lhere)
  endif

! if routine already defined as write routine, do nothing
  do rtn_i = 1, write_routines_nb
   if (routine_ind == write_routines_index (rtn_i) ) then
       return
   endif
  enddo

  write_routines_nb = write_routines_nb + 1
  call alloc ('write_routines_index', write_routines_index, write_routines_nb)
  write_routines_index (write_routines_nb) = routine_ind

  write(6,'(4a)') trim(lhere),': ', trim(routine_name),' defined as writing routine'

 end subroutine routine_write

! ===================================================================================
  subroutine writing_routines
! -----------------------------------------------------------------------------------
! Description   : call writing routines at each block
!
! Created       : J. Toulouse, 20 Oct 2005
! -----------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'writing_routines'
  integer rtn_i

! begin
# if defined (DEBUG)
   call routine_enter (lhere)
# endif

! for MPI: write only for main process (idtask=0)
  if (idtask /= 0) then
   return
  endif

  do rtn_i = 1, write_routines_nb
    call exe_by_address_0 (routines(write_routines_index(rtn_i))%address)
  enddo !rtn_i

# if defined (DEBUG)
   call routine_exit (lhere)
# endif

 end subroutine writing_routines

! ===================================================================================
  subroutine reinit_writing_routines
! -----------------------------------------------------------------------------------
! Description   : reinitialize array of writing routines
!
! Created       : J. Toulouse, 04 Mar 2006
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'reinit_writing_routines'

! begin

  write_routines_nb = 0
  call release ('write_routines_index', write_routines_index)

 end subroutine reinit_writing_routines

end module write_mod
