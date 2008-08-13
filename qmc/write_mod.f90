module write_mod

  use basic_tools_mod
  use constants_mod
  use variables_mod
  use objects_mod

  integer                      :: routines_write_block_nb = 0
  integer, allocatable         :: routines_write_block_index (:)
  integer                      :: routines_write_final_nb = 0
  integer, allocatable         :: routines_write_final_index (:)

  contains

! ===================================================================================
  subroutine routine_write_block_request (routine_name)
! -----------------------------------------------------------------------------------
! Description   : requesting a writing routine to be called at each block
! Description   : store index of routine
!
! Created       : J. Toulouse, 04 Mar 2006
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! local
  character(len=max_string_len_rout), save :: lhere = 'routine_write_block_request'
  integer routine_ind, rtn_i


! begin
! index of routine
  routine_ind = routine_index (routine_name)
  if (routine_ind == 0) then
    call die (lhere, 'routine >'+trim(routine_name)+'< is not catalogued.')
  endif

! if routine already defined as write routine, do nothing
  do rtn_i = 1, routines_write_block_nb
   if (routine_ind == routines_write_block_index (rtn_i) ) then
       return
   endif
  enddo

  routines_write_block_nb = routines_write_block_nb + 1
  call alloc ('routines_write_block_index', routines_write_block_index, routines_write_block_nb)
  routines_write_block_index (routines_write_block_nb) = routine_ind

!  write(6,'(4a)') trim(lhere),': ', trim(routine_name),' defined as writing routine'

 end subroutine routine_write_block_request

! ===================================================================================
  subroutine routine_write_final_request (routine_name)
! -----------------------------------------------------------------------------------
! Description   : requesting a writing routine to be called at the end of MC run
! Description   : store index of routine
!
! Created       : J. Toulouse, 23 Jul 2008
! -----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! local
  character(len=max_string_len_rout), save :: lhere = 'routine_write_final_request'
  integer routine_ind
  integer rtn_i


! begin
! index of routine
  routine_ind = routine_index (routine_name)
  if (routine_ind == 0) then
     write(6,'(4a)') trim(lhere),': routine ', trim(routine_name), ' is not catalogued.'
    call die (lhere)
  endif

! if routine already defined as writing routine, do nothing
  do rtn_i = 1, routines_write_final_nb
   if (routine_ind == routines_write_final_index (rtn_i)) then
       return
   endif
  enddo

  routines_write_final_nb = routines_write_final_nb + 1
  call alloc ('routines_write_final_index', routines_write_final_index, routines_write_final_nb)
  routines_write_final_index (routines_write_final_nb) = routine_ind

!  write(6,'(4a)') trim(lhere),': ', trim(routine_name),' defined as writing routine'

 end subroutine routine_write_final_request

! ===================================================================================
  subroutine routines_write_block
! -----------------------------------------------------------------------------------
! Description   : call writing routines at each block
!
! Created       : J. Toulouse, 20 Oct 2005
! -----------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'routines_write_block'
  integer rtn_i

! begin
# if defined (DEBUG)
   call routine_enter (lhere)
# endif

! for MPI: write only for main process (idtask=0)
  if (idtask /= 0) then
   return
  endif

  do rtn_i = 1, routines_write_block_nb
    call exe_by_address_0 (routines(routines_write_block_index(rtn_i))%address)
  enddo !rtn_i

# if defined (DEBUG)
   call routine_exit (lhere)
# endif

 end subroutine routines_write_block

! ===================================================================================
  subroutine routines_write_final
! -----------------------------------------------------------------------------------
! Description   : call writing routines at the end of a MC run
!
! Created       : J. Toulouse, 23 Jul 2008
! -----------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'routines_write_final'
  integer rtn_i

! begin
# if defined (DEBUG)
   call routine_enter (lhere)
# endif

! for MPI: write only for main process (idtask=0)
  if (idtask /= 0) then
   return
  endif

  do rtn_i = 1, routines_write_final_nb
    call exe_by_address_0 (routines(routines_write_final_index(rtn_i))%address)
  enddo !rtn_i

# if defined (DEBUG)
   call routine_exit (lhere)
# endif

 end subroutine routines_write_final

! ===================================================================================
  subroutine reinit_routines_write_block
! -----------------------------------------------------------------------------------
! Description   : reinitialize array of writing routines
!
! Created       : J. Toulouse, 04 Mar 2006
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'reinit_routines_write_block'

! begin
  routines_write_block_nb = 0
  call release ('routines_write_block_index', routines_write_block_index)

 end subroutine reinit_routines_write_block

! ===================================================================================
  subroutine reinit_routines_write_final
! -----------------------------------------------------------------------------------
! Description   : reinitialize array of writing routines
!
! Created       : J. Toulouse, 23 Jul 2008
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'reinit_routines_write_final'

! begin
  routines_write_final_nb = 0
  call release ('routines_write_final_index', routines_write_final_index)

 end subroutine reinit_routines_write_final

end module write_mod
