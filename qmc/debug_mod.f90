module debug_mod

  use all_tools_mod

  logical                                     :: l_track_NaN = .false.

  contains

!===========================================================================
  subroutine debug_menu
!---------------------------------------------------------------------------
! Description : menu for debug options
!
! Created     : J. Toulouse, 16 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'debug_menu'
  integer rtn_i, routine_ind
  integer                                  :: debug_routines_list_nb = 0
# if defined (PATHSCALE)
   character(len=max_string_len) :: debug_routines_list (max_string_array_len) ! for pathscale compiler
# else
   character(len=max_string_len), allocatable  :: debug_routines_list (:)
# endif

! begin

# if !defined (DEBUG)
    call die (lhere, 'the code must be compiled in DEBUG mode for debug options.')
# endif

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))

  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for debug menu:'
   write(6,'(a)') 'debug'
   write(6,'(a)') ' trace = [boolean] active tracing of nodes and routines (default=false)'
   write(6,'(a)') ' routines ... end : list of routines to debug'
   write(6,'(a)') ' print_all_objects = [bool] print all objects (default=false)'
!   write(6,'(a)') ' track_NaN = [boolean] track NaN (default=false)'
   write(6,'(a)') 'end'

  case ('trace')
   call get_next_value (l_trace)

  case ('routines')
   call get_next_value_list ('debug_routines_list', debug_routines_list, debug_routines_list_nb)

  case ('print_all_objects')
   call get_next_value (l_print_all_objects)

!  case ('track_NaN')
!   call get_next_value (l_track_NaN)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines

  do rtn_i = 1, debug_routines_list_nb
    routine_ind = routine_catalog_and_index (debug_routines_list(rtn_i))
    routines(routine_ind)%debug = .true.
  enddo

  end subroutine debug_menu

end module debug_mod
