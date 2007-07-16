module print_mod

  use basic_tools_mod
  use constants_mod
  use variables_mod
  use objects_mod
  use parser_tools_mod

  integer                      :: objects_print_at_each_block_nb = 0
  integer, allocatable         :: objects_print_at_each_block_index (:)

  contains

!===========================================================================
  subroutine print_menu
!---------------------------------------------------------------------------
! Description : menu for printing objects
!
! Created     : J. Toulouse, 29 Nov 2005
! Modified    : J. Toulouse, 15 Nov 2006, print at each block
!---------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'print_menu'
  integer objects_to_print_now_nb
  character(len=max_string_len), allocatable :: objects_to_print_now (:)
  integer objects_to_print_block_nb
  character(len=max_string_len), allocatable :: objects_to_print_block (:)
  integer obj_i
  
! begin
  objects_to_print_now_nb = 0
  objects_to_print_block_nb = 0

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))

  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for print menu:'
   write(6,'(a)') 'print'
   write(6,'(a)') '  now ... end: list of objects to print now'
   write(6,'(a)') '  block ... end: list of objects to print at each block'
   write(6,'(a)') 'end'
   write(6,*)


  case ('now')
   call get_next_value_list ('objects_to_print_now', objects_to_print_now, objects_to_print_now_nb)

!  provide and print objects now
   do obj_i = 1, objects_to_print_now_nb
    call object_provide (objects_to_print_now (obj_i))
    call object_write (objects_to_print_now (obj_i))
   enddo

  case ('block')
   call get_next_value_list ('objects_to_print_block', objects_to_print_block, objects_to_print_block_nb)

!  request printing of objects at each block
   do obj_i = 1, objects_to_print_block_nb
    call object_print_at_each_block_request (objects_to_print_block (obj_i))
   enddo

  case ('end')
   exit

  case default
   write(6,'(3a)') trim(lhere),': unknown keyword = ',trim(word)
   call die (lhere)
  end select

  enddo ! end loop over menu lines

  end subroutine print_menu

! ===================================================================================
  subroutine object_print_at_each_block_request (object_name)
! ----------------------------------------------------------------------------------
! Description   : request for printing an object at each block
!
! Created       : J. Toulouse, 15 Nov 2006
! -----------------------------------------------------------------------------------
  implicit none

! input 
  character(len=*), intent(in) :: object_name
  
! local
  character(len=max_string_len_rout), save :: lhere = 'object_print_at_each_block_request'
  integer object_ind
  

! begin
  call object_add_once_and_index (object_name, object_ind)

  call append_once (objects_print_at_each_block_index, object_ind)
  objects_print_at_each_block_nb = size(objects_print_at_each_block_index)

  write(6,'(4a)') trim(lhere),': object ', trim(object_name),' will be printed at each block.'

 end subroutine object_print_at_each_block_request

! ===================================================================================
  subroutine objects_print_at_each_block
! -----------------------------------------------------------------------------------
! Description   : print requested objects at each block
!
! Created       : J. Toulouse, 15 Nov 2006
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'objects_print_at_each_block'
  integer obj_i
  
! begin

  do obj_i = 1, objects_print_at_each_block_nb
    call object_provide_by_index (objects_print_at_each_block_index (obj_i))
    call object_write_by_index (objects_print_at_each_block_index (obj_i))
  enddo

 end subroutine objects_print_at_each_block

! ===================================================================================
  subroutine reinit_objects_print_at_each_block
! -----------------------------------------------------------------------------------
! Description   : reinitialize array of requested objects to be printed at each block
!
! Created       : J. Toulouse, 15 Nov 2006
! -----------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'reinit_objects_print_at_each_block'

! begin
  objects_print_at_each_block_nb = 0
  call release ('objects_print_at_each_block_index', objects_print_at_each_block_index)

 end subroutine reinit_objects_print_at_each_block

end module print_mod
