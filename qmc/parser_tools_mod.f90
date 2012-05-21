module parser_tools_mod

  use basic_tools_mod
  use constants_mod
  use objects_mod

  character(len=50000)                      :: current_line = ''
  integer                                   :: position_in_current_line = 0
  integer                                   :: unit_input = 5
  character(len=max_string_len)             :: word
  logical                                   :: l_echo = .false.

!---------------------------------------------------------------------------
  interface get_next_value
!---------------------------------------------------------------------------
  module procedure                  &
            get_next_value_string,  &
            get_next_value_integer, &
            get_next_value_double,    &
            get_next_value_logical
  end interface

!---------------------------------------------------------------------------
  interface read_next_line_list
!---------------------------------------------------------------------------
   module procedure               &
            read_next_line_list_integer, &
            read_next_line_list_double
  end interface

!---------------------------------------------------------------------------
  interface get_next_value_list_noalloc
!---------------------------------------------------------------------------
   module procedure               &
            get_next_value_list_noalloc_integer, &
            get_next_value_list_noalloc_double
  end interface

!---------------------------------------------------------------------------
  interface get_next_value_list
!---------------------------------------------------------------------------
   module procedure               &
            get_next_value_list_string,  &
            get_next_value_list_integer, &
            get_next_value_list_double
  end interface

!---------------------------------------------------------------------------
  interface get_next_value_list_object
!---------------------------------------------------------------------------
   module procedure               &
            get_next_value_list_object_string,  &
            get_next_value_list_object_integer, &
            get_next_value_list_object_double
  end interface

  contains

!---------------------------------------------------------------------------
  subroutine read_next_line (iostat)
!---------------------------------------------------------------------------
! Description : read next line in input file
!
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! output
  integer, intent(out) :: iostat

! local

! begin
  iostat = 0
  position_in_current_line = 0


! read current line
  read(unit_input,'(a)',iostat=iostat) current_line

! no next line found
  if(iostat < 0) return

  if (l_echo) then
   write(6,'(2a)') 'input> ',trim(current_line)
  endif

! convert to lower case
  call upplow (current_line)

!  write(6,*)trim(lhere),': next line found'
!  write(6,*)trim(lhere),': line >',trim(current_line),'<'

  end subroutine read_next_line

!===========================================================================
  subroutine get_next_command (word)
!---------------------------------------------------------------------------
! Description : find next command in input file
!
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! output
  character(len=max_string_len), intent(out) :: word

! local
  integer iostat
  integer current_line_length
  character (len=1) current_char
  integer first_char_in_word_index, last_char_in_word_index
  logical first_char_in_word_found, last_char_in_word_found

! begin

! added to handle the case where there is no 'end' after the optimization section
  call read_next_line (iostat)
  if (iostat < 0) then
   word = 'exit'
   return
  endif

! find first character of word
  first_char_in_word_found = .false.
  do

    current_line_length = len(trim(current_line))

!   loop over characters in current line
    do while (position_in_current_line < current_line_length)

     position_in_current_line = position_in_current_line + 1

     current_char = current_line(position_in_current_line:position_in_current_line)

!    ignore blank characters
     if(current_char == ' ')  cycle

!    ignore remaining of line after ! (comments)
     if(current_char == '!')  exit

     first_char_in_word_index = position_in_current_line
     first_char_in_word_found = .true.
     exit

    enddo

    if (first_char_in_word_found) exit

    call read_next_line (iostat)
    if (iostat < 0) then
     word = 'exit'
     return
    endif

  enddo

! find last character of word
  last_char_in_word_found = .false.

!  loop over characters in current line
    do while (position_in_current_line < current_line_length)

     position_in_current_line = position_in_current_line + 1

     current_char = current_line(position_in_current_line:position_in_current_line)

     if (position_in_current_line == current_line_length) then
      last_char_in_word_index = position_in_current_line
      last_char_in_word_found = .true.
      exit
     endif

     if (current_char == ' ' .or. current_char == '!') then
      last_char_in_word_index = position_in_current_line - 1
      last_char_in_word_found = .true.
      exit
     endif

    enddo

!  word
   word = current_line(first_char_in_word_index:last_char_in_word_index)

  end subroutine get_next_command

!===========================================================================
  subroutine get_next_word (word)
!---------------------------------------------------------------------------
! Description : find next word in input file
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! output
  character(len=max_string_len), intent(out) :: word

! local
  character(len=max_string_len_rout), save :: lhere = 'get_next_word'
  integer iostat
  integer current_line_length
  character(len=1) current_char
  integer first_char_in_word_index, last_char_in_word_index
  logical first_char_in_word_found, last_char_in_word_found

! begin

! find first character of word
  first_char_in_word_found = .false.
  do

    current_line_length = len(trim(current_line))

!   loop over characters in current line
    do while (position_in_current_line < current_line_length)

     position_in_current_line = position_in_current_line + 1

     current_char = current_line(position_in_current_line:position_in_current_line)

!     write(6,*) trim(here),': current_char=',current_char

!    ignore blank characters
     if(current_char == ' ')  cycle

!    ignore '='
     if(current_char == '=')  cycle

!    ignore remaining of line after ! (comments)
     if(current_char == '!')  exit

     first_char_in_word_index = position_in_current_line
     first_char_in_word_found = .true.
     exit

    enddo

    if (first_char_in_word_found) exit

    call read_next_line (iostat)
    if(iostat < 0) then
     call die(lhere,'unexpected end of file')
    endif

  enddo

! find last character of word
  last_char_in_word_found = .false.
  last_char_in_word_index = first_char_in_word_index

!  loop over characters in current line
    do while (position_in_current_line < current_line_length)

     position_in_current_line = position_in_current_line + 1

     current_char = current_line(position_in_current_line:position_in_current_line)

     if(position_in_current_line == current_line_length) then
      last_char_in_word_index = position_in_current_line
      last_char_in_word_found = .true.
      exit
     endif

     if(current_char == ' ' .or. current_char == '=' .or. current_char == '!') then
      last_char_in_word_index = position_in_current_line - 1
      last_char_in_word_found = .true.
      exit
     endif

    enddo

!  word
   word = current_line(first_char_in_word_index:last_char_in_word_index)

!  write(6,*) trim(lhere),': next word found >',trim(word),'<'

  end subroutine get_next_word

!===========================================================================
  function words_number_in_string (string)
!---------------------------------------------------------------------------
! Description : returns the number of words in a string
!
! Created     : J. Toulouse, 08 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: string

! output
  integer :: words_number_in_string

! local
!  character(len=max_string_len_rout), save :: lhere = 'words_number_in_string'
  integer string_length, position_in_string
  character(len=1) current_char
  integer first_char_in_word_index, last_char_in_word_index
  logical first_char_in_word_found, last_char_in_word_found

! begin
  words_number_in_string = 0
  position_in_string = 0

! loop over succesive words until end of string
  do

!   find first character of word
    first_char_in_word_found = .false.
    do
    
      string_length = len(trim(string))
    
!     loop over characters in string
      do while (position_in_string < string_length)
    
       position_in_string = position_in_string + 1
    
       current_char = string(position_in_string:position_in_string)
    
!       write(6,*) trim(here),': current_char=',current_char
    
!      ignore blank characters
       if(current_char == ' ')  cycle
    
!      ignore '='
       if(current_char == '=')  cycle
    
!      ignore remaining of line after ! (comments)
       if(current_char == '!')  exit
    
       first_char_in_word_index = position_in_string
       first_char_in_word_found = .true.
       exit
    
      enddo
    
      if (first_char_in_word_found) exit
    
!     end of string is reached
      return
    
    enddo
    
!   find last character of word
    last_char_in_word_found = .false.
    last_char_in_word_index = first_char_in_word_index
    
!    loop over characters in current line
      do while (position_in_string < string_length)
    
       position_in_string = position_in_string + 1
    
       current_char = string(position_in_string:position_in_string)
    
       if(position_in_string == string_length) then
        last_char_in_word_index = position_in_string
        last_char_in_word_found = .true.
        exit
       endif
    
       if(current_char == ' ' .or. current_char == '=' .or. current_char == '!') then
        last_char_in_word_index = position_in_string - 1
        last_char_in_word_found = .true.
        exit
       endif
    
      enddo
    
!    word
!     word = string(first_char_in_word_index:last_char_in_word_index)
!    write(6,*) trim(lhere),': next word found >',trim(word),'<'
     words_number_in_string = words_number_in_string + 1
     
  enddo ! end loop over words

  end function words_number_in_string

!===========================================================================
  subroutine read_next_word_in_line (word)
!---------------------------------------------------------------------------
! Description : read next word in the line
! Description : returns 'end_of_line' when end of the line is reached
!
! Created     : J. Toulouse, 07 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! output
  character(len=max_string_len), intent(out) :: word

! local
!  character(len=max_string_len_rout), save :: lhere = 'read_next_word_in_line'
  integer current_line_length
  character(len=1) current_char
  integer first_char_in_word_index, last_char_in_word_index
  logical first_char_in_word_found, last_char_in_word_found

! begin

! find first character of word
  first_char_in_word_found = .false.
  do

    current_line_length = len(trim(current_line))

!   loop over characters in current line
    do while (position_in_current_line < current_line_length)

     position_in_current_line = position_in_current_line + 1

     current_char = current_line(position_in_current_line:position_in_current_line)

!     write(6,*) trim(here),': current_char=',current_char

!    ignore blank characters
     if(current_char == ' ')  cycle

!    ignore '='
     if(current_char == '=')  cycle

!    ignore remaining of line after ! (comments)
     if(current_char == '!')  exit

     first_char_in_word_index = position_in_current_line
     first_char_in_word_found = .true.
     exit

    enddo

    if (first_char_in_word_found) exit

!   end of line is reached
    word = 'end_of_line'
    return

  enddo

! find last character of word
  last_char_in_word_found = .false.
  last_char_in_word_index = first_char_in_word_index

!  loop over characters in current line
    do while (position_in_current_line < current_line_length)

     position_in_current_line = position_in_current_line + 1

     current_char = current_line(position_in_current_line:position_in_current_line)

     if(position_in_current_line == current_line_length) then
      last_char_in_word_index = position_in_current_line
      last_char_in_word_found = .true.
      exit
     endif

     if(current_char == ' ' .or. current_char == '=' .or. current_char == '!') then
      last_char_in_word_index = position_in_current_line - 1
      last_char_in_word_found = .true.
      exit
     endif

    enddo

!  word
   word = current_line(first_char_in_word_index:last_char_in_word_index)

!  write(6,*) trim(lhere),': next word found >',trim(word),'<'

  end subroutine read_next_word_in_line

!===========================================================================
  subroutine get_next_value_string (value)
!---------------------------------------------------------------------------
! Description : find next value in line
!
! Created     : J. Toulouse, 5 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! output
  character(len=*), intent(out) :: value

! local
  character(len=max_string_len) value_string

! begin
  call get_next_word (value_string)
  value = value_string

  end subroutine get_next_value_string

!===========================================================================
  subroutine get_next_value_integer (value)
!---------------------------------------------------------------------------
! Description : find next value in line
!
! Created     : J. Toulouse, 5 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! output
  integer, intent(out) :: value

! local
  character(len=max_string_len) value_string

! begin
  call get_next_word (value_string)
  call cnvint(value_string, value)

  end subroutine get_next_value_integer

!---------------------------------------------------------------------------
  subroutine get_next_value_double (value)
!---------------------------------------------------------------------------
! Description : find next value in line
!
! Created     : J. Toulouse, 5 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! output
  real(dp), intent(out) :: value

! local
  character(len=max_string_len) value_string

! begin
  call get_next_word (value_string)
  call cnvdbl(value_string, value)

  end subroutine get_next_value_double

!---------------------------------------------------------------------------
  subroutine get_next_value_logical (value)
!---------------------------------------------------------------------------
! Description : find next value in line
!
! Created     : J. Toulouse, 22 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! output
  logical, intent(out) :: value

! local
  character(len=max_string_len_rout), save :: lhere = 'get_next_value_logical'
  character(len=max_string_len) value_string

! begin
  call get_next_word (value_string)

  if (trim(value_string) == 'true' .or. trim(value_string) == '.true.' ) then
    value = .true.
  elseif (trim(value_string) == 'false' .or. trim(value_string) == '.false.' ) then
    value = .false.
  else
    call die (lhere, 'value >'+trim(value_string)+'< must be "true" or "false".')
  endif

  end subroutine get_next_value_logical

!===========================================================================
  subroutine get_next_value_list_string (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : find next list of values (not a catalogued object)
!
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent (in)              :: value_list_name

! output
# if defined (PATHSCALE)
   character(len=*), intent(out) :: value_list (max_string_array_len) ! for pathscale compiler
# else
   character(len=*), allocatable, intent(out) :: value_list (:)
# endif
  integer, intent(out)                       :: value_list_nb

! local
  character(len=max_string_len_rout), save :: lhere = 'get_next_value_list_string'
  character(len=max_string_len) value_string

! begin
  value_string = ''
  value_list_nb = 0

! loop over succesive words until 'end' is found
  do

    call get_next_word (value_string)

    if (value_string == 'end') exit

    value_list_nb = value_list_nb + 1

# if !defined (PATHSCALE)
   call alloc (value_list_name, value_list, value_list_nb) ! commented out for pathscale compiler
# endif
    value_list (value_list_nb) = trim(value_string)
!    write(6,*) 'value_string=',value_string

  enddo ! end loop

  if (value_list_nb == 0 ) then
    call die (lhere, 'no values found.')
  endif

  end subroutine get_next_value_list_string

!===========================================================================
  subroutine get_next_value_list_integer (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : find next list of values (not a catalogued object)
!
! Created     : J. Toulouse, 25 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: value_list_name

! output
# if defined (PATHSCALE)
   integer, intent(out) :: value_list (max_int_array_len) ! for pathscale compiler
# else
   integer, allocatable, intent(out) :: value_list (:)
# endif
  integer, intent(out)              :: value_list_nb

! local
# if defined (PATHSCALE)
   character(len=max_string_len):: value_list_string (max_string_array_len) ! for pathscale compiler
# else
   character(len=max_string_len), allocatable :: value_list_string (:)
# endif
  integer i

! begin
  call get_next_value_list_string (value_list_name, value_list_string, value_list_nb)

# if !defined (PATHSCALE)
  call alloc (value_list_name, value_list, value_list_nb)
# endif

  do i = 1, value_list_nb
    call cnvint (value_list_string(i), value_list(i))
  enddo

  end subroutine get_next_value_list_integer

!===========================================================================
  subroutine get_next_value_list_double (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : find next list of values (not a catalogued object)
!
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: value_list_name

! output
# if defined (PATHSCALE)
   real(dp), intent(out)              :: value_list (max_double_array_len) ! for pathscale compiler
# else
   real(dp), allocatable, intent(out) :: value_list (:)
# endif
  integer, intent(out)              :: value_list_nb

! local
# if defined (PATHSCALE)
   character(len=max_string_len) :: value_list_string (max_string_array_len) ! for pathscale compiler
# else
   character(len=max_string_len), allocatable :: value_list_string (:)
# endif
  integer i

! begin
  call get_next_value_list_string (value_list_name, value_list_string, value_list_nb)

# if !defined (PATHSCALE)
  call alloc (value_list_name, value_list, value_list_nb)
# endif

  do i = 1, value_list_nb
    call cnvdbl (value_list_string(i), value_list(i))
  enddo

  end subroutine get_next_value_list_double

!===========================================================================
  subroutine read_next_line_list_string (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : read next line and store it in an array
!
! Created     : J. Toulouse, 07 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: value_list_name

! output
  character (len=*), allocatable, intent(out)  :: value_list (:)
  integer, intent(out)  :: value_list_nb

! local
  character(len=max_string_len_rout), save :: lhere = 'read_next_line_list_string'
  character(len=max_string_len) value_string
  integer iostat

! begin
  value_string = ''
  value_list_nb = 0

! read line
  call read_next_line (iostat)
  if(iostat < 0) then
    call die (lhere, 'error when reading line') 
  endif

! loop over succesive words until end of line
  do

    call read_next_word_in_line (value_string)

    if (value_string == 'end_of_line') exit

    value_list_nb = value_list_nb + 1

    call alloc (value_list_name, value_list, value_list_nb)
    value_list (value_list_nb) = trim(value_string)
!    write(6,*) 'value_string=',value_string

  enddo ! end loop

  if (value_list_nb == 0 ) then
    call die (lhere, 'no values found.')
  endif

  end subroutine read_next_line_list_string

!===========================================================================
  subroutine read_next_line_list_integer (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : read next line and store it in an array
!
! Created     : J. Toulouse, 07 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: value_list_name

! output
  integer, allocatable, intent(out)  :: value_list (:)
  integer, intent(out)  :: value_list_nb

! local
  character(len=max_string_len), allocatable :: value_list_string (:)
  integer i

! begin
  call read_next_line_list_string (value_list_name, value_list_string, value_list_nb)

  call alloc (value_list_name, value_list, value_list_nb)

  do i = 1, value_list_nb
    call cnvint (value_list_string(i), value_list(i))
  enddo

  end subroutine read_next_line_list_integer

!===========================================================================
  subroutine read_next_line_list_double (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : read next line and store it in an array
!
! Created     : J. Toulouse, 07 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: value_list_name

! output
  real(dp), allocatable, intent(out) :: value_list (:)
  integer, intent(out) :: value_list_nb

! local
  character(len=max_string_len), allocatable :: value_list_string (:)
  integer i

! begin
  call read_next_line_list_string (value_list_name, value_list_string, value_list_nb)

  call alloc (value_list_name, value_list, value_list_nb)

  do i = 1, value_list_nb
    call cnvdbl (value_list_string(i), value_list(i))
  enddo

  end subroutine read_next_line_list_double

!===========================================================================
  subroutine get_next_value_list_noalloc_integer (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : find next list of values (not a catalogued object, not allocated)
!
! Created     : J. Toulouse, 07 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: value_list_name

! output
   integer, intent(out) :: value_list (:)
  integer, intent(out)  :: value_list_nb

! local
  character(len=max_string_len), allocatable :: value_list_string (:)
  integer i

! begin
  call get_next_value_list_string (value_list_name, value_list_string, value_list_nb)

  do i = 1, value_list_nb
    call cnvint (value_list_string(i), value_list(i))
  enddo

  end subroutine get_next_value_list_noalloc_integer

!===========================================================================
  subroutine get_next_value_list_noalloc_double (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : find next list of values (not a catalogued object, not allocated)
!
! Created     : J. Toulouse, 07 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: value_list_name

! output
  real(dp), intent(out) :: value_list (:)
  integer, intent(out)  :: value_list_nb

! local
  character(len=max_string_len), allocatable :: value_list_string (:)
  integer i

! begin
  call get_next_value_list_string (value_list_name, value_list_string, value_list_nb)

  do i = 1, value_list_nb
    call cnvdbl (value_list_string(i), value_list(i))
  enddo

  end subroutine get_next_value_list_noalloc_double

!===========================================================================
  subroutine get_next_value_list_object_string (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : find next list of values and allocate it as a catalogued object
!
! Created     : J. Toulouse, 10 Mar 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent (in)              :: value_list_name

! output
# if defined (PATHSCALE)
   character(len=*), intent(out) :: value_list (max_string_array_len) ! for pathscale compiler
# else
   character(len=*), allocatable, intent(out) :: value_list (:)
# endif
  integer, intent(out)                       :: value_list_nb

! local
  character(len=max_string_len_rout), save :: lhere = 'get_next_value_list_object_string'
  character(len=max_string_len) value_string

! begin
  value_string = ''
  value_list_nb = 0

! loop over succesive words until 'end' is found
  do

    call get_next_word (value_string)

    if (value_string == 'end') exit

    value_list_nb = value_list_nb + 1

# if !defined (PATHSCALE)
   call object_alloc (value_list_name, value_list, value_list_nb) ! commented out for pathscale compiler
# endif
    value_list (value_list_nb) = trim(value_string)

  enddo ! end loop

  if (value_list_nb == 0 ) then
    call die (lhere, 'no values found.')
  endif

  call object_modified (value_list_name)

  end subroutine get_next_value_list_object_string

!===========================================================================
  subroutine get_next_value_list_object_integer (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : find next list of values and allocate it as a catalogued object
!
! Created     : J. Toulouse, 10 Mar 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: value_list_name

! output
# if defined (PATHSCALE)
   integer, intent(out) :: value_list (max_int_array_len) ! for pathscale compiler
# else
   integer, allocatable, intent(out) :: value_list (:)
# endif
  integer, intent(out)              :: value_list_nb

! local
# if defined (PATHSCALE)
   character(len=max_string_len):: value_list_string (max_string_array_len) ! for pathscale compiler
# else
   character(len=max_string_len), allocatable :: value_list_string (:)
# endif
  integer i

! begin
  call get_next_value_list_string (value_list_name, value_list_string, value_list_nb)

# if !defined (PATHSCALE)
  call object_alloc (value_list_name, value_list, value_list_nb)
# endif

  do i = 1, value_list_nb
    call cnvint (value_list_string(i), value_list(i))
  enddo

  call object_modified (value_list_name)

  end subroutine get_next_value_list_object_integer

!===========================================================================
  subroutine get_next_value_list_object_double (value_list_name, value_list, value_list_nb)
!---------------------------------------------------------------------------
! Description : find next list of values and allocate it as a catalogued object
!
! Created     : J. Toulouse, 10 Mar 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: value_list_name

! output
# if defined (PATHSCALE)
   real(dp), intent(out)              :: value_list (max_double_array_len) ! for pathscale compiler
# else
   real(dp), allocatable, intent(out) :: value_list (:)
# endif
  integer, intent(out)              :: value_list_nb

! local
# if defined (PATHSCALE)
   character(len=max_string_len) :: value_list_string (max_string_array_len) ! for pathscale compiler
# else
   character(len=max_string_len), allocatable :: value_list_string (:)
# endif
  integer i

! begin
  call get_next_value_list_string (value_list_name, value_list_string, value_list_nb)

# if !defined (PATHSCALE)
  call object_alloc (value_list_name, value_list, value_list_nb)
# endif

  do i = 1, value_list_nb
    call cnvdbl (value_list_string(i), value_list(i))
  enddo

  call object_modified (value_list_name)

  end subroutine get_next_value_list_object_double

!===========================================================================
  subroutine read_up_to_end
!---------------------------------------------------------------------------
! Description : read lines in input files up to next 'end' keyword or end of file
!
! Created     : J. Toulouse, 5 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len) string_temp !, lowstring
!  integer iostat
  logical l_echo_save

  l_echo_save = l_echo
  l_echo = .false.

  do

! read current line
!  read(5,'(a)',iostat=iostat) string

  call get_next_command (string_temp)

! no end found
!  if(iostat < 0) then
!   write(6,*) trim(lhere),': no "end" keyword found in input file'
!   call die(here, 'no "end" keyword found in input file')
!   exit
!  endif

! convert to lower case
!  lowstring = string
!  call upplow (lowstring)

  if(trim(string_temp) == 'end') then
   if (l_echo_save) then
    write(6,'(2a)') 'input> ',trim(string_temp)
   endif
   exit
  endif

  if(trim(string_temp) == 'exit') then
   exit
  endif

  enddo

  l_echo = l_echo_save

  end subroutine read_up_to_end

end module parser_tools_mod
