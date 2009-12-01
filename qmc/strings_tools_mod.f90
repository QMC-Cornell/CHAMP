module strings_tools_mod

  use basic_tools_mod
  use constants_mod

!===============================================================
  interface string
!---------------------------------------------------------------
   module procedure string_integer
  end interface string

!===============================================================
  interface string2
!---------------------------------------------------------------
   module procedure string2_integer
  end interface string2

!===============================================================
  interface operator (+)
!---------------------------------------------------------------
    module procedure &
      string_string_cat, &
      integer_string_cat, &
      string_integer_cat, &
      double_string_cat, &
      string_double_cat
  end interface

  contains

! ==================================================================================
  function string_string_cat (string_1, string_2) result(string)
! ----------------------------------------------------------------------------------
! Description : concatenates two strings
!
! Created     : J. Toulouse, 21 Mar 2007
! ----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)  :: string_1, string_2

! output
  character(len=max_string_len) :: string


! begin
  string = ''
  string = trim(string_1)//trim(string_2)

  end function string_string_cat

! ==================================================================================
  function integer_string_cat (intg, string) result(string_cat)
! ----------------------------------------------------------------------------------
! Description : concatenates one integer and one string
!
! Created     : J. Toulouse, 18 Apr 2007
! ----------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: intg
  character(len=*), intent(in)  :: string

! output
  character(len=max_string_len) :: string_cat

! local
  character(len=max_string_len) intg_string

! begin
  string_cat = ''
  write (intg_string,'(i5)') intg
  string_cat = trim(intg_string)//trim(string)

  end function integer_string_cat

! ==================================================================================
  function string_integer_cat (string, intg) result(string_cat)
! ----------------------------------------------------------------------------------
! Description : concatenates one integer and one string
!
! Created     : J. Toulouse, 18 Apr 2007
! ----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)  :: string
  integer, intent(in) :: intg

! output
  character(len=max_string_len) :: string_cat

! local
  character(len=max_string_len) intg_string

! begin
  string_cat = ''
  write (intg_string,'(i5)') intg
  string_cat = trim(string)//trim(intg_string)

  end function string_integer_cat

! ==================================================================================
  function double_string_cat (real, string) result(string_cat)
! ----------------------------------------------------------------------------------
! Description : concatenates one real and one string
!
! Created     : J. Toulouse, 18 Apr 2007
! ----------------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in) :: real
  character(len=*), intent(in)  :: string

! output
  character(len=max_string_len) :: string_cat

! local
  character(len=max_string_len) real_string

! begin
  string_cat = ''
  write (real_string,'(es15.8)') real
  string_cat = trim(real_string)//trim(string)

  end function double_string_cat

! ==================================================================================
  function string_double_cat (string, real) result(string_cat)
! ----------------------------------------------------------------------------------
! Description : concatenates one real and one string
!
! Created     : J. Toulouse, 18 Apr 2007
! ----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)  :: string
  real(dp), intent(in) :: real

! output
  character(len=max_string_len) :: string_cat

! local
  character(len=max_string_len) real_string

! begin
  string_cat = ''
  write (real_string,'(es15.8)') real
  string_cat = trim(string)//trim(real_string)

  end function string_double_cat

! ==================================================================================
  function string_integer (intg) result(string)
! ----------------------------------------------------------------------------------
! Description : convert to string a integer
!
! Created     : J. Toulouse, 21 May 2007
! ----------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: intg

! output
  character(len=max_string_len) :: string


! begin
  string = ''
  write (string,'(i1)') intg
  string = trim(string)

  end function string_integer

! ==================================================================================
  function string2_integer (intg) result(string)
! ----------------------------------------------------------------------------------
! Description : convert to string a integer
!
! Created     : J. Toulouse, 21 May 2007
! ----------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: intg

! output
  character(len=max_string_len) :: string

! begin
  string = ''
  write (string,'(i3)') intg
  string = trim(string)

  end function string2_integer

! ==================================================================================
  function string_to_integer (string)
! ----------------------------------------------------------------------------------
! Description : convert a string into an integer
!
! Created     : J. Toulouse, 26 Jul 2007
! ----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)  :: string

! output
  integer string_to_integer


! begin
  call cnvint (string, string_to_integer)

  end function string_to_integer

! ==================================================================================
  function is_string_integer (string)
! ----------------------------------------------------------------------------------
! Description : is string an integer?
!
! Created     : J. Toulouse, 03 Mar 2009
! ----------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)  :: string

! output
  logical is_string_integer

! local
  character(len=1) c
  integer i 

! begin
  is_string_integer = .true.

  do i = 1, len(trim(string))
   c = string(i:i)
   if (c /= '0' .and. c /= '1' .and. c /= '2' .and. &
       c /= '3' .and. c /= '4' .and. c /= '5' .and. &
       c /= '6' .and. c /= '7' .and. c /= '8' .and. &
       c /= '9') then
      is_string_integer = .false.
      exit
    endif
  enddo

  end function is_string_integer

end module strings_tools_mod
