module files_tools_mod

  use basic_tools_mod
  use strings_tools_mod
  use constants_mod
  use objects_mod

! Declaration of global variables and default values
  integer                                     :: files_nb = 0
  character (len=max_string_len), allocatable :: file_names (:)
  integer, allocatable                        :: file_units (:)

  contains

!===========================================================================
  subroutine catalog_file (filename, file_unit)
!---------------------------------------------------------------------------
! Description : catalog a file and attribute it an unique file unit
! Description : if file_unit <= 0
!
! Created     : J. Toulouse, 21 Mar 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: filename

! input/output
  integer, intent(inout) :: file_unit

! local
  character (len=max_string_len_rout), save :: lhere = 'catalog_file'
  integer file_i

! begin

! check if file name is already catalogued
  do file_i = 1, files_nb
   if (trim(filename) == trim(file_names (file_i)) ) then
!    if file_unit > 0, check that file_unit agrees with the previously attributed one
     if (file_unit > 0 .and. file_unit /= file_units (file_i)) then
      write(6,'(4a,i3)') trim(lhere),': the file ',trim(filename),' is asked to be catalgued with file unit ',file_unit
      write(6,'(3a,i3)') trim(lhere),': but it has already been catalogued with file unit ',file_units (file_i)
      call die (lhere, 'a file is catalogued twice with different file units')
     endif
     file_unit = file_units (file_i)
     return
   endif
  enddo

!  check if file unit is already catalogued
   if (file_unit > 0) then
    do file_i = 1, files_nb
     if (file_unit == file_units (file_i) .and. trim(filename) /= trim(file_names (file_i))) then
      write(6,'(4a,i3)') trim(lhere),': the file ',trim(filename),' is asked to be catalgued with file unit ',file_unit
      write(6,'(3a,i3)') trim(lhere),': but the file unit is already associated to the file',trim(file_names (file_i))
      call die (lhere, 'two different files are catalogued with the same file unit')
     endif
    enddo
   endif

! add new file name
  files_nb = files_nb + 1
  call object_alloc ('file_names', file_names, files_nb)
  file_names (files_nb) = filename

! if needed, attribute a unique file unit to the new file (file units start at 40)
  if (file_unit <= 0) then
   file_unit = max_element (file_units) + 1
   file_unit = max(40,file_unit)
   if (file_unit > 99)  then
      call die (lhere, ' largest possible file unit 99 reached')
   endif
  endif
  call object_alloc ('file_units', file_units, files_nb)
  file_units (files_nb) = file_unit

!  write(6,'(3a,i3)') 'Catalog new file ',trim(filename),' with file unit ',file_unit

  end subroutine catalog_file

!===========================================================================
  subroutine open_file_or_die (filename, file_unit)
!---------------------------------------------------------------------------
! Description : opens file and attributes an unique unit if file_unit = 0
!
! Created     : J. Toulouse, 21 Mar 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: filename

! input/output
  integer, intent(inout) :: file_unit

! local
  character (len=max_string_len_rout), save :: lhere = 'open_file_or_die'
  integer iostat

! begin

! catalog file
  call catalog_file (filename, file_unit)

! checkings
!  call file_exist_or_die (filename)
  call file_not_opened_or_die (filename)
  call require (lhere, 'file_unit > 0', file_unit > 0) !fp

! open file
  open (file=trim(filename), unit=file_unit, iostat=iostat)
  if (iostat /= 0) then
   call die (lhere, 'error during opening file >'+trim(filename)+'<')
  endif

  end subroutine open_file_or_die

!===========================================================================
  subroutine file_not_opened_or_die (filename)
!---------------------------------------------------------------------------
! Description : check if file is opened, die otherwise
!
! Created     : J. Toulouse, 21 Mar 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: filename

! local
  character(len=max_string_len_rout), save :: lhere = 'file_not_opened_or_die'
  logical :: opened

! begin

  inquire(file=trim(filename), opened=opened)

  if (opened) then
   call die (lhere, 'try to open the file '+trim(filename)+' that is already opened.')
  endif

  end subroutine file_not_opened_or_die

!===========================================================================
  function file_exist (filename) result(exist)
!---------------------------------------------------------------------------
! Description : check is file exists
!
! Created     : J. Toulouse, 22 Mar 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: filename

! local
  character(len=max_string_len_rout), save :: lhere = 'file_exist'
  logical :: exist

! begin
  inquire(file=trim(filename), exist=exist)

  end function file_exist

!===========================================================================
  subroutine file_exist_or_die (filename)
!---------------------------------------------------------------------------
! Description : check is file exists, die otherwise
!
! Created     : J. Toulouse, 21 Mar 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character (len=*), intent(in) :: filename

! local
  character(len=max_string_len_rout), save :: lhere = 'file_exist_or_die'
  logical :: exist

! begin
  if (.not. file_exist(filename)) then
   call die (lhere, 'file >'+trim(filename)+'< does not exist.')
  endif

  end subroutine file_exist_or_die


!===========================================================================
  function lines_number_in_file (file_unit) result(lines_number)
!--------------------------------------------------------------------------
! Description : returns the number of lines in a open file of unit file_unit
!
! Created     : J. Toulouse, 21 Mar 2007
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in)  :: file_unit

! output
  integer              :: lines_number

! local
  character(len=max_string_len_rout), save :: lhere = 'lines_number_in_file'
  character(len=max_string_len) record
  integer iostat

! initilization
  iostat = 0
  lines_number = 0

! rewind file
  rewind file_unit

! loop over file's lines
  do

!  read current line
   read(file_unit,'(a)',iostat=iostat) record

   if (iostat < 0) exit

   lines_number = lines_number + 1

  enddo

! rewind file
  rewind file_unit

  write(6,'(2a,i2,a,i10,a)') trim(lhere),': the open file of unit ', file_unit,' has ',lines_number,' lines.'

  end function lines_number_in_file

end module files_tools_mod
