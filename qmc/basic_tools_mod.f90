module basic_tools_mod

  use constants_mod
  use variables_mod
  use types_mod
  use mpi_mod

! Interfaces

!===============================================================
  interface die
!---------------------------------------------------------------
   module procedure die_message, &
                    die_nomessage, &
                    die_basic
  end interface die

!===============================================================
   interface alloc
!---------------------------------------------------------------
   module procedure alloc_integer_1, &
                    alloc_integer_2, &
                    alloc_integer_3, &
                    alloc_integer_row_1, &
                    alloc_double_1,  &
                    alloc_double_2,  &
                    alloc_double_3,  &
                    alloc_double_4,  &
                    alloc_double_5,  &
                    alloc_double_row_1, &
                    alloc_complex_1,  &
                    alloc_complex_2,  &
                    alloc_complex_3,  &
                    alloc_complex_4,  &
                    alloc_logical_1,  &
                    alloc_logical_2,  &
                    alloc_logical_3,  &
                    alloc_string_1,   &
                    alloc_string_row_1

   end interface alloc

!===============================================================
   interface alloc_range
!---------------------------------------------------------------
   module procedure alloc_range_double_1, &
                    alloc_range_double_2, &
                    alloc_range_double_3
   end interface alloc_range

!===============================================================
   interface release
!---------------------------------------------------------------
   module procedure release_integer_1, &
                    release_integer_2,  &
                    release_integer_3,  &
                    release_integer_row_1,  &
                    release_double_1,  &
                    release_double_2,  &
                    release_double_3,  &
                    release_double_4,  &
                    release_double_5,  &
                    release_double_row_1,  &
                    release_complex_1,  &
                    release_complex_2,  &
                    release_complex_3,  &
                    release_complex_4,  &
                    release_logical_1,  &
                    release_logical_2,  &
                    release_logical_3,  &
                    release_string_1,   &
                    release_string_row_1

   end interface release

!===============================================================
   interface flatten
!---------------------------------------------------------------
   module procedure flatten_integer_2   , &
                    flatten_double_2

   end interface flatten

!===============================================================
   interface unflatten
!---------------------------------------------------------------
   module procedure unflatten_integer_2   , &
                    unflatten_double_2

   end interface unflatten

!===============================================================
  interface copy
!---------------------------------------------------------------
   module procedure copy_integer_1, &
                    copy_double_1

  end interface copy

!===============================================================
  interface move
!---------------------------------------------------------------
   module procedure move_integer_1, &
                    move_double_1

  end interface move

!===============================================================
  interface append
!---------------------------------------------------------------
   module procedure append_integer_0_to_1, &
                    append_integer_1_to_1, &
                    append_double_0_to_1,  &
                    append_double_1_to_1, &
                    append_string_0_to_1

  end interface append

!===============================================================
  interface append_once
!---------------------------------------------------------------
   module procedure append_once_integer
!                    append_once_string ! commented out for pathscale compiler
  end interface append_once

!===============================================================
  interface last
!---------------------------------------------------------------
   module procedure last_integer
  end interface last

!===============================================================
  interface swap
!---------------------------------------------------------------
   module procedure swap_integer, &
                    swap_double,  &
                    swap_double1
  end interface swap

!===============================================================
  interface sort
!---------------------------------------------------------------
   module procedure sort_integer1, &
                    sort_integer1_integer1, &
                    sort_integer1_integer1_double1, &
                    sort_double1_double1_double2
  end interface sort

!===============================================================
  interface is_sorted
!---------------------------------------------------------------
   module procedure is_sorted_integer1
  end interface is_sorted

!===============================================================
  interface sort_first_dim
!---------------------------------------------------------------
   module procedure sort_first_dim_integer
  end interface sort_first_dim

!===============================================================
  interface remove_elt_in_array
!---------------------------------------------------------------
   module procedure remove_elt_in_array_integer, &
                    remove_elt_in_array_double_1, &
                    remove_elt_in_array_double_3

  end interface remove_elt_in_array

!===============================================================
  interface is_equal_or_die
!---------------------------------------------------------------
   module procedure is_equal_or_die_double_1, &
                    is_equal_or_die_double_2, &
                    is_equal_or_die_double_3
  end interface is_equal_or_die

!===============================================================
  interface array_is_zero
!---------------------------------------------------------------
   module procedure array_is_zero_double

  end interface array_is_zero

!===============================================================
  interface arrays_equal
!---------------------------------------------------------------
   module procedure arrays_equal_logical,  &
                    arrays_equal_integer,  &
                    arrays_equal_double_1, &
                    arrays_equal_double_2
  end interface arrays_equal

!===============================================================
  interface elt_in_array
!---------------------------------------------------------------
   module procedure elt_in_array_integer, &
                    elt_in_array_string

  end interface elt_in_array


!===============================================================
  interface write_array
!---------------------------------------------------------------
   module procedure write_array_double_1
  end interface write_array

!===============================================================
  interface mysize
!---------------------------------------------------------------
   module procedure mysize_integer_1, &
                    mysize_double_1, &
                    mysize_string_1
  end interface mysize

!===============================================================
  interface last_element
!---------------------------------------------------------------
   module procedure last_element_integer
  end interface last_element

!===============================================================
  interface max_element
!---------------------------------------------------------------
   module procedure max_element_integer
  end interface max_element

!===============================================================
  interface array_norm
!---------------------------------------------------------------
   module procedure array_norm_double_1
  end interface array_norm

  contains

!===========================================================================
  subroutine die_basic
!---------------------------------------------------------------------------
! Description : stop the program
!
! Created     : J. Toulouse, 15 Oct 2005
!---------------------------------------------------------------------------
  implicit none


  write(6,*)
  write(6,'(a)') 'The program died.'
  stop 'The program died.'

  end subroutine die_basic

!===========================================================================
  subroutine die_nomessage (routine_name)
!---------------------------------------------------------------------------
! Description : stop the program with no error message
!
! Created     : J. Toulouse, 15 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: routine_name

! local
!  integer ierr

  write(6,*)
!  write(6,'(a,i3)') 'idtask=',idtask
  write(6,'(3a)') 'The program died because of a fatal error in routine ',trim(routine_name),'.'

! MPI beg -------------------------------------------------------------------
!# if defined (MPI)
!   call mpi_finalize (ierr)
!   if (ierr /= 0) then
!    write (6,'(2a)') trim(lhere),': error in mpi_finalize'
!   endif
!# endif
! MPI end -------------------------------------------------------------------

  stop 'The program died.'

  end subroutine die_nomessage

!===========================================================================
  subroutine die_message (routine_name, message)
!---------------------------------------------------------------------------
! Description : stop the program with an error message
!
! Created     : J. Toulouse, 15 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent (in) :: routine_name, message

  write(6,*)
  write(6,'(3a)') trim(routine_name),': ', trim(message)
  write(6,'(3a)') 'The program died because of a fatal error in routine ',trim(routine_name),'.'
  stop 'The program died.'

  end subroutine die_message

!===========================================================================
!  function locate_string_in_array (array, string) result(result) ! commented out for pathscale compiler
!---------------------------------------------------------------------------
! Description : return index of string in array
! Created     : J. Toulouse, 15 Oct 2005
!---------------------------------------------------------------------------
!  implicit none
!
!! input
!  character (len=*), allocatable :: array(:)
!  character (len=*) string
!
!! output
!  integer result
!
!! local
!  character(len=max_string_len_rout) lhere
!  integer i
!
!! begin
!  lhere = 'locate_string_in_array'
!
!  do i = 1, size(array)
!    if (string == array(i) ) then
!       result = i
!       return
!    endif
!   enddo
!
!  write(6,*) trim(lhere),': string not found in array'
!  call die(lhere)
!
!  return
!  end function locate_string_in_array

!===========================================================================
  subroutine alloc_string_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : Allocation or resize a array of string of 1 dimension
!
! Created     : J. Toulouse, 05 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1

! output
!  character(len=max_string_len), allocatable  :: object (:) ! comment out for lbasis
  character(len=*), allocatable  :: object (:) ! is this syntax ok?

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_string_1'
  integer i, object_dim, dim_min, all_err
  character(len=max_string_len), allocatable  :: object_temp (:)

! begin
! allocate object if not already allocated
  if(.not. allocated(object)) then
   allocate (object(dim1), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die(lhere,'allocation failed')
   endif

!  initialize array
   object(:) = ''

   else

! resize object if already allocated with different dimension
   object_dim = size(object)

   if(object_dim /= dim1) then

    dim_min =  min(object_dim, dim1)

    allocate (object_temp (dim_min))

    do i = 1, dim_min
     object_temp(i) = object(i)
    enddo

    call release (object_name, object)

    allocate (object (dim1))
    object(:) = ''

    do i = 1, dim_min
     object(i) = object_temp(i)
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_string_1

!===========================================================================
  subroutine alloc_string_row_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 11 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1

! output
  type (type_string_row), allocatable, intent(out)   :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_string_row_1'
  integer i, object_dim, dim_min, all_err
  type (type_string_row), allocatable          :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die (lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim = size(object)

   if(object_dim /= dim1) then

    dim_min =  min(object_dim, dim1)

    allocate (object_temp (dim_min))

    do i = 1, dim_min
     object_temp(i) = object(i)
    enddo

    call release (object_name, object)

    allocate (object (dim1))

    do i = 1, dim_min
     object(i) = object_temp(i)
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_string_row_1

!===========================================================================
  subroutine alloc_integer_1 (object_name, object, dim1) ! new
!  recursive subroutine alloc_integer_1 (object_name, object, dim1) ! old
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1

! output
  integer, allocatable          :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_integer_1'
  integer all_err
  integer i, object_dim, dim_min
  integer, allocatable          :: object_temp (:)

! begin

! allocate object if not already allocated
  if (.not. allocated(object)) then

   allocate (object(dim1), stat = all_err)

   if (all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die (lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:) = 0

   else

! resize object if already allocated with different dimension
   object_dim = size(object)

   if (object_dim /= dim1) then

    dim_min =  min(object_dim, dim1)

!    call alloc ('object_temp', object_temp, dim_min) ! old
    allocate (object_temp(dim_min)) ! new

    do i = 1, dim_min
     object_temp(i) = object(i)
    enddo

    call release (object_name, object)

!    call alloc (object_name, object, dim1) ! old
    allocate (object(dim1)) ! new
    object(:) = 0 !new

    do i = 1, dim_min
     object(i) = object_temp(i)
    enddo

!    call release ('object_temp', object_temp) ! old
    deallocate (object_temp) ! new

! this below does not work because reshape cannot reallocate an array with a larger size
!    object = reshape (object, (/ dim1 /))
!    if (object_dim < dim1) then
!      object (object_dim+1:dim1) = 0
!   endif

   endif

  endif

  end subroutine alloc_integer_1

!===========================================================================
  subroutine alloc_integer_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2

! output
  integer, allocatable          :: object (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_integer_2'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  integer, allocatable          :: object_temp (:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:,:) = 0

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)

   if (object_dim1 /= dim1 .or. object_dim2 /= dim2) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)

    allocate (object_temp(dim_min1, dim_min2))

    do i = 1, dim_min1
     do j = 1, dim_min2
      object_temp(i,j) = object(i,j)
     enddo
    enddo

    call release (object_name, object)

    allocate (object(dim1, dim2))
    object(:,:) = 0

    do i = 1, dim_min1
     do j = 1, dim_min2
      object(i,j) = object_temp(i,j)
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_integer_2

!===========================================================================
  subroutine alloc_integer_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2, dim3

! output
  integer, allocatable          :: object (:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_integer_3'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  integer k, object_dim3, dim_min3
  integer, allocatable          :: object_temp (:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2,dim3), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2, dim3
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:,:,:) = 0

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)

    allocate (object_temp(dim_min1, dim_min2, dim_min3))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       object_temp(i,j,k) = object(i,j,k)
      enddo
     enddo
    enddo

    call release (object_name, object)

    allocate (object(dim1, dim2, dim3))
    object(:,:,:) = 0

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       object(i,j,k) = object_temp(i,j,k)
      enddo
     enddo
    enddo

    deallocate(object_temp)

   endif

  endif

  end subroutine alloc_integer_3

!===========================================================================
  subroutine alloc_integer_row_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 15 Dec 2004
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1

! output
  type (type_integer_row), allocatable, intent(out)   :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_integer_row_1'
  integer all_err
  integer i, object_dim, dim_min
  type (type_integer_row), allocatable          :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die (lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim = size(object)

   if(object_dim /= dim1) then

    dim_min =  min(object_dim, dim1)

    allocate (object_temp(dim_min))

    do i = 1, dim_min
     object_temp(i) = object(i)
    enddo

    call release (object_name, object)

    allocate (object(dim1))

    do i = 1, dim_min
     object(i) = object_temp(i)
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_integer_row_1

!===========================================================================
  subroutine alloc_double_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1

! output
  real(dp) , allocatable          :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_double_1'
  integer all_err
  integer i, object_dim, dim_min
  real(dp) , allocatable          :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:) = 0.d0

   else

! resize object if already allocated with different dimension
   object_dim = size(object)

   if(object_dim /= dim1) then

    dim_min =  min(object_dim, dim1)

    allocate (object_temp (dim_min))

    do i = 1, dim_min
     object_temp(i) = object(i)
    enddo

    call release (object_name, object)

    allocate (object(dim1))
    object(:) = 0.d0

    do i = 1, dim_min
     object(i) = object_temp(i)
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_double_1

!===========================================================================
  subroutine alloc_double_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2

! output
  real(dp) , allocatable          :: object (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_double_2'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  real(dp) , allocatable          :: object_temp (:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:,:) = 0.d0

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)

    allocate (object_temp (dim_min1, dim_min2))

    do i = 1, dim_min1
     do j = 1, dim_min2
      object_temp(i,j) = object(i,j)
     enddo
    enddo

    call release (object_name, object)

    allocate (object (dim1, dim2))
    object(:,:) = 0.d0

    do i = 1, dim_min1
     do j = 1, dim_min2
      object(i,j) = object_temp(i,j)
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_double_2

!===========================================================================
  subroutine alloc_double_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2, dim3

! output
  real(dp) , allocatable          :: object (:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_double_3'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  integer k, object_dim3, dim_min3
  real(dp) , allocatable          :: object_temp (:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2,dim3), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2, dim3
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:,:,:) = 0.d0

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)

    allocate (object_temp (dim_min1, dim_min2, dim_min3))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       object_temp(i,j,k) = object(i,j,k)
      enddo
     enddo
    enddo

    call release (object_name, object)

    allocate (object (dim1, dim2, dim3))
    object(:,:,:) = 0.d0

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       object(i,j,k) = object_temp(i,j,k)
      enddo
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_double_3

!===========================================================================
  subroutine alloc_double_4 (object_name, object, dim1, dim2, dim3, dim4)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2, dim3, dim4

! output
  real(dp) , allocatable          :: object (:,:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_double_4'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  integer k, object_dim3, dim_min3
  integer l, object_dim4, dim_min4
  real(dp) , allocatable          :: object_temp (:,:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2,dim3,dim4), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2, dim3, dim4
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:,:,:,:) = 0.d0

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)
   object_dim4 = size(object,4)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3 .or. object_dim4 /= dim4) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)
    dim_min4 =  min(object_dim4, dim4)

    allocate (object_temp (dim_min1, dim_min2, dim_min3, dim_min4))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       do l = 1, dim_min4
        object_temp(i,j,k,l) = object(i,j,k,l)
       enddo
      enddo
     enddo
    enddo

    call release (object_name, object)

    allocate (object (dim1, dim2, dim3, dim4))
    object(:,:,:,:) = 0.d0

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       do l = 1, dim_min4
        object(i,j,k,l) = object_temp(i,j,k,l)
       enddo
      enddo
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_double_4

!===========================================================================
  subroutine alloc_double_5 (object_name, object, dim1, dim2, dim3, dim4, dim5)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 26 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2, dim3, dim4, dim5

! output
  real(dp) , allocatable          :: object (:,:,:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_double_5'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  integer k, object_dim3, dim_min3
  integer l, object_dim4, dim_min4
  integer n, object_dim5, dim_min5
  real(dp) , allocatable          :: object_temp (:,:,:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2,dim3,dim4,dim5), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2, dim3, dim4, dim5
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:,:,:,:,:) = 0.d0

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)
   object_dim4 = size(object,4)
   object_dim5 = size(object,5)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3 .or. object_dim4 /= dim4 .or. object_dim5 /= dim5) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)
    dim_min4 =  min(object_dim4, dim4)
    dim_min5 =  min(object_dim5, dim5)

    allocate (object_temp (dim_min1, dim_min2, dim_min3, dim_min4, dim_min5))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       do l = 1, dim_min4
        do n = 1, dim_min5
         object_temp(i,j,k,l,n) = object(i,j,k,l,n)
        enddo
       enddo
      enddo
     enddo
    enddo

    call release (object_name, object)

    allocate (object (dim1, dim2, dim3, dim4, dim5))
    object(:,:,:,:,:) = 0.d0

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       do l = 1, dim_min4
        do n = 1, dim_min5
         object(i,j,k,l,n) = object_temp(i,j,k,l,n)
        enddo
       enddo
      enddo
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_double_5

!===========================================================================
  subroutine alloc_double_row_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 15 Dec 2004
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)        :: object_name
  integer, intent(in)                 :: dim1

! output
  type (type_real_row), allocatable, intent(out)   :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_double_row_1'
  integer all_err
  integer i, object_dim, dim_min
  type (type_real_row), allocatable          :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die (lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim = size(object)

   if(object_dim /= dim1) then

    dim_min =  min(object_dim, dim1)

    allocate (object_temp (dim_min))

    do i = 1, dim_min
     object_temp(i) = object(i)
    enddo

    call release (object_name, object)

    allocate (object (dim1))

    do i = 1, dim_min
     object(i) = object_temp(i)
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_double_row_1

!===========================================================================
  subroutine alloc_complex_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1

! output
  complex(dpc) , allocatable          :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_complex_1'
  integer all_err
  integer i, object_dim, dim_min
  complex(dpc) , allocatable          :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die(lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim = size(object)

   if(object_dim /= dim1) then

    dim_min =  min(object_dim, dim1)

    allocate (object_temp (dim_min))

    do i = 1, dim_min
     object_temp(i) = object(i)
    enddo

    call release (object_name, object)

    allocate (object (dim1))

    do i = 1, dim_min
     object(i) = object_temp(i)
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_complex_1

!===========================================================================
  subroutine alloc_complex_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2

! output
  complex(dpc) , allocatable          :: object (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_complex_2'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  complex(dpc) , allocatable          :: object_temp (:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2
    call die(lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)

    allocate (object_temp (dim_min1, dim_min2))

    do i = 1, dim_min1
     do j = 1, dim_min2
      object_temp(i,j) = object(i,j)
     enddo
    enddo

    call release (object_name, object)

    allocate (object (dim1, dim2))

    do i = 1, dim_min1
     do j = 1, dim_min2
      object(i,j) = object_temp(i,j)
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_complex_2

!===========================================================================
  subroutine alloc_complex_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2, dim3

! output
  complex(dpc) , allocatable          :: object (:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_complex_3'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  integer k, object_dim3, dim_min3
  complex(dpc) , allocatable          :: object_temp (:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2,dim3), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2, dim3
    call die(lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)

    allocate (object_temp (dim_min1, dim_min2, dim_min3))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       object_temp(i,j,k) = object(i,j,k)
      enddo
     enddo
    enddo

    call release (object_name, object)

    allocate (object (dim1, dim2, dim3))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       object(i,j,k) = object_temp(i,j,k)
      enddo
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_complex_3

!===========================================================================
  subroutine alloc_complex_4 (object_name, object, dim1, dim2, dim3, dim4)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2, dim3, dim4

! output
  complex(dpc) , allocatable          :: object (:,:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_complex_4'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  integer k, object_dim3, dim_min3
  integer l, object_dim4, dim_min4
  complex(dpc) , allocatable          :: object_temp (:,:,:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2,dim3,dim4), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2, dim3, dim4
    call die(lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)
   object_dim4 = size(object,4)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3 .or. object_dim4 /= dim4) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)
    dim_min4 =  min(object_dim4, dim4)

    allocate (object_temp (dim_min1, dim_min2, dim_min3, dim_min4))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       do l = 1, dim_min4
        object_temp(i,j,k,l) = object(i,j,k,l)
       enddo
      enddo
     enddo
    enddo

    call release (object_name, object)

    allocate (object (dim1, dim2, dim3, dim4))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       do l = 1, dim_min4
        object(i,j,k,l) = object_temp(i,j,k,l)
       enddo
      enddo
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_complex_4

!===========================================================================
  subroutine alloc_logical_1 (object_name, object, dim1)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1

! output
  logical , allocatable          :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_logical_1'
  integer all_err
  integer i, object_dim, dim_min
  logical , allocatable          :: object_temp (:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1
    call die(lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim = size(object)

   if(object_dim /= dim1) then

    dim_min =  min(object_dim, dim1)

    allocate (object_temp (dim_min))

    do i = 1, dim_min
     object_temp(i) = object(i)
    enddo

    call release (object_name, object)

    allocate (object (dim1))

    do i = 1, dim_min
     object(i) = object_temp(i)
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_logical_1

!===========================================================================
  subroutine alloc_logical_2 (object_name, object, dim1, dim2)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2

! output
  logical , allocatable          :: object (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_logical_2'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  logical , allocatable          :: object_temp (:,:)

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(dim1,dim2), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2
    call die(lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)

    allocate (object_temp (dim_min1, dim_min2))

    do i = 1, dim_min1
     do j = 1, dim_min2
      object_temp(i,j) = object(i,j)
     enddo
    enddo

    call release (object_name, object)

    allocate (object (dim1, dim2))

    do i = 1, dim_min1
     do j = 1, dim_min2
      object(i,j) = object_temp(i,j)
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_logical_2

!===========================================================================
  subroutine alloc_logical_3 (object_name, object, dim1, dim2, dim3)
!---------------------------------------------------------------------------
! Description : allocate a object
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)                    :: object_name
  integer                             :: dim1, dim2, dim3

! output
  logical , allocatable               :: object (:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_logical_3'
  integer all_err
  integer i, object_dim1, dim_min1
  integer j, object_dim2, dim_min2
  integer k, object_dim3, dim_min3
  logical, allocatable          :: object_temp (:,:,:)

! begin

! allocate object if not already allocated
  if (.not. allocated(object)) then

   allocate (object(dim1,dim2,dim3), stat = all_err)

   if (all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': dimensions are ', dim1, dim2, dim3
    call die(lhere,'allocation failed')
   endif

   else

! resize object if already allocated with different dimension
   object_dim1 = size(object,1)
   object_dim2 = size(object,2)
   object_dim3 = size(object,3)

   if(object_dim1 /= dim1 .or. object_dim2 /= dim2 .or. object_dim3 /= dim3) then

    dim_min1 =  min(object_dim1, dim1)
    dim_min2 =  min(object_dim2, dim2)
    dim_min3 =  min(object_dim3, dim3)

    allocate (object_temp (dim_min1, dim_min2, dim_min3))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       object_temp(i,j,k) = object(i,j,k)
      enddo
     enddo
    enddo

    call release (object_name, object)

    allocate (object (dim1, dim2, dim3))

    do i = 1, dim_min1
     do j = 1, dim_min2
      do k = 1, dim_min3
       object(i,j,k) = object_temp(i,j,k)
      enddo
     enddo
    enddo

    deallocate (object_temp)

   endif

  endif

  end subroutine alloc_logical_3

!===========================================================================
  subroutine alloc_range_double_1 (object_name, object, lower1, upper1)
!---------------------------------------------------------------------------
! Description : allocate a object with its range
!
! Created     : J. Toulouse, 05 Jun 2010
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)     :: object_name
  integer, intent(in)              :: lower1, upper1

! output
  real(dp) , allocatable           :: object (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_range_double_1'
  integer all_err, object_lower1, object_upper1

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(lower1:upper1), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
    write(6,*) trim(lhere),': range is ', lower1, upper1
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:) = 0.d0

   else

!  reallocate object if already allocated with different dimension
   object_lower1 = lbound(object, 1)
   object_upper1 = ubound(object, 1)

   if (object_lower1 /= lower1 .or. object_upper1 /= upper1) then

    call release (object_name, object)

    allocate (object(lower1:upper1), stat = all_err)

    if(all_err /= 0) then
     write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
     write(6,*) trim(lhere),': range is ', lower1, upper1
     call die(lhere,'allocation failed')
    endif

!   initialize array to 0
    object(:) = 0.d0

   endif

  endif

  end subroutine alloc_range_double_1

!===========================================================================
  subroutine alloc_range_double_2 (object_name, object, lower1, upper1, lower2, upper2)
!---------------------------------------------------------------------------
! Description : allocate a object with its range
!
! Created     : J. Toulouse, 05 Jun 2010
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)     :: object_name
  integer, intent(in)              :: lower1, upper1, lower2, upper2

! output
  real(dp) , allocatable           :: object (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_range_double_2'
  integer all_err, object_lower1, object_upper1, object_lower2, object_upper2

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(lower1:upper1,lower2:upper2), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
     write(6,*) trim(lhere),': ranges are ', lower1, upper1, lower2, upper2
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:,:) = 0.d0

   else

!  reallocate object if already allocated with different dimension
   object_lower1 = lbound(object, 1)
   object_upper1 = ubound(object, 1)
   object_lower2 = lbound(object, 2)
   object_upper2 = ubound(object, 2)

   if (object_lower1 /= lower1 .or. object_upper1 /= upper1 .or. object_lower2 /= lower2 .or. object_upper2 /= upper2) then

    call release (object_name, object)

    allocate (object(lower1:upper1,lower2:upper2), stat = all_err)

    if(all_err /= 0) then
     write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
     write(6,*) trim(lhere),': ranges are ', lower1, upper1, lower2, upper2
     call die(lhere,'allocation failed')
    endif

!   initialize array to 0
    object(:,:) = 0.d0

   endif

  endif

  end subroutine alloc_range_double_2

!===========================================================================
  subroutine alloc_range_double_3 (object_name, object, lower1, upper1, lower2, upper2, lower3, upper3)
!---------------------------------------------------------------------------
! Description : allocate a object with its range
!
! Created     : J. Toulouse, 05 Jun 2010
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)     :: object_name
  integer, intent(in)              :: lower1, upper1, lower2, upper2, lower3, upper3

! output
  real(dp) , allocatable           :: object (:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'alloc_range_double_3'
  integer all_err, object_lower1, object_upper1, object_lower2, object_upper2, object_lower3, object_upper3

! begin

! allocate object if not already allocated
  if(.not. allocated(object)) then

   allocate (object(lower1:upper1,lower2:upper2,lower3:upper3), stat = all_err)

   if(all_err /= 0) then
    write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
     write(6,*) trim(lhere),': ranges are ', lower1, upper1, lower2, upper2, lower3, upper3
    call die(lhere,'allocation failed')
   endif

!  initialize array to 0
   object(:,:,:) = 0.d0

   else

!  reallocate object if already allocated with different dimension
   object_lower1 = lbound(object, 1)
   object_upper1 = ubound(object, 1)
   object_lower2 = lbound(object, 2)
   object_upper2 = ubound(object, 2)
   object_lower3 = lbound(object, 3)
   object_upper3 = ubound(object, 3)

   if (object_lower1 /= lower1 .or. object_upper1 /= upper1 .or. object_lower2 /= lower2 .or. object_upper2 /= upper2  &
       .or. object_lower3 /= lower3 .or. object_upper3 /= upper3) then

    call release (object_name, object)

   allocate (object(lower1:upper1,lower2:upper2,lower3:upper3), stat = all_err)

    if(all_err /= 0) then
     write(6,*) trim(lhere),': allocation for object ', trim(object_name),' failed'
     write(6,*) trim(lhere),': ranges are ', lower1, upper1, lower2, upper2, lower3, upper3
     call die(lhere,'allocation failed')
    endif

!   initialize array to 0
    object(:,:,:) = 0.d0

   endif

  endif

  end subroutine alloc_range_double_3

!===========================================================================
  subroutine release_integer_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  integer, allocatable     :: object(:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_integer_1'
  integer all_err


! build object if not allocated
  if (allocated(object)) then
   deallocate(object, stat = all_err)
!   write(6,'(4a,i3)') trim(lhere),': deallocate object ', trim(object_name),' => size(object)=',size(object)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die (lhere, 'deallocation failed')
   endif
  endif


  end subroutine release_integer_1

!===========================================================================
  subroutine release_integer_2 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  integer, allocatable    :: object(:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_integer_2'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_integer_2

!===========================================================================
  subroutine release_integer_3 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  integer, allocatable    :: object(:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_integer_3'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_integer_3

!===========================================================================
  subroutine release_integer_row_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                          :: object_name

! input/output
  type (type_integer_row), allocatable, intent(inout)   :: object(:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_integer_row_1'
  integer all_err


! build object if not allocated
  if (allocated(object)) then
   deallocate(object, stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die (lhere, 'deallocation failed')
   endif
  endif

  end subroutine release_integer_row_1

!===========================================================================
  subroutine release_double_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  real(dp), allocatable    :: object(:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_double_1'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_double_1

!===========================================================================
  subroutine release_double_2 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  real(dp), allocatable    :: object(:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_double_2'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_double_2

!===========================================================================
  subroutine release_double_3 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  real(dp), allocatable    :: object(:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_double_3'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_double_3

!===========================================================================
  subroutine release_double_4 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  real(dp), allocatable    :: object(:,:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_double_4'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_double_4

!===========================================================================
  subroutine release_double_5 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 26 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)            :: object_name
  real(dp), allocatable    :: object(:,:,:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_double_5'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_double_5

!===========================================================================
  subroutine release_double_row_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                       :: object_name

! input/output
  type (type_real_row), allocatable, intent(inout)   :: object(:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_double_row_1'
  integer all_err


! build object if not allocated
  if (allocated(object)) then
   deallocate(object, stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die (lhere, 'deallocation failed')
   endif
  endif

  end subroutine release_double_row_1

!===========================================================================
  subroutine release_complex_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  complex(dpc), allocatable    :: object(:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_complex_1'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_complex_1

!===========================================================================
  subroutine release_complex_2 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  complex(dpc), allocatable    :: object(:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_complex_2'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_complex_2

!===========================================================================
  subroutine release_complex_3 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  complex(dpc), allocatable    :: object(:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_complex_3'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_complex_3

!===========================================================================
  subroutine release_complex_4 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  complex(dpc), allocatable    :: object(:,:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_complex_4'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_complex_4

!===========================================================================
  subroutine release_logical_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  logical, allocatable    :: object(:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_logical_1'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_logical_1

!===========================================================================
  subroutine release_logical_2 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  logical, allocatable    :: object(:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_logical_2'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_logical_2

!===========================================================================
  subroutine release_logical_3 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character*(*)            :: object_name
  logical, allocatable    :: object(:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_logical_3'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_logical_3

!===========================================================================
  subroutine release_string_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*)            :: object_name
  character(len=max_string_len), allocatable    :: object(:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_string_1'
  integer all_err

! begin

! build object if not allocated
  if(allocated(object)) then
   deallocate(object,stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die(lhere,'deallocation failed')
   endif
  endif

  end subroutine release_string_1

!===========================================================================
  subroutine release_string_row_1 (object_name, object)
!---------------------------------------------------------------------------
! Description : deallocate object
!
! Created     : J. Toulouse, 11 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in)                       :: object_name

! input/output
  type (type_string_row), allocatable, intent(inout)   :: object(:)

! local
  character(len=max_string_len_rout), save :: lhere = 'release_string_row_1'
  integer all_err


! build object if not allocated
  if (allocated(object)) then
   deallocate(object, stat = all_err)
   if(all_err /= 0) then
    write(6,*) trim(lhere),': deallocation for object ',trim(object_name),' failed'
    call die (lhere, 'deallocation failed')
   endif
  endif

  end subroutine release_string_row_1

!===========================================================================
  subroutine flatten_integer_2 (object_1, object_2, dim1, dim2)
!---------------------------------------------------------------------------
! Description : 2D array -> 1D array
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in)     :: object_2 (:,:)
  integer, intent(in)     :: dim1, dim2

! output
  integer, intent(out)    :: object_1 (:)

! local
  integer i, j

! begin
  do i = 1, dim1
    do j = 1, dim2
     object_1 ((i-1)*dim2 + j) = object_2 (i,j)
    enddo
  enddo

  end subroutine flatten_integer_2

!===========================================================================
  subroutine flatten_double_2 (object_1, object_2, dim1, dim2)
!---------------------------------------------------------------------------
! Description : 2D array -> 1D array
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)        :: object_2 (:,:)
  integer, intent(in)         :: dim1, dim2

! output
  real(dp), intent(out)       :: object_1 (:)

! local
  integer i, j

! begin
  do i = 1, dim1
    do j = 1, dim2
     object_1 ((i-1)*dim2 + j) = object_2 (i,j)
    enddo
  enddo

  end subroutine flatten_double_2

!===========================================================================
  subroutine unflatten_integer_2 (object_1, object_2, dim1, dim2)
!---------------------------------------------------------------------------
! Description : 1D array -> 2D array
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in)     :: object_1 (:)
  integer, intent(in)     :: dim1, dim2

! output
  integer, intent(out)    :: object_2 (:,:)

! local
  integer i, j

  do i = 1, dim1
    do j = 1, dim2
     object_2 (i,j) = object_1 ((i-1)*dim2 + j)
    enddo
  enddo

  end subroutine unflatten_integer_2

!===========================================================================
  subroutine unflatten_double_2 (object_1, object_2, dim1, dim2)
!---------------------------------------------------------------------------
! Description : 1D array -> 2D array
!
! Created     : J. Toulouse, 19 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)                 :: object_1 (:)
  integer, intent(in)                  :: dim1, dim2

! output
  real(dp), intent(out)                :: object_2 (:,:)

! local
  integer i, j

! begin
  do i = 1, dim1
    do j = 1, dim2
     object_2 (i,j) = object_1 ((i-1)*dim2 + j)
    enddo
  enddo

  end subroutine unflatten_double_2

!===========================================================================
  subroutine copy_integer_1 (array1, array2)
!---------------------------------------------------------------------------
! Description : copy array1 into array2, allocating array2 if necesssary
!
! Created     : J. Toulouse, 16 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  integer, allocatable, intent(in) :: array1 (:)

! input/output
  integer, allocatable, intent(inout) :: array2 (:)


! begin
  if (mysize(array1) /= 0) then
   call alloc ('array2', array2, mysize(array1))
   array2 (:) = array1 (:)
  endif

  end subroutine copy_integer_1

!===========================================================================
  subroutine copy_double_1 (array1, array2)
!---------------------------------------------------------------------------
! Description : copy array1 into array2, allocating array2 if necesssary
!
! Created     : J. Toulouse, 16 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), allocatable, intent(in) :: array1 (:)

! input/output
  real(dp), allocatable, intent(inout) :: array2 (:)


! begin
  if (mysize(array1) /= 0) then
   call alloc ('array2', array2, mysize(array1))
   array2 (:) = array1 (:)
  endif

  end subroutine copy_double_1

!===========================================================================
  subroutine move_integer_1 (array1, array2)
!---------------------------------------------------------------------------
! Description : move array1 into array2
!
! Created     : J. Toulouse, 27 Mar 2010
!---------------------------------------------------------------------------
  implicit none

! input
  integer, allocatable, intent(in) :: array1 (:)

! input/output
  integer, allocatable, intent(inout) :: array2 (:)


! begin
  call copy (array1, array2)
  call release ('array1', array1)

  end subroutine move_integer_1

!===========================================================================
  subroutine move_double_1 (array1, array2)
!---------------------------------------------------------------------------
! Description : move array1 into array2
!
! Created     : J. Toulouse, 27 Mar 2010
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), allocatable, intent(in) :: array1 (:)

! input/output
  real(dp), allocatable, intent(inout) :: array2 (:)


! begin
  call copy (array1, array2)
  call release ('array1', array1)

  end subroutine move_double_1

!===========================================================================
  subroutine append_integer_0_to_1 (array, intg)
!---------------------------------------------------------------------------
! Description : append a integer into an array
!
! Created     : J. Toulouse, 18 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: intg

! input/output
  integer, allocatable, intent(inout) :: array(:)

! local
  integer array_nb

! begin
  array_nb = mysize(array)
  array_nb = array_nb + 1
  call alloc ('array', array, array_nb)
  array (array_nb) = intg

  end subroutine append_integer_0_to_1

!===========================================================================
  subroutine append_double_0_to_1 (array, double)
!---------------------------------------------------------------------------
! Description : append a real(dp) number into an array
!
! Created     : J. Toulouse, 13 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in) :: double

! input/output
  real(dp), allocatable, intent(inout) :: array (:)

! local
  integer array_nb

! begin
  array_nb = mysize(array)
  array_nb = array_nb + 1
  call alloc ('array', array, array_nb)
  array (array_nb) = double

  end subroutine append_double_0_to_1

!===========================================================================
  subroutine append_string_0_to_1 (array, string)
!---------------------------------------------------------------------------
! Description : append a string into an array
!
! Created     : J. Toulouse, 11 Apr 2009
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: string

! input/output
  character(len=*), allocatable, intent(inout) :: array (:)

! local
  integer array_nb !, i

! begin
!  write(6,*) "string=",string
  array_nb = mysize(array)
  array_nb = array_nb + 1
  call alloc ('array', array, array_nb)
  array (array_nb) = string

!  do i=1,array_nb
!  write(6,*) "array >",trim(array(i)),'<'
!  enddo

  end subroutine append_string_0_to_1

!===========================================================================
  subroutine append_integer_1_to_1 (array1, array2)
!---------------------------------------------------------------------------
! Description : append a array of integers into an array
!
! Created     : J. Toulouse, 13 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  integer, allocatable, intent(in) :: array2 (:)

! input/output
  integer, allocatable, intent(inout) :: array1 (:)

! local
  integer i

! begin
  do i = 1, mysize(array2)
    call append (array1, array2(i))
  enddo

  end subroutine append_integer_1_to_1

!===========================================================================
  subroutine append_double_1_to_1 (array1, array2)
!---------------------------------------------------------------------------
! Description : append a array of integers into an array
!
! Created     : J. Toulouse, 13 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), allocatable, intent(in) :: array2 (:)

! input/output
  real(dp), allocatable, intent(inout) :: array1 (:)

! local
  integer i

! begin
  do i = 1, mysize(array2)
    call append (array1, array2(i))
  enddo

  end subroutine append_double_1_to_1

!!===========================================================================
!  subroutine append_once_string (array, string) ! commented out for pathscale compiler
!!---------------------------------------------------------------------------
!! Description : append only once a string to an array
!!
!! Created     : J. Toulouse, 15 Oct 2005
!---------------------------------------------------------------------------
!  implicit none
!
!! input
!  character(len=*) string
!
!! input/output
!  character(len=*), allocatable :: array(:)
!
!! local
!!  character(len=max_string_len_rout) lhere
!  integer i, array_nb
!
!! begin
!!  lhere = 'append_once_string'
!
!! this is needed since size(array) is not 0 if array is not allocated but has been previously allocated
!  if (.not. allocated(array)) then
!   array_nb = 0
!  else
!   array_nb = size(array)
!  endif
!
!! if string is already in array, do nothing
!  if (array_nb /= 0 ) then
!   do i = 1, array_nb
!    if (string == array(i) ) then
!       return
!    endif
!   enddo
!  endif
!
!! add string to array
!  array_nb = array_nb + 1
!  call alloc ('array', array, array_nb)
!  array (array_nb) = string
!
!  end subroutine append_once_string

!===========================================================================
  subroutine append_once_integer (array, intg)
!---------------------------------------------------------------------------
! Description : append only once an integer to an array
!
! Created     : J. Toulouse, 15 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: intg

! input/output
  integer, intent(inout), allocatable :: array (:)

! local
  integer i, array_nb

! begin

  if (.not. allocated(array)) then
   array_nb = 0
  else
   array_nb = size(array)
  endif


! if intg is already in array, do nothing
  if (array_nb /= 0 ) then
   do i = 1, array_nb
    if (intg == array(i) ) then
       return
    endif
   enddo
  endif

! add intg to array
  array_nb = array_nb + 1
  call alloc ('array', array, array_nb)
  array (array_nb) = intg

  end subroutine append_once_integer

!===========================================================================
  function last_integer (array) result(result)
!---------------------------------------------------------------------------
! Description : return last element of an array
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  integer, allocatable, intent(in) :: array (:)

! output
  integer :: result

! local
  character(len=max_string_len_rout), save :: lhere = 'last_integer'
  integer array_nb

! begin
  array_nb = mysize (array)

  if (array_nb == 0) then
   write (6,'(2a)') trim(lhere), ': size (array) = 0'
   call die (lhere)
  endif

  result = array (array_nb)

  end function last_integer

!===========================================================================
  subroutine swap_integer (integer1, integer2)
!---------------------------------------------------------------------------
! Description : swap two integers
!
! Created     : J. Toulouse, 12 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input/output
  integer, intent(inout) :: integer1, integer2

! local
  integer temp

! begin
  temp = integer1
  integer1 = integer2
  integer2 = temp

  end subroutine swap_integer

!===========================================================================
  subroutine swap_double (double1, double2)
!---------------------------------------------------------------------------
! Description : swap two real(dp) numbers
!
! Created     : J. Toulouse, 12 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input/output
  real(dp), intent(inout) :: double1, double2

! local
  real(dp) temp

! begin
  temp = double1
  double1 = double2
  double2 = temp

  end subroutine swap_double

!===========================================================================
  subroutine swap_double1 (double1, double2)
!---------------------------------------------------------------------------
! Description : swap two real(dp) arrays
!
! Created     : J. Toulouse, 30 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input/output
  real(dp), intent(inout) :: double1 (:), double2 (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'swap_double1'
  real(dp), allocatable :: temp (:)
  integer array1_nb, array2_nb

! begin
  array1_nb = size (double1)
  array2_nb = size (double2)

  if (array1_nb /= array2_nb) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': the two arrays are of different sizes: array1_nb=',array1_nb,' /= array2_nb=',array2_nb
   call die (lhere)
  endif

  call alloc ('temp', temp, array1_nb)
  temp (:)  = double1 (:)
  double1 (:) = double2 (:)
  double2 (:) = temp (:)

  end subroutine swap_double1

!===========================================================================
  subroutine sort_integer1 (array)
!---------------------------------------------------------------------------
! Description : sort in increasing order an array of integer
!
! Created     : J. Toulouse, 28 Oct 2005
!---------------------------------------------------------------------------
  implicit none

! input/output
  integer, intent(inout), allocatable :: array(:)

! local
  integer i, j, array_nb

! begin
  array_nb = size(array)

  do i = 1, array_nb
    do j = i+1, array_nb
       if (array(j) < array(i)) then
        call swap (array(i), array(j))
       endif
    enddo
  enddo

  end subroutine sort_integer1

!===========================================================================
  subroutine sort_integer1_integer1 (array1, array2)
!---------------------------------------------------------------------------
! Description : sort in increasing order array1, and then array2
!
! Created     : J. Toulouse, 12 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input/output
  integer, intent(inout) :: array1 (:)
  integer, intent(inout) :: array2 (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'sort_integer1_integer1'
  integer i, j, array1_nb, array2_nb

! begin
  array1_nb = size(array1)
  array2_nb = size(array2)

  if (array1_nb /= array2_nb) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': the two arrays are of different sizes: array1_nb=',array1_nb,' /= array2_nb=',array2_nb
   call die (lhere)
  endif

  do i = 1, array1_nb
    do j = i+1, array1_nb

!       order using array1
       if (array1(j) < array1(i)) then
        call swap (array1(i), array1(j))
        call swap (array2(i), array2(j))
       endif

!       if array1(j) == array1(i), order using array2
       if (array1(j) == array1(i)) then
         if (array2(j) < array2(i)) then
           call swap (array1(i), array1(j))
           call swap (array2(i), array2(j))
         endif
       endif

    enddo
  enddo

  end subroutine sort_integer1_integer1

!===========================================================================
  subroutine sort_integer1_integer1_double1 (array1, array2, array3)
!---------------------------------------------------------------------------
! Description : sort in increasing order array1, then array2, and then array3
!
! Created     : J. Toulouse, 12 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input/output
  integer, intent(inout) :: array1 (:)
  integer, intent(inout) :: array2 (:)
  real(dp), intent(inout) :: array3 (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'sort_integer1_integer1_double1'
  integer i, j, array1_nb, array2_nb, array3_nb

! begin
  array1_nb = size(array1)
  array2_nb = size(array2)
  array3_nb = size(array3)

  if (array1_nb /= array2_nb .or. array1_nb /= array3_nb .or. array2_nb /= array3_nb) then
   write(6,'(2a,i8,a,i8,a,i8)') trim(lhere),': arrays are of different sizes: array1_nb=',array1_nb,', array2_nb=',array2_nb,', array3_nb=',array3_nb
   call die (lhere)
  endif

  do i = 1, array1_nb
    do j = i+1, array1_nb

!       order using array1
       if (array1(j) < array1(i)) then
        call swap (array1(i), array1(j))
        call swap (array2(i), array2(j))
        call swap (array3(i), array3(j))
       endif

!       if array1(j) == array1(i), order using array2
       if (array1(j) == array1(i)) then
         if (array2(j) < array2(i)) then
           call swap (array1(i), array1(j))
           call swap (array2(i), array2(j))
           call swap (array3(i), array3(j))
         endif

!        if array2(j) == array2(i), order using array3
         if (array2(j) == array2(i)) then
         if (array3(j) < array3(i)) then
           call swap (array1(i), array1(j))
           call swap (array2(i), array2(j))
           call swap (array3(i), array3(j))
         endif
       endif
       endif

    enddo
  enddo

  end subroutine sort_integer1_integer1_double1

!===========================================================================
  subroutine sort_double1_double1_double2 (array1, array2, array3)
!---------------------------------------------------------------------------
! Description : sort in increasing order array1, then array2 (and not array3)
!
! Created     : J. Toulouse, 30 Jan 2007
!---------------------------------------------------------------------------
  implicit none

! input/output
  real(dp), intent(inout) :: array1 (:)
  real(dp), intent(inout) :: array2 (:)
  real(dp), intent(inout) :: array3 (:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'sort_double1_double1_double2'
  integer i, j, array1_nb, array2_nb, array3_nb

! begin
  array1_nb = size(array1)
  array2_nb = size(array2)
  array3_nb = size(array3,2)

  if (array1_nb /= array2_nb .or. array1_nb /= array3_nb .or. array2_nb /= array3_nb) then
   write(6,'(2a,i8,a,i8,a,i8)') trim(lhere),': arrays are of different sizes: array1_nb=',array1_nb,', array2_nb=',array2_nb,', array3_nb=',array3_nb
   call die (lhere)
  endif

  do i = 1, array1_nb
    do j = i+1, array1_nb

!       order using array1
       if (array1(j) < array1(i)) then
        call swap (array1(i), array1(j))
        call swap (array2(i), array2(j))
        call swap (array3(:,i), array3(:,j))
       endif

!       if array1(j) == array1(i), order using array2
       if (array1(j) == array1(i)) then
         if (array2(j) < array2(i)) then
           call swap (array1(i), array1(j))
           call swap (array2(i), array2(j))
           call swap (array3(:,i), array3(:,j))
         endif

       endif

    enddo
  enddo

  end subroutine sort_double1_double1_double2

!===========================================================================
  function is_sorted_integer1 (array)
!---------------------------------------------------------------------------
! Description : is array sorted in increasing order?
!
! Created     : J. Toulouse, 04 Jul 2007
!---------------------------------------------------------------------------
  implicit none

! input/output
  integer, intent(inout) :: array(:)
  logical is_sorted_integer1

! local
  integer i, j, array_nb

! begin
  is_sorted_integer1 = .true.

  array_nb = size(array)
  do i = 1, array_nb
    do j = i+1, array_nb
       if (array(j) < array(i)) then
        is_sorted_integer1 = .false.
        return
       endif
    enddo
  enddo

  end function is_sorted_integer1

!===========================================================================
  subroutine sort_first_dim_integer (array)
!---------------------------------------------------------------------------
! Description : sort in increasing order the first dim of an array of integer
!
! Created     : J. Toulouse, 29 Nov 2005
!---------------------------------------------------------------------------
  implicit none

! input/output
  integer         , allocatable :: array(:,:)

! local
  integer i, j, k, dim1, dim2, temp

! begin

  dim1 = size(array,1)
  dim2 = size(array,2)

  do k = 1, dim2
  do i = 1, dim1
    do j = i+1, dim1
       if (array(j,k) < array(i,k)) then
        temp = array(i,k)
        array(i,k) = array(j,k)
        array(j,k) = temp
       endif
    enddo
  enddo
  enddo

  end subroutine sort_first_dim_integer

! ==============================================================================
  subroutine sort_and_sign (array, sign)
! ------------------------------------------------------------------------------
! Description   : sort array and determine the signature of the sort permutation
!
! Created       : J. Toulouse, 07 Dec 2005
! ------------------------------------------------------------------------------
  implicit none

! input/output
  integer array (:)

! output
  integer sign

! local
  integer array_nb
  integer i, j, temp

  array_nb = size(array)
  sign = 1


  do i = 1, array_nb
    do j = i+1, array_nb

!      swap two elements
       if (array(j) < array(i)) then
        temp = array(i)
        array(i) = array(j)
        array(j) = temp
        sign = sign * (-1)
       endif

    enddo
  enddo

  end subroutine sort_and_sign

! ==============================================================================
  subroutine remove_elt_in_array_integer (array, elt_position)
! ------------------------------------------------------------------------------
! Description   : remove element at position elt_position in array
!
! Created       : J. Toulouse, 15 Dec 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: elt_position

! input/output
  integer, allocatable, intent(inout) :: array (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'remove_elt_in_array_integer'
  integer array_nb

! begin
  array_nb = mysize(array)

  if (array_nb == 0) then
   write(6,'(2a)') trim(lhere),': size(array) = 0'
   call die (lhere)
  endif

  if (elt_position <= 0) then
   write(6,'(2a,i8,a)') trim(lhere),': elt_position=',elt_position,' <= 0'
   call die (lhere)
  endif

  if (elt_position > array_nb) then
   write(6,'(3a,i8,a,i8)') 'Fatal error in routine ',trim(lhere),': elt_position=',elt_position,' >  size(array)=', array_nb
   call die (lhere)
  endif

  array (elt_position:array_nb) = eoshift (array(elt_position:array_nb), shift = 1)

  call alloc ('array', array, array_nb-1)

  end subroutine remove_elt_in_array_integer

! ==============================================================================
  subroutine remove_elt_in_array_double_1 (array, elt_position)
! ------------------------------------------------------------------------------
! Description   : remove element at position elt_position in array
!
! Created       : J. Toulouse, 15 Dec 2006
! ------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: elt_position

! input/output
  real(dp), allocatable, intent(inout) :: array (:)

! local
  character(len=max_string_len_rout), save :: lhere = 'remove_elt_in_array_double_1'
  integer array_nb

! begin
  array_nb = mysize(array)

  if (array_nb == 0) then
   write(6,'(2a)') trim(lhere),': size(array) = 0'
   call die (lhere)
  endif

  if (elt_position <= 0) then
   write(6,'(2a,i8,a)') trim(lhere),': elt_position=',elt_position,' <= 0'
   call die (lhere)
  endif

  if (elt_position > array_nb) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': elt_position=',elt_position,' >  size(array)=', array_nb
   call die (lhere)
  endif

  array (elt_position:array_nb) = eoshift (array(elt_position:array_nb), shift = 1)

  call alloc ('array', array, array_nb-1)

  end subroutine remove_elt_in_array_double_1

! ==================================================================================
  subroutine remove_elt_in_array_double_3 (array, elt_position)
! ----------------------------------------------------------------------------------
! Description  : remove element at position elt_position in last dimension of array
!
! Created      : J. Toulouse, 15 Dec 2006
! ----------------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in) :: elt_position

! input/output
  real(dp), allocatable, intent(inout) :: array (:,:,:)

! local
  character(len=max_string_len_rout), save :: lhere = 'remove_elt_in_array_double_3'
  integer array_nb1, array_nb2, array_nb3

! begin
  array_nb1 = size (array,1)
  array_nb2 = size (array,2)
  array_nb3 = size (array,3)

  if (array_nb3 == 0) then
   write(6,'(2a)') trim(lhere),': size (array,3) = 0'
   call die (lhere)
  endif

  if (elt_position <= 0) then
   write(6,'(2a,i8,a)') trim(lhere),': elt_position=',elt_position,' <= 0'
   call die (lhere)
  endif

  if (elt_position > array_nb3) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': elt_position=',elt_position,' >  size(array,3)=', array_nb3
   call die (lhere)
  endif

  array (:,:,elt_position:array_nb3) = eoshift (array(:,:,elt_position:array_nb3), shift = 1, dim = 3)
  call alloc ('array', array, array_nb1, array_nb2, array_nb3-1)

  end subroutine remove_elt_in_array_double_3

! ==============================================================================
  subroutine replace_elt_in_array (array, elt_old, elt_new)
! ------------------------------------------------------------------------------
! Description   : replace elt_old into elt_new in an array
! Description   : check if elt_old is present only once
!
! Created       : J. Toulouse, 07 Dec 2005
! ------------------------------------------------------------------------------
  implicit none

! input
  integer array (:)
  integer elt_old, elt_new

! local
  character(len=max_string_len_rout), save :: lhere = 'change_elt_in_array'
  integer array_nb, change_nb
  integer i

! begin

  array_nb = size(array)

  change_nb = 0

  do i = 1, array_nb

   if (array (i) == elt_old) then
     array (i) = elt_new
     change_nb = change_nb + 1
   endif

  enddo

  if (change_nb /= 1) then
   write(6,'(2a,i3,a,i3,a)') trim(lhere),': element ', elt_old, ' has been found ', change_nb, ' times in the array.'
   write(6,*) trim(lhere),': array=',array
   call die (lhere)
  endif

  end subroutine replace_elt_in_array

!===========================================================================
  function array_is_zero_double (array) result(result)
!---------------------------------------------------------------------------
! Description : test if array is zero
!
! Created     : J. Toulouse, 16 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)  :: array (:)

! output
  logical :: result

! local
  integer i

! begin

  do i = 1, size(array)
    if (array (i) /= 0.d0) then
     result = .false.
     return
    endif
  enddo

  result = .true.
  return
  end function array_is_zero_double

!===========================================================================
  function arrays_equal_integer (array1, array2) result(result)
!---------------------------------------------------------------------------
! Description : test if two arrays are equal
!
! Created     : J. Toulouse, 29 Nov 2005
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in)  :: array1(:), array2(:)

! output
  logical :: result

! local
  integer i, array1_nb, array2_nb

! begin
  array1_nb = size(array1)
  array2_nb = size(array2)

  if (array1_nb /= array2_nb) then
    result = .false.
    return
  endif

  result = .true.

  do i = 1, array1_nb
    if (array1(i) /= array2(i)) then
     result = .false.
     exit
    endif
  enddo

  return
  end function arrays_equal_integer

!===========================================================================
  function arrays_equal_double_1 (array1, array2) result(result)
!---------------------------------------------------------------------------
! Description : test if two arrays are equal
!
! Created     : J. Toulouse, 12 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)  :: array1(:)
  real(dp), intent(in)  :: array2(:)

! output
  logical :: result

! local
  integer i, array1_nb, array2_nb

! begin
  array1_nb = size(array1)
  array2_nb = size(array2)

  if (array1_nb /= array2_nb) then
    result = .false.
    return
  endif

  result = .true.

  do i = 1, array1_nb
    if (array1(i) /= array2(i)) then
     result = .false.
     exit
    endif
  enddo

  return
  end function arrays_equal_double_1

!===========================================================================
  function arrays_equal_double_2 (array1, array2) result(result)
!---------------------------------------------------------------------------
! Description : test if two arrays are equal
!
! Created     : J. Toulouse, 12 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)  :: array1(:,:)
  real(dp), intent(in)  :: array2(:,:)

! output
  logical :: result

! local
  integer i, array1_nb1, array2_nb1
  integer j, array1_nb2, array2_nb2

! begin
  array1_nb1 = size(array1,1)
  array1_nb2 = size(array1,2)
  array2_nb1 = size(array2,1)
  array2_nb2 = size(array2,2)

  if (array1_nb1 /= array2_nb1 .or. array1_nb2 /= array2_nb2) then
    result = .false.
    return
  endif

  result = .true.

  do i = 1, array1_nb1
   do j = 1, array1_nb2
    if (array1(i,j) /= array2(i,j)) then
     result = .false.
     exit
    endif
   enddo
  enddo

  return
  end function arrays_equal_double_2

!===========================================================================
  function arrays_equal_logical (array1, array2) result(result)
!---------------------------------------------------------------------------
! Description : test if two arrays are equal
!
! Created     : J. Toulouse, 29 Nov 2005
!---------------------------------------------------------------------------
  implicit none

! input
  logical, intent(in)  :: array1(:)
  logical, intent(in)  :: array2(:)

! output
  logical  :: result

! local
  integer i, array1_nb, array2_nb

! begin
  array1_nb = size(array1)
  array2_nb = size(array2)

  if (array1_nb /= array2_nb) then
    result = .false.
    return
  endif

  result = .true.

  do i = 1, array1_nb
    if (array1(i) .neqv. array2(i)) then
     result = .false.
     exit
    endif
  enddo

  return
  end function arrays_equal_logical

!===========================================================================
  function elt_in_array_integer (array, intg) result(result)
!---------------------------------------------------------------------------
! Description : test intg is in array
!
! Created     : J. Toulouse, 15 Dec 2006
!---------------------------------------------------------------------------
  implicit none

! input
  integer, intent(in)  :: array (:)
  integer, intent(in)  :: intg

! output
  logical  :: result

! local
  integer i

! begin
  do i = 1, size (array)
   if (array (i) == intg) then
    result = .true.
    return
   endif
  enddo

  result = .false.
  return

  end function elt_in_array_integer

!===========================================================================
  function elt_in_array_string (array, stg) result(result)
!---------------------------------------------------------------------------
! Description : test if stg is in array
!
! Created     : J. Toulouse, 28 Mar 2008
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=max_string_len), intent(in)  :: array (:)
  character(len=*), intent(in)  :: stg

! output
  logical  :: result

! local
  integer i

! begin
  do i = 1, size (array)
   if (trim(array (i)) == trim(stg)) then
    result = .true.
    return
   endif
  enddo

  result = .false.
  return

  end function elt_in_array_string

!===========================================================================
  subroutine is_equal_or_die_double_1 (array1, array2, tol, print_message)
!---------------------------------------------------------------------------
! Description : test if two arrays are equal within a tolerance
!
! Created     : J. Toulouse, 29 Nov 2005
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp), intent(in)  :: array1(:), array2(:)
  real(dp), intent(in)  :: tol
  logical, optional, intent(in) :: print_message

! local
  character(len=max_string_len_rout), save :: lhere = 'is_equal_or_die_double_1'
  integer i, array1_nb, array2_nb

! begin
  if (present(print_message)) then
   if (print_message) then
     write(6,'(2a)') trim(lhere),': checking equality'
   endif
  endif

  array1_nb = size(array1)
  array2_nb = size(array2)

  if (array1_nb /= array2_nb) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': array1_nb=',array1_nb,' /= array2_nb=',array2_nb
   call die (lhere)
  endif

  do i = 1, array1_nb
    if (dabs(array1(i) -array2(i)) > tol) then
      write(6,'(2a,i8)') trim(lhere),': element # ',i
      write(6,'(2a,es15.8,a,es15.8)') trim(lhere),': array1=',array1(i),' /= array2=',array2(i)
      write(6,'(2a,es15.8)') trim(lhere),': tolerance=',tol
      call die (lhere)
    endif
  enddo

  end subroutine is_equal_or_die_double_1

!===========================================================================
  subroutine is_equal_or_die_double_2 (array1, array2, tol, print_message)
!---------------------------------------------------------------------------
! Description : test if two arrays are equal within a tolerance
!
! Created     : J. Toulouse, 29 Nov 2005
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp)  :: array1(:,:)
  real(dp)  :: array2(:,:)
  real(dp)  :: tol
  logical, optional :: print_message

! local
  character(len=max_string_len_rout), save :: lhere = 'is_equal_or_die_double_2'
  integer i, array1_dim1, array2_dim1
  integer j, array1_dim2, array2_dim2

! begin
  if (present(print_message)) then
   if (print_message) then
     write(6,'(2a)') trim(lhere),': checking equality'
   endif
  endif

  array1_dim1 = size(array1,1)
  array2_dim1 = size(array2,1)
  array1_dim2 = size(array1,2)
  array2_dim2 = size(array2,2)

  if (array1_dim1 /= array2_dim1) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': array1_dim1=',array1_dim1,' /= array2_dim1=',array2_dim1
   call die (lhere)
  endif
  if (array1_dim2 /= array2_dim2) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': array1_dim2=',array1_dim2,' /= array2_dim2=',array2_dim2
   call die (lhere)
  endif

  do i = 1, array1_dim1
   do j = 1, array1_dim2
    if (dabs(array1(i,j) -array2(i,j)) > tol) then
      write(6,'(2a,i8,i8)') trim(lhere),': element # ',i,j
      write(6,'(2a,es15.8,a,es15.8)') trim(lhere),': array1=',array1(i,j),' /= array2=',array2(i,j)
      write(6,'(2a,es15.8)') trim(lhere),': tolerance=',tol
      call die (lhere)
    endif
  enddo
  enddo

  end subroutine is_equal_or_die_double_2

!===========================================================================
  subroutine is_equal_or_die_double_3 (array1, array2, tol, print_message)
!---------------------------------------------------------------------------
! Description : test if two arrays are equal within a tolerance
!
! Created     : J. Toulouse, 29 Nov 2005
!---------------------------------------------------------------------------
  implicit none

! input
  real(dp)  :: array1(:,:,:)
  real(dp)  :: array2(:,:,:)
  real(dp)  :: tol
  logical, optional :: print_message

! local
  character(len=max_string_len_rout), save :: lhere = 'is_equal_or_die_double_3'
  integer i, array1_dim1, array2_dim1
  integer j, array1_dim2, array2_dim2
  integer k, array1_dim3, array2_dim3

! begin
  if (present(print_message)) then
   if (print_message) then
     write(6,'(2a)') trim(lhere),': checking equality'
   endif
  endif

  array1_dim1 = size(array1,1)
  array2_dim1 = size(array2,1)
  array1_dim2 = size(array1,2)
  array2_dim2 = size(array2,2)
  array1_dim3 = size(array1,3)
  array2_dim3 = size(array2,3)

  if (array1_dim1 /= array2_dim1) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': array1_dim1=',array1_dim1,' /= array2_dim1=',array2_dim1
   call die (lhere)
  endif
  if (array1_dim2 /= array2_dim2) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': array1_dim2=',array1_dim2,' /= array2_dim2=',array2_dim2
   call die (lhere)
  endif
  if (array1_dim3 /= array2_dim3) then
   write(6,'(2a,i8,a,i8)') trim(lhere),': array1_dim3=',array1_dim3,' /= array2_dim3=',array2_dim3
   call die (lhere)
  endif

  do i = 1, array1_dim1
   do j = 1, array1_dim2
    do k = 1, array1_dim3
    if (dabs(array1(i,j,k) -array2(i,j,k)) > tol) then
      write(6,'(2a,i8,i8,i8)') trim(lhere),': element # ',i,j,k
      write(6,'(2a,es15.8,a,es15.8)') trim(lhere),': array1=',array1(i,j,k),' /= array2=',array2(i,j,k)
      write(6,'(2a,es15.8)') trim(lhere),': tolerance=',tol
      call die (lhere)
    endif
  enddo
  enddo
  enddo

  end subroutine is_equal_or_die_double_3

!===========================================================================
  subroutine write_array_double_1 (array_name, array)
!---------------------------------------------------------------------------
! Description : write array
!
! Created     : J. Toulouse, 06 Dec 2005
!---------------------------------------------------------------------------
  implicit none

! input
  character(len=*) array_name
  real(dp)  :: array(:)

! local
  character(len=max_string_len_rout), save :: lhere = 'wirte_array_double_1'
  integer i, array_nb

! begin

  array_nb = size(array)

  if (array_nb <= 0) then
   write(6,*) trim(lhere),': array_nb=',array_nb,' <= 0'
   call die (lhere)
  endif

  do i = 1, array_nb
    write(6,'(a,a,i3,a,es15.8)') trim (array_name), '(',i,')=', array(i)
  enddo

  end subroutine write_array_double_1

! ================================================================================
  function cpu_to_string (rseconds)
! --------------------------------------------------------------------------------
! Description : converts "rsecs" into a string "..h..m..s..c"
!
! Created     : F. Colonna, 06 Oct 2000
! Modified    : J. Toulouse, 22 Jan 2009
! --------------------------------------------------------------------------------
  implicit none

! input/ouput
  character(len=20) :: cpu_to_string
  real(dp), intent(in) :: rseconds

! local
  character(len=20) :: string
  real(dp)          :: hours, minutes, seconds, cents
  integer           :: ihours, iminutes, iseconds, icents

! begin
  hours    = rseconds/3600.d0
  ihours   = int(hours)

  minutes  = (hours - ihours)*60.d0
  iminutes = int(minutes)

  seconds  = (minutes - iminutes)*60d0
  iseconds = int(seconds)

  cents    = (seconds-iseconds)*100.d0
  icents   = int(cents)

  write(string,'(i4.4,a,i2.2,a,i2.2,a,i2.2,a)') ihours,'h',iminutes,'m',iseconds,'s',icents,'c'

  cpu_to_string = string

  end function cpu_to_string

! ================================================================================
  subroutine print_cpu_time
! --------------------------------------------------------------------------------
! Description : print current cpu time elapsed since start of the execution
!
! Created     : J. Toulouse, 16 Dec 2005
! --------------------------------------------------------------------------------
  implicit none

! local:
  real(dp)              :: cpu_duration, cpu_time_now
  character(len=max_string_len) :: cpu_string

! begin
  call cpu_time (cpu_time_now)

  cpu_duration = cpu_time_now - cpu_time_start
  cpu_string = cpu_to_string (cpu_duration)
!  write(6,'(a,f,a,a)') 'CPU time=',cpu_duration,' s = ',trim(cpu_string)

  if (idtask == 0) then
   write(6,'(2a)') 'Total CPU time is ',trim(cpu_string)
  endif

  end subroutine print_cpu_time

! ================================================================================
  subroutine print_cpu_time_in_seconds (message)
! --------------------------------------------------------------------------------
! Description : print current cpu time elapsed since start of the execution
! Description : and since last check, in second
!
! Created     : J. Toulouse, 27 Mar 2007
! --------------------------------------------------------------------------------
  implicit none

! local:
  real(dp)  :: cpu_duration, cpu_time_now, cpu_duration_last
  character(len=*) :: message

! begin
  call cpu_time (cpu_time_now)
  cpu_duration = cpu_time_now - cpu_time_start
  cpu_duration_last = cpu_time_now - cpu_time_last
  cpu_time_last = cpu_time_now

  if (idtask == 0) then
   write(6,'(2a,f11.2,a,f11.2,a)') trim(message),' (total CPU time is',cpu_duration,' s, CPU time since last check is',cpu_duration_last,' s)'
  endif

  end subroutine print_cpu_time_in_seconds

! ================================================================================
  subroutine get_date (date_nice)
! --------------------------------------------------------------------------------
! Description : return system current date in nice format
!
! Created     : J. Toulouse, 10 Jul 2007
! --------------------------------------------------------------------------------
  implicit none

! output
  character (len=*), intent(out) :: date_nice

! local:
  character (len=max_string_len) :: date, time, zone

! begin
  call date_and_time (date, time, zone)
  write (date_nice,'(a4,a,a2,a,a2,x,a2,a,a2,a,a2,a,a5,a)') date(1:4),'-',date(5:6),'-',date(7:8),time(1:2),':',time(3:4),':',time(4:6),' (',zone(1:5),')'

  end subroutine get_date

! ================================================================================
  subroutine require ( lhere, string, bool)
! --------------------------------------------------------------------------------
! Description :  die if bool = false
!
! Created     : J. Toulouse, 17 Feb 2006
! --------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), intent(in) :: lhere !fp
  character(len=*), intent(in) :: string
  logical, intent(in)          :: bool

  if (.not. bool) then
    write(6,'(6a)') 'require: in routine ',trim(lhere),', condition ', trim(string), ' is false.' !fp
    call die (lhere)
  endif

  end subroutine require

! ! ================================================================================
!   subroutine is_a_number_or_die (array_name, array)
! ! --------------------------------------------------------------------------------
! ! Description : die if one element of array is not a number
! !
! ! Created     : J. Toulouse, 24 Oct 2006
! ! --------------------------------------------------------------------------------
!   implicit none

! ! input
!   character(len=*), intent(in) :: array_name
!   real(dp)             :: array (:)

! ! local
!   character(len=max_string_len_rout), save :: lhere = 'is_a_number_or_die'
!   integer i

! ! comment this out for xlf
! !  do i = 1, size(array)
! !   if (array(i) == "NaN") then
! !    write(6,'(6a,i3,a,f,a)') trim(lhere),': in routine ',trim(here),', ', trim(array_name), ' (',i,')=',array(i), ' is not a number!'
! !    call die (lhere)
! !   endif
! !  enddo

!   end subroutine is_a_number_or_die

!===========================================================================
  function mysize_integer_1 (array) result(result)
!---------------------------------------------------------------------------
! Description : return size of array if allocated, and 0 if not allocated
! Description : this is needed since size(array) is not 0 if array is not allocated but has been previously allocated
!
! Created     : J. Toulouse, 24 Dec 2006
! --------------------------------------------------------------------------------
  implicit none

! input
  integer, allocatable, intent(in) :: array (:)

! output
  integer :: result


! begin

  if (.not. allocated (array)) then
   result = 0
  else
   result = size(array)
  endif

  return
  end function mysize_integer_1

!===========================================================================
  function mysize_double_1 (array) result(result)
!---------------------------------------------------------------------------
! Description : return size of array if allocated, and 0 if not allocated
! Description : this is needed since size(array) is not 0 if array is not allocated but has been previously allocated
!
! Created     : J. Toulouse, 24 Dec 2006
! --------------------------------------------------------------------------------
  implicit none

! input
  real(dp), allocatable, intent(in) :: array (:)

! output
  integer :: result


! begin

  if (.not. allocated (array)) then
   result = 0
  else
   result = size(array)
  endif

  return
  end function mysize_double_1

!===========================================================================
  function mysize_string_1 (array) result(result)
!---------------------------------------------------------------------------
! Description : return size of array if allocated, and 0 if not allocated
! Description : this is needed since size(array) is not 0 if array is not allocated but has been previously allocated
!
! Created     : J. Toulouse, 11 Apr 2009
! --------------------------------------------------------------------------------
  implicit none

! input
  character(len=*), allocatable, intent(in) :: array (:)

! output
  integer :: result

! begin

  if (.not. allocated (array)) then
   result = 0
  else
   result = size(array)
  endif

  return
  end function mysize_string_1


!===========================================================================
  function last_element_integer (array) result(result)
!---------------------------------------------------------------------------
! Description : return last element of array, or 0 if array is empty
!
! Created     : J. Toulouse, 22 Mar 2007
! --------------------------------------------------------------------------------
  implicit none

! input
  integer, allocatable, intent(in) :: array (:)

! output
  integer :: result


! begin
  if (.not. allocated (array)) then
   result = 0
  else
   result = array (size(array))
  endif

  return
  end function last_element_integer

!===========================================================================
  function max_element_integer (array) result(result)
!---------------------------------------------------------------------------
! Description : return the max element of array, or 0 if array is empty
!
! Created     : J. Toulouse, 22 Mar 2007
! --------------------------------------------------------------------------------
  implicit none

! input
  integer, allocatable, intent(in) :: array (:)

! output
  integer :: result


! begin
  if (.not. allocated (array)) then
   result = 0
  else
   result = maxval(array)
  endif

  return
  end function max_element_integer

!===========================================================================
  function array_norm_double_1 (array) result(result)
!---------------------------------------------------------------------------
! Description : return the cartesian norm of array
!
! Created     : J. Toulouse, 25 Mar 2007
! --------------------------------------------------------------------------------
  implicit none

! input
  real(dp), allocatable, intent(in) :: array (:)

! output
  real(dp) :: result

! local
  integer i

! begin
  result = 0.d0

  do i = 1, mysize(array)
    result = result + array (i) **2
  enddo

  result = dsqrt(result)

  return
  end function array_norm_double_1

end module basic_tools_mod
