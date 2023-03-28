module hash_table_mod
!Written by Tyler Anderson. 
!A simple hash table (an array of singly linked lists), used to hash determinants.
  implicit none
  
  type bucket
    integer, allocatable :: key(:)
    integer :: val
    type(bucket), pointer :: next => null()
  end type

contains

function equal_arrays(a1, a2)
  implicit none

  logical :: equal_arrays
  integer :: a1(:), a2(:), i

  if (size(a1).eq.size(a2)) then
    equal_arrays = all(a1.eq.a2) 
  else
    equal_arrays = .false.
  endif
end function

function hash_array(a) result(hash)
!hash function written by Junhao Li
  implicit none

  integer, intent(in) :: a(:)
  integer :: hash, i

  hash = 0
  do i=1, size(a)
    hash = xor(hash, a(i) + 1640531527 + ishft(hash, 6) + ishft(hash, -2))
  enddo
end function

subroutine hash_table_lookup(ht, key, bkt)
  implicit none

  type(bucket), pointer :: bkt
  type(bucket), target :: ht(:)
  integer :: key(:), ib

  ib = modulo(hash_array(key), size(ht))+1
  bkt => ht(ib)
  do while (allocated(bkt%key))
    if (equal_arrays(bkt%key, key)) return
    bkt => bkt%next
  enddo
end subroutine

subroutine hash_table_add(ht, key, val)
  implicit none

  type(bucket), pointer :: bkt
  type(bucket), target :: ht(:)
  integer :: key(:), val

  call hash_table_lookup(ht, key, bkt)
  if (.not. allocated(bkt%key)) then
      allocate(bkt%key(size(key)))
      allocate(bkt%next)
      bkt%key = key
      bkt%val = val
  endif
end subroutine

function hash_table_get(ht, key, success)
  implicit none

  integer :: key(:)
  type(bucket), target :: ht(:)
  logical :: success
  type(bucket), pointer :: bkt
  integer hash_table_get

  call hash_table_lookup(ht, key, bkt)
  success = allocated(bkt%key)
  hash_table_get = bkt%val
  if (.not.success) hash_table_get = 0
end function

subroutine hash_table_print(ht)
  implicit none

  type(bucket), pointer :: bkt
  type(bucket), target :: ht(:)
  integer :: i

  do i=1, size(ht)
    bkt => ht(i)
    write(6,*) "LIST ", i
    do while (allocated(bkt%key))
      write(6,*) bkt%key(:)
      bkt => bkt%next
    enddo
  enddo
end subroutine hash_table_print

end module hash_table_mod
