module types_mod

!  integer, parameter :: i8b = SELECTED_INT_KIND(18)
!  integer, parameter :: i4b = SELECTED_INT_KIND(9)
!  integer, parameter :: i2b = SELECTED_INT_KIND(4)
!  integer, parameter :: i1b = SELECTED_INT_KIND(2)
!  integer, parameter :: LGT = KIND(.true.)

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))


! some derived types:
  type type_string_row
   character(len=100), allocatable  :: row (:)
  end type type_string_row

  type type_integer_row
   integer, allocatable      :: row (:)
  end type type_integer_row

  type type_real_row
   real(dp), allocatable     :: row (:)
  end type type_real_row

end module types_mod
