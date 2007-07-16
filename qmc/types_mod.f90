module types_mod

! 
!  integer, parameter :: I4B = SELECTED_INT_KIND(9)
!  integer, parameter :: I2B = SELECTED_INT_KIND(4)
!  integer, parameter :: I1B = SELECTED_INT_KIND(2)
!  integer, parameter :: LGT = KIND(.true.)
  
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  

! some derived types:
  type type_integer_row
   integer, allocatable              :: row (:)
  end type type_integer_row

  type type_real_row
   real(dp), allocatable     :: row (:)
  end type type_real_row

end module types_mod
