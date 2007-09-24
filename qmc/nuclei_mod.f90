module nuclei_mod

  use all_tools_mod

! Declaration of global variables and default values
  real(dp), allocatable     :: dist_nn  (:,:)

  contains

! ==============================================================================
  subroutine dist_nn_bld
! ------------------------------------------------------------------------------
! Description   : distance |Ri - Rj| between two nuclei
!
! Created       : J. Toulouse, 26 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer  cent_i, cent_j, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('dist_nn')

   call object_needed ('ncent')
   call object_needed ('ndim')
   call object_needed ('cent')

   return

  endif

! allocations
  call object_alloc ('dist_nn', dist_nn, ncent, ncent)
  dist_nn (:,:) = 0.d0

  do cent_i = 1, ncent
    do cent_j = cent_i+1, ncent
      do dim_i = 1, ndim
        dist_nn (cent_i, cent_j) = dist_nn (cent_i, cent_j) + (cent(dim_i,cent_i) - cent(dim_i,cent_j))**2
      enddo ! dim_i
      dist_nn (cent_i, cent_j) = dsqrt (dist_nn (cent_i, cent_j))
      dist_nn (cent_j, cent_i) = dist_nn (cent_i, cent_j)
    enddo ! cent_j
  enddo ! cent_i

  end subroutine dist_nn_bld

end module nuclei_mod
