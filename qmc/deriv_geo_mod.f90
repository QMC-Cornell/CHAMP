module deriv_geo_mod

  use all_tools_mod
  use forces_pulay_mod

! Declaration of global variables and default values
  integer                                :: param_geo_nb = 0
  real(dp), allocatable                  :: dpsi_geo (:)
  real(dp), allocatable                  :: deloc_geo (:)

  contains

! ==============================================================================
  subroutine param_geo_nb_bld
! ------------------------------------------------------------------------------
! Description   : number of geometrical parameters
!
! Created       : J. Toulouse, 26 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local

! header
  if (header_exe) then

   call object_create ('param_geo_nb', param_geo_nb_index)

   call object_needed ('forces_nb')

   return

  endif

! begin
  param_geo_nb = forces_nb

  end subroutine param_geo_nb_bld

! ==============================================================================
  subroutine dpsi_geo_bld
! ------------------------------------------------------------------------------
! Description   : logarithmic derivatives of Psi with respect to geometrical parameters
!
! Created       : J. Toulouse, 06 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_geo_i, force_i

! header
  if (header_exe) then

   call object_create ('dpsi_geo', dpsi_geo_index)

   call object_needed ('param_geo_nb')
   call object_needed ('forces_nb')
   call object_needed ('dpsi_rn')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsi_geo', dpsi_geo, param_geo_nb)

  do param_geo_i = 1, param_geo_nb
   force_i = param_geo_i
   dpsi_geo (param_geo_i) = dpsi_rn (force_i)
  enddo

  end subroutine dpsi_geo_bld

! ==============================================================================
  subroutine deloc_geo_bld
! ------------------------------------------------------------------------------
! Description   : derivatives of local energy with respect to geometrical parameters
!
! Created       : J. Toulouse, 06 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_geo_i, force_i

! header
  if (header_exe) then

   call object_create ('deloc_geo', deloc_geo_index)

   call object_needed ('param_geo_nb')
   call object_needed ('forces_nb')
   call object_needed ('deloc_rn_num')

   return

  endif

! begin

! allocations
  call object_alloc ('deloc_geo', deloc_geo, param_geo_nb)

  do param_geo_i = 1, param_geo_nb
   force_i = param_geo_i
   deloc_geo (param_geo_i) = deloc_rn_num (force_i)
  enddo

  end subroutine deloc_geo_bld

end module deriv_geo_mod
