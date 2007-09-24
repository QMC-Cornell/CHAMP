module deriv_geo_mod

  use all_tools_mod
  use forces_mod

! Declaration of global variables and default values
  integer                        :: param_geo_nb = 0

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

end module deriv_geo_mod
