module eloc_mod

  use all_tools_mod
  use montecarlo_mod
  use psi_mod

! Declaration of global variables and default values
  real(dp)                       :: eloc_kin
  real(dp)                       :: eloc_kin_av
  real(dp)                       :: eloc_kin_av_err
  real(dp)                       :: eloc_pot
  real(dp)                       :: eloc_pot_loc
  real(dp)                       :: eloc_pot_nloc
  real(dp)                       :: pe_en
  real(dp)                       :: pe_ee
  real(dp)                       :: eloc_pot_en
  real(dp)                       :: eloc_pot_en_av
  real(dp)                       :: eloc_pot_ee
  real(dp)                       :: eloc_pot_ee_av

  contains

! ==============================================================================
  subroutine eloc_kin_bld
! ------------------------------------------------------------------------------
! Description   : Kinetic local energy
!
! Created       : J. Toulouse, 05 Dec 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dim_i, elec_i
  real(dp) sum_grd_psi_over_psi_square

! header
  if (header_exe) then

   call object_create ('eloc_kin')

   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('sum_lap_lnpsi')
   call object_needed ('grd_psi_over_psi_wlk')

   return

  endif

! begin

! allocations
  call object_associate ('eloc_kin', eloc_kin)
  call object_associate ('eloc_kin_av', eloc_kin_av)
  call object_associate ('eloc_kin_av_err', eloc_kin_av_err)

  sum_grd_psi_over_psi_square = 0.d0

  do dim_i = 1, ndim
    do elec_i = 1, nelec
       sum_grd_psi_over_psi_square = sum_grd_psi_over_psi_square + grd_psi_over_psi_wlk (dim_i, elec_i, 1)**2
    enddo
  enddo

  eloc_kin =  -0.5d0 * (sum_lap_lnpsi + sum_grd_psi_over_psi_square)

 end subroutine eloc_kin_bld

! ==============================================================================
  subroutine eloc_bld
! ------------------------------------------------------------------------------
! Description   : total local energy
!
! Created       : J. Toulouse, 05 Dec 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! header
  if (header_exe) then

!   call object_create ('eloc')

   call object_needed ('eloc_kin')
   call object_needed ('eloc_pot')

   return

  endif

! begin
!  write(6,*) trim(here),': entering'

! allocations
  call object_associate ('eloc', eloc)
!  call object_associate ('eloc_av', eloc_av)
!  call object_associate ('eloc_av_err', eloc_av_err)

  eloc = eloc_kin + eloc_pot

!  write(6,*) trim(here),': eloc=', eloc

 end subroutine eloc_bld

! ==============================================================================
  subroutine eloc_pot_en_bld
! ------------------------------------------------------------------------------
! Description   : electron-nucleus local potential energy
!
! Created       : J. Toulouse, 18 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! header
  if (header_exe) then

   call object_create ('eloc_pot_en')
   call object_average_define ('eloc_pot_en', 'eloc_pot_en_av')

   call object_needed ('pe_en')

   return

  endif

! begin

! allocations
  call object_associate ('eloc_pot_en', eloc_pot_en)
  call object_associate ('eloc_pot_en_av', eloc_pot_en_av)

  eloc_pot_en = pe_en

 end subroutine eloc_pot_en_bld

! ==============================================================================
  subroutine eloc_pot_ee_bld
! ------------------------------------------------------------------------------
! Description   : electron-electron local potential energy
!
! Created       : J. Toulouse, 18 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! header
  if (header_exe) then

   call object_create ('eloc_pot_ee')
   call object_average_define ('eloc_pot_ee', 'eloc_pot_ee_av')

   call object_needed ('pe_ee')

   return

  endif

! begin

! allocations
  call object_associate ('eloc_pot_ee', eloc_pot_ee)
  call object_associate ('eloc_pot_ee_av', eloc_pot_ee_av)

  eloc_pot_ee = pe_ee

 end subroutine eloc_pot_ee_bld

end module eloc_mod
