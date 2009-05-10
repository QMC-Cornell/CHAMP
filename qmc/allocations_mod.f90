module allocations_mod

  use all_tools_mod

  contains

! ===================================================================================
  subroutine common_allocations
! -----------------------------------------------------------------------------------
! Description   : allocate arrays used in FIT, VMC and DMC
!
! Created       : J. Toulouse, 09 May 2009
! -----------------------------------------------------------------------------------
  use orbitals_mod, only: orb_tot_nb
  include 'modules.h'
  implicit none
  include 'commons.h'

  call object_provide ('nelec')
  call object_provide ('orb_tot_nb')

  call alloc ('orb', orb, nelec, orb_tot_nb)
  call alloc ('dorb', dorb, 3, nelec, orb_tot_nb)
  call alloc ('ddorb', ddorb, nelec, orb_tot_nb)
  call alloc ('orbe', orbe, orb_tot_nb)
  call alloc ('dorbe', dorbe, 3, orb_tot_nb)
  call alloc ('ddorbe', ddorbe, orb_tot_nb)

  if(ibasis.eq.3) then
   call alloc ('cdorb', cdorb, 3, orb_tot_nb)
   call alloc ('cddorb', cddorb, orb_tot_nb)
  endif

 end subroutine common_allocations

end module allocations_mod
