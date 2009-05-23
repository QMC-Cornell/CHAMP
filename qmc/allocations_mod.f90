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

  call object_provide ('nelec')
  call object_provide ('orb_tot_nb')

  call alloc ('orb', orb, nelec, orb_tot_nb)
  call alloc ('dorb', dorb, 3, nelec, orb_tot_nb)
  call alloc ('ddorb', ddorb, nelec, orb_tot_nb)
  call alloc ('orbe', orbe, orb_tot_nb)
  call alloc ('dorbe', dorbe, 3, orb_tot_nb)
  call alloc ('ddorbe', ddorbe, orb_tot_nb)

  call object_provide ('nupdn_square')
  call object_provide ('ndetupdn')
  call alloc ('slmui', slmui, nupdn_square, ndetupdn)
  call alloc ('slmdi', slmdi, nupdn_square, ndetupdn)
  call alloc ('fpu', fpu, 3, nupdn_square, ndetupdn)
  call alloc ('fpd', fpd, 3, nupdn_square, ndetupdn)
  call alloc ('fppu', fppu, nupdn_square, ndetupdn)
  call alloc ('fppd', fppd, nupdn_square, ndetupdn)
  call alloc ('detu', detu, ndetupdn)
  call alloc ('detd', detd, ndetupdn)
  call alloc ('ddeti_deti', ddeti_deti, 3, nelec, ndetupdn)
  call alloc ('d2edeti_deti', d2edeti_deti, nelec, ndetupdn)

  call object_provide ('nparmd')
  call alloc ('deti_det', deti_det, nparmd)
  call alloc ('ddeti_det', ddeti_det, 3, nelec, nparmd)
  call alloc ('d2deti_det', d2deti_det, nparmd)
  call alloc ('detij_det', detij_det, nparmd, nparmd)

  if(ibasis.eq.3) then
   call alloc ('cdorb', cdorb, 3, orb_tot_nb)
   call alloc ('cddorb', cddorb, orb_tot_nb)
  endif

 end subroutine common_allocations

end module allocations_mod
