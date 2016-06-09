module allocations_mod

  use all_tools_mod

  contains

! ===================================================================================
  subroutine common_allocations
! -----------------------------------------------------------------------------------
! Description   : allocate arrays used in FIT, VMC and DMC
! To improved   : allocate only when necessary!
!
! Created       : J. Toulouse, 09 May 2009
! -----------------------------------------------------------------------------------
  use orbitals_mod, only: orb_tot_nb
  use all_modules_mod
  implicit none

  call object_provide ('nelec')
  call object_provide ('nelec_pair')
  call object_provide ('ncent')
  call object_provide ('nbasis')
  call object_provide ('orb_tot_nb')
  call object_provide ('nctype')
  call object_provide ('nparm')
  call object_provide ('nparmjs')

  call alloc ('phin', phin, nbasis, nelec)
  call alloc ('dphin', dphin, 3, nbasis, nelec)
  call alloc ('d2phin', d2phin, nbasis, nelec)

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

  call alloc ('slmin', slmin, nupdn_square, ndetupdn)
  call alloc ('detn', detn, ndetupdn)
  call alloc ('ddeti_detin', ddeti_detin, 3, nelec, ndetupdn)
  call alloc ('d2edeti_detin', d2edeti_detin, nelec, ndetupdn)

  call object_provide ('nparmd')
  call alloc ('deti_det', deti_det, nparmd)
  call alloc ('ddeti_det', ddeti_det, 3, nelec, nparmd)
  call alloc ('d2deti_det', d2deti_det, nparmd)
  call alloc ('detij_det', detij_det, nparmd, nparmd)
! only for pseudo?
  call alloc ('deti_new', deti_new, nparmd)

  call alloc ('rshift', rshift, 3, nelec, ncent)
  call alloc ('rvec_en', rvec_en, 3, nelec, ncent)
  call alloc ('r_en', r_en, nelec, ncent)
  call alloc ('rvec_ee', rvec_ee, 3, nelec_pair)
  call alloc ('r_ee', r_ee, nelec_pair)
  call alloc ('pot_ee', pot_ee, nelec)

  call alloc ('rshift_sav', rshift_sav, 3, ncent)
  call alloc ('rvec_en_sav', rvec_en_sav, 3, ncent)
  call alloc ('r_en_sav', r_en_sav, ncent)
  call alloc ('rvec_ee_sav', rvec_ee_sav, 3, nelec)
  call alloc ('r_ee_sav', r_ee_sav, nelec)

  call alloc ('d1d2a', d1d2a, nctype)
  call alloc ('d2d2a', d2d2a, nctype)
  call alloc ('d1d2b', d1d2b, 2)
  call alloc ('d2d2b', d2d2b, 2)
  call alloc ('didk', didk, nparmjs)

  call alloc ('vd', vd, 3, nelec)
  call alloc ('vj', vj, 3, nelec)
  call alloc ('div_vo', div_vo, nelec)

  call alloc ('ekineo', ekineo, nelec)
  call alloc ('ekinen', ekinen, nelec)

! jastrow
  call alloc ('fso', fso, nelec, nelec)
  call alloc ('fijo', fijo, 3, nelec, nelec)
  call alloc ('d2ijo', d2ijo, nelec, nelec)
  call alloc ('fjo', fjo, 3, nelec)
  call alloc ('fsn', fsn, nelec, nelec)
  call alloc ('fijn', fijn, 3, nelec, nelec)
  call alloc ('d2ijn', d2ijn, nelec, nelec)
  call alloc ('fjn', fjn, 3, nelec)

  call alloc ('gvalue', gvalue, nparmjs)
  call alloc ('g', g, 3, nelec, nparmjs)
  call alloc ('d2g', d2g, nparmjs)
  call alloc ('go', go, nelec, nelec, nparmjs)

! optimization
  call alloc ('denergy', denergy, nparm)

  call alloc ('try', try, NRAD)
  call alloc ('suc', suc, NRAD)
  call alloc ('trunfb', trunfb, NRAD)
  call alloc ('rprob', rprob, NRAD)
  call alloc ('rprobup', rprobup, NRAD)
  call alloc ('rprobdn', rprobdn, NRAD)
  call alloc ('ekin', ekin, NRAD)
  call alloc ('ekin2', ekin2, NRAD)

  if(ibasis.eq.3) then
   call alloc ('cdorb', cdorb, 3, orb_tot_nb)
   call alloc ('cddorb', cddorb, orb_tot_nb)
  endif

! for pseudo
  if(nloc > 0) then
    call alloc ('vps', vps, nelec, ncent, MPS_L)
  endif

 end subroutine common_allocations

end module allocations_mod
