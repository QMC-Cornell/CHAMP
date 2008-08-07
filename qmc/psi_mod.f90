module psi_mod

  use all_tools_mod
  use orbitals_mod
  use montecarlo_mod
  use determinants_mod

! Declaration of global variables and default values
  logical                        :: l_maximize_psi2 = .false.
  integer                        :: max_psi2_trial_nb = 10
  real(dp), allocatable          :: grd_det_unq_dn (:,:,:)
  real(dp), allocatable          :: grd_det_unq_up (:,:,:)
  real(dp), allocatable          :: lap_det_unq_dn (:,:)
  real(dp), allocatable          :: lap_det_unq_up (:,:)
  real(dp), allocatable          :: grd_psid_over_psid (:,:)
  real(dp), allocatable          :: lap_psid_over_psid (:)
  real(dp), allocatable          :: lap_lnpsid (:)
  real(dp)                       :: sum_lap_lnpsid
  real(dp)                       :: sum_lap_lnj
  real(dp)                       :: sum_lap_lnpsi
  real(dp), allocatable          :: grd_psi_over_psi_wlk (:,:,:)
  real(dp), allocatable          :: grd_psi_over_psi_old (:,:)
  real(dp), allocatable          :: grd_psi_over_psi_sq_wlk (:,:)
  real(dp), allocatable          :: lap_psi_over_psi_wlk (:,:)
  real(dp), allocatable          :: div_grd_psi_over_psi_wlk (:,:)

  contains

!===========================================================================
  subroutine wavefunction_menu
!---------------------------------------------------------------------------
! Description : menu for wavefunction
!
! Created     : J. Toulouse, 06 Apr 2008
!---------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'wavefunction_menu'

! begin
  write(6,*) 
  write(6,'(a)') 'Beginning of wavefunction menu ---------------------------------------------------------------------------'


! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for wavefunction menu:'
   write(6,'(a)') ' wavefunction'
   write(6,'(a)') '  maximization ... end: menu for maximimization of wave function square'
   write(6,'(a)') ' end'
   write(6,*)

  case ('maximization')
   call psi2_maximization_menu

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines
  write(6,'(a)') 'End of wavefunction menu ---------------------------------------------------------------------------------'

  end subroutine wavefunction_menu

!===========================================================================
  subroutine psi2_maximization_menu
!---------------------------------------------------------------------------
! Description : menu for maximization of wave function square
!
! Created     : J. Toulouse, 06 Apr 2008
!---------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'psi2_maximization_menu'

! begin
  l_maximize_psi2 = .true.
  write(6,'(a)') 'Requesting maximization of wave function square.'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for maximization menu:'
   write(6,'(a)') ' maximization'
   write(6,'(a)') ' max_psi2_trial_nb = [real] : number of trial of maximization (default=10)'
   write(6,'(a)') ' end'

  case ('max_psi2_trial_nb')
   call get_next_value (max_psi2_trial_nb)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  call maximize_psi2

  end subroutine psi2_maximization_menu

! ==============================================================================
  subroutine grd_det_unq_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of unique spin up and down determinants
!
! Created       : J. Toulouse, 05 Dec 2005
! Modified      : J. Toulouse, 04 Jul 2007: unique determinants
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_unq_up_i, det_unq_dn_i, dim_i, i, j

! header
  if (header_exe) then

   call object_create ('grd_det_unq_up')
   call object_create ('grd_det_unq_dn')


   call object_needed ('ndim')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('dorb')
   call object_needed ('detu')
   call object_needed ('detd')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_det_unq_up', grd_det_unq_up, ndim, nup, ndetup)
  call object_alloc ('grd_det_unq_dn', grd_det_unq_dn, ndim, ndn, ndetdn)

  grd_det_unq_up = 0.d0
  grd_det_unq_dn = 0.d0

! loop over unique spin-up determinants
  do det_unq_up_i = 1, ndetup
     do dim_i = 1, ndim
      do i = 1, nup
       do j= 1, nup
        grd_det_unq_up (dim_i, i, det_unq_up_i) = grd_det_unq_up (dim_i, i, det_unq_up_i)  +  &
          slater_mat_trans_inv_up (i, j, det_unq_up_i) * dorb (dim_i, i, det_unq_orb_lab_srt_up (j, det_unq_up_i))
       enddo
       grd_det_unq_up (dim_i, i, det_unq_up_i) = grd_det_unq_up (dim_i, i, det_unq_up_i) * detu (det_unq_up_i)
      enddo
     enddo ! dim_i
  enddo ! det_unq_up_i

! loop over unique spin-down determinants
  do det_unq_dn_i = 1, ndetdn
     do dim_i = 1, ndim
      do i = 1, ndn
       do j = 1, ndn
        grd_det_unq_dn (dim_i, i, det_unq_dn_i) = grd_det_unq_dn (dim_i, i, det_unq_dn_i)  +  &
          slater_mat_trans_inv_dn (i, j, det_unq_dn_i) * dorb (dim_i, nup + i, det_unq_orb_lab_srt_dn (j, det_unq_dn_i))
       enddo
       grd_det_unq_dn (dim_i, i, det_unq_dn_i) = grd_det_unq_dn (dim_i, i, det_unq_dn_i) * detd (det_unq_dn_i)
      enddo
     enddo ! dim_i
  enddo ! det_unq_dn_i

  end subroutine grd_det_unq_bld

! ==============================================================================
  subroutine lap_det_unq_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of unique spin-up and down determinants
!
! Created       : J. Toulouse, 05 Dec 2005
! Modified      : J. Toulouse, 04 Jul 2007: unique determinants
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_unq_up_i, det_unq_dn_i, dim_i, i, j

! header
  if (header_exe) then

   call object_create ('lap_det_unq_up')
   call object_create ('lap_det_unq_dn')

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('ddorb')
   call object_needed ('detu')
   call object_needed ('detd')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_det_unq_up', lap_det_unq_up, nup, ndetup)
  call object_alloc ('lap_det_unq_dn', lap_det_unq_dn, ndn, ndetdn)

  lap_det_unq_up = 0.d0
  lap_det_unq_dn = 0.d0

! loop over unique spin-up determinants
  do det_unq_up_i = 1, ndetup
      do i = 1, nup
       do j = 1, nup
        lap_det_unq_up (i, det_unq_up_i) = lap_det_unq_up (i, det_unq_up_i)  +  &
          slater_mat_trans_inv_up (i, j, det_unq_up_i) * ddorb (i, det_unq_orb_lab_srt_up (j, det_unq_up_i))
       enddo
        lap_det_unq_up (i, det_unq_up_i) = lap_det_unq_up (i, det_unq_up_i) * detu (det_unq_up_i)
      enddo
  enddo ! det_unq_up_i

! loop over unique spin-down determinants
  do det_unq_dn_i = 1, ndetdn
      do i = 1, ndn
       do j = 1, ndn
        lap_det_unq_dn (i, det_unq_dn_i) = lap_det_unq_dn (i, det_unq_dn_i)  +  &
          slater_mat_trans_inv_dn (i, j, det_unq_dn_i) * ddorb (nup + i, det_unq_orb_lab_srt_dn (j, det_unq_dn_i))
       enddo
        lap_det_unq_dn (i, det_unq_dn_i) = lap_det_unq_dn (i, det_unq_dn_i) * detd (det_unq_dn_i)
      enddo
  enddo ! det_unq_dn_i

  end subroutine lap_det_unq_bld

! ==============================================================================
  subroutine grd_psid_over_psid_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of determinant part of the wave function psid over psid
! Description   : (Gradient psid)/psid
!
! Created       : J. Toulouse, 05 Dec 2005
! Modified      : J. Toulouse, 04 Jul 2007: unique determinants
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_i, det_unq_up_i, det_unq_dn_i
  integer dim_i
  integer elec_i, elec_up_i, elec_dn_i
  integer csf_i, det_in_csf_i
  real(dp) coefficient

! header
  if (header_exe) then

   call object_create ('grd_psid_over_psid')

   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('csf_coef')
   call object_needed ('cdet_in_csf')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('grd_det_unq_up')
   call object_needed ('grd_det_unq_dn')
   call object_needed ('psido')


   return

  endif

! begin

! allocations
  call object_alloc ('grd_psid_over_psid', grd_psid_over_psid, ndim, nelec)

  grd_psid_over_psid = 0.d0

  do csf_i = 1, ncsf
   do det_in_csf_i = 1, ndet_in_csf (csf_i)

    det_i = iwdet_in_csf (det_in_csf_i, csf_i)
    det_unq_up_i = det_to_det_unq_up (det_i)
    det_unq_dn_i = det_to_det_unq_dn (det_i)

    coefficient = csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)

    do dim_i = 1, ndim

    elec_i = 0
    do elec_up_i = 1, nup
     elec_i = elec_i + 1
     grd_psid_over_psid (dim_i, elec_i) = grd_psid_over_psid (dim_i, elec_i) + coefficient * grd_det_unq_up (dim_i, elec_up_i, det_unq_up_i) * detd (det_unq_dn_i)
    enddo

    do elec_dn_i = 1, ndn
     elec_i = elec_i + 1
     grd_psid_over_psid (dim_i, elec_i) = grd_psid_over_psid (dim_i, elec_i) + coefficient * detu (det_unq_up_i) * grd_det_unq_dn (dim_i, elec_dn_i, det_unq_dn_i)
    enddo

    enddo ! dim_i

   enddo ! det_in_csf_i
  enddo ! csf_i

  grd_psid_over_psid = grd_psid_over_psid / psido

 end subroutine grd_psid_over_psid_bld

! ==============================================================================
  subroutine lap_psid_over_psid_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of determinant part of the wave function psid over psid
! Description   : (Laplacian psid)/psid
!
! Created       : J. Toulouse, 05 Dec 2005
! Modified      : J. Toulouse, 04 Jul 2007: unique determinants
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_i, det_unq_up_i, det_unq_dn_i
  integer elec_i, elec_up_i, elec_dn_i
  integer csf_i, det_in_csf_i
  real(dp) coefficient

! header
  if (header_exe) then

   call object_create ('lap_psid_over_psid')

   call object_needed ('nelec')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('csf_coef')
   call object_needed ('cdet_in_csf')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('lap_det_unq_up')
   call object_needed ('lap_det_unq_dn')
   call object_needed ('psido')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_psid_over_psid', lap_psid_over_psid, nelec)

  lap_psid_over_psid = 0.d0

  do csf_i = 1, ncsf
   do det_in_csf_i = 1, ndet_in_csf (csf_i)

    det_i = iwdet_in_csf (det_in_csf_i, csf_i)
    det_unq_up_i = det_to_det_unq_up (det_i)
    det_unq_dn_i = det_to_det_unq_dn (det_i)

    coefficient = csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)

    elec_i = 0
    do elec_up_i = 1, nup
     elec_i = elec_i + 1
     lap_psid_over_psid (elec_i) = lap_psid_over_psid (elec_i) + coefficient * lap_det_unq_up (elec_up_i, det_unq_up_i) * detd (det_unq_dn_i)
    enddo

    do elec_dn_i = 1, ndn
     elec_i = elec_i + 1
     lap_psid_over_psid (elec_i) = lap_psid_over_psid (elec_i) + coefficient * detu (det_unq_up_i) * lap_det_unq_dn (elec_dn_i, det_unq_dn_i)
    enddo


   enddo ! det_in_csf_i
  enddo ! csf_i

  lap_psid_over_psid = lap_psid_over_psid / psido

 end subroutine lap_psid_over_psid_bld

! ==============================================================================
  subroutine lap_lnpsid_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian ( ln (psid) )
!
! Created       : J. Toulouse, 05 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i
  integer elec_i
  real(dp) sum

! header
  if (header_exe) then

   call object_create ('lap_lnpsid')

   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('lap_psid_over_psid')
   call object_needed ('grd_psid_over_psid')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_lnpsid', lap_lnpsid, nelec)

  do elec_i = 1, nelec
    sum = 0
    do dim_i = 1, ndim
     sum = sum +  grd_psid_over_psid (dim_i, elec_i)**2
    enddo
     lap_lnpsid (elec_i) = lap_psid_over_psid (elec_i) - sum
  enddo

 end subroutine lap_lnpsid_bld

! ==============================================================================
  subroutine sum_lap_lnpsid_bld
! ------------------------------------------------------------------------------
! Description   : Sum_i lap_lnpsid (i)
!
! Created       : J. Toulouse, 05 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer elec_i

! header
  if (header_exe) then

   call object_create ('sum_lap_lnpsid')

   call object_needed ('nelec')
   call object_needed ('lap_lnpsid')

   return

  endif

! begin

! allocations
  call object_associate ('sum_lap_lnpsid', sum_lap_lnpsid)

  sum_lap_lnpsid = 0.d0

  do elec_i = 1, nelec
     sum_lap_lnpsid = sum_lap_lnpsid + lap_lnpsid (elec_i)
  enddo

 end subroutine sum_lap_lnpsid_bld

! ==============================================================================
  subroutine grd_psi_over_psi_old_bld
! ------------------------------------------------------------------------------
! Description   : (Gradient psi)/psi
!
! Created       : J. Toulouse, 05 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i
  integer elec_i

! header
  if (header_exe) then

   call object_create ('grd_psi_over_psi_old')

   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('grd_psid_over_psid')
   call object_needed ('vj')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_psi_over_psi_old', grd_psi_over_psi_old, ndim, nelec)

  do dim_i = 1, ndim
    do elec_i = 1, nelec
      grd_psi_over_psi_old (dim_i, elec_i) = grd_psid_over_psid (dim_i, elec_i) + vj (dim_i, elec_i)
    enddo
  enddo

 end subroutine grd_psi_over_psi_old_bld

! ==============================================================================
  subroutine grd_psi_over_psi_wlk_bld
! ------------------------------------------------------------------------------
! Description   : (Gradient psi)/psi
!
! Created       : J. Toulouse, 06 Mar 2006
! Modified      : J. Toulouse, 17 Nov 2006, walkers
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i
  integer elec_i

! header
  if (header_exe) then

   call object_create ('grd_psi_over_psi_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('ndim')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_psi_over_psi_wlk', grd_psi_over_psi_wlk, ndim, nelec, nwalk)

  if (index(mode, 'vmc') /= 0) then
   call object_provide_by_index (grd_psi_over_psi_wlk_bld_index, vold_index)
   grd_psi_over_psi_wlk (:,:,1) = vold (1:ndim, 1:nelec)

  elseif (index(mode, 'dmc') /= 0) then
   call object_provide_by_index (grd_psi_over_psi_wlk_bld_index, voldw_index)
   grd_psi_over_psi_wlk (:,:,:) = voldw (1:ndim, 1:nelec, 1:nwalk, 1)

  else
   write(6,'(4a)') trim(here),': mode=',trim(mode),' should contain either vmc or dmc.'
   call die (here)

  endif

! call object_provide  ('grd_psi_over_psi_old')
! call is_equal_or_die (grd_psi_over_psi, grd_psi_over_psi_old, 10.d-10, .true.)

 end subroutine grd_psi_over_psi_wlk_bld

! ==============================================================================
  subroutine grd_psi_over_psi_sq_wlk_bld
! ------------------------------------------------------------------------------
! Description   : [(Gradient Psi)/Psi]^2
!
! Created       : J. Toulouse, 06 Mar 2006
! Modified      : J. Toulouse, 17 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, elec_i, walk_i

! header
  if (header_exe) then

   call object_create ('grd_psi_over_psi_sq_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('grd_psi_over_psi_wlk')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_psi_over_psi_sq_wlk', grd_psi_over_psi_sq_wlk, nelec, nwalk)

  grd_psi_over_psi_sq_wlk (:,:) = 0.d0

  do walk_i = 1, nwalk
   do elec_i = 1, nelec
    do dim_i = 1, ndim
      grd_psi_over_psi_sq_wlk (elec_i, walk_i) = grd_psi_over_psi_sq_wlk (elec_i, walk_i) + grd_psi_over_psi_wlk (dim_i, elec_i, walk_i)**2
    enddo
   enddo
  enddo

 end subroutine grd_psi_over_psi_sq_wlk_bld

! ==============================================================================
  subroutine sum_lap_lnpsi_bld
! ------------------------------------------------------------------------------
! Description   : Sum_i lap_lnpsi (i)
!
! Created       : J. Toulouse, 05 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('sum_lap_lnpsi')

   call object_needed ('sum_lap_lnpsid')
   call object_needed ('sum_lap_lnj')

   return

  endif

! begin

! allocations
  call object_associate ('sum_lap_lnpsi', sum_lap_lnpsi)

  sum_lap_lnpsi = sum_lap_lnpsid + sum_lap_lnj

 end subroutine sum_lap_lnpsi_bld

! ==============================================================================
  subroutine div_grd_psi_over_psi_wlk_bld
! ------------------------------------------------------------------------------
! Description   : div ( (grad Psi)/Psi )
!
! Created       : J. Toulouse, 06 Mar 2006
! Modified      : J. Toulouse, 17 Nov 2006, walkers
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('div_grd_psi_over_psi_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')

   return

  endif

! begin
! allocations
  call object_alloc ('div_grd_psi_over_psi_wlk', div_grd_psi_over_psi_wlk, nelec, nwalk)

  if (index(mode, 'vmc') /= 0) then
   call object_provide_by_index (div_grd_psi_over_psi_wlk_bld_index, div_vo_index)
   div_grd_psi_over_psi_wlk (:,1) = div_vo (1:nelec)

  elseif (index(mode, 'dmc') /= 0) then
   call object_provide_by_index (div_grd_psi_over_psi_wlk_bld_index, div_vow_index)
   div_grd_psi_over_psi_wlk (:,:) = div_vow (1:nelec, 1:nwalk)

  else
   write(6,'(4a)') trim(here),': mode=',trim(mode),' should contain either vmc or dmc.'
   call die (here)

  endif

 end subroutine div_grd_psi_over_psi_wlk_bld

! ==============================================================================
  subroutine lap_psi_over_psi_wlk_bld
! ------------------------------------------------------------------------------
! Description   : (Laplacian Psi) / Psi
!
! Created       : J. Toulouse, 06 Mar 2006
! Modified      : J. Toulouse, 17 Nov 2006, walkers
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('lap_psi_over_psi_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('div_grd_psi_over_psi_wlk')
   call object_needed ('grd_psi_over_psi_sq_wlk')

   return

  endif

! begin
! allocations
  call object_alloc ('lap_psi_over_psi_wlk', lap_psi_over_psi_wlk, nelec, nwalk)

  lap_psi_over_psi_wlk (:,:) = div_grd_psi_over_psi_wlk (:,:) + grd_psi_over_psi_sq_wlk (:,:)

!  call object_provide ('lap_psid_over_psid')
!  call is_equal_or_die (lap_psi_over_psi, lap_psid_over_psid, 10.d-10, .true.) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 end subroutine lap_psi_over_psi_wlk_bld

! ==============================================================================
  subroutine maximize_psi2
! ------------------------------------------------------------------------------
! Description   : find the maximum of the wave function square psi2
!
! Created       : J. Toulouse, 06 Apr 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'find_max_of_psi'
  real(dp) psi2, psi2_tol, psi2_best
  integer elec_dim_nb, vectex_nb, vertex_i, elec_i, trial_i, trial_best_i
  integer itmax, iter, iter_all
  real(dp), allocatable  :: coord_elec_vec (:), coord_elec_vec_best (:)
  real(dp), allocatable  :: coord_elec_vec_simplex (:,:)
  real(dp), allocatable  :: psi2_simplex (:)
  logical converged 

! begin
  write(6,*) 
  write(6,'(a)') 'Search for maximum of the wave function...'

! initialization
  psi2_best = 0.d0

! needed objects
  call object_provide ('ndim')
  call object_provide ('nelec')

! allocations
  elec_dim_nb = ndim * nelec
  vectex_nb = elec_dim_nb + 1
  call alloc ('coord_elec_vec', coord_elec_vec, elec_dim_nb)
  call alloc ('coord_elec_vec_best', coord_elec_vec_best, elec_dim_nb)
  call alloc ('coord_elec_vec_simplex', coord_elec_vec_simplex, vectex_nb, elec_dim_nb)
  call alloc ('psi2_simplex', psi2_simplex, vectex_nb)

  write(6,'(a,i5)') 'Number of trial of maximization = ', max_psi2_trial_nb
  do trial_i = 1, max_psi2_trial_nb

!  initialize simplex
   write(6,*)
   write(6,'(a,i3)') 'maximization trial # ', trial_i
   write(6,'(a)') 'initialize simplex:'
   do vertex_i = 1, vectex_nb
     write(6,'(a,i3)') 'vertex # ',vertex_i
     call mc_configs_read
     do elec_i = 1, nelec
       write(6,'(a,i3, 3f)') 'electron # ',elec_i, xold (1:ndim, elec_i)
     enddo ! elec_i
     call flatten  (coord_elec_vec, xold, ndim, nelec)
     coord_elec_vec_simplex (vertex_i, :) = coord_elec_vec (:)
     psi2_simplex (vertex_i) = psi2_eval (coord_elec_vec)
     write(6,'(a,f)') 'wave function square = ', psi2_simplex (vertex_i)
   enddo ! vertex_i
   
!  simplex algorithm
   write(6,*)
   write(6,'(a)') 'Start maximization by simplex algorithm...'
   converged  = .false.
   itmax = 100
   iter_all = 0
   psi2_tol = 1.d-8
   write(6,'(a,f)') 'convergence threshold on wave function square = ',psi2_tol
   do
     call amoeba(coord_elec_vec_simplex,psi2_simplex,vectex_nb,elec_dim_nb,elec_dim_nb,psi2_tol,minus_psi2_eval,iter,itmax,converged)
   
     iter_all = iter_all + iter
  
     write(6,*)
     write(6,'(a,i10)') 'Iteration # ',iter_all
   
     coord_elec_vec = coord_elec_vec_simplex (1,:)
     call unflatten  (coord_elec_vec, xold, ndim, nelec)
     do elec_i = 1, nelec
       write(6,'(a,i3, 3f)') 'electron # ',elec_i, xold (1:ndim, elec_i)
     enddo ! elec_i
     psi2 = psi2_eval (coord_elec_vec)
     write(6,'(a,f)') 'wave function square = ', psi2
  

     if (converged) exit
   enddo

   if (psi2 > psi2_best) then
       trial_best_i = trial_i
       coord_elec_vec_best (:) = coord_elec_vec (:)
       psi2_best = psi2
       write(6,'(a)') 'convergence reached, best so far.'
   else
       write(6,'(a)') 'convergence reached.'
   endif

  enddo ! trial_i

  write(6,*)
  write(6,'(a)') 'Final results:'
  write(6,'(a,i10)') 'the best maximum was found at maximization trial # ', trial_best_i
  call unflatten  (coord_elec_vec_best, xold, ndim, nelec)
  do elec_i = 1, nelec
     write(6,'(a,i3, 3f)') 'electron # ',elec_i, xold (1:ndim, elec_i)
  enddo ! elec_i
  write(6,'(a,f)') 'wave function square = ', psi2_best

 end subroutine maximize_psi2

! ==============================================================================
  function psi2_eval (coord_elec_vec)
! ------------------------------------------------------------------------------
! Description   : returns value of wave function square  for electron coordinates coord_elec
!
! Created       : J. Toulouse, 06 Apr 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  real(dp), intent(in) :: coord_elec_vec (ndim*nelec)

! output
  real(dp) :: psi2_eval

! begin
  call unflatten  (coord_elec_vec, xold, ndim, nelec)
  call hpsi (xold,psido,psijo,vold,div_vo,d2o,peo,peio,eold(1),denergy,1)
  psi2_eval = (dexp(psijo)*psido)**2

!  write(6,*) 'eold=',eold(1)

  end function psi2_eval

! ==============================================================================
  function minus_psi2_eval (coord_elec_vec)
! ------------------------------------------------------------------------------
! Description   : returns  (- psi2) for minimization
!
! Created       : J. Toulouse, 08 Apr 2008
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! input
  real(dp), intent(in) :: coord_elec_vec (ndim*nelec)

! output
  real(dp) :: minus_psi2_eval

  minus_psi2_eval = - psi2_eval (coord_elec_vec)

  end function minus_psi2_eval

end module psi_mod
