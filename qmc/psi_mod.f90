module psi_mod

  use all_tools_mod
  use orbitals_mod
  use montecarlo_mod
  use determinants_mod

! Declaration of global variables and default values
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
   call object_provide_in_node_by_index (grd_psi_over_psi_wlk_bld_index, vold_index)
   grd_psi_over_psi_wlk (:,:,1) = vold (1:ndim, 1:nelec)

  elseif (index(mode, 'dmc') /= 0) then
   call object_provide_in_node_by_index (grd_psi_over_psi_wlk_bld_index, voldw_index)
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
   call object_provide_in_node_by_index (div_grd_psi_over_psi_wlk_bld_index, div_vo_index)
   div_grd_psi_over_psi_wlk (:,1) = div_vo (1:nelec)

  elseif (index(mode, 'dmc') /= 0) then
   call object_provide_in_node_by_index (div_grd_psi_over_psi_wlk_bld_index, div_vow_index)
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

end module psi_mod
