module csfs_mod

  use all_tools_mod
  use orbitals_mod
  use determinants_mod

! declaration of global variables and default values
  type (type_integer_row), allocatable :: det_unq_up_in_csf (:)
  type (type_integer_row), allocatable :: det_unq_dn_in_csf (:)
  type (type_real_row),    allocatable :: cdet_unq_in_csf (:)
  real(dp), allocatable                :: csf_prefac (:)

  real(dp)                  :: wfdet_ovlp
  real(dp), allocatable     :: dens_mat_wfdet_up (:,:)
  real(dp), allocatable     :: dens_mat_wfdet_dn (:,:)
  real(dp), allocatable     :: dens_mat_wfdet (:,:)

  real(dp), allocatable     :: csfs_ovlp (:,:)
  real(dp), allocatable     :: csfs_wfdet_ovlp (:)

  contains

! ==============================================================================
  subroutine det_unq_in_csf_bld
! ------------------------------------------------------------------------------
! Description : each csf as a sorted list of unique determinants
!
! Created     : J. Toulouse, 14 Dec 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer csf_i, det_i, det_in_csf_i
  integer det_unq_up_cur, det_unq_dn_cur
  real(dp) cdet_unq_cur

! header
  if (header_exe) then

   call object_create ('det_unq_in_csf_nb')
   call object_create ('det_unq_up_in_csf')
   call object_create ('det_unq_dn_in_csf')
   call object_create ('cdet_unq_in_csf')
   call object_create ('csf_prefac')

   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('det_unq_up_in_csf', det_unq_up_in_csf, ncsf)
  call object_alloc ('det_unq_dn_in_csf', det_unq_dn_in_csf, ncsf)
  call object_alloc ('cdet_unq_in_csf', cdet_unq_in_csf, ncsf)
  call object_alloc ('csf_prefac', csf_prefac, ncsf)

  do csf_i = 1, ncsf

   do det_in_csf_i = 1, ndet_in_csf (csf_i)

!    index of original determinant
     det_i = iwdet_in_csf (det_in_csf_i, csf_i)

!    indexes of corresponding unique spin-up and spin-dn determinants
     det_unq_up_cur = det_to_det_unq_up (det_i)
     det_unq_dn_cur = det_to_det_unq_dn (det_i)

!    coefficients of unique determinants, with sign
     cdet_unq_cur = cdet_in_csf (det_in_csf_i, csf_i)

     call append (det_unq_up_in_csf (csf_i)%row, det_unq_up_cur)
     call append (det_unq_dn_in_csf (csf_i)%row, det_unq_dn_cur)
     call append (cdet_unq_in_csf (csf_i)%row, cdet_unq_cur)

   enddo ! det_in_csf_i

!  sort determinants in csf
   call sort (det_unq_up_in_csf (csf_i)%row (:), det_unq_dn_in_csf (csf_i)%row (:), cdet_unq_in_csf (csf_i)%row (:))

!  prefactor of csf = coefficient of first determinant
   csf_prefac (csf_i) = cdet_unq_in_csf (csf_i)%row (1)

!  normalize csf so that the coefficient of first determinant is 1
   cdet_unq_in_csf (csf_i)%row (:) = cdet_unq_in_csf (csf_i)%row (:) / csf_prefac (csf_i)

  enddo ! csf_i

!  write(6,'(2a)') trim(here), ': csfs in terms of sorted unique determinants:'
!  do csf_i = 1, ncsf
!   write(6,'(2a,i4,a)') trim(here), ': csf # ', csf_i,' :'
!   write(6,'(2a,100i3)') trim(here), ': det_unq_up_in_csf =', det_unq_up_in_csf (csf_i)%row (:)
!   write(6,'(2a,100i3)') trim(here), ': det_unq_dn_in_csf =', det_unq_dn_in_csf (csf_i)%row (:)
!   write(6,'(2a,100f7.3)')  trim(here), ': cdet_unq_in_csf =', cdet_unq_in_csf (csf_i)%row (:)
!  enddo

  end subroutine det_unq_in_csf_bld

! ==============================================================================
  subroutine wfdet_ovlp_bld_old
! ------------------------------------------------------------------------------
! Description   : overlap of the determinantal part of the wave function
! Description   : assuming orthnomality of orbitals
!
! Created       : J. Toulouse, 29 Nov 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer csf_i, csf_j
  integer det_i, det_j, det_unq_up_i, det_unq_dn_i, det_unq_up_j, det_unq_dn_j
  integer det_in_csf_i, det_in_csf_j

! header
  if (header_exe) then

   call object_create ('wfdet_ovlp')

   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('csf_coef')
   call object_needed ('orb_occ_in_det_unq_up')
   call object_needed ('orb_occ_in_det_unq_dn')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_associate ('wfdet_ovlp', wfdet_ovlp)

  wfdet_ovlp = 0.d0

  do csf_i = 1, ncsf
   do det_in_csf_i = 1, ndet_in_csf (csf_i)
    det_i = iwdet_in_csf (det_in_csf_i, csf_i)
    det_unq_up_i = det_to_det_unq_up (det_i)
    det_unq_dn_i = det_to_det_unq_dn (det_i)

     do csf_j = 1, ncsf
      do det_in_csf_j = 1, ndet_in_csf (csf_j)
       det_j = iwdet_in_csf (det_in_csf_j, csf_j)
       det_unq_up_j = det_to_det_unq_up (det_j)
       det_unq_dn_j = det_to_det_unq_dn (det_j)

        if (arrays_equal (orb_occ_in_det_unq_up (:, det_unq_up_i), orb_occ_in_det_unq_up (:, det_unq_up_j)) .and. &
            arrays_equal (orb_occ_in_det_unq_dn (:, det_unq_dn_i), orb_occ_in_det_unq_dn (:, det_unq_dn_j))) then

          wfdet_ovlp = wfdet_ovlp + csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i) * csf_coef (csf_j, 1) * cdet_in_csf (det_in_csf_j, csf_j)

        endif

       enddo ! det_in_csf_j
      enddo ! csf_j
    enddo ! det_in_csf_i
   enddo ! csf_i

!   write(6,*) trim(here),': wfdet_ovlp=',wfdet_ovlp

  end subroutine wfdet_ovlp_bld_old

! ==============================================================================
  subroutine csfs_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap of the CSFs  (without Jastrow)
! Description   : assuming orthonormality of orbitals
!
! Created       : J. Toulouse, 31 Mar 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer csf_i, csf_j
  integer det_i, det_j, det_unq_up_i, det_unq_dn_i, det_unq_up_j, det_unq_dn_j
  integer det_in_csf_i, det_in_csf_j

! header
  if (header_exe) then

   call object_create ('csfs_ovlp')

   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('orb_occ_in_det_unq_up')
   call object_needed ('orb_occ_in_det_unq_dn')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('csfs_ovlp', csfs_ovlp, ncsf, ncsf)

  csfs_ovlp (:,:) = 0.d0

  do csf_i = 1, ncsf
   do det_in_csf_i = 1, ndet_in_csf (csf_i)
    det_i = iwdet_in_csf (det_in_csf_i, csf_i)
    det_unq_up_i = det_to_det_unq_up (det_i)
    det_unq_dn_i = det_to_det_unq_dn (det_i)

     do csf_j = 1, ncsf
      do det_in_csf_j = 1, ndet_in_csf (csf_j)
       det_j = iwdet_in_csf (det_in_csf_j, csf_j)
       det_unq_up_j = det_to_det_unq_up (det_j)
       det_unq_dn_j = det_to_det_unq_dn (det_j)

        if (arrays_equal (orb_occ_in_det_unq_up (:, det_unq_up_i), orb_occ_in_det_unq_up (:, det_unq_up_j)) .and. &
            arrays_equal (orb_occ_in_det_unq_dn (:, det_unq_dn_i), orb_occ_in_det_unq_dn (:, det_unq_dn_j))) then

          csfs_ovlp (csf_i, csf_j) = csfs_ovlp (csf_i, csf_j) + cdet_in_csf (det_in_csf_i, csf_i) * cdet_in_csf (det_in_csf_j, csf_j)

        endif

       enddo ! det_in_csf_j
      enddo ! csf_j
    enddo ! det_in_csf_i
   enddo ! csf_i


  end subroutine csfs_ovlp_bld

! ==============================================================================
  subroutine csfs_wfdet_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap of the CSFs  with the multi-configuration wave function (without Jastrow)
! Description   : assuming orthonormality of orbitals
!
! Created       : J. Toulouse, 31 Mar 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer csf_i, csf_j

! header
  if (header_exe) then

   call object_create ('csfs_wfdet_ovlp')

   call object_needed ('ncsf')
   call object_needed ('csfs_ovlp')
   call object_needed ('csf_coef')

   return

  endif

! begin

! allocations
  call object_alloc ('csfs_wfdet_ovlp', csfs_wfdet_ovlp, ncsf)

  csfs_wfdet_ovlp (:) = 0.d0

  do csf_i = 1, ncsf
     do csf_j = 1, ncsf
       csfs_wfdet_ovlp (csf_i) = csfs_wfdet_ovlp (csf_i) + csf_coef (csf_j, 1) * csfs_ovlp (csf_i, csf_j)
     enddo ! csf_j
  enddo ! csf_i


  end subroutine csfs_wfdet_ovlp_bld

! ==============================================================================
  subroutine wfdet_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : overlap of the determinantal part of the wave function
! Description   : assuming orthnomality of orbitals
!
! Created       : J. Toulouse, 31 Mar 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer csf_i

! header
  if (header_exe) then

   call object_create ('wfdet_ovlp')

   call object_needed ('ncsf')
   call object_needed ('csfs_wfdet_ovlp')
   call object_needed ('csf_coef')

   return

  endif

! begin

! allocations
  call object_associate ('wfdet_ovlp', wfdet_ovlp)

  wfdet_ovlp = 0.d0

  do csf_i = 1, ncsf
      wfdet_ovlp = wfdet_ovlp + csf_coef (csf_i, 1) * csfs_wfdet_ovlp (csf_i)
  enddo ! csf_i

  end subroutine wfdet_ovlp_bld

! ==============================================================================
  subroutine dens_mat_wfdet_bld
! ------------------------------------------------------------------------------
! Description   : density matrix of determinantal part of wave function
! Description   : assuming orthonormality of orbitals
!
! Created       : J. Toulouse, 29 Nov 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_i, orb_j
  integer csf_i, csf_j
  integer det_i, det_unq_up_i, det_unq_dn_i
  integer det_j, det_unq_up_j, det_unq_dn_j
  integer det_in_csf_i, det_in_csf_j
  real(dp) dens_mat_wfdet_up_trace
  real(dp) dens_mat_wfdet_dn_trace
  real(dp) dens_mat_wfdet_trace

! header
  if (header_exe) then

   call object_create ('dens_mat_wfdet_up')
   call object_create ('dens_mat_wfdet_dn')
   call object_create ('dens_mat_wfdet')

   call object_needed ('nelec')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('orb_tot_nb')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('csf_coef')
   call object_needed ('orb_occ_in_adet_unq_up')
   call object_needed ('orb_occ_in_adet_unq_dn')
   call object_needed ('sgn_adet_unq_dn')
   call object_needed ('sgn_adet_unq_up')
   call object_needed ('wfdet_ovlp')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('dens_mat_wfdet_up', dens_mat_wfdet_up, orb_tot_nb, orb_tot_nb)
  call object_alloc ('dens_mat_wfdet_dn', dens_mat_wfdet_dn, orb_tot_nb, orb_tot_nb)
  call object_alloc ('dens_mat_wfdet', dens_mat_wfdet, orb_tot_nb, orb_tot_nb)

  dens_mat_wfdet_up = 0.d0
  dens_mat_wfdet_dn = 0.d0

  do orb_i = 1, orb_tot_nb
  do orb_j = 1, orb_tot_nb

  do csf_i = 1, ncsf
   do det_in_csf_i = 1, ndet_in_csf (csf_i)
    det_i = iwdet_in_csf (det_in_csf_i, csf_i)
    det_unq_up_i = det_to_det_unq_up (det_i)
    det_unq_dn_i = det_to_det_unq_dn (det_i)

     do csf_j = 1, ncsf
      do det_in_csf_j = 1, ndet_in_csf (csf_j)
       det_j = iwdet_in_csf (det_in_csf_j, csf_j)
       det_unq_up_j = det_to_det_unq_up (det_j)
       det_unq_dn_j = det_to_det_unq_dn (det_j)

!       spin-up contribution
        if (arrays_equal (orb_occ_in_adet_unq_up (:, orb_i, det_unq_up_i), orb_occ_in_adet_unq_up (:, orb_j, det_unq_up_j)) .and. &
            arrays_equal (orb_occ_in_det_unq_dn (:, det_unq_dn_i), orb_occ_in_det_unq_dn (:, det_unq_dn_j))) then

          dens_mat_wfdet_up (orb_i, orb_j) = dens_mat_wfdet_up (orb_i, orb_j) +                        &
                                          sgn_adet_unq_up (orb_i, det_unq_up_i) * sgn_adet_unq_up (orb_j, det_unq_up_j)      &
                                        * csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)      &
                                        * csf_coef (csf_j, 1) * cdet_in_csf (det_in_csf_j, csf_j)

        endif

!       spin-down contribution
        if (arrays_equal (orb_occ_in_adet_unq_dn (:, orb_i, det_unq_dn_i), orb_occ_in_adet_unq_dn (:, orb_j, det_unq_dn_j)) .and. &
            arrays_equal (orb_occ_in_det_unq_up (:, det_unq_up_i), orb_occ_in_det_unq_up (:, det_unq_up_j))) then

          dens_mat_wfdet_dn (orb_i, orb_j) = dens_mat_wfdet_dn (orb_i, orb_j) +                       &
                                          sgn_adet_unq_dn (orb_i, det_unq_dn_i) * sgn_adet_unq_dn (orb_j, det_unq_dn_j)     &
                                        * csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)     &
                                        * csf_coef (csf_j, 1) * cdet_in_csf (det_in_csf_j, csf_j)

        endif

       enddo ! det_in_csf_j
      enddo ! csf_j
    enddo ! det_in_csf_i
   enddo ! csf_i

   enddo ! orb_j
   enddo ! orb_i

   dens_mat_wfdet_up = dens_mat_wfdet_up/wfdet_ovlp
   dens_mat_wfdet_dn = dens_mat_wfdet_dn/wfdet_ovlp
   dens_mat_wfdet = dens_mat_wfdet_up + dens_mat_wfdet_dn

   write(6,*) trim(here),': dens_mat_wfdet=',dens_mat_wfdet

!  check trace up
   dens_mat_wfdet_up_trace = 0.d0
   do orb_i = 1, orb_tot_nb
    dens_mat_wfdet_up_trace = dens_mat_wfdet_up_trace + dens_mat_wfdet_up (orb_i, orb_i)
   enddo
   write(6,*) trim(here),': dens_mat_wfdet_up_trace=',dens_mat_wfdet_up_trace

   if (dabs(dens_mat_wfdet_up_trace-nup) > 10d-6) then
    write(6,*) trim(here),': dens_mat_wfdet_up_trace=',dens_mat_wfdet_up_trace,' /= nup=',nup
    call die (here)
   endif

!  check trace down
   dens_mat_wfdet_dn_trace = 0.d0
   do orb_i = 1, orb_tot_nb
    dens_mat_wfdet_dn_trace = dens_mat_wfdet_dn_trace + dens_mat_wfdet_dn (orb_i, orb_i)
   enddo
   write(6,*) trim(here),': dens_mat_wfdet_dn_trace=',dens_mat_wfdet_dn_trace

   if (dabs(dens_mat_wfdet_dn_trace-ndn) > 10d-6) then
    write(6,*) trim(here),': dens_mat_wfdet_dn_trace=',dens_mat_wfdet_dn_trace,' /= ndn=',ndn
    call die (here)
   endif

!  check trace
   dens_mat_wfdet_trace = 0.d0
   do orb_i = 1, orb_tot_nb
    dens_mat_wfdet_trace = dens_mat_wfdet_trace + dens_mat_wfdet (orb_i, orb_i)
   enddo
   write(6,*) trim(here),': dens_mat_wfdet_trace=',dens_mat_wfdet_trace

   if (dabs(dens_mat_wfdet_trace-nelec) > 10d-6) then
    write(6,*) trim(here),': dens_mat_wfdet_trace=',dens_mat_wfdet_trace,' /= nelec=',nelec
    call die (here)
   endif

  end subroutine dens_mat_wfdet_bld

end module csfs_mod
