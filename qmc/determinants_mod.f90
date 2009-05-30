module determinants_mod

  use all_tools_mod
  use orbitals_mod
  use basis_mod

! Declaration of global variables and default values
!  integer, allocatable              :: det_unq_orb_lab_up (:,:)
!  integer, allocatable              :: det_unq_orb_lab_dn (:,:)
!  integer, allocatable              :: det_orb_lab_srt_up (:,:)
!  integer, allocatable              :: det_orb_lab_srt_dn (:,:)
!  integer, allocatable              :: det_orb_lab_srt_sgn_up (:)
!  integer, allocatable              :: det_orb_lab_srt_sgn_dn (:)

!  integer                           :: det_unq_up_nb = 0
!  integer                           :: det_unq_dn_nb = 0
!  integer, allocatable              :: det_unq_orb_lab_srt_up (:,:)
!  integer, allocatable              :: det_unq_orb_lab_srt_dn (:,:)
!  integer, allocatable              :: det_unq_orb_lab_srt_sgn_up (:)
!  integer, allocatable              :: det_unq_orb_lab_srt_sgn_dn (:)
!  integer, allocatable              :: det_to_det_unq_up (:)
!  integer, allocatable              :: det_to_det_unq_dn (:)
!  integer, allocatable              :: det_to_det_unq_sgn_up (:)
!  integer, allocatable              :: det_to_det_unq_sgn_dn (:)
!  integer, allocatable              :: det_unq_up_to_det (:)
!  integer, allocatable              :: det_unq_dn_to_det (:)

  real(dp), allocatable             :: slater_mat_trans_inv_up (:,:,:)
  real(dp), allocatable             :: slater_mat_trans_inv_dn (:,:,:)

  logical, allocatable              :: orb_occ_in_adet_unq_up (:,:,:)
  logical, allocatable              :: orb_occ_in_adet_unq_dn (:,:,:)
  integer, allocatable              :: sgn_adet_unq_up (:,:)
  integer, allocatable              :: sgn_adet_unq_dn (:,:)

  contains

!===========================================================================
  subroutine determinants_rd
!---------------------------------------------------------------------------
! Description : read Slater determinants
!
! Created     : J. Toulouse, 07 Apr 2009
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'orb_coefs_rd'
  character(len=500) line
  integer iostat, elec_i, det_i

! begin
  call object_provide ('nelec')

  det_i = 0

  do
   read(unit_input,'(a)',iostat=iostat) line
!   write(6,*) 'line >',trim(line),'<'

!  no next line found
   if(iostat < 0) then
     call die (lhere, 'error while reading determinants')
   endif

!  convert to lower case
   call upplow (line)

!  exit when 'end' is read
   if (index(line,'end') /= 0) exit

   det_i = det_i + 1
   if (det_i > MDET) then
       call die (lhere, ' det_i='+det_i+' > MDET='+MDET)
   endif

   call alloc ('iworbd', iworbd, nelec, det_i)
   read(line,*,iostat=iostat) (iworbd(elec_i,det_i),elec_i=1,nelec)
   if(iostat < 0) then
     call die (lhere, 'error while reading determinants')
   endif
  enddo

  ndet = det_i
  call object_modified ('ndet')
!!!  call object_modified ('iworbd') iworbd is sorted after reading csfs

  end subroutine determinants_rd

!! ==============================================================================
!  subroutine det_orb_lab_srt_bld
!! ------------------------------------------------------------------------------
!! Description   : determinants represented as list of orbital labels sorted
!! Description   : plus associated sign
!!
!! Created       : J. Toulouse, 28 Oct 2005
!! ------------------------------------------------------------------------------
!  implicit none
!
!! local
!  integer det_i
!
!! header
!  if (header_exe) then
!
!   call object_create ('det_orb_lab_srt_up')
!   call object_create ('det_orb_lab_srt_dn')
!   call object_create ('det_orb_lab_srt_sgn_up')
!   call object_create ('det_orb_lab_srt_sgn_dn')
!
!   call object_needed ('ndet')
!   call object_needed ('nup')
!   call object_needed ('ndn')
!   call object_needed ('det_orb_lab_up')
!   call object_needed ('det_orb_lab_dn')
!
!   return
!
!  endif
!
!! begin
!
!! allocations
!  call object_alloc ('det_orb_lab_srt_up', det_orb_lab_srt_up, nup, ndet)
!  call object_alloc ('det_orb_lab_srt_dn', det_orb_lab_srt_dn, ndn, ndet)
!  call object_alloc ('det_orb_lab_srt_sgn_up', det_orb_lab_srt_sgn_up, ndet)
!  call object_alloc ('det_orb_lab_srt_sgn_dn', det_orb_lab_srt_sgn_dn, ndet)
!
!  det_orb_lab_srt_up = det_orb_lab_up
!  det_orb_lab_srt_dn = det_orb_lab_dn
!
!  do det_i = 1, ndet
!     call sort_and_sign (det_orb_lab_srt_up (:, det_i), det_orb_lab_srt_sgn_up (det_i))
!     call sort_and_sign (det_orb_lab_srt_dn (:, det_i), det_orb_lab_srt_sgn_dn (det_i))
!  enddo ! det_i
!
!! Imposing that the determinants must be given already sorted to avoid taking care of permutations elsewehere
!  do det_i = 1, ndet
!     if (.not. arrays_equal (det_orb_lab_up (:, det_i), det_orb_lab_srt_up (:, det_i))) then
!        write(6,'(2a,i,a)') trim(here),' spin-up determinant # ', det_i,' is not sorted:'
!        write(6,'(2a,100i)') trim(here),':', det_orb_lab_up (:, det_i)
!        call die (here)
!     endif
!     if (.not. arrays_equal (det_orb_lab_dn (:, det_i), det_orb_lab_srt_dn (:, det_i))) then
!        write(6,'(2a,i,a)') trim(here),' spin-up determinant # ', det_i,' is not sorted:'
!        write(6,'(2a,100i)') trim(here),':', det_orb_lab_dn (:, det_i)
!        call die (here)
!     endif
!  enddo
!
!  end subroutine det_orb_lab_srt_bld

!! ==============================================================================
!  subroutine det_unq_orb_lab_srt_bld
!! ------------------------------------------------------------------------------
!! Description   : unique determinants represented as list of orbital labels sorted
!! Description   : plus associated sign
!!
!! Created       : J. Toulouse, 14 Dec 2006
!! ------------------------------------------------------------------------------
!  implicit none
!
!! local
!  integer det_i
!  integer det_unq_cur_up, det_unq_sgn_cur_up, det_unq_up_k
!  integer det_unq_cur_dn, det_unq_sgn_cur_dn, det_unq_dn_k
!
!! header
!  if (header_exe) then
!
!   call object_create ('det_unq_up_nb')
!   call object_create ('det_unq_dn_nb')
!   call object_create ('det_unq_orb_lab_srt_up')
!   call object_create ('det_unq_orb_lab_srt_dn')
!   call object_create ('det_unq_orb_lab_srt_sgn_up')
!   call object_create ('det_unq_orb_lab_srt_sgn_dn')
!   call object_create ('det_to_det_unq_up')
!   call object_create ('det_to_det_unq_dn')
!   call object_create ('det_to_det_unq_sgn_up')
!   call object_create ('det_to_det_unq_sgn_dn')
!   call object_create ('det_unq_up_to_det')
!   call object_create ('det_unq_dn_to_det')
!
!   call object_needed ('ndet')
!   call object_needed ('nup')
!   call object_needed ('ndn')
!   call object_needed ('det_orb_lab_srt_up')
!   call object_needed ('det_orb_lab_srt_dn')
!   call object_needed ('det_orb_lab_srt_sgn_up')
!   call object_needed ('det_orb_lab_srt_sgn_dn')
!
!   return
!
!  endif
!
!! begin
!
!! allocations
!  call object_alloc ('det_to_det_unq_up', det_to_det_unq_up, ndet)
!  call object_alloc ('det_to_det_unq_dn', det_to_det_unq_dn, ndet)
!  call object_alloc ('det_to_det_unq_sgn_up', det_to_det_unq_sgn_up, ndet)
!  call object_alloc ('det_to_det_unq_sgn_dn', det_to_det_unq_sgn_dn, ndet)
!
!! initialization
!  det_unq_up_nb = 0
!  det_unq_dn_nb = 0
!
!  do det_i = 1, ndet
!
!!   spin-up determinants
!!   check if current determinant has already been encountered
!    det_unq_cur_up = 0
!    do det_unq_up_k = 1, det_unq_up_nb
!      if (arrays_equal (det_orb_lab_srt_up (:, det_i), det_unq_orb_lab_srt_up (:, det_unq_up_k))) then
!         det_unq_cur_up = det_unq_up_k
!         det_unq_sgn_cur_up = det_orb_lab_srt_sgn_up (det_i) * det_unq_orb_lab_srt_sgn_up (det_unq_up_k)
!         exit
!      endif
!    enddo
!
!!   if current determinant is a new determinant, add it to the list of unique determinants
!    if (det_unq_cur_up == 0) then
!       det_unq_up_nb = det_unq_up_nb + 1
!       det_unq_cur_up = det_unq_up_nb
!       det_unq_sgn_cur_up = 1
!       call object_alloc ('det_unq_orb_lab_srt_up', det_unq_orb_lab_srt_up, nup, det_unq_up_nb)
!       call object_alloc ('det_unq_orb_lab_srt_sgn_up', det_unq_orb_lab_srt_sgn_up, det_unq_up_nb)
!       call object_alloc ('det_unq_up_to_det', det_unq_up_to_det, det_unq_up_nb)
!       det_unq_orb_lab_srt_up (:, det_unq_up_nb) = det_orb_lab_srt_up (:, det_i)
!       det_unq_orb_lab_srt_sgn_up (det_unq_up_nb) = det_orb_lab_srt_sgn_up (det_i)
!
!!      from unique determinants to original determinants
!       det_unq_up_to_det (det_unq_up_nb) = det_i
!
!     endif
!
!!    from original determinants to unique determinants
!     det_to_det_unq_up (det_i) = det_unq_cur_up
!     det_to_det_unq_sgn_up (det_i) = det_unq_sgn_cur_up
!
!
!!   spin-dn determinants
!!   check if current determinant has already been encountered
!    det_unq_cur_dn = 0
!    do det_unq_dn_k = 1, det_unq_dn_nb
!      if (arrays_equal (det_orb_lab_srt_dn (:, det_i), det_unq_orb_lab_srt_dn (:, det_unq_dn_k))) then
!         det_unq_cur_dn = det_unq_dn_k
!         det_unq_sgn_cur_dn = det_orb_lab_srt_sgn_dn (det_i) * det_unq_orb_lab_srt_sgn_dn (det_unq_dn_k)
!         exit
!      endif
!    enddo
!
!!   if current determinant is a new determinant, add it to the list of unique determinants
!    if (det_unq_cur_dn == 0) then
!       det_unq_dn_nb = det_unq_dn_nb + 1
!       det_unq_cur_dn = det_unq_dn_nb
!       det_unq_sgn_cur_dn = 1
!       call object_alloc ('det_unq_orb_lab_srt_dn', det_unq_orb_lab_srt_dn, ndn, det_unq_dn_nb)
!       call object_alloc ('det_unq_orb_lab_srt_sgn_dn', det_unq_orb_lab_srt_sgn_dn, det_unq_dn_nb)
!       call object_alloc ('det_unq_dn_to_det', det_unq_dn_to_det, det_unq_dn_nb)
!       det_unq_orb_lab_srt_dn (:, det_unq_dn_nb) = det_orb_lab_srt_dn (:, det_i)
!       det_unq_orb_lab_srt_sgn_dn (det_unq_dn_nb) =  det_orb_lab_srt_sgn_dn (det_i)
!
!!      from unique determinants to original determinants
!       det_unq_dn_to_det (det_unq_dn_nb) = det_i
!
!     endif
!
!!    from original determinants to unique determinants
!     det_to_det_unq_dn (det_i) = det_unq_cur_dn
!     det_to_det_unq_sgn_dn (det_i) = det_unq_sgn_cur_dn
!
!  enddo ! det_i
!
!!  write(6,'(2a,i4)') trim(here),': number of unique spin-up determinants in ground-state wave function =', det_unq_up_nb
!
!!  do det_unq_up_k = 1, det_unq_up_nb
!!   write(6,'(2a,i3,a,20i3)') trim(here),': determinant unique # ',det_unq_up_k,' : ',det_unq_orb_lab_srt_up (:, det_unq_up_k)
!!  enddo
!
!!  write(6,'(2a,i4)') trim(here),': number of unique spin-dn determinants in ground-state wave function =', det_unq_dn_nb
!
!!  do det_unq_dn_k = 1, det_unq_dn_nb
!!   write(6,'(2a,i3,a,20i3)') trim(here),': determinant unique # ',det_unq_dn_k,' : ',det_unq_orb_lab_srt_dn (:, det_unq_dn_k)
!!  enddo
!
!!  write(6,'(2a)') trim(here),': correspondance with original determinants:'
!
!!  do det_i = 1, ndet
!!   write(6,'(2a,i3,a,i3,a,i3)') trim(here),': orginal determinant ',det_i,' -> unique spin-up and spin-dn determinants ',det_to_det_unq_up (det_i) * det_to_det_unq_sgn_up (det_i), ' |', det_to_det_unq_dn (det_i) * det_to_det_unq_sgn_dn (det_i)
!!  enddo
!
!  end subroutine det_unq_orb_lab_srt_bld

! ==============================================================================
  subroutine slater_mat_trans_inv_bld
! ------------------------------------------------------------------------------
! Description   : Build inverse of transpose spin-up and down slater matrices
! Description   : in two-dimensional format
!
! Created       : J. Toulouse, 12 Oct 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer det_unq_up_i, det_unq_dn_i, orb_i
  integer elec_up_i, elec_dn_i, mat_i

! header
  if (header_exe) then

   call object_create ('slater_mat_trans_inv_up')
   call object_create ('slater_mat_trans_inv_dn')

   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('slmui')
   call object_needed ('slmdi')

   return

  endif

! begin

! allocations
  call object_alloc ('slater_mat_trans_inv_up', slater_mat_trans_inv_up, nup, nup, ndetup)
  call object_alloc ('slater_mat_trans_inv_dn', slater_mat_trans_inv_dn, ndn, ndn, ndetdn)

! build slater_mat_trans_inv_up
  do det_unq_up_i = 1, ndetup
    mat_i = 0
    do orb_i = 1, nup
      do elec_up_i = 1, nup
        mat_i = mat_i + 1
        slater_mat_trans_inv_up (orb_i, elec_up_i, det_unq_up_i) = slmui (mat_i, det_unq_up_i)
      enddo
    enddo
  enddo

! build slater_mat_trans_inv_dn
  do det_unq_dn_i = 1, ndetdn
    mat_i = 0
    do orb_i = 1, ndn
      do elec_dn_i = 1, ndn
        mat_i = mat_i + 1
        slater_mat_trans_inv_dn (orb_i, elec_dn_i, det_unq_dn_i) = slmdi (mat_i, det_unq_dn_i)
      enddo
    enddo
  enddo

  end subroutine slater_mat_trans_inv_bld

! ==============================================================================
  subroutine orb_occ_in_adet_unq_bld
! ------------------------------------------------------------------------------
! Description   : orbital occupations in unique determinants where the orbital i
! Description   : was been anihilated a_i | det >
!
! Created       : J. Toulouse, 29 Nov 2005
! Modified      : J. Toulouse, 04 Jul 2007: unique detertimants
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer det_unq_up_i,  det_unq_dn_i, orb_i

! header
  if (header_exe) then

   call object_create ('orb_occ_in_adet_unq_up')
   call object_create ('orb_occ_in_adet_unq_dn')
   call object_create ('sgn_adet_unq_up')
   call object_create ('sgn_adet_unq_dn')

   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('orb_tot_nb')
   call object_needed ('orb_occ_in_det_unq_up')
   call object_needed ('orb_occ_in_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('orb_occ_in_adet_unq_up', orb_occ_in_adet_unq_up, orb_tot_nb, orb_tot_nb, ndetup)
  call object_alloc ('orb_occ_in_adet_unq_dn', orb_occ_in_adet_unq_dn, orb_tot_nb, orb_tot_nb, ndetdn)
  call object_alloc ('sgn_adet_unq_up', sgn_adet_unq_up, orb_tot_nb, ndetup)
  call object_alloc ('sgn_adet_unq_dn', sgn_adet_unq_dn, orb_tot_nb, ndetdn)

  orb_occ_in_adet_unq_up = .false.
  orb_occ_in_adet_unq_dn = .false.
  sgn_adet_unq_up = 0
  sgn_adet_unq_dn = 0

  do orb_i = 1, orb_tot_nb

!   loop over unique spin-up determinants
    do det_unq_up_i = 1, ndetup

     orb_occ_in_adet_unq_up (:, orb_i, det_unq_up_i) = orb_occ_in_det_unq_up (:, det_unq_up_i)

!    anihilate orbital orb_i in det up
     if (orb_occ_in_det_unq_up (orb_i, det_unq_up_i)) then
       orb_occ_in_adet_unq_up (orb_i, orb_i, det_unq_up_i) = .false.
       sgn_adet_unq_up (orb_i, det_unq_up_i) = (-1)**(orb_i - 1)
     endif

    enddo ! det_unq_up_i


!   loop over unique spin-dn determinants
    do det_unq_dn_i = 1, ndetdn

     orb_occ_in_adet_unq_dn (:, orb_i, det_unq_dn_i) = orb_occ_in_det_unq_dn (:, det_unq_dn_i)

!    anihilate orbital orb_i in det dn
     if (orb_occ_in_det_unq_dn (orb_i, det_unq_dn_i)) then
       orb_occ_in_adet_unq_dn (orb_i, orb_i, det_unq_dn_i) = .false.
       sgn_adet_unq_dn (orb_i, det_unq_dn_i) = (-1)**(orb_i - 1)
     endif

    enddo ! det_unq_dn_i

   enddo ! orb_i

  end subroutine orb_occ_in_adet_unq_bld

end module determinants_mod
