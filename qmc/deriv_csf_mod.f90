module deriv_csf_mod

  use all_tools_mod
  use montecarlo_mod
  use eloc_mod
  use determinants_mod

! Declaration of global variables and default values
  integer                                :: csf_pairs_nb
  integer,  allocatable                  :: csf_pairs (:,:)

!  real(dp)                       :: det1_det
  real(dp), allocatable          :: dpsi_csf (:)
  real(dp), allocatable          :: dpsi_csf_av (:)

  real(dp), allocatable          :: dpsi_csf_dpsi_csf (:)
  real(dp), allocatable          :: dpsi_csf_dpsi_csf_av (:)
  real(dp), allocatable          :: dpsi_csf_dpsi_csf_covar (:,:)
  real(dp), allocatable          :: dpsi_csf_sq (:)
  real(dp), allocatable          :: dpsi_csf_sq_av (:)
  real(dp), allocatable          :: dpsi_csf_eloc (:)
  real(dp), allocatable          :: dpsi_csf_eloc_av (:)
  real(dp), allocatable          :: dpsi_csf_sq_eloc (:)
  real(dp), allocatable          :: dpsi_csf_sq_eloc_av (:)
  real(dp), allocatable          :: deloc_csf    (:)
  real(dp), allocatable          :: deloc_csf_av (:)
  real(dp), allocatable          :: dpsi_csf_deloc_csf (:)
  real(dp), allocatable          :: dpsi_csf_deloc_csf_av (:)
  real(dp), allocatable          :: e_csf (:)
  real(dp), allocatable          :: delta_e_csf (:)

  contains

! ==============================================================================
  subroutine csf_pairs_bld
! ------------------------------------------------------------------------------
! Description   : pairs of csfs parameters
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j, pair

! header
  if (header_exe) then

   call object_create ('csf_pairs_nb')
   call object_create ('csf_pairs')

   call object_needed ('nparmcsf')

   return

  endif

! begin
  csf_pairs_nb = nparmcsf * (nparmcsf + 1) / 2

  call object_alloc ('csf_pairs', csf_pairs, nparmcsf, nparmcsf)
  pair = 0
  do i = 1, nparmcsf
   do j = i, nparmcsf
      pair = pair + 1
      csf_pairs (i,j) = pair
      csf_pairs (j,i) = pair
   enddo
  enddo

  end subroutine csf_pairs_bld

! ==============================================================================
  subroutine dpsi_csf_bld
! ------------------------------------------------------------------------------
! Description   :  Logarithmic derivatives of Psi with respect to
! Description   :  csf parameters
! Description   :  d ln Psi / d j = (1/Psi) * d Psi / d j
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i

! header
  if (header_exe) then

   call object_create ('dpsi_csf')

   call object_needed ('nparmcsf')
   call object_needed ('deti_det')
!temp   call object_needed ('csf_coef') !!   commented because causes circular dependencies in ovlp_ovlp_fn
!   call object_needed ('iwcsf')  !!
!   call object_needed ('det1_det')  !!

   return

  endif

! begin
  call object_alloc ('dpsi_csf', dpsi_csf, nparmcsf)
  call object_alloc ('dpsi_csf_av', dpsi_csf_av, nparmcsf)

  do i = 1, nparmcsf
    dpsi_csf (i) = deti_det (i)

! rotation
!    dpsi_csf (i) = deti_det (i+1) -  csf_coef(i+1,1)/(1.d0 + csf_coef(1,1)) * (1.d0 + deti_det(1))
  enddo

  end subroutine dpsi_csf_bld

! ==============================================================================
  subroutine deloc_csf_bld
! ------------------------------------------------------------------------------
! Description   :  derivative of local energy Psi with respect to
! Description   :  csf  parameters to optimize
! Description   :  d eloc / d c
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i

! header
  if (header_exe) then

   call object_create ('deloc_csf')

   call object_needed ('nparmcsf')
   call object_needed ('denergy')
   call object_needed ('csf_coef') !!

   return

  endif

! begin
  call object_alloc ('deloc_csf', deloc_csf, nparmcsf)
  call object_alloc ('deloc_csf_av', deloc_csf_av, nparmcsf)

! warning: this is terribly fragile!
  do i = 1, nparmcsf
   deloc_csf (i) = denergy (i)

! rotation
!   deloc_csf (i) = denergy (i+1) - csf_coef(i+1,1)/(1.d0 + csf_coef(1,1)) * denergy (1)
  enddo

  end subroutine deloc_csf_bld

! ==============================================================================
  subroutine dpsi_csf_sq_bld
! ------------------------------------------------------------------------------
! Description   : square of dpsi_csf
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_csf_sq')

   call object_needed ('nparmcsf')
   call object_needed ('dpsi_csf')

   return

  endif

! begin
! allocations
  call object_alloc ('dpsi_csf_sq', dpsi_csf_sq, nparmcsf)
  call object_alloc ('dpsi_csf_sq_av', dpsi_csf_sq_av, nparmcsf)

  dpsi_csf_sq (:) = dpsi_csf (:)**2

 end subroutine dpsi_csf_sq_bld

! ==============================================================================
  subroutine dpsi_csf_dpsi_csf_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_csf_dpsi_csf = dpsi_csf (i) * dpsi_csf (j)
!
! Created       : J. Toulouse, 27 Mar 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j

! header
  if (header_exe) then

   call object_create ('dpsi_csf_dpsi_csf')

   call object_needed ('dpsi_csf')
   call object_needed ('nparmcsf')
   call object_needed ('csf_pairs_nb')
   call object_needed ('csf_pairs')

   return

  endif

! begin
  call object_alloc ('dpsi_csf_dpsi_csf', dpsi_csf_dpsi_csf, csf_pairs_nb)
  call object_alloc ('dpsi_csf_dpsi_csf_av', dpsi_csf_dpsi_csf_av, csf_pairs_nb)

  do i = 1, nparmcsf
   do j = i, nparmcsf
    dpsi_csf_dpsi_csf (csf_pairs(i,j))= dpsi_csf (i) *  dpsi_csf (j)
   enddo
  enddo

  end subroutine dpsi_csf_dpsi_csf_bld

! ==============================================================================
  subroutine dpsi_csf_dpsi_csf_covar_bld
! ------------------------------------------------------------------------------
! Description   :  < dpsi_csf (i) * dpsi_csf (j) > - <dpsi_csf (i)> <dpsi_csf (j)>
!
! Created       : J. Toulouse, 27 Mar 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j

! header
  if (header_exe) then

   call object_create ('dpsi_csf_dpsi_csf_covar')

   call object_needed ('nparmcsf')
   call object_needed ('csf_pairs')
   call object_needed ('dpsi_csf_av')
   call object_needed ('dpsi_csf_dpsi_csf_av')

   return

  endif

! begin
  call object_alloc ('dpsi_csf_dpsi_csf_covar', dpsi_csf_dpsi_csf_covar, nparmcsf, nparmcsf)

  do i = 1, nparmcsf
   do j = 1, nparmcsf
    dpsi_csf_dpsi_csf_covar (i,j)= dpsi_csf_dpsi_csf_av (csf_pairs(i,j)) - dpsi_csf_av (i) * dpsi_csf_av (j)
   enddo
  enddo

  end subroutine dpsi_csf_dpsi_csf_covar_bld

! ==============================================================================
  subroutine dpsi_csf_eloc_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_csf_eloc = dpsi_csf * eloc
!
! Created       : J. Toulouse, 18 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_csf_eloc')

   call object_needed ('nparmj')
   call object_needed ('dpsi_csf')
   call object_needed ('eloc')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsi_csf_eloc', dpsi_csf_eloc, nparmj)
  call object_alloc ('dpsi_csf_eloc_av', dpsi_csf_eloc_av, nparmj)

  dpsi_csf_eloc = dpsi_csf * eloc

 end subroutine dpsi_csf_eloc_bld

! ==============================================================================
  subroutine dpsi_csf_sq_eloc_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_csf_sq_eloc = dpsi_csf_sq * eloc
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_csf_sq_eloc')

   call object_needed ('nparmcsf')
   call object_needed ('dpsi_csf_sq')
   call object_needed ('eloc')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsi_csf_sq_eloc', dpsi_csf_sq_eloc, nparmcsf)
  call object_alloc ('dpsi_csf_sq_eloc_av', dpsi_csf_sq_eloc_av, nparmcsf)

  dpsi_csf_sq_eloc = dpsi_csf_sq * eloc

 end subroutine dpsi_csf_sq_eloc_bld

! ==============================================================================
  subroutine dpsi_csf_deloc_csf_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_csf * deloc_csf
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_csf_deloc_csf')

   call object_needed ('nparmcsf')
   call object_needed ('dpsi_csf')
   call object_needed ('deloc_csf')

   return

  endif

! begin
! allocations
  call object_alloc ('dpsi_csf_deloc_csf', dpsi_csf_deloc_csf, nparmcsf)
  call object_alloc ('dpsi_csf_deloc_csf_av', dpsi_csf_deloc_csf_av, nparmcsf)

  dpsi_csf_deloc_csf (:) = dpsi_csf (:) * deloc_csf (:)

  end subroutine dpsi_csf_deloc_csf_bld

! ==============================================================================
  subroutine e_csf_bld
! ------------------------------------------------------------------------------
! Description   :  energy of derivatives of the Jastrow  (diagonal of Hamiltonian in linear method)
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('e_csf')

   call object_needed ('nparmcsf')
   call object_needed ('dpsi_csf_sq_av')
   call object_needed ('dpsi_csf_sq_eloc_av')
   call object_needed ('dpsi_csf_deloc_csf_av')

   return

  endif

! begin

! allocations
  call object_alloc ('e_csf', e_csf, nparmcsf)

  e_csf (:) = (dpsi_csf_sq_eloc_av (:) + dpsi_csf_deloc_csf_av (:) ) / dpsi_csf_sq_av (:)

 end subroutine e_csf_bld

! ==============================================================================
  subroutine delta_e_csf_bld
! ------------------------------------------------------------------------------
! Description   : denominator energy differences for optimization of Jastrow
! Description   : with perturbation method
!
! Created       : J. Toulouse, 17 Feb 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('delta_e_csf')

   call object_needed ('nparmcsf')
   call object_needed ('e_csf')
   call object_needed ('eloc_av')

   return

  endif

! begin

! allocations
  call object_alloc ('delta_e_csf', delta_e_csf, nparmcsf)

  delta_e_csf (:)  = e_csf (:) - eloc_av

 end subroutine delta_e_csf_bld

end module deriv_csf_mod

