module deriv_jas_mod

  use all_tools_mod
  use montecarlo_mod
  use eloc_mod

! Declaration of global variables and default values
  integer                                :: jas_pairs_nb
  integer,  allocatable                  :: jas_pairs (:,:)

  real(dp), allocatable          :: dpsi_jas (:)
  real(dp), allocatable          :: dpsi_jas_av (:)
  real(dp), allocatable          :: d2gvalue (:,:)
  real(dp), allocatable          :: d2gvalue_av (:,:)
  real(dp), allocatable          :: d2psi_jas (:)
  real(dp), allocatable          :: d2psi_jas_av (:)
  real(dp), allocatable          :: d2eloc_jas (:)
  real(dp), allocatable          :: d2eloc_jas_av (:)

  real(dp), allocatable          :: dpsi_jas_sq (:)
  real(dp), allocatable          :: dpsi_jas_sq_av (:)
  real(dp), allocatable          :: dpsi_jas_eloc (:)
  real(dp), allocatable          :: dpsi_jas_eloc_av (:)
  real(dp), allocatable          :: dpsi_jas_sq_eloc (:)
  real(dp), allocatable          :: dpsi_jas_sq_eloc_av (:)
  real(dp), allocatable          :: deloc_jas    (:)
  real(dp), allocatable          :: deloc_jas_av (:)
  real(dp), allocatable          :: dpsi_jas_deloc_jas (:)
  real(dp), allocatable          :: dpsi_jas_deloc_jas_av (:)
  real(dp), allocatable          :: e_jas (:)
  real(dp), allocatable          :: delta_e_jas (:)

  contains

! ==============================================================================
  subroutine jas_pairs_bld
! ------------------------------------------------------------------------------
! Description   : pairs of Jastrow paramaters
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j, pair

! header
  if (header_exe) then

   call object_create ('jas_pairs_nb')
   call object_create ('jas_pairs')

   call object_needed ('nparmj')

   return

  endif

! begin
  jas_pairs_nb = nparmj * (nparmj + 1) / 2

  call object_alloc ('jas_pairs', jas_pairs, nparmj, nparmj)
  pair = 0
  do i = 1, nparmj
   do j = i, nparmj
      pair = pair + 1
      jas_pairs (i,j) = pair
      jas_pairs (j,i) = pair
   enddo
  enddo

  end subroutine jas_pairs_bld

! ==============================================================================
  subroutine dpsi_jas_bld
! ------------------------------------------------------------------------------
! Description   :  Logarithmic derivatives of Psi with respect to
! Description   :  jastrow parameters
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

   call object_create ('dpsi_jas')

   call object_needed ('nparmj')
   call object_needed ('gvalue')

   return

  endif

! begin
  call object_alloc ('dpsi_jas', dpsi_jas, nparmj)
  call object_alloc ('dpsi_jas_av', dpsi_jas_av, nparmj)

  do i = 1, nparmj
   dpsi_jas (i) = gvalue (i)
  enddo

  end subroutine dpsi_jas_bld

! ==============================================================================
  subroutine d2psi_jas_bld
! ------------------------------------------------------------------------------
! Description   :  second-derivative of psi wrt to Jastrow parameter divided by psi
! Description   :  (1/psi) * d^2 psi / d c^ 2
!
! Created       : J. Toulouse, 19 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, it, isb, j, iparm, jparm, nparm0, pair

! header
  if (header_exe) then

   call object_create ('d2psi_jas')
   call object_create ('d2gvalue')

   call object_needed ('nparmj')
   call object_needed ('jas_pairs_nb')
   call object_needed ('gvalue')
   call object_needed ('nctype')
   call object_needed ('nparma')
   call object_needed ('d2d2a')
   call object_needed ('d1d2a')
   call object_needed ('iwjasa')
   call object_needed ('nparmb')
   call object_needed ('nspin2b')
   call object_needed ('iwjasb')
   call object_needed ('d2d2b')
   call object_needed ('d1d2b')

   return

  endif

! begin
  call object_alloc ('d2psi_jas', d2psi_jas, jas_pairs_nb)
  call object_alloc ('d2psi_jas_av', d2psi_jas_av, jas_pairs_nb)
  call object_alloc ('d2gvalue', d2gvalue, nparmj, nparmj)
  call object_alloc ('d2gvalue_av', d2gvalue_av, nparmj, nparmj)

  d2gvalue = 0.d0

  nparm0 = 0
  do it = 1, nctype
   nparm0 = nparm0 + npointa(it)
   do i = 1, nparma(it)
    iparm = nparm0 + i
    if (iwjasa(i,it) == 2) then
     d2gvalue(iparm,iparm) = d2d2a(it)
     do j = 1, nparma(it)
      if (iwjasa(j,it) == 1) then
       jparm = nparm0 + j
       if (jparm > iparm) then
        d2gvalue(jparm,iparm) = d1d2a(it)
       else
        d2gvalue(iparm,jparm) = d1d2a(it)
       endif
      endif
     enddo ! j
    endif
   enddo ! i
  enddo ! it

  nparm0=npointa(nctype)+nparma(nctype)
  do isb = 1, nspin2b
   if (isb == 2) nparm0=nparm0+nparmb(1)
    do i=1, nparmb(isb)
     iparm=nparm0+i
      if (iwjasb(i,isb) == 2) then
       d2gvalue(iparm,iparm) = d2d2b(isb)
       do j=1,nparmb(isb)
        if (iwjasb(j,isb) == 1) then
          jparm=nparm0+j
          if (jparm > iparm) then
           d2gvalue(jparm,iparm) = d1d2b(isb)
          else
           d2gvalue(iparm,jparm) = d1d2b(isb)
          endif
        endif
       enddo ! j
      endif
     enddo ! i
    enddo ! isb

! gvalue (i,j) with i<j correct
! but gvalue (i,j) with i>j incorrect!

  pair = 0
  do i = 1, nparmj
   do j = i, nparmj
     pair = pair + 1
     d2psi_jas (pair) = d2gvalue (j, i) + gvalue (i) *  gvalue (j)
   enddo
  enddo

  end subroutine d2psi_jas_bld

! ==============================================================================
  subroutine deloc_jas_bld
! ------------------------------------------------------------------------------
! Description   :  derivative of local energy Psi with respect to
! Description   :  Jastrow  parameters to optimize
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

   call object_create ('deloc_jas')

   call object_needed ('nparmcsf')
   call object_needed ('nparmj')
   call object_needed ('denergy')

   return

  endif

! begin
  call object_alloc ('deloc_jas', deloc_jas, nparmj)
  call object_alloc ('deloc_jas_av', deloc_jas_av, nparmj)

! warning: this is ugly and terribly fragile!
  do i = 1, nparmj
   deloc_jas (i) = denergy (nparmcsf+i)
  enddo

  end subroutine deloc_jas_bld

! ==============================================================================
  subroutine d2eloc_jas_bld
! ------------------------------------------------------------------------------
! Description   :  second-order derivative of local energy Psi with respect to
! Description   :  Jastrow  parameters to optimize
! Description   :  d^2 eloc / d c^2
!
! Description   :  E_{L,ij} = (H Psi_ij)/Psi - (Psi_ij/Psi) E_L - (Psi_i/Psi) E_{L,j} - (Psi_j/Psi) E_{L,i}
! Description   :  or
! Description   :  E_{L,ij} = (-1/2) [ \nabla^2 f_ij + 2 \nabla f_ij . V + 2 \nabla f_i . \nabla f_j] + potential part

! Description   :  f_ij = 0 except for a1,a2,b1,b2
! Description   :  Without pseudopotential and a1,a2,b1,b2, E_{L,ij} =  - \nabla f_i . \nabla f_j
! Description   :                                                    =  -  g(i) . g(j)

! Created       : J. Toulouse, 16 May 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none
!
! local
  integer i,j, elec_i, dim_i
!
! header
  if (header_exe) then

   call object_create ('d2eloc_jas')

   call object_needed ('nparmj')
   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('jas_pairs_nb')
   call object_needed ('jas_pairs')
   call object_needed ('g')

   return

  endif

! begin
  call object_alloc ('d2eloc_jas', d2eloc_jas, jas_pairs_nb)
  call object_alloc ('d2eloc_jas_av', d2eloc_jas_av, jas_pairs_nb)

  d2eloc_jas (:) = 0.d0

  do i = 1, nparmj
   do j = i, nparmj
     do elec_i = 1, nelec
       do dim_i = 1, ndim
         d2eloc_jas (jas_pairs (i,j)) = d2eloc_jas (jas_pairs (i,j)) - g(dim_i,elec_i,i) * g(dim_i,elec_i,j)
       enddo
     enddo
   enddo
  enddo

  end subroutine d2eloc_jas_bld

! ==============================================================================
  subroutine dpsi_jas_sq_bld
! ------------------------------------------------------------------------------
! Description   : square of dpsi_jas
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_jas_sq')

   call object_needed ('nparmj')
   call object_needed ('dpsi_jas')

   return

  endif

! begin
! allocations
  call object_alloc ('dpsi_jas_sq', dpsi_jas_sq, nparmj)
  call object_alloc ('dpsi_jas_sq_av', dpsi_jas_sq_av, nparmj)

  dpsi_jas_sq (:) = dpsi_jas (:)**2

 end subroutine dpsi_jas_sq_bld

! ==============================================================================
  subroutine dpsi_jas_eloc_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_jas_eloc = dpsi_jas * eloc
!
! Created       : J. Toulouse, 18 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_jas_eloc')

   call object_needed ('nparmj')
   call object_needed ('dpsi_jas')
   call object_needed ('eloc')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsi_jas_eloc', dpsi_jas_eloc, nparmj)
  call object_alloc ('dpsi_jas_eloc_av', dpsi_jas_eloc_av, nparmj)

  dpsi_jas_eloc = dpsi_jas * eloc

 end subroutine dpsi_jas_eloc_bld

! ==============================================================================
  subroutine dpsi_jas_sq_eloc_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_jas_sq_eloc = dpsi_jas_sq * eloc
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_jas_sq_eloc')

   call object_needed ('nparmj')
   call object_needed ('dpsi_jas_sq')
   call object_needed ('eloc')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsi_jas_sq_eloc', dpsi_jas_sq_eloc, nparmj)
  call object_alloc ('dpsi_jas_sq_eloc_av', dpsi_jas_sq_eloc_av, nparmj)

  dpsi_jas_sq_eloc = dpsi_jas_sq * eloc

 end subroutine dpsi_jas_sq_eloc_bld

! ==============================================================================
  subroutine dpsi_jas_deloc_jas_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_jas * deloc_jas
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local

! header
  if (header_exe) then

   call object_create ('dpsi_jas_deloc_jas')

   call object_needed ('nparmj')
   call object_needed ('dpsi_jas')
   call object_needed ('deloc_jas')

   return

  endif

! begin
! allocations
  call object_alloc ('dpsi_jas_deloc_jas', dpsi_jas_deloc_jas, nparmj)
  call object_alloc ('dpsi_jas_deloc_jas_av', dpsi_jas_deloc_jas_av, nparmj)

  dpsi_jas_deloc_jas (:) = dpsi_jas (:) * deloc_jas (:)

  end subroutine dpsi_jas_deloc_jas_bld

! ==============================================================================
  subroutine e_jas_bld
! ------------------------------------------------------------------------------
! Description   :  energy of derivatives of the Jastrow
! Description   :  (diagonal of Hamiltonian in linear method)
!
! Created       : J. Toulouse, 17 Feb 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('e_jas')

   call object_needed ('nparmj')
   call object_needed ('dpsi_jas_sq_av')
   call object_needed ('dpsi_jas_sq_eloc_av')
   call object_needed ('dpsi_jas_deloc_jas_av')

   call object_needed ('dpsi_jas_av')
   call object_needed ('deloc_jas_av')

   return

  endif

! begin

! allocations
  call object_alloc ('e_jas', e_jas, nparmj)

! energy with semiorthogonal basis
  e_jas (:) = ( dpsi_jas_sq_eloc_av (:) + dpsi_jas_deloc_jas_av (:)               &
                - dpsi_jas_av (:) * dpsi_jas_eloc_av (:)                          &
                - dpsi_jas_av (:) * ( dpsi_jas_eloc_av (:) + deloc_jas_av (:) )   &
                + dpsi_jas_av (:) * dpsi_jas_av (:) * eloc_av )                   &
               / (dpsi_jas_sq_av (:) - dpsi_jas_av (:) **2 )

 end subroutine e_jas_bld

! ==============================================================================
  subroutine delta_e_jas_bld
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

   call object_create ('delta_e_jas')

   call object_needed ('nparmj')
   call object_needed ('e_jas')
   call object_needed ('eloc_av')

   return

  endif

! begin

! allocations
  call object_alloc ('delta_e_jas', delta_e_jas, nparmj)

  delta_e_jas (:)  = e_jas (:) - eloc_av

  call object_write (here, 'delta_e_jas')

 end subroutine delta_e_jas_bld

end module deriv_jas_mod

