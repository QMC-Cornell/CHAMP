module derivatives_fast_mod

  use all_tools_mod
  use electrons_mod
  use psi_mod

! Declaration of global variables and default values
  real(dp),allocatable           :: lap_det_over_det_fast(:)
  real(dp),allocatable           :: lap_det_over_det_fast_av(:)
  real(dp),allocatable           :: lap_det_over_det_fast_av_err(:)
  real(dp)                       :: lap_det_over_det_legacy
  real(dp),allocatable           :: grd_det_over_det_fast(:,:,:)
  real(dp),allocatable           :: grd_det_over_det_legacy(:,:)

  contains

! ==============================================================================
  subroutine grd_det_over_det_fast_bld
! ------------------------------------------------------------------------------
! Description   : Grad(Psi)/Psi
!
! Created       : B. Mussard, Jan 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer dim_i,elec_i,elec_k

! header
  if (header_exe) then

   call object_create ('grd_det_over_det_fast')

   call object_needed('ndim')
   call object_needed('nelec')
   !call object_needed('nwalk')
   call object_needed('dorb')
   call object_needed('det_unq_orb_lab_srt_up')
   call object_needed('det_unq_orb_lab_srt_dn')
   call object_needed('slater_mat_trans_inv_up')
   call object_needed('slater_mat_trans_inv_dn')
   call object_needed('grd_det_over_det_legacy')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_det_over_det_fast', grd_det_over_det_fast, ndim, nelec, 1)

  ! for all dimensions and electrons:
  ! Grad(DET)/DET = -1/2 tr(A^-1.B)
  ! where A is Slater determinant
  ! and   B_ij=Grad_(dim,i)(\phi_j(r_i))
  grd_det_over_det_fast = 0.d0
  do dim_i=1,ndim
    do elec_i=1,nup
      do elec_k=1,nup
        grd_det_over_det_fast(dim_i,elec_i,1)=grd_det_over_det_fast(dim_i,elec_i,1) &
                                           & +slater_mat_trans_inv_up(elec_i,elec_k,1)*dorb(dim_i,elec_i,det_unq_orb_lab_srt_up(elec_k,1))
      enddo
    enddo
    do elec_i=1,ndn
      do elec_k=1,ndn
        grd_det_over_det_fast(dim_i,nup+elec_i,1)=grd_det_over_det_fast(dim_i,nup+elec_i,1) &
                                               & +slater_mat_trans_inv_dn(elec_i,elec_k,1)*dorb(dim_i,nup+elec_i,det_unq_orb_lab_srt_dn(elec_k,1))
      enddo
    enddo
  enddo

  call is_equal_or_die(grd_det_over_det_legacy, grd_det_over_det_fast(:,:,1), 10.d-10, .false.)

  end subroutine grd_det_over_det_fast_bld

! ==============================================================================
  subroutine lap_det_over_det_fast_bld
! ------------------------------------------------------------------------------
! Description   : Kinetic local energy
!
! Created       : B. Mussard, Jan 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer elec_i,orb_j

! header
  if (header_exe) then

   call object_create ('lap_det_over_det_fast')

   call object_needed('nup')
   call object_needed('ndn')
   call object_needed('nwalk')
   call object_needed('ddorb')
   call object_needed('det_unq_orb_lab_srt_up')
   call object_needed('det_unq_orb_lab_srt_dn')
   call object_needed('slater_mat_trans_inv_up')
   call object_needed('slater_mat_trans_inv_dn')
   call object_needed('lap_det_over_det_legacy')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_det_over_det_fast', lap_det_over_det_fast, 1)
  call object_alloc ('lap_det_over_det_fast_av', lap_det_over_det_fast_av, 1)
  call object_alloc ('lap_det_over_det_fast_av_err', lap_det_over_det_fast_av_err, 1)
  call object_associate ('lap_det_over_det_fast', lap_det_over_det_fast,1)
  call object_associate ('lap_det_over_det_fast_av', lap_det_over_det_fast_av,1)
  call object_associate ('lap_det_over_det_fast_av_err', lap_det_over_det_fast_av_err,1)

  ! T(DET)/DET = -1/2 tr(A^-1.B)
  ! where A is Slater determinant
  ! and   B_ij=\laplacian\phi_j(r_i)
  lap_det_over_det_fast = 0.d0
  do elec_i=1,nup
    do orb_j=1,nup
      lap_det_over_det_fast(1)=lap_det_over_det_fast(1) &
                            & +slater_mat_trans_inv_up(elec_i,orb_j,1)*ddorb(elec_i,det_unq_orb_lab_srt_up(orb_j,1))
    enddo
  enddo
  do elec_i=1,ndn
    do orb_j=1,ndn
      lap_det_over_det_fast(1)=lap_det_over_det_fast(1) &
                            & +slater_mat_trans_inv_dn(elec_i,orb_j,1)*ddorb(nup+elec_i,det_unq_orb_lab_srt_dn(orb_j,1))
    enddo
  enddo
  lap_det_over_det_fast(1)=-0.5d0*lap_det_over_det_fast(1)

  if (abs(lap_det_over_det_fast(1)-lap_det_over_det_legacy).ge.1.d-10) then
    write(6,'(a,f15.10)') 'is_equal_or_die_scalar:   checking equality ',lap_det_over_det_fast(1)-lap_det_over_det_legacy
  endif

  end subroutine lap_det_over_det_fast_bld

end module derivatives_fast_mod
