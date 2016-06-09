module forces_pulay_mod

  use all_tools_mod
  use basis_mod
  use orbitals_mod
  use determinants_mod
  use jastrow_mod
  use montecarlo_mod
  use eloc_mod

! Declaration of global variables and default values
  integer, allocatable                   :: force_to_bas_nb (:)
  type (type_integer_row), allocatable   :: force_to_bas (:)
  real(dp), allocatable                  :: dphin_norm_drn (:,:,:)
  real(dp), allocatable                  :: dorb_drn (:,:,:)
  real(dp), allocatable                  :: ddet_drn_unq_up (:,:)
  real(dp), allocatable                  :: ddet_drn_unq_dn (:,:)
  real(dp), allocatable                  :: ddet_drn_col_unq_up (:,:,:)
  real(dp), allocatable                  :: ddet_drn_col_unq_dn (:,:,:)
  real(dp), allocatable                  :: dpsid_drn (:)
  real(dp), allocatable                  :: dpsid_rn (:)
  real(dp), allocatable                  :: dpsij_drn (:)
  real(dp), allocatable                  :: dpsij_rn (:)
  real(dp), allocatable                  :: dpsi_rn (:)
  real(dp), allocatable                  :: dpsi_rn_num (:)
  real(dp), allocatable                  :: dpsi_rn_av (:)
  real(dp), allocatable                  :: dpsi_rn_av_var (:)
  real(dp), allocatable                  :: dpsi_rn_eloc (:)
  real(dp), allocatable                  :: dpsi_rn_eloc_av (:)
  real(dp), allocatable                  :: dpsi_rn_eloc_av_var (:)
  real(dp), allocatable                  :: dpsi_rn_av_eloc_av_covar (:)
  real(dp), allocatable                  :: dpsi_rn_eloc_av_eloc_av_covar (:)
  real(dp), allocatable                  :: dpsi_rn_eloc_av_dpsi_rn_av_covar (:)
  real(dp), allocatable                  :: deloc_rn_num (:)
  real(dp), allocatable                  :: deloc_rn_num_av (:)
  real(dp), allocatable                  :: forces_pulay_av (:)
  real(dp), allocatable                  :: forces_pulay_av_var (:)
  real(dp), allocatable                  :: forces_pulay_av_err (:)

  contains

! ==============================================================================
  subroutine force_to_bas_bld
! ------------------------------------------------------------------------------
! Description   : correspondance array from a force component to the invoved basis functions
!
! Created       : J. Toulouse, 27 Jul 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, cent_i, bas_i

! header
  if (header_exe) then

   call object_create ('force_to_bas_nb')
   call object_create ('force_to_bas')

   call object_needed ('forces_nb')
   call object_needed ('forces_cent')
   call object_needed ('nbasis')
   call object_needed ('basis_fns_cent')

   return

  endif

! begin

! allocation
  call object_alloc ('force_to_bas_nb', force_to_bas_nb, forces_nb)
  call object_alloc ('force_to_bas', force_to_bas, forces_nb)

  do force_i = 1, forces_nb
    cent_i = forces_cent (force_i)
    do bas_i = 1, nbasis
      if (cent_i == basis_fns_cent (bas_i)) then
         force_to_bas_nb (force_i) = force_to_bas_nb (force_i) + 1
         call append (force_to_bas (force_i)%row, bas_i)
      endif
    enddo ! bas_i
  enddo ! force_i

  end subroutine force_to_bas_bld

! ==============================================================================
  subroutine dphin_norm_drn_bld
! ------------------------------------------------------------------------------
! Description   : derivatives of the normalized basis functions wrt to nuclear coordinate
! Description   : this is simply the opposite of the derivatives wrt electron coordinate
!
! Created       : J. Toulouse, 27 Jul 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer bas_i, elec_i, dim_i

! header
  if (header_exe) then

   call object_create ('dphin_norm_drn')

   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('dphin')
   call object_needed ('norm_basis')

   return

  endif

! begin

! allocation
  call object_alloc ('dphin_norm_drn', dphin_norm_drn, ndim, nelec, nbasis)

  do bas_i = 1, nbasis
    do elec_i = 1, nelec
      do dim_i = 1, ndim
        dphin_norm_drn (dim_i, elec_i, bas_i) = - norm_basis (bas_i) * dphin (dim_i, bas_i, elec_i)
      enddo ! elec_i
    enddo ! dim_i
  enddo ! bas_i

  end subroutine dphin_norm_drn_bld

! ===============================================================================
  subroutine dorb_drn_bld
! -------------------------------------------------------------------------------
! Description   : derivatives of occupied orbitals wrt to nuclear coordinates
!
! Created       : J. Toulouse, 27 Jul 2008
! -------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, dim_i, bas_i, orb_i, force_to_bas_i

! header
  if (header_exe) then

   call object_create ('dorb_drn')

   call object_needed ('forces_nb')
   call object_needed ('orb_occ_last_in_wf_lab')
   call object_needed ('nelec')
   call object_needed ('force_to_bas_nb')
   call object_needed ('force_to_bas')
!   call object_needed ('coef_orb_on_norm_basis')
   call object_needed ('coef')
!   call object_needed ('dphin_norm_drn')
   call object_needed ('dphin')

   return

  endif

! begin

! allocation
  call object_alloc ('dorb_drn', dorb_drn, nelec, orb_occ_last_in_wf_lab, forces_nb)
  dorb_drn (:,:,:) = 0.d0

  do force_i = 1, forces_nb
   dim_i = forces_direct (force_i)
   do orb_i = 1, orb_occ_last_in_wf_lab
    do force_to_bas_i = 1, force_to_bas_nb (force_i)
     bas_i = force_to_bas (force_i)%row (force_to_bas_i)
!     dorb_drn (1:nelec, orb_i, force_i) = dorb_drn (1:nelec, orb_i, force_i) + coef_orb_on_norm_basis (bas_i, orb_i, 1) * dphin_norm_drn (dim_i, 1:nelec, bas_i)
     dorb_drn (1:nelec, orb_i, force_i) = dorb_drn (1:nelec, orb_i, force_i) - coef (bas_i, orb_i, 1) * dphin (dim_i, bas_i, 1:nelec)      
    enddo ! force_to_bas_i
   enddo ! orb_i
  enddo ! force_i

  end subroutine dorb_drn_bld

! ==============================================================================
  subroutine ddet_drn_unq_bld
! ------------------------------------------------------------------------------
! Description   : derivatives of unique spin-up and spin-down determinants
! Description   : with respect to nuclear coordinates
! Description   : by updating reference determinant via the Sherman-Morrison formula
!
! Created       : J. Toulouse, 27 Jul 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer det_unq_up_i, det_unq_dn_i
  integer orb_i, col_i, i
  integer force_i
  real(dp) factor_up, factor_dn

! header
  if (header_exe) then

   call object_create ('ddet_drn_unq_up')
   call object_create ('ddet_drn_unq_dn')
   call object_create ('ddet_drn_col_unq_up')
   call object_create ('ddet_drn_col_unq_dn')

   call object_needed ('forces_nb')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('dorb_drn')
   call object_needed ('detu')
   call object_needed ('detd')

   return

  endif

! begin

! allocations
  call object_alloc ('ddet_drn_unq_up', ddet_drn_unq_up, ndetup, forces_nb)
  call object_alloc ('ddet_drn_unq_dn', ddet_drn_unq_dn, ndetdn, forces_nb)
  call object_alloc ('ddet_drn_col_unq_up', ddet_drn_col_unq_up, nup, ndetup, forces_nb)
  call object_alloc ('ddet_drn_col_unq_dn', ddet_drn_col_unq_dn, ndn, ndetdn, forces_nb)

! loop over force components
  do force_i = 1, forces_nb

!  loop over unique spin-up determinants
   do det_unq_up_i = 1, ndetup

     ddet_drn_unq_up (det_unq_up_i, force_i) = 0.d0

!    loop over columns (=orbitals) of determinant
     do col_i = 1, nup

       orb_i = det_unq_orb_lab_srt_up (col_i, det_unq_up_i)

       factor_up = 0.d0
       do i = 1, nup
        factor_up = factor_up + slater_mat_trans_inv_up (i, col_i, det_unq_up_i) * dorb_drn (i, orb_i, force_i)
       enddo

       ddet_drn_col_unq_up (col_i, det_unq_up_i, force_i) = factor_up * detu (det_unq_up_i)
       ddet_drn_unq_up (det_unq_up_i, force_i) = ddet_drn_unq_up (det_unq_up_i, force_i) + ddet_drn_col_unq_up (col_i, det_unq_up_i, force_i)

     enddo ! col_i

   enddo ! det_unq_up_i

!  loop over unique spin-dn determinants
   do det_unq_dn_i = 1, ndetdn

     ddet_drn_unq_dn (det_unq_dn_i, force_i) = 0.d0

!    loop over columns (=orbitals) of determinant
     do col_i = 1, ndn

       orb_i = det_unq_orb_lab_srt_dn (col_i, det_unq_dn_i)

       factor_dn = 0.d0
       do i = 1, ndn
        factor_dn = factor_dn + slater_mat_trans_inv_dn (i, col_i, det_unq_dn_i) * dorb_drn (nup + i, orb_i, force_i)
       enddo

       ddet_drn_col_unq_dn (col_i, det_unq_dn_i, force_i) = factor_dn * detd (det_unq_dn_i)
       ddet_drn_unq_dn (det_unq_dn_i, force_i) = ddet_drn_unq_dn (det_unq_dn_i, force_i) + ddet_drn_col_unq_dn (col_i, det_unq_dn_i, force_i)

     enddo ! col_i

   enddo ! det_unq_dn_i

  enddo ! force_i

  end subroutine ddet_drn_unq_bld

! ==============================================================================
  subroutine dpsid_rn_bld
! ------------------------------------------------------------------------------
! Description   : logarithmic derivatives of determinantal part of Psi with respect to nuclear coordinates
!
! Created       : J. Toulouse, 27 Jul 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i
  integer csf_i, det_in_csf_i, det_i
  integer det_unq_up_i, det_unq_dn_i

! header
  if (header_exe) then

   call object_create ('dpsid_rn')
   call object_create ('dpsid_drn')

   call object_needed ('forces_nb')
   call object_needed ('ncsf')
   call object_needed ('ndet')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('csf_coef')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('ddet_drn_unq_up')
   call object_needed ('ddet_drn_unq_dn')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('psi_det')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsid_rn', dpsid_rn, forces_nb)
  call object_alloc ('dpsid_drn', dpsid_drn, forces_nb)

  dpsid_drn (:) = 0.d0

! loop over force components
  do force_i = 1, forces_nb

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

        det_i = iwdet_in_csf (det_in_csf_i, csf_i)
        det_unq_up_i = det_to_det_unq_up (det_i)
        det_unq_dn_i = det_to_det_unq_dn (det_i)

        dpsid_drn (force_i) = dpsid_drn (force_i) + csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i) *  &
        (ddet_drn_unq_up (det_unq_up_i, force_i) * detd (det_unq_dn_i) + detu (det_unq_up_i) * ddet_drn_unq_dn (det_unq_dn_i, force_i))

     enddo ! det_in_csf_i
   enddo ! csf_i

    dpsid_rn (force_i) = dpsid_drn (force_i) / psi_det

  enddo ! force_i

  end subroutine dpsid_rn_bld

! ==============================================================================
  subroutine dpsij_rn_bld
! ------------------------------------------------------------------------------
! Description   : logarithmic derivatives of Jastrow part of Psi with respect to nuclear coordinates
!
! Created       : J. Toulouse, 28 Jul 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsij_rn')

   call object_needed ('forces_nb')
   call object_needed ('dfen_drn')
   call object_needed ('dfeen_drn')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsij_rn', dpsij_rn, forces_nb)

  dpsij_rn (:) = dfen_drn (:) + dfeen_drn (:)

  end subroutine dpsij_rn_bld

! ==============================================================================
  subroutine dpsi_rn_bld
! ------------------------------------------------------------------------------
! Description   : logarithmic derivatives of Psi with respect to nuclear coordinates
!
! Created       : J. Toulouse, 27 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! header
  if (header_exe) then

   call object_create ('dpsi_rn')
   call object_average_define ('dpsi_rn', 'dpsi_rn_av')
   call object_variance_define ('dpsi_rn_av', 'dpsi_rn_av_var')
   call object_covariance_define ('dpsi_rn_av', 'eloc_av', 'dpsi_rn_av_eloc_av_covar')

   call object_needed ('forces_nb')
   call object_needed ('dpsij_rn')
   call object_needed ('dpsid_rn')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsi_rn', dpsi_rn, forces_nb)
  call object_alloc ('dpsi_rn_av', dpsi_rn_av, forces_nb)
  call object_alloc ('dpsi_rn_av_var', dpsi_rn_av_var, forces_nb)
  call object_alloc ('dpsi_rn_av_eloc_av_covar', dpsi_rn_av_eloc_av_covar, forces_nb)

  dpsi_rn (:) = dpsij_rn (:) + dpsid_rn (:)

  end subroutine dpsi_rn_bld

! ==============================================================================
  subroutine dpsi_rn_num_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_rn calculating by numerical differentiation for checkings
!
! Created       : J. Toulouse, 28 Jul 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, cent_i, dim_i
  real(dp) :: cent_save, d, psi1, psi2, dpsi

! header
  if (header_exe) then

   call object_create ('dpsi_rn_num')

   call object_needed ('forces_nb')
   call object_needed ('forces_cent')
   call object_needed ('forces_direct')
   call object_needed ('cent')
   call object_needed ('psijo')
   call object_needed ('psid_det')

   return

  endif

! begin
  call object_alloc ('dpsi_rn_num', dpsi_rn_num, forces_nb)

! calculating wave function
!  call hpsi(xold,psi_det,psijo,vold,div_vo,d2o,peo,peio,eold(1),denergy,1)
!  call object_modified_by_index (eold_index)
!  call object_modified_by_index (eloc_index)
!  call object_modified_by_index (psi_det_index)
!  call object_modified_by_index (denergy_index)
!  call object_modified_by_index (vold_index)
!  call object_modified_by_index (div_vo_index)
  psi1 = dexp(psijo)*psi_det

! loop over force components
  do force_i = 1, forces_nb
   cent_i = forces_cent (force_i)
   dim_i = forces_direct (force_i)

!  moving nuclear coordinates
   d = 0.000001d0
   cent_save = cent (dim_i, cent_i)
   cent (dim_i, cent_i) = cent (dim_i, cent_i) + d
   call object_modified ('cent')

!  recalculating wave function
   call hpsi(xold,psi_det,psijo,vold,div_vo,d2o,peo,peio,eold(1),denergy,1)
   call object_modified_by_index (eold_index)
   call object_modified_by_index (eloc_index)
   call object_modified_by_index (psi_det_index)
   call object_modified_by_index (denergy_index)
   call object_modified_by_index (vold_index)
   call object_modified_by_index (div_vo_index)
   psi2 = dexp(psijo)*psi_det

!  calculating numerical derivative
   dpsi = (psi2 - psi1)/d
   dpsi_rn_num (force_i) = dpsi /psi1
   
!  restoring the nuclear coordinates
   cent (dim_i, cent_i) = cent_save
   call object_modified ('cent')

  enddo ! force_i
  
  end subroutine dpsi_rn_num_bld

! ==============================================================================
  subroutine dpsi_rn_eloc_bld
! ------------------------------------------------------------------------------
! Description   : dpsi_rn * eloc 
!
! Created       : J. Toulouse, 27 Jul 2008
! ------------------------------------------------------------------------------
  implicit none


! header
  if (header_exe) then

   call object_create ('dpsi_rn_eloc')
   call object_average_define ('dpsi_rn_eloc', 'dpsi_rn_eloc_av')
   call object_variance_define ('dpsi_rn_eloc_av', 'dpsi_rn_eloc_av_var')
   call object_covariance_define ('dpsi_rn_eloc_av', 'eloc_av', 'dpsi_rn_eloc_av_eloc_av_covar')
   call object_covariance_define ('dpsi_rn_eloc_av', 'dpsi_rn_av', 'dpsi_rn_eloc_av_dpsi_rn_av_covar')

   call object_needed ('forces_nb')
   call object_needed ('dpsi_rn')
   call object_needed ('eloc')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsi_rn_eloc', dpsi_rn_eloc, forces_nb)
  call object_alloc ('dpsi_rn_eloc_av', dpsi_rn_eloc_av, forces_nb)
  call object_alloc ('dpsi_rn_eloc_av_var', dpsi_rn_eloc_av_var, forces_nb)
  call object_alloc ('dpsi_rn_eloc_av_eloc_av_covar', dpsi_rn_eloc_av_eloc_av_covar, forces_nb)
  call object_alloc ('dpsi_rn_eloc_av_dpsi_rn_av_covar', dpsi_rn_eloc_av_dpsi_rn_av_covar, forces_nb)

  dpsi_rn_eloc (:) = dpsi_rn (:) * eloc

  end subroutine dpsi_rn_eloc_bld

! ==============================================================================
  subroutine deloc_rn_num_bld
! ------------------------------------------------------------------------------
! Description   : deloc_rn calculating by numerical differentiation for checkings
!
! Created       : J. Toulouse, 04 Aug 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer force_i, cent_i, dim_i
  real(dp) :: cent_save, d, eloc1, eloc2

! header
  if (header_exe) then

   call object_create ('deloc_rn_num')
   call object_average_define ('deloc_rn_num', 'deloc_rn_num_av')

   call object_needed ('forces_nb')
   call object_needed ('forces_cent')
   call object_needed ('forces_direct')
   call object_needed ('cent')
   call object_needed ('eloc')

   return

  endif

! begin
  call object_alloc ('deloc_rn_num', deloc_rn_num, forces_nb)
  call object_alloc ('deloc_rn_num_av', deloc_rn_num_av, forces_nb)

! current local energy
  eloc1 = eloc

! loop over force components
  do force_i = 1, forces_nb
   cent_i = forces_cent (force_i)
   dim_i = forces_direct (force_i)

!  moving nuclear coordinates
   d = 0.000001d0
   cent_save = cent (dim_i, cent_i)
   cent (dim_i, cent_i) = cent (dim_i, cent_i) + d
   call object_modified ('cent')

!  recalculating wave function and local energy (including nuclear potential energy!)
   call pot_nn(cent,znuc,iwctype,ncent,pecent)
   call hpsi(xold,psi_det,psijo,vold,div_vo,d2o,peo,peio,eold(1),denergy,1)
   call object_modified_by_index (eold_index)
   call object_modified_by_index (eloc_index)
   call object_modified_by_index (psi_det_index)
   call object_modified_by_index (denergy_index)
   call object_modified_by_index (vold_index)
   call object_modified_by_index (div_vo_index)
   eloc2 = eold(1)

!  calculating numerical derivative
   deloc_rn_num (force_i) = (eloc2 - eloc1)/d
   
!  restoring the nuclear coordinates and nuclear potential energy
   cent (dim_i, cent_i) = cent_save
   call object_modified ('cent')
   call pot_nn(cent,znuc,iwctype,ncent,pecent)

  enddo ! force_i
  
  end subroutine deloc_rn_num_bld

! ==============================================================================
  subroutine forces_pulay_av_bld
! ------------------------------------------------------------------------------
! Description   : Averaged Pulay correction to forces
!
! Created       : J. Toulouse, 27 Jul 2008
! ------------------------------------------------------------------------------
  implicit none


! header
  if (header_exe) then

   call object_create ('forces_pulay_av')

   call object_needed ('forces_nb')
   call object_needed ('dpsi_rn_eloc_av')
   call object_needed ('dpsi_rn_av')
   call object_needed ('eloc_av')

   return

  endif

! begin

! allocations
  call object_alloc ('forces_pulay_av', forces_pulay_av, forces_nb)

  forces_pulay_av (:) = -2.d0 * (dpsi_rn_eloc_av (:) - dpsi_rn_av (:) * eloc_av)

  end subroutine forces_pulay_av_bld

! ==============================================================================
  subroutine forces_pulay_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of averaged Pulay correction to forces
!
! Created       : J. Toulouse, 27 Jul 2008
! ------------------------------------------------------------------------------
  implicit none


! header
  if (header_exe) then

   call object_create ('forces_pulay_av_var')
   call object_error_define_from_variance ('forces_pulay_av_var', 'forces_pulay_av_err')

   call object_needed ('forces_nb')
   call object_needed ('dpsi_rn_eloc_av_var')
   call object_needed ('eloc_av')
   call object_needed ('dpsi_rn_av_var')
   call object_needed ('dpsi_rn_av')
   call object_needed ('eloc_av_var')
   call object_needed ('dpsi_rn_eloc_av_dpsi_rn_av_covar')
   call object_needed ('dpsi_rn_eloc_av_eloc_av_covar')
   call object_needed ('dpsi_rn_av_eloc_av_covar')

   return

  endif

! begin

! allocations
  call object_alloc ('forces_pulay_av_var', forces_pulay_av_var, forces_nb)
  call object_alloc ('forces_pulay_av_err', forces_pulay_av_err, forces_nb)

  forces_pulay_av_var (:) = 4.d0*dpsi_rn_eloc_av_var (:) + 4.d0*(eloc_av**2)*dpsi_rn_av_var (:) &
                          + 4.d0*(dpsi_rn_av (:)**2)*eloc_av_var - 8.d0*eloc_av*dpsi_rn_eloc_av_dpsi_rn_av_covar (:) &
                          - 8.d0*dpsi_rn_av (:)*dpsi_rn_eloc_av_eloc_av_covar (:) +8.d0*dpsi_rn_av (:)*eloc_av*dpsi_rn_av_eloc_av_covar (:)

  end subroutine forces_pulay_av_var_bld

end module forces_pulay_mod
