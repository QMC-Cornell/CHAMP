module deriv_orb_mod

  use all_tools_mod
  use montecarlo_mod
  use orbitals_mod
  use determinants_mod
  use csfs_mod
  use eloc_mod

! Declaration of global variables and default values

  integer                                :: param_orb_nb = 0
  integer                                :: single_ex_nb = 0
  integer                                :: param_orb_nb_old = 0  ! temp
  integer                                :: single_ex_nb_old = 0  ! temp
  integer                                :: deriv_orb_pairs_nb = 0

! new
  integer                                :: csf_ex_unq_nb = 0
  integer, allocatable                   :: csf_ex_unq_ref (:)

  integer, allocatable                   :: det_ex_unq_ref_up (:)
  integer, allocatable                   :: det_ex_unq_ref_dn (:)
  integer, allocatable                   :: det_unq_in_csf_ex_unq_nb (:)
  type (type_integer_row), allocatable   :: det_unq_up_in_csf_ex_unq (:)
  type (type_integer_row), allocatable   :: det_unq_dn_in_csf_ex_unq (:)
  type (type_real_row), allocatable      :: cdet_unq_in_csf_ex_unq (:)

  integer, allocatable                   :: csf_unq_in_wf_ex_nb (:)
  type (type_integer_row), allocatable   :: csf_unq_in_wf_ex (:)
  type (type_integer_row), allocatable   :: csf_unq_ref_in_wf_ex (:)
  type (type_real_row), allocatable      :: csf_unq_prefac_in_wf_ex (:)

  integer, allocatable                   :: csf_unq_in_dpsi_orb_nb (:)
  type (type_integer_row), allocatable   :: csf_unq_in_dpsi_orb (:)
  type (type_integer_row), allocatable   :: csf_unq_ref_in_dpsi_orb (:)
  type (type_real_row), allocatable      :: csf_unq_prefac_in_dpsi_orb (:)

  integer, allocatable                   :: ex_orb_ind (:)
  integer, allocatable                   :: ex_orb_ind_rev (:)
  integer, allocatable                   :: ex_orb_1st_lab (:)
  integer, allocatable                   :: ex_orb_2nd_lab (:)
  integer, allocatable                   :: ex_orb_ind_old (:) ! temp
  integer, allocatable                   :: ex_orb_ind_rev_old (:) ! temp
  integer, allocatable                   :: ex_orb_1st_lab_old (:) ! temp
  integer, allocatable                   :: ex_orb_2nd_lab_old (:) ! temp

  logical, allocatable                   :: is_ex_act_act (:)
  integer, allocatable                   :: ex_act_act (:)
  logical, allocatable                   :: ex_up_is_zero (:,:)
  logical, allocatable                   :: ex_dn_is_zero (:,:)
  integer, allocatable                   :: det_ex_orb_lab_up (:,:,:)
  integer, allocatable                   :: det_ex_orb_lab_dn (:,:,:)
  integer, allocatable                   :: det_ex_orb_lab_srt_up (:,:,:)
  integer, allocatable                   :: det_ex_orb_lab_srt_dn (:,:,:)
  integer, allocatable                   :: det_ex_orb_lab_srt_sgn_up (:,:)
  integer, allocatable                   :: det_ex_orb_lab_srt_sgn_dn (:,:)
  logical, allocatable                   :: is_det_ex_up (:,:)
  logical, allocatable                   :: is_det_ex_dn (:,:)
  integer, allocatable                   :: iwdet_ex_up (:,:)
  integer, allocatable                   :: iwdet_ex_dn (:,:)
  integer                                :: det_ex_unq_up_nb
  integer                                :: det_ex_unq_dn_nb
  integer, allocatable                   :: det_ex_unq_up_orb_1st_pos (:)
  integer, allocatable                   :: det_ex_unq_dn_orb_1st_pos (:)
  integer, allocatable                   :: det_ex_unq_up_orb_2nd_lab (:)
  integer, allocatable                   :: det_ex_unq_dn_orb_2nd_lab (:)
  integer, allocatable                   :: iwdet_ex_ref_up (:)
  integer, allocatable                   :: iwdet_ex_ref_dn (:)
  integer, allocatable                   :: det_ex_unq_orb_lab_up (:,:)
  integer, allocatable                   :: det_ex_unq_orb_lab_dn (:,:)
  integer, allocatable                   :: det_ex_unq_orb_lab_srt_up (:,:)
  integer, allocatable                   :: det_ex_unq_orb_lab_srt_dn (:,:)
  integer, allocatable                   :: det_ex_unq_orb_lab_srt_sgn_up (:)
  integer, allocatable                   :: det_ex_unq_orb_lab_srt_sgn_dn (:)
  real(dp), allocatable                  :: det_ex_unq_up (:)
  real(dp), allocatable                  :: det_ex_unq_dn (:)
  real(dp), allocatable                  :: det_ex_unq_up_2 (:)
  real(dp), allocatable                  :: det_ex_unq_dn_2 (:)
  real(dp), allocatable                  :: det_ex_up (:,:)
  real(dp), allocatable                  :: det_ex_dn (:,:)
  real(dp), allocatable                  :: det_ex (:,:)
  real(dp), allocatable                  :: psid_ex (:)
  real(dp), allocatable                  :: slater_mat_ex_trans_up (:,:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_dn (:,:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_up (:, :, :)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_dn (:, :, :)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_up_2 (:,:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_dn_2 (:,:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_up_3 (:,:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_dn_3 (:,:,:)
  real(dp), allocatable                  :: mat_flat_up (:,:)
  real(dp), allocatable                  :: mat_flat_dn (:,:)
  real(dp), allocatable                  :: grd_det_ex_unq_dn (:,:,:)
  real(dp), allocatable                  :: grd_det_ex_unq_up (:,:,:)
  real(dp), allocatable                  :: lap_det_ex_unq_dn (:,:)
  real(dp), allocatable                  :: lap_det_ex_unq_up (:,:)
  real(dp), allocatable                  :: grd_det_ex_dn (:,:,:,:)
  real(dp), allocatable                  :: grd_det_ex_up (:,:,:,:)
  real(dp), allocatable                  :: lap_det_ex_dn (:,:,:)
  real(dp), allocatable                  :: lap_det_ex_up (:,:,:)
  real(dp), allocatable                  :: grd_psid_ex_over_psid (:,:,:)
  real(dp), allocatable                  :: lap_psid_ex_over_psid (:,:)
  real(dp), allocatable                  :: lap_lnpsid_ex (:,:)
  real(dp), allocatable                  :: sum_lap_lnpsid_ex (:)
  real(dp), allocatable                  :: sum_lap_lnj_ex (:)
  real(dp), allocatable                  :: sum_lap_lnpsi_ex (:)
  real(dp), allocatable                  :: grd_psi_ex_over_psi (:,:,:)
  real(dp), allocatable                  :: eloc_kin_ex (:)
  real(dp), allocatable                  :: eloc_pot_nloc_ex (:)
  real(dp), allocatable                  :: eloc_pot_ex (:)
  real(dp), allocatable                  :: eloc_ex (:)
  real(dp), allocatable                  :: vpot_ex (:,:)
  real(dp), allocatable                  :: vpsp_ex (:)
  real(dp), allocatable                  :: psid_ex_in_x (:)
  integer                                :: electron

  real(dp), allocatable                  :: dpsi_orb(:)
  real(dp), allocatable                  :: dpsi_orb_test(:)

  real(dp), allocatable                  :: deloc_orb (:)

  real(dp), allocatable                  :: delta_eps (:)

  contains

! ==============================================================================
  subroutine single_ex_wf_bld
! ------------------------------------------------------------------------------
! Description   : build single orbital excitation list for general wave function
! Description   : ex_orb_ind : index of direct excitation Eij
! Description   : ex_orb_ind_rev : index of reverse excitation Eji (for active-active excitations)
!
! Created       : J. Toulouse, 12 Oct 2005
! Revised       : J. Toulouse, 24 Oct 2005: open shells
! Revised       : J. Toulouse, 23 Mar 2006: active-active excitations for non-CASSCF wave functions
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_i, elec_up_i
  integer orb_opt_cls_i, orb_opt_opn_i, orb_opt_act_i, orb_opt_act_j, orb_opt_vir_i
  integer orb_opt_cls, orb_opt_opn, orb_opt_act, orb_opt_vir, orb_opt_act_lab_i, orb_opt_act_lab_j

! header
  if (header_exe) then

   call object_create ('param_orb_nb')
   call object_create ('single_ex_nb')
   call object_create ('deriv_orb_pairs_nb')
   call object_create ('ex_orb_ind')
   call object_create ('ex_orb_ind_rev')
   call object_create ('ex_orb_1st_lab')
   call object_create ('ex_orb_2nd_lab')

   call object_needed ('orb_opt_nb')
   call object_needed ('orb_opt_lab')
   call object_needed ('orb_opt_occ_nb')
   call object_needed ('orb_opt_occ_lab')
   call object_needed ('orb_opt_cls_nb')
   call object_needed ('orb_opt_cls_lab')
   call object_needed ('orb_opt_opn_nb')
   call object_needed ('orb_opt_opn_lab')
   call object_needed ('orb_opt_act_nb')
   call object_needed ('orb_opt_act_lab')
   call object_needed ('orb_opt_vir_nb')
   call object_needed ('orb_opt_vir_lab')
   call object_needed ('orb_sym_lab')
   call object_needed ('orb_ex_forbidden')

   return

  endif

! begin
  write(6,'(a)') ' Constructing list of single orbital excitations based on orbital occupancies alone...'

  single_ex_nb = 0
  param_orb_nb = 0

! single-determinant wave function
  if (ndet == 1) then

! closed -> open
  do orb_opt_cls_i = 1, orb_opt_cls_nb
     orb_opt_cls = orb_opt_cls_lab (orb_opt_cls_i)
    do orb_opt_opn_i = 1, orb_opt_opn_nb
       orb_opt_opn = orb_opt_opn_lab (orb_opt_opn_i)
       if (orb_sym_lab (orb_opt_cls) == orb_sym_lab (orb_opt_opn) .and. .not. orb_ex_forbidden (orb_opt_cls, orb_opt_opn)) then
       param_orb_nb = param_orb_nb + 1
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
       call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_ind (param_orb_nb) = single_ex_nb
       ex_orb_ind_rev (param_orb_nb) = 0
       ex_orb_1st_lab (single_ex_nb) = orb_opt_cls
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_opn
       endif
     enddo ! orb_opt_opn_i
  enddo ! orb_opt_cls_i

! closed -> virtual
  do orb_opt_cls_i = 1, orb_opt_cls_nb
       orb_opt_cls = orb_opt_cls_lab (orb_opt_cls_i)
    do orb_opt_vir_i = 1, orb_opt_vir_nb
       orb_opt_vir = orb_opt_vir_lab (orb_opt_vir_i)
       if (orb_sym_lab (orb_opt_cls) == orb_sym_lab (orb_opt_vir) .and. .not. orb_ex_forbidden (orb_opt_cls, orb_opt_vir)) then
       param_orb_nb = param_orb_nb + 1
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
       call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_ind (param_orb_nb) = single_ex_nb
       ex_orb_ind_rev (param_orb_nb) = 0
       ex_orb_1st_lab (single_ex_nb)  = orb_opt_cls
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_vir
       endif
     enddo ! orb_opt_vir_i
  enddo ! orb_opt_cls_i

! open -> virtual
  do orb_opt_opn_i = 1, orb_opt_opn_nb
       orb_opt_opn = orb_opt_opn_lab (orb_opt_opn_i)
    do orb_opt_vir_i = 1, orb_opt_vir_nb
       orb_opt_vir = orb_opt_vir_lab (orb_opt_vir_i)
       if (orb_sym_lab (orb_opt_opn) == orb_sym_lab (orb_opt_vir) .and. .not. orb_ex_forbidden (orb_opt_opn, orb_opt_vir)) then
       param_orb_nb = param_orb_nb + 1
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
       call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_ind (param_orb_nb) = single_ex_nb
       ex_orb_ind_rev (param_orb_nb) = 0
       ex_orb_1st_lab (single_ex_nb)  = orb_opt_opn
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_vir
       endif
     enddo ! orb_opt_vir_i
  enddo ! orb_opt_opn_i

  else

! multi-determinantal wave function

! closed -> active
  do orb_opt_cls_i = 1, orb_opt_cls_nb
    orb_opt_cls = orb_opt_cls_lab (orb_opt_cls_i)
    do orb_opt_act_i = 1, orb_opt_act_nb
       orb_opt_act = orb_opt_act_lab (orb_opt_act_i)
       if (orb_sym_lab (orb_opt_cls) == orb_sym_lab (orb_opt_act) .and. .not. orb_ex_forbidden (orb_opt_cls, orb_opt_act)) then
       param_orb_nb = param_orb_nb + 1
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
       call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_ind (param_orb_nb) = single_ex_nb
       ex_orb_ind_rev (param_orb_nb) = 0
       ex_orb_1st_lab (single_ex_nb) = orb_opt_cls
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_act
       endif
     enddo ! orb_opt_act_i
  enddo ! orb_opt_cls_i

! closed -> virtual
  do orb_opt_cls_i = 1, orb_opt_cls_nb
    orb_opt_cls = orb_opt_cls_lab (orb_opt_cls_i)
    do orb_opt_vir_i = 1, orb_opt_vir_nb
       orb_opt_vir = orb_opt_vir_lab (orb_opt_vir_i)
       if (orb_sym_lab (orb_opt_cls) == orb_sym_lab (orb_opt_vir) .and. .not. orb_ex_forbidden (orb_opt_cls, orb_opt_vir)) then
       param_orb_nb = param_orb_nb + 1
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
       call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_ind (param_orb_nb) = single_ex_nb
       ex_orb_ind_rev (param_orb_nb) = 0
       ex_orb_1st_lab (single_ex_nb) = orb_opt_cls
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_vir
       endif
     enddo ! orb_opt_vir_i
  enddo ! orb_opt_cls_i

! active -> active (these excitations are redundant for CASSCF wave function)
  if (.not. l_casscf) then
   do orb_opt_act_i = 1, orb_opt_act_nb
    orb_opt_act_lab_i = orb_opt_act_lab (orb_opt_act_i)
    do orb_opt_act_j = orb_opt_act_i+1, orb_opt_act_nb
       orb_opt_act_lab_j = orb_opt_act_lab (orb_opt_act_j)
       if (orb_sym_lab (orb_opt_act_lab_i) == orb_sym_lab (orb_opt_act_lab_j) .and. .not. orb_ex_forbidden (orb_opt_act_lab_i, orb_opt_act_lab_j)) then

       param_orb_nb = param_orb_nb + 1
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
       call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_ind (param_orb_nb) = single_ex_nb
       ex_orb_1st_lab (single_ex_nb) = orb_opt_act_lab_i
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_act_lab_j

! reverse excitation (j->i)
!       ex_orb_ind_rev (param_orb_nb) = 0 ! desactivate reverse excitation
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_ind_rev (param_orb_nb) = single_ex_nb
       ex_orb_1st_lab (single_ex_nb) = orb_opt_act_lab_j
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_act_lab_i

       endif
     enddo ! orb_opt_act_j
   enddo ! orb_opt_act_i
  endif

! active -> virtual
  do orb_opt_act_i = 1, orb_opt_act_nb
    orb_opt_act = orb_opt_act_lab (orb_opt_act_i)
    do orb_opt_vir_i = 1, orb_opt_vir_nb
       orb_opt_vir = orb_opt_vir_lab (orb_opt_vir_i)
       if (orb_sym_lab (orb_opt_act) == orb_sym_lab (orb_opt_vir) .and. .not. orb_ex_forbidden (orb_opt_act, orb_opt_vir)) then
       param_orb_nb = param_orb_nb + 1
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
       call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_ind (param_orb_nb) = single_ex_nb
       ex_orb_ind_rev (param_orb_nb) = 0
       ex_orb_1st_lab (single_ex_nb) = orb_opt_act
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_vir
       endif
     enddo ! orb_opt_act_i
  enddo ! orb_opt_act_i

  endif

  write(6,'(a,i10)') ' Number of single orbital excitations = ', single_ex_nb
  write(6,'(a,i10)') ' Number of orbital derivatives        = ', param_orb_nb

! excitation pairs number
  deriv_orb_pairs_nb = param_orb_nb * (param_orb_nb + 1) / 2

! temporary: new routine
  if (l_check_redundant_orbital_derivative) then

! temporary: save created objects for comparison with new routine
  param_orb_nb_old = param_orb_nb
  single_ex_nb_old = single_ex_nb
  call copy (ex_orb_ind, ex_orb_ind_old)
  call copy (ex_orb_ind_rev, ex_orb_ind_rev_old)
  call copy (ex_orb_1st_lab, ex_orb_1st_lab_old)
  call copy (ex_orb_2nd_lab, ex_orb_2nd_lab_old)

! temporary: call new routine checking for redundancies
  call node_exe ('single_ex_wf_bld_2')

  endif

 end subroutine single_ex_wf_bld

! ==============================================================================
  subroutine single_ex_wf_bld_2
! ------------------------------------------------------------------------------
! Description   : Build single orbital excitation list for general wave function
! Description   : taking into account some (but not all) redundancies
! Description   : For large CASCF wave functions, this routine does not seem to recognize redundant active-> active excitations
! Description   : recognize redundant active-> active excitations -> has to be improved
! Description   : This routine is to replace single_ex_wf_bld
! Description   : For now, the routine is only used as replacement for single_ex_wf_bld
! Description   : information gathered on excitations in this routine should be used in the program
! Description   : in lieu of single_ex_det_bld
!
! Created       : J. Toulouse, 25 Oct 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer orb_opt_i, orb_opt_j, orb_opt_lab_i, orb_opt_lab_j
  integer ex_dir_rev
  logical dpsi_orb_is_zero, ex_is_zero, csf_ex_unq_is_zero, dpsi_orb_is_redundant

  integer det_ex_unq_up_k, det_ex_unq_dn_k
  integer det_unq_cur_up, det_ex_unq_sgn_cur_up
  integer det_unq_cur_dn, det_ex_unq_sgn_cur_dn
  real(dp) cdet_unq_cur
  integer, allocatable :: det_ex_cur_orb_lab_up (:)
  integer, allocatable :: det_ex_cur_orb_lab_srt_up (:)
  integer, allocatable :: det_ex_cur_orb_lab_dn (:)
  integer, allocatable :: det_ex_cur_orb_lab_srt_dn (:)
  integer det_ex_cur_orb_lab_srt_sgn_up, det_ex_cur_orb_lab_srt_sgn_dn

  integer det_in_csf_i, det_unq_up_i, det_unq_dn_i
  integer det_unq_up_k, det_unq_dn_k

  integer det_unq_in_csf_ex_cur_nb
  integer, allocatable :: det_unq_up_in_csf_ex_cur (:)
  integer, allocatable :: det_unq_dn_in_csf_ex_cur (:)
  real(dp), allocatable :: cdet_unq_in_csf_ex_cur (:)
  integer det_unq_in_csf_ex_cur_temp_nb, det_unq_in_csf_ex_cur_i
  integer, allocatable :: det_unq_up_in_csf_ex_cur_temp (:)
  integer, allocatable :: det_unq_dn_in_csf_ex_cur_temp (:)
  real(dp), allocatable :: cdet_unq_in_csf_ex_cur_temp (:)

  integer csf_i, csf_k, csf_ex_unq_i, csf_cur
  integer csf_opt_k, csf_opt_k_lab
  integer csf_unq_in_wf_ex_cur_nb
  integer, allocatable :: csf_unq_in_wf_ex_cur (:)
  integer, allocatable :: csf_unq_ref_in_wf_ex_cur (:)
  real(dp), allocatable :: csf_unq_prefac_in_wf_ex_cur (:)
  real(dp) :: csf_ex_cur_prefac

  integer single_ex_added_nb
  integer ex_i, ex_cur
  integer ex_orb_ind_cur, ex_orb_ind_rev_cur
  integer dpsi_orb_i
  integer csf_unq_in_dpsi_orb_cur_nb, csf_unq_in_dpsi_orb_cur_i
  integer, allocatable :: csf_unq_in_dpsi_orb_cur (:)
  integer, allocatable :: csf_unq_ref_in_dpsi_orb_cur (:)
  real(dp), allocatable :: csf_unq_prefac_in_dpsi_orb_cur (:)
  integer csf_unq_in_dpsi_orb_cur_temp_nb
  integer, allocatable :: csf_unq_in_dpsi_orb_cur_temp (:)
  integer, allocatable :: csf_unq_ref_in_dpsi_orb_cur_temp (:)
  real(dp), allocatable :: csf_unq_prefac_in_dpsi_orb_cur_temp (:)

! tenp
   integer ex_j, orb_1st, orb_2nd
   logical ex_found

! header
  if (header_exe) then

!   call object_create ('det_ex_unq_up_nb')
!   call object_create ('det_ex_unq_ref_up')
!   call object_create ('det_ex_unq_orb_lab_srt_up')
!   call object_create ('det_ex_unq_orb_lab_srt_sgn_up')
!   call object_create ('det_ex_unq_dn_nb')
!   call object_create ('det_ex_unq_ref_dn')
!   call object_create ('det_ex_unq_orb_lab_srt_dn')
!   call object_create ('det_ex_unq_orb_lab_srt_sgn_dn')
!   call object_create ('csf_ex_unq_nb')
!   call object_create ('csf_ex_unq_ref')
!   call object_create ('det_unq_in_csf_ex_unq_nb')
!   call object_create ('det_unq_up_in_csf_ex_unq')
!   call object_create ('det_unq_dn_in_csf_ex_unq')
!   call object_create ('cdet_unq_in_csf_ex_unq')

!   call object_create ('single_ex_nb')
!   call object_create ('csf_unq_in_wf_ex_nb')
!   call object_create ('csf_unq_in_wf_ex')
!   call object_create ('csf_unq_ref_in_wf_ex')
!   call object_create ('csf_unq_prefac_in_wf_ex')
!   call object_create ('ex_orb_ind')
!   call object_create ('ex_orb_ind_rev')
!   call object_create ('ex_orb_1st_lab')
!   call object_create ('ex_orb_2nd_lab')

!   call object_create ('param_orb_nb')
!   call object_create ('deriv_orb_pairs_nb')
!   call object_create ('csf_unq_in_dpsi_orb_nb')
!   call object_create ('csf_unq_in_dpsi_orb')
!   call object_create ('csf_unq_ref_in_dpsi_orb')
!   call object_create ('csf_unq_prefac_in_dpsi_orb')

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ncsf')
   call object_needed ('ndet')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('nparmcsf')
   call object_needed ('iwcsf')
   call object_needed ('orb_opt_nb')
   call object_needed ('orb_opt_lab')
   call object_needed ('orb_sym_lab')
   call object_needed ('orb_occ_in_wf')
   call object_needed ('orb_cls_in_wf')
   call object_needed ('orb_act_in_wf')
   call object_needed ('orb_opn_in_wf')
   call object_needed ('orb_vir_in_wf')
   call object_needed ('orb_ex_forbidden')

   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('orb_occ_in_wf')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')

   call object_needed ('det_unq_up_in_csf')
   call object_needed ('det_unq_dn_in_csf')
   call object_needed ('cdet_unq_in_csf')

   return

  endif

! begin
  write(6,'(a)') ' Constructing list of single orbital excitations checking for additional redundancies...'

! allocations
  call alloc ('det_ex_cur_orb_lab_up', det_ex_cur_orb_lab_up, nup)
  call alloc ('det_ex_cur_orb_lab_srt_up', det_ex_cur_orb_lab_srt_up, nup)
  call alloc ('det_ex_cur_orb_lab_dn', det_ex_cur_orb_lab_dn, ndn)
  call alloc ('det_ex_cur_orb_lab_srt_dn', det_ex_cur_orb_lab_srt_dn, ndn)

! initializations
  single_ex_nb = 0
  param_orb_nb = 0
  ex_orb_ind_rev_cur = 0

! loops over orbitals
  do orb_opt_i = 1, orb_opt_nb
   orb_opt_lab_i = orb_opt_lab (orb_opt_i)

   do orb_opt_j = orb_opt_i+1, orb_opt_nb
    orb_opt_lab_j = orb_opt_lab (orb_opt_j)

!    write(6,'(2a,i4,a,i4)') trim(here),': considering orbital excitation ', orb_opt_lab_i, ' -> ', orb_opt_lab_j
    dpsi_orb_is_zero = .true.

!   skip user-specified forbidden excitations
    if (orb_ex_forbidden (orb_opt_lab_i, orb_opt_lab_j)) then
      cycle
    endif

!   excitations between orbitals of same symmetry only
    if (orb_sym_lab (orb_opt_lab_i) /= orb_sym_lab (orb_opt_lab_j)) then
      cycle
    endif

!   skip obvious redundancies:
!   skip closed -> closed excitations
    if (orb_cls_in_wf (orb_opt_lab_i) .and. orb_cls_in_wf (orb_opt_lab_j)) then
      cycle
    endif

!   skip virtual -> virtual excitations
    if (orb_vir_in_wf (orb_opt_lab_i) .and. orb_vir_in_wf (orb_opt_lab_j)) then
      cycle
    endif

!   skip open -> open excitations for single-determinant wave functions
    if (ndet == 1 .and. (orb_opn_in_wf (orb_opt_lab_i) .and. orb_opn_in_wf (orb_opt_lab_j))) then
      cycle
    endif

!   skip active -> active excitations for CASSCF wave functions
    if (l_casscf .and. (orb_act_in_wf (orb_opt_lab_i) .and. orb_act_in_wf (orb_opt_lab_j))) then
      cycle
    endif

!   loop over the direct excitation i->j and the reverse excitation j->i
    do ex_dir_rev = 1, 2

    ex_is_zero = .true.

!   swap orbitals for reverse excitation j->i
    if (ex_dir_rev == 2) then
     call swap (orb_opt_lab_i, orb_opt_lab_j)
    endif

!   skip if orbital orb_opt_lab_i is not occupied in the wave function
    if (.not. orb_occ_in_wf (orb_opt_lab_i)) then
      cycle
    endif

!   skip if orbital orb_opt_lab_j is closed in the wave function
    if (orb_cls_in_wf (orb_opt_lab_j)) then
      cycle
    endif

!   temporary: only active -> active excitations
!    if ( .not. (orb_act_in_wf (orb_opt_lab_i) .and.  orb_act_in_wf (orb_opt_lab_j))) then
!      cycle
!    endif

!    write(6,'(2a)') trim(here),': --------------------------------------------------------------------------------------'
!    write(6,'(2a,i4,a,i4)') trim(here),': considering orbital excitation ', orb_opt_lab_i, ' -> ', orb_opt_lab_j


!   construct excitations from ground-state wave-function
    do csf_i = 1, ncsf

!     write(6,'(2a,i4)') trim(here),': consider CSF # ',csf_i

     csf_ex_unq_is_zero = .true.

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

!      indexes unique spin-up and spin-dn determinants
       det_unq_up_i = det_unq_up_in_csf (csf_i)%row (det_in_csf_i)
       det_unq_dn_i = det_unq_dn_in_csf (csf_i)%row (det_in_csf_i)

!      spin-up excited determinants:
!      if orb_opt_lab_i is occupied in this deterninant and orb_opt_lab_j is not occupied
       if (elt_in_array (det_unq_orb_lab_srt_up (:, det_unq_up_i), orb_opt_lab_i) .and.  &
           .not. elt_in_array (det_unq_orb_lab_srt_up (:, det_unq_up_i), orb_opt_lab_j)) then

        csf_ex_unq_is_zero = .false.

!       build current excited determinant
        det_ex_cur_orb_lab_up (:) = det_unq_orb_lab_srt_up (:, det_unq_up_i)
        call replace_elt_in_array (det_ex_cur_orb_lab_up (:), orb_opt_lab_i, orb_opt_lab_j)
!        write(6,'(2a,100i4)') trim(here),': construct excited spin-up determinant:',det_ex_cur_orb_lab_up

!       sort current excited determinant
        det_ex_cur_orb_lab_srt_up (:) = det_ex_cur_orb_lab_up (:)
        call sort_and_sign (det_ex_cur_orb_lab_srt_up (:), det_ex_cur_orb_lab_srt_sgn_up)

!       check if current excited determinant is a determinant in the ground-state wave function
        det_unq_cur_up = 0
        do det_unq_up_k = 1, ndetup
         if (arrays_equal (det_ex_cur_orb_lab_srt_up (:), det_unq_orb_lab_srt_up (:, det_unq_up_k))) then
           det_unq_cur_up = det_unq_up_k
           det_ex_unq_sgn_cur_up = det_ex_cur_orb_lab_srt_sgn_up
!           write(6,'(2a)') trim(here), ': spin-up determinant is a ground-state spin-up determinant'
          exit
         endif
        enddo

!       check if current excited determinant is a excited determinant already encountered
        if (det_unq_cur_up == 0) then
         do det_ex_unq_up_k = 1, det_ex_unq_up_nb
          if (arrays_equal (det_ex_cur_orb_lab_srt_up (:), det_ex_unq_orb_lab_srt_up (:, det_ex_unq_up_k))) then
           det_unq_cur_up =  ndetup + det_ex_unq_up_k
           det_ex_unq_sgn_cur_up = det_ex_cur_orb_lab_srt_sgn_up * det_ex_unq_orb_lab_srt_sgn_up (det_ex_unq_up_k)
!           write(6,'(2a)') trim(here), ': spin-up determinant is a excited spin-up determinant already encountered'
           exit
          endif
         enddo
        endif

!       if current excited determinant is a new determinant, add it to the list of unique excited determinants
        if (det_unq_cur_up == 0) then
         det_ex_unq_up_nb = det_ex_unq_up_nb + 1
         det_unq_cur_up = ndetup + det_ex_unq_up_nb
         det_ex_unq_sgn_cur_up = 1
         call object_alloc ('det_ex_unq_ref_up', det_ex_unq_ref_up, det_ex_unq_up_nb)
         call object_alloc ('det_ex_unq_orb_lab_srt_up', det_ex_unq_orb_lab_srt_up, nup, det_ex_unq_up_nb)
         call object_alloc ('det_ex_unq_orb_lab_srt_sgn_up', det_ex_unq_orb_lab_srt_sgn_up, det_ex_unq_up_nb)
         det_ex_unq_ref_up (det_ex_unq_up_nb) = det_unq_up_i
         det_ex_unq_orb_lab_srt_up (:, det_ex_unq_up_nb) = det_ex_cur_orb_lab_srt_up (:)
         det_ex_unq_orb_lab_srt_sgn_up (det_ex_unq_up_nb) = det_ex_cur_orb_lab_srt_sgn_up
!         write(6,'(2a,100i4)') trim(here), ': add new excited spin-up determinant:', det_ex_cur_orb_lab_srt_up (:)
        endif

!       compute coefficient of current excited determinant with sign
        cdet_unq_cur = cdet_unq_in_csf (csf_i)%row (det_in_csf_i) * det_ex_unq_sgn_cur_up

!       add this excited determinant to the current excited csf
        call append (det_unq_up_in_csf_ex_cur, det_unq_cur_up)
        call append (det_unq_dn_in_csf_ex_cur, det_unq_dn_i)
        call append (cdet_unq_in_csf_ex_cur, cdet_unq_cur)

       endif ! spin-up determinants

!      spin-dn excited determinants:
!      if orb_opt_lab_i is occupied in this deterninant and orb_opt_lab_j is not occupied
       if (elt_in_array (det_unq_orb_lab_srt_dn (:, det_unq_dn_i), orb_opt_lab_i) .and.  &
           .not. elt_in_array (det_unq_orb_lab_srt_dn (:, det_unq_dn_i), orb_opt_lab_j)) then

        csf_ex_unq_is_zero = .false.

!       build current excited determinant
        det_ex_cur_orb_lab_dn (:) = det_unq_orb_lab_srt_dn (:, det_unq_dn_i)
        call replace_elt_in_array (det_ex_cur_orb_lab_dn (:), orb_opt_lab_i, orb_opt_lab_j)
!        write(6,'(2a,100i4)') trim(here),': construct excited spin-up determinant:',det_ex_cur_orb_lab_dn

!       sort current excited determinant
        det_ex_cur_orb_lab_srt_dn (:) = det_ex_cur_orb_lab_dn (:)
        call sort_and_sign (det_ex_cur_orb_lab_srt_dn (:), det_ex_cur_orb_lab_srt_sgn_dn)

!       check if current excited determinant is a determinant in the ground-state wave function
        det_unq_cur_dn = 0
        do det_unq_dn_k = 1, ndetdn
         if (arrays_equal (det_ex_cur_orb_lab_srt_dn (:), det_unq_orb_lab_srt_dn (:, det_unq_dn_k))) then
           det_unq_cur_dn = det_unq_dn_k
           det_ex_unq_sgn_cur_dn = det_ex_cur_orb_lab_srt_sgn_dn
!           write(6,'(2a)') trim(here), ': spin-dn determinant is a ground-state spin-dn determinant'
          exit
         endif
        enddo

!       check if current excited determinant is a excited determinant already encountered
        if (det_unq_cur_dn == 0) then
         do det_ex_unq_dn_k = 1, det_ex_unq_dn_nb
          if (arrays_equal (det_ex_cur_orb_lab_srt_dn (:), det_ex_unq_orb_lab_srt_dn (:, det_ex_unq_dn_k))) then
           det_unq_cur_dn =  ndetdn + det_ex_unq_dn_k
           det_ex_unq_sgn_cur_dn = det_ex_cur_orb_lab_srt_sgn_dn * det_ex_unq_orb_lab_srt_sgn_dn (det_ex_unq_dn_k)
!           write(6,'(2a)') trim(here), ': spin-dn determinant is a excited spin-dn determinant already encountered'
           exit
          endif
         enddo
        endif

!       if current excited determinant is a new determinant, add it to the list of unique excited determinants
        if (det_unq_cur_dn == 0) then
         det_ex_unq_dn_nb = det_ex_unq_dn_nb + 1
         det_unq_cur_dn = ndetdn + det_ex_unq_dn_nb
         det_ex_unq_sgn_cur_dn = 1
         call object_alloc ('det_ex_unq_ref_dn', det_ex_unq_ref_dn, det_ex_unq_dn_nb)
         call object_alloc ('det_ex_unq_orb_lab_srt_dn', det_ex_unq_orb_lab_srt_dn, ndn, det_ex_unq_dn_nb)
         call object_alloc ('det_ex_unq_orb_lab_srt_sgn_dn', det_ex_unq_orb_lab_srt_sgn_dn, det_ex_unq_dn_nb)
         det_ex_unq_ref_dn (det_ex_unq_dn_nb) = det_unq_dn_i
         det_ex_unq_orb_lab_srt_dn (:, det_ex_unq_dn_nb) = det_ex_cur_orb_lab_srt_dn (:)
         det_ex_unq_orb_lab_srt_sgn_dn (det_ex_unq_dn_nb) = det_ex_cur_orb_lab_srt_sgn_dn
!         write(6,'(2a,100i4)') trim(here), ': add new excited spin-dn determinant:', det_ex_cur_orb_lab_srt_dn (:)
        endif

!       compute coefficient of current excited determinant with sign
        cdet_unq_cur = cdet_unq_in_csf (csf_i)%row (det_in_csf_i) * det_ex_unq_sgn_cur_dn

!       add this excited determinant to the current excited csf
        call append (det_unq_up_in_csf_ex_cur, det_unq_up_i)
        call append (det_unq_dn_in_csf_ex_cur, det_unq_cur_dn)
        call append (cdet_unq_in_csf_ex_cur, cdet_unq_cur)

       endif ! spin-dn determinants

     enddo ! det_in_csf_i

!    current excited csf
     if (.not. csf_ex_unq_is_zero) then

      ex_is_zero = .false.

!     sort only unique determinants in current excited csf
      call sort (det_unq_up_in_csf_ex_cur (:), det_unq_dn_in_csf_ex_cur (:), cdet_unq_in_csf_ex_cur (:))
      det_unq_in_csf_ex_cur_nb = mysize (det_unq_up_in_csf_ex_cur)

!     combine unique determinants in current excited csf
      det_unq_in_csf_ex_cur_temp_nb = det_unq_in_csf_ex_cur_nb
      call copy (det_unq_up_in_csf_ex_cur, det_unq_up_in_csf_ex_cur_temp)
      call copy (det_unq_dn_in_csf_ex_cur, det_unq_dn_in_csf_ex_cur_temp)
      call copy (cdet_unq_in_csf_ex_cur, cdet_unq_in_csf_ex_cur_temp)
      call release ('det_unq_up_in_csf_ex_cur', det_unq_up_in_csf_ex_cur)
      call release ('det_unq_dn_in_csf_ex_cur', det_unq_dn_in_csf_ex_cur)
      call release ('cdet_unq_in_csf_ex_cur', cdet_unq_in_csf_ex_cur)

      det_unq_in_csf_ex_cur_nb = 1
      call append (det_unq_up_in_csf_ex_cur, det_unq_up_in_csf_ex_cur_temp (1))
      call append (det_unq_dn_in_csf_ex_cur, det_unq_dn_in_csf_ex_cur_temp (1))
      call append (cdet_unq_in_csf_ex_cur, cdet_unq_in_csf_ex_cur_temp (1))

      do det_unq_in_csf_ex_cur_i = 2, det_unq_in_csf_ex_cur_temp_nb
        if (det_unq_up_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i) == det_unq_up_in_csf_ex_cur (det_unq_in_csf_ex_cur_nb) .and.  &
            det_unq_up_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i) == det_unq_dn_in_csf_ex_cur (det_unq_in_csf_ex_cur_nb)) then
          cdet_unq_in_csf_ex_cur (det_unq_in_csf_ex_cur_nb) = cdet_unq_in_csf_ex_cur (det_unq_in_csf_ex_cur_nb)  &
                                                              + cdet_unq_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i)
        else
         det_unq_in_csf_ex_cur_nb = det_unq_in_csf_ex_cur_nb + 1
         call append (det_unq_up_in_csf_ex_cur, det_unq_up_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i))
         call append (det_unq_dn_in_csf_ex_cur, det_unq_dn_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i))
         call append (cdet_unq_in_csf_ex_cur, cdet_unq_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i))
        endif
      enddo


!     prefactor of csf = coefficient of first determinant
      csf_ex_cur_prefac = cdet_unq_in_csf_ex_cur (1)

!     normalize csf so that the coefficient of first determinant is 1
      cdet_unq_in_csf_ex_cur (:) = cdet_unq_in_csf_ex_cur (:) / csf_ex_cur_prefac

!      write(6,'(2a)') trim(here), ': construct excited csf:'
!      write(6,'(2a,100i4)') trim(here),': det_unq_up_in_csf_ex_cur=',det_unq_up_in_csf_ex_cur (:)
!      write(6,'(2a,100i4)') trim(here),': det_unq_dn_in_csf_ex_cur=',det_unq_dn_in_csf_ex_cur (:)
!      write(6,'(2a,100f7.3)')  trim(here),': cdet_unq_in_csf_ex_cur=', cdet_unq_in_csf_ex_cur (:)

!     check if current excited csf is a csf in the ground-state wave function
      csf_cur = 0
      do csf_k = 1, ncsf
       if (arrays_equal (det_unq_up_in_csf_ex_cur (:), det_unq_up_in_csf (csf_k)%row (:)) .and. &
           arrays_equal (det_unq_dn_in_csf_ex_cur (:), det_unq_dn_in_csf (csf_k)%row (:)) .and. &
           arrays_equal (cdet_unq_in_csf_ex_cur (:), cdet_unq_in_csf (csf_k)%row (:))) then
           csf_cur = csf_k
!           write(6,'(2a)') trim(here),': excited csf is ground-state csf'
       endif
      enddo ! csf_k

!     check if current excited csf is a linear combinaison of ground-state csfs
!     if (csf_cur == 0) then
!     endif

!     check if current excited csf is a excited csf already encountered
      if (csf_cur == 0) then
       do csf_ex_unq_i = 1, csf_ex_unq_nb
        if (arrays_equal (det_unq_up_in_csf_ex_cur (:), det_unq_up_in_csf_ex_unq (csf_ex_unq_i)%row (:)) .and. &
            arrays_equal (det_unq_dn_in_csf_ex_cur (:), det_unq_dn_in_csf_ex_unq (csf_ex_unq_i)%row (:)) .and. &
            arrays_equal (cdet_unq_in_csf_ex_cur (:), cdet_unq_in_csf_ex_unq (csf_ex_unq_i)%row (:))) then
            csf_cur = ncsf + csf_ex_unq_i
!            write(6,'(2a)') trim(here),': excited csf is excited csf already encountered'
         endif
       enddo ! csf_ex_unq_i
      endif


!    if current excited csf is a new excited csf, add it to the list of unique excited csf
     if (csf_cur == 0) then
      csf_ex_unq_nb = csf_ex_unq_nb + 1
      csf_cur = ncsf + csf_ex_unq_nb

!     reference csf
      call object_alloc ('csf_ex_unq_ref', csf_ex_unq_ref, csf_ex_unq_nb)
      csf_ex_unq_ref (csf_ex_unq_nb) = csf_i

!     list of excited determinants in excited csf
      call object_alloc ('det_unq_in_csf_ex_unq_nb', det_unq_in_csf_ex_unq_nb, csf_ex_unq_nb)
      call object_alloc ('det_unq_up_in_csf_ex_unq', det_unq_up_in_csf_ex_unq, csf_ex_unq_nb)
      call object_alloc ('det_unq_dn_in_csf_ex_unq', det_unq_dn_in_csf_ex_unq, csf_ex_unq_nb)
      call object_alloc ('cdet_unq_in_csf_ex_unq', cdet_unq_in_csf_ex_unq, csf_ex_unq_nb)
      det_unq_in_csf_ex_unq_nb (csf_ex_unq_nb) = det_unq_in_csf_ex_cur_nb
      call copy (det_unq_up_in_csf_ex_cur, det_unq_up_in_csf_ex_unq (csf_ex_unq_nb)%row)
      call copy (det_unq_dn_in_csf_ex_cur, det_unq_dn_in_csf_ex_unq (csf_ex_unq_nb)%row)
      call copy (cdet_unq_in_csf_ex_cur, cdet_unq_in_csf_ex_unq (csf_ex_unq_nb)%row)
!      write(6,'(2a)') trim(here),': add new excited csf'
     endif ! csf_cur == 0

!    release temporary arrays
     call release ('det_unq_up_in_csf_ex_cur', det_unq_up_in_csf_ex_cur)
     call release ('det_unq_dn_in_csf_ex_cur', det_unq_dn_in_csf_ex_cur)
     call release ('cdet_unq_in_csf_ex_cur', cdet_unq_in_csf_ex_cur)

!    add this csf to the current singly-excited wave function
     call append (csf_unq_in_wf_ex_cur, csf_cur)
     call append (csf_unq_ref_in_wf_ex_cur, csf_i)
     call append (csf_unq_prefac_in_wf_ex_cur, csf_ex_cur_prefac)
     csf_unq_in_wf_ex_cur_nb = mysize (csf_unq_in_wf_ex_cur)

     endif ! if (.not. csf_ex_unq_is_zero)

    enddo ! ncsf

!   current singly-excited wave function
    if (.not. ex_is_zero) then

!    sort excited csf in current singly-excited wave function
     call sort (csf_unq_in_wf_ex_cur (:), csf_unq_ref_in_wf_ex_cur (:), csf_unq_prefac_in_wf_ex_cur (:))
!     write(6,'(2a)') trim(here), ': construct singly-excited wave function:'
!     write(6,'(2a,100i4)') trim(here),': csf_unq_in_wf_ex_cur=',csf_unq_in_wf_ex_cur (:)
!     write(6,'(2a,100i4)') trim(here),': csf_unq_ref_in_wf_ex_cur=',csf_unq_ref_in_wf_ex_cur (:)
!     write(6,'(2a,100f7.3)') trim(here),': csf_unq_prefac_in_wf_ex_cur=',csf_unq_prefac_in_wf_ex_cur (:)

!    check if current singly-excited wave function is a singly-excited wave function already encountered
!    (do not check for equalities between csf coefficients)
     ex_cur = 0
     do ex_i = 1, single_ex_nb
      if (arrays_equal (csf_unq_in_wf_ex_cur (:), csf_unq_in_wf_ex (ex_i)%row (:)) .and.            &
          arrays_equal (csf_unq_ref_in_wf_ex_cur (:), csf_unq_ref_in_wf_ex (ex_i)%row (:)) .and.    &
          arrays_equal (csf_unq_prefac_in_wf_ex_cur (:), csf_unq_prefac_in_wf_ex (ex_i)%row (:))) then
         ex_cur = ex_i
 !        write(6,'(2a)') trim(here), ': singly-excited wave function is a singly-excited wave function already encountered'
      endif
     enddo ! csf_ex_unq_i

!    if current singly-excited wave function is a new singly-excited wave function, add it to the list of unique singly-excited wave functions
     if (ex_cur == 0) then
      single_ex_nb = single_ex_nb + 1
      single_ex_added_nb = single_ex_added_nb + 1
      ex_cur = single_ex_nb
      call object_alloc ('csf_unq_in_wf_ex_nb', csf_unq_in_wf_ex_nb, single_ex_nb)
      call object_alloc ('csf_unq_in_wf_ex', csf_unq_in_wf_ex, single_ex_nb)
      call object_alloc ('csf_unq_ref_in_wf_ex', csf_unq_ref_in_wf_ex, single_ex_nb)
      call object_alloc ('csf_unq_prefac_in_wf_ex', csf_unq_prefac_in_wf_ex, single_ex_nb)
      csf_unq_in_wf_ex_nb (single_ex_nb) = csf_unq_in_wf_ex_cur_nb
      call copy (csf_unq_in_wf_ex_cur, csf_unq_in_wf_ex (single_ex_nb)%row)
      call copy (csf_unq_ref_in_wf_ex_cur, csf_unq_ref_in_wf_ex (single_ex_nb)%row)
      call copy (csf_unq_prefac_in_wf_ex_cur, csf_unq_prefac_in_wf_ex (single_ex_nb)%row)

      call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
      call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
      ex_orb_1st_lab (single_ex_nb) = orb_opt_lab_i
      ex_orb_2nd_lab (single_ex_nb) = orb_opt_lab_j

 !     write(6,'(2a)') trim(here), ': add new singly-excited wave function'
     endif ! ex_cur == 0

!     build current orbital derivative: E(i->j) - E(j->i)
      dpsi_orb_is_zero = .false.

!     direct excitation i->j
      if (ex_dir_rev == 1) then

        call append (csf_unq_in_dpsi_orb_cur, csf_unq_in_wf_ex_cur)
        call append (csf_unq_ref_in_dpsi_orb_cur, csf_unq_ref_in_wf_ex_cur)
        call append (csf_unq_prefac_in_dpsi_orb_cur, csf_unq_prefac_in_wf_ex_cur)
        csf_unq_in_dpsi_orb_cur_nb = mysize (csf_unq_in_dpsi_orb_cur)

        ex_orb_ind_cur = ex_cur

 !    reverse excitation j->i
      else

!       opposite sign
        csf_unq_prefac_in_wf_ex_cur (:) = - csf_unq_prefac_in_wf_ex_cur (:)

        call append (csf_unq_in_dpsi_orb_cur, csf_unq_in_wf_ex_cur)
        call append (csf_unq_ref_in_dpsi_orb_cur, csf_unq_ref_in_wf_ex_cur)
        call append (csf_unq_prefac_in_dpsi_orb_cur, csf_unq_prefac_in_wf_ex_cur)
        csf_unq_in_dpsi_orb_cur_nb = mysize (csf_unq_in_dpsi_orb_cur)

        ex_orb_ind_rev_cur = ex_cur

      endif ! if (ex_dir_rev == 1)

!    release temporary arrays
     call release ('csf_unq_in_wf_ex_cur', csf_unq_in_wf_ex_cur)
     call release ('csf_unq_ref_in_wf_ex_cur', csf_unq_ref_in_wf_ex_cur)
     call release ('csf_unq_prefac_in_wf_ex_cur', csf_unq_prefac_in_wf_ex_cur)

    endif ! if (.not. ex_is_zero) then

    enddo ! ex_dir_rev


    if (.not. dpsi_orb_is_zero) then

    dpsi_orb_is_redundant = .false.

!   sort csf in current orbital derivative
    call sort (csf_unq_in_dpsi_orb_cur (:), csf_unq_ref_in_dpsi_orb_cur (:), csf_unq_prefac_in_dpsi_orb_cur (:))

!   combine unique csfs in current orbital derivative
    csf_unq_in_dpsi_orb_cur_temp_nb = csf_unq_in_dpsi_orb_cur_nb
    call copy (csf_unq_in_dpsi_orb_cur, csf_unq_in_dpsi_orb_cur_temp)
    call copy (csf_unq_ref_in_dpsi_orb_cur, csf_unq_ref_in_dpsi_orb_cur_temp)
    call copy (csf_unq_prefac_in_dpsi_orb_cur, csf_unq_prefac_in_dpsi_orb_cur_temp)
    call release ('csf_unq_in_dpsi_orb_cur', csf_unq_in_dpsi_orb_cur)
    call release ('csf_unq_ref_in_dpsi_orb_cur', csf_unq_ref_in_dpsi_orb_cur)
    call release ('csf_unq_prefac_in_dpsi_orb_cur', csf_unq_prefac_in_dpsi_orb_cur)

    csf_unq_in_dpsi_orb_cur_nb = 1
    call append (csf_unq_in_dpsi_orb_cur, csf_unq_in_dpsi_orb_cur_temp (1))
    call append (csf_unq_ref_in_dpsi_orb_cur, csf_unq_ref_in_dpsi_orb_cur_temp (1))
    call append (csf_unq_prefac_in_dpsi_orb_cur, csf_unq_prefac_in_dpsi_orb_cur_temp (1))

    do csf_unq_in_dpsi_orb_cur_i = 2, csf_unq_in_dpsi_orb_cur_temp_nb
      if (csf_unq_in_dpsi_orb_cur_temp (csf_unq_in_dpsi_orb_cur_i) == csf_unq_in_dpsi_orb_cur (csf_unq_in_dpsi_orb_cur_nb) .and.  &
          csf_unq_ref_in_dpsi_orb_cur_temp (csf_unq_in_dpsi_orb_cur_i) == csf_unq_ref_in_dpsi_orb_cur (csf_unq_in_dpsi_orb_cur_nb)) then
        csf_unq_prefac_in_dpsi_orb_cur (csf_unq_in_dpsi_orb_cur_nb) = csf_unq_prefac_in_dpsi_orb_cur (csf_unq_in_dpsi_orb_cur_nb)  &
                                                            + csf_unq_prefac_in_dpsi_orb_cur_temp (csf_unq_in_dpsi_orb_cur_i)
      else
       csf_unq_in_dpsi_orb_cur_nb = csf_unq_in_dpsi_orb_cur_nb + 1
       call append (csf_unq_in_dpsi_orb_cur, csf_unq_in_dpsi_orb_cur_temp (csf_unq_in_dpsi_orb_cur_i))
       call append (csf_unq_ref_in_dpsi_orb_cur, csf_unq_ref_in_dpsi_orb_cur_temp (csf_unq_in_dpsi_orb_cur_i))
       call append (csf_unq_prefac_in_dpsi_orb_cur, csf_unq_prefac_in_dpsi_orb_cur_temp (csf_unq_in_dpsi_orb_cur_i))
      endif
    enddo

!    write(6,'(2a)') trim(here), ': construct orbital derivative:'
!    write(6,'(2a,100i4)') trim(here),': csf_unq_in_dpsi_orb_cur=',csf_unq_in_dpsi_orb_cur (:)
!    write(6,'(2a,100i4)') trim(here),': csf_unq_ref_in_dpsi_orb_cur=',csf_unq_ref_in_dpsi_orb_cur (:)
!    write(6,'(2a,100f7.3)') trim(here),': csf_unq_prefac_in_dpsi_orb_cur=',csf_unq_prefac_in_dpsi_orb_cur (:)

!   remove csfs with zero prefactors
    csf_unq_in_dpsi_orb_cur_i = 1
    do while (csf_unq_in_dpsi_orb_cur_i <= csf_unq_in_dpsi_orb_cur_nb)
      if (csf_unq_prefac_in_dpsi_orb_cur (csf_unq_in_dpsi_orb_cur_i) == 0.d0) then
        call remove_elt_in_array (csf_unq_in_dpsi_orb_cur, csf_unq_in_dpsi_orb_cur_i)
        call remove_elt_in_array (csf_unq_ref_in_dpsi_orb_cur, csf_unq_in_dpsi_orb_cur_i)
        call remove_elt_in_array (csf_unq_prefac_in_dpsi_orb_cur, csf_unq_in_dpsi_orb_cur_i)
        csf_unq_in_dpsi_orb_cur_nb = mysize (csf_unq_in_dpsi_orb_cur)
      else
        csf_unq_in_dpsi_orb_cur_i = csf_unq_in_dpsi_orb_cur_i + 1
      endif
    enddo

!   if all csf prefactors of current orbital derivative are zero, this orbital derivative is redundant
    if (csf_unq_in_dpsi_orb_cur_nb == 0.d0) then
        dpsi_orb_is_redundant = .true.
!        write(6, '(2a)') trim(here),': orbital derivative is redundant (zero)'
    endif

!   if current orbital derivative is a linear combinaison of optimized ground-state csfs, then
!   this orbital derivative is redundant.
!   (Even if not optimized, the first csf is also considered as part of the optimization space.
!   This is necessary so that for a CASSCF wave function all active->active excitations are redundant.
!   I am sure if one should not include also the first csf is the optimization for this to be really correct,
!   as for instance danish MCSCF people do).
    if (.not. dpsi_orb_is_redundant .and. l_opt_csf) then
      dpsi_orb_is_redundant = .true.
      do csf_unq_in_dpsi_orb_cur_i = 1, csf_unq_in_dpsi_orb_cur_nb
        if (.not. (elt_in_array (iwcsf (1:nparmcsf), csf_unq_in_dpsi_orb_cur (csf_unq_in_dpsi_orb_cur_i)) .or. &
                  csf_unq_in_dpsi_orb_cur (csf_unq_in_dpsi_orb_cur_i) == 1)) then
            dpsi_orb_is_redundant = .false.
            exit
        endif
      enddo
      if (dpsi_orb_is_redundant) then
!        write(6, '(2a)') trim(here),': orbital derivative is redundant with linear combinaison of optimized ground-state csfs'
      endif
    endif

!   check if current orbital derivative is an orbital derivative already encountered
    if (.not. dpsi_orb_is_redundant) then
     do dpsi_orb_i = 1, param_orb_nb
      if (arrays_equal (csf_unq_in_dpsi_orb_cur (:), csf_unq_in_dpsi_orb (dpsi_orb_i)%row (:)) .and.            &
          (csf_unq_in_dpsi_orb_cur_nb == 1 .or.                                                           &
          (arrays_equal (csf_unq_ref_in_dpsi_orb_cur (:), csf_unq_ref_in_dpsi_orb (dpsi_orb_i)%row (:)) .and.   &
          arrays_equal (csf_unq_prefac_in_dpsi_orb_cur (:), csf_unq_prefac_in_dpsi_orb (dpsi_orb_i)%row (:))))) then
         dpsi_orb_is_redundant = .true.
!         write(6, '(2a)') trim(here),': orbital derivative is redundant with another orbital derivative'
      endif
     enddo ! dpsi_orb_i
    endif

!   if current orbital derivative is redundant, remove added unique singly-excited wave functions
    if (dpsi_orb_is_redundant) then
      single_ex_nb = single_ex_nb - single_ex_added_nb
      call object_alloc ('csf_unq_in_wf_ex_nb', csf_unq_in_wf_ex_nb, single_ex_nb)
      call object_alloc ('csf_unq_in_wf_ex', csf_unq_in_wf_ex, single_ex_nb)
      call object_alloc ('csf_unq_ref_in_wf_ex', csf_unq_ref_in_wf_ex, single_ex_nb)
      call object_alloc ('csf_unq_prefac_in_wf_ex', csf_unq_prefac_in_wf_ex, single_ex_nb)
      call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
      call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
!      write(6,'(2a)') trim(here), ': remove singly-excited wave function'
    endif

!   if current orbital derivative is a not redundant, add it to the list of orbital derivatives
    if (.not. dpsi_orb_is_redundant) then
      param_orb_nb = param_orb_nb + 1
      call object_alloc ('csf_unq_in_dpsi_orb_nb', csf_unq_in_dpsi_orb_nb, param_orb_nb)
      call object_alloc ('csf_unq_in_dpsi_orb', csf_unq_in_dpsi_orb, param_orb_nb)
      call object_alloc ('csf_unq_ref_in_dpsi_orb', csf_unq_ref_in_dpsi_orb, param_orb_nb)
      call object_alloc ('csf_unq_prefac_in_dpsi_orb', csf_unq_prefac_in_dpsi_orb, param_orb_nb)
      csf_unq_in_dpsi_orb_nb (param_orb_nb) = csf_unq_in_dpsi_orb_cur_nb
      call copy (csf_unq_in_dpsi_orb_cur, csf_unq_in_dpsi_orb(param_orb_nb)%row)
      call copy (csf_unq_ref_in_dpsi_orb_cur, csf_unq_ref_in_dpsi_orb(param_orb_nb)%row)
      call copy (csf_unq_prefac_in_dpsi_orb_cur, csf_unq_prefac_in_dpsi_orb(param_orb_nb)%row)

      call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
      call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
      if (ex_orb_ind_cur == 0) then
       write(6,'(2a)') trim(here),': orbital derivative with vanishing direct excitation'
       write(6,'(2a)') trim(here),': this case is not yet implemented'
       write(6,'(2a,i4)') trim(here),': the index of reverse excitation is ',ex_orb_ind_rev_cur
       call die (here)
      endif
      ex_orb_ind (param_orb_nb) = ex_orb_ind_cur
      ex_orb_ind_rev (param_orb_nb) = ex_orb_ind_rev_cur

!      write(6,'(2a)') trim(here),': add new orbital derivative'
    endif ! if (.not. dpsi_orb_is_redundant)

!   release temporary arrays
    call release ('csf_unq_in_dpsi_orb_cur', csf_unq_in_dpsi_orb_cur)
    call release ('csf_unq_ref_in_dpsi_orb_cur', csf_unq_ref_in_dpsi_orb_cur)
    call release ('csf_unq_prefac_in_dpsi_orb_cur', csf_unq_prefac_in_dpsi_orb_cur)

    endif ! if (.not. dpsi_orb_is_zero)

!   reinitializations
    single_ex_added_nb = 0
    ex_orb_ind_cur = 0
    ex_orb_ind_rev_cur = 0

!  reset label of first orbital
   orb_opt_lab_i = orb_opt_lab_j

   enddo ! orb_opt_j

  enddo ! orb_opt_i

  write(6,'(a,i10)') ' Number of single orbital excitations = ', single_ex_nb
  write(6,'(a,i10)') ' Number of orbital derivatives        = ', param_orb_nb

! excitation pairs number
  deriv_orb_pairs_nb = param_orb_nb * (param_orb_nb + 1) / 2

! check with object from old routine
!  write(6,'(2a,i10)') trim(here),': old number of single excitations for orbital optimization: ', single_ex_nb_old
!  write(6,'(2a,i10)') trim(here),': old number of orbital derivatives: ', param_orb_nb_old
!
!  do ex_i = 1, single_ex_nb_old
!   orb_1st = ex_orb_1st_lab_old (ex_i)
!   orb_2nd = ex_orb_2nd_lab_old (ex_i)
!   ex_found = .false.
!   do ex_j = 1, single_ex_nb
!    if (orb_1st == ex_orb_1st_lab (ex_j) .and. orb_2nd == ex_orb_2nd_lab (ex_j)) then
!      ex_found = .true.
!     exit
!    endif
!   enddo
!   if (.not. ex_found) then
!    write(6,'(2a,i4,a,i4,a)') trim(here),': excitation ',orb_1st, ' ->', orb_2nd, ' is no longer present'
!   endif
!  enddo

 end subroutine single_ex_wf_bld_2

! ==============================================================================
  subroutine single_ex_det_bld
! ------------------------------------------------------------------------------
! Description   : build single orbital excitation information for determinants
! Description   : to be replaced by single_ex_wf_bld_2
!
! Created       : J. Toulouse, 25 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer ex_i, orb_1st, orb_2nd
  integer det_unq_up_i, det_unq_up_k, det_unq_dn_i, det_unq_dn_k

! header
  if (header_exe) then

   call object_create ('det_ex_unq_up_nb')
   call object_create ('det_ex_unq_dn_nb')
   call object_create ('det_ex_unq_up_orb_1st_pos')
   call object_create ('det_ex_unq_dn_orb_1st_pos')
   call object_create ('det_ex_unq_up_orb_2nd_lab')
   call object_create ('det_ex_unq_dn_orb_2nd_lab')
   call object_create ('det_ex_orb_lab_srt_up')
   call object_create ('det_ex_orb_lab_srt_dn')
   call object_create ('det_ex_orb_lab_srt_sgn_up')
   call object_create ('det_ex_orb_lab_srt_sgn_dn')
   call object_create ('det_ex_unq_orb_lab_up')
   call object_create ('det_ex_unq_orb_lab_dn')
   call object_create ('det_ex_unq_orb_lab_srt_up')
   call object_create ('det_ex_unq_orb_lab_srt_dn')
   call object_create ('det_ex_unq_orb_lab_srt_sgn_up')
   call object_create ('det_ex_unq_orb_lab_srt_sgn_dn')
   call object_create ('ex_up_is_zero')
   call object_create ('ex_dn_is_zero')
   call object_create ('is_det_ex_up')
   call object_create ('is_det_ex_dn')
   call object_create ('iwdet_ex_up')
   call object_create ('iwdet_ex_dn')
   call object_create ('iwdet_ex_ref_up')
   call object_create ('iwdet_ex_ref_dn')

   call object_needed ('single_ex_nb')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('ex_orb_1st_lab')
   call object_needed ('ex_orb_2nd_lab')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('orb_occ_in_det_unq_up')
   call object_needed ('orb_occ_in_det_unq_dn')
   call object_needed ('orb_pos_in_det_unq_up')
   call object_needed ('orb_pos_in_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('ex_up_is_zero', ex_up_is_zero, single_ex_nb, ndetup)
  call object_alloc ('ex_dn_is_zero', ex_dn_is_zero, single_ex_nb, ndetdn)
  call object_alloc ('det_ex_orb_lab_up', det_ex_orb_lab_up, nup, single_ex_nb, ndetup)
  call object_alloc ('det_ex_orb_lab_dn', det_ex_orb_lab_dn, ndn, single_ex_nb, ndetdn)
  call object_alloc ('det_ex_orb_lab_srt_up', det_ex_orb_lab_srt_up, nup, single_ex_nb, ndetup)
  call object_alloc ('det_ex_orb_lab_srt_dn', det_ex_orb_lab_srt_dn, ndn, single_ex_nb, ndetdn)
  call object_alloc ('det_ex_orb_lab_srt_sgn_up', det_ex_orb_lab_srt_sgn_up, single_ex_nb, ndetup)
  call object_alloc ('det_ex_orb_lab_srt_sgn_dn', det_ex_orb_lab_srt_sgn_dn, single_ex_nb, ndetdn)
  call object_alloc ('is_det_ex_up', is_det_ex_up, single_ex_nb, ndetup)
  call object_alloc ('is_det_ex_dn', is_det_ex_dn, single_ex_nb, ndetdn)
  call object_alloc ('iwdet_ex_up', iwdet_ex_up, single_ex_nb, ndetup)
  call object_alloc ('iwdet_ex_dn', iwdet_ex_dn, single_ex_nb, ndetdn)


  ex_up_is_zero = .true.
  ex_dn_is_zero = .true.

  det_ex_orb_lab_up = 0
  det_ex_orb_lab_dn = 0
  det_ex_orb_lab_srt_up = 0
  det_ex_orb_lab_srt_dn = 0

  det_ex_unq_up_nb = 0
  det_ex_unq_dn_nb = 0
  is_det_ex_up = .true.
  is_det_ex_dn = .true.
  iwdet_ex_dn = 0
  iwdet_ex_up = 0

! loop over single orbital excitations
  do ex_i = 1, single_ex_nb

     orb_1st = ex_orb_1st_lab (ex_i)
     orb_2nd = ex_orb_2nd_lab (ex_i)

!    loop over unique spin-up determinants
     do det_unq_up_i = 1, ndetup

     if (orb_occ_in_det_unq_up (orb_1st, det_unq_up_i) .and. .not. orb_occ_in_det_unq_up (orb_2nd, det_unq_up_i)) then
       ex_up_is_zero (ex_i, det_unq_up_i) = .false.

!      build current excited determinant
       det_ex_orb_lab_up (:, ex_i, det_unq_up_i) = det_unq_orb_lab_srt_up (:, det_unq_up_i)
       call replace_elt_in_array (det_ex_orb_lab_up (:, ex_i, det_unq_up_i), orb_1st, orb_2nd)
       det_ex_orb_lab_srt_up (:, ex_i, det_unq_up_i) = det_ex_orb_lab_up (:, ex_i, det_unq_up_i)
       call sort_and_sign (det_ex_orb_lab_srt_up (:, ex_i, det_unq_up_i), det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i))

!      check if current excited determinant is a determinant in the ground-state wave function
       do det_unq_up_k = 1, ndetup
        if (arrays_equal (det_ex_orb_lab_srt_up (:, ex_i, det_unq_up_i), det_unq_orb_lab_srt_up (:, det_unq_up_k))) then
          is_det_ex_up (ex_i, det_unq_up_i) = .false.
          iwdet_ex_up (ex_i, det_unq_up_i) = det_unq_up_k
          exit
        endif
       enddo

!      check if current excited determinant is a excited determinant already encountered
       do det_unq_up_k = 1, det_ex_unq_up_nb
        if (arrays_equal (det_ex_orb_lab_srt_up (:, ex_i, det_unq_up_i), det_ex_unq_orb_lab_srt_up (:, det_unq_up_k))) then
          iwdet_ex_up (ex_i, det_unq_up_i) =  det_unq_up_k
          exit
        endif
       enddo

!      if current excited determinant is a new determinant, add it to the list of unique excited determinants
       if (iwdet_ex_up (ex_i, det_unq_up_i) == 0) then
         det_ex_unq_up_nb = det_ex_unq_up_nb + 1
         iwdet_ex_up (ex_i, det_unq_up_i) = det_ex_unq_up_nb
         call object_alloc ('iwdet_ex_ref_up', iwdet_ex_ref_up, det_ex_unq_up_nb)
         call object_alloc ('det_ex_unq_orb_lab_up', det_ex_unq_orb_lab_up, nup, det_ex_unq_up_nb)
         call object_alloc ('det_ex_unq_orb_lab_srt_up', det_ex_unq_orb_lab_srt_up, nup, det_ex_unq_up_nb)
         call object_alloc ('det_ex_unq_orb_lab_srt_sgn_up', det_ex_unq_orb_lab_srt_sgn_up, det_ex_unq_up_nb)
         call object_alloc ('det_ex_unq_up_orb_1st_pos', det_ex_unq_up_orb_1st_pos, det_ex_unq_up_nb)
         call object_alloc ('det_ex_unq_up_orb_2nd_lab', det_ex_unq_up_orb_2nd_lab, det_ex_unq_up_nb)
         iwdet_ex_ref_up (det_ex_unq_up_nb) = det_unq_up_i
         det_ex_unq_orb_lab_up (:, det_ex_unq_up_nb ) = det_ex_orb_lab_up (:, ex_i, det_unq_up_i)
         det_ex_unq_orb_lab_srt_up (:, det_ex_unq_up_nb ) = det_ex_orb_lab_srt_up (:, ex_i, det_unq_up_i)
         det_ex_unq_orb_lab_srt_sgn_up (det_ex_unq_up_nb) = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i)
         det_ex_unq_up_orb_1st_pos (det_ex_unq_up_nb) = orb_pos_in_det_unq_up (orb_1st, det_unq_up_i)
         det_ex_unq_up_orb_2nd_lab (det_ex_unq_up_nb) = orb_2nd
       endif

     endif

     enddo ! det_unq_up_i

!    loop over unique spin-down determinants
     do det_unq_dn_i = 1, ndetdn

     if (orb_occ_in_det_unq_dn (orb_1st, det_unq_dn_i) .and. .not. orb_occ_in_det_unq_dn (orb_2nd, det_unq_dn_i)) then
       ex_dn_is_zero (ex_i, det_unq_dn_i) = .false.

!      build current excited determinant
       det_ex_orb_lab_dn (:, ex_i, det_unq_dn_i) = det_unq_orb_lab_srt_dn (:, det_unq_dn_i)
       call replace_elt_in_array (det_ex_orb_lab_dn (:, ex_i, det_unq_dn_i), orb_1st, orb_2nd)
       det_ex_orb_lab_srt_dn (:, ex_i, det_unq_dn_i) = det_ex_orb_lab_dn (:, ex_i, det_unq_dn_i)
       call sort_and_sign (det_ex_orb_lab_srt_dn (:, ex_i, det_unq_dn_i), det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i))

!      check if current excited determinant is a determinant in the ground-state wave function
       do det_unq_dn_k = 1, ndetdn
        if (arrays_equal (det_ex_orb_lab_srt_dn (:, ex_i, det_unq_dn_i), det_unq_orb_lab_srt_dn (:, det_unq_dn_k))) then
          is_det_ex_dn (ex_i, det_unq_dn_i) = .false.
          iwdet_ex_dn (ex_i, det_unq_dn_i) = det_unq_dn_k
          exit
        endif
       enddo

!      check if current excited determinant is a excited determinant already encountered
       do det_unq_dn_k = 1, det_ex_unq_dn_nb
        if (arrays_equal (det_ex_orb_lab_srt_dn (:, ex_i, det_unq_dn_i), det_ex_unq_orb_lab_srt_dn (:, det_unq_dn_k))) then
          iwdet_ex_dn (ex_i, det_unq_dn_i) =  det_unq_dn_k
          exit
        endif
       enddo


!      if current excited determinant is a new determinant, add it to the list of unique excited determinants
       if (iwdet_ex_dn (ex_i, det_unq_dn_i) == 0) then
         det_ex_unq_dn_nb = det_ex_unq_dn_nb + 1
         iwdet_ex_dn (ex_i, det_unq_dn_i) = det_ex_unq_dn_nb
         call object_alloc ('iwdet_ex_ref_dn', iwdet_ex_ref_dn, det_ex_unq_dn_nb)
         call object_alloc ('det_ex_unq_orb_lab_dn', det_ex_unq_orb_lab_dn, ndn, det_ex_unq_dn_nb)
         call object_alloc ('det_ex_unq_orb_lab_srt_dn', det_ex_unq_orb_lab_srt_dn, ndn, det_ex_unq_dn_nb)
         call object_alloc ('det_ex_unq_orb_lab_srt_sgn_dn', det_ex_unq_orb_lab_srt_sgn_dn, det_ex_unq_dn_nb)
         call object_alloc ('det_ex_unq_dn_orb_1st_pos', det_ex_unq_dn_orb_1st_pos, det_ex_unq_dn_nb)
         call object_alloc ('det_ex_unq_dn_orb_2nd_lab', det_ex_unq_dn_orb_2nd_lab, det_ex_unq_dn_nb)
         iwdet_ex_ref_dn (det_ex_unq_dn_nb) = det_unq_dn_i
         det_ex_unq_orb_lab_dn (:, det_ex_unq_dn_nb) = det_ex_orb_lab_dn (:, ex_i, det_unq_dn_i)
         det_ex_unq_orb_lab_srt_dn (:, det_ex_unq_dn_nb) = det_ex_orb_lab_srt_dn (:, ex_i, det_unq_dn_i)
         det_ex_unq_orb_lab_srt_sgn_dn (det_ex_unq_dn_nb) = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i)
         det_ex_unq_dn_orb_1st_pos (det_ex_unq_dn_nb) = orb_pos_in_det_unq_dn (orb_1st, det_unq_dn_i)
         det_ex_unq_dn_orb_2nd_lab (det_ex_unq_dn_nb) = orb_2nd
       endif

     endif

     enddo ! det_unq_dn_i

  enddo ! ex_i

  write(6,'(a,i10)') ' Number of unique spin-up   excited determinants = ',det_ex_unq_up_nb
  write(6,'(a,i10)') ' Number of unique spin-down excited determinants = ',det_ex_unq_dn_nb

 end subroutine single_ex_det_bld

! ==============================================================================
  subroutine det_ex_unq_bld
! ------------------------------------------------------------------------------
! Description   :  Calculate value of unque singly-excited spin-up and spin-down determinants
! Description   :  by updating reference determinant via the Sherman-Morrison formula
!
! Created       : J. Toulouse, 27 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_ex_unq_up_i, det_ex_unq_dn_i
  integer det_unq_ref_up, det_unq_ref_dn
  integer i
  integer col_up, orb_up
  integer col_dn, orb_dn
  real(dp) factor_up, factor_dn

! header
  if (header_exe) then

   call object_create ('det_ex_unq_up', det_ex_unq_up_index)
   call object_create ('det_ex_unq_dn', det_ex_unq_dn_index)

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('det_ex_unq_up_nb')
   call object_needed ('det_ex_unq_dn_nb')
   call object_needed ('ex_up_is_zero')
   call object_needed ('ex_dn_is_zero')
   call object_needed ('iwdet_ex_ref_up')
   call object_needed ('iwdet_ex_ref_dn')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('det_ex_unq_up_orb_1st_pos')
   call object_needed ('det_ex_unq_dn_orb_1st_pos')
   call object_needed ('det_ex_unq_up_orb_2nd_lab')
   call object_needed ('det_ex_unq_dn_orb_2nd_lab')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('orb')
   call object_needed ('detu')
   call object_needed ('detd')

   return

  endif

! begin

! allocations
  call object_alloc ('det_ex_unq_up', det_ex_unq_up, det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_dn', det_ex_unq_dn, det_ex_unq_dn_nb)

! spin up determinants
  do det_ex_unq_up_i = 1, det_ex_unq_up_nb

      col_up  = det_ex_unq_up_orb_1st_pos (det_ex_unq_up_i)
      orb_up  = det_ex_unq_up_orb_2nd_lab (det_ex_unq_up_i)
      det_unq_ref_up = iwdet_ex_ref_up (det_ex_unq_up_i)

      factor_up = 0.d0
      do i = 1, nup
       factor_up = factor_up + slater_mat_trans_inv_up (i, col_up, det_unq_ref_up) * orb (i, orb_up)
      enddo

      det_ex_unq_up (det_ex_unq_up_i) = factor_up * detu (det_unq_ref_up)

   enddo ! det_ex_unq_up_i

! spin down determinants
  do det_ex_unq_dn_i = 1, det_ex_unq_dn_nb

      col_dn  = det_ex_unq_dn_orb_1st_pos (det_ex_unq_dn_i)
      orb_dn  = det_ex_unq_dn_orb_2nd_lab (det_ex_unq_dn_i)
      det_unq_ref_dn = iwdet_ex_ref_dn (det_ex_unq_dn_i)

      factor_dn = 0.d0
      do i = 1, ndn
       factor_dn = factor_dn + slater_mat_trans_inv_dn (i, col_dn, det_unq_ref_dn) * orb (nup + i, orb_dn)
      enddo

      det_ex_unq_dn (det_ex_unq_dn_i) = factor_dn * detd (det_unq_ref_dn)

   enddo ! det_ex_unq_dn_i

  end subroutine det_ex_unq_bld

! ==============================================================================
  subroutine det_ex_bld
! ------------------------------------------------------------------------------
! Description   : all excited determinants
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer ex_i, ex_rev_i
  integer csf_i, det_in_csf_i, det_i, det_unq_up_i, det_unq_dn_i
  integer iwdet, sgn

! header
  if (header_exe) then

   call object_create ('det_ex_up')
   call object_create ('det_ex_dn')
   call object_create ('det_ex')

   call object_needed ('single_ex_nb')
   call object_needed ('ncsf')
   call object_needed ('ndet')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('is_det_ex_up')
   call object_needed ('is_det_ex_dn')
   call object_needed ('iwdet_ex_up')
   call object_needed ('iwdet_ex_dn')
   call object_needed ('det_ex_unq_up')
   call object_needed ('det_ex_unq_dn')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('det_ex_orb_lab_srt_sgn_up')
   call object_needed ('det_ex_orb_lab_srt_sgn_dn')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('det_ex_up', det_ex_up, single_ex_nb, ndet)
  call object_alloc ('det_ex_dn', det_ex_dn, single_ex_nb, ndet)
  call object_alloc ('det_ex', det_ex, single_ex_nb, ndet)

! loop over single orbital excitations
  do ex_i = 1, single_ex_nb

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

        det_i = iwdet_in_csf (det_in_csf_i, csf_i)
        det_unq_up_i = det_to_det_unq_up (det_i)
        det_unq_dn_i = det_to_det_unq_dn (det_i)

!       spin-up excited determinants:

!       current excited determinant is zero
        if (ex_up_is_zero (ex_i, det_unq_up_i)) then
         det_ex_up (ex_i, det_i) = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif (.not. is_det_ex_up (ex_i, det_unq_up_i)) then
         iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i)
         det_ex_up (ex_i, det_i) = sgn * detu (iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i) * det_ex_unq_orb_lab_srt_sgn_up (iwdet)
         det_ex_up (ex_i, det_i) = sgn * det_ex_unq_up (iwdet)
        endif

!       spin-down excited determinants:

!       current excited determinant is zero

        if (ex_dn_is_zero (ex_i, det_unq_dn_i)) then
         det_ex_dn (ex_i, det_i) = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif ( .not. is_det_ex_dn (ex_i, det_unq_dn_i) ) then
         iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i)
         det_ex_dn (ex_i, det_i) = sgn * detu (iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i) * det_ex_unq_orb_lab_srt_sgn_dn (iwdet)
         det_ex_dn (ex_i, det_i) = sgn * det_ex_unq_dn (iwdet)
        endif

        det_ex (ex_i, det_i) = det_ex_up (ex_i, det_i) * detd (det_unq_dn_i) + detu (det_unq_up_i) * det_ex_dn (ex_i, det_i)

     enddo ! det_in_csf_i
   enddo ! csf_i

  enddo ! ex_i

  end subroutine det_ex_bld

! ==============================================================================
  subroutine dpsi_orb_bld
! ------------------------------------------------------------------------------
! Description   : Logarithm derivatives of Psi with respect to orbital rotations kij
! Description   : d ln Psi / d kij = (1/Psi) * d Psi / d kij
! Description   :                  = (1/Psi) * Psi (i->j)
! Description   : where Psi is the reference wave function
! Description   : and Psi (i->j) = E(i->j) Psi is a singly excited wave function
!
! Created       : J. Toulouse, 12 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer ex_i, ex_rev_i
  integer csf_i, det_in_csf_i, det_i, dorb_i
  real(dp) factor_up, factor_dn
  real(dp) det

! header
  if (header_exe) then

   call object_create ('dpsi_orb')
   call object_create ('psid_ex')

   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_ind_rev')
   call object_needed ('ncsf')
   call object_needed ('ndet')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('det_ex')
   call object_needed ('psido')

   return

  endif

! begin

! allocations
  call object_alloc ('psid_ex', psid_ex, param_orb_nb)
  call object_alloc ('dpsi_orb', dpsi_orb, param_orb_nb)

  psid_ex = 0.d0

  do dorb_i = 1, param_orb_nb

   ex_i = ex_orb_ind (dorb_i)
   ex_rev_i = ex_orb_ind_rev (dorb_i)

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

        det_i = iwdet_in_csf (det_in_csf_i, csf_i)

        det = det_ex (ex_i, det_i)

!       reverse excitation for non casscf wave function
        if (.not. l_casscf .and. ex_rev_i /= 0) then
!          (Eij-Eji) Det
           det = det - det_ex (ex_rev_i, det_i)
        endif

        psid_ex (dorb_i) = psid_ex (dorb_i) + csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i) * det

     enddo ! det_in_csf_i
   enddo ! csf_i

    dpsi_orb (dorb_i) = psid_ex (dorb_i) / psido

 enddo ! ex_i


!    do ex_i = 1, single_ex_nb
!     write(6,'(a,a,i,a,f)') trim(here),': ex_i=',ex_i,' dpsi_orb=', dpsi_orb (ex_i)
!    enddo

! test for helium
!     write(6,*) trim(here),': ---------------------------------------------'
!     write(6,*) trim(here),': test for helium:'
!     write(6,*) trim(here),': factor_up=',orb(1,2)/orb(1,1)
!     write(6,*) trim(here),': factor_dn=',orb(2,2)/orb(2,1)
!     write(6,*) trim(here),': ---------------------------------------------'
!    stop

! test for lithium
!     write(6,*) trim(here),': ---------------------------------------------'
!     write(6,*) trim(here),': test for lithium:'
!     psi_0 = (orb(1,1)*orb(2,2)-orb(1,2)*orb(2,1))*orb(3,1)
!     psi_12_dn = (orb(1,1)*orb(2,2)-orb(1,2)*orb(2,1))*orb(3,2)
!     psi_13_up = (orb(1,3)*orb(2,2)-orb(1,2)*orb(2,3))*orb(3,1)
!     psi_13_dn = (orb(1,1)*orb(2,2)-orb(1,2)*orb(2,1))*orb(3,3)
!     dpsi_orb (1) = psi_12_dn/psi_0
!     dpsi_orb (2) = (psi_13_up + psi_13_dn)/psi_0
!
!    do ex_i = 1, single_ex_nb
!    write(6,'(a,a,i,a,f)') trim(here),': ex_i=',ex_i,' dpsi_orb=', dpsi_orb (ex_i)
!    enddo
!     write(6,*) trim(here),': ---------------------------------------------'
!    stop

!  write(6,*) trim(here),': exiting'
!  call object_provide ('dpsi_orb_test')
!  dpsi_orb = dpsi_orb_test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  call is_equal_or_die (dpsi_orb, dpsi_orb_test, 1.d-8)

  end subroutine dpsi_orb_bld

! ==============================================================================
  subroutine slater_mat_ex_trans_inv_bld
! ------------------------------------------------------------------------------
! Description   : inverse of transpose of Slater matrix corresponding to
! Description   : to excited determinants by inversion
! Description   : suprisingly, this is faster than using the Sherman-Morison formula!
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_i, orb_i, elec_i
  integer i, j, k, l

! header
  if (header_exe) then

   call object_create ('slater_mat_ex_trans_inv_up')
   call object_create ('slater_mat_ex_trans_inv_dn')
!   call object_create ('det_ex_unq_up')
!   call object_create ('det_ex_unq_dn')

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('det_ex_unq_up_nb')
   call object_needed ('det_ex_unq_dn_nb')
   call object_needed ('orb')

   return

  endif

! begin

! allocations
  call object_alloc ('slater_mat_ex_trans_up', slater_mat_ex_trans_up, nup, nup, det_ex_unq_up_nb)
  call object_alloc ('slater_mat_ex_trans_dn', slater_mat_ex_trans_dn, ndn, ndn, det_ex_unq_dn_nb)
  call object_alloc ('slater_mat_ex_trans_inv_up', slater_mat_ex_trans_inv_up, nup, nup, det_ex_unq_up_nb)
  call object_alloc ('slater_mat_ex_trans_inv_dn', slater_mat_ex_trans_inv_dn, ndn, ndn, det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_up', det_ex_unq_up, det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_dn', det_ex_unq_dn, det_ex_unq_dn_nb)

! spin-up
  do det_i = 1, det_ex_unq_up_nb
    do orb_i = 1, nup
      do elec_i = 1, nup
        slater_mat_ex_trans_up (orb_i, elec_i, det_i) = orb (elec_i, det_ex_unq_orb_lab_up (orb_i, det_i))
      enddo
    enddo
  enddo

! spin-dn
  do det_i = 1, det_ex_unq_dn_nb
    do orb_i = 1, ndn
      do elec_i = 1, ndn
        slater_mat_ex_trans_dn (orb_i, elec_i, det_i) = orb (nup + elec_i, det_ex_unq_orb_lab_dn (orb_i, det_i))
      enddo
    enddo
  enddo


!!! WAS
!!$  call object_provide ('slater_mat_ex_trans_inv_up_2')
!!$  call object_provide ('slater_mat_ex_trans_inv_dn_2')
!!$  slater_mat_ex_trans_up= slater_mat_ex_trans_inv_up_2
!!$  slater_mat_ex_trans_dn= slater_mat_ex_trans_inv_dn_2
!!$  return
!!!!



  call object_alloc ('mat_flat_up', mat_flat_up, nup*nup, det_ex_unq_up_nb)
  call object_alloc ('mat_flat_dn', mat_flat_dn, ndn*ndn, det_ex_unq_dn_nb)



  do det_i = 1, det_ex_unq_up_nb
   call flatten (mat_flat_up (:,det_i), slater_mat_ex_trans_up (:,:,det_i), nup, nup)
!debug WAS
!!$   write(*,*) "slater_mat_ex_trans_up (:,:,det_i)", slater_mat_ex_trans_up (:,:,det_i)
!!$   write(*,*) "doing deti ", det_i, det_ex_unq_up_nb
!!$   write(*,*) "mat_falt_up",  mat_flat_up (:, det_i)
!!$   write(*,*)  "det_ex_unq_up_nb",  det_ex_unq_up_nb
!
   call matinv (mat_flat_up (:,det_i), nup, det_ex_unq_up (det_i))
   call unflatten (mat_flat_up (:,det_i), slater_mat_ex_trans_inv_up (:,:,det_i), nup, nup)
  enddo

  do det_i = 1, det_ex_unq_dn_nb
   call flatten (mat_flat_dn (:,det_i), slater_mat_ex_trans_dn (:,:,det_i), ndn, ndn)
   call matinv (mat_flat_dn (:,det_i), ndn, det_ex_unq_dn (det_i))
   call unflatten (mat_flat_dn (:,det_i), slater_mat_ex_trans_inv_dn (:,:,det_i), ndn, ndn)
  enddo

   call object_modified_by_index (det_ex_unq_up_index)
   call object_modified_by_index (det_ex_unq_dn_index)

! tests!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   call object_provide ('slater_mat_ex_trans_inv_up_2')
!   call object_provide ('slater_mat_ex_trans_inv_dn_2')
!
!   call object_write_2 ('slater_mat_ex_trans_inv_up', 'slater_mat_ex_trans_inv_up_2')
!   call object_write_2 ('slater_mat_ex_trans_inv_dn', 'slater_mat_ex_trans_inv_dn_2')
!   call is_equal_or_die (slater_mat_ex_trans_inv_up, slater_mat_ex_trans_inv_up_2, 1.d-5)
!   call is_equal_or_die (slater_mat_ex_trans_inv_dn, slater_mat_ex_trans_inv_dn_2, 1.d-5)
!
!   call object_provide ('slater_mat_ex_trans_inv_up_3')
!   call object_provide ('slater_mat_ex_trans_inv_dn_3')
!
!   call object_write_2 ('slater_mat_ex_trans_inv_up', 'slater_mat_ex_trans_inv_up_3')
!   call object_write_2 ('slater_mat_ex_trans_inv_dn', 'slater_mat_ex_trans_inv_dn_3')
!   call is_equal_or_die (slater_mat_ex_trans_inv_up, slater_mat_ex_trans_inv_up_3, 1.d-5)
!   call is_equal_or_die (slater_mat_ex_trans_inv_dn, slater_mat_ex_trans_inv_dn_3, 1.d-5)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine slater_mat_ex_trans_inv_bld

! ==============================================================================
  subroutine slater_mat_ex_trans_inv_2_bld
! ------------------------------------------------------------------------------
! Description   :  Inverse of transpose of Slater matrix corresponding to
! Description   :  to excited determinants using the Sherman-Morison formula
!
! Created       : J. Toulouse, 27 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_up_i, det_dn_i
  integer det_ref_up, det_ref_dn
  integer i, j, k, l
  integer col_up, orb_up
  integer col_dn, orb_dn
  real(dp) factor_up, factor_dn
  real(dp) sum_up, sum_dn

! header
  if (header_exe) then

   call object_create ('slater_mat_ex_trans_inv_up_2')
   call object_create ('slater_mat_ex_trans_inv_dn_2')

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('det_ex_unq_up_nb')
   call object_needed ('det_ex_unq_dn_nb')
   call object_needed ('iwdet_ex_ref_up')
   call object_needed ('iwdet_ex_ref_dn')
   call object_needed ('det_ex_unq_up_orb_1st_pos')
   call object_needed ('det_ex_unq_dn_orb_1st_pos')
   call object_needed ('det_ex_unq_up_orb_2nd_lab')
   call object_needed ('det_ex_unq_dn_orb_2nd_lab')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('orb')

   return

  endif

! begin

! allocations
  call object_alloc ('slater_mat_ex_trans_inv_up_2', slater_mat_ex_trans_inv_up_2, nup, nup, det_ex_unq_up_nb)
  call object_alloc ('slater_mat_ex_trans_inv_dn_2', slater_mat_ex_trans_inv_dn_2, ndn, ndn, det_ex_unq_dn_nb)

! spin up
  do det_up_i = 1, det_ex_unq_up_nb

      col_up  = det_ex_unq_up_orb_1st_pos (det_up_i)
      orb_up  = det_ex_unq_up_orb_2nd_lab (det_up_i)
      det_ref_up = iwdet_ex_ref_up (det_up_i)

      factor_up = 0.d0
      do i = 1, nup
       factor_up = factor_up + slater_mat_trans_inv_up (i, col_up, det_ref_up) * orb (i, orb_up)
      enddo

      i = col_up

      do k = 1, nup
        do j = 1 , nup

         sum_up = 0.d0
         do l = 1, nup
          sum_up = sum_up + slater_mat_trans_inv_up (l, j, det_ref_up) * orb (l, orb_up)
         enddo

         if (i == j) then
          sum_up = sum_up - 1.d0
         endif

         slater_mat_ex_trans_inv_up_2 (k, j, det_up_i) = slater_mat_trans_inv_up (k, j, det_ref_up) &
                                 - slater_mat_trans_inv_up (k, i, det_ref_up) * sum_up/factor_up

        enddo ! j
      enddo ! k

   enddo ! det_up_i

! spin down
  do det_dn_i = 1, det_ex_unq_dn_nb

      col_dn  = det_ex_unq_dn_orb_1st_pos (det_dn_i)
      orb_dn  = det_ex_unq_dn_orb_2nd_lab (det_dn_i)
      det_ref_dn = iwdet_ex_ref_dn (det_dn_i)

      factor_dn = 0.d0
      do i = 1, ndn
       factor_dn = factor_dn + slater_mat_trans_inv_dn (i, col_dn, det_ref_dn) * orb (nup + i, orb_dn)
      enddo

      i = col_dn

      do k = 1, ndn
        do j = 1 , ndn

         sum_dn = 0.d0
         do l = 1, ndn
          sum_dn = sum_dn + slater_mat_trans_inv_dn (l, j, det_ref_dn) * orb (nup + l, orb_dn)
         enddo

         if (i == j) then
          sum_dn = sum_dn - 1.d0
         endif

         slater_mat_ex_trans_inv_dn_2 (k, j, det_dn_i) = slater_mat_trans_inv_dn (k, j, det_ref_dn) &
                                 - slater_mat_trans_inv_dn (k, i, det_ref_dn) * sum_dn/factor_dn

        enddo ! j
      enddo ! k

   enddo ! det_dn_i


!  Post-conditions
!   call object_provide ('slater_mat_ex_trans_inv_up')
!   call object_provide ('slater_mat_ex_trans_inv_dn')

!   call object_write_2 ('slater_mat_ex_trans_inv_up', 'slater_mat_ex_trans_inv_up_2')
!   call object_write_2 ('slater_mat_ex_trans_inv_dn', 'slater_mat_ex_trans_inv_dn_2')
!!!   call is_equal_or_die (slater_mat_ex_trans_inv_up, slater_mat_ex_trans_inv_up_2, 1.d-8)
!!!   call is_equal_or_die (slater_mat_ex_trans_inv_dn, slater_mat_ex_trans_inv_dn_2, 1.d-8)

  end subroutine slater_mat_ex_trans_inv_2_bld

! ==============================================================================
  subroutine slater_mat_ex_trans_inv_3_bld
! ------------------------------------------------------------------------------
! Description   :  inverse of transpose of Slater matrix corresponding to
! Description   :  to excited determinants using the Sherman-Morison formula
!
! Created       : J. Toulouse, 27 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_up_i, det_dn_i
  integer det_ref_up, det_ref_dn
  integer i, j, k, l
  integer col_up, orb_up
  integer col_dn, orb_dn
  real(dp) factor_up_inv, factor_dn_inv
  real(dp), allocatable :: ratio_up (:), ratio_dn (:)

! header
  if (header_exe) then

!   call object_create ('slater_mat_ex_trans_inv_up')
!   call object_create ('slater_mat_ex_trans_inv_dn')

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('det_ex_unq_up_nb')
   call object_needed ('det_ex_unq_dn_nb')
   call object_needed ('iwdet_ex_ref_up')
   call object_needed ('iwdet_ex_ref_dn')
   call object_needed ('det_ex_unq_up_orb_1st_pos')
   call object_needed ('det_ex_unq_dn_orb_1st_pos')
   call object_needed ('det_ex_unq_up_orb_2nd_lab')
   call object_needed ('det_ex_unq_dn_orb_2nd_lab')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('orb')

   return

  endif

! begin

! allocations
  call object_alloc ('slater_mat_ex_trans_inv_up', slater_mat_ex_trans_inv_up, nup, nup, det_ex_unq_up_nb)
  call object_alloc ('slater_mat_ex_trans_inv_dn', slater_mat_ex_trans_inv_dn, ndn, ndn, det_ex_unq_dn_nb)
  call alloc ('ratio_up', ratio_up, nup)
  call alloc ('ratio_dn', ratio_dn, ndn)

! spin up
  do det_up_i = 1, det_ex_unq_up_nb

      col_up  = det_ex_unq_up_orb_1st_pos (det_up_i)
      orb_up  = det_ex_unq_up_orb_2nd_lab (det_up_i)
      det_ref_up = iwdet_ex_ref_up (det_up_i)

      do l = 1, nup
       ratio_up (l) = 0.d0
       do i = 1, nup
         ratio_up (l) = ratio_up (l) + slater_mat_trans_inv_up (i, l, det_ref_up) * orb (i, orb_up)
       enddo  ! i
      enddo ! l

      factor_up_inv = 1.d0/ratio_up (col_up)
      ratio_up (col_up) = ratio_up (col_up) - 1.d0

      do k = 1, nup
        do j = 1 , nup
         slater_mat_ex_trans_inv_up (k, j, det_up_i) = slater_mat_trans_inv_up (k, j, det_ref_up) &
                                 - slater_mat_trans_inv_up (k, col_up, det_ref_up) * ratio_up(j)* factor_up_inv
        enddo ! j
      enddo ! k

   enddo ! det_up_i

! spin dn
  do det_dn_i = 1, det_ex_unq_dn_nb

      col_dn  = det_ex_unq_dn_orb_1st_pos (det_dn_i)
      orb_dn  = det_ex_unq_dn_orb_2nd_lab (det_dn_i)
      det_ref_dn = iwdet_ex_ref_dn (det_dn_i)

      do l = 1, ndn
       ratio_dn (l) = 0.d0
       do i = 1, ndn
         ratio_dn (l) = ratio_dn (l) + slater_mat_trans_inv_dn (i, l, det_ref_dn) * orb (nup + i, orb_dn)
       enddo  ! i
      enddo ! l

      factor_dn_inv = 1.d0/ratio_dn (col_dn)
      ratio_dn (col_dn) = ratio_dn (col_dn) - 1.d0

      do k = 1, ndn
        do j = 1 , ndn
         slater_mat_ex_trans_inv_dn (k, j, det_dn_i) = slater_mat_trans_inv_dn (k, j, det_ref_dn) &
                                 - slater_mat_trans_inv_dn (k, col_dn, det_ref_dn) * ratio_dn(j)* factor_dn_inv
        enddo ! j
      enddo ! k

   enddo ! det_dn_i

!  Post-conditions
!   call object_provide ('slater_mat_ex_trans_inv_up')
!   call object_provide ('slater_mat_ex_trans_inv_dn')

!   call object_write_2 ('slater_mat_ex_trans_inv_up', 'slater_mat_ex_trans_inv_up_2')
!   call object_write_2 ('slater_mat_ex_trans_inv_dn', 'slater_mat_ex_trans_inv_dn_2')
!!!   call is_equal_or_die (slater_mat_ex_trans_inv_up, slater_mat_ex_trans_inv_up_2, 1.d-8)
!!!   call is_equal_or_die (slater_mat_ex_trans_inv_dn, slater_mat_ex_trans_inv_dn_2, 1.d-8)

  end subroutine slater_mat_ex_trans_inv_3_bld

! ==============================================================================
  subroutine grd_det_ex_unq_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of unique spin-up and down excited determinants
!
! Created       : J. Toulouse, 08 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_i, dim_i, i, j

! header
  if (header_exe) then

   call object_create ('grd_det_ex_unq_up')
   call object_create ('grd_det_ex_unq_dn')

   call object_needed ('ndim')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('det_ex_unq_up_nb')
   call object_needed ('det_ex_unq_dn_nb')
   call object_needed ('slater_mat_ex_trans_inv_up')
   call object_needed ('slater_mat_ex_trans_inv_dn')
   call object_needed ('det_ex_unq_orb_lab_up')
   call object_needed ('det_ex_unq_orb_lab_dn')
   call object_needed ('dorb')
   call object_needed ('det_ex_unq_up')
   call object_needed ('det_ex_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_det_ex_unq_up', grd_det_ex_unq_up, ndim, nup, det_ex_unq_up_nb)
  call object_alloc ('grd_det_ex_unq_dn', grd_det_ex_unq_dn, ndim, ndn, det_ex_unq_dn_nb)

  grd_det_ex_unq_up = 0.d0
  grd_det_ex_unq_dn = 0.d0

! spin up determinants
  do det_i = 1, det_ex_unq_up_nb
   do dim_i = 1, ndim
    do i = 1, nup
     do j = 1, nup
      grd_det_ex_unq_up (dim_i, i, det_i) = grd_det_ex_unq_up (dim_i, i, det_i)  +  &
        slater_mat_ex_trans_inv_up (i, j, det_i) * dorb (dim_i, i, det_ex_unq_orb_lab_up (j, det_i))
     enddo
        grd_det_ex_unq_up (dim_i, i, det_i) = grd_det_ex_unq_up (dim_i, i, det_i) * det_ex_unq_up (det_i)
    enddo

   enddo ! dim_i
  enddo ! det_i

! spin down determinants
  do det_i = 1, det_ex_unq_dn_nb
   do dim_i = 1, ndim
    do i = 1, ndn
     do j = 1, ndn
      grd_det_ex_unq_dn (dim_i, i, det_i) = grd_det_ex_unq_dn (dim_i, i, det_i)  +  &
        slater_mat_ex_trans_inv_dn (i, j, det_i) * dorb (dim_i, nup + i, det_ex_unq_orb_lab_dn (j, det_i))
     enddo
        grd_det_ex_unq_dn (dim_i, i, det_i) = grd_det_ex_unq_dn (dim_i, i, det_i) * det_ex_unq_dn (det_i)
    enddo
   enddo ! dim_i
  enddo ! det_i

! write(6,*) trim(here), ': grd_det_unq_ex_up=',grd_det_ex_unq_up
! write(6,*) trim(here), ': grd_det_unq_ex_dn=',grd_det_ex_unq_dn
!
! write(6,*) trim(here), ': slater_mat_ex_trans_inv_up(1,1,1)=',slater_mat_ex_trans_inv_up(1,1,1)
! write(6,*) trim(here), ': det_ex_unq_up(1)=',det_ex_unq_up(1)

!  write(6,*) trim(here), ': grd_det_unq_ex_up(1,1,1)=',grd_det_ex_unq_up(1,1,1)
!  write(6,*) trim(here), ': dorb(1,1,2)=',dorb(1,1,2)
!  write(6,*) trim(here), ': grd_det_unq_ex_up(2,1,1)=',grd_det_ex_unq_up(2,1,1)
!  write(6,*) trim(here), ': dorb(2,1,2)=',dorb(2,1,2)
!  write(6,*) trim(here), ': grd_det_unq_ex_up(3,1,1)=',grd_det_ex_unq_up(3,1,1)
!  write(6,*) trim(here), ': dorb(3,1,2)=',dorb(3,1,2)
!
!  write(6,*) trim(here), ': grd_det_unq_ex_up(1,1,2)=',grd_det_ex_unq_up(1,1,2)
!  write(6,*) trim(here), ': dorb(1,1,3)=',dorb(1,1,3)

  end subroutine grd_det_ex_unq_bld

! ==============================================================================
  subroutine lap_det_ex_unq_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of unique spin-up and down excited determinants
!
! Created       : J. Toulouse, 08 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_i, dim_i, i, j

! header
  if (header_exe) then

   call object_create ('lap_det_ex_unq_up')
   call object_create ('lap_det_ex_unq_dn')

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('det_ex_unq_up_nb')
   call object_needed ('det_ex_unq_dn_nb')
   call object_needed ('slater_mat_ex_trans_inv_up')
   call object_needed ('slater_mat_ex_trans_inv_dn')
   call object_needed ('det_ex_unq_orb_lab_up')
   call object_needed ('det_ex_unq_orb_lab_dn')
   call object_needed ('ddorb')
   call object_needed ('det_ex_unq_up')
   call object_needed ('det_ex_unq_dn')


   return

  endif

! begin

! allocations
  call object_alloc ('lap_det_ex_unq_up', lap_det_ex_unq_up, nup, det_ex_unq_up_nb)
  call object_alloc ('lap_det_ex_unq_dn', lap_det_ex_unq_dn, ndn, det_ex_unq_dn_nb)

  lap_det_ex_unq_up = 0.d0
  lap_det_ex_unq_dn = 0.d0

! spin up determinants
   do det_i = 1, det_ex_unq_up_nb
      do i = 1, nup
       do j = 1, nup
        lap_det_ex_unq_up (i, det_i) = lap_det_ex_unq_up (i, det_i)  +  &
          slater_mat_ex_trans_inv_up (i, j, det_i) * ddorb (i, det_ex_unq_orb_lab_up (j, det_i))
       enddo
        lap_det_ex_unq_up (i, det_i) = lap_det_ex_unq_up (i, det_i)  * det_ex_unq_up (det_i)
      enddo
    enddo  ! det_i

! spin down determinants
   do det_i = 1, det_ex_unq_dn_nb
      do i = 1, ndn
       do j = 1, ndn
        lap_det_ex_unq_dn (i, det_i) = lap_det_ex_unq_dn (i, det_i)  +  &
          slater_mat_ex_trans_inv_dn (i, j, det_i) * ddorb (nup + i, det_ex_unq_orb_lab_dn (j, det_i))
       enddo
        lap_det_ex_unq_dn (i, det_i) = lap_det_ex_unq_dn (i, det_i)  * det_ex_unq_dn (det_i)
      enddo
    enddo ! det_i

  end subroutine lap_det_ex_unq_bld

! ==============================================================================
  subroutine grd_det_ex_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of all excited determinants
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_i
  integer dim_i
  integer ex_i
  integer iwdet
  integer sgn
  integer csf_i, det_in_csf_i
  integer det_unq_up_i, det_unq_dn_i

! header
  if (header_exe) then

   call object_create ('grd_det_ex_up')
   call object_create ('grd_det_ex_dn')

   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('grd_det_unq_up')
   call object_needed ('grd_det_unq_dn')
   call object_needed ('grd_det_ex_unq_up')
   call object_needed ('grd_det_ex_unq_dn')
   call object_needed ('det_ex_orb_lab_srt_sgn_up')
   call object_needed ('det_ex_orb_lab_srt_sgn_dn')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_det_ex_up', grd_det_ex_up, ndim, nup, single_ex_nb, ndet)
  call object_alloc ('grd_det_ex_dn', grd_det_ex_dn, ndim, ndn, single_ex_nb, ndet)

  do ex_i = 1, single_ex_nb

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

        det_i = iwdet_in_csf (det_in_csf_i, csf_i)
        det_unq_up_i = det_to_det_unq_up (det_i)
        det_unq_dn_i = det_to_det_unq_dn (det_i)

!       spin-up excited determinants:

!       current excited determinant is zero
        if (ex_up_is_zero (ex_i, det_unq_up_i)) then
         grd_det_ex_up (:,:,ex_i, det_i) = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif (.not. is_det_ex_up (ex_i, det_unq_up_i)) then
         iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i)
         grd_det_ex_up (:, :, ex_i, det_i) = sgn * grd_det_unq_up (:,:,iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i) * det_ex_unq_orb_lab_srt_sgn_up (iwdet)
         grd_det_ex_up (:,:, ex_i, det_i) = sgn * grd_det_ex_unq_up (:,:,iwdet)
        endif

!       spin-down excited determinants:

!       current excited determinant is zero
        if (ex_dn_is_zero (ex_i, det_unq_dn_i)) then
         grd_det_ex_dn (:,:,ex_i, det_i) = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif (.not. is_det_ex_dn (ex_i, det_unq_dn_i)) then
         iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i)
         grd_det_ex_dn (:, :, ex_i, det_i) = sgn * grd_det_unq_dn (:,:,iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i) * det_ex_unq_orb_lab_srt_sgn_dn (iwdet)
         grd_det_ex_dn (:,:, ex_i, det_i) = sgn * grd_det_ex_unq_dn (:,:,iwdet)
        endif

     enddo ! det_in_csf_i
   enddo ! csf_i

 enddo ! ex_i

 end subroutine grd_det_ex_bld

! ==============================================================================
  subroutine lap_det_ex_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of all excited determinants
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer det_i
  integer csf_i, det_in_csf_i
  integer ex_i, iwdet, sgn
  integer det_unq_up_i, det_unq_dn_i

! header
  if (header_exe) then

   call object_create ('lap_det_ex_up')
   call object_create ('lap_det_ex_dn')

   call object_needed ('nelec')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('lap_det_unq_up')
   call object_needed ('lap_det_unq_dn')
   call object_needed ('lap_det_ex_unq_up')
   call object_needed ('lap_det_ex_unq_dn')
   call object_needed ('det_ex_orb_lab_srt_sgn_up')
   call object_needed ('det_ex_orb_lab_srt_sgn_dn')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_det_ex_up', lap_det_ex_up, nup, single_ex_nb, ndet)
  call object_alloc ('lap_det_ex_dn', lap_det_ex_dn, ndn, single_ex_nb, ndet)

  do ex_i = 1, single_ex_nb

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

        det_i = iwdet_in_csf (det_in_csf_i, csf_i)
        det_unq_up_i = det_to_det_unq_up (det_i)
        det_unq_dn_i = det_to_det_unq_dn (det_i)

!       spin-up excited determinants:

!       current excited determinant is zero
        if (ex_up_is_zero (ex_i, det_unq_up_i)) then
         lap_det_ex_up (:,ex_i, det_i) = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif (.not. is_det_ex_up (ex_i, det_unq_up_i)) then
         iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i)
         lap_det_ex_up (:, ex_i, det_i) = sgn * lap_det_unq_up (:,iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i) * det_ex_unq_orb_lab_srt_sgn_up (iwdet)
         lap_det_ex_up (:, ex_i, det_i) = sgn * lap_det_ex_unq_up (:,iwdet)
        endif

!       spin-down excited determinants:

!       current excited determinant is zero
        if (ex_dn_is_zero (ex_i, det_unq_dn_i)) then
         lap_det_ex_dn (:,ex_i, det_i) = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif (.not. is_det_ex_dn (ex_i, det_unq_dn_i)) then
         iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i)
         lap_det_ex_dn ( :, ex_i, det_i) = sgn * lap_det_unq_dn (:,iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i) * det_ex_unq_orb_lab_srt_sgn_dn (iwdet)
         lap_det_ex_dn (:, ex_i, det_i) = sgn * lap_det_ex_unq_dn (:,iwdet)
        endif

     enddo ! det_in_csf_i
   enddo ! csf_i

 enddo ! ex_i

 end subroutine lap_det_ex_bld

! ==============================================================================
  subroutine grd_psid_ex_over_psid_bld
! ------------------------------------------------------------------------------
! Description   : (Gradient psid excited)/psid excited
!
! Created       : J. Toulouse, 08 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer ex_i, ex_rev_i
  integer dorb_i, det_i, dim_i, det_unq_up_i, det_unq_dn_i
  integer elec_i, elec_up_i, elec_dn_i
  integer csf_i, det_in_csf_i
  real(dp) coefficient, grd_det

! header
  if (header_exe) then

   call object_create ('grd_psid_ex_over_psid')

   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_ind_rev')
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
   call object_needed ('det_ex_up')
   call object_needed ('det_ex_dn')
   call object_needed ('grd_det_unq_up')
   call object_needed ('grd_det_unq_dn')
   call object_needed ('grd_det_ex_up')
   call object_needed ('grd_det_ex_dn')
   call object_needed ('psid_ex')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_psid_ex_over_psid', grd_psid_ex_over_psid, ndim, nelec, param_orb_nb)

  grd_psid_ex_over_psid = 0.d0

  do dorb_i = 1, param_orb_nb

    ex_i = ex_orb_ind (dorb_i)
    ex_rev_i = ex_orb_ind_rev (dorb_i)

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

!        Eij Det
         grd_det = grd_det_ex_up (dim_i, elec_up_i, ex_i, det_i) * detd (det_unq_dn_i) +  &
                   grd_det_unq_up (dim_i, elec_up_i, det_unq_up_i) * det_ex_dn (ex_i, det_i)

         if (.not. l_casscf .and. ex_rev_i /= 0) then
!         (Eij-Eji) Det
          grd_det = grd_det - ( grd_det_ex_up (dim_i, elec_up_i, ex_rev_i, det_i) * detd (det_unq_dn_i) +  &
                                grd_det_unq_up (dim_i, elec_up_i, det_unq_up_i) * det_ex_dn (ex_rev_i, det_i) )
          endif

         grd_psid_ex_over_psid (dim_i, elec_i, dorb_i) = grd_psid_ex_over_psid (dim_i, elec_i, dorb_i) + coefficient * grd_det

        enddo

        do elec_dn_i = 1, ndn
          elec_i = elec_i + 1

!         Eij Det
          grd_det = det_ex_up (ex_i, det_i) * grd_det_unq_dn (dim_i, elec_dn_i, det_unq_dn_i)  +        &
                            detu (det_unq_up_i) *  grd_det_ex_dn (dim_i, elec_dn_i, ex_i, det_i)

         if (.not. l_casscf .and. ex_rev_i /= 0) then
!         (Eij-Eji) Det
          grd_det =  grd_det                                                                         &
                    - ( det_ex_up (ex_rev_i, det_i) * grd_det_unq_dn (dim_i, elec_dn_i, det_unq_dn_i)  +        &
                            detu (det_unq_up_i) *  grd_det_ex_dn (dim_i, elec_dn_i, ex_rev_i, det_i) )
          endif

          grd_psid_ex_over_psid (dim_i, elec_i, dorb_i) = grd_psid_ex_over_psid (dim_i, elec_i, dorb_i) + coefficient * grd_det
        enddo

       enddo ! dim_i

     enddo ! det_in_csf_i
   enddo ! csf_i

  grd_psid_ex_over_psid (:,:,dorb_i) = grd_psid_ex_over_psid (:,:,dorb_i)/psid_ex (dorb_i)

 enddo ! dorb_i

 end subroutine grd_psid_ex_over_psid_bld

! ==============================================================================
  subroutine lap_psid_ex_over_psid_bld
! ------------------------------------------------------------------------------
! Description   : (Laplacian psid excited)/psid excited
!
! Created       : J. Toulouse, 08 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer ex_i, ex_rev_i
  integer dorb_i, det_i, det_unq_up_i, det_unq_dn_i
  integer elec_i, elec_up_i, elec_dn_i
  integer csf_i, det_in_csf_i
  real(dp) coefficient, lap_det

! header
  if (header_exe) then

   call object_create ('lap_psid_ex_over_psid')

   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_ind_rev')
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
   call object_needed ('det_ex_up')
   call object_needed ('det_ex_dn')
   call object_needed ('lap_det_unq_up')
   call object_needed ('lap_det_unq_dn')
   call object_needed ('lap_det_ex_up')
   call object_needed ('lap_det_ex_dn')
   call object_needed ('psid_ex')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_psid_ex_over_psid', lap_psid_ex_over_psid, nelec, param_orb_nb)

  lap_psid_ex_over_psid (:,:) = 0.d0

  do dorb_i = 1, param_orb_nb

   ex_i = ex_orb_ind (dorb_i)
   ex_rev_i = ex_orb_ind_rev (dorb_i)

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

      det_i = iwdet_in_csf (det_in_csf_i, csf_i)
      det_unq_up_i = det_to_det_unq_up (det_i)
      det_unq_dn_i = det_to_det_unq_dn (det_i)

      coefficient = csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)

        elec_i = 0
        do elec_up_i = 1, nup
         elec_i = elec_i + 1

!          Eij Det
           lap_det = lap_det_ex_up (elec_up_i, ex_i, det_i) * detd (det_unq_dn_i) +  &
                            lap_det_unq_up (elec_up_i, det_unq_up_i) * det_ex_dn (ex_i, det_i)

         if (.not. l_casscf .and. ex_rev_i /= 0) then
!          (Eij-Eji) Det
           lap_det = lap_det                      &
                    - (lap_det_ex_up (elec_up_i, ex_rev_i, det_i) * detd (det_unq_dn_i) +  &
                            lap_det_unq_up (elec_up_i, det_unq_up_i) * det_ex_dn (ex_rev_i, det_i) )
         endif

         lap_psid_ex_over_psid (elec_i, dorb_i) = lap_psid_ex_over_psid (elec_i, dorb_i) + coefficient * lap_det
        enddo

        do elec_dn_i = 1, ndn
          elec_i = elec_i + 1

!          Eij Det
           lap_det = det_ex_up (ex_i, det_i) * lap_det_unq_dn (elec_dn_i, det_unq_dn_i)  +        &
                            detu (det_unq_up_i) * lap_det_ex_dn (elec_dn_i, ex_i, det_i)

          if (.not. l_casscf .and. ex_rev_i /= 0) then
           lap_det = lap_det                                     &
                    - (det_ex_up (ex_rev_i, det_i) * lap_det_unq_dn (elec_dn_i, det_unq_dn_i)  +        &
                            detu (det_unq_up_i) * lap_det_ex_dn (elec_dn_i, ex_rev_i, det_i))
          endif

          lap_psid_ex_over_psid (elec_i, dorb_i) = lap_psid_ex_over_psid (elec_i, dorb_i) + coefficient * lap_det
        enddo

     enddo ! det_in_csf_i
   enddo ! csf_i

  lap_psid_ex_over_psid (:,dorb_i) = lap_psid_ex_over_psid (:,dorb_i)/psid_ex (dorb_i)

 enddo ! dorb_i

 end subroutine lap_psid_ex_over_psid_bld

! ==============================================================================
  subroutine lap_lnpsid_ex_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian ( ln (psid_ex) )
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, elec_i, ex_i
  real(dp) sum

! header
  if (header_exe) then

   call object_create ('lap_lnpsid_ex')

   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('param_orb_nb')
   call object_needed ('lap_psid_ex_over_psid')
   call object_needed ('grd_psid_ex_over_psid')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_lnpsid_ex', lap_lnpsid_ex, nelec, param_orb_nb)

  do ex_i = 1, param_orb_nb
   do elec_i = 1, nelec
    sum = 0.d0
    do dim_i = 1, ndim
     sum = sum + grd_psid_ex_over_psid (dim_i, elec_i, ex_i)**2
    enddo
    lap_lnpsid_ex (elec_i, ex_i) = lap_psid_ex_over_psid (elec_i, ex_i) - sum
   enddo
  enddo

 end subroutine lap_lnpsid_ex_bld

! ==============================================================================
  subroutine sum_lap_lnpsid_ex_bld
! ------------------------------------------------------------------------------
! Description   : Sum_i lap_lnpsid (i)
!
! Created       : J. Toulouse, 05 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer elec_i, ex_i

! header
  if (header_exe) then

   call object_create ('sum_lap_lnpsid_ex')

   call object_needed ('param_orb_nb')
   call object_needed ('nelec')
   call object_needed ('lap_lnpsid_ex')

   return

  endif

! begin

! allocations
  call object_alloc ('sum_lap_lnpsid_ex', sum_lap_lnpsid_ex, param_orb_nb)

  sum_lap_lnpsid_ex (:) = 0.d0

  do ex_i = 1, param_orb_nb
   do elec_i = 1, nelec
     sum_lap_lnpsid_ex (ex_i) = sum_lap_lnpsid_ex (ex_i) + lap_lnpsid_ex (elec_i, ex_i)
   enddo
  enddo

 end subroutine sum_lap_lnpsid_ex_bld

! ==============================================================================
  subroutine grd_psi_ex_over_psi_bld
! ------------------------------------------------------------------------------
! Description   : (Gradient excited psi)/ excited psi
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer ex_i

! header
  if (header_exe) then

   call object_create ('grd_psi_ex_over_psi')

   call object_needed ('param_orb_nb')
   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('grd_psid_ex_over_psid')
   call object_needed ('vj')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_psi_ex_over_psi', grd_psi_ex_over_psi, ndim, nelec, param_orb_nb)

  do ex_i = 1, param_orb_nb
     grd_psi_ex_over_psi (1:ndim, 1:nelec, ex_i) = grd_psid_ex_over_psid (1:ndim, 1:nelec, ex_i) + vj (1:ndim, 1:nelec)
  enddo

 end subroutine grd_psi_ex_over_psi_bld

! ==============================================================================
  subroutine sum_lap_lnpsi_ex_bld
! ------------------------------------------------------------------------------
! Description   : Sum_i lap_lnpsi (i) : determinant part + Jastrow
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('sum_lap_lnpsi_ex')

   call object_needed ('param_orb_nb')
   call object_needed ('sum_lap_lnpsid_ex')
   call object_needed ('sum_lap_lnj')

   return

  endif

! begin

! allocations
  call object_alloc ('sum_lap_lnpsi_ex', sum_lap_lnpsi_ex, param_orb_nb)

  sum_lap_lnpsi_ex (:) = sum_lap_lnpsid_ex (:) + sum_lap_lnj

 end subroutine sum_lap_lnpsi_ex_bld

! ==============================================================================
  subroutine eloc_kin_ex_bld
! ------------------------------------------------------------------------------
! Description   : Kinetic local energy of excited wave function
!
! Created       : J. Toulouse, 05 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer dim_i, elec_i, ex_i
  real(dp) sum_grd_psi_over_psi_square

! header
  if (header_exe) then

   call object_create ('eloc_kin_ex')

   call object_needed ('param_orb_nb')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('sum_lap_lnpsi_ex')
   call object_needed ('grd_psi_ex_over_psi')

   return

  endif

! begin

! allocations
  call object_alloc ('eloc_kin_ex', eloc_kin_ex, param_orb_nb)

  do ex_i = 1, param_orb_nb

  sum_grd_psi_over_psi_square = 0.d0

  do dim_i = 1, ndim
    do elec_i = 1, nelec
       sum_grd_psi_over_psi_square = sum_grd_psi_over_psi_square + grd_psi_ex_over_psi (dim_i, elec_i, ex_i)**2
    enddo
  enddo

  eloc_kin_ex (ex_i) =  -0.5d0 * (sum_lap_lnpsi_ex (ex_i) + sum_grd_psi_over_psi_square)

  enddo ! ex_i

 end subroutine eloc_kin_ex_bld

! ==============================================================================
  subroutine psid_ex_in_x_bld
! ------------------------------------------------------------------------------
! Description   :  Calculate value of determinantal part of single-excited wave function
! Description   :  in a grid point of the angular quadrature grid for pseudopotential
! Description   :  for one electron, update done by Shermann-Morison formula
!
! Created       : J. Toulouse, 16 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer ex_i, ex_rev_i
  integer det_i, dorb_i, det_unq_up_i, det_unq_dn_i
  integer j, iel_up, iel_dn
  integer csf_i
  integer det_in_csf_i
  integer iwdet
  real(dp) factor_up, factor_dn
  real(dp) sgn
  real(dp) coefficient
  real(dp) det_ex_up_in_x, det_ex_dn_in_x
  real(dp) det_ex_rev_up_in_x, det_ex_rev_dn_in_x

! header
  if (header_exe) then

   call object_create ('psid_ex_in_x', psid_ex_in_x_index)

   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_ind_rev')
   call object_needed ('electron')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('csf_coef')
   call object_needed ('cdet_in_csf')
   call object_needed ('ex_up_is_zero')
   call object_needed ('ex_dn_is_zero')
   call object_needed ('iwdet_ex_up')
   call object_needed ('iwdet_ex_dn')
   call object_needed ('det_ex_orb_lab_srt_sgn_up')
   call object_needed ('det_ex_orb_lab_srt_sgn_dn')
   call object_needed ('slater_mat_ex_trans_inv_up')
   call object_needed ('slater_mat_ex_trans_inv_dn')
   call object_needed ('orbe')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('detn')
   call object_needed ('det_ex_up')
   call object_needed ('det_ex_dn')
   call object_needed ('det_ex_unq_up')
   call object_needed ('det_ex_unq_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('psid_ex_in_x', psid_ex_in_x, param_orb_nb)

  psid_ex_in_x = 0.d0

  do dorb_i = 1, param_orb_nb

   ex_i = ex_orb_ind (dorb_i)
   ex_rev_i = ex_orb_ind_rev (dorb_i)

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

        det_i = iwdet_in_csf (det_in_csf_i, csf_i)
        det_unq_up_i = det_to_det_unq_up (det_i)
        det_unq_dn_i = det_to_det_unq_dn (det_i)

        coefficient = csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)

!       electron is spin up
        if (electron <= nup) then

        iel_up = electron

!       current excited determinant is zero
        if (ex_up_is_zero (ex_i, det_unq_up_i)) then
         det_ex_up_in_x = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif ( .not. is_det_ex_up (ex_i, det_unq_up_i) ) then
         iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i)
         det_ex_up_in_x = sgn * detn (iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i) * det_ex_unq_orb_lab_srt_sgn_up (iwdet)

         factor_up = 0.d0
         do j = 1, nup
         factor_up = factor_up + slater_mat_ex_trans_inv_up (iel_up, j, iwdet) * orbe (det_ex_unq_orb_lab_up (j, iwdet))
         enddo

         det_ex_up_in_x  = sgn * factor_up * det_ex_unq_up (iwdet)

        endif

         psid_ex_in_x (dorb_i) = psid_ex_in_x (dorb_i) + coefficient * (det_ex_up_in_x * detd (det_unq_dn_i) + detn(det_unq_up_i) * det_ex_dn (ex_i, det_i))


!       reverse excitation for active-active non-CASSCF excitations
        if (.not. l_casscf .and. ex_rev_i /=0 ) then

!       current excited determinant is zero
        if (ex_up_is_zero (ex_rev_i, det_unq_up_i)) then
         det_ex_rev_up_in_x = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif (.not. is_det_ex_up (ex_rev_i, det_unq_up_i)) then
         iwdet = iwdet_ex_up (ex_rev_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_rev_i, det_unq_up_i)
         det_ex_rev_up_in_x = sgn * detn (iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_up (ex_rev_i, det_unq_up_i)
         sgn = det_ex_orb_lab_srt_sgn_up (ex_rev_i, det_unq_up_i) * det_ex_unq_orb_lab_srt_sgn_up (iwdet)

         factor_up = 0.d0
         do j = 1, nup
         factor_up = factor_up +  &
          slater_mat_ex_trans_inv_up (iel_up, j, iwdet) * orbe (det_ex_unq_orb_lab_up (j, iwdet))
         enddo

         det_ex_rev_up_in_x  = sgn * factor_up * det_ex_unq_up (iwdet)

        endif

         psid_ex_in_x (dorb_i) = psid_ex_in_x (dorb_i) - coefficient *   &
                   (det_ex_rev_up_in_x * detd (det_unq_dn_i) + detn(det_unq_up_i) * det_ex_dn (ex_rev_i, det_i))

        endif ! if .not. l_cascsf


!       electron is spin down
        else
        iel_dn = electron - nup


!       current excited determinant is zero
        if (ex_dn_is_zero (ex_i, det_unq_dn_i)) then
         det_ex_dn_in_x = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif ( .not. is_det_ex_dn (ex_i, det_unq_dn_i) ) then
         iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i)
         det_ex_dn_in_x = sgn * detn (iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i) * det_ex_unq_orb_lab_srt_sgn_dn (iwdet)

         factor_dn = 0.d0
         do j = 1, ndn
         factor_dn = factor_dn +  &
          slater_mat_ex_trans_inv_dn (iel_dn, j, iwdet) * orbe (det_ex_unq_orb_lab_dn (j, iwdet))
         enddo

         det_ex_dn_in_x  = sgn * factor_dn * det_ex_unq_dn (iwdet)

        endif

        psid_ex_in_x (dorb_i) = psid_ex_in_x (dorb_i) + coefficient * (det_ex_up (ex_i, det_i) * detn (det_unq_dn_i) + detu(det_unq_up_i) * det_ex_dn_in_x)

!       reverse excitation for active-active non-CASSCF excitations
        if (.not. l_casscf .and. ex_rev_i /=0 ) then

!       current excited determinant is zero
        if (ex_dn_is_zero (ex_rev_i, det_unq_dn_i)) then
         det_ex_rev_dn_in_x = 0.d0

!       current excited determinant is a determinant in the ground-state wave function
        elseif ( .not. is_det_ex_dn (ex_rev_i, det_unq_dn_i) ) then
         iwdet = iwdet_ex_dn (ex_rev_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_rev_i, det_unq_dn_i)
         det_ex_rev_dn_in_x = sgn * detn (iwdet)

!       current excited determinant is a true excited determinant
        else
         iwdet = iwdet_ex_dn (ex_rev_i, det_unq_dn_i)
         sgn = det_ex_orb_lab_srt_sgn_dn (ex_rev_i, det_unq_dn_i) * det_ex_unq_orb_lab_srt_sgn_dn (iwdet)

         factor_dn = 0.d0
         do j = 1, ndn
         factor_dn = factor_dn + slater_mat_ex_trans_inv_dn (iel_dn, j, iwdet) * orbe (det_ex_unq_orb_lab_dn (j, iwdet))
         enddo

         det_ex_rev_dn_in_x  = sgn * factor_dn * det_ex_unq_dn (iwdet)

        endif

        psid_ex_in_x (dorb_i) = psid_ex_in_x (dorb_i) - coefficient * (det_ex_up (ex_rev_i, det_i) * detn (det_unq_dn_i) + detu(det_unq_up_i) * det_ex_rev_dn_in_x)

        endif ! if .not. l_cascsf

        endif ! electron

     enddo ! det_in_csf_i
   enddo ! csf_i

 enddo ! dorb_i

 end subroutine psid_ex_in_x_bld

! ==============================================================================
  subroutine eloc_pot_nloc_ex_bld
! ------------------------------------------------------------------------------
! Description   : nonlocal pseudopoential contribution to local energy for excited wave functions
! Description   : eloc_pot_nloc_ex = vpsp_ex / psid_ex
! Description   : where vpsp_ex = (V_nonloc D_ex J)/J and psid_ex = D_ex = determinantal part of excited wave functions
!
! Created       : J. Toulouse, 16 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer ex_i

! header
  if (header_exe) then

   call object_create ('eloc_pot_nloc_ex')

   call object_needed ('param_orb_nb')
   call object_needed ('vpsp_ex')
   call object_needed ('psid_ex')

   return

  endif

! begin

! allocations
  call object_alloc ('eloc_pot_nloc_ex', eloc_pot_nloc_ex, param_orb_nb)

  eloc_pot_nloc_ex = vpsp_ex / psid_ex

 end subroutine eloc_pot_nloc_ex_bld

! ==============================================================================
  subroutine eloc_pot_ex_bld
! ------------------------------------------------------------------------------
! Description   : total local potential energy for excited wave functions
! Description   : with or without (nonlocal) pseudopotential
!
! Created       : J. Toulouse, 15 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'eloc_pot_ex_bld'

! header
  if (header_exe) then

   call object_create ('eloc_pot_ex')

   call object_needed ('param_orb_nb')
   call object_needed ('nloc')

   return

  endif

! begin

! allocations
  call object_alloc ('eloc_pot_ex', eloc_pot_ex, param_orb_nb)


! without a pseudopotential
  if (nloc <= 0) then

   call object_provide_by_index (eloc_pot_ex_bld_index, eloc_pot_index)
   eloc_pot_ex = eloc_pot

! with a pseudopotential
  else

   call object_provide_by_index (eloc_pot_ex_bld_index, eloc_pot_loc_index)
   call object_provide_by_index (eloc_pot_ex_bld_index, eloc_pot_nloc_ex_index)
   eloc_pot_ex = eloc_pot_loc + eloc_pot_nloc_ex

  endif

 end subroutine eloc_pot_ex_bld

! ==============================================================================
  subroutine eloc_ex_bld
! ------------------------------------------------------------------------------
! Description   : Total local energy of excited wave functions (without weight)
! Description   : (H Psi_ex)/Psi_ex
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('eloc_ex')

   call object_needed ('param_orb_nb')
   call object_needed ('eloc_kin_ex')
   call object_needed ('eloc_pot_ex')

   return

  endif

! begin

! allocations
  call object_alloc ('eloc_ex', eloc_ex, param_orb_nb)

  eloc_ex (:) = eloc_kin_ex (:) + eloc_pot_ex

 end subroutine eloc_ex_bld

! ==============================================================================
  subroutine deloc_orb_bld
! ------------------------------------------------------------------------------
! Description   :  derivative of local energy wrt orbital rotation parameters kappa
! Description   :  dEloc/dkappa = [ (H dPsi/dkappa)/(dPsi/dkappa) - Eloc ] dlnPsi/dkappa
!
! Created       : J. Toulouse, 14 Jan 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! header
  if (header_exe) then

   call object_create ('deloc_orb')

   call object_needed ('param_orb_nb')
   call object_needed ('eloc_ex')
   call object_needed ('eloc')
   call object_needed ('dpsi_orb')

   return

  endif

! begin

! allocations
  call object_alloc ('deloc_orb', deloc_orb, param_orb_nb)

  deloc_orb (:) = (eloc_ex (:) - eloc) * dpsi_orb (:)

 end subroutine deloc_orb_bld

! ==============================================================================
  subroutine delta_eps_bld
! ------------------------------------------------------------------------------
! Description  : orbital energy difference for single excitations
!
! Created      : J. Toulouse, 04 Nov 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer ex_i
  integer orb_1st, orb_2nd
  real(dp) eps_1st, eps_2nd


! header
  if (header_exe) then

   call object_create ('delta_eps')

   call object_needed ('param_orb_nb')
   call object_needed ('orb_energies')
   call object_needed ('ex_orb_1st_lab')
   call object_needed ('ex_orb_2nd_lab')
   call object_needed ('ndet')

   return

  endif

! begin

! allocations
  call object_alloc ('delta_eps', delta_eps, param_orb_nb)

  do ex_i = 1, param_orb_nb
     orb_1st = ex_orb_1st_lab (ex_i)
     orb_2nd = ex_orb_2nd_lab (ex_i)
     eps_1st   = orb_energies (orb_1st)
     eps_2nd   = orb_energies (orb_2nd)
     delta_eps (ex_i) = eps_2nd - eps_1st
  enddo ! ex_i


! correction of closed -> active and active -> virtual excitation energies for MCSCF wave functions
  if (ndet > 1) then
    call object_provide ('orb_cls_in_wf')
    call object_provide ('orb_act_in_wf')
    call object_provide ('dens_mat_wfdet')
    do ex_i = 1, param_orb_nb
     orb_1st = ex_orb_1st_lab (ex_i)
     orb_2nd = ex_orb_2nd_lab (ex_i)

     if ((orb_cls_in_wf (orb_1st) .and. orb_act_in_wf (orb_2nd)) .or.  &
        ((orb_act_in_wf (orb_1st) .and. orb_vir_in_wf (orb_2nd)))) then

     delta_eps (ex_i) = delta_eps (ex_i) +  &
       lambda*( dens_mat_wfdet(orb_2nd,orb_2nd) + 2 - dens_mat_wfdet(orb_1st,orb_1st) )/2.d0

     endif
    enddo ! ex_i
  endif


  end subroutine delta_eps_bld

end module deriv_orb_mod
