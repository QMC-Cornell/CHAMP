module deriv_orb_mod

  use deriv_fast_mod
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

  integer                                :: orb_mix_lab_nb
  logical, allocatable                   :: orb_mix_lab (:)

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
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_up (:, :, :)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_dn (:, :, :)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_up_2 (:,:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_dn_2 (:,:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_up_3 (:,:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_dn_3 (:,:,:)
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

  logical                                :: l_slater_mat_ex_trans_inv_sm = .true.

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
  use all_modules_mod
  implicit none

! local
  integer  ex_i, orb_i, dorb_i
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
   call object_create ('orb_mix_lab')
   call object_create ('orb_mix_lab_nb')

   call object_needed ('orb_tot_nb')
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

  else ! ndet \= 1

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

!      direct excitation i->j
       param_orb_nb = param_orb_nb + 1
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
       call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_ind (param_orb_nb) = single_ex_nb
       ex_orb_1st_lab (single_ex_nb) = orb_opt_act_lab_i
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_act_lab_j

!      reverse excitation j->i
       single_ex_nb = single_ex_nb + 1
       call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
       call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
       ex_orb_1st_lab (single_ex_nb) = orb_opt_act_lab_j
       ex_orb_2nd_lab (single_ex_nb) = orb_opt_act_lab_i

!      if orthogonality constraint is imposed, then the reverse excitation is combined with the direct excitation in one orbital derivative
       if (l_active_orb_ortho_constraint) then
        ex_orb_ind_rev (param_orb_nb) = single_ex_nb
!       ex_orb_ind_rev (param_orb_nb) = 0 ! desactivate reverse excitation

!      if orthogonality constraint is not imposed, then the reverse excitation is a new independent orbtital derivative
       else
        param_orb_nb = param_orb_nb + 1
        call object_alloc ('ex_orb_ind', ex_orb_ind, param_orb_nb)
        call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
        ex_orb_ind (param_orb_nb) = single_ex_nb
        ex_orb_ind_rev (param_orb_nb) = 0
       endif

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

! determine (occupied or vitual) orbitals mixed in orbital optimization
  call object_alloc ('orb_mix_lab', orb_mix_lab, orb_tot_nb)
  orb_mix_lab (:) = .false.
  do ex_i = 1, single_ex_nb
   orb_mix_lab (ex_orb_1st_lab (ex_i)) = .true.
   orb_mix_lab (ex_orb_2nd_lab (ex_i)) = .true.
  enddo ! ex_i
  orb_mix_lab_nb = 0
  do orb_i = 1, orb_tot_nb
    if (orb_mix_lab (orb_i)) then
      orb_mix_lab_nb = orb_mix_lab_nb + 1
    endif
  enddo

  write(6,'(a,i10)') ' Number of single orbital excitations = ', single_ex_nb
  write(6,'(a,i10)') ' Number of orbital derivatives        = ', param_orb_nb
  write(6,'(a,i10)') ' Number of orbitals involved          = ', orb_mix_lab_nb


! excitation pairs number
  deriv_orb_pairs_nb = param_orb_nb * (param_orb_nb + 1) / 2

! temporary: new routine
  if (l_check_redundant_orbital_derivative) then

!  temporary: save created objects for comparison with new routine
   param_orb_nb_old = param_orb_nb
   single_ex_nb_old = single_ex_nb
   call copy (ex_orb_ind, ex_orb_ind_old)
   call copy (ex_orb_ind_rev, ex_orb_ind_rev_old)
   call copy (ex_orb_1st_lab, ex_orb_1st_lab_old)
   call copy (ex_orb_2nd_lab, ex_orb_2nd_lab_old)

!  temporary: call new routine checking for redundancies
   call node_exe ('single_ex_wf_bld_2')

  else
  endif

! write(6,*)
! do dorb_i = 1, param_orb_nb
!   write(6,'(''dorb_i, ex_orb_ind(dorb_i), ex_orb_ind_rev(dorb_i), ex_orb_1st_lab (ex_orb_ind(dorb_i)), ex_orb_2nd_lab (ex_orb_ind(dorb_i)), ex_orb_1st_lab (ex_orb_ind_rev(dorb_i)), ex_orb_2nd_lab (ex_orb_ind_rev(dorb_i)),'', 9i5)') dorb_i, ex_orb_ind(dorb_i), ex_orb_ind_rev(dorb_i), ex_orb_1st_lab (ex_orb_ind(dorb_i)), ex_orb_2nd_lab (ex_orb_ind(dorb_i)), ex_orb_1st_lab (ex_orb_ind_rev(dorb_i)), ex_orb_2nd_lab (ex_orb_ind_rev(dorb_i))
! enddo

  if (l_print_orbital_excitations .and. .not. l_check_redundant_orbital_derivative) then
  write(6,*)
  write(6,'(a)') ' Orbital excitations:'
   do dorb_i = 1, param_orb_nb
    write(6,'(a,i4,a,i4,a,i4,a,i4)') ' orbital parameter # ',dorb_i,' corresponds to single excitation # ',ex_orb_ind (dorb_i),' : ',ex_orb_1st_lab (ex_orb_ind (dorb_i)),' -> ', ex_orb_2nd_lab (ex_orb_ind (dorb_i))
    if (ex_orb_ind_rev (dorb_i) /= 0) then
     write(6,'(a,i4,a,i4,a,i4)') '                             and reverse single excitation # ',ex_orb_ind_rev (dorb_i),' : ',ex_orb_1st_lab (ex_orb_ind_rev (dorb_i)),' -> ', ex_orb_2nd_lab (ex_orb_ind_rev (dorb_i))
    endif
   enddo ! dorb_i
  write(6,*)
  endif

 end subroutine single_ex_wf_bld

! ==============================================================================
  subroutine single_ex_wf_bld_2
! ------------------------------------------------------------------------------
! Description   : Build single orbital excitation list for general wave function
! Description   : taking into account some (but not all) redundancies
! Description   : For large CASCF wave functions, this routine does not seem to recognize
! Description   : recognize some redundant active-> active excitations -> has to be improved
! Description   : This routine is to replace single_ex_wf_bld
! Description   : For now, the routine is only used as replacement for single_ex_wf_bld
! Description   : information gathered on excitations in this routine should be used in the program
! Description   : in lieu of single_ex_det_bld
!
! Created       : J. Toulouse, 25 Oct 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  use hash_table_mod
  implicit none

! local
  integer orb_opt_i, orb_opt_j, orb_opt_lab_i, orb_opt_lab_j, dorb_i, orb_i
  integer ex_dir_rev, ex_dir_rev_nonortho, ex_dir_rev_nb, ex_dir_rev_nonortho_nb
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
  integer det_i, det_ex_sgn_up, det_ex_sgn_dn

  integer det_unq_in_csf_ex_cur_nb
  integer, allocatable :: det_unq_up_in_csf_ex_cur (:)
  integer, allocatable :: det_unq_dn_in_csf_ex_cur (:)
  real(dp), allocatable :: cdet_unq_in_csf_ex_cur (:)
  integer det_unq_in_csf_ex_cur_temp_nb, det_unq_in_csf_ex_cur_i
  integer, allocatable :: det_unq_up_in_csf_ex_cur_temp (:)
  integer, allocatable :: det_unq_dn_in_csf_ex_cur_temp (:)
  real(dp), allocatable :: cdet_unq_in_csf_ex_cur_temp (:)

  integer csf_i, csf_k, csf_ex_unq_i, csf_cur
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

  integer :: single_ex_max 
  integer :: det_ex_up_max 
  integer :: det_ex_dn_max 
  integer :: csf_ex_unq_max

  type(bucket), allocatable :: up_dets(:), dn_dets(:), csfs(:)
  logical :: success

! tenp
 !  integer ex_j, orb_1st, orb_2nd
!   logical ex_found

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

! Most arrays will be adjusted at the end of the routine and have the following max sizes:
  single_ex_max  = nup*(orb_tot_nb-nup)+ndn*(orb_tot_nb-ndn)
  det_ex_up_max  = ndetup*single_ex_max
  det_ex_dn_max  = ndetdn*single_ex_max
  csf_ex_unq_max = ncsf*single_ex_max

  call object_alloc ('det_ex_unq_ref_up', det_ex_unq_ref_up, det_ex_up_max)
  call object_alloc ('det_ex_unq_orb_lab_srt_up', det_ex_unq_orb_lab_srt_up, nup, det_ex_up_max)
  call object_alloc ('det_ex_unq_orb_lab_srt_sgn_up', det_ex_unq_orb_lab_srt_sgn_up, det_ex_up_max)

  call object_alloc ('det_ex_unq_ref_dn', det_ex_unq_ref_dn, det_ex_dn_max)
  call object_alloc ('det_ex_unq_orb_lab_srt_dn', det_ex_unq_orb_lab_srt_dn, ndn, det_ex_dn_max)
  call object_alloc ('det_ex_unq_orb_lab_srt_sgn_dn', det_ex_unq_orb_lab_srt_sgn_dn, det_ex_dn_max)

  call object_alloc ('csf_ex_unq_ref', csf_ex_unq_ref, csf_ex_unq_max)
  call object_alloc ('det_unq_in_csf_ex_unq_nb', det_unq_in_csf_ex_unq_nb, csf_ex_unq_max)
  call object_alloc ('det_unq_up_in_csf_ex_unq', det_unq_up_in_csf_ex_unq, csf_ex_unq_max)
  call object_alloc ('det_unq_dn_in_csf_ex_unq', det_unq_dn_in_csf_ex_unq, csf_ex_unq_max)
  call object_alloc ('cdet_unq_in_csf_ex_unq', cdet_unq_in_csf_ex_unq, csf_ex_unq_max)


  !Use hash table to speed up finding redundant orbital params TA
  allocate(up_dets(det_ex_up_max+ndetup))
  do det_i=1, ndetup
    call hash_table_add(up_dets, det_unq_orb_lab_srt_up(:,det_i), det_i)
  enddo

  allocate(dn_dets(det_ex_dn_max+ndetdn))
  do det_i=1, ndetdn
    call hash_table_add(dn_dets, det_unq_orb_lab_srt_dn(:,det_i), det_i)
  enddo

  if (l_active_orb_ortho_constraint) then
   ex_dir_rev_nb = 2
   ex_dir_rev_nonortho_nb = 1
  else
   ex_dir_rev_nb = 1
   ex_dir_rev_nonortho_nb = 2
  endif

! initializations
  single_ex_nb = 0
  param_orb_nb = 0
  ex_orb_ind_rev_cur = 0
  single_ex_added_nb = 0

! loops over orbitals
  do orb_opt_i = 1, orb_opt_nb
   orb_opt_lab_i = orb_opt_lab (orb_opt_i)

   do orb_opt_j = orb_opt_i+1, orb_opt_nb
    orb_opt_lab_j = orb_opt_lab (orb_opt_j)

!    write(6,'(2a,i4,a,i4)') trim(here),': considering orbital excitation ', orb_opt_lab_i, ' -> ', orb_opt_lab_j

!   skip user-specified forbidden excitations
    if (orb_ex_forbidden (orb_opt_lab_i, orb_opt_lab_j)) then
      cycle
    endif

!   skip excitations between orbitals of different symmetry
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
!   one loop or the other is used depending whether the orthogonality constraint is imposed or not
    do ex_dir_rev_nonortho = 1, ex_dir_rev_nonortho_nb

    dpsi_orb_is_zero = .true.

    do ex_dir_rev = 1, ex_dir_rev_nb

    ex_is_zero = .true.

!   swap orbitals for reverse excitation j->i
    if (ex_dir_rev == 2 .or. ex_dir_rev_nonortho == 2) then
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

!    write(6,'(2a,i4)') trim(here),': consider CSF # ',csf_i

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

        det_unq_cur_up = hash_table_get(up_dets, det_ex_cur_orb_lab_srt_up(:), success) !TA
        if (.not.(det_unq_cur_up.eq.0)) then
          if (det_unq_cur_up.le.ndetup) then
            det_ex_unq_sgn_cur_up = det_ex_cur_orb_lab_srt_sgn_up
          else
            det_ex_unq_sgn_cur_up = det_ex_cur_orb_lab_srt_sgn_up * det_ex_unq_orb_lab_srt_sgn_up (det_unq_cur_up - ndetup)
          endif
        endif

!       if current excited determinant is a new determinant, add it to the list of unique excited determinants
        if (det_unq_cur_up == 0) then
         det_ex_unq_up_nb = det_ex_unq_up_nb + 1
         det_unq_cur_up = ndetup + det_ex_unq_up_nb
         det_ex_unq_sgn_cur_up = 1
         det_ex_unq_ref_up (det_ex_unq_up_nb) = det_unq_up_i
         det_ex_unq_orb_lab_srt_up (:, det_ex_unq_up_nb) = det_ex_cur_orb_lab_srt_up (:)
         det_ex_unq_orb_lab_srt_sgn_up (det_ex_unq_up_nb) = det_ex_cur_orb_lab_srt_sgn_up
         call hash_table_add(up_dets, det_ex_cur_orb_lab_srt_up(:), det_unq_cur_up) !TA
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

        det_unq_cur_dn = hash_table_get(dn_dets, det_ex_cur_orb_lab_srt_dn(:), success) !TA
        if (.not.(det_unq_cur_dn.eq.0)) then
          if (det_unq_cur_dn.le.ndetdn) then
            det_ex_unq_sgn_cur_dn = det_ex_cur_orb_lab_srt_sgn_dn
          else
            det_ex_unq_sgn_cur_dn = det_ex_cur_orb_lab_srt_sgn_dn * det_ex_unq_orb_lab_srt_sgn_dn (det_unq_cur_dn - ndetdn)
          endif
        endif

!       if current excited determinant is a new determinant, add it to the list of unique excited determinants

        if (det_unq_cur_dn == 0) then
         det_ex_unq_dn_nb = det_ex_unq_dn_nb + 1
         det_unq_cur_dn = ndetdn + det_ex_unq_dn_nb
         det_ex_unq_sgn_cur_dn = 1
         det_ex_unq_ref_dn (det_ex_unq_dn_nb) = det_unq_dn_i
         det_ex_unq_orb_lab_srt_dn (:, det_ex_unq_dn_nb) = det_ex_cur_orb_lab_srt_dn (:)
         det_ex_unq_orb_lab_srt_sgn_dn (det_ex_unq_dn_nb) = det_ex_cur_orb_lab_srt_sgn_dn
         call hash_table_add(dn_dets, det_ex_cur_orb_lab_srt_dn(:), det_unq_cur_dn) !TA
!        write(6,'(2a,100i4)') trim(here), ': add new excited spin-dn determinant:', det_ex_cur_orb_lab_srt_dn (:)
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
!FIXED BUG - should be det_unq_dn_in_csf_ex_cur_temp - TA
!        if (det_unq_up_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i) == det_unq_up_in_csf_ex_cur (det_unq_in_csf_ex_cur_nb) .and.  &
!            det_unq_up_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i) == det_unq_dn_in_csf_ex_cur (det_unq_in_csf_ex_cur_nb)) then
        if (det_unq_up_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i) == det_unq_up_in_csf_ex_cur (det_unq_in_csf_ex_cur_nb) .and.  &
            det_unq_dn_in_csf_ex_cur_temp (det_unq_in_csf_ex_cur_i) == det_unq_dn_in_csf_ex_cur (det_unq_in_csf_ex_cur_nb)) then
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
!     write(6,'(2a)') trim(here), ': construct excited csf:'
!     write(6,'(2a,100i4)') trim(here),': det_unq_up_in_csf_ex_cur=',det_unq_up_in_csf_ex_cur (:)
!     write(6,'(2a,100i4)') trim(here),': det_unq_dn_in_csf_ex_cur=',det_unq_dn_in_csf_ex_cur (:)
!     write(6,'(2a,100f7.3)')  trim(here),': cdet_unq_in_csf_ex_cur=', cdet_unq_in_csf_ex_cur (:)

!     check if current excited csf is a csf in the ground-state wave function
      csf_cur = 0
      do csf_k = 1, ncsf
       if (arrays_equal (det_unq_up_in_csf_ex_cur (:), det_unq_up_in_csf (csf_k)%row (:)) .and. &
           arrays_equal (det_unq_dn_in_csf_ex_cur (:), det_unq_dn_in_csf (csf_k)%row (:)) .and. &
           arrays_equal (cdet_unq_in_csf_ex_cur (:), cdet_unq_in_csf (csf_k)%row (:))) then
           csf_cur = csf_k
!          write(6,'(2a)') trim(here),': excited csf is ground-state csf'
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
!           write(6,'(2a)') trim(here),': excited csf is excited csf already encountered'
         endif
       enddo ! csf_ex_unq_i
      endif
!    if current excited csf is a new excited csf, add it to the list of unique excited csf
     if (csf_cur == 0) then
      csf_ex_unq_nb = csf_ex_unq_nb + 1
      csf_cur = ncsf + csf_ex_unq_nb

      csf_ex_unq_ref (csf_ex_unq_nb) = csf_i

      det_unq_in_csf_ex_unq_nb (csf_ex_unq_nb) = det_unq_in_csf_ex_cur_nb
      call copy (det_unq_up_in_csf_ex_cur, det_unq_up_in_csf_ex_unq (csf_ex_unq_nb)%row)
      call copy (det_unq_dn_in_csf_ex_cur, det_unq_dn_in_csf_ex_unq (csf_ex_unq_nb)%row)
      call copy (cdet_unq_in_csf_ex_cur, cdet_unq_in_csf_ex_unq (csf_ex_unq_nb)%row)
!     write(6,'(2a)') trim(here),': add new excited csf'
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
!    write(6,'(2a)') trim(here), ': construct singly-excited wave function:'
!    write(6,'(2a,100i4)') trim(here),': csf_unq_in_wf_ex_cur=',csf_unq_in_wf_ex_cur (:)
!    write(6,'(2a,100i4)') trim(here),': csf_unq_ref_in_wf_ex_cur=',csf_unq_ref_in_wf_ex_cur (:)
!    write(6,'(2a,100f7.3)') trim(here),': csf_unq_prefac_in_wf_ex_cur=',csf_unq_prefac_in_wf_ex_cur (:)

!    check if current singly-excited wave function is a singly-excited wave function already encountered
!    (do not check for equalities between csf coefficients)
     ex_cur = 0
     do ex_i = 1, single_ex_nb
       if (arrays_equal (csf_unq_in_wf_ex_cur (:), csf_unq_in_wf_ex (ex_i)%row (:)) .and.            &
           arrays_equal (csf_unq_ref_in_wf_ex_cur (:), csf_unq_ref_in_wf_ex (ex_i)%row (:)) .and.    &
           arrays_equal (csf_unq_prefac_in_wf_ex_cur (:), csf_unq_prefac_in_wf_ex (ex_i)%row (:))) then
          ex_cur = ex_i
!         write(6,'(2a)') trim(here), ': singly-excited wave function is a singly-excited wave function already encountered'
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

!    write(6,'(2a)') trim(here), ': add new singly-excited wave function'
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

!   write(6,'(2a)') trim(here), ': construct orbital derivative:'
!   write(6,'(2a,100i4)') trim(here),': csf_unq_in_dpsi_orb_cur=',csf_unq_in_dpsi_orb_cur (:)
!   write(6,'(2a,100i4)') trim(here),': csf_unq_ref_in_dpsi_orb_cur=',csf_unq_ref_in_dpsi_orb_cur (:)
!   write(6,'(2a,100f7.3)') trim(here),': csf_unq_prefac_in_dpsi_orb_cur=',csf_unq_prefac_in_dpsi_orb_cur (:)

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
!       write(6, '(2a)') trim(here),': orbital derivative is redundant (zero)'
    endif

!   if current orbital derivative is a linear combinaison of optimized ground-state csfs, then
!   this orbital derivative is redundant.
!   (Even if not optimized, the first csf is also considered as part of the optimization space.
!   This is necessary so that for a CASSCF wave function all active->active excitations are redundant.
!   I am not sure if one should not include also the first csf is the optimization for this to be really correct,
!   as for instance Danish MCSCF people do).
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
!       write(6, '(2a)') trim(here),': orbital derivative is redundant with linear combinaison of optimized ground-state csfs'
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
!        write(6, '(2a)') trim(here),': orbital derivative is redundant with another orbital derivative'
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
!     write(6,'(2a)') trim(here), ': remove singly-excited wave function'
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

      if (ex_orb_ind_cur /= 0) then
        ex_orb_ind (param_orb_nb) = ex_orb_ind_cur
        ex_orb_ind_rev (param_orb_nb) = ex_orb_ind_rev_cur
      else

!       write(6,'(2a,i4,a,i4)') trim(here),': orbital excitation ', orb_opt_lab_j, ' -> ', orb_opt_lab_i
!       write(6,'(2a)') trim(here),': orbital derivative with vanishing direct excitation but not vanishing reverse excitation'
!       write(6,'(2a)') trim(here),': this case is not yet implemented'
!       write(6,'(2a,i4)') trim(here),': the index of reverse excitation is ',ex_orb_ind_rev_cur
!       call die (here)

!       if direct excitation i->j is zero, swap direct and reverse excitations
!       this simply corresponds to changing the sign of the derivative (so the optimized parameter will get the opposite sign)
!       in particular, this can happen if the orbitals are in an unexpected order, e.g.  i<j with i virtual and j occupied
!       Warning: I did not check this very carefully
        ex_orb_ind (param_orb_nb) = ex_orb_ind_rev_cur
        ex_orb_ind_rev (param_orb_nb) = ex_orb_ind_cur
      endif

!     write(6,'(2a)') trim(here),': add new orbital derivative'
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

    enddo ! ex_dir_rev_nonortho

!   reset label of first orbital
    orb_opt_lab_i = orb_opt_lab_j

   enddo ! orb_opt_j

  enddo ! orb_opt_i

  !resize to actual size
  call object_alloc ('det_ex_unq_ref_up', det_ex_unq_ref_up, det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_orb_lab_srt_up', det_ex_unq_orb_lab_srt_up, nup, det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_orb_lab_srt_sgn_up', det_ex_unq_orb_lab_srt_sgn_up, det_ex_unq_up_nb)

  call object_alloc ('det_ex_unq_ref_dn', det_ex_unq_ref_dn, det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_orb_lab_srt_dn', det_ex_unq_orb_lab_srt_dn, ndn, det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_orb_lab_srt_sgn_dn', det_ex_unq_orb_lab_srt_sgn_dn, det_ex_unq_dn_nb)

  call object_alloc ('csf_ex_unq_ref', csf_ex_unq_ref, csf_ex_unq_nb)
  call object_alloc ('det_unq_in_csf_ex_unq_nb', det_unq_in_csf_ex_unq_nb, csf_ex_unq_nb)
  call object_alloc ('det_unq_up_in_csf_ex_unq', det_unq_up_in_csf_ex_unq, csf_ex_unq_nb)
  call object_alloc ('det_unq_dn_in_csf_ex_unq', det_unq_dn_in_csf_ex_unq, csf_ex_unq_nb)
  call object_alloc ('cdet_unq_in_csf_ex_unq', cdet_unq_in_csf_ex_unq, csf_ex_unq_nb)

! determine (occupied or vitual) orbitals mixed in orbital optimization
  call object_alloc ('orb_mix_lab', orb_mix_lab, orb_tot_nb)
  orb_mix_lab (:) = .false.
  do ex_i = 1, single_ex_nb
   orb_mix_lab (ex_orb_1st_lab (ex_i)) = .true.
   orb_mix_lab (ex_orb_2nd_lab (ex_i)) = .true.
  enddo ! ex_i
  orb_mix_lab_nb = 0
  do orb_i = 1, orb_tot_nb
    if (orb_mix_lab (orb_i)) then
      orb_mix_lab_nb = orb_mix_lab_nb + 1
    endif
  enddo

  write(6,'(a,i10)') ' Number of single orbital excitations = ', single_ex_nb
  write(6,'(a,i10)') ' Number of orbital derivatives        = ', param_orb_nb
  write(6,'(a,i10)') ' Number of orbitals involved          = ', orb_mix_lab_nb

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

  if (l_print_orbital_excitations) then
  write(6,*)
  write(6,'(a)') ' Orbital excitations:'
   do dorb_i = 1, param_orb_nb
    write(6,'(a,i4,a,i4,a,i4,a,i4)') ' orbital parameter # ',dorb_i,' corresponds to single excitation # ',ex_orb_ind (dorb_i),' : ',ex_orb_1st_lab (ex_orb_ind (dorb_i)),' -> ', ex_orb_2nd_lab (ex_orb_ind (dorb_i))
    if (ex_orb_ind_rev (dorb_i) /= 0) then
     write(6,'(a,i4,a,i4,a,i4)') '                             and reverse single excitation # ',ex_orb_ind_rev (dorb_i),' : ',ex_orb_1st_lab (ex_orb_ind_rev (dorb_i)),' -> ', ex_orb_2nd_lab (ex_orb_ind_rev (dorb_i))
    endif
   enddo ! dorb_i
  write(6,*)
  endif

 end subroutine single_ex_wf_bld_2

! ==============================================================================
  subroutine single_ex_det_bld
! ------------------------------------------------------------------------------
! Description   : build single orbital excitation information for determinants
! Description   : to be replaced by single_ex_wf_bld_2
!
! Created       : J. Toulouse, 25 Oct 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  use hash_table_mod
  implicit none

! local
  integer ex_i, orb_1st, orb_2nd
  integer det_unq_up_i, det_unq_up_k, det_unq_dn_i, det_unq_dn_k
  integer :: det_ex_unq_up_max
  integer :: det_ex_unq_dn_max
  integer :: det_i
  type(bucket), allocatable :: up_dets(:), dn_dets(:), csfs(:)
  logical :: success

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

  det_ex_unq_up_max = single_ex_nb*ndetup
  det_ex_unq_dn_max = single_ex_nb*ndetdn

  call object_alloc ('iwdet_ex_ref_up', iwdet_ex_ref_up, det_ex_unq_up_max)
  call object_alloc ('det_ex_unq_orb_lab_up', det_ex_unq_orb_lab_up, nup, det_ex_unq_up_max)
  call object_alloc ('det_ex_unq_orb_lab_srt_up', det_ex_unq_orb_lab_srt_up, nup, det_ex_unq_up_max)
  call object_alloc ('det_ex_unq_orb_lab_srt_sgn_up', det_ex_unq_orb_lab_srt_sgn_up, det_ex_unq_up_max)
  call object_alloc ('det_ex_unq_up_orb_1st_pos', det_ex_unq_up_orb_1st_pos, det_ex_unq_up_max)
  call object_alloc ('det_ex_unq_up_orb_2nd_lab', det_ex_unq_up_orb_2nd_lab, det_ex_unq_up_max)

  call object_alloc ('iwdet_ex_ref_dn', iwdet_ex_ref_dn, det_ex_unq_dn_max)
  call object_alloc ('det_ex_unq_orb_lab_dn', det_ex_unq_orb_lab_dn, ndn, det_ex_unq_dn_max)
  call object_alloc ('det_ex_unq_orb_lab_srt_dn', det_ex_unq_orb_lab_srt_dn, ndn, det_ex_unq_dn_max)
  call object_alloc ('det_ex_unq_orb_lab_srt_sgn_dn', det_ex_unq_orb_lab_srt_sgn_dn, det_ex_unq_dn_max)
  call object_alloc ('det_ex_unq_dn_orb_1st_pos', det_ex_unq_dn_orb_1st_pos, det_ex_unq_dn_max)
  call object_alloc ('det_ex_unq_dn_orb_2nd_lab', det_ex_unq_dn_orb_2nd_lab, det_ex_unq_dn_max)

  allocate(up_dets(det_ex_unq_up_max+ndetup))
  do det_i=1, ndetup
    call hash_table_add(up_dets, det_unq_orb_lab_srt_up(:,det_i), det_i)
  enddo

  allocate(dn_dets(det_ex_unq_dn_max+ndetdn))
  do det_i=1, ndetdn
    call hash_table_add(dn_dets, det_unq_orb_lab_srt_dn(:,det_i), det_i)
  enddo

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

       det_unq_up_k = hash_table_get(up_dets, det_ex_orb_lab_srt_up (:, ex_i, det_unq_up_i), success) !TA
       if (.not.(det_unq_up_k.eq.0)) then
         if (det_unq_up_k.le.ndetup) then
           is_det_ex_up (ex_i, det_unq_up_i) = .false.
           iwdet_ex_up (ex_i, det_unq_up_i) = det_unq_up_k
         else
           iwdet_ex_up (ex_i, det_unq_up_i) = det_unq_up_k - ndetup
         endif
       endif

!      if current excited determinant is a new determinant, add it to the list of unique excited determinants
       if (iwdet_ex_up (ex_i, det_unq_up_i) == 0) then
         det_ex_unq_up_nb = det_ex_unq_up_nb + 1
         iwdet_ex_up (ex_i, det_unq_up_i) = det_ex_unq_up_nb
         iwdet_ex_ref_up (det_ex_unq_up_nb) = det_unq_up_i
         det_ex_unq_orb_lab_up (:, det_ex_unq_up_nb) = det_ex_orb_lab_up (:, ex_i, det_unq_up_i)
         det_ex_unq_orb_lab_srt_up (:, det_ex_unq_up_nb) = det_ex_orb_lab_srt_up (:, ex_i, det_unq_up_i)
         det_ex_unq_orb_lab_srt_sgn_up (det_ex_unq_up_nb) = det_ex_orb_lab_srt_sgn_up (ex_i, det_unq_up_i)
         det_ex_unq_up_orb_1st_pos (det_ex_unq_up_nb) = orb_pos_in_det_unq_up (orb_1st, det_unq_up_i)
         det_ex_unq_up_orb_2nd_lab (det_ex_unq_up_nb) = orb_2nd
         call hash_table_add(up_dets, det_ex_orb_lab_srt_up(:, ex_i, det_unq_up_i), ndetup + det_ex_unq_up_nb)
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

       det_unq_dn_k = hash_table_get(dn_dets, det_ex_orb_lab_srt_dn (:, ex_i, det_unq_dn_i), success) !TA
       if (.not.(det_unq_dn_k.eq.0)) then
         if (det_unq_dn_k.le.ndetdn) then
           is_det_ex_dn (ex_i, det_unq_dn_i) = .false.
           iwdet_ex_dn (ex_i, det_unq_dn_i) = det_unq_dn_k
         else
           iwdet_ex_dn (ex_i, det_unq_dn_i) = det_unq_dn_k - ndetdn
         endif
       endif

!      if current excited determinant is a new determinant, add it to the list of unique excited determinants
       if (iwdet_ex_dn (ex_i, det_unq_dn_i) == 0) then
         det_ex_unq_dn_nb = det_ex_unq_dn_nb + 1
         iwdet_ex_dn (ex_i, det_unq_dn_i) = det_ex_unq_dn_nb
         iwdet_ex_ref_dn (det_ex_unq_dn_nb) = det_unq_dn_i
         det_ex_unq_orb_lab_dn (:, det_ex_unq_dn_nb) = det_ex_orb_lab_dn (:, ex_i, det_unq_dn_i)
         det_ex_unq_orb_lab_srt_dn (:, det_ex_unq_dn_nb) = det_ex_orb_lab_srt_dn (:, ex_i, det_unq_dn_i)
         det_ex_unq_orb_lab_srt_sgn_dn (det_ex_unq_dn_nb) = det_ex_orb_lab_srt_sgn_dn (ex_i, det_unq_dn_i)
         det_ex_unq_dn_orb_1st_pos (det_ex_unq_dn_nb) = orb_pos_in_det_unq_dn (orb_1st, det_unq_dn_i)
         det_ex_unq_dn_orb_2nd_lab (det_ex_unq_dn_nb) = orb_2nd
         call hash_table_add(dn_dets, det_ex_orb_lab_srt_dn(:, ex_i, det_unq_dn_i), ndetdn + det_ex_unq_dn_nb)
       endif

     endif

     enddo ! det_unq_dn_i

  enddo ! ex_i

! resize to actual size
  call object_alloc ('iwdet_ex_ref_up', iwdet_ex_ref_up, det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_orb_lab_up', det_ex_unq_orb_lab_up, nup, det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_orb_lab_srt_up', det_ex_unq_orb_lab_srt_up, nup, det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_orb_lab_srt_sgn_up', det_ex_unq_orb_lab_srt_sgn_up, det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_up_orb_1st_pos', det_ex_unq_up_orb_1st_pos, det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_up_orb_2nd_lab', det_ex_unq_up_orb_2nd_lab, det_ex_unq_up_nb)

  call object_alloc ('iwdet_ex_ref_dn', iwdet_ex_ref_dn, det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_orb_lab_dn', det_ex_unq_orb_lab_dn, ndn, det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_orb_lab_srt_dn', det_ex_unq_orb_lab_srt_dn, ndn, det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_orb_lab_srt_sgn_dn', det_ex_unq_orb_lab_srt_sgn_dn, det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_dn_orb_1st_pos', det_ex_unq_dn_orb_1st_pos, det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_dn_orb_2nd_lab', det_ex_unq_dn_orb_2nd_lab, det_ex_unq_dn_nb)

! Warning: tmp debug printout
! do ex_i=1,single_ex_nb
!   do det_unq_up_i=1,ndetup
!     write(6,'(''ex_i, det_unq_up_i, iwdet_ex_up(ex_i,det_unq_up_i)'',3i5)') ex_i, det_unq_up_i, iwdet_ex_up(ex_i,det_unq_up_i)
!   enddo
! enddo

! do ex_i=1,single_ex_nb
!   do det_unq_dn_i=1,ndetdn
!     write(6,'(''ex_i, det_unq_dn_i, iwdet_ex_dn(ex_i,det_unq_dn_i)'',3i5)') ex_i, det_unq_dn_i, iwdet_ex_dn(ex_i,det_unq_dn_i)
!   enddo
! enddo

! write(6,'(''det_ex_unq_up_nb='',i5)') det_ex_unq_up_nb
! write(6,'(''iwdet_ex_ref_up'',10000i5)') iwdet_ex_ref_up
! write(6,'(''iwdet_ex_ref_dn'',10000i5)') iwdet_ex_ref_dn

! write(6,'(''det_ex_unq_orb_lab_up'',10000i5)') det_ex_unq_orb_lab_up
! write(6,'(''det_ex_unq_orb_lab_dn'',10000i5)') det_ex_unq_orb_lab_dn

! write(6,'(''det_ex_unq_up_orb_1st_pos'',10000i5)') det_ex_unq_up_orb_1st_pos
! write(6,'(''det_ex_unq_up_orb_2nd_lab'',10000i5)') det_ex_unq_up_orb_2nd_lab
! write(6,'(''det_ex_unq_dn_orb_1st_pos'',10000i5)') det_ex_unq_dn_orb_1st_pos
! write(6,'(''det_ex_unq_dn_orb_2nd_lab'',10000i5)') det_ex_unq_dn_orb_2nd_lab

! write(6,'(''param_orb_nb='',i5)') param_orb_nb
! write(6,'(''ex_orb_ind'',10000i5)') ex_orb_ind
! write(6,'(''ex_orb_ind_rev'',10000i5)') ex_orb_ind_rev

! write(6,'(''single_ex_nb='',i5)') single_ex_nb
! write(6,'(''ex_orb_1st_lab'',10000i5)') ex_orb_1st_lab
! write(6,'(''ex_orb_2nd_lab'',10000i5)') ex_orb_2nd_lab
! End tmp debug printout

  write(6,'(a,i10)') ' Number of unique spin-up   excited determinants = ',det_ex_unq_up_nb
  write(6,'(a,i10)') ' Number of unique spin-down excited determinants = ',det_ex_unq_dn_nb

  end subroutine single_ex_det_bld

!==============================================================================
  subroutine dpsi_orb_bld
!------------------------------------------------------------------------------
!Description   : Calculate derivative of determinantal psi with respect to 
!Description   : orbitals.
!
!Created       : Tyler Anderson, Fall 2020
!------------------------------------------------------------------------------
    implicit none

    integer     :: iparam, i, j, k, iup, idn, jup, jdn
    real(dp)    :: dpsi, diff

    if (header_exe) then
      call object_create ('dpsi_orb')
      call object_needed ('param_orb_nb')
      call object_needed ('ex_orb_ind')
      call object_needed ('ex_orb_ind_rev')
      call object_needed ('ex_orb_1st_lab')
      call object_needed ('ex_orb_2nd_lab')
      call object_needed ('gup')
      call object_needed ('gdn')
      call object_needed ('atup')
      call object_needed ('atdn')
      return
    endif

    call object_alloc ('dpsi_orb', dpsi_orb, param_orb_nb)

    do iparam=1, param_orb_nb
      i = ex_orb_1st_lab(ex_orb_ind(iparam))
      j = ex_orb_2nd_lab(ex_orb_ind(iparam))
      
      dpsi = 0

      iup = ioccup(i)
      if (iup.gt.0) then
        do k=1,nup
          dpsi = dpsi + gup(iup,k)*atup(k,j)
        enddo
      endif

      idn = ioccdn(i)
      if (idn.gt.0) then
        do k=1,ndn
          dpsi = dpsi + gdn(idn,k)*atdn(k,j)
        enddo
      endif

      if ((.not. l_casscf) .and. (ex_orb_ind_rev(iparam) /= 0)) then

      jup = ioccup(j)
      if (jup.gt.0) then
        do k=1,nup
          dpsi = dpsi - gup(jup,k)*atup(k,i)
        enddo
      endif

      jdn = ioccdn(j)
      if (jdn.gt.0) then
        do k=1,ndn
          dpsi = dpsi - gdn(jdn,k)*atdn(k,i)
        enddo
      endif

      endif

      dpsi_orb(iparam) = dpsi
    enddo
  end subroutine

!==============================================================================
  subroutine deloc_orb_bld
!------------------------------------------------------------------------------
!Description   : Calculate derivative of local energy with respect to 
!Description   : orbitals.
!
!Created       : Tyler Anderson, Fall 2020
!------------------------------------------------------------------------------

    implicit none

    real(dp)    :: deloc
    integer     :: i, j, k, iparam, iup, jup, idn, jdn

    if (header_exe) then
      call object_create ('deloc_orb')
      call object_needed ('dgup')
      call object_needed ('dgdn')
      call object_needed ('gup')
      call object_needed ('gdn')
      call object_needed ('atup')
      call object_needed ('atdn')
      call object_needed ('elocal_orb_up')
      call object_needed ('elocal_orb_dn')
      return
    endif

    call object_alloc ('deloc_orb', deloc_orb, param_orb_nb)

    do iparam=1, param_orb_nb
      i = ex_orb_1st_lab(ex_orb_ind(iparam))
      j = ex_orb_2nd_lab(ex_orb_ind(iparam))

      deloc = 0

      iup = ioccup(i)
      if (iup.gt.0) then
        do k=1,nup
          deloc = &
          deloc + gup(iup,k)*elocal_orb_up(k,j) + dgup(iup,k)*atup(k,j)
        enddo
      endif

      idn = ioccdn(i)
      if (idn.gt.0) then
        do k=1,ndn
          deloc = &
          deloc + gdn(idn,k)*elocal_orb_dn(k,j) + dgdn(idn,k)*atdn(k,j)
        enddo
      endif

      if ((.not. l_casscf) .and. (ex_orb_ind_rev(iparam) /= 0)) then

      jup = ioccup(j)
      if (jup.gt.0) then
        do k=1,nup
          deloc = &
          deloc - gup(jup,k)*elocal_orb_up(k,i) - dgup(jup,k)*atup(k,i)
        enddo
      endif

      jdn = ioccdn(j)
      if (jdn.gt.0) then
        do k=1,ndn
          deloc = &
          deloc - gdn(jdn,k)*elocal_orb_dn(k,i) - dgdn(jdn,k)*atdn(k,i)
        enddo
      endif

      endif

      deloc_orb(iparam) = deloc
    enddo
  end subroutine

! ==============================================================================
  subroutine delta_eps_bld
! ------------------------------------------------------------------------------
! Description  : orbital energy difference for single excitations
!
! Created      : J. Toulouse, 04 Nov 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

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
