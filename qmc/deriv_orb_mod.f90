module deriv_orb_mod

  use all_tools_mod
  use montecarlo_mod
  use orbitals_mod
  use determinants_mod
  use csfs_mod
  use eloc_mod

! Declaration of global variables and default values

! single_ex_wf_bld_check 
  integer                                :: param_orb_nb = 0
  integer                                :: single_ex_nb = 0
  integer                                :: deriv_orb_pairs_nb = 0
  integer, allocatable                   :: ex_orb_ind (:)
  integer, allocatable                   :: ex_orb_ind_rev (:)
  integer, allocatable                   :: ex_orb_1st_lab (:)
  integer, allocatable                   :: ex_orb_2nd_lab (:)
  integer                                :: orb_mix_lab_nb
  logical, allocatable                   :: orb_mix_lab (:)

  integer                                :: det_ex_unq_up_nb
  integer                                :: det_ex_unq_dn_nb
  integer, allocatable                   :: det_ex_unq_up_orb_1st_pos (:)
  integer, allocatable                   :: det_ex_unq_dn_orb_1st_pos (:)
  integer, allocatable                   :: det_ex_unq_up_orb_2nd_lab (:)
  integer, allocatable                   :: det_ex_unq_dn_orb_2nd_lab (:)
  integer, allocatable                   :: det_ex_unq_orb_lab_up (:,:)
  integer, allocatable                   :: det_ex_unq_orb_lab_dn (:,:)
  integer, allocatable                   :: det_ex_unq_sgn_up (:,:)
  integer, allocatable                   :: det_ex_unq_sgn_dn (:,:)
  integer, allocatable                   :: iwdet_ex_up (:,:)
  integer, allocatable                   :: iwdet_ex_dn (:,:)
  integer, allocatable                   :: iwdet_ex_ref_up (:)
  integer, allocatable                   :: iwdet_ex_ref_dn (:)

! det_ex_unq_bld
  real(dp), allocatable                  :: det_ex_unq_up (:)
  real(dp), allocatable                  :: det_ex_unq_dn (:)

! det_ex_bld
  real(dp), allocatable                  :: det_ex_up (:,:)
  real(dp), allocatable                  :: det_ex_dn (:,:)
  real(dp), allocatable                  :: det_ex (:,:)

! dpsi_orb_bld
  real(dp), allocatable                  :: psid_ex (:)
  real(dp), allocatable                  :: dpsi_orb(:)
  real(dp), allocatable                  :: dcsf_orb(:,:)

! slater_mat_ex_trans_inv_bld
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_up (:, :, :)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_dn (:, :, :)
  logical                                :: l_slater_mat_ex_trans_inv_sm = .true.

! slater_mat_ex_trans_inv_bld_2
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_up_2 (:,:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_inv_dn_2 (:,:,:)

! grd
  real(dp), allocatable                  :: grd_det_ex_unq_up (:,:,:)
  real(dp), allocatable                  :: grd_det_ex_unq_dn (:,:,:)
  real(dp), allocatable                  :: grd_det_ex_dn (:,:,:,:)
  real(dp), allocatable                  :: grd_det_ex_up (:,:,:,:)
  real(dp), allocatable                  :: grd_psid_ex_over_psid (:,:,:)
  real(dp), allocatable                  :: grd_psi_ex_over_psi (:,:,:)

! lap
  real(dp), allocatable                  :: lap_det_ex_unq_up (:,:)
  real(dp), allocatable                  :: lap_det_ex_unq_dn (:,:)
  real(dp), allocatable                  :: lap_det_ex_dn (:,:,:)
  real(dp), allocatable                  :: lap_det_ex_up (:,:,:)
  real(dp), allocatable                  :: lap_psid_ex_over_psid (:,:)
  real(dp), allocatable                  :: lap_lnpsid_ex (:,:)
  real(dp), allocatable                  :: sum_lap_lnpsid_ex (:)
  real(dp), allocatable                  :: sum_lap_lnpsi_ex (:)

! eloc
  real(dp), allocatable                  :: eloc_kin_ex (:)
  real(dp), allocatable                  :: eloc_pot_nloc_ex (:)
  real(dp), allocatable                  :: eloc_pot_ex (:)
  real(dp), allocatable                  :: eloc_ex (:)
  real(dp), allocatable                  :: vpot_ex (:,:)
  real(dp), allocatable                  :: vpsp_ex (:)

! psid_ex_in_x_bld
  real(dp), allocatable                  :: psid_ex_in_x (:)
  integer                                :: electron

! deloc_orb_bld
  real(dp), allocatable                  :: deloc_orb (:)

! delta_eps_bld
  real(dp), allocatable                  :: delta_eps (:)

! double_ex_det_bld
  integer                                :: double_ex_nb
  integer                                :: det_ex2_unq_up_nb
  integer                                :: det_ex2_unq_dn_nb
  integer, allocatable                   :: det_ex2_unq_up_orb_info (:,:)
  integer, allocatable                   :: det_ex2_unq_dn_orb_info (:,:)

! det_ex2_unq_bld
  real(dp), allocatable                  :: det_ex2_unq_up (:)
  real(dp), allocatable                  :: det_ex2_unq_dn (:)

! det_ex2_bld
  real(dp), allocatable                  :: det_ex2_up (:,:)
  real(dp), allocatable                  :: det_ex2_dn (:,:)
  real(dp), allocatable                  :: det_ex2 (:,:)

! d2psi_orb_bld
  real(dp), allocatable                  :: psid_ex2 (:)
  real(dp), allocatable                  :: d2psi_orb(:)
  real(dp), allocatable                  :: d2csf_orb(:,:)

  contains

! ==============================================================================
  subroutine single_ex_wf_bld
! ------------------------------------------------------------------------------
! Description   : This routines finds single orbital excitations
! Description   : and unique orbital derivatives.
! Description   :
! Description   : For each pair of orbitals, if the excitation is trivially allowed
! Description   : (in terms of the classes of both orbitals, symmetry, etc...)
! Description   : -  the unique determinants spawned by the excitation are constructed (not calculated)
! Description   : -  the unique CSFs are constructed
! Description   :    (an excitation can spawn excited determinant that do not form a new excited CSF)
! Description   : -  the unique WF is constructed
! Description   :    (similarly, an excitation that spawns excited CSF can form an excited WF that is not new)
! Description   :
! Description   : A unique orbital derivative is a derivative that, 
! Description   : when applied to a CI expansion, yields a unique excited wavefunction.
! Description   :
! Description   : Information on all unique excited determinants is gathered for use outside of the routine.
!
! Created       : B. Mussard, July 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'single_ex_wf_bld_check'
  integer :: i, j, ireinit
  integer :: orb_opt_lab_i, orb_opt_lab_j
  integer :: ex_i, ex_j
  integer :: csf_i, det_in_csf_i, det_i, orb_i, dpsi_i
  integer :: csf_k, dpsi_k
  integer :: this_csf, this_ex
  integer :: det_unq_up_i, det_unq_dn_i
  integer :: det_unq_up_k, det_unq_dn_k
  integer :: iwdet
  integer :: ex_dir_rev_nonortho, ex_dir_rev_nonortho_nb
  integer :: ex_dir_rev, ex_dir_rev_nb
  integer :: ex_orb_ind_tmp,ex_orb_ind_rev_tmp
  logical :: dpsi_is_redundant
  integer :: single_ex_added_nb
  real(dp) :: tA, tB, tAintermed, tBintermed, tSort, tDet, tCSFa,tCSFb,tCSFc, tWF, tDPSI, tBegin, tEnd
  integer :: single_ex_max,det_ex_max,csf_ex_unq_max,det_in_csf_max,csf_in_wf_max,csf_in_dpsi_max,wf_unq_max,dpsi_unq_max

! Unique DET (most object are global, only the following are local)
  integer, allocatable :: det_orb_lab_up (:)
  integer, allocatable :: det_orb_lab_dn (:)
  integer, allocatable :: det_orb_lab_srt_up (:)
  integer, allocatable :: det_orb_lab_srt_dn (:)
  integer, allocatable :: det_ex_to_calc_sgn_up (:)
  integer, allocatable :: det_ex_to_calc_sgn_dn (:)

! Unique CSF
!   an "CSF=\sum_k^n c_k Dkup Dkdn" is the data of
!    - the number of determinant "n",
!    - the "Dkup" UP determinant,
!    - the "Dkdn" DOWN determinants
!    - the coefficients "c_k"
! We need an unsorted and a sorted version, as well as the unique data:
  integer              :: det_in_csf_nb_tmp
  integer, allocatable :: det_in_csf_up_tmp(:)
  integer, allocatable :: det_in_csf_dn_tmp(:)
  real(8), allocatable :: det_in_csf_coef_tmp(:)
  integer              :: det_in_csf_nb
  integer, allocatable :: det_in_csf_up(:)
  integer, allocatable :: det_in_csf_dn(:)
  real(8), allocatable :: det_in_csf_coef(:)

  integer              :: csf_ex_unq_nb
  integer, allocatable :: csf_ex_unq_ref(:)
  integer, allocatable :: csf_ex_unq_size(:)

! Unique WF
!   a "WF=\sum_I^n c_I CSF_I" is the data of
!    - the number of CSF "n",
!    - the "CSF_I"
!    - the reference of these CSFs,
!    - the coefficients "c_I"
! We need the whole data, as well as the unique data:
  integer              :: csf_in_wf_nb
  integer, allocatable :: csf_in_wf(:)
  integer, allocatable :: csf_in_wf_ref(:)
  real(8), allocatable :: csf_in_wf_coef(:)
  integer,                allocatable :: wf_unq_nb(:)
  type(type_integer_row), allocatable :: wf_unq_csf(:)
  type(type_integer_row), allocatable :: wf_unq_ref(:)
  type(type_real_row)   , allocatable :: wf_unq_coef(:)

! Unique DPSI
!   a derivative DPSI
!   contains the excited WF and the reverse excited WF
!   so : it is determined by the same data as the WF
! We need an unsorted and a sorted version, as well as the unique data:
  integer              :: csf_in_dpsi_nb_tmp
  integer, allocatable :: csf_in_dpsi_tmp(:)
  integer, allocatable :: csf_in_dpsi_ref_tmp(:)
  real(8), allocatable :: csf_in_dpsi_coef_tmp(:)
  integer              :: csf_in_dpsi_nb
  integer, allocatable :: csf_in_dpsi(:)
  integer, allocatable :: csf_in_dpsi_ref(:)
  real(8), allocatable :: csf_in_dpsi_coef(:)

  integer, allocatable                   :: dpsi_unq_nb (:)
  type (type_integer_row), allocatable   :: dpsi_unq_csf (:)
  type (type_integer_row), allocatable   :: dpsi_unq_ref (:)
  type (type_real_row),    allocatable   :: dpsi_unq_coef (:)


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

   call object_create ('det_ex_unq_up_nb')
   call object_create ('det_ex_unq_dn_nb')
   call object_create ('det_ex_unq_up_orb_1st_pos')
   call object_create ('det_ex_unq_dn_orb_1st_pos')
   call object_create ('det_ex_unq_up_orb_2nd_lab')
   call object_create ('det_ex_unq_dn_orb_2nd_lab')
   call object_create ('det_ex_unq_orb_lab_up')
   call object_create ('det_ex_unq_orb_lab_dn')

   call object_create ('det_ex_unq_sgn_up')
   call object_create ('det_ex_unq_sgn_dn')
   call object_create ('iwdet_ex_up')
   call object_create ('iwdet_ex_dn')
   call object_create ('iwdet_ex_ref_up')
   call object_create ('iwdet_ex_ref_dn')

   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('det_unq_up_in_csf')
   call object_needed ('det_unq_dn_in_csf')
   call object_needed ('cdet_unq_in_csf')

   call object_needed ('orb_tot_nb')
   call object_needed ('orb_opt_lab')
   call object_needed ('orb_ex_forbidden')
   call object_needed ('orb_sym_lab')
   call object_needed ('orb_cls_in_wf')
   call object_needed ('orb_vir_in_wf')
   call object_needed ('orb_opn_in_wf')
   call object_needed ('orb_act_in_wf')
   call object_needed ('orb_occ_in_wf')
   call object_needed ('orb_occ_in_det_unq_up')
   call object_needed ('orb_occ_in_det_unq_dn')

   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('cdet_in_csf')
   call object_needed ('iwcsf')

   return
  endif

  write(6,'(a)') ' Constructing list of single orbital excitations checking for redundancies...'

! allocations
! Most arrays will be adjusted at the end of the routine and have the following max sizes:
  single_ex_max  =orb_tot_nb*orb_tot_nb
  det_ex_max     =(ndetup+ndetdn)*single_ex_max
  csf_ex_unq_max =2*det_ex_max*ndetup
  det_in_csf_max =50
  csf_in_wf_max  =ncsf
  csf_in_dpsi_max=2*ncsf
  wf_unq_max     =single_ex_max
  dpsi_unq_max   =single_ex_max
! This suppress the need to resize an array to often.

! FOR PART 1: UNIQUE DETS
  det_ex_unq_up_nb = 0
  det_ex_unq_dn_nb = 0
  call object_alloc ('iwdet_ex_up',               iwdet_ex_up,                 single_ex_max, ndetup)
  call object_alloc ('iwdet_ex_dn',               iwdet_ex_dn,                 single_ex_max, ndetdn)
  call object_alloc ('det_ex_unq_sgn_up',         det_ex_unq_sgn_up,           single_ex_max, ndetup)
  call object_alloc ('det_ex_unq_sgn_dn',         det_ex_unq_sgn_dn,           single_ex_max, ndetdn)
  call object_alloc ('iwdet_ex_ref_up',           iwdet_ex_ref_up,             det_ex_max)
  call object_alloc ('iwdet_ex_ref_dn',           iwdet_ex_ref_dn,             det_ex_max)
  call object_alloc ('det_ex_unq_orb_lab_up',     det_ex_unq_orb_lab_up, nup,  det_ex_max)
  call object_alloc ('det_ex_unq_orb_lab_dn',     det_ex_unq_orb_lab_dn, ndn,  det_ex_max)
  call object_alloc ('det_unq_orb_lab_srt_up',    det_unq_orb_lab_srt_up, nup, ndetup+det_ex_max)
  call object_alloc ('det_unq_orb_lab_srt_dn',    det_unq_orb_lab_srt_dn, ndn, ndetdn+det_ex_max)
  call object_alloc ('det_ex_unq_up_orb_1st_pos', det_ex_unq_up_orb_1st_pos,   det_ex_max)
  call object_alloc ('det_ex_unq_dn_orb_1st_pos', det_ex_unq_dn_orb_1st_pos,   det_ex_max)
  call object_alloc ('det_ex_unq_up_orb_2nd_lab', det_ex_unq_up_orb_2nd_lab,   det_ex_max)
  call object_alloc ('det_ex_unq_dn_orb_2nd_lab', det_ex_unq_dn_orb_2nd_lab,   det_ex_max)
! local arrays (not objects)
  call alloc ('det_orb_lab_up',     det_orb_lab_up,     nup)
  call alloc ('det_orb_lab_dn',     det_orb_lab_dn,     ndn)
  call alloc ('det_orb_lab_srt_up', det_orb_lab_srt_up, nup)
  call alloc ('det_orb_lab_srt_dn', det_orb_lab_srt_dn, ndn)
  call alloc ('det_ex_to_calc_sgn_up',  det_ex_to_calc_sgn_up,  det_ex_max)
  call alloc ('det_ex_to_calc_sgn_dn',  det_ex_to_calc_sgn_dn,  det_ex_max)

! FOR PART 2: UNIQUE CSFS
  det_in_csf_nb_tmp=0
  call alloc('det_in_csf_up_tmp',   det_in_csf_up_tmp,   det_in_csf_max)
  call alloc('det_in_csf_dn_tmp',   det_in_csf_dn_tmp,   det_in_csf_max)
  call alloc('det_in_csf_coef_tmp', det_in_csf_coef_tmp, det_in_csf_max)
  det_in_csf_nb=0
  call alloc('det_in_csf_up',       det_in_csf_up,       det_in_csf_max)
  call alloc('det_in_csf_dn',       det_in_csf_dn,       det_in_csf_max)
  call alloc('det_in_csf_coef',     det_in_csf_coef,     det_in_csf_max)

  csf_ex_unq_nb=0
  call alloc('csf_ex_unq_ref',           csf_ex_unq_ref,    csf_ex_unq_max)
  call alloc('csf_ex_unq_size',          csf_ex_unq_size,   csf_ex_unq_max)
  call object_alloc('det_unq_up_in_csf', det_unq_up_in_csf, ncsf+csf_ex_unq_max)
  call object_alloc('det_unq_dn_in_csf', det_unq_dn_in_csf, ncsf+csf_ex_unq_max)
  call object_alloc('cdet_unq_in_csf'  , cdet_unq_in_csf,   ncsf+csf_ex_unq_max)

! FOR PART 3: UNIQUE WFS
  csf_in_wf_nb=0
  call alloc('csf_in_wf',      csf_in_wf,      csf_in_wf_max)
  call alloc('csf_in_wf_ref',  csf_in_wf_ref,  csf_in_wf_max)
  call alloc('csf_in_wf_coef', csf_in_wf_coef, csf_in_wf_max)

  call alloc('wf_unq_nb',      wf_unq_nb,      wf_unq_max)
  call alloc('wf_unq_csf',     wf_unq_csf,     wf_unq_max)
  call alloc('wf_unq_ref',     wf_unq_ref,     wf_unq_max)
  call alloc('wf_unq_coef',    wf_unq_coef,    wf_unq_max)

! FOR PART 4: UNIQUE DPSI
  csf_in_dpsi_nb_tmp=0
  call alloc('csf_in_dpsi_tmp',      csf_in_dpsi_tmp,      csf_in_dpsi_max)
  call alloc('csf_in_dpsi_ref_tmp',  csf_in_dpsi_ref_tmp,  csf_in_dpsi_max)
  call alloc('csf_in_dpsi_coef_tmp', csf_in_dpsi_coef_tmp, csf_in_dpsi_max)
  csf_in_dpsi_nb=0
  call alloc('csf_in_dpsi',          csf_in_dpsi,          csf_in_dpsi_max)
  call alloc('csf_in_dpsi_ref',      csf_in_dpsi_ref,      csf_in_dpsi_max)
  call alloc('csf_in_dpsi_coef',     csf_in_dpsi_coef,     csf_in_dpsi_max)

  call alloc('dpsi_unq_nb',          dpsi_unq_nb,          dpsi_unq_max)
  call alloc('dpsi_unq_csf',         dpsi_unq_csf,         dpsi_unq_max)
  call alloc('dpsi_unq_ref',         dpsi_unq_ref,         dpsi_unq_max)
  call alloc('dpsi_unq_coef',        dpsi_unq_coef,        dpsi_unq_max)

! RELEVANT INFORMATION
  single_ex_nb     = 0
  param_orb_nb     = 0
  call object_alloc ('ex_orb_ind'    , ex_orb_ind    , single_ex_max)
  call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, single_ex_max)
  call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_max)
  call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_max)

! two scenarii
  if (l_active_orb_ortho_constraint) then
    ex_dir_rev_nb          = 2
    ex_dir_rev_nonortho_nb = 1
  else
    ex_dir_rev_nb          = 1
    ex_dir_rev_nonortho_nb = 2
  endif

  tSort=0.d0
  tDet=0.d0
  tCSFa=0.d0
  tCSFb=0.d0
  tCSFc=0.d0
  tWF=0.d0
  tDPSI=0.d0
  call cpu_time(tA)
  tBegin=tA

! loop over all pairs of orbitals
  ex_i=0
  do i=1, orb_tot_nb
  orb_opt_lab_i=orb_opt_lab(i)
  do j=i+1, orb_tot_nb
  orb_opt_lab_j=orb_opt_lab(j)

! loop over the direct  excitation i->j
!       and the reverse excitation j->i
! (one loop or the other is used depending
!  whether the orthogonality constraint is imposed or not)
  do ex_dir_rev_nonortho = 1, ex_dir_rev_nonortho_nb
  single_ex_added_nb = 0
  ex_orb_ind_tmp     = 0
  ex_orb_ind_rev_tmp = 0
  !DPSI_tmp initialization
  csf_in_dpsi_nb_tmp     = 0
  csf_in_dpsi_tmp        = csf_in_dpsi_max+1
  csf_in_dpsi_ref_tmp    = csf_in_dpsi_max+1
  csf_in_dpsi_coef_tmp   = csf_in_dpsi_max+1
  do ex_dir_rev = 1, ex_dir_rev_nb

!   swap orbitals for reverse excitation j- >i
    if (ex_dir_rev ==  2 .or. ex_dir_rev_nonortho ==  2) call swap (orb_opt_lab_i, orb_opt_lab_j)

!   skip for numerous reasons
    if ((orb_ex_forbidden (orb_opt_lab_i, orb_opt_lab_j)) &
    .or. (orb_sym_lab (orb_opt_lab_i) /= orb_sym_lab (orb_opt_lab_j)) &
    .or. (orb_cls_in_wf (orb_opt_lab_i) .and. orb_cls_in_wf (orb_opt_lab_j)) &
    .or. (orb_vir_in_wf (orb_opt_lab_i) .and. orb_vir_in_wf (orb_opt_lab_j)) &
    .or. (ndet == 1 .and. (orb_opn_in_wf (orb_opt_lab_i) .and. orb_opn_in_wf (orb_opt_lab_j))) &
    .or. (l_casscf .and. (orb_act_in_wf (orb_opt_lab_i) .and. orb_act_in_wf (orb_opt_lab_j))) &
    .or. (.not. orb_occ_in_wf (orb_opt_lab_i)) &
    .or. (orb_cls_in_wf (orb_opt_lab_j))) then
      cycle
    endif

    ex_i=ex_i+1
    if (ex_i>single_ex_max) call die(lhere,"Recompile with a superior value of 'single_ex_max' in deriv_orb_mod.f90 (A)")

!1\ ================================================================================= 
!   Gather information on unique excited determinants
!   ================================================================================= 
!1a\loop over unique spin-up determinants
    do det_unq_up_i = 1, ndetup
    if (orb_occ_in_det_unq_up (orb_opt_lab_i, det_unq_up_i) &
      .and. .not. orb_occ_in_det_unq_up (orb_opt_lab_j, det_unq_up_i)) then

!     build current excited determinant
      det_orb_lab_up = det_unq_orb_lab_srt_up (:, det_unq_up_i)
      call replace_elt_in_array (det_orb_lab_up, orb_opt_lab_i, orb_opt_lab_j)
      det_orb_lab_srt_up = det_orb_lab_up
      call sort_and_sign (det_orb_lab_srt_up, det_ex_unq_sgn_up (ex_i, det_unq_up_i))

!     check if current excited determinant is an already known determinant
      do det_unq_up_k = 1, ndetup+det_ex_unq_up_nb
       if (arrays_equal (det_orb_lab_srt_up, det_unq_orb_lab_srt_up (:, det_unq_up_k))) then
         iwdet_ex_up (ex_i, det_unq_up_i) = det_unq_up_k
         if (det_unq_up_k.gt.ndetup) then
           det_ex_unq_sgn_up (ex_i, det_unq_up_i)=det_ex_unq_sgn_up (ex_i, det_unq_up_i)*det_ex_to_calc_sgn_up(det_unq_up_k-ndetup)
         endif
         exit
       endif
      enddo

!     if current excited determinant is a new determinant, add it to the list of unique determinants
      if (iwdet_ex_up (ex_i, det_unq_up_i) == 0) then
        det_ex_unq_up_nb = det_ex_unq_up_nb + 1
        if (det_ex_unq_up_nb>det_ex_max) call die(lhere,"Recompile with a superior value of 'det_ex_max' in deriv_orb_mod.f90 (B)")
        iwdet_ex_up (ex_i, det_unq_up_i) = ndetup+det_ex_unq_up_nb
        iwdet_ex_ref_up (det_ex_unq_up_nb) = det_unq_up_i
        det_ex_unq_orb_lab_up (:, det_ex_unq_up_nb ) = det_orb_lab_up
        det_unq_orb_lab_srt_up (:, ndetup+det_ex_unq_up_nb ) = det_orb_lab_srt_up
        det_ex_to_calc_sgn_up (det_ex_unq_up_nb) = det_ex_unq_sgn_up (ex_i, det_unq_up_i)
        det_ex_unq_sgn_up (ex_i, det_unq_up_i)=det_ex_unq_sgn_up (ex_i, det_unq_up_i)*det_ex_to_calc_sgn_up(det_ex_unq_up_nb)
        det_ex_unq_up_orb_1st_pos (det_ex_unq_up_nb) = orb_pos_in_det_unq_up (orb_opt_lab_i, det_unq_up_i)
        det_ex_unq_up_orb_2nd_lab (det_ex_unq_up_nb) = orb_opt_lab_j
      endif
    endif
    enddo ! det_unq_up_i

!1b\loop over unique spin-down determinants
    do det_unq_dn_i = 1, ndetdn
    if (orb_occ_in_det_unq_dn (orb_opt_lab_i, det_unq_dn_i) &
      .and. .not. orb_occ_in_det_unq_dn (orb_opt_lab_j, det_unq_dn_i)) then

!     build current excited determinant
      det_orb_lab_dn = det_unq_orb_lab_srt_dn (:, det_unq_dn_i)
      call replace_elt_in_array (det_orb_lab_dn, orb_opt_lab_i, orb_opt_lab_j)
      det_orb_lab_srt_dn = det_orb_lab_dn
      call sort_and_sign (det_orb_lab_srt_dn, det_ex_unq_sgn_dn (ex_i, det_unq_dn_i))

!     check if current excited determinant is an already known determinant
      do det_unq_dn_k = 1, ndetdn+det_ex_unq_dn_nb
       if (arrays_equal (det_orb_lab_srt_dn, det_unq_orb_lab_srt_dn (:, det_unq_dn_k))) then
         iwdet_ex_dn (ex_i, det_unq_dn_i) = det_unq_dn_k
         if (det_unq_dn_k.gt.ndetdn) then
           det_ex_unq_sgn_dn (ex_i, det_unq_dn_i)=det_ex_unq_sgn_dn (ex_i, det_unq_dn_i)*det_ex_to_calc_sgn_dn(det_unq_dn_k-ndetdn)
         endif
         exit
       endif
      enddo

!     if current excited determinant is a new determinant, add it to the list of unique determinants
      if (iwdet_ex_dn (ex_i, det_unq_dn_i) == 0) then
        det_ex_unq_dn_nb = det_ex_unq_dn_nb + 1
        if (det_ex_unq_dn_nb>det_ex_max) call die(lhere,"Recompile with a superior value of 'det_ex_max' in deriv_orb_mod.f90 (C)")
        iwdet_ex_dn (ex_i, det_unq_dn_i) = ndetdn+det_ex_unq_dn_nb
        iwdet_ex_ref_dn (det_ex_unq_dn_nb) = det_unq_dn_i
        det_ex_unq_orb_lab_dn (:, det_ex_unq_dn_nb) = det_orb_lab_dn
        det_unq_orb_lab_srt_dn (:, ndetdn+det_ex_unq_dn_nb) = det_orb_lab_srt_dn
        det_ex_to_calc_sgn_dn (det_ex_unq_dn_nb) = det_ex_unq_sgn_dn (ex_i, det_unq_dn_i)
        det_ex_unq_sgn_dn (ex_i, det_unq_dn_i)=det_ex_unq_sgn_dn (ex_i, det_unq_dn_i)*det_ex_to_calc_sgn_dn(det_ex_unq_dn_nb)
        det_ex_unq_dn_orb_1st_pos (det_ex_unq_dn_nb) = orb_pos_in_det_unq_dn (orb_opt_lab_i, det_unq_dn_i)
        det_ex_unq_dn_orb_2nd_lab (det_ex_unq_dn_nb) = orb_opt_lab_j
      endif
      if (l_triplet) det_ex_unq_sgn_dn (ex_i, det_unq_dn_i ) = -det_ex_unq_sgn_dn (ex_i, det_unq_dn_i )

    endif
    enddo ! det_unq_dn_i
  call cpu_time(tB)
  tDet=tDet+tB-tA
  tA=tB


!2\ ================================================================================= 
!   Gather and sort information on the excited CSF
!   (in the process: gather information on unique CSF)
!   ================================================================================= 
    !CSF_in_WF initialization
    csf_in_wf_nb  =0
    csf_in_wf     =csf_in_wf_max+1
    csf_in_wf_ref =csf_in_wf_max+1
    csf_in_wf_coef=csf_in_wf_max+1

    do csf_i = 1, ncsf

!2a\  Gather the unsorted information on CSFs
!     (loop through determinants in the reference CSF)
      !DET_in_CSF_tmp initialization
      det_in_csf_nb_tmp  =0
      det_in_csf_up_tmp  =det_in_csf_max+1
      det_in_csf_dn_tmp  =det_in_csf_max+1
      det_in_csf_coef_tmp=det_in_csf_max+1
      write(6,'(''csf_i. ndet_in_csf (csf_i), det_in_csf_nb_tmp'',9i5)') csf_i, ndet_in_csf (csf_i), det_in_csf_nb_tmp
      do det_in_csf_i = 1, ndet_in_csf (csf_i)
        det_i=iwdet_in_csf(det_in_csf_i, csf_i)
        det_unq_up_i = det_to_det_unq_up (det_i)
        det_unq_dn_i = det_to_det_unq_dn (det_i)

        iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
        if (iwdet.ne.0) then
          det_in_csf_nb_tmp=det_in_csf_nb_tmp+1
          if (det_in_csf_nb_tmp>det_in_csf_max) then
             write(6,'(''1csf_i, ndet_in_csf (csf_i), det_in_csf_i, det_in_csf_nb_tmp, det_in_csf_max='',9i5)') csf_i, ndet_in_csf (csf_i), det_in_csf_i, det_in_csf_nb_tmp, det_in_csf_max 
             call systemflush(6)
             call die(lhere,"Recompile with a superior value of 'det_in_csf_max' in deriv_orb_mod.f90 (D)")
          endif
          det_in_csf_up_tmp  (det_in_csf_nb_tmp)=iwdet
          det_in_csf_dn_tmp  (det_in_csf_nb_tmp)=det_unq_dn_i
          det_in_csf_coef_tmp(det_in_csf_nb_tmp)=cdet_in_csf(det_in_csf_i, csf_i)*det_ex_unq_sgn_up(ex_i, det_unq_up_i)
        endif
        iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
        if (iwdet.ne.0) then
          det_in_csf_nb_tmp=det_in_csf_nb_tmp+1
          if (det_in_csf_nb_tmp>det_in_csf_max) then
             write(6,'(''2csf_i, det_in_csf_i, det_in_csf_nb_tmp, det_in_csf_max='',9i5)') csf_i, det_in_csf_i, det_in_csf_nb_tmp, det_in_csf_max
             call systemflush(6)
             call die(lhere,"Recompile with a superior value of 'det_in_csf_max' in deriv_orb_mod.f90 (E)")
          endif
          det_in_csf_up_tmp  (det_in_csf_nb_tmp)=det_unq_up_i
          det_in_csf_dn_tmp  (det_in_csf_nb_tmp)=iwdet
          det_in_csf_coef_tmp(det_in_csf_nb_tmp)=cdet_in_csf(det_in_csf_i, csf_i)*det_ex_unq_sgn_dn(ex_i, det_unq_dn_i)
        endif
      enddo ! det_in_csf_i


!2b\  Sort the CSF information
!     and fold repetitions (non unique couple of UP and DOWN determinants) together
      if (det_in_csf_nb_tmp.ne.0) then
!       Here we do resize the arrays to make the sort quicker
        call alloc('det_in_csf_up_tmp',det_in_csf_up_tmp,det_in_csf_nb_tmp)
        call alloc('det_in_csf_dn_tmp',det_in_csf_dn_tmp,det_in_csf_nb_tmp)
        call alloc('det_in_csf_coef_tmp',det_in_csf_coef_tmp,det_in_csf_nb_tmp)
  call cpu_time(tAintermed)
        call sort(det_in_csf_up_tmp, det_in_csf_dn_tmp, det_in_csf_coef_tmp)
  call cpu_time(tBintermed)
  tSort=tSort+tBintermed-tAintermed
        call alloc('det_in_csf_up_tmp',det_in_csf_up_tmp,det_in_csf_max)
        call alloc('det_in_csf_dn_tmp',det_in_csf_dn_tmp,det_in_csf_max)
        call alloc('det_in_csf_coef_tmp',det_in_csf_coef_tmp,det_in_csf_max)

        det_in_csf_nb  =1
        det_in_csf_up  =0
        det_in_csf_dn  =0
        det_in_csf_coef=0 
        det_in_csf_up(1)  =det_in_csf_up_tmp(1)
        det_in_csf_dn(1)  =det_in_csf_dn_tmp(1)
        det_in_csf_coef(1)=det_in_csf_coef_tmp(1)
        do csf_k=2, det_in_csf_nb_tmp
          if ((det_in_csf_up_tmp(csf_k)==det_in_csf_up(det_in_csf_nb)).and.&
              (det_in_csf_dn_tmp(csf_k)==det_in_csf_dn(det_in_csf_nb))) then
            det_in_csf_coef(det_in_csf_nb)=det_in_csf_coef(det_in_csf_nb)+det_in_csf_coef_tmp(csf_k)
          elseif (det_in_csf_up_tmp(csf_k).ne.det_in_csf_max+1) then
            det_in_csf_nb=det_in_csf_nb+1
            if (det_in_csf_nb>det_in_csf_max) then
               write(6,'(''3csf_i, det_in_csf_i, det_in_csf_nb_tmp, det_in_csf_max='',9i5)') csf_i, det_in_csf_i, det_in_csf_nb_tmp, det_in_csf_max
               call systemflush(6) 
               call die(lhere,"Recompile with a superior value of 'det_in_csf_max' in deriv_orb_mod.f90 (F)")
            endif
            det_in_csf_up  (det_in_csf_nb)=det_in_csf_up_tmp(csf_k)
            det_in_csf_dn  (det_in_csf_nb)=det_in_csf_dn_tmp(csf_k)
            det_in_csf_coef(det_in_csf_nb)=det_in_csf_coef_tmp(csf_k)
          endif
        enddo
  call cpu_time(tB)
  tCSFa=tCSFa+tB-tA
  tA=tB

!2c\    Check if this properly sorted data on CSF is an already known CSF...
        this_csf=0
        do csf_k=1, ncsf+csf_ex_unq_nb  
          if (arrays_equal(det_in_csf_up(:det_in_csf_nb) , det_unq_up_in_csf (csf_k)%row (:)).and.&
              arrays_equal(det_in_csf_dn(:det_in_csf_nb) , det_unq_dn_in_csf (csf_k)%row (:)).and.&
              arrays_equal(det_in_csf_coef(:det_in_csf_nb),cdet_unq_in_csf   (csf_k)%row (:))) then
            this_csf=csf_k
          endif
        enddo
  call cpu_time(tB)
  tCSFb=tCSFb+tB-tA
  tA=tB
!       ... and if this CSF is unknown, add it
        if (this_csf==0) then
          csf_ex_unq_nb=csf_ex_unq_nb+1
          if (csf_ex_unq_nb>csf_ex_unq_max) call die(lhere,"Recompile with a superior value of 'csf_ex_unq_max' in deriv_orb_mod.f90 (G)")
          this_csf=ncsf+csf_ex_unq_nb
          csf_ex_unq_ref(csf_ex_unq_nb) =csf_i
          csf_ex_unq_size(csf_ex_unq_nb)=det_in_csf_nb
          det_unq_up_in_csf(ncsf+csf_ex_unq_nb)%row=det_in_csf_up(:det_in_csf_nb)
          det_unq_dn_in_csf(ncsf+csf_ex_unq_nb)%row=det_in_csf_dn(:det_in_csf_nb)
          cdet_unq_in_csf  (ncsf+csf_ex_unq_nb)%row=det_in_csf_coef(:det_in_csf_nb)
        endif

!3a\    Gather this CSF into the information on the excited WF
        csf_in_wf_nb=csf_in_wf_nb+1
        if (csf_in_wf_nb>csf_in_wf_max) call die(lhere,"Recompile with a superior value of 'csf_in_wf_max' in deriv_orb_mod.f90 (H)")
        csf_in_wf(csf_in_wf_nb)     =this_csf
        csf_in_wf_ref(csf_in_wf_nb) =csf_i
        csf_in_wf_coef(csf_in_wf_nb)=det_in_csf_coef(1)
  call cpu_time(tB)
  tCSFc=tCSFc+tB-tA
  tA=tB

      endif ! csf_ex != 0
    enddo ! csf_i

!3\ ================================================================================= 
!   Sort information on the excited WF
!   (in the process: gather information on unique WF)
!   ================================================================================= 
    if (csf_in_wf_nb.ne.0) then
!3b\  Sort the WF information
!     Here we do resize the arrays to make the sort quicker
      call alloc('csf_in_wf',csf_in_wf,csf_in_wf_nb)
      call alloc('csf_in_wf_ref',csf_in_wf_ref,csf_in_wf_nb)
      call alloc('csf_in_wf_coef',csf_in_wf_coef,csf_in_wf_nb)
  call cpu_time(tAintermed)
      call sort(csf_in_wf, csf_in_wf_ref, csf_in_wf_coef)
  call cpu_time(tBintermed)
  tSort=tSort+tBintermed-tAintermed
      call alloc('csf_in_wf',csf_in_wf,csf_in_wf_max)
      call alloc('csf_in_wf_ref',csf_in_wf_ref,csf_in_wf_max)
      call alloc('csf_in_wf_coef',csf_in_wf_coef,csf_in_wf_max)

!3c\  Check if this properly sorted data on WF is an already known WF...
      this_ex=0
      do ex_j = 1, single_ex_nb
        if (arrays_equal (csf_in_wf      (:csf_in_wf_nb) , wf_unq_csf (ex_j)%row (:)) .and. &
            arrays_equal (csf_in_wf_ref  (:csf_in_wf_nb) , wf_unq_ref (ex_j)%row (:)) .and. &
            arrays_equal (csf_in_wf_coef (:csf_in_wf_nb) , wf_unq_coef(ex_j)%row (:))) then
          this_ex = ex_j
        endif
      enddo
!     ...if this WF is unknown, add it
      if (this_ex == 0) then
        single_ex_nb = single_ex_nb + 1
        if (single_ex_nb>single_ex_max) call die(lhere,"Recompile with a superior value of 'single_ex_max' in deriv_orb_mod.f90 (I)")
        single_ex_added_nb = single_ex_added_nb + 1
        this_ex = single_ex_nb
        wf_unq_nb  (single_ex_nb)   = csf_in_wf_nb
        call copy_portion(csf_in_wf     , wf_unq_csf (single_ex_nb)%row, csf_in_wf_nb)
        call copy_portion(csf_in_wf_ref , wf_unq_ref (single_ex_nb)%row, csf_in_wf_nb)
        call copy_portion(csf_in_wf_coef, wf_unq_coef(single_ex_nb)%row, csf_in_wf_nb)
        ex_orb_ind     (single_ex_nb) = single_ex_nb
        ex_orb_1st_lab (single_ex_nb) = orb_opt_lab_i
        ex_orb_2nd_lab (single_ex_nb) = orb_opt_lab_j
      endif
  call cpu_time(tB)
  tWF=tWF+tB-tA
  tA=tB

!4\ ================================================================================= 
!   Gather and sort information on the derivative DPSI
!   (in the process: gather information on unique versus redundant DPSI)
!   ================================================================================= 
!4a\  Gather information
      if (ex_dir_rev==1) then
        ex_orb_ind_tmp = this_ex
        ex_orb_ind_rev_tmp = 0
      else
        ex_orb_ind_rev_tmp =  this_ex
        ex_orb_ind_tmp = 0
        csf_in_wf_coef(:csf_in_wf_nb)=-csf_in_wf_coef(:csf_in_wf_nb)
      endif
      if (csf_in_dpsi_nb_tmp+csf_in_wf_nb>csf_in_dpsi_max) call die(lhere,"Recompile with a superior value of 'csf_in_dpsi_max' in deriv_orb_mod.f90 (J)")
      csf_in_dpsi_tmp     (csf_in_dpsi_nb_tmp+1:csf_in_dpsi_nb_tmp+csf_in_wf_nb)=csf_in_wf     (:csf_in_wf_nb)
      csf_in_dpsi_ref_tmp (csf_in_dpsi_nb_tmp+1:csf_in_dpsi_nb_tmp+csf_in_wf_nb)=csf_in_wf_ref (:csf_in_wf_nb)
      csf_in_dpsi_coef_tmp(csf_in_dpsi_nb_tmp+1:csf_in_dpsi_nb_tmp+csf_in_wf_nb)=csf_in_wf_coef(:csf_in_wf_nb)
      csf_in_dpsi_nb_tmp=csf_in_dpsi_nb_tmp+csf_in_wf_nb

    endif ! wf_ex != 0
  enddo ! ex_dir_rev

!4b\Sort the information
!   and fold repetitions (non unique couple of CSF and REF) together
  if (csf_in_dpsi_nb_tmp.ne.0) then
    dpsi_is_redundant=.false.

!   Here we do resize the arrays to make the sort quicker
    call alloc('csf_in_dpsi_tmp',csf_in_dpsi_tmp,csf_in_dpsi_nb_tmp)
    call alloc('csf_in_dpsi_ref_tmp',csf_in_dpsi_ref_tmp,csf_in_dpsi_nb_tmp)
    call alloc('csf_in_dpsi_coef_tmp',csf_in_dpsi_coef_tmp,csf_in_dpsi_nb_tmp)
 call cpu_time(tAintermed)
    call sort(csf_in_dpsi_tmp,csf_in_dpsi_ref_tmp,csf_in_dpsi_coef_tmp)
 call cpu_time(tBintermed)
 tSort=tSort+tBintermed-tAintermed
    call alloc('csf_in_dpsi_tmp',csf_in_dpsi_tmp,csf_in_dpsi_max)
    call alloc('csf_in_dpsi_ref_tmp',csf_in_dpsi_ref_tmp,csf_in_dpsi_max)
    call alloc('csf_in_dpsi_coef_tmp',csf_in_dpsi_coef_tmp,csf_in_dpsi_max)

    csf_in_dpsi_nb=1
    csf_in_dpsi(1)     =csf_in_dpsi_tmp(1)
    csf_in_dpsi_ref(1) =csf_in_dpsi_ref_tmp(1)
    csf_in_dpsi_coef(1)=csf_in_dpsi_coef_tmp(1)
    do dpsi_k=2,csf_in_dpsi_nb_tmp
      if ((csf_in_dpsi_tmp(dpsi_k)    ==csf_in_dpsi(csf_in_dpsi_nb)) .and.  &
          (csf_in_dpsi_ref_tmp(dpsi_k)==csf_in_dpsi_ref(csf_in_dpsi_nb))) then
        csf_in_dpsi_coef(csf_in_dpsi_nb)=csf_in_dpsi_coef(csf_in_dpsi_nb)+csf_in_dpsi_coef_tmp(dpsi_k)
      elseif (csf_in_dpsi_tmp(dpsi_k).ne.csf_in_dpsi_max+1) then
        csf_in_dpsi_nb=csf_in_dpsi_nb + 1
        if (csf_in_dpsi_nb>csf_in_dpsi_max) call die(lhere,"Recompile with a superior value of 'csf_in_dpsi_max' in deriv_orb_mod.f90 (K)")
        csf_in_dpsi(csf_in_dpsi_nb)     =csf_in_dpsi_tmp(dpsi_k)
        csf_in_dpsi_ref(csf_in_dpsi_nb) =csf_in_dpsi_ref_tmp(dpsi_k)
        csf_in_dpsi_coef(csf_in_dpsi_nb)=csf_in_dpsi_coef_tmp(dpsi_k)
      endif
    enddo

!4c\Check if this DPSI is redundant...
    !if all first coefficients are zero
    dpsi_is_redundant=.true.
    do dpsi_k=1,csf_in_dpsi_nb
      if (csf_in_dpsi_coef(dpsi_k).ne.0.d0) then
        dpsi_is_redundant=.false.
        exit
      endif
    enddo

    !if all CSF are from the ground-state
    if ((.not. dpsi_is_redundant).and.(l_opt_csf)) then
      dpsi_is_redundant=.true.
      do dpsi_k=1,csf_in_dpsi_nb
        if (((.not.(elt_in_array(iwcsf, csf_in_dpsi(dpsi_k)))) &
             .and. (csf_in_dpsi_coef(dpsi_k).ne.0.d0))             &
           .or. (csf_in_dpsi_nb==1)) then
          dpsi_is_redundant=.false.
        endif
      enddo
    endif

    !if is it a DPSI already known
    if (.not. dpsi_is_redundant) then
     do dpsi_i = 1, param_orb_nb
      if ((arrays_equal (csf_in_dpsi (:csf_in_dpsi_nb),     dpsi_unq_csf (dpsi_i)%row (:))) &
       .and.( &
          (arrays_equal (csf_in_dpsi_ref (:csf_in_dpsi_nb), dpsi_unq_ref (dpsi_i)%row (:)) .and. & 
           arrays_equal (csf_in_dpsi_coef(:csf_in_dpsi_nb), dpsi_unq_coef(dpsi_i)%row (:))) &
          .or.(csf_in_dpsi_nb == 1))  ) then
         dpsi_is_redundant = .true.
      endif
     enddo
    endif

    !if redundant: remove added unique singly-excited wave functions
    if (dpsi_is_redundant) then
      single_ex_nb = single_ex_nb - single_ex_added_nb
      !NEW
      !wf_unq_nb     (single_ex_nb+1:)=0
      !do ireinit=single_ex_nb+1,wf_unq_max
      !  wf_unq_csf (ireinit)%row=0
      !  wf_unq_ref (ireinit)%row=0
      !  wf_unq_coef(ireinit)%row=0
      !enddo
      !ex_orb_1st_lab(single_ex_nb+1:)=0
      !ex_orb_2nd_lab(single_ex_nb+1:)=0
      !OLD
      !call object_alloc ('wf_unq_nb',   wf_unq_nb,   single_ex_nb)
      !call object_alloc ('wf_unq_csf',  wf_unq_csf,  single_ex_nb)
      !call object_alloc ('wf_unq_ref',  wf_unq_ref,  single_ex_nb)
      !call object_alloc ('wf_unq_coef', wf_unq_coef, single_ex_nb)
      !call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
      !call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)
    endif

!4d\...or unique
    if (.not. dpsi_is_redundant) then
      param_orb_nb = param_orb_nb + 1
      if (param_orb_nb>single_ex_max) call die(lhere,"Recompile with a superior value of 'mamax' in deriv_orb_mod.f90 (L)")
      dpsi_unq_nb  (param_orb_nb)    =csf_in_dpsi_nb
      call copy_portion(csf_in_dpsi     , dpsi_unq_csf (param_orb_nb)%row, csf_in_dpsi_nb)
      call copy_portion(csf_in_dpsi_ref , dpsi_unq_ref (param_orb_nb)%row, csf_in_dpsi_nb)
      call copy_portion(csf_in_dpsi_coef, dpsi_unq_coef(param_orb_nb)%row, csf_in_dpsi_nb)

      if (ex_orb_ind_tmp /= 0) then
       ex_orb_ind    (param_orb_nb) = ex_orb_ind_tmp
       ex_orb_ind_rev(param_orb_nb) = ex_orb_ind_rev_tmp
      else
!      Warning: I did not check this very carefully
       ex_orb_ind    (param_orb_nb) = ex_orb_ind_rev_tmp
       ex_orb_ind_rev(param_orb_nb) = ex_orb_ind_tmp
      endif
    endif
  call cpu_time(tB)
  tDPSI=tDPSI+tB-tA
  tA=tB

  endif ! dpsi_ex != 0
  enddo ! ex_dir_rev_nonortho

  orb_opt_lab_i = orb_opt_lab_j
  enddo !orb_opt_lab_i
  enddo !orb_opt_lab_j
  call cpu_time(tEnd)
  write(6,'(a)') '   Timing of the different steps:'
  write(6,'(a,f8.3)')  '     Det:   ',tDet
  write(6,'(a,f8.3)')  '     CSFa:  ',tCSFa
  write(6,'(a,f8.3,a)')'     CSFb:  ',tCSFb,'  <- usually the most demanding'
  write(6,'(a,f8.3)')  '     CSFc:  ',tCSFc
  write(6,'(a,f8.3)')  '     WF:    ',tWF
  write(6,'(a,f8.3)')  '     DPSI:  ',tDPSI
  write(6,'(a,f8.3)')  '     Total: ',tEnd-tBegin
  write(6,'(a,f8.3)')  '     Sort:  ',tSort
  write(6,'(a)') '   Max of the array dimensions versus actual dimensions (and ratio):'
  write(6,'(a,i8,i8,i8)') '     single_ex:  ',single_ex_max  ,single_ex_nb    ,single_ex_max  /single_ex_nb
  write(6,'(a,i8,i8,i8)') '     det_ex:     ',det_ex_max     ,det_ex_unq_up_nb,det_ex_max     /det_ex_unq_up_nb
  write(6,'(a,i8,i8,i8)') '     csf_ex_unq: ',csf_ex_unq_max ,csf_ex_unq_nb   ,csf_ex_unq_max /csf_ex_unq_nb
  write(6,'(a,i8,i8,i8)') '     det_in_csf: ',det_in_csf_max ,det_in_csf_nb   ,det_in_csf_max /det_in_csf_nb
  write(6,'(a,i8,i8,i8)') '     csf_in_wf:  ',csf_in_wf_max  ,csf_in_wf_nb    ,csf_in_wf_max  /csf_in_wf_nb
  write(6,'(a,i8,i8,i8)') '     csf_in_dpsi:',csf_in_dpsi_max,csf_in_dpsi_nb  ,csf_in_dpsi_max/csf_in_dpsi_nb
  write(6,'(a,i8,i8,i8)') '     wf_unq:     ',wf_unq_max     ,wf_unq_nb(1)    ,wf_unq_max     /wf_unq_nb(1)
  write(6,'(a,i8,i8,i8)') '     dpsi_unq:   ',dpsi_unq_max   ,dpsi_unq_nb(1)  ,dpsi_unq_max   /dpsi_unq_nb(1)

! determine orbitals mixed in orbital optimization 
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

! excitation pairs number
  deriv_orb_pairs_nb = param_orb_nb * (param_orb_nb + 1) / 2

  write(6,'(a,i10)') ' Number of single orbital excitations = ', single_ex_nb
  write(6,'(a,i10)') ' Number of orbital derivatives        = ', param_orb_nb
  write(6,'(a,i10)') ' Number of orbitals involved          = ', orb_mix_lab_nb
  write(6,'(a,i10)') ' Number of unique spin-up   excited determinants = ', det_ex_unq_up_nb
  write(6,'(a,i10)') ' Number of unique spin-down excited determinants = ', det_ex_unq_dn_nb

! adjust the arrays to their real size
  call object_alloc ('iwdet_ex_up',               iwdet_ex_up,                 single_ex_nb, ndetup)
  call object_alloc ('iwdet_ex_dn',               iwdet_ex_dn,                 single_ex_nb, ndetdn)
  call object_alloc ('det_ex_unq_sgn_up',         det_ex_unq_sgn_up,           single_ex_nb, ndetup)
  call object_alloc ('det_ex_unq_sgn_dn',         det_ex_unq_sgn_dn,           single_ex_nb, ndetdn)
  call object_alloc ('iwdet_ex_ref_up',           iwdet_ex_ref_up,             det_ex_unq_up_nb)
  call object_alloc ('iwdet_ex_ref_dn',           iwdet_ex_ref_dn,             det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_orb_lab_up',     det_ex_unq_orb_lab_up, nup,  det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_orb_lab_dn',     det_ex_unq_orb_lab_dn, ndn,  det_ex_unq_up_nb)
  call object_alloc ('det_unq_orb_lab_srt_up',    det_unq_orb_lab_srt_up, nup, ndetup+det_ex_unq_up_nb)
  call object_alloc ('det_unq_orb_lab_srt_dn',    det_unq_orb_lab_srt_dn, ndn, ndetdn+det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_up_orb_1st_pos', det_ex_unq_up_orb_1st_pos,   det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_up_orb_2nd_lab', det_ex_unq_up_orb_2nd_lab,   det_ex_unq_up_nb)
  call object_alloc ('det_ex_unq_dn_orb_1st_pos', det_ex_unq_dn_orb_1st_pos,   det_ex_unq_dn_nb)
  call object_alloc ('det_ex_unq_dn_orb_2nd_lab', det_ex_unq_dn_orb_2nd_lab,   det_ex_unq_dn_nb)

  call object_alloc('det_unq_up_in_csf', det_unq_up_in_csf, ncsf+csf_ex_unq_nb)
  call object_alloc('det_unq_dn_in_csf', det_unq_dn_in_csf, ncsf+csf_ex_unq_nb)
  call object_alloc('cdet_unq_in_csf'  , cdet_unq_in_csf,   ncsf+csf_ex_unq_nb)

  call object_alloc ('ex_orb_ind',     ex_orb_ind,     param_orb_nb)
  call object_alloc ('ex_orb_ind_rev', ex_orb_ind_rev, param_orb_nb)
  call object_alloc ('ex_orb_1st_lab', ex_orb_1st_lab, single_ex_nb)
  call object_alloc ('ex_orb_2nd_lab', ex_orb_2nd_lab, single_ex_nb)

  call release ('det_orb_lab_up',     det_orb_lab_up)
  call release ('det_orb_lab_dn',     det_orb_lab_dn)
  call release ('det_orb_lab_srt_up', det_orb_lab_srt_up)
  call release ('det_orb_lab_srt_dn', det_orb_lab_srt_dn)
  call release ('det_ex_to_calc_sgn_up',  det_ex_to_calc_sgn_up)
  call release ('det_ex_to_calc_sgn_dn',  det_ex_to_calc_sgn_dn)

  call release('det_in_csf_up_tmp',   det_in_csf_up_tmp)
  call release('det_in_csf_dn_tmp',   det_in_csf_dn_tmp)
  call release('det_in_csf_coef_tmp', det_in_csf_coef_tmp)
  call release('det_in_csf_up',       det_in_csf_up)
  call release('det_in_csf_dn',       det_in_csf_dn)
  call release('det_in_csf_coef',     det_in_csf_coef)

  call release('csf_ex_unq_ref',           csf_ex_unq_ref)
  call release('csf_ex_unq_size',          csf_ex_unq_size)
  
  call release('csf_in_wf',      csf_in_wf)
  call release('csf_in_wf_ref',  csf_in_wf_ref)
  call release('csf_in_wf_coef', csf_in_wf_coef)

  call release('wf_unq_nb',      wf_unq_nb)
  call release('wf_unq_csf',     wf_unq_csf)
  call release('wf_unq_ref',     wf_unq_ref)
  call release('wf_unq_coef',    wf_unq_coef)

  call release('csf_in_dpsi_tmp',      csf_in_dpsi_tmp)
  call release('csf_in_dpsi_ref_tmp',  csf_in_dpsi_ref_tmp)
  call release('csf_in_dpsi_coef_tmp', csf_in_dpsi_coef_tmp)
  call release('csf_in_dpsi',          csf_in_dpsi)
  call release('csf_in_dpsi_ref',      csf_in_dpsi_ref)
  call release('csf_in_dpsi_coef',     csf_in_dpsi_coef)

  call release('dpsi_unq_nb',          dpsi_unq_nb)
  call release('dpsi_unq_csf',         dpsi_unq_csf)
  call release('dpsi_unq_ref',         dpsi_unq_ref)
  call release('dpsi_unq_coef',        dpsi_unq_coef)

  !stop

  end subroutine single_ex_wf_bld

! ==============================================================================
  subroutine det_ex_unq_bld
! ------------------------------------------------------------------------------
! Description   :  Calculate value of unque singly-excited spin-up and spin-down determinants
! Description   :  by updating reference determinant via the Sherman-Morrison formula
!
! Created       : J. Toulouse, 27 Oct 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer det_ex_unq_up_i, det_ex_unq_dn_i
  integer iwdet_ref
  integer i
  integer col_up, orb_up
  integer col_dn, orb_dn
  real(dp) factor_up, factor_dn

! header
  if (header_exe) then

   call object_create ('det_ex_unq_up')
   call object_create ('det_ex_unq_dn')

   call object_needed ('det_ex_unq_up_nb')
   call object_needed ('det_ex_unq_dn_nb')
   call object_needed ('det_ex_unq_up_orb_1st_pos')
   call object_needed ('det_ex_unq_dn_orb_1st_pos')
   call object_needed ('det_ex_unq_up_orb_2nd_lab')
   call object_needed ('det_ex_unq_dn_orb_2nd_lab')
   call object_needed ('iwdet_ex_ref_up')
   call object_needed ('iwdet_ex_ref_dn')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('orb')
   call object_needed ('detu')
   call object_needed ('detd')

   return

  endif

! begin

! allocations
  call alloc ('det_ex_unq_up', det_ex_unq_up, det_ex_unq_up_nb)
  call alloc ('det_ex_unq_dn', det_ex_unq_dn, det_ex_unq_dn_nb)

! spin up determinants
  do det_ex_unq_up_i = 1, det_ex_unq_up_nb

    col_up  = det_ex_unq_up_orb_1st_pos (det_ex_unq_up_i)
    orb_up  = det_ex_unq_up_orb_2nd_lab (det_ex_unq_up_i)
    iwdet_ref = iwdet_ex_ref_up (det_ex_unq_up_i)

    factor_up = 0.d0
    do i = 1, nup
     factor_up = factor_up + slater_mat_trans_inv_up (i, col_up, iwdet_ref) * orb (i, orb_up)
    enddo

    det_ex_unq_up(det_ex_unq_up_i) = factor_up * detu (iwdet_ref)

    !write(6,*) '>det_ex_unq_up',det_ex_unq_up(det_ex_unq_up_i)
  enddo ! det_ex_unq_up_i

! spin down determinants
  do det_ex_unq_dn_i = 1, det_ex_unq_dn_nb

    col_dn  = det_ex_unq_dn_orb_1st_pos (det_ex_unq_dn_i)
    orb_dn  = det_ex_unq_dn_orb_2nd_lab (det_ex_unq_dn_i)
    iwdet_ref = iwdet_ex_ref_dn (det_ex_unq_dn_i)

    factor_dn = 0.d0
    do i = 1, ndn
     factor_dn = factor_dn + slater_mat_trans_inv_dn (i, col_dn, iwdet_ref) * orb (nup + i, orb_dn)
    enddo

    det_ex_unq_dn(det_ex_unq_dn_i) = factor_dn * detd (iwdet_ref)

    !write(6,*) '>det_ex_unq_dn',det_ex_unq_dn(det_ex_unq_dn_i)
  enddo ! det_ex_unq_dn_i

  end subroutine det_ex_unq_bld

! ==============================================================================
  subroutine det_ex_bld
! ------------------------------------------------------------------------------
! Description   : all excited determinants
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer ex_i
  integer csf_i, det_in_csf_i, det_i, det_unq_up_i, det_unq_dn_i
  integer iwdet, sgn

! header
  if (header_exe) then

   call object_create ('det_ex_up')
   call object_create ('det_ex_dn')
   call object_create ('det_ex')

   call object_needed ('single_ex_nb')
   call object_needed ('ndet')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('iwdet_ex_up')
   call object_needed ('iwdet_ex_dn')
   call object_needed ('det_ex_unq_sgn_up')
   call object_needed ('det_ex_unq_sgn_dn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('det_ex_unq_up')
   call object_needed ('det_ex_unq_dn')

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
        iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
        sgn = det_ex_unq_sgn_up(ex_i,det_unq_up_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetup)) then
          det_ex_up (ex_i, det_i) = sgn * detu (iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetup+det_ex_unq_up_nb)) then
          det_ex_up (ex_i, det_i) = sgn * det_ex_unq_up (iwdet-ndetup)
        endif


!       spin-down excited determinants:
        iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
        sgn = det_ex_unq_sgn_dn(ex_i,det_unq_dn_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetdn)) then
          det_ex_dn (ex_i, det_i) = sgn * detd (iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetdn+det_ex_unq_dn_nb)) then
          det_ex_dn (ex_i, det_i) = sgn * det_ex_unq_dn (iwdet-ndetdn)
        endif

        det_ex (ex_i, det_i) = det_ex_up (ex_i, det_i) * detd (det_unq_dn_i) + detu (det_unq_up_i) * det_ex_dn (ex_i, det_i)

     enddo ! det_in_csf_i
   enddo ! csf_i

  enddo ! ex_i

  !write(6,*) '>det_ex',det_ex

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
  use all_modules_mod
  implicit none

! local
  integer ex_i, ex_rev_i
  integer csf_i, det_in_csf_i, det_i, dorb_i
!  real(dp) factor_up, factor_dn
  real(dp) detex
  real(dp), allocatable                  :: dpsi_orb_test(:)

! header
  if (header_exe) then

   call object_create ('dpsi_orb')
   call object_create ('dcsf_orb')
   call object_create ('psid_ex')

   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_ind_rev')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_ex')
   call object_needed ('cdet_in_csf')
   call object_needed ('csf_coef')
   call object_needed ('psi_det')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsi_orb', dpsi_orb, param_orb_nb)
  call object_alloc ('dcsf_orb', dcsf_orb, ncsf, param_orb_nb)
  call object_alloc ('psid_ex', psid_ex, param_orb_nb)

  psid_ex = 0.d0
  dcsf_orb =  0.d0

  do dorb_i = 1, param_orb_nb
    ex_i = ex_orb_ind (dorb_i)
    ex_rev_i = ex_orb_ind_rev (dorb_i)

    do csf_i = 1, ncsf
      do det_in_csf_i = 1, ndet_in_csf (csf_i)
        det_i = iwdet_in_csf (det_in_csf_i, csf_i)

        detex = det_ex (ex_i, det_i)
!       reverse excitation for non casscf wave function: (Eij-Eji) Det
        if (.not. l_casscf .and. ex_rev_i /= 0) then
          detex = detex - det_ex (ex_rev_i, det_i)
        endif

        dcsf_orb(csf_i,dorb_i)=dcsf_orb(csf_i,dorb_i) + cdet_in_csf (det_in_csf_i, csf_i) * detex
      enddo ! det_in_csf_i

      psid_ex (dorb_i) = psid_ex (dorb_i) + csf_coef (csf_i, 1) * dcsf_orb(csf_i,dorb_i)
      dcsf_orb(csf_i,dorb_i)=dcsf_orb(csf_i,dorb_i) /  psi_det

    enddo ! csf_i
    dpsi_orb (dorb_i) = psid_ex (dorb_i) / psi_det
    !write(6,*) '>dpsi_orb',dpsi_orb(dorb_i)
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
! Description   :  inverse of transpose of Slater matrix corresponding to
! Description   :  to excited determinants using the Sherman-Morison formula with O(nup^2 + ndn^2) scaling
! Description   :  or with inversion from scratch (seems faster for small electron numbers)
! Description   :  memory could be saved by using either removing slater_mat_ex_trans_up and slater_mat_ex_trans_dn
!
! Created       : J. Toulouse, 27 Oct 2005
! Modified      : J. Toulouse, 22 Apr 2015: merge Sherman-Morison and inversion from scratch in a single subroutine
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer det_up_i, det_dn_i 
  integer iwdet_ref
  integer i, j, k, l
  integer col_up, orb_up
  integer col_dn, orb_dn
  real(dp) factor_up_inv, factor_dn_inv
  real(dp), allocatable :: ratio_up (:), ratio_dn (:)

  integer det_i, orb_i, elec_i
  real(dp), allocatable                  :: slater_mat_ex_trans_up (:,:)
  real(dp), allocatable                  :: slater_mat_ex_trans_dn (:,:)
  real(dp), allocatable                  :: mat_flat_up (:)
  real(dp), allocatable                  :: mat_flat_dn (:)

! header
  if (header_exe) then

   call object_create ('slater_mat_ex_trans_inv_up')
   call object_create ('slater_mat_ex_trans_inv_dn')

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('det_ex_unq_up_nb')
   call object_needed ('det_ex_unq_dn_nb')
   call object_needed ('det_ex_unq_up_orb_1st_pos')
   call object_needed ('det_ex_unq_dn_orb_1st_pos')
   call object_needed ('det_ex_unq_up_orb_2nd_lab')
   call object_needed ('det_ex_unq_dn_orb_2nd_lab')
   call object_needed ('iwdet_ex_ref_up')
   call object_needed ('iwdet_ex_ref_dn')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('orb')

   call object_needed ('det_ex_unq_orb_lab_up')
   call object_needed ('det_ex_unq_orb_lab_dn')
   call object_needed ('det_ex_unq_up')
   call object_needed ('det_ex_unq_dn')
   call object_needed ('detu')
   call object_needed ('detd')

   return

  endif

! allocations
  call object_alloc ('slater_mat_ex_trans_inv_up', slater_mat_ex_trans_inv_up, nup, nup, det_ex_unq_up_nb)
  call object_alloc ('slater_mat_ex_trans_inv_dn', slater_mat_ex_trans_inv_dn, ndn, ndn, det_ex_unq_dn_nb)

! calculate inverses with Sherman-Morison
  if (l_slater_mat_ex_trans_inv_sm) then

!  spin up
   call alloc ('ratio_up', ratio_up, nup)
   do det_up_i = 1, det_ex_unq_up_nb
   
       col_up  = det_ex_unq_up_orb_1st_pos (det_up_i)
       orb_up  = det_ex_unq_up_orb_2nd_lab (det_up_i)
       iwdet_ref = iwdet_ex_ref_up (det_up_i)
   
       do l = 1, nup
        ratio_up (l) = 0.d0
        do i = 1, nup
          ratio_up (l) = ratio_up (l) + slater_mat_trans_inv_up (i, l, iwdet_ref) * orb (i, orb_up)
        enddo  ! i
       enddo ! l
   
       factor_up_inv = 1.d0/ratio_up (col_up)
       ratio_up (col_up) = ratio_up (col_up) - 1.d0
   
       do k = 1, nup
         do j = 1 , nup
          slater_mat_ex_trans_inv_up (k, j, det_up_i) = slater_mat_trans_inv_up (k, j, iwdet_ref) &
                                  - slater_mat_trans_inv_up (k, col_up, iwdet_ref) * ratio_up(j)* factor_up_inv
         enddo ! j
       enddo ! k
   
    enddo ! det_up_i
    call release ('ratio_up', ratio_up)
   
!  spin dn
   call alloc ('ratio_dn', ratio_dn, ndn)
   do det_dn_i = 1, det_ex_unq_dn_nb
   
       col_dn  = det_ex_unq_dn_orb_1st_pos (det_dn_i)
       orb_dn  = det_ex_unq_dn_orb_2nd_lab (det_dn_i)
       iwdet_ref = iwdet_ex_ref_dn (det_dn_i)
   
       do l = 1, ndn
        ratio_dn (l) = 0.d0
        do i = 1, ndn
          ratio_dn (l) = ratio_dn (l) + slater_mat_trans_inv_dn (i, l, iwdet_ref) * orb (nup + i, orb_dn)
        enddo  ! i
       enddo ! l
   
       factor_dn_inv = 1.d0/ratio_dn (col_dn)
       ratio_dn (col_dn) = ratio_dn (col_dn) - 1.d0
   
       do k = 1, ndn
         do j = 1 , ndn
          slater_mat_ex_trans_inv_dn (k, j, det_dn_i) = slater_mat_trans_inv_dn (k, j, iwdet_ref) &
                                  - slater_mat_trans_inv_dn (k, col_dn, iwdet_ref) * ratio_dn(j)* factor_dn_inv
         enddo ! j
       enddo ! k
   
    enddo ! det_dn_i
    call release ('ratio_dn', ratio_dn)

! calculate inverses from scratch
  else

!  spin-up
   call alloc ('slater_mat_ex_trans_up', slater_mat_ex_trans_up, nup, nup)
   call alloc ('mat_flat_up', mat_flat_up, nup*nup)
   do det_i = 1, det_ex_unq_up_nb
     do orb_i = 1, nup
       do elec_i = 1, nup
         slater_mat_ex_trans_up (orb_i, elec_i) = orb (elec_i, det_ex_unq_orb_lab_up (orb_i, det_i))
       enddo
     enddo
    call flatten (mat_flat_up (:), slater_mat_ex_trans_up (:,:), nup, nup)
    call matinv (mat_flat_up (:), nup, det_ex_unq_up(det_i))
    call unflatten (mat_flat_up (:), slater_mat_ex_trans_inv_up (:,:, det_i), nup, nup)
   enddo
   call release ('slater_mat_ex_trans_up', slater_mat_ex_trans_up)
   call release ('mat_flat_up', mat_flat_up)
   call object_modified_by_index (detu_index)
   
!  spin-dn
   call alloc ('slater_mat_ex_trans_dn', slater_mat_ex_trans_dn, ndn, ndn)
   call alloc ('mat_flat_dn', mat_flat_dn, ndn*ndn)
   do det_i = 1, det_ex_unq_dn_nb
     do orb_i = 1, ndn
       do elec_i = 1, ndn
         slater_mat_ex_trans_dn (orb_i, elec_i) = orb (nup + elec_i, det_ex_unq_orb_lab_dn (orb_i, det_i))
       enddo
     enddo
    call flatten (mat_flat_dn (:), slater_mat_ex_trans_dn (:,:), ndn, ndn)
    call matinv (mat_flat_dn (:), ndn, det_ex_unq_dn(det_i))
    call unflatten (mat_flat_dn (:), slater_mat_ex_trans_inv_dn (:,:, det_i), ndn, ndn)
   enddo
   call release ('slater_mat_ex_trans_dn', slater_mat_ex_trans_dn)
   call release ('mat_flat_dn', mat_flat_dn)
   call object_modified_by_index (detd_index)

  endif

  do det_i = 1, det_ex_unq_up_nb
  !write(6,*) '>slater_mat_ex_trans_inv_up',slater_mat_ex_trans_inv_up(:,:,det_i)
  enddo
  do det_i = 1, det_ex_unq_dn_nb
  !write(6,*) '>slater_mat_ex_trans_inv_dn',slater_mat_ex_trans_inv_dn(:,:,det_i)
  enddo

  end subroutine slater_mat_ex_trans_inv_bld

! ==============================================================================
  subroutine slater_mat_ex_trans_inv_2_bld
! ------------------------------------------------------------------------------
! Description   :  Inverse of transpose of Slater matrix corresponding to
! Description   :  to excited determinants using the Sherman-Morison formula
! Description   :  Bad implementation with O(nup^3 + ndn^3) scaling
!
! Created       : J. Toulouse, 27 Oct 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer det_up_i, det_dn_i
  integer iwdet_ref
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
   call object_needed ('det_ex_unq_up_orb_1st_pos')
   call object_needed ('det_ex_unq_dn_orb_1st_pos')
   call object_needed ('det_ex_unq_up_orb_2nd_lab')
   call object_needed ('det_ex_unq_dn_orb_2nd_lab')
   call object_needed ('iwdet_ex_ref_up')
   call object_needed ('iwdet_ex_ref_dn')
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
      iwdet_ref = iwdet_ex_ref_up (det_up_i)

      factor_up = 0.d0
      do i = 1, nup
       factor_up = factor_up + slater_mat_trans_inv_up (i, col_up, iwdet_ref) * orb (i, orb_up)
      enddo

      i = col_up

      do k = 1, nup
        do j = 1 , nup

         sum_up = 0.d0
         do l = 1, nup
          sum_up = sum_up + slater_mat_trans_inv_up (l, j, iwdet_ref) * orb (l, orb_up)
         enddo

         if (i == j) then
          sum_up = sum_up - 1.d0
         endif

         slater_mat_ex_trans_inv_up_2 (k, j, det_up_i) = slater_mat_trans_inv_up (k, j, iwdet_ref) &
                                 - slater_mat_trans_inv_up (k, i, iwdet_ref) * sum_up/factor_up

        enddo ! j
      enddo ! k

   enddo ! det_up_i

! spin down
  do det_dn_i = 1, det_ex_unq_dn_nb

      col_dn  = det_ex_unq_dn_orb_1st_pos (det_dn_i)
      orb_dn  = det_ex_unq_dn_orb_2nd_lab (det_dn_i)
      iwdet_ref = iwdet_ex_ref_dn (det_dn_i)

      factor_dn = 0.d0
      do i = 1, ndn
       factor_dn = factor_dn + slater_mat_trans_inv_dn (i, col_dn, iwdet_ref) * orb (nup + i, orb_dn)
      enddo

      i = col_dn

      do k = 1, ndn
        do j = 1 , ndn

         sum_dn = 0.d0
         do l = 1, ndn
          sum_dn = sum_dn + slater_mat_trans_inv_dn (l, j, iwdet_ref) * orb (nup + l, orb_dn)
         enddo

         if (i == j) then
          sum_dn = sum_dn - 1.d0
         endif

         slater_mat_ex_trans_inv_dn_2 (k, j, det_dn_i) = slater_mat_trans_inv_dn (k, j, iwdet_ref) &
                                 - slater_mat_trans_inv_dn (k, i, iwdet_ref) * sum_dn/factor_dn

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
  subroutine grd_det_ex_unq_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of unique spin-up and down excited determinants
!
! Created       : J. Toulouse, 08 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

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
   call object_needed ('dorb')
   call object_needed ('det_ex_unq_orb_lab_up')
   call object_needed ('det_ex_unq_orb_lab_dn')
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
        grd_det_ex_unq_up (dim_i, i, det_i) = grd_det_ex_unq_up (dim_i, i, det_i) * det_ex_unq_up(det_i)
        !write(6,*) '>grd_det_ex_unq_ex_up',grd_det_ex_unq_up(dim_i, i, det_i)

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
        grd_det_ex_unq_dn (dim_i, i, det_i) = grd_det_ex_unq_dn (dim_i, i, det_i) * det_ex_unq_dn(det_i)
        !write(6,*) '>grd_det_ex_unq_ex_dn',grd_det_ex_unq_dn(dim_i, i, det_i)
      enddo
    enddo ! dim_i
  enddo ! det_i

  end subroutine grd_det_ex_unq_bld

! ==============================================================================
  subroutine lap_det_ex_unq_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of unique spin-up and down excited determinants
!
! Created       : J. Toulouse, 08 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer det_i, i, j

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
   call object_needed ('ddorb')
   call object_needed ('det_ex_unq_orb_lab_up')
   call object_needed ('det_ex_unq_orb_lab_dn')
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
      lap_det_ex_unq_up (i, det_i) = lap_det_ex_unq_up (i, det_i) * det_ex_unq_up(det_i)
      !write(6,*) '>lap_det_ex_unq_up',lap_det_ex_unq_up(i,det_i)
    enddo
  enddo  ! det_i

! spin down determinants
  do det_i = 1, det_ex_unq_dn_nb
    do i = 1, ndn
      do j = 1, ndn
        lap_det_ex_unq_dn (i, det_i) = lap_det_ex_unq_dn (i, det_i)  +  &
         slater_mat_ex_trans_inv_dn (i, j, det_i) * ddorb (nup + i, det_ex_unq_orb_lab_dn (j, det_i))
      enddo
      lap_det_ex_unq_dn (i, det_i) = lap_det_ex_unq_dn (i, det_i)  * det_ex_unq_dn(det_i)
      !write(6,*) '>lap_det_ex_unq_dn',lap_det_ex_unq_dn(i,det_i)
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
  use all_modules_mod
  implicit none

! local
  integer det_i
  integer ex_i
  integer iwdet
  integer sgn
  integer csf_i, det_in_csf_i
  integer det_unq_up_i, det_unq_dn_i

! header
  if (header_exe) then

   call object_create ('grd_det_ex_up')
   call object_create ('grd_det_ex_dn')

   call object_needed ('ndim')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('single_ex_nb')
   call object_needed ('ndet')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('iwdet_ex_up')
   call object_needed ('iwdet_ex_dn')
   call object_needed ('det_ex_unq_sgn_up')
   call object_needed ('det_ex_unq_sgn_dn')
   call object_needed ('grd_det_unq_up')
   call object_needed ('grd_det_unq_dn')
   call object_needed ('grd_det_ex_unq_up')
   call object_needed ('grd_det_ex_unq_dn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')

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
        iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
        sgn = det_ex_unq_sgn_up (ex_i, det_unq_up_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetup)) then
          grd_det_ex_up (:, :, ex_i, det_i) = sgn * grd_det_unq_up (:,:,iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetup+det_ex_unq_up_nb)) then
          grd_det_ex_up (:, :, ex_i, det_i) = sgn * grd_det_ex_unq_up (:,:,iwdet-ndetup)
        endif

!       spin-down excited determinants:
        iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
        sgn = det_ex_unq_sgn_dn (ex_i, det_unq_dn_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetdn)) then
          grd_det_ex_dn (:, :, ex_i, det_i) = sgn * grd_det_unq_dn (:,:,iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetdn+det_ex_unq_dn_nb)) then
          grd_det_ex_dn (:, :, ex_i, det_i) = sgn * grd_det_ex_unq_dn (:,:,iwdet-ndetdn)
        endif

     enddo ! det_in_csf_i
   enddo ! csf_i

   !write(6,*) '>grd_det_ex_up',grd_det_ex_up(:,:,ex_i,:)
   !write(6,*) '>grd_det_ex_dn',grd_det_ex_dn(:,:,ex_i,:)

 enddo ! ex_i

 end subroutine grd_det_ex_bld

! ==============================================================================
  subroutine lap_det_ex_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of all excited determinants
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer det_i
  integer csf_i, det_in_csf_i
  integer ex_i, iwdet, sgn
  integer det_unq_up_i, det_unq_dn_i

! header
  if (header_exe) then

   call object_create ('lap_det_ex_up')
   call object_create ('lap_det_ex_dn')

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('single_ex_nb')
   call object_needed ('ndet')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('iwdet_ex_up')
   call object_needed ('iwdet_ex_dn')
   call object_needed ('det_ex_unq_sgn_up')
   call object_needed ('det_ex_unq_sgn_dn')
   call object_needed ('lap_det_unq_up')
   call object_needed ('lap_det_unq_dn')
   call object_needed ('lap_det_ex_unq_up')
   call object_needed ('lap_det_ex_unq_dn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')

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
        iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
        sgn = det_ex_unq_sgn_up (ex_i, det_unq_up_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetup)) then
          lap_det_ex_up (:, ex_i, det_i) = sgn * lap_det_unq_up (:,iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetup+det_ex_unq_up_nb)) then
          lap_det_ex_up (:, ex_i, det_i) = sgn * lap_det_ex_unq_up (:,iwdet-ndetup)
        endif

!       spin-down excited determinants:
        iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
        sgn = det_ex_unq_sgn_dn (ex_i, det_unq_dn_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetdn)) then
          lap_det_ex_dn (:, ex_i, det_i) = sgn * lap_det_unq_dn (:,iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetdn+det_ex_unq_dn_nb)) then
          lap_det_ex_dn (:, ex_i, det_i) = sgn * lap_det_ex_unq_dn (:,iwdet-ndetdn)
        endif

     enddo ! det_in_csf_i
   enddo ! csf_i

   !write(6,*) '>lap_det_ex_up',lap_det_ex_up(:,ex_i,:)
   !write(6,*) '>lap_det_ex_dn',lap_det_ex_dn(:,ex_i,:)

 enddo ! ex_i

 end subroutine lap_det_ex_bld

! ==============================================================================
  subroutine grd_psid_ex_over_psid_bld
! ------------------------------------------------------------------------------
! Description   : (Gradient psid excited)/psid excited
!
! Created       : J. Toulouse, 08 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer ex_i, ex_rev_i
  integer dorb_i, det_i, dim_i, det_unq_up_i, det_unq_dn_i
  integer elec_i, elec_up_i, elec_dn_i
  integer csf_i, det_in_csf_i
  real(dp) coefficient, grd_det

! header
  if (header_exe) then

   call object_create ('grd_psid_ex_over_psid')

   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_ind_rev')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('csf_coef')
   call object_needed ('cdet_in_csf')
   call object_needed ('nup')
   call object_needed ('ndn')
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

        enddo ! elec_up

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
        enddo ! elec_dn

       enddo ! dim_i


     enddo ! det_in_csf_i
   enddo ! csf_i

   grd_psid_ex_over_psid (:,:,dorb_i) = grd_psid_ex_over_psid (:,:,dorb_i)/psid_ex (dorb_i)

   !write(6,*) '>grd_psid_ex_over_psid',grd_psid_ex_over_psid(:,:,dorb_i)
 enddo ! dorb_i

 end subroutine grd_psid_ex_over_psid_bld

! ==============================================================================
  subroutine lap_psid_ex_over_psid_bld
! ------------------------------------------------------------------------------
! Description   : (Laplacian psid excited)/psid excited
!
! Created       : J. Toulouse, 08 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer ex_i, ex_rev_i
  integer dorb_i, det_i, det_unq_up_i, det_unq_dn_i
  integer elec_i, elec_up_i, elec_dn_i
  integer csf_i, det_in_csf_i
  real(dp) coefficient, lap_det

! header
  if (header_exe) then

   call object_create ('lap_psid_ex_over_psid')

   call object_needed ('nelec')
   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_ind_rev')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('csf_coef')
   call object_needed ('cdet_in_csf')
   call object_needed ('nup')
   call object_needed ('ndn')
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

   !write(6,*) '>lap_psid_ex_over_psid',lap_psid_ex_over_psid(:,dorb_i)
 enddo ! dorb_i


 end subroutine lap_psid_ex_over_psid_bld

! ==============================================================================
  subroutine lap_lnpsid_ex_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian ( ln (psid_ex) )
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer dim_i, elec_i, ex_i
  real(dp) sum

! header
  if (header_exe) then

   call object_create ('lap_lnpsid_ex')

   call object_needed ('nelec')
   call object_needed ('param_orb_nb')
   call object_needed ('ndim')
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
    !write(6,*) '>lap_lnpsid_ex',lap_lnpsid_ex(:,ex_i)
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
  use all_modules_mod
  implicit none

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
   !write(6,*) '>sum_lap_lnpsid_ex',sum_lap_lnpsid_ex(ex_i)
  enddo

 end subroutine sum_lap_lnpsid_ex_bld

! ==============================================================================
  subroutine grd_psi_ex_over_psi_bld
! ------------------------------------------------------------------------------
! Description   : (Gradient excited psi)/ excited psi
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer ex_i

! header
  if (header_exe) then

   call object_create ('grd_psi_ex_over_psi')

   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('param_orb_nb')
   call object_needed ('grd_psid_ex_over_psid')
   call object_needed ('vj')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_psi_ex_over_psi', grd_psi_ex_over_psi, ndim, nelec, param_orb_nb)

  do ex_i = 1, param_orb_nb
     grd_psi_ex_over_psi (1:ndim, 1:nelec, ex_i) = grd_psid_ex_over_psid (1:ndim, 1:nelec, ex_i) + vj (1:ndim, 1:nelec)
     !write(6,*) '>grd_psi_ex_over_psi',grd_psi_ex_over_psi (1:ndim, 1:nelec, ex_i)
  enddo

 end subroutine grd_psi_ex_over_psi_bld

! ==============================================================================
  subroutine sum_lap_lnpsi_ex_bld
! ------------------------------------------------------------------------------
! Description   : Sum_i lap_lnpsi (i) : determinant part + Jastrow
!
! Created       : J. Toulouse, 09 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

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
  !write(6,*) '>sum_lap_lnpsi_ex',sum_lap_lnpsi_ex

 end subroutine sum_lap_lnpsi_ex_bld

! ==============================================================================
  subroutine eloc_kin_ex_bld
! ------------------------------------------------------------------------------
! Description   : Kinetic local energy of excited wave function
!
! Created       : J. Toulouse, 05 Dec 2005
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer dim_i, elec_i, ex_i
  real(dp) sum_grd_psi_over_psi_square

! header
  if (header_exe) then

   call object_create ('eloc_kin_ex')

   call object_needed ('param_orb_nb')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('grd_psi_ex_over_psi')
   call object_needed ('sum_lap_lnpsi_ex')

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
  !write(6,*) '>eloc_kin_ex',eloc_kin_ex (ex_i)

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
  use all_modules_mod
  implicit none

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
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('csf_coef')
   call object_needed ('cdet_in_csf')
   call object_needed ('electron')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('iwdet_ex_up')
   call object_needed ('iwdet_ex_dn')
   call object_needed ('det_ex_unq_sgn_up')
   call object_needed ('det_ex_unq_sgn_dn')
   call object_needed ('detn')
   call object_needed ('orbe')
   call object_needed ('slater_mat_ex_trans_inv_up')
   call object_needed ('slater_mat_ex_trans_inv_dn')
   call object_needed ('detu')
   call object_needed ('detd')
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
        det_ex_up_in_x = 0.d0
        det_ex_rev_up_in_x = 0.d0
        det_ex_dn_in_x = 0.d0
        det_ex_rev_dn_in_x = 0.d0

        iwdet = iwdet_ex_up (ex_i, det_unq_up_i)
        sgn = det_ex_unq_sgn_up (ex_i, det_unq_up_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetup)) then
         det_ex_up_in_x = sgn * detn (iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetup+det_ex_unq_up_nb)) then

         factor_up = 0.d0
         do j = 1, nup
         factor_up = factor_up + slater_mat_ex_trans_inv_up (iel_up, j, iwdet-ndetup) * orbe (det_ex_unq_orb_lab_up (j, iwdet-ndetup))
         enddo

         det_ex_up_in_x  = sgn * factor_up * det_ex_unq_up(iwdet-ndetup)

        endif

         psid_ex_in_x (dorb_i) = psid_ex_in_x (dorb_i) + coefficient * (det_ex_up_in_x * detd (det_unq_dn_i) + detn(det_unq_up_i) * det_ex_dn (ex_i, det_i))


!       reverse excitation for active-active non-CASSCF excitations
        if (.not. l_casscf .and. ex_rev_i /=0 ) then

        iwdet = iwdet_ex_up (ex_rev_i, det_unq_up_i)
        sgn = det_ex_unq_sgn_up (ex_rev_i, det_unq_up_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetup)) then
         det_ex_rev_up_in_x = sgn * detn (iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetup+det_ex_unq_up_nb)) then

         factor_up = 0.d0
         do j = 1, nup
         factor_up = factor_up +  &
          slater_mat_ex_trans_inv_up (iel_up, j, iwdet-ndetup) * orbe (det_ex_unq_orb_lab_up (j, iwdet-ndetup))
         enddo

         det_ex_rev_up_in_x  = sgn * factor_up * det_ex_unq_up(iwdet-ndetup)

        endif

         psid_ex_in_x (dorb_i) = psid_ex_in_x (dorb_i) - coefficient *   &
                   (det_ex_rev_up_in_x * detd (det_unq_dn_i) + detn(det_unq_up_i) * det_ex_dn (ex_rev_i, det_i))

        endif ! if .not. l_casscf


!       electron is spin down
        else
        iel_dn = electron - nup


        iwdet = iwdet_ex_dn (ex_i, det_unq_dn_i)
        sgn = det_ex_unq_sgn_dn (ex_i, det_unq_dn_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetdn)) then
         det_ex_dn_in_x = sgn * detn (iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetdn+det_ex_unq_dn_nb)) then

         factor_dn = 0.d0
         do j = 1, ndn
         factor_dn = factor_dn +  &
          slater_mat_ex_trans_inv_dn (iel_dn, j, iwdet-ndetdn) * orbe (det_ex_unq_orb_lab_dn (j, iwdet-ndetdn))
         enddo

         det_ex_dn_in_x  = sgn * factor_dn * det_ex_unq_dn(iwdet-ndetdn)

        endif

        psid_ex_in_x (dorb_i) = psid_ex_in_x (dorb_i) + coefficient * (det_ex_up (ex_i, det_i) * detn (det_unq_dn_i) + detu(det_unq_up_i) * det_ex_dn_in_x)

!       reverse excitation for active-active non-CASSCF excitations
        if (.not. l_casscf .and. ex_rev_i /=0 ) then

        iwdet = iwdet_ex_dn (ex_rev_i, det_unq_dn_i)
        sgn = det_ex_unq_sgn_dn (ex_rev_i, det_unq_dn_i)
        if ((iwdet.ne.0).and.(iwdet.le.ndetdn)) then
         det_ex_rev_dn_in_x = sgn * detn (iwdet)
        elseif ((iwdet.ne.0).and.(iwdet.le.ndetdn+det_ex_unq_dn_nb)) then

         factor_dn = 0.d0
         do j = 1, ndn
         factor_dn = factor_dn + slater_mat_ex_trans_inv_dn (iel_dn, j, iwdet-ndetdn) * orbe (det_ex_unq_orb_lab_dn (j, iwdet-ndetdn))
         enddo

         det_ex_rev_dn_in_x  = sgn * factor_dn * det_ex_unq_dn(iwdet-ndetdn)

        endif

        psid_ex_in_x (dorb_i) = psid_ex_in_x (dorb_i) - coefficient * (det_ex_up (ex_rev_i, det_i) * detn (det_unq_dn_i) + detu(det_unq_up_i) * det_ex_rev_dn_in_x)

        endif ! if .not. l_casscf

        endif ! electron

     enddo ! det_in_csf_i
   enddo ! csf_i

 enddo ! dorb_i

 !write(6,*) '>psid_ex_in_x',psid_ex_in_x

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
  use all_modules_mod
  implicit none

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
  use all_modules_mod
  implicit none


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
  use all_modules_mod
  implicit none

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
  !write(6,*) '>eloc_ex',eloc_ex

 end subroutine eloc_ex_bld

! ==============================================================================
  subroutine deloc_orb_bld
! ------------------------------------------------------------------------------
! Description   :  derivative of local energy wrt orbital rotation parameters kappa
! Description   :  dEloc/dkappa = [ (H dPsi/dkappa)/(dPsi/dkappa) - Eloc ] dlnPsi/dkappa
!
! Created       : J. Toulouse, 14 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

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
   call object_needed ('ex_orb_1st_lab')
   call object_needed ('ex_orb_2nd_lab')
   call object_needed ('orb_energies')
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

! ==============================================================================
  subroutine double_ex_det_bld
! ------------------------------------------------------------------------------
! Description   : 
!
! Created       : B. Mussard, July 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, ex_i, ex_j, ex_ij, orbi_1st, orbi_2nd, orbj_1st, orbj_2nd
  integer det_unq_up_i, det_unq_up_k, det_unq_dn_i, det_unq_dn_k
  integer, allocatable :: det_orb_lab_up (:)
  integer, allocatable :: det_orb_lab_dn (:)
  integer, allocatable :: det_orb_lab_srt_up (:)
  integer, allocatable :: det_orb_lab_srt_dn (:)
  integer, allocatable :: det_ex_to_calc_sgn_up (:)
  integer, allocatable :: det_ex_to_calc_sgn_dn (:)
  integer :: det_ex2_max

! header
  if (header_exe) then

   call object_create ('double_ex_nb')
   call object_create ('det_ex2_unq_up_nb')
   call object_create ('det_ex2_unq_dn_nb')
   call object_create ('det_ex2_unq_up_orb_info')
   call object_create ('det_ex2_unq_dn_orb_info')

   call object_needed ('single_ex_nb')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('ex_orb_1st_lab')
   call object_needed ('ex_orb_2nd_lab')
   call object_needed ('orb_occ_in_det_unq_up')
   call object_needed ('orb_occ_in_det_unq_dn')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('orb_pos_in_det_unq_up')
   call object_needed ('orb_pos_in_det_unq_dn')
   call object_needed ('det_ex_unq_sgn_up')
   call object_needed ('det_ex_unq_sgn_dn')
   call object_needed ('iwdet_ex_up')
   call object_needed ('iwdet_ex_dn')

   return

  endif

! begin

  double_ex_nb=single_ex_nb*single_ex_nb
  det_ex2_max =(ndetup+ndetdn)*double_ex_nb

! allocations
  call object_alloc ('iwdet_ex_up', iwdet_ex_up, single_ex_nb+double_ex_nb, ndetup)
  call object_alloc ('iwdet_ex_dn', iwdet_ex_dn, single_ex_nb+double_ex_nb, ndetdn)
  call object_alloc ('det_ex_unq_sgn_up', det_ex_unq_sgn_up, single_ex_nb+double_ex_nb, ndetup)
  call object_alloc ('det_ex_unq_sgn_dn', det_ex_unq_sgn_dn, single_ex_nb+double_ex_nb, ndetdn)
  det_ex2_unq_up_nb = 0
  det_ex2_unq_dn_nb = 0

! these arrays will be adjusted at the end of the routine
  call object_alloc ('iwdet_ex_ref_up', iwdet_ex_ref_up, det_ex_unq_up_nb+det_ex2_max)
  call object_alloc ('iwdet_ex_ref_dn', iwdet_ex_ref_dn, det_ex_unq_dn_nb+det_ex2_max)
  call object_alloc ('det_ex_unq_orb_lab_up', det_ex_unq_orb_lab_up, nup, det_ex_unq_up_nb+det_ex2_max)
  call object_alloc ('det_ex_unq_orb_lab_dn', det_ex_unq_orb_lab_dn, ndn, det_ex_unq_up_nb+det_ex2_max)
  call object_alloc ('det_unq_orb_lab_srt_up', det_unq_orb_lab_srt_up, nup, ndetup+det_ex_unq_up_nb+det_ex2_max)
  call object_alloc ('det_unq_orb_lab_srt_dn', det_unq_orb_lab_srt_dn, ndn, ndetdn+det_ex_unq_dn_nb+det_ex2_max)
  call object_alloc ('det_ex2_unq_up_orb_info', det_ex2_unq_up_orb_info, det_ex2_max,4)
  call object_alloc ('det_ex2_unq_dn_orb_info', det_ex2_unq_dn_orb_info, det_ex2_max,4)

! local arrays (not objects)
  call alloc ('det_orb_lab_up', det_orb_lab_up, nup)
  call alloc ('det_orb_lab_dn', det_orb_lab_dn, ndn)
  call alloc ('det_orb_lab_srt_up', det_orb_lab_srt_up, nup)
  call alloc ('det_orb_lab_srt_dn', det_orb_lab_srt_dn, ndn)
  call alloc ('det_ex_to_calc_sgn_up',  det_ex_to_calc_sgn_up,  det_ex2_max)
  call alloc ('det_ex_to_calc_sgn_dn',  det_ex_to_calc_sgn_dn,  det_ex2_max)

  ! double loop over single orbital excitations
  ex_ij=0
  do ex_i = 1, single_ex_nb
    orbi_1st = ex_orb_1st_lab (ex_i)
    orbi_2nd = ex_orb_2nd_lab (ex_i)
    do ex_j = 1, single_ex_nb
      orbj_1st = ex_orb_1st_lab (ex_j)
      orbj_2nd = ex_orb_2nd_lab (ex_j)
      ex_ij=ex_ij+1

      if ((orbi_1st.eq.orbj_1st) .or. (orbi_2nd.eq.orbj_2nd)) then
        cycle
      endif

!     loop over unique spin-up determinants
      do det_unq_up_i = 1, ndetup
      if (           orb_occ_in_det_unq_up (orbi_1st, det_unq_up_i) &
        &.and.       orb_occ_in_det_unq_up (orbj_1st, det_unq_up_i) &
        &.and. .not. orb_occ_in_det_unq_up (orbi_2nd, det_unq_up_i) &
        &.and. .not. orb_occ_in_det_unq_up (orbj_2nd, det_unq_up_i)) then

!       build current excited determinant
        det_orb_lab_up = det_unq_orb_lab_srt_up (:, det_unq_up_i)
        call replace_elt_in_array (det_orb_lab_up, orbi_1st, orbi_2nd)
        call replace_elt_in_array (det_orb_lab_up, orbj_1st, orbj_2nd)
        det_orb_lab_srt_up = det_orb_lab_up
        call sort_and_sign (det_orb_lab_srt_up, det_ex_unq_sgn_up (single_ex_nb+ex_ij, det_unq_up_i))

!       check if current excited determinant is an already known determinant
        do det_unq_up_k = 1, ndetup+det_ex_unq_up_nb+det_ex2_unq_up_nb
         if (arrays_equal (det_orb_lab_srt_up, det_unq_orb_lab_srt_up (:, det_unq_up_k))) then
           iwdet_ex_up (single_ex_nb+ex_ij, det_unq_up_i) = det_unq_up_k
           if (det_unq_up_k.gt.ndetup+det_ex_unq_up_nb) then
             det_ex_unq_sgn_up (single_ex_nb+ex_ij, det_unq_up_i)=det_ex_unq_sgn_up (single_ex_nb+ex_ij,det_unq_up_i)*det_ex_to_calc_sgn_up(det_unq_up_k-ndetup-det_ex_unq_up_nb)
           endif
           exit
         endif
        enddo

!       if current excited determinant is a new determinant, add it to the list of unique determinants
        if (iwdet_ex_up (single_ex_nb+ex_ij, det_unq_up_i) == 0) then
          det_ex2_unq_up_nb = det_ex2_unq_up_nb + 1
          iwdet_ex_up (single_ex_nb+ex_ij, det_unq_up_i) = ndetup+det_ex_unq_up_nb+det_ex2_unq_up_nb
          iwdet_ex_ref_up (det_ex_unq_up_nb+det_ex2_unq_up_nb) = det_unq_up_i
          det_ex_unq_orb_lab_up (:, det_ex_unq_up_nb+det_ex2_unq_up_nb ) = det_orb_lab_up
          det_unq_orb_lab_srt_up (:, ndetup+det_ex_unq_up_nb+det_ex2_unq_up_nb ) = det_orb_lab_srt_up
          det_ex_to_calc_sgn_up (det_ex2_unq_up_nb) = det_ex_unq_sgn_up (single_ex_nb+ex_ij,det_unq_up_i)
          det_ex_unq_sgn_up (single_ex_nb+ex_ij, det_unq_up_i)=det_ex_unq_sgn_up (single_ex_nb+ex_ij,det_unq_up_i)*det_ex_to_calc_sgn_up(det_ex2_unq_up_nb)
          det_ex2_unq_up_orb_info(det_ex2_unq_up_nb,1) = orbi_1st!orb_pos_in_det_unq_up (orbi_1st, det_unq_up_i)
          det_ex2_unq_up_orb_info(det_ex2_unq_up_nb,2) = orbj_1st!orb_pos_in_det_unq_up (orbj_1st, det_unq_up_i)
          det_ex2_unq_up_orb_info(det_ex2_unq_up_nb,3) = orbi_2nd
          det_ex2_unq_up_orb_info(det_ex2_unq_up_nb,4) = orbj_2nd
        endif

      endif

      enddo ! det_unq_up_i

!     loop over unique spin-down determinants
      do det_unq_dn_i = 1, ndetdn

      if (orb_occ_in_det_unq_dn (orbi_1st, det_unq_dn_i) &
        &.and. orb_occ_in_det_unq_dn (orbj_1st, det_unq_dn_i) &
        &.and. .not. orb_occ_in_det_unq_dn (orbi_2nd, det_unq_dn_i) &
        &.and. .not. orb_occ_in_det_unq_dn (orbj_2nd, det_unq_dn_i)) then

!       build current excited determinant
        det_orb_lab_dn = det_unq_orb_lab_srt_dn (:, det_unq_dn_i)
        call replace_elt_in_array (det_orb_lab_dn, orbi_1st, orbi_2nd)
        call replace_elt_in_array (det_orb_lab_dn, orbj_1st, orbj_2nd)
        det_orb_lab_srt_dn = det_orb_lab_dn
        call sort_and_sign (det_orb_lab_srt_dn, det_ex_unq_sgn_dn (single_ex_nb+ex_ij, det_unq_dn_i))

!       check if current excited determinant is an already known determinant
        do det_unq_dn_k = 1, ndetdn+det_ex_unq_dn_nb+det_ex2_unq_dn_nb
         if (arrays_equal (det_orb_lab_srt_dn, det_unq_orb_lab_srt_dn (:, det_unq_dn_k))) then
           iwdet_ex_dn (single_ex_nb+ex_ij, det_unq_dn_i) = det_unq_dn_k
           if (det_unq_dn_k.gt.ndetdn+det_ex_unq_dn_nb) then
             det_ex_unq_sgn_dn (single_ex_nb+ex_ij, det_unq_dn_i)=det_ex_unq_sgn_dn (single_ex_nb+ex_ij,det_unq_dn_i)*det_ex_to_calc_sgn_dn(det_unq_dn_k-ndetdn-det_ex_unq_dn_nb)
           endif
           exit
         endif
        enddo

!       if current excited determinant is a new determinant, add it to the list of excited determinants
        if (iwdet_ex_dn (single_ex_nb+ex_ij, det_unq_dn_i) == 0) then
          det_ex2_unq_dn_nb = det_ex2_unq_dn_nb + 1
          iwdet_ex_dn (single_ex_nb+ex_ij, det_unq_dn_i) = ndetdn+det_ex_unq_dn_nb+det_ex2_unq_dn_nb
          iwdet_ex_ref_dn (det_ex_unq_dn_nb+det_ex2_unq_dn_nb) = det_unq_dn_i
          det_ex_unq_orb_lab_dn (:, det_ex_unq_dn_nb+det_ex2_unq_dn_nb) = det_orb_lab_dn
          det_unq_orb_lab_srt_dn (:, ndetdn+det_ex_unq_dn_nb+det_ex2_unq_dn_nb) = det_orb_lab_srt_dn
          det_ex_to_calc_sgn_dn (det_ex2_unq_dn_nb) = det_ex_unq_sgn_dn (single_ex_nb+ex_ij,det_unq_dn_i)
          det_ex_unq_sgn_dn (single_ex_nb+ex_ij, det_unq_dn_i)=det_ex_unq_sgn_dn (single_ex_nb+ex_ij,det_unq_dn_i)*det_ex_to_calc_sgn_dn(det_ex2_unq_dn_nb)
          det_ex2_unq_dn_orb_info (det_ex2_unq_dn_nb,1) = orbi_1st!orb_pos_in_det_unq_dn (orbi_1st, det_unq_dn_i)
          det_ex2_unq_dn_orb_info (det_ex2_unq_dn_nb,2) = orbj_1st!orb_pos_in_det_unq_dn (orbj_1st, det_unq_dn_i)
          det_ex2_unq_dn_orb_info (det_ex2_unq_dn_nb,3) = orbi_2nd
          det_ex2_unq_dn_orb_info (det_ex2_unq_dn_nb,4) = orbj_2nd
        endif

      endif

     enddo ! det_unq_dn_i

    enddo ! ex_j
  enddo ! ex_i

! adjust the arrays to their real size
  call object_alloc ('iwdet_ex_ref_up', iwdet_ex_ref_up, det_ex_unq_up_nb+det_ex2_unq_up_nb)
  call object_alloc ('iwdet_ex_ref_dn', iwdet_ex_ref_dn, det_ex_unq_dn_nb+det_ex2_unq_dn_nb)
  call object_alloc ('det_ex_unq_orb_lab_up', det_ex_unq_orb_lab_up, nup, det_ex_unq_up_nb+det_ex2_unq_up_nb)
  call object_alloc ('det_ex_unq_orb_lab_dn', det_ex_unq_orb_lab_dn, ndn, det_ex_unq_up_nb+det_ex2_unq_up_nb)
  call object_alloc ('det_unq_orb_lab_srt_up', det_unq_orb_lab_srt_up, nup, ndetup+det_ex_unq_up_nb+det_ex2_unq_up_nb)
  call object_alloc ('det_unq_orb_lab_srt_dn', det_unq_orb_lab_srt_dn, ndn, ndetup+det_ex_unq_up_nb+det_ex2_unq_dn_nb)
  call object_alloc ('det_ex2_unq_up_orb_info', det_ex2_unq_up_orb_info, det_ex2_unq_up_nb,4)
  call object_alloc ('det_ex2_unq_dn_orb_info', det_ex2_unq_dn_orb_info, det_ex2_unq_dn_nb,4)

  call release ('det_orb_lab_up', det_orb_lab_up)
  call release ('det_orb_lab_dn', det_orb_lab_dn)
  call release ('det_orb_lab_srt_up', det_orb_lab_srt_up)
  call release ('det_orb_lab_srt_dn', det_orb_lab_srt_dn)
  call release ('det_ex_to_calc_sgn_up',  det_ex_to_calc_sgn_up)
  call release ('det_ex_to_calc_sgn_dn',  det_ex_to_calc_sgn_dn)

  write(6,'(a,i10)') ' Number of unique spin-up   doubly excited determinants = ',det_ex2_unq_up_nb
  write(6,'(a,i10)') ' Number of unique spin-down doubly excited determinants = ',det_ex2_unq_dn_nb

  !write(6,*) '>det_ex2_unq_up_nb         ' ,det_ex2_unq_up_nb
  !write(6,*) '>det_ex2_unq_dn_nb         ' ,det_ex2_unq_dn_nb
  !write(6,*) '>det_ex2_unq_up_orb_info   ' ,det_ex2_unq_up_orb_info
  !write(6,*) '>det_ex2_unq_dn_orb_info   ' ,det_ex2_unq_dn_orb_info
  !write(6,*) '>det_ex_unq_sgn_up         ' ,det_ex_unq_sgn_up
  !write(6,*) '>det_ex_unq_sgn_dn         ' ,det_ex_unq_sgn_dn
  !do ex_i=1,ndetup
  !write(6,*) '>iwdet_ex_up              ' ,iwdet_ex_up(:,ex_i)
  !enddo
  !do ex_i=1,ndetdn
  !write(6,*) '>iwdet_ex_dn              ' ,iwdet_ex_dn(:,ex_i)
  !enddo
  !write(6,*) '>iwdet_ex_ref_up          ' ,iwdet_ex_ref_up
  !write(6,*) '>iwdet_ex_ref_dn          ' ,iwdet_ex_ref_dn

  end subroutine double_ex_det_bld

! ==============================================================================
  subroutine det_ex2_unq_bld
! ------------------------------------------------------------------------------
! Description   :  
!
! Created       : B. Mussard, July 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer det_ex2_unq_up_i, det_ex2_unq_dn_i
  integer p,q,r,s,pq,rs,ps,rq
  integer iwdet_ref,ex_i,done
  real(8) :: det_pq,det_rs,det_ps,det_rq

  integer :: i
  real(8),allocatable :: p_mat(:,:),q_mat(:,:),work(:,:)

! header
  if (header_exe) then

   call object_create ('det_ex2_unq_up')
   call object_create ('det_ex2_unq_dn')

   call object_needed ('det_ex2_unq_up_nb')
   call object_needed ('det_ex2_unq_dn_nb')
   call object_needed ('det_ex2_unq_up_orb_info')
   call object_needed ('det_ex2_unq_dn_orb_info')
   call object_needed ('iwdet_ex_ref_up')
   call object_needed ('iwdet_ex_ref_dn')
   call object_needed ('single_ex_nb')
   call object_needed ('ex_orb_1st_lab')
   call object_needed ('ex_orb_2nd_lab')
   call object_needed ('det_ex_up')
   call object_needed ('det_ex_dn')
   call object_needed ('detu')
   call object_needed ('detd')

   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('orb')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('slater_mat_ex_trans_inv_up')
   call object_needed ('slater_mat_ex_trans_inv_dn')

   return

  endif

! begin

! allocations
  call alloc ('det_ex2_unq_up', det_ex2_unq_up, det_ex2_unq_up_nb)
  call alloc ('det_ex2_unq_dn', det_ex2_unq_dn, det_ex2_unq_dn_nb)

! spin up determinants
  do det_ex2_unq_up_i = 1, det_ex2_unq_up_nb
    p=det_ex2_unq_up_orb_info(det_ex2_unq_up_i,1)
    r=det_ex2_unq_up_orb_info(det_ex2_unq_up_i,2)
    q=det_ex2_unq_up_orb_info(det_ex2_unq_up_i,3)
    s=det_ex2_unq_up_orb_info(det_ex2_unq_up_i,4)
    iwdet_ref = iwdet_ex_ref_up (det_ex_unq_up_nb+det_ex2_unq_up_i)

    if (itest.eq.0) then
      done=0
      do ex_i=1,single_ex_nb
        if      (ex_orb_1st_lab(ex_i).eq.p &
            .and.ex_orb_2nd_lab(ex_i).eq.q) then
          det_pq=det_ex_up(ex_i,iwdet_ref)
          done=done+1
        elseif  (ex_orb_1st_lab(ex_i).eq.r &
            .and.ex_orb_2nd_lab(ex_i).eq.s) then
          det_rs=det_ex_up(ex_i,iwdet_ref)
          done=done+1
        elseif  (ex_orb_1st_lab(ex_i).eq.p &
            .and.ex_orb_2nd_lab(ex_i).eq.s) then
          det_ps=det_ex_up(ex_i,iwdet_ref)
          done=done+1
        elseif  (ex_orb_1st_lab(ex_i).eq.r &
            .and.ex_orb_2nd_lab(ex_i).eq.q) then
          det_rq=det_ex_up(ex_i,iwdet_ref)
          done=done+1
        endif
      enddo
      if (done.ne.4) call die('done.ne.4','')
      det_ex2_unq_up(det_ex2_unq_up_i)=(det_pq*det_rs-det_ps*det_rq)/detu(iwdet_ref)

    elseif(itest.eq.1) then
      call alloc('p_mat',p_mat,2,nup)
      p_mat=0
      p_mat(1,p)=1
      p_mat(2,r)=1
      call alloc('q_mat',q_mat,nup,2)
      q_mat=0
      do i=1,nup
        q_mat(i,1)=orb(i,q)
        q_mat(i,2)=orb(i,s)
      enddo
      call alloc('work' ,work ,2,2)
      work=0
      if (iwdet_ref.le.ndetup) then
        work=matmul(matmul(p_mat,transpose(slater_mat_trans_inv_up(:,:,iwdet_ref))),q_mat)
      elseif (iwdet_ref.le.ndetup+det_ex_unq_up_nb) then
        work=matmul(matmul(p_mat,transpose(slater_mat_ex_trans_inv_up(:,:,iwdet_ref-ndetup))),q_mat)
      elseif (iwdet_ref.le.ndetup+det_ex_unq_up_nb+det_ex2_unq_up_nb) then
        call die('','')
      endif
      det_ex2_unq_up(det_ex2_unq_up_i)=(work(1,1)*work(2,2)-work(1,2)*work(2,1))*detu(iwdet_ref)

    endif
    !write(6,'(a,f20.8)') '>det_ex2_unq_up',det_ex2_unq_up(det_ex2_unq_up_i)

  enddo ! det_ex_unq_up_i

! spin down determinants
  do det_ex2_unq_dn_i = 1, det_ex2_unq_dn_nb
    p=det_ex2_unq_dn_orb_info(det_ex2_unq_dn_i,1)
    r=det_ex2_unq_dn_orb_info(det_ex2_unq_dn_i,2)
    q=det_ex2_unq_dn_orb_info(det_ex2_unq_dn_i,3)
    s=det_ex2_unq_dn_orb_info(det_ex2_unq_dn_i,4)
    iwdet_ref = iwdet_ex_ref_dn (det_ex_unq_dn_nb+det_ex2_unq_dn_i)

    if (itest.eq.0) then
      done=0
      do ex_i=1,single_ex_nb
        if      (ex_orb_1st_lab(ex_i).eq.p &
            .and.ex_orb_2nd_lab(ex_i).eq.q) then
          det_pq=det_ex_dn(ex_i,iwdet_ref)
          done=done+1
        elseif  (ex_orb_1st_lab(ex_i).eq.r &
            .and.ex_orb_2nd_lab(ex_i).eq.s) then
          det_rs=det_ex_dn(ex_i,iwdet_ref)
          done=done+1
        elseif  (ex_orb_1st_lab(ex_i).eq.p &
            .and.ex_orb_2nd_lab(ex_i).eq.s) then
          det_ps=det_ex_dn(ex_i,iwdet_ref)
          done=done+1
        elseif  (ex_orb_1st_lab(ex_i).eq.r &
            .and.ex_orb_2nd_lab(ex_i).eq.q) then
          det_rq=det_ex_dn(ex_i,iwdet_ref)
          done=done+1
        endif
      enddo
      if (done.ne.4) call die('done.ne.4','')
      det_ex2_unq_dn(det_ex2_unq_dn_i)=(det_pq*det_rs-det_ps*det_rq)/detd(iwdet_ref)

    elseif(itest.eq.1) then
      call alloc('p_mat',p_mat,2,ndn)
      p_mat=0
      p_mat(1,p)=1
      p_mat(2,r)=1
      call alloc('q_mat',q_mat,ndn,2)
      q_mat=0
      do i=1,ndn
        q_mat(i,1)=orb(nup+i,q)
        q_mat(i,2)=orb(nup+i,s)
      enddo
      call alloc('work' ,work ,2,2)
      work=0
      if (iwdet_ref.le.ndetdn) then
        work=matmul(matmul(p_mat,transpose(slater_mat_trans_inv_dn(:,:,iwdet_ref))),q_mat)
      elseif (iwdet_ref.le.ndetdn+det_ex_unq_dn_nb) then
        work=matmul(matmul(p_mat,transpose(slater_mat_ex_trans_inv_dn(:,:,iwdet_ref-ndetdn))),q_mat)
      elseif (iwdet_ref.le.ndetdn+det_ex_unq_dn_nb+det_ex2_unq_dn_nb) then
        call die('','')
      endif
      det_ex2_unq_dn(det_ex2_unq_dn_i)=(work(1,1)*work(2,2)-work(1,2)*work(2,1))*detd(iwdet_ref)

    endif
    !write(6,'(a,f20.8)') '>det_ex2_unq_dn',det_ex2_unq_dn(det_ex2_unq_dn_i)

  enddo ! det_ex_unq_dn_i

  end subroutine det_ex2_unq_bld

! ==============================================================================
  subroutine det_ex2_bld
! ------------------------------------------------------------------------------
! Description   : all doubly-excited determinants
!
! Created       : B. Mussard, July 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer ex_i,ex_j,ex_ij
  integer csf_i, det_in_csf_i, det_i, det_unq_up_i, det_unq_dn_i
  integer iwdet, sgn

! header
  if (header_exe) then

   call object_create ('det_ex2')
   call object_create ('det_ex2_up')
   call object_create ('det_ex2_dn')

   call object_needed ('double_ex_nb')
   call object_needed ('ndet')
   call object_needed ('single_ex_nb')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('iwdet_ex_up')
   call object_needed ('iwdet_ex_dn')
   call object_needed ('det_ex_unq_sgn_up')
   call object_needed ('det_ex_unq_sgn_dn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('det_ex_unq_up_nb')
   call object_needed ('det_ex_unq_dn_nb')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('det_ex_unq_up')
   call object_needed ('det_ex_unq_dn')
   call object_needed ('det_ex2_unq_up')
   call object_needed ('det_ex2_unq_dn')
   call object_needed ('det_ex_up')
   call object_needed ('det_ex_dn')

   return

  endif

! begin

! allocations
  call object_alloc ('det_ex2_up', det_ex2_up, double_ex_nb, ndet)
  call object_alloc ('det_ex2_dn', det_ex2_dn, double_ex_nb, ndet)
  call object_alloc ('det_ex2', det_ex2, double_ex_nb, ndet)

! loop over single orbital excitations
  ex_ij=0
  do ex_i = 1, single_ex_nb
    do ex_j = 1, single_ex_nb
      ex_ij=ex_ij+1
      do csf_i = 1, ncsf
        do det_in_csf_i = 1, ndet_in_csf (csf_i)
          det_i = iwdet_in_csf (det_in_csf_i, csf_i)
          det_unq_up_i = det_to_det_unq_up (det_i)
          det_unq_dn_i = det_to_det_unq_dn (det_i)

!         spin-up excited determinants:
          iwdet = iwdet_ex_up (single_ex_nb+ex_ij, det_unq_up_i)
          sgn = det_ex_unq_sgn_up (single_ex_nb+ex_ij, det_unq_up_i)
          if ((iwdet.ne.0).and.(iwdet.le.ndetup)) then
            det_ex2_up (ex_ij, det_i) = sgn * detu (iwdet)
          elseif ((iwdet.ne.0).and.(iwdet.le.ndetup+det_ex_unq_up_nb)) then
            det_ex2_up (ex_ij, det_i) = sgn * det_ex_unq_up (iwdet-ndetup)
          elseif ((iwdet.ne.0).and.(iwdet.le.ndetup+det_ex_unq_up_nb+det_ex2_unq_up_nb)) then
            det_ex2_up (ex_ij, det_i) = sgn * det_ex2_unq_up (iwdet-ndetup-det_ex_unq_up_nb)
          endif

!         spin-down excited determinants:
          iwdet = iwdet_ex_dn (single_ex_nb+ex_ij, det_unq_dn_i)
          sgn = det_ex_unq_sgn_dn (single_ex_nb+ex_ij, det_unq_dn_i)
          if ((iwdet.ne.0).and.(iwdet.le.ndetdn)) then
            det_ex2_dn (ex_ij, det_i) = sgn * detd (iwdet)
          elseif ((iwdet.ne.0).and.(iwdet.le.ndetdn+det_ex_unq_dn_nb)) then
            det_ex2_dn (ex_ij, det_i) = sgn * det_ex_unq_dn (iwdet-ndetdn)
          elseif ((iwdet.ne.0).and.(iwdet.le.ndetdn+det_ex_unq_dn_nb+det_ex2_unq_dn_nb)) then
            det_ex2_dn (ex_ij, det_i) = sgn * det_ex2_unq_dn (iwdet-ndetdn-det_ex_unq_dn_nb)
          endif

          det_ex2(ex_ij, det_i) = det_ex2_up (ex_ij, det_i) * detd (det_unq_dn_i) &
                                + det_ex_up (ex_i, det_i)   * det_ex_dn (ex_j, det_i)   &
                                + det_ex_up (ex_j, det_i)   * det_ex_dn (ex_i, det_i)   &
                                + detu(det_unq_up_i)        * det_ex2_dn (ex_ij, det_i)
        enddo ! det_in_csf_i
      enddo ! csf_i
    enddo ! ex_j
  enddo ! ex_i

  !write(6,*) '>det_ex2',det_ex2

  end subroutine det_ex2_bld

!===========================================================================
  subroutine  d2psi_orb_bld
!---------------------------------------------------------------------------
! Description : 
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer ex_i,ex_j,ex_rev_i,ex_rev_j
  integer ex_ij,ex_rev_i_j,ex_i_rev_j,ex_rev_i_rev_j
  integer ex_ji,ex_rev_j_i,ex_j_rev_i,ex_rev_j_rev_i
  integer csf_i, det_in_csf_i, det_i, dorb_i,dorb_j,dorb_ij
  real(dp) detex2

! header
  if (header_exe) then

   call object_create ('psid_ex2')
   call object_create ('d2psi_orb')
   call object_create ('d2csf_orb')

   call object_needed ('deriv_orb_pairs_nb')
   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_ind_rev')
   call object_needed ('single_ex_nb')
   call object_needed ('ncsf')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('det_ex2')
   call object_needed ('cdet_in_csf')
   call object_needed ('csf_coef')
   call object_needed ('psi_det')

   return

  endif

! begin

! allocations
  call object_alloc ('psid_ex2', psid_ex2, deriv_orb_pairs_nb)
  call object_alloc ('d2psi_orb', d2psi_orb, deriv_orb_pairs_nb)
  call object_alloc ('d2csf_orb', d2csf_orb, ncsf, deriv_orb_pairs_nb)

  d2csf_orb =  0.d0
  psid_ex2 =  0.d0
  dorb_ij=0
  do dorb_i = 1, param_orb_nb
    do dorb_j = 1, dorb_i
      dorb_ij=dorb_ij+1

      ex_i = ex_orb_ind (dorb_i)
      ex_rev_i = ex_orb_ind_rev (dorb_i)
      ex_j = ex_orb_ind (dorb_j)
      ex_rev_j = ex_orb_ind_rev (dorb_j)

      ex_ij         =(ex_i-1)    *single_ex_nb+ex_j
      ex_rev_i_j    =(ex_rev_i-1)*single_ex_nb+ex_j
      ex_i_rev_j    =(ex_i-1)    *single_ex_nb+ex_rev_j
      ex_rev_i_rev_j=(ex_rev_i-1)*single_ex_nb+ex_rev_j
      ex_ji         =(ex_j-1)    *single_ex_nb+ex_i
      ex_rev_j_i    =(ex_rev_j-1)*single_ex_nb+ex_i
      ex_j_rev_i    =(ex_j-1)    *single_ex_nb+ex_rev_i
      ex_rev_j_rev_i=(ex_rev_j-1)*single_ex_nb+ex_rev_i

      do csf_i = 1, ncsf
        do det_in_csf_i = 1, ndet_in_csf (csf_i)
          det_i = iwdet_in_csf (det_in_csf_i, csf_i)
          detex2 = det_ex2 (ex_ij, det_i)+det_ex2 (ex_ji, det_i)
          if (.not. l_casscf .and. ex_rev_i /= 0) then
            detex2 =  detex2-det_ex2 (ex_rev_i_j, det_i)-det_ex2 (ex_j_rev_i, det_i)
          endif
          if (.not. l_casscf .and. ex_rev_j /= 0) then
            detex2 =  detex2-det_ex2 (ex_i_rev_j, det_i)-det_ex2 (ex_rev_j_i, det_i)
          endif
          if (.not. l_casscf .and. ex_rev_i /= 0 .and. ex_rev_j /= 0) then
            detex2 =  detex2+det_ex2 (ex_rev_i_rev_j, det_i)+det_ex2 (ex_rev_j_rev_i, det_i)
          endif
          d2csf_orb(csf_i,dorb_ij)=d2csf_orb(csf_i,dorb_ij) + cdet_in_csf (det_in_csf_i, csf_i) *  0.5d0*detex2
        enddo ! det_in_csf_i

        psid_ex2 (dorb_ij) = psid_ex2 (dorb_ij) + csf_coef (csf_i, 1) * d2csf_orb(csf_i,dorb_ij)
        d2csf_orb(csf_i,dorb_ij)=d2csf_orb(csf_i,dorb_ij) /  psi_det
      enddo ! csf_i
      d2psi_orb(dorb_ij) = psid_ex2(dorb_ij) / psi_det

    !write(6,*) '>d2psi_orb',d2psi_orb(dorb_ij)
    enddo ! dorb_j
  enddo ! dorb_i

  end subroutine  d2psi_orb_bld

end module deriv_orb_mod
