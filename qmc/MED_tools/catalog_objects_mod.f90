module catalog_objects_mod

  use basic_tools_mod
  use constants_mod
  use variables_mod
  use objects_mod
  use indexes_mod

  contains

! ==============================================================================
  subroutine catalog_objects
! ------------------------------------------------------------------------------
! Description   : pre-catalog some objects often modified to speed up
!
! Created       : J. Toulouse, 05 Feb 2006
! ------------------------------------------------------------------------------
  implicit none

! begin
  call catalog_one_object ('xold', xold_index)
  call catalog_one_object ('xoldw', xoldw_index)
  call catalog_one_object ('denergy', denergy_index)
!  call catalog_one_object ('det_ex_unq_up', det_ex_unq_up_index)
!  call catalog_one_object ('det_ex_unq_dn', det_ex_unq_dn_index)
  call catalog_one_object ('gvalue', gvalue_index)
  call catalog_one_object ('g', g_index)
  call catalog_one_object ('d2d2a', d2d2a_index)
  call catalog_one_object ('d1d2a', d1d2a_index)
  call catalog_one_object ('d1d2b', d1d2b_index)
  call catalog_one_object ('d2d2b', d2d2b_index)
  call catalog_one_object ('vj', vj_index)
  call catalog_one_object ('sum_lap_lnj', sum_lap_lnj_index)
  call catalog_one_object ('orb', orb_index)
  call catalog_one_object ('dorb', dorb_index)
  call catalog_one_object ('ddorb', ddorb_index)
  call catalog_one_object ('deti_det', deti_det_index)
  call catalog_one_object ('psi_det', psi_det_index)
  call catalog_one_object ('psi_jas', psi_jas_index)
  call catalog_one_object ('psijo', psijo_index)
  call catalog_one_object ('eloc', eloc_index)
  call catalog_one_object ('slmui', slmui_index)
  call catalog_one_object ('slmdi', slmdi_index)
  call catalog_one_object ('detu', detu_index)
  call catalog_one_object ('detd', detd_index)
  call catalog_one_object ('electron', electron_index)
  call catalog_one_object ('vpsp_ex', vpsp_ex_index)
  call catalog_one_object ('dvpsp_exp', dvpsp_exp_index) !fp
  call catalog_one_object ('orbe', orbe_index)
  call catalog_one_object ('detn', detn_index)
  call catalog_one_object ('slmin', slmin_index) !fp
  call catalog_one_object ('eloc_vmc', eloc_vmc_index)
  call catalog_one_object ('eloc_dmc', eloc_dmc_index)
  call catalog_one_object ('eloc_bav', eloc_bav_index)
  call catalog_one_object ('eloc_av', eloc_av_index)
  call catalog_one_object ('eloc_av_err', eloc_av_err_index)
  call catalog_one_object ('eloc_pot_nloc', eloc_pot_nloc_index)
  call catalog_one_object ('vold', vold_index)
  call catalog_one_object ('voldw', voldw_index)
  call catalog_one_object ('div_vo', div_vo_index)
  call catalog_one_object ('div_vow', div_vow_index)
  call catalog_one_object ('dpsi_csf', dpsi_csf_index)
  call catalog_one_object ('dpsi_jas', dpsi_jas_index)
  call catalog_one_object ('dpsi_orb', dpsi_orb_index)
  call catalog_one_object ('dpsi_csf_av', dpsi_csf_av_index)
  call catalog_one_object ('dpsi_jas_av', dpsi_jas_av_index)
  call catalog_one_object ('dpsi_orb_av', dpsi_orb_av_index)
  call catalog_one_object ('deloc_csf', deloc_csf_index)
  call catalog_one_object ('deloc_jas', deloc_jas_index)
  call catalog_one_object ('d2eloc_jas', d2eloc_jas_index)
  call catalog_one_object ('deloc_orb', deloc_orb_index)
  call catalog_one_object ('deloc_av', deloc_av_index)
  call catalog_one_object ('eloc_pot', eloc_pot_index)
  call catalog_one_object ('eloc_pot_loc', eloc_pot_loc_index)
  call catalog_one_object ('eloc_pot_nloc_ex', eloc_pot_nloc_ex_index)
  call catalog_one_object ('eloc_pot_nloc_exp', eloc_pot_nloc_exp_index) !fp
  call catalog_one_object ('pe_ee', pe_ee_index)
  call catalog_one_object ('pe_en', pe_en_index)
  call catalog_one_object ('jas_pairs_nb', jas_pairs_nb_index)
  call catalog_one_object ('d2psi_jas', d2psi_jas_index)
  call catalog_one_object ('opt_nwt_nb', opt_nwt_nb_index)
  call catalog_one_object ('opt_lin_nb', opt_lin_nb_index)
  call catalog_one_object ('opt_ptb_nb', opt_ptb_nb_index)
  call catalog_one_object ('grad_nwt', grad_nwt_index)
  call catalog_one_object ('grad', grad_index)
  call catalog_one_object ('grad_ptb', grad_ptb_index)
  call catalog_one_object ('det1_det', det1_det_index)
  call catalog_one_object ('intra_sp_histo_av', intra_sp_histo_av_index)
  call catalog_one_object ('intra_sp_histo_av_err', intra_sp_histo_av_err_index)
  call catalog_one_object ('intra_sp_zv1_av', intra_sp_zv1_av_index)
  call catalog_one_object ('intra_sp_zv1_av_err', intra_sp_zv1_av_err_index)
  call catalog_one_object ('intra_sp_zv2_av', intra_sp_zv2_av_index)
  call catalog_one_object ('intra_sp_zv2_av_err', intra_sp_zv2_av_err_index)
  call catalog_one_object ('intra_sp_zv3_av', intra_sp_zv3_av_index)
  call catalog_one_object ('intra_sp_zv3_av_err', intra_sp_zv3_av_err_index)
  call catalog_one_object ('intra_sp_zv4_av', intra_sp_zv4_av_index)
  call catalog_one_object ('intra_sp_zv4_av_err', intra_sp_zv4_av_err_index)
  call catalog_one_object ('intra_sp_zv5_av', intra_sp_zv5_av_index)
  call catalog_one_object ('intra_sp_zv5_av_err', intra_sp_zv5_av_err_index)
  call catalog_one_object ('intra_sp_zvzb1_av', intra_sp_zvzb1_av_index)
  call catalog_one_object ('intra_sp_zvzb1_av_err', intra_sp_zvzb1_av_err_index)
  call catalog_one_object ('intra_sp_zvzb3_av', intra_sp_zvzb3_av_index)
  call catalog_one_object ('intra_sp_zvzb3_av_err', intra_sp_zvzb3_av_err_index)
  call catalog_one_object ('intra_sp_zvzb4_av', intra_sp_zvzb4_av_index)
  call catalog_one_object ('intra_sp_zvzb4_av_err', intra_sp_zvzb4_av_err_index)
  call catalog_one_object ('intra_sp_zv1zb3_av', intra_sp_zv1zb3_av_index)
  call catalog_one_object ('intra_sp_zv1zb3_av_err', intra_sp_zv1zb3_av_err_index)
  call catalog_one_object ('intra_sp_zvzb5_av', intra_sp_zvzb5_av_index)
  call catalog_one_object ('intra_sp_zvzb5_av_err', intra_sp_zvzb5_av_err_index)
  call catalog_one_object ('hess_sor', hess_sor_index)
  call catalog_one_object ('hess_lzr', hess_lzr_index)
  call catalog_one_object ('hess_uf', hess_uf_index)
  call catalog_one_object ('hess_tu', hess_tu_index)
  call catalog_one_object ('hess_lin', hess_lin_index)
  call catalog_one_object ('nparmcsf', nparmcsf_index)
  call catalog_one_object ('nparmj', nparmj_index)
  call catalog_one_object ('param_orb_nb', param_orb_nb_index)
  call catalog_one_object ('delta_jas_nwt', delta_jas_nwt_index)
  call catalog_one_object ('delta_jas_lin', delta_jas_lin_index)
  call catalog_one_object ('delta_jas_ptb', delta_jas_ptb_index)
  call catalog_one_object ('delta_csf_nwt', delta_csf_nwt_index)
  call catalog_one_object ('delta_csf_lin', delta_csf_lin_index)
  call catalog_one_object ('delta_csf_ptb', delta_csf_ptb_index)
  call catalog_one_object ('delta_coef_ex_nwt', delta_coef_ex_nwt_index)
  call catalog_one_object ('delta_coef_ex_lin', delta_coef_ex_lin_index)
  call catalog_one_object ('delta_coef_ex_ptb', delta_coef_ex_ptb_index)
  call catalog_one_object ('wt', wt_index)
  call catalog_one_object ('fprod', fprod_index)
  call catalog_one_object ('wgcum', wgcum_index)
  call catalog_one_object ('eold', eold_index)
  call catalog_one_object ('eoldw', eoldw_index)
  call catalog_one_object ('dist_ee_min', dist_ee_min_index)
  call catalog_one_object ('dist_ee_max', dist_ee_max_index)
  call catalog_one_object ('nwalk', nwalk_index)
  call catalog_one_object ('phin', phin_index)
  call catalog_one_object ('dphin', dphin_index)
  call catalog_one_object ('d2phin', d2phin_index)
  call catalog_one_object ('r_en', r_en_index)
  call catalog_one_object ('rvec_en', rvec_en_index)
  call catalog_one_object ('coef', coef_index)
  call catalog_one_object ('coef_orb_on_norm_basis', coef_orb_on_norm_basis_index)
  call catalog_one_object ('coef_orb_on_ortho_basis', coef_orb_on_ortho_basis_index)
  call catalog_one_object ('add_diag_mult_exp', add_diag_mult_exp_index)
  call catalog_one_object ('dpsi_av', dpsi_av_index)

!! pjas WAS
  call catalog_one_object ('psid_pjas', psid_pjas_index)
  call catalog_one_object ('vd', vd_index)

  call catalog_one_object ('cos_star_en', cos_star_en_index)
  call catalog_one_object ('grad_cos_star_en', grad_cos_star_en_index)

  call catalog_one_object ('sin_star_en', sin_star_en_index)
  call catalog_one_object ('grad_sin_star_en', grad_sin_star_en_index)


  call catalog_one_object ('star_en', star_en_index)
  call catalog_one_object ('grad_star_en', grad_star_en_index)

  call catalog_one_object ('cos_star_ee', cos_star_ee_index)
  call catalog_one_object ('grad_cos_star_ee', grad_cos_star_ee_index)

  call catalog_one_object ('sin_star_ee', sin_star_ee_index)
  call catalog_one_object ('grad_sin_star_ee', grad_sin_star_ee_index)

  call catalog_one_object ('star_ee', star_ee_index)
  call catalog_one_object ('grad_star_ee', grad_star_ee_index)

  call catalog_one_object ('dvpsp_pjas', dvpsp_pjas_index)

  call catalog_one_object ('dpsi_pjasee', dpsi_pjasee_index)
  call catalog_one_object ('grad_dpsi_pjasee',grad_dpsi_pjasee_index)
  call catalog_one_object ('lap_dpsi_pjasee',lap_dpsi_pjasee_index)

  call catalog_one_object ('dpsi_pjasen', dpsi_pjasen_index)
  call catalog_one_object ('grad_dpsi_pjasen',grad_dpsi_pjasen_index)
  call catalog_one_object ('lap_dpsi_pjasen',lap_dpsi_pjasen_index)
!!!
!! solid
!  call catalog_one_object ('ngvec_orb', ngvec_orb_index)
!!!

  call catalog_one_object ('current_walker', current_walker_index)
  call catalog_one_object ('current_walker_weight', current_walker_weight_index)
  call catalog_one_object ('total_iterations_block_nb', total_iterations_block_nb_index)
  call catalog_one_object ('total_iterations_nb', total_iterations_nb_index)
  call catalog_one_object ('walker_weights_sum_block', walker_weights_sum_block_index)
  call catalog_one_object ('walker_weights_sum', walker_weights_sum_index)
  call catalog_one_object ('walker_weights_sq_sum_block', walker_weights_sq_sum_block_index)
  call catalog_one_object ('walker_weights_sq_sum', walker_weights_sq_sum_index)

 end subroutine catalog_objects

! ==============================================================================
  subroutine catalog_one_object (object_name, object_ind)
! ------------------------------------------------------------------------------
! Description   : pre-catalog a object and store its index
!
! Created       : J. Toulouse, 05 Feb 2006
! ------------------------------------------------------------------------------
  implicit none

! i/o
  character(len=*), intent(in) :: object_name
  integer, intent(inout) :: object_ind

! local
  character(len=max_string_len_rout) :: lhere = 'catalog_one_object'

! begin

! make sure that object_ind = 0
  if (object_ind /= 0) then
   write(6,*) trim(lhere), ': object_name=',trim(object_name)
   write(6,*) trim(lhere), ': calling with object_ind=',object_ind, '/=0'
   call die (lhere)
  endif

! index of object
  object_ind = object_index (object_name)

! if object not found, catalogue it
  if (object_ind == 0) then
     call object_add (object_name)
     object_ind = objects_nb
  endif

 end subroutine catalog_one_object


end module catalog_objects_mod
