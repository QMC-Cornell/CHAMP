module catalog_routines_mod

  use all_tools_mod
  use nodes_mod
  use routines_mod
  use basis_mod
  use orbitals_mod
  use jastrow_mod
  use optimization_mod
  use opt_nwt_mod
  use opt_lin_mod
  use opt_ptb_mod
  use opt_ovlp_fn_mod
  use montecarlo_mod
  use determinants_mod
  use eloc_mod
  use grid_mod
  use deriv_mod
  use deriv_jas_mod
  use density_mod
  use intracule_mod
  use extracule_mod
  use periodic_jastrow_mod  !WAS
  use nuclei_mod
  use forces_mod
  use dipole_moment_mod
  use backflow_mod, only: xi_bld, d_xi_j_beta_d_r_i_alpha_bld, d_xi_j_beta_d_r_i_alpha_num_bld, lap_i_xi_j_beta_bld, lap_i_xi_j_beta_num_bld, &
       & xi_ee_bld, eta_bld, asymp_eta_bld, c_param_eta_bld, d_eta_d_r_i_alpha_second_bld, d_eta_d_r_j_alpha_first_bld, d_xi_ee_j_beta_d_r_i_alpha_bld, d_xi_ee_j_beta_d_r_i_alpha_num_bld, d_eta_d_r_ji_bld, lap_i_xi_ee_j_beta_bld,  lap_i_xi_ee_j_beta_num_bld, lap_i_eta_second_bld, lap_j_eta_first_bld, d2_eta_d_r_ji_2_bld, &
       & xi_en_bld, mu_bld, asymp_mu_pseudo_bld, d_param_mu_bld, d_mu_d_r_j_inuc_bld, d_mu_d_r_j_alpha_bld, d_xi_en_j_beta_d_r_i_alpha_bld, d_xi_en_j_beta_d_r_i_alpha_num_bld, lap_i_xi_en_j_beta_bld, lap_i_xi_en_j_beta_num_bld, lap_j_mu_bld, d2_mu_d_r_j_inuc_2_bld, &
       & xi_een_phi_bld, phi_bld, a_param_phi_bld, a_param_phi_parallel_cond_bld, a_param_phi_anti_parallel_cond_bld, dep_a_param_phi_anti_parallel_bld, dep_a_param_phi_parallel_bld, nb_a_param_phi_bld, d_xi_een_phi_j_beta_d_r_i_alpha_bld, d_xi_een_phi_j_beta_d_r_i_alpha_num_bld, d_phi_d_r_j_alpha_first_bld, d_phi_d_r_i_alpha_second_bld, d_phi_d_r_ji_bld, d_phi_d_r_j_inuc_bld, lap_i_xi_een_phi_j_beta_bld, lap_i_xi_een_phi_j_beta_num_bld, lap_j_phi_first_bld, lap_i_phi_second_bld, d2_phi_d_r_ji_2_bld, d2_phi_d_r_j_inuc_d_r_ji_bld, d2_phi_d_r_j_inuc_2_bld, &
       & xi_een_theta_bld, theta_bld, b_param_theta_bld, nb_b_param_theta_bld, dep_b_param_theta_bld, b_param_theta_cond_bld, d_xi_een_theta_j_beta_d_r_i_alpha_bld, d_xi_een_theta_j_beta_d_r_i_alpha_num_bld, d_theta_d_r_j_alpha_first_bld, d_theta_d_r_i_alpha_second_bld, d_theta_d_r_ji_bld, d_theta_d_r_j_inuc_bld, lap_i_xi_een_theta_j_beta_bld, lap_i_xi_een_theta_j_beta_num_bld, lap_j_theta_first_bld, lap_i_theta_second_bld, d2_theta_d_r_ji_2_bld, d2_theta_d_r_j_inuc_d_r_ji_bld, d2_theta_d_r_j_inuc_2_bld, &
       & smooth_cutoff_g_bld, d_smooth_cutoff_g_d_r_bld, d2_smooth_cutoff_g_d_r_2_bld, &
       & d_scaled_dist_en_wlk_d_r_bld, d2_scaled_dist_en_wlk_d_r_2_bld, d_scaled_dist_ee_wlk_d_r_bld, d2_scaled_dist_ee_wlk_d_r_2_bld, d_scaled_dist_een_ee_wlk_d_r_bld, d_scaled_dist_een_en_wlk_d_r_bld, d2_scaled_dist_een_ee_wlk_d_r_2_bld, d2_scaled_dist_een_en_wlk_d_r_2_bld, scaled_dist_een_ee_wlk_bld, scaled_dist_een_en_wlk_bld, scaled_dist_ee_wlk_bld, scaled_dist_en_wlk_bld, asymp_scaled_dist_two_body_bld!fp
  use backflow_basis_mod, only: orbitals_bf_bld, orbital_coeff_bf_bld, d_orbitals_d_x_beta_bf_bld, &
       & normalized_bas_funcs_bf_bld, bas_funcs_bf_bld, bas_funcs_normalization_bf_bld, d_normalized_bas_funcs_d_x_beta_bf_bld, d_bas_funcs_d_x_beta_bf_bld, d2_bas_funcs_d_x_beta_d_x_gamma_bf_bld, d2_bas_funcs_d_x_beta_d_x_gamma_bf_num_bld,&
       & slat_bas_funcs_bf_bld, slat_bas_funcs_normalization_bf_bld, slat_bas_funcs_radial_bf_bld, d_slat_bas_funcs_d_x_beta_bf_bld, d_slat_bas_funcs_radial_d_x_beta_bf_bld, d_slat_bas_funcs_radial_d_x_bf_bld, d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf_bld, d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf_bld, d2_slat_bas_funcs_radial_d_x2_bf_bld,&
       & bas_funcs_3d_theta_bf_bld, bas_funcs_3d_phi_bf_bld, d_bas_funcs_3d_theta_d_x_beta_bf_bld, d_bas_funcs_3d_phi_d_x_beta_bf_bld, d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf_bld, d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf_bld

  contains

!===============================================================================
  subroutine build_tree
! ------------------------------------------------------------------------------
! Description   : build production tree of nodes
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none


! begin
  call catalog_all_nodes_and_routines
  call execute_node_headers

 end subroutine build_tree

! ==============================================================================
  subroutine catalog_all_nodes_and_routines
! ------------------------------------------------------------------------------
! Description   : catalog all catalog_nodes
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none


! begin

! basic
  call catalog_one_node ('coord_elec_bld', coord_elec_bld, coord_elec_bld_index)
  call catalog_one_node ('coord_elec_wlk_bld', coord_elec_wlk_bld, coord_elec_wlk_bld_index)
  call catalog_one_node ('eloc_wlk_bld', eloc_wlk_bld, eloc_wlk_bld_index)
  call catalog_one_node ('eloc_test_bld', eloc_test_bld)
  call catalog_one_node ('eloc_test2_bld', eloc_test2_bld, eloc_test2_bld_index)
  call catalog_one_node ('eloc_test3_bld', eloc_test3_bld, eloc_test3_bld_index)
  call catalog_one_node ('eloc_test4_bld', eloc_test4_bld, eloc_test4_bld_index)
  call catalog_one_node ('eloc_wlk_test_bld', eloc_wlk_test_bld, eloc_wlk_test_bld_index)
  call catalog_one_node ('eloc_wlk_test2_bld', eloc_wlk_test2_bld, eloc_wlk_test2_bld_index)
  call catalog_one_node ('sigma_bld', sigma_bld)
  call catalog_one_node ('eloc_bld', eloc_bld, eloc_bld_index)
  call catalog_one_node ('eloc_vmc_bld', eloc_vmc_bld)
  call catalog_one_node ('eloc_var_bld', eloc_var_bld)
  call catalog_one_node ('eloc_sq_bld', eloc_sq_bld)
  call catalog_one_node ('dist_e_wlk_bld', dist_e_wlk_bld)
  call catalog_one_node ('vec_ee_xyz_wlk_bld', vec_ee_xyz_wlk_bld)
  call catalog_one_node ('dist_ee_wlk_bld', dist_ee_wlk_bld)
  call catalog_one_node ('spin_two_elec_wlk_bld', spin_two_elec_wlk_bld)
  call catalog_one_node ('walker_weights_bld', walker_weights_bld, walker_weights_bld_index)
  call catalog_one_node ('vec_en_xyz_wlk_bld', vec_en_xyz_wlk_bld)
  call catalog_one_node ('dist_en_wlk_bld', dist_en_wlk_bld)
  call catalog_one_node ('dist_en_scaled_wlk_bld', dist_en_scaled_wlk_bld)
  call catalog_one_node ('dist_en_scaled_deriv1_wlk_bld', dist_en_scaled_deriv1_wlk_bld)
  call catalog_one_node ('dist_ee_scaled2_wlk_bld', dist_ee_scaled2_wlk_bld)
  call catalog_one_node ('dist_en_scaled2_wlk_bld', dist_en_scaled2_wlk_bld)
  call catalog_one_node ('dist_en_scaled2_deriv1_wlk_bld', dist_en_scaled2_deriv1_wlk_bld)
  call catalog_one_node ('dist_ee_scaled2_deriv1_wlk_bld', dist_ee_scaled2_deriv1_wlk_bld)
  call catalog_one_node ('grd_dist_en_bld', grd_dist_en_bld)
  call catalog_one_node ('lap_dist_en_bld', lap_dist_en_bld)
  call catalog_one_node ('elec_nb_close_atom_wlk_bld', elec_nb_close_atom_wlk_bld)
  call catalog_one_node ('elec_spin_nb_close_atom_wlk_bld', elec_spin_nb_close_atom_wlk_bld)
  call catalog_one_node ('frac_elec_nb_close_atom_input_wlk_bld', frac_elec_nb_close_atom_input_wlk_bld)
  call catalog_one_node ('frac_elec_spin_nb_close_atom_input_wlk_bld', frac_elec_spin_nb_close_atom_input_wlk_bld)
  call catalog_one_node ('elec_nb_close_atom_encount_bld', elec_nb_close_atom_encount_bld)
  call catalog_one_node ('elec_spin_nb_close_atom_encount_bld', elec_spin_nb_close_atom_encount_bld)
  call catalog_one_node ('eloc_pot_en_bld', eloc_pot_en_bld)
  call catalog_one_node ('eloc_pot_ee_bld', eloc_pot_ee_bld)
  call catalog_one_node ('eloc_pot_ee_zv_bld', eloc_pot_ee_zv_bld)
  call catalog_one_node ('eloc_kin_jas_bld', eloc_kin_jas_bld)
  call catalog_one_node ('eloc_kin_jas_pot_ee_av_bld', eloc_kin_jas_pot_ee_av_bld)
  call catalog_one_node ('dist_nn_bld', dist_nn_bld)
  call catalog_one_node ('mass_nucl_bld', mass_nucl_bld)
  call catalog_one_node ('mass_nucl_total_bld', mass_nucl_total_bld)
  call catalog_one_node ('mass_nucl_center_bld', mass_nucl_center_bld)

! basis
  call catalog_one_node ('basis_ovlp_bld', basis_ovlp_bld)
  call catalog_one_node ('dbasis_ovlp_dz_bld', dbasis_ovlp_dz_bld)
  call catalog_one_node ('dbasis_ovlp_dz_in_eig_basis_bld', dbasis_ovlp_dz_in_eig_basis_bld)
  call catalog_one_node ('basis_ovlp_eig_bld', basis_ovlp_eig_bld)
  call catalog_one_node ('basis_ovlp_12_bld', basis_ovlp_12_bld)
  call catalog_one_node ('basis_ovlp_m12_bld', basis_ovlp_m12_bld)
  call catalog_one_node ('dbasis_ovlp_m12_dz_bld', dbasis_ovlp_m12_dz_bld)
  call catalog_one_node ('norm_basis_bld', norm_basis_bld)
  call catalog_one_node ('phin_norm_bld', phin_norm_bld)
  call catalog_one_node ('phin_ortho_bld', phin_ortho_bld)
  call catalog_one_node ('dphin_norm_dz_bld', dphin_norm_dz_bld)
  call catalog_one_node ('dphin_ortho_dz_bld', dphin_ortho_dz_bld)
  call catalog_one_node ('grd_dphin_norm_dz_bld', grd_dphin_norm_dz_bld)
  call catalog_one_node ('grd_dphin_ortho_dz_bld', grd_dphin_ortho_dz_bld)
  call catalog_one_node ('lap_dphin_norm_dz_bld', lap_dphin_norm_dz_bld)
  call catalog_one_node ('lap_dphin_ortho_dz_bld', lap_dphin_ortho_dz_bld)
  call catalog_one_node ('exp_opt_lab_read_bld', exp_opt_lab_read_bld)
  call catalog_one_node ('exp_opt_lab_bld', exp_opt_lab_bld)

! orbitals
  call catalog_one_node ('orb_occupations_bld', orb_occupations_bld)
  call catalog_one_node ('orb_ovlp_bld', orb_ovlp_bld)
  call catalog_one_node ('orb_occ_ovlp_bld', orb_occ_ovlp_bld)
  call catalog_one_node ('orb_occ_ovlp_inv_bld', orb_occ_ovlp_inv_bld)
  call catalog_one_node ('orb_cls_ovlp_inv_bld', orb_cls_ovlp_inv_bld)
  call catalog_one_node ('orb_cls_ovlp_bld', orb_cls_ovlp_bld)
  call catalog_one_node ('orb_act_ovlp_bld', orb_act_ovlp_bld)
  call catalog_one_node ('orb_vir_ovlp_bld', orb_vir_ovlp_bld)
  call catalog_one_node ('orb_cls_ovlp_eig_bld', orb_cls_ovlp_eig_bld)
  call catalog_one_node ('orb_act_ovlp_eig_bld', orb_act_ovlp_eig_bld)
  call catalog_one_node ('orb_vir_ovlp_eig_bld', orb_vir_ovlp_eig_bld)
  call catalog_one_node ('orb_cls_ovlp_m12_bld', orb_cls_ovlp_m12_bld)
  call catalog_one_node ('orb_act_ovlp_m12_bld', orb_act_ovlp_m12_bld)
  call catalog_one_node ('orb_vir_ovlp_m12_bld', orb_vir_ovlp_m12_bld)
  call catalog_one_node ('orb_sym_lab_default_bld', orb_sym_lab_default_bld)
  call catalog_one_node ('orb_opt_lab_bld', orb_opt_lab_bld)
  call catalog_one_node ('orb_ex_forbidden_bld', orb_ex_forbidden_bld)
  call catalog_one_node ('orb_optimized_bld', orb_optimized_bld)
  call catalog_one_node ('lo_bld', lo_bld)
  call catalog_one_node ('orb_on_x_bld', orb_on_x_bld)
  call catalog_one_node ('grid_on_x_bld', grid_on_x_bld)

! determinants
  call catalog_one_node ('det_unq_orb_lab_srt_bld', det_unq_orb_lab_srt_bld)
!  call catalog_one_node ('det_orb_lab_srt_bld', det_orb_lab_srt_bld)
!  call catalog_one_node ('det_unq_orb_lab_srt_bld', det_unq_orb_lab_srt_bld)
  call catalog_one_node ('orb_occ_in_adet_unq_bld', orb_occ_in_adet_unq_bld)

! csf
  call catalog_one_node ('wfdet_ovlp_bld', wfdet_ovlp_bld)
  call catalog_one_node ('det_unq_in_csf_bld', det_unq_in_csf_bld)
  call catalog_one_node ('csfs_ovlp_bld', csfs_ovlp_bld)
  call catalog_one_node ('csfs_wfdet_ovlp_bld', csfs_wfdet_ovlp_bld)
  call catalog_one_node ('dens_mat_wfdet_bld', dens_mat_wfdet_bld)
  call catalog_one_node ('vb_weights_hiberty_bld', vb_weights_hiberty_bld)
  call catalog_one_node ('first_csf_over_psid_bld', first_csf_over_psid_bld)
  call catalog_one_node ('csf_over_psid_bld', csf_over_psid_bld)
  call catalog_one_node ('product_csf_over_psid_bld', product_csf_over_psid_bld)
  call catalog_one_node ('vb_weights_chirgwin_coulson_bld', vb_weights_chirgwin_coulson_bld)
  call catalog_one_node ('vb_weights_lowdin_bav_bld', vb_weights_lowdin_bav_bld)
  call catalog_one_node ('vb_weights_lowdin_like_bav_bld', vb_weights_lowdin_like_bav_bld)
  call catalog_one_node ('iwcsf_bld', iwcsf_bld)

! general optimization method
  call catalog_one_node ('delta_param_bld', delta_param_bld)
  call catalog_one_node ('gradient_bld', gradient_bld, gradient_bld_index)
  call catalog_one_node ('gradient_energy_bld', gradient_energy_bld, gradient_energy_bld_index)
  call catalog_one_node ('gradient_energy_blk_bld', gradient_energy_blk_bld)
  call catalog_one_node ('gradient_variance_bld', gradient_variance_bld)
  call catalog_one_node ('gradient_norm_bld', gradient_norm_bld)
  call catalog_one_node ('hessian_variance_lm_bld', hessian_variance_lm_bld)
  call catalog_one_node ('hessian_variance_lmcov_bld', hessian_variance_lmcov_bld)
  call catalog_one_node ('hessian_variance_lin_bld', hessian_variance_lin_bld)
  call catalog_one_node ('hessian_variance_bld', hessian_variance_bld, hessian_variance_bld_index)

! newton optimization method
  call catalog_one_node ('hess_nwt_energy_bld', hess_nwt_energy_bld, hess_nwt_energy_bld_index)
  call catalog_one_node ('hess_nwt_bld', hess_nwt_bld, hess_nwt_bld_index)
  call catalog_one_node ('hess_nwt_eig_bld', hess_nwt_eig_bld)
  call catalog_one_node ('hess_nwt_norm_bld', hess_nwt_norm_bld)
  call catalog_one_node ('hess_nwt_stab_bld', hess_nwt_stab_bld)
  call catalog_one_node ('hess_nwt_stab_inv_bld', hess_nwt_stab_inv_bld)
  call catalog_one_node ('delta_nwt_bld', delta_nwt_bld)
  call catalog_one_node ('sh_sor_bld', sh_sor_bld)
  call catalog_one_node ('g_sor_bld', g_sor_bld)
  call catalog_one_node ('hess_sor_bld', hess_sor_bld)
  call catalog_one_node ('hess_lzr_bld', hess_lzr_bld)
  call catalog_one_node ('hess_uf_bld', hess_uf_bld)
  call catalog_one_node ('hess_tu_bld', hess_tu_bld)
  call catalog_one_node ('hess_lin_bld', hess_lin_bld)

! linear optimization method
  call catalog_one_node ('ovlp_lin_bld', ovlp_lin_bld)
  call catalog_one_node ('ovlp_lin_renorm_bld', ovlp_lin_renorm_bld)
  call catalog_one_node ('ham_lin_energy_bld', ham_lin_energy_bld)
!  call catalog_one_node ('ham_lin_energy_av_bld', ham_lin_energy_av_bld)
  call catalog_one_node ('ham_lin_variance_bld', ham_lin_variance_bld)
  call catalog_one_node ('ham_lin_bld', ham_lin_bld)
  call catalog_one_node ('ham_lin_renorm_bld', ham_lin_renorm_bld)
  call catalog_one_node ('ham_lin_renorm_stab_bld', ham_lin_renorm_stab_bld)
  call catalog_one_node ('delta_lin_bld', delta_lin_bld)
  call catalog_one_node ('psi_lin_norm_sq_bld', psi_lin_norm_sq_bld)
  call catalog_one_node ('ham_eigval_av_bld', ham_eigval_av_bld)

! perturbative optimization method
  call catalog_one_node ('delta_ptb_bld', delta_ptb_bld)
  call catalog_one_node ('e_ptb_bld', e_ptb_bld)
  call catalog_one_node ('delta_e_ptb_bld', delta_e_ptb_bld)

! overlap FN wavefunction method
  call catalog_one_node ('dpsi_over_jas2_bld', dpsi_over_jas2_bld)
  call catalog_one_node ('dpsi2_over_jas2_bld', dpsi2_over_jas2_bld)
  call catalog_one_node ('first_csf_over_jas2_bld', first_csf_over_jas2_bld)
  call catalog_one_node ('first_csf2_over_jas2_bld', first_csf2_over_jas2_bld)
  call catalog_one_node ('one_over_jas2_bld', one_over_jas2_bld)
  call catalog_one_node ('delta_ovlp_fn_linear_bld', delta_ovlp_fn_linear_bld)
  call catalog_one_node ('delta_ovlp_fn_bld', delta_ovlp_fn_bld)
  call catalog_one_node ('ovlp_trial_fn_bld', ovlp_trial_fn_bld)
  call catalog_one_node ('ovlp_trial_fn_over_ovlp_trial_bld', ovlp_trial_fn_over_ovlp_trial_bld)
  call catalog_one_node ('gradient_ovlp_bld', gradient_ovlp_bld)
  call catalog_one_node ('gradient_ovlp_norm_bld', gradient_ovlp_norm_bld)
  call catalog_one_node ('ovlp_fn_bld', ovlp_fn_bld)
  call catalog_one_node ('delta_ovlp_fn_nonsemiortho_bld', delta_ovlp_fn_nonsemiortho_bld)

! jastrow derivatives
  call catalog_one_node ('jas_pairs_bld', jas_pairs_bld)
  call catalog_one_node ('dpsi_jas_bld', dpsi_jas_bld)
  call catalog_one_node ('dpsi_jas_sq_bld', dpsi_jas_sq_bld)
  call catalog_one_node ('dpsi_jas_eloc_bld', dpsi_jas_eloc_bld)
  call catalog_one_node ('dpsi_jas_sq_eloc_bld', dpsi_jas_sq_eloc_bld)
  call catalog_one_node ('dpsi_jas_deloc_jas_bld', dpsi_jas_deloc_jas_bld)
  call catalog_one_node ('d2psi_jas_bld', d2psi_jas_bld)
  call catalog_one_node ('deloc_jas_bld', deloc_jas_bld)
  call catalog_one_node ('d2eloc_jas_bld', d2eloc_jas_bld)
  call catalog_one_node ('e_jas_bld', e_jas_bld)
  call catalog_one_node ('delta_e_jas_bld', delta_e_jas_bld)

! csf derivatives
  call catalog_one_node ('csf_pairs_bld', csf_pairs_bld)
  call catalog_one_node ('dpsi_csf_bld', dpsi_csf_bld)
  call catalog_one_node ('dpsi_csf_dpsi_csf_bld', dpsi_csf_dpsi_csf_bld)
  call catalog_one_node ('dpsi_csf_dpsi_csf_covar_bld', dpsi_csf_dpsi_csf_covar_bld)
  call catalog_one_node ('dpsi_csf_sq_bld', dpsi_csf_sq_bld)
  call catalog_one_node ('dpsi_csf_sq_eloc_bld', dpsi_csf_sq_eloc_bld)
  call catalog_one_node ('dpsi_csf_deloc_csf_bld', dpsi_csf_deloc_csf_bld)
  call catalog_one_node ('dpsi_csf_eloc_bld', dpsi_csf_eloc_bld)
  call catalog_one_node ('deloc_csf_bld', deloc_csf_bld)
  call catalog_one_node ('e_csf_bld', e_csf_bld)
  call catalog_one_node ('delta_e_csf_bld', delta_e_csf_bld)
  call catalog_one_node ('delta_csf_rot_bld', delta_csf_rot_bld)

! general derivatives
  call catalog_one_node ('param_nb_bld', param_nb_bld, param_nb_bld_index)
  call catalog_one_node ('dpsi_bld', dpsi_bld, dpsi_bld_index)
  call catalog_one_node ('dpsi_sq_bld', dpsi_sq_bld)
  call catalog_one_node ('dpsi_var_bld', dpsi_var_bld)
  call catalog_one_node ('dpsi_sq_covar_bld', dpsi_sq_covar_bld)
  call catalog_one_node ('dpsi_dpsi_bld', dpsi_dpsi_bld)
  call catalog_one_node ('dpsi_dpsi_covar_bld', dpsi_dpsi_covar_bld)
  call catalog_one_node ('dpsi_dpsi_covar_inv_bld', dpsi_dpsi_covar_inv_bld)
  call catalog_one_node ('dpsi_dpsi_uwcovar_bld', dpsi_dpsi_uwcovar_bld)
  call catalog_one_node ('dpsi_dpsi_uwcovar_inv_bld', dpsi_dpsi_uwcovar_inv_bld)
  call catalog_one_node ('dpsi_dpsi_eloc_sq_bld', dpsi_dpsi_eloc_sq_bld)
  call catalog_one_node ('dpsi_deloc_bld', dpsi_deloc_bld)
  call catalog_one_node ('dpsi_deloc_covar_bld', dpsi_deloc_covar_bld)
  call catalog_one_node ('deloc_bld', deloc_bld, deloc_bld_index)
  call catalog_one_node ('deloc_sq_bld', deloc_sq_bld)
  call catalog_one_node ('deloc_var_bld', deloc_var_bld)
  call catalog_one_node ('deloc_av_abs_max_bld', deloc_av_abs_max_bld)
  call catalog_one_node ('dpsi_sq_eloc_bld', dpsi_sq_eloc_bld)
  call catalog_one_node ('dpsi_eloc_sq_bld',dpsi_eloc_sq_bld)
  call catalog_one_node ('dpsi_eloc_sq_covar_bld',dpsi_eloc_sq_covar_bld)
  call catalog_one_node ('deloc_eloc_bld', deloc_eloc_bld)
  call catalog_one_node ('deloc_eloc_covar_bld', deloc_eloc_covar_bld)
  call catalog_one_node ('deloc_deloc_bld', deloc_deloc_bld)
  call catalog_one_node ('deloc_deloc_covar_bld', deloc_deloc_covar_bld)
  call catalog_one_node ('deloc_deloc_blk_covar_bld', deloc_deloc_blk_covar_bld)
  call catalog_one_node ('deloc_deloc_covar_inv_bld', deloc_deloc_covar_inv_bld)
  call catalog_one_node ('deloc_deloc_blk_covar_inv_bld', deloc_deloc_blk_covar_inv_bld)
  call catalog_one_node ('dpsi_dpsi_eloc_bld', dpsi_dpsi_eloc_bld)
  call catalog_one_node ('dpsi_dpsi_c_eloc_av_bld', dpsi_dpsi_c_eloc_av_bld)
  call catalog_one_node ('dpsi_eloc_bld', dpsi_eloc_bld)
  call catalog_one_node ('dpsi_eloc_covar_bld', dpsi_eloc_covar_bld)
  call catalog_one_node ('dpsi_eloc_covar_var_bld', dpsi_eloc_covar_var_bld)
  call catalog_one_node ('dpsi_eloc_covar_deloc_av_bld', dpsi_eloc_covar_deloc_av_bld)
  call catalog_one_node ('dpsi_eloc_covar_deloc_av_var_bld', dpsi_eloc_covar_deloc_av_var_bld)
  call catalog_one_node ('dpsi_eloc_blk_covar_bld', dpsi_eloc_blk_covar_bld)
  call catalog_one_node ('d2psi_bld', d2psi_bld, d2psi_bld_index)
  call catalog_one_node ('d2psi_eloc_bld', d2psi_eloc_bld)
  call catalog_one_node ('dpsi_deloc_eloc_bld', dpsi_deloc_eloc_bld)
  call catalog_one_node ('dpsi_deloc_eloc_tc_av_bld', dpsi_deloc_eloc_tc_av_bld)
  call catalog_one_node ('dpsi_dpsi_eloc_eloc_qc_av_bld', dpsi_dpsi_eloc_eloc_qc_av_bld)

! orbital derivatives
  call catalog_one_node ('single_ex_wf_bld', single_ex_wf_bld)
  call catalog_one_node ('single_ex_wf_bld_2', single_ex_wf_bld_2)
  call catalog_one_node ('single_ex_det_bld', single_ex_det_bld)
  call catalog_one_node ('slater_mat_trans_inv_bld', slater_mat_trans_inv_bld)

!  call catalog_one_node ('single_ex_det_test_bld', single_ex_det_test_bld)
!  call catalog_one_node ('det_ij_test_bld', det_ij_test_bld)
!  call catalog_one_node ('dpsi_orb_test_bld', dpsi_orb_test_bld)

  call catalog_one_node ('det_ex_unq_bld', det_ex_unq_bld)
  call catalog_one_node ('det_ex_bld', det_ex_bld)
  call catalog_one_node ('dpsi_orb_bld', dpsi_orb_bld)

!  call catalog_one_node ('pot_efp_bld', pot_efp_bld)

!  call catalog_one_node ('pot_efp_over_delta_eps_bld', pot_efp_over_delta_eps_bld)
!  call catalog_one_node ('delta_coef_pefp_bld', delta_coef_pefp_bld)
  call catalog_one_node ('delta_eps_bld', delta_eps_bld)
!  call catalog_one_node ('grad_orb_bld', grad_orb_bld)
!  call catalog_one_node ('grad_orb_norm_bld', grad_orb_norm_bld)
!  call catalog_one_node ('hess_orb_pefp_bld', hess_orb_pefp_bld)
!  call catalog_one_node ('hess_orb_pefp_stab_bld', hess_orb_pefp_stab_bld)
!  call catalog_one_node ('hess_orb_pefp_inv_bld', hess_orb_pefp_inv_bld)
!  call catalog_one_node ('hess_orb_pefp_stab_inv_bld', hess_orb_pefp_stab_inv_bld)
!  call catalog_one_node ('hess_orb_pefp_eig_abs_bld', hess_orb_pefp_eig_abs_bld)

  call catalog_one_node ('grd_det_unq_bld', grd_det_unq_bld)
  call catalog_one_node ('grd_psid_over_psid_bld', grd_psid_over_psid_bld)
  call catalog_one_node ('grd_psi_over_psi_wlk_bld', grd_psi_over_psi_wlk_bld, grd_psi_over_psi_wlk_bld_index)
  call catalog_one_node ('div_grd_psi_over_psi_wlk_bld', div_grd_psi_over_psi_wlk_bld, div_grd_psi_over_psi_wlk_bld_index)
  call catalog_one_node ('grd_psi_over_psi_sq_wlk_bld', grd_psi_over_psi_sq_wlk_bld)
  call catalog_one_node ('grd_psi_over_psi_old_bld', grd_psi_over_psi_old_bld)
  call catalog_one_node ('lap_psi_over_psi_wlk_bld', lap_psi_over_psi_wlk_bld)
  call catalog_one_node ('lap_det_unq_bld', lap_det_unq_bld)
  call catalog_one_node ('lap_psid_over_psid_bld', lap_psid_over_psid_bld)
  call catalog_one_node ('lap_lnpsid_bld', lap_lnpsid_bld)
  call catalog_one_node ('sum_lap_lnpsid_bld', sum_lap_lnpsid_bld)
  call catalog_one_node ('sum_lap_lnpsi_bld', sum_lap_lnpsi_bld)
  call catalog_one_node ('eloc_kin_bld', eloc_kin_bld)
!  call catalog_one_node ('eloc_bld', eloc_bld)

  call catalog_one_node ('slater_mat_ex_trans_inv_bld', slater_mat_ex_trans_inv_bld)
  call catalog_one_node ('slater_mat_ex_trans_inv_2_bld', slater_mat_ex_trans_inv_2_bld)
  call catalog_one_node ('grd_det_ex_unq_bld', grd_det_ex_unq_bld)
  call catalog_one_node ('grd_det_ex_bld', grd_det_ex_bld)
  call catalog_one_node ('lap_det_ex_unq_bld', lap_det_ex_unq_bld)
  call catalog_one_node ('lap_det_ex_bld', lap_det_ex_bld)
  call catalog_one_node ('grd_psid_ex_over_psid_bld', grd_psid_ex_over_psid_bld)
  call catalog_one_node ('grd_psi_ex_over_psi_bld', grd_psi_ex_over_psi_bld)
  call catalog_one_node ('lap_psid_ex_over_psid_bld', lap_psid_ex_over_psid_bld)
  call catalog_one_node ('lap_lnpsid_ex_bld', lap_lnpsid_ex_bld)
  call catalog_one_node ('sum_lap_lnpsid_ex_bld', sum_lap_lnpsid_ex_bld)
  call catalog_one_node ('sum_lap_lnpsi_ex_bld', sum_lap_lnpsi_ex_bld)
  call catalog_one_node ('eloc_kin_ex_bld', eloc_kin_ex_bld)
  call catalog_one_node ('eloc_pot_ex_bld', eloc_pot_ex_bld, eloc_pot_ex_bld_index)
  call catalog_one_node ('eloc_ex_bld', eloc_ex_bld)
  call catalog_one_node ('psid_ex_in_x_bld', psid_ex_in_x_bld)
  call catalog_one_node ('eloc_pot_nloc_ex_bld', eloc_pot_nloc_ex_bld)
  call catalog_one_node ('deloc_orb_bld', deloc_orb_bld)
  call catalog_one_node ('delta_mat_rot_1st_order_bld', delta_mat_rot_1st_order_bld)
  call catalog_one_node ('delta_mat_rot_real_bld', delta_mat_rot_real_bld)
  call catalog_one_node ('delta_coef_pw_bld', delta_coef_pw_bld)

! exponent derivatives
  call catalog_one_node ('slater_mat_trans_inv_in_x_bld', slater_mat_trans_inv_in_x_bld)
  call catalog_one_node ('basis_fns_cent_bld', basis_fns_cent_bld)
  call catalog_one_node ('basis_fns_type_bld', basis_fns_type_bld)
  call catalog_one_node ('basis_fns_name_bld', basis_fns_name_bld)
  call catalog_one_node ('param_exp_nb_bld', param_exp_nb_bld)
  call catalog_one_node ('orbital_depends_on_opt_exp_bld', orbital_depends_on_opt_exp_bld)
  call catalog_one_node ('dnorm_basis_dz_bld', dnorm_basis_dz_bld)
  call catalog_one_node ('dphin_dz_bld', dphin_dz_bld)
  call catalog_one_node ('grd_dphin_dz_bld', grd_dphin_dz_bld)
  call catalog_one_node ('lap_dphin_dz_bld', lap_dphin_dz_bld)
  call catalog_one_node ('dorb_dexp_bld', dorb_dexp_bld, dorb_dexp_bld_index)
  call catalog_one_node ('grd_dorb_dexp_bld', grd_dorb_dexp_bld, grd_dorb_dexp_bld_index)
  call catalog_one_node ('lap_dorb_dexp_bld', lap_dorb_dexp_bld, lap_dorb_dexp_bld_index)
  call catalog_one_node ('ddet_dexp_unq_bld', ddet_dexp_unq_bld)
  call catalog_one_node ('grd_ddet_dexp_unq_bld', grd_ddet_dexp_unq_bld)
  call catalog_one_node ('lap_ddet_dexp_unq_bld', lap_ddet_dexp_unq_bld)
  call catalog_one_node ('dpsi_exp_bld', dpsi_exp_bld)
  call catalog_one_node ('dpsid_exp_in_x_bld', dpsid_exp_in_x_bld) !fp
  call catalog_one_node ('ddet_dexp_unq_in_x_bld', ddet_dexp_unq_in_x_bld) !fp
  call catalog_one_node ('dpsi_lnexp_bld', dpsi_lnexp_bld)
  call catalog_one_node ('grd_dpsid_exp_over_dpsid_exp_bld', grd_dpsid_exp_over_dpsid_exp_bld)
  call catalog_one_node ('lap_dpsid_exp_over_dpsid_exp_bld', lap_dpsid_exp_over_dpsid_exp_bld)
  call catalog_one_node ('lap_ln_dpsid_exp_bld', lap_ln_dpsid_exp_bld)
  call catalog_one_node ('sum_lap_ln_dpsid_exp_bld', sum_lap_ln_dpsid_exp_bld)
  call catalog_one_node ('grd_dpsi_exp_over_dpsi_exp_bld', grd_dpsi_exp_over_dpsi_exp_bld)
  call catalog_one_node ('sum_lap_ln_dpsi_exp_bld', sum_lap_ln_dpsi_exp_bld)
  call catalog_one_node ('eloc_kin_exp_bld', eloc_kin_exp_bld)
  call catalog_one_node ('eloc_exp_bld', eloc_exp_bld)
  call catalog_one_node ('eloc_pot_nloc_exp_bld', eloc_pot_nloc_exp_bld) !fp
  call catalog_one_node ('eloc_pot_exp_bld', eloc_pot_exp_bld, eloc_pot_exp_bld_index) !fp
  call catalog_one_node ('deloc_exp_bld', deloc_exp_bld)
  call catalog_one_node ('deloc_lnexp_bld', deloc_lnexp_bld)
  call catalog_one_node ('slater_mat_exp_trans_inv_bld', slater_mat_exp_trans_inv_bld)

! geometry derivatives
  call catalog_one_node ('param_geo_nb_bld', param_geo_nb_bld)
  call catalog_one_node ('dpsi_geo_bld', dpsi_geo_bld)
  call catalog_one_node ('deloc_geo_bld', deloc_geo_bld)

! jastrow
  call catalog_one_node ('fen_bld', fen_bld)
  call catalog_one_node ('feen_bld', feen_bld)
  call catalog_one_node ('dfen_drn_bld', dfen_drn_bld)
  call catalog_one_node ('dfeen_drn_bld', dfeen_drn_bld)

! dipole moment
  call catalog_one_node ('dipole_moment_origin_bld', dipole_moment_origin_bld)
  call catalog_one_node ('dipole_moment_nucl_bld', dipole_moment_nucl_bld)
  call catalog_one_node ('dipole_moment_bld', dipole_moment_bld)
  call catalog_one_node ('dipole_moment_dpsi_bld', dipole_moment_dpsi_bld)
  call catalog_one_node ('dipole_moment_dpsi_dpsi_bld', dipole_moment_dpsi_dpsi_bld)
  call catalog_one_node ('dipole_moment_deloc_bld', dipole_moment_deloc_bld)
  call catalog_one_node ('dipole_moment_deloc_covar_bld', dipole_moment_deloc_covar_bld)
  call catalog_one_node ('dipole_moment_deloc_blk_covar_bld', dipole_moment_deloc_blk_covar_bld)
  call catalog_one_node ('dipole_moment_zv_coef_bld', dipole_moment_zv_coef_bld)
  call catalog_one_node ('dipole_moment_zv_coef_blk_bld', dipole_moment_zv_coef_blk_bld)
  call catalog_one_node ('dipole_moment_zv_av_bld', dipole_moment_zv_av_bld)
  call catalog_one_node ('dipole_moment_zv_bav_bld', dipole_moment_zv_bav_bld)
  call catalog_one_node ('dipole_moment_zv_bav_sq_bld', dipole_moment_zv_bav_sq_bld)
  call catalog_one_node ('dipole_moment_zv_av_var_bld', dipole_moment_zv_av_var_bld)
  call catalog_one_node ('dipole_moment_delta_zb_av_bld', dipole_moment_delta_zb_av_bld)
  call catalog_one_node ('dipole_moment_delta_zb_av_var_bld', dipole_moment_delta_zb_av_var_bld)
  call catalog_one_node ('dipole_moment_zvzb_av_bld', dipole_moment_zvzb_av_bld)
  call catalog_one_node ('dipole_moment_zvzb_bav_bld', dipole_moment_zvzb_bav_bld)
  call catalog_one_node ('dipole_moment_zvzb_av_var_bld', dipole_moment_zvzb_av_var_bld)
  call catalog_one_node ('dipole_moment_lin_av_bld', dipole_moment_lin_av_bld)

! grids
  call catalog_one_node ('grid_r_bld', grid_r_bld)
  call catalog_one_node ('grid_xy_z_bld', grid_xy_z_bld)
  call catalog_one_node ('grid_xyz_bld', grid_xyz_bld)

! density
  call catalog_one_node ('dens_bld', dens_bld)
  call catalog_one_node ('dens_histo_bld', dens_histo_bld)
  call catalog_one_node ('dens_zv1_bld', dens_zv1_bld)
  call catalog_one_node ('dens_zv5_bld', dens_zv5_bld)
  call catalog_one_node ('dens_q1_bld', dens_q1_bld)
  call catalog_one_node ('dens_q1_eloc_bld', dens_q1_eloc_bld)
  call catalog_one_node ('dens_zb1_av_bld', dens_zb1_av_bld)
  call catalog_one_node ('dens_zvzb1_av_bld', dens_zvzb1_av_bld)
  call catalog_one_node ('dens_q5_bld', dens_q5_bld)
  call catalog_one_node ('dens_q5_eloc_bld', dens_q5_eloc_bld)
  call catalog_one_node ('dens_zb5_av_bld', dens_zb5_av_bld)
  call catalog_one_node ('dens_zvzb5_av_bld', dens_zvzb5_av_bld)
  call catalog_one_node ('dens_xy_z_bld', dens_xy_z_bld)
  call catalog_one_node ('dens_xy_z_histo_bld', dens_xy_z_histo_bld)
  call catalog_one_node ('dens_3d_bld', dens_3d_bld)
  call catalog_one_node ('dens_3d_histo_bld', dens_3d_histo_bld)
  call catalog_one_node ('dens_3d_zv1_bld', dens_3d_zv1_bld)
  call catalog_one_node ('dens_3d_zv2_bld', dens_3d_zv2_bld)
  call catalog_one_node ('dens_3d_zv2_linear_av_bld', dens_3d_zv2_linear_av_bld)
  call catalog_one_node ('dens_3d_zv2_linear_coef_bld', dens_3d_zv2_linear_coef_bld)
  call catalog_one_node ('dens_3d_zv2_deloc_covar_bld', dens_3d_zv2_deloc_covar_bld)
  call catalog_one_node ('dens_3d_zv2_deloc_bld', dens_3d_zv2_deloc_bld)
  call catalog_one_node ('dens_3d_zv2_linear_av_var_bld', dens_3d_zv2_linear_av_var_bld)
  call catalog_one_node ('dens_fourier_bld', dens_fourier_bld)
  call catalog_one_node ('dens_fourier_real_bld', dens_fourier_real_bld)

! intracule
  call catalog_one_node ('intra_bld', intra_bld)
  call catalog_one_node ('intra_4pir2_bld', intra_4pir2_bld)
!  call catalog_one_node ('intra_histo_bld', intra_histo_bld)
!  call catalog_one_node ('intra_zv1_bld', intra_zv1_bld)
!  call catalog_one_node ('intra_zv2_bld', intra_zv2_bld)
  call catalog_one_node ('intra_sp_bld', intra_sp_bld, intra_sp_bld_index)
  call catalog_one_node ('intra_sp_4pir2_bld', intra_sp_4pir2_bld)
  call catalog_one_node ('intra_sp_histo_bld', intra_sp_histo_bld)
  call catalog_one_node ('intra_sp_zv1_bld', intra_sp_zv1_bld)
  call catalog_one_node ('intra_sp_zv2_bld', intra_sp_zv2_bld)
  call catalog_one_node ('intra_sp_zv3_bld', intra_sp_zv3_bld)
  call catalog_one_node ('intra_sp_zv4_bld', intra_sp_zv4_bld)
  call catalog_one_node ('intra_sp_zv5_bld', intra_sp_zv5_bld)
!  call catalog_one_node ('intra_sp_q1_eloc_bld', intra_sp_q1_eloc_bld)
!  call catalog_one_node ('intra_sp_q1_bld', intra_sp_q1_bld)
!  call catalog_one_node ('intra_sp_zb1_av_bld', intra_sp_zb1_av_bld)
!  call catalog_one_node ('intra_sp_zvzb1_av_bld', intra_sp_zvzb1_av_bld)
  call catalog_one_node ('intra_sp_zvzb1_bld', intra_sp_zvzb1_bld)
  call catalog_one_node ('intra_sp_zvzb3_bld', intra_sp_zvzb3_bld)
  call catalog_one_node ('intra_sp_zvzb4_bld', intra_sp_zvzb4_bld)
  call catalog_one_node ('intra_sp_zv1zb3_bld', intra_sp_zv1zb3_bld)
!  call catalog_one_node ('intra_sp_q5_eloc_wlk_bld', intra_sp_q5_eloc_wlk_bld)
!  call catalog_one_node ('intra_sp_q5_wlk_bld', intra_sp_q5_wlk_bld)
!  call catalog_one_node ('intra_sp_zb5_av_bld', intra_sp_zb5_av_bld)
!  call catalog_one_node ('intra_sp_zvzb5_av_bld', intra_sp_zvzb5_av_bld)
  call catalog_one_node ('intra_sp_zvzb5_bld', intra_sp_zvzb5_bld)
  call catalog_one_node ('intra_3d_bld', intra_3d_bld)
  call catalog_one_node ('intra_3d_histo_bld', intra_3d_histo_bld)
  call catalog_one_node ('intra_3d_zv1_bld', intra_3d_zv1_bld)
  call catalog_one_node ('intra_3d_zv2_bld', intra_3d_zv2_bld)

! extracule
  call catalog_one_node ('grid_extra_bld', grid_extra_bld)
  call catalog_one_node ('extracule_bld', extracule_bld)
  call catalog_one_node ('extracule_4pir2_bld', extracule_4pir2_bld)
  call catalog_one_node ('extracule_zv1_bld', extracule_zv1_bld)

! forces
  call catalog_one_node ('forces_list_bld', forces_list_bld)
  call catalog_one_node ('forces_nn_bld', forces_nn_bld)
  call catalog_one_node ('forces_bare_bld', forces_bare_bld)
  call catalog_one_node ('forces_zv_bld', forces_zv_bld)
  call catalog_one_node ('forces_zvzb_av_bld', forces_zvzb_av_bld)
  call catalog_one_node ('forces_q_bld', forces_q_bld)
  call catalog_one_node ('forces_q_eloc_bld', forces_q_eloc_bld)
  call catalog_one_node ('forces_zvzb_av_var_bld', forces_zvzb_av_var_bld)
  call catalog_one_node ('forces_zv_linear_coef_bld', forces_zv_linear_coef_bld)
  call catalog_one_node ('forces_zv_deloc_bld', forces_zv_deloc_bld)
  call catalog_one_node ('forces_zv_deloc_covar_bld', forces_zv_deloc_covar_bld)
  call catalog_one_node ('forces_zv_linear_av_bld', forces_zv_linear_av_bld)
  call catalog_one_node ('forces_zv_linear_av_var_bld', forces_zv_linear_av_var_bld)
  call catalog_one_node ('forces_zv_linear_var_bld', forces_zv_linear_var_bld)
  call catalog_one_node ('forces_zv_sq_bld', forces_zv_sq_bld)
  call catalog_one_node ('forces_zv_var_bld', forces_zv_var_bld)
  call catalog_one_node ('forces_zv_pulay_av_bld', forces_zv_pulay_av_bld)
  call catalog_one_node ('forces_zv_pulay_av_var_bld', forces_zv_pulay_av_var_bld)
  call catalog_one_node ('forces_zv_deriv_bld', forces_zv_deriv_bld)
  call catalog_one_node ('forces_zv_deriv_sq_bld', forces_zv_deriv_sq_bld)
  call catalog_one_node ('forces_zv_deriv_var_bld', forces_zv_deriv_var_bld)
  call catalog_one_node ('forces_zv_deriv_pulay_av_bld', forces_zv_deriv_pulay_av_bld)
  call catalog_one_node ('forces_zv_deriv_pulay_av_var_bld', forces_zv_deriv_pulay_av_var_bld)
  call catalog_one_node ('forces_zv_deriv_deloc_bld', forces_zv_deriv_deloc_bld)
  call catalog_one_node ('forces_zv_deriv_deloc_covar_bld', forces_zv_deriv_deloc_covar_bld)
  call catalog_one_node ('forces_zv_deriv_linear_coef_bld', forces_zv_deriv_linear_coef_bld)
  call catalog_one_node ('forces_zv_deriv_linear_var_bld', forces_zv_deriv_linear_var_bld)
  call catalog_one_node ('forces_zv_deriv_linear_av_bld', forces_zv_deriv_linear_av_bld)
  call catalog_one_node ('forces_zv_deriv_linear_av_var_bld', forces_zv_deriv_linear_av_var_bld)

!backflow fp
  call catalog_one_node('xi_bld', xi_bld)
  call catalog_one_node('d_xi_j_beta_d_r_i_alpha_bld', d_xi_j_beta_d_r_i_alpha_bld)
  call catalog_one_node('d_xi_j_beta_d_r_i_alpha_num_bld', d_xi_j_beta_d_r_i_alpha_num_bld)
  call catalog_one_node('lap_i_xi_j_beta_bld', lap_i_xi_j_beta_bld)
  call catalog_one_node('lap_i_xi_j_beta_num_bld', lap_i_xi_j_beta_num_bld)

  call catalog_one_node('xi_en_bld', xi_en_bld)
  call catalog_one_node('mu_bld', mu_bld)
  call catalog_one_node('asymp_mu_pseudo_bld', asymp_mu_pseudo_bld)
  call catalog_one_node('d_mu_d_r_j_inuc_bld', d_mu_d_r_j_inuc_bld)
  call catalog_one_node('d_mu_d_r_j_alpha_bld', d_mu_d_r_j_alpha_bld)
  call catalog_one_node('d_param_mu_bld',  d_param_mu_bld)
  call catalog_one_node('d_xi_en_j_beta_d_r_i_alpha_bld', d_xi_en_j_beta_d_r_i_alpha_bld)
  call catalog_one_node('d_xi_en_j_beta_d_r_i_alpha_num_bld', d_xi_en_j_beta_d_r_i_alpha_num_bld)
  call catalog_one_node('lap_i_xi_en_j_beta_bld', lap_i_xi_en_j_beta_bld)
  call catalog_one_node('lap_i_xi_en_j_beta_num_bld', lap_i_xi_en_j_beta_num_bld)
  call catalog_one_node('lap_j_mu_bld', lap_j_mu_bld)
  call catalog_one_node('d2_mu_d_r_j_inuc_2_bld', d2_mu_d_r_j_inuc_2_bld)

  call catalog_one_node('xi_ee_bld', xi_ee_bld)
  call catalog_one_node('eta_bld', eta_bld)
  call catalog_one_node('asymp_eta_bld', asymp_eta_bld)
  call catalog_one_node('c_param_eta_bld', c_param_eta_bld)
  call catalog_one_node('d_xi_ee_j_beta_d_r_i_alpha_num_bld', d_xi_ee_j_beta_d_r_i_alpha_num_bld)
  call catalog_one_node('d_xi_ee_j_beta_d_r_i_alpha_bld', d_xi_ee_j_beta_d_r_i_alpha_bld)
  call catalog_one_node('d_eta_d_r_j_alpha_first_bld',d_eta_d_r_j_alpha_first_bld)
  call catalog_one_node('d_eta_d_r_i_alpha_second_bld',d_eta_d_r_i_alpha_second_bld)
  call catalog_one_node('d_eta_d_r_ji_bld', d_eta_d_r_ji_bld)
  call catalog_one_node('d2_eta_d_r_ji_2_bld', d2_eta_d_r_ji_2_bld)
  call catalog_one_node('lap_i_xi_ee_j_beta_bld',lap_i_xi_ee_j_beta_bld)
  call catalog_one_node('lap_i_xi_ee_j_beta_num_bld',lap_i_xi_ee_j_beta_num_bld)
  call catalog_one_node('lap_i_eta_second_bld', lap_i_eta_second_bld)
  call catalog_one_node('lap_j_eta_first_bld', lap_j_eta_first_bld)

  call catalog_one_node('xi_een_phi_bld', xi_een_phi_bld)
  call catalog_one_node('a_param_phi_bld', a_param_phi_bld)
  call catalog_one_node('a_param_phi_parallel_cond_bld', a_param_phi_parallel_cond_bld)
  call catalog_one_node('a_param_phi_anti_parallel_cond_bld', a_param_phi_anti_parallel_cond_bld)
  call catalog_one_node('dep_a_param_phi_parallel_bld', dep_a_param_phi_parallel_bld)
  call catalog_one_node('dep_a_param_phi_anti_parallel_bld', dep_a_param_phi_anti_parallel_bld)
  call catalog_one_node('nb_a_param_phi_bld', nb_a_param_phi_bld)
  call catalog_one_node('phi_bld', phi_bld)
  call catalog_one_node('d_phi_d_r_j_inuc_bld',d_phi_d_r_j_inuc_bld)
  call catalog_one_node('d_xi_een_phi_j_beta_d_r_i_alpha_bld',d_xi_een_phi_j_beta_d_r_i_alpha_bld)
  call catalog_one_node('d_xi_een_phi_j_beta_d_r_i_alpha_num_bld',d_xi_een_phi_j_beta_d_r_i_alpha_num_bld)
  call catalog_one_node('d_phi_d_r_ji_bld',d_phi_d_r_ji_bld)
  call catalog_one_node('d_phi_d_r_j_alpha_first_bld',d_phi_d_r_j_alpha_first_bld)
  call catalog_one_node('d_phi_d_r_i_alpha_second_bld',d_phi_d_r_i_alpha_second_bld)
  call catalog_one_node('lap_i_xi_een_phi_j_beta_bld', lap_i_xi_een_phi_j_beta_bld)
  call catalog_one_node('lap_i_xi_een_phi_j_beta_num_bld', lap_i_xi_een_phi_j_beta_num_bld)
  call catalog_one_node('lap_j_phi_first_bld', lap_j_phi_first_bld)
  call catalog_one_node('lap_i_phi_second_bld', lap_i_phi_second_bld)
  call catalog_one_node('d2_phi_d_r_j_inuc_d_r_ji_bld', d2_phi_d_r_j_inuc_d_r_ji_bld)
  call catalog_one_node('d2_phi_d_r_ji_2_bld', d2_phi_d_r_ji_2_bld)
  call catalog_one_node('d2_phi_d_r_j_inuc_2_bld', d2_phi_d_r_j_inuc_2_bld)

  call catalog_one_node('xi_een_theta_bld', xi_een_theta_bld)
  call catalog_one_node('b_param_theta_bld', b_param_theta_bld)
  call catalog_one_node('b_param_theta_cond_bld', b_param_theta_cond_bld)
  call catalog_one_node('nb_b_param_theta_bld', nb_b_param_theta_bld)
  call catalog_one_node('theta_bld', theta_bld)
  call catalog_one_node('dep_b_param_theta_bld', dep_b_param_theta_bld)
  call catalog_one_node('d_theta_d_r_j_inuc_bld',d_theta_d_r_j_inuc_bld)
  call catalog_one_node('d_xi_een_theta_j_beta_d_r_i_alpha_bld',d_xi_een_theta_j_beta_d_r_i_alpha_bld)
  call catalog_one_node('d_xi_een_theta_j_beta_d_r_i_alpha_num_bld',d_xi_een_theta_j_beta_d_r_i_alpha_num_bld)
  call catalog_one_node('d_theta_d_r_ji_bld',d_theta_d_r_ji_bld)
  call catalog_one_node('d_theta_d_r_j_alpha_first_bld',d_theta_d_r_j_alpha_first_bld)
  call catalog_one_node('d_theta_d_r_i_alpha_second_bld',d_theta_d_r_i_alpha_second_bld)
  call catalog_one_node('lap_i_xi_een_theta_j_beta_bld', lap_i_xi_een_theta_j_beta_bld)
  call catalog_one_node('lap_i_xi_een_theta_j_beta_num_bld', lap_i_xi_een_theta_j_beta_num_bld)
  call catalog_one_node('lap_j_theta_first_bld', lap_j_theta_first_bld)
  call catalog_one_node('lap_i_theta_second_bld', lap_i_theta_second_bld)
  call catalog_one_node('d2_theta_d_r_j_inuc_d_r_ji_bld', d2_theta_d_r_j_inuc_d_r_ji_bld)
  call catalog_one_node('d2_theta_d_r_ji_2_bld', d2_theta_d_r_ji_2_bld)
  call catalog_one_node('d2_theta_d_r_j_inuc_2_bld', d2_theta_d_r_j_inuc_2_bld)

  call catalog_one_node('scaled_dist_een_ee_wlk_bld', scaled_dist_een_ee_wlk_bld)
  call catalog_one_node('scaled_dist_een_en_wlk_bld', scaled_dist_een_en_wlk_bld)
  call catalog_one_node('scaled_dist_ee_wlk_bld', scaled_dist_ee_wlk_bld)
  call catalog_one_node('scaled_dist_en_wlk_bld', scaled_dist_en_wlk_bld)
  call catalog_one_node('asymp_scaled_dist_two_body_bld', asymp_scaled_dist_two_body_bld)
  call catalog_one_node('d_scaled_dist_en_wlk_d_r_bld', d_scaled_dist_en_wlk_d_r_bld)
  call catalog_one_node('d2_scaled_dist_en_wlk_d_r_2_bld', d2_scaled_dist_en_wlk_d_r_2_bld)
  call catalog_one_node('d_scaled_dist_ee_wlk_d_r_bld', d_scaled_dist_ee_wlk_d_r_bld)
  call catalog_one_node('d2_scaled_dist_ee_wlk_d_r_2_bld', d2_scaled_dist_ee_wlk_d_r_2_bld)
  call catalog_one_node('d_scaled_dist_een_en_wlk_d_r_bld', d_scaled_dist_een_en_wlk_d_r_bld)
  call catalog_one_node('d_scaled_dist_een_ee_wlk_d_r_bld', d_scaled_dist_een_ee_wlk_d_r_bld)
  call catalog_one_node('d2_scaled_dist_een_en_wlk_d_r_2_bld', d2_scaled_dist_een_en_wlk_d_r_2_bld)
  call catalog_one_node('d2_scaled_dist_een_ee_wlk_d_r_2_bld', d2_scaled_dist_een_ee_wlk_d_r_2_bld)

  call catalog_one_node('smooth_cutoff_g_bld',smooth_cutoff_g_bld)
  call catalog_one_node('d_smooth_cutoff_g_d_r_bld',d_smooth_cutoff_g_d_r_bld)
  call catalog_one_node('d2_smooth_cutoff_g_d_r_2_bld',d2_smooth_cutoff_g_d_r_2_bld)

  !backflow_basis
  call catalog_one_node('orbitals_bf_bld', orbitals_bf_bld)
  call catalog_one_node('orbital_coeff_bf_bld', orbital_coeff_bf_bld)
  call catalog_one_node('d_orbitals_d_x_beta_bf_bld', d_orbitals_d_x_beta_bf_bld)

  call catalog_one_node('normalized_bas_funcs_bf_bld', normalized_bas_funcs_bf_bld)
  call catalog_one_node('bas_funcs_bf_bld', bas_funcs_bf_bld)
  call catalog_one_node('bas_funcs_normalization_bf_bld', bas_funcs_normalization_bf_bld)
  call catalog_one_node('d_normalized_bas_funcs_d_x_beta_bf_bld', d_normalized_bas_funcs_d_x_beta_bf_bld)
  call catalog_one_node('d_bas_funcs_d_x_beta_bf_bld', d_bas_funcs_d_x_beta_bf_bld)
  call catalog_one_node('d2_bas_funcs_d_x_beta_d_x_gamma_bf_bld', d2_bas_funcs_d_x_beta_d_x_gamma_bf_bld)
  call catalog_one_node('d2_bas_funcs_d_x_beta_d_x_gamma_bf_num_bld', d2_bas_funcs_d_x_beta_d_x_gamma_bf_num_bld)

  call catalog_one_node('slat_bas_funcs_bf_bld', slat_bas_funcs_bf_bld)
  call catalog_one_node('slat_bas_funcs_normalization_bf_bld', slat_bas_funcs_normalization_bf_bld)
  call catalog_one_node('slat_bas_funcs_radial_bf_bld', slat_bas_funcs_radial_bf_bld)
  call catalog_one_node('d_slat_bas_funcs_d_x_beta_bf_bld', d_slat_bas_funcs_d_x_beta_bf_bld)
  call catalog_one_node('d_slat_bas_funcs_radial_d_x_beta_bf_bld', d_slat_bas_funcs_radial_d_x_beta_bf_bld)
  call catalog_one_node('d_slat_bas_funcs_radial_d_x_bf_bld', d_slat_bas_funcs_radial_d_x_bf_bld)
  call catalog_one_node('d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf_bld', d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf_bld)
  call catalog_one_node('d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf_bld', d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf_bld)
  call catalog_one_node('d2_slat_bas_funcs_radial_d_x2_bf_bld', d2_slat_bas_funcs_radial_d_x2_bf_bld)

  call catalog_one_node('bas_funcs_3d_theta_bf_bld', bas_funcs_3d_theta_bf_bld)
  call catalog_one_node('bas_funcs_3d_phi_bf_bld', bas_funcs_3d_phi_bf_bld)
  call catalog_one_node('d_bas_funcs_3d_theta_d_x_beta_bf_bld', d_bas_funcs_3d_theta_d_x_beta_bf_bld)
  call catalog_one_node('d_bas_funcs_3d_phi_d_x_beta_bf_bld', d_bas_funcs_3d_phi_d_x_beta_bf_bld)
  call catalog_one_node('d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf_bld', d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf_bld)
  call catalog_one_node('d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf_bld', d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf_bld)

! forces pulay
  call catalog_one_node ('force_to_bas_bld', force_to_bas_bld)
  call catalog_one_node ('dphin_norm_drn_bld', dphin_norm_drn_bld)
  call catalog_one_node ('dorb_drn_bld', dorb_drn_bld)
  call catalog_one_node ('ddet_drn_unq_bld', ddet_drn_unq_bld)
  call catalog_one_node ('dpsid_rn_bld', dpsid_rn_bld)
  call catalog_one_node ('dpsij_rn_bld', dpsij_rn_bld)
  call catalog_one_node ('dpsi_rn_bld', dpsi_rn_bld)
  call catalog_one_node ('dpsi_rn_num_bld', dpsi_rn_num_bld)
  call catalog_one_node ('dpsi_rn_eloc_bld', dpsi_rn_eloc_bld)
  call catalog_one_node ('deloc_rn_num_bld', deloc_rn_num_bld)
  call catalog_one_node ('forces_pulay_av_bld', forces_pulay_av_bld)
  call catalog_one_node ('forces_pulay_av_var_bld', forces_pulay_av_var_bld)

! average routines (not nodes!)
!  call catalog_one_routine ('dpsi_orb_av_bld', dpsi_orb_av_bld)
!  call catalog_one_routine ('dpsi_orb_dpsi_orb_av_bld', dpsi_orb_dpsi_orb_av_bld)
!  call catalog_one_routine ('dpsi_orb_eloc_av_bld', dpsi_orb_eloc_av_bld)
!  call catalog_one_routine ('dpsi_orb_sq_eloc_ex_av_bld', dpsi_orb_sq_eloc_ex_av_bld)
!  call catalog_one_routine ('dpsi_orb_sq_av_bld', dpsi_orb_sq_av_bld)

! writing routines (not nodes!)
  call catalog_one_routine ('intra_wrt', intra_wrt)
  call catalog_one_routine ('intra_3d_wrt', intra_3d_wrt)
  call catalog_one_routine ('dens_wrt', dens_wrt)
  call catalog_one_routine ('dens_xy_z_wrt', dens_xy_z_wrt)
  call catalog_one_routine ('dens_3d_wrt', dens_3d_wrt)
  call catalog_one_routine ('dens_fourier_wrt', dens_fourier_wrt)
  call catalog_one_routine ('extracule_wrt', extracule_wrt)
  call catalog_one_routine ('dipole_moment_wrt', dipole_moment_wrt)
  call catalog_one_routine ('forces_wrt', forces_wrt)
  call catalog_one_routine ('vb_weights_wrt', vb_weights_wrt)

!!pjas WAS
  call catalog_one_node ('pjas_init_bld', pjas_init_bld)
  call catalog_one_node ('deloc_pjasen_bld',deloc_pjasen_bld)
!  call catalog_one_node ('d2psi_pjase_bld',  d2psi_pjase_bld)
  call catalog_one_node ('deloc_pjasee_bld',deloc_pjasee_bld)
  call catalog_one_node ('dpsi_pjas_bld',  dpsi_pjas_bld)
  call catalog_one_node ('deloc_pjas_bld',deloc_pjas_bld)
  call catalog_one_node ('deloc_pot_nloc_pjas_bld',deloc_pot_nloc_pjas_bld)

 end subroutine catalog_all_nodes_and_routines

end module catalog_routines_mod
