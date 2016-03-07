module catalog_averages_mod

  use basic_tools_mod
  use constants_mod
  use variables_mod
  use objects_mod
  use average_mod

  contains

! ===================================================================================
  subroutine define_averages_and_errors
! -----------------------------------------------------------------------------------
! Description   : define all averages and errors
!
! Created       : J. Toulouse, 15 Jan 2006
! -----------------------------------------------------------------------------------
  implicit none

! begin
  call object_error_define ('e_loc_av','e_loc_av_err')
  call object_average_define ('eloc_sq','eloc_sq_av')
  call object_average_define ('eloc_kin','eloc_kin_av')
  call object_error_define ('eloc_kin_av','eloc_kin_av_err')

! general derivatives
  call object_average_define ('dpsi_sq', 'dpsi_sq_av')
  call object_error_define ('dpsi_dpsi_av', 'dpsi_dpsi_av_err')
  call object_average_define ('dpsi_deloc', 'dpsi_deloc_av')
  call object_average_define ('dpsi_sq_eloc', 'dpsi_sq_eloc_av')
  call object_average_define ('dpsi_dpsi_eloc', 'dpsi_dpsi_eloc_av')
  call object_average_define ('d2psi', 'd2psi_av')
  call object_average_define ('d2psi_eloc', 'd2psi_eloc_av')
  call object_average_define ('dpsi_orb', 'dpsi_orb_av')
  call object_error_define   ('dpsi_orb_av', 'dpsi_orb_err')
  call object_average_define ('dpsi_orb_dpsi_orb','dpsi_orb_dpsi_orb_av')
  call object_average_define ('dpsi_orb_eloc', 'dpsi_orb_eloc_av')

! orbital derivatives
  call object_average_define ('dpsi_orb_sq', 'dpsi_orb_sq_av')
  call object_average_define ('dpsi_orb_sq_eloc_ex','dpsi_orb_sq_eloc_ex_av')

  call object_average_define ('dpsi_orb_sq_eloc','dpsi_orb_sq_eloc_av')
  call object_average_define ('dpsi_orb_dpsi_orb_eloc_ex','dpsi_orb_dpsi_orb_eloc_ex_av')
  call object_error_define ('e_ex','e_ex_err')
  call object_error_define ('e_ex_2','e_ex_2_err')
  call object_average_define ('dpsi_orb_deloc_orb','dpsi_orb_deloc_orb_av')
  call object_error_define   ('dpsi_orb_deloc_orb_av','dpsi_orb_deloc_orb_av_err')
  call object_error_define ('dpsi_orb_sq_eloc_av','dpsi_orb_sq_eloc_av_err')
  call object_error_define ('grad_orb_norm','grad_orb_norm_err')
  call object_error_define ('grad_orb', 'grad_orb_err')
  call object_error_define ('delta_e_orb', 'delta_e_orb_err')

! jastrow derivatives
  call object_average_define ('dpsi_jas', 'dpsi_jas_av')
  call object_average_define ('d2psi_jas', 'd2psi_jas_av')
  call object_average_define ('d2gvalue', 'd2gvalue_av')
  call object_average_define ('dpsi_jas_sq', 'dpsi_jas_sq_av')
  call object_average_define ('dpsi_jas_eloc', 'dpsi_jas_eloc_av')
  call object_average_define ('dpsi_jas_sq_eloc', 'dpsi_jas_sq_eloc_av')
  call object_average_define ('dpsi_jas_deloc_jas', 'dpsi_jas_deloc_jas_av')
  call object_average_define ('deloc_jas', 'deloc_jas_av')
  call object_average_define ('d2eloc_jas', 'd2eloc_jas_av')

! csf derivatives
  call object_average_define ('dpsi_csf', 'dpsi_csf_av')
  call object_average_define ('dpsi_csf_dpsi_csf', 'dpsi_csf_dpsi_csf_av')
  call object_average_define ('dpsi_csf_sq', 'dpsi_csf_sq_av')
  call object_average_define ('dpsi_csf_sq_eloc', 'dpsi_csf_sq_eloc_av')
  call object_average_define ('dpsi_csf_deloc_csf', 'dpsi_csf_deloc_csf_av')
  call object_average_define ('deloc_csf', 'deloc_csf_av')
  call object_average_define ('dpsi_csf_eloc', 'dpsi_csf_eloc_av')

! newton method
  call object_error_define ('hess_nwt', 'hess_nwt_err')
  call object_error_define ('hess_nwt_norm', 'hess_nwt_norm_err')
  call object_error_define ('gradient_norm', 'gradient_norm_err')

! linear method
  call object_error_define ('ham_lin', 'ham_lin_err')

! perturbative method

! intracule
  call object_average_define ('intra_sp_histo', 'intra_sp_histo_av')
  call object_error_define ('intra_sp_histo_av', 'intra_sp_histo_av_err')
  call object_average_define ('intra_sp_zv1', 'intra_sp_zv1_av')
  call object_error_define ('intra_sp_zv1_av', 'intra_sp_zv1_av_err')
  call object_average_define ('intra_sp_zv2', 'intra_sp_zv2_av')
  call object_error_define ('intra_sp_zv2_av', 'intra_sp_zv2_av_err')
  call object_average_define ('intra_sp_zv3', 'intra_sp_zv3_av')
  call object_error_define ('intra_sp_zv3_av', 'intra_sp_zv3_av_err')
  call object_average_define ('intra_sp_zv4', 'intra_sp_zv4_av')
  call object_error_define ('intra_sp_zv4_av', 'intra_sp_zv4_av_err')
  call object_average_define ('intra_sp_zvzb3', 'intra_sp_zvzb3_av')
  call object_error_define ('intra_sp_zvzb3_av', 'intra_sp_zvzb3_av_err')
  call object_average_define ('intra_sp_zvzb4', 'intra_sp_zvzb4_av')
  call object_error_define ('intra_sp_zvzb4_av', 'intra_sp_zvzb4_av_err')
  call object_average_define ('intra_sp_zv1zb3', 'intra_sp_zv1zb3_av')
  call object_error_define ('intra_sp_zv1zb3_av', 'intra_sp_zv1zb3_av_err')
  call object_average_define ('intra_3d_histo', 'intra_3d_histo_av')
  call object_error_define ('intra_3d_histo_av', 'intra_3d_histo_av_err')
  call object_average_define ('intra_3d_zv1', 'intra_3d_zv1_av')
  call object_error_define ('intra_3d_zv1_av', 'intra_3d_zv1_av_err')
  call object_average_define ('intra_3d_zv2', 'intra_3d_zv2_av')
  call object_error_define ('intra_3d_zv2_av', 'intra_3d_zv2_av_err')

! extracule
  call object_average_define ('extracule_zv1', 'extracule_zv1_av')
  call object_error_define ('extracule_zv1_av', 'extracule_zv1_av_err')

 end subroutine define_averages_and_errors

end module catalog_averages_mod
