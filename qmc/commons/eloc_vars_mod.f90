module eloc_vars_mod
  use types_mod
  implicit none
  real(dp)                       :: eloc
  real(dp)                       :: eloc_vmc
  real(dp)                       :: eloc_dmc
  real(dp)                       :: eloc_av
  real(dp)                       :: eloc_av_err
  real(dp)                       :: eloc_bav
  real(dp)                       :: eloc_av_var
  real(dp)                       :: eloc_sq
  real(dp)                       :: eloc_sq_av
  real(dp)                       :: eloc_var
  real(dp)                       :: sigma
  real(dp)                       :: error_sigma
  real(dp)                       :: eloc_kin
  real(dp)                       :: eloc_kin_av
  real(dp)                       :: eloc_kin_av_err
  real(dp)                       :: eloc_pot
  real(dp)                       :: eloc_pot_loc
  real(dp)                       :: eloc_pot_nloc
  real(dp)                       :: pe_en
  real(dp)                       :: pe_ee
  real(dp)                       :: eloc_pot_en
  real(dp)                       :: eloc_pot_en_av
  real(dp)                       :: eloc_pot_ee
  real(dp)                       :: eloc_pot_ee_av
  real(dp)                       :: eloc_pot_ee_av_err
  real(dp)                       :: eloc_pot_ee_zv
  real(dp)                       :: eloc_pot_ee_zv_av
  real(dp)                       :: eloc_pot_ee_zv_av_err
  real(dp)                       :: eloc_kin_jas
  real(dp)                       :: eloc_kin_jas_av
  real(dp)                       :: eloc_kin_jas_av_err
  real(dp)                       :: eloc_kin_jas_pot_ee_av
  real(dp)                       :: eloc_kin_jas_pot_ee_av_err
end module eloc_vars_mod 
