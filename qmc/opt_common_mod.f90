module opt_common_mod

  use all_tools_mod

! Declaration of global variables and default values
  real(dp)                :: energy_sav
  real(dp)                :: energy_err_sav
  real(dp)                :: energy_sigma_sav
  real(dp)                :: error_sigma_sav
  real(dp)                :: ene_var_sav

  real(dp)                        :: psi_lin_var_norm = 0.d0
  real(dp)                        :: psi_lin_var_norm_max = 10.d0

end module opt_common_mod
