module contrldmc_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: taueff(:)
 double precision tau,rttau,rtrttau,tautot
 double precision taunow ! needed for tmoves
 double precision pow_rewt
 integer :: nfprod,idmc,ipq,itau_eff,itau_integ=1,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
 logical tmoves
 logical :: l_modified_adrift, l_method_of_images, l_tau_diffusion
 logical :: l_psi_approx, l_psi_approx_norm_old

end module contrldmc_mod
