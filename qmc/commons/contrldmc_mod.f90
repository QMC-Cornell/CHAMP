module contrldmc_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: taueff(:)
 double precision tau,rttau,rtrttau,tautot
 double precision taunow ! needed for tmoves
 integer :: nfprod,idmc,ipq,itau_eff,itau_integ=1,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
 logical tmoves
 logical l_improved_gf

end module contrldmc_mod
