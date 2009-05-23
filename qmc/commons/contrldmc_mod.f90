module contrldmc_mod

 use constants_mod
 implicit none
 save

 double precision tau,rttau,taueff(MFORCE),tautot
 integer nfprod,idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e

end module contrldmc_mod
