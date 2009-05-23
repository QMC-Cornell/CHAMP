module config_dmc_mod

 use constants_mod
 implicit none
 save

 double precision xoldw(3,MELEC,MWALK,MFORCE),voldw(3,MELEC,MWALK,MFORCE)
 double precision psidow(MWALK,MFORCE),psijow(MWALK,MFORCE),peow(MWALK,MFORCE),peiow(MWALK,MFORCE),d2ow(MWALK,MFORCE)

end module config_dmc_mod
