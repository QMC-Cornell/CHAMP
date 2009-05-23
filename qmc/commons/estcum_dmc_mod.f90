module estcum_dmc_mod 

 use constants_mod
 implicit none
 save

 double precision wcum,w_acc_cum,wfcum,wgcum(MFORCE),wg_acc_cum,wdcum
 double precision wgdcum, wcum1,w_acc_cum1,wfcum1,wgcum1(MFORCE),wg_acc_cum1     
 double precision wdcum1, ecum,efcum,egcum(MFORCE),ecum1,efcum1,egcum1(MFORCE)
 double precision ei1cum,ei2cum,ei3cum, pecum(MFORCE),peicum(MFORCE),tpbcum(MFORCE),tjfcum(MFORCE),r2cum
 double precision ricum,taucum(MFORCE)

end module estcum_dmc_mod 
