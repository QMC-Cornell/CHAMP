module estcum_dmc_mod 

 use constants_mod
 implicit none
 save

 double precision, allocatable :: wgcum(:),wgcum1(:),egcum(:),egcum1(:),pecum(:),peicum(:),tpbcum(:),tjfcum(:),taucum(:)
 double precision, allocatable :: ovlp_ovlp_fn_cum(:,:)
 double precision wcum,w_acc_cum,wfcum,wg_acc_cum,wdcum
 double precision wgdcum, wcum1,w_acc_cum1,wfcum1,wg_acc_cum1     
 double precision wdcum1, ecum,efcum,ecum1,efcum1
 double precision ei1cum,ei2cum,ei3cum,r1cum,r2cum,r3cum,r4cum
 double precision ricum

end module estcum_dmc_mod 
