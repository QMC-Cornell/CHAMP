module estsum_dmc_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: wgsum(:),wsum1(:),wgsum1(:),egsum(:),esum1(:),egsum1(:)
 double precision, allocatable :: pesum(:),peisum(:),tpbsum(:),tjfsum(:),tausum(:) 
 double precision, allocatable :: ovlp_ovlp_fn_sum(:,:)
 double precision wsum,w_acc_sum,wfsum,wg_acc_sum,wdsum   
 double precision wgdsum, w_acc_sum1,wfsum1,wg_acc_sum1
 double precision wdsum1, esum,efsum,efsum1
 double precision ei1sum,ei2sum,ei3sum,r1sum,r2sum,r3sum,r4sum
 double precision risum

end module estsum_dmc_mod
