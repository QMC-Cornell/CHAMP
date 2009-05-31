module periodic_mod

 use constants_mod
 implicit none
 save

 integer igvec(3,NGVEC_BIGX),igmult(NGNORM_BIGX),igvec_sim(3,NGVEC_SIM_BIGX)
 integer, allocatable :: ireal_imag(:)
 integer igmult_sim(NGNORM_SIM_BIGX),kvec(3,MKPTS),k_inv(MKPTS),nband(MKPTS)
 integer ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
 integer ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
 integer ng1d(3),ng1d_sim(3),npoly,ncoef,np,isrange
 
 double precision rlatt(3,3),glatt(3,3),rlatt_sim(3,3),glatt_sim(3,3)
 double precision rlatt_inv(3,3),glatt_inv(3,3),rlatt_sim_inv(3,3),glatt_sim_inv(3,3)
 double precision cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
 double precision gvec(3,NGVEC_BIGX),gnorm(NGNORM_BIGX)
 double precision gvec_sim(3,NGVEC_SIM_BIGX),gnorm_sim(NGNORM_SIM_BIGX)
 double precision rkvec_shift(3),rkvec(3,MKPTS),rknorm(MKPTS),znuc_sum,znuc2_sum,vcell,vcell_sim

end module periodic_mod
