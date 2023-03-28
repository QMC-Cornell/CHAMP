module periodic_mod

 use constants_mod
 implicit none
 save

 integer, allocatable :: igvec(:,:),igmult(:),igvec_sim(:,:)
 integer, allocatable :: ireal_imag(:)
 integer, allocatable :: igmult_sim(:)
 integer, allocatable :: kvec(:,:),k_inv(:),nband(:)
 integer ngnorm,ngvec,ngnorm_sim,ngvec_sim,ngnorm_orb,ngvec_orb,nkvec
 integer ngnorm_big,ngvec_big,ngnorm_sim_big,ngvec_sim_big
 integer, allocatable :: ng1d(:),ng1d_sim(:)
 integer npoly,ncoef,np,isrange
 
 double precision, allocatable :: rlatt(:,:),glatt(:,:),rlatt_sim(:,:),glatt_sim(:,:)
 double precision, allocatable :: rlatt_inv(:,:),glatt_inv(:,:),rlatt_sim_inv(:,:),glatt_sim_inv(:,:)
 double precision cutr,cutr_sim,cutg,cutg_sim,cutg_big,cutg_sim_big
 double precision, allocatable :: gvec(:,:),gnorm(:)
 double precision, allocatable :: gvec_sim(:,:),gnorm_sim(:)
 double precision znuc_sum,znuc2_sum,vcell,vcell_sim
 double precision, allocatable :: rkvec_shift(:),rkvec(:,:),rknorm(:)

end module periodic_mod
