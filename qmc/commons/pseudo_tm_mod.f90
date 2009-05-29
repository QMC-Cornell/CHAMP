module pseudo_tm_mod

 implicit none
 save 

 double precision, allocatable :: rmax_coul(:),rmax_nloc(:),exp_h_ps(:),r0_ps(:),vpseudo(:,:,:),d2pot(:,:,:)
 integer, allocatable :: igrid_ps(:), nr_ps(:)

end module pseudo_tm_mod
