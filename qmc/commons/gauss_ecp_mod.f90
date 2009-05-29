module gauss_ecp_mod

 implicit none
 save
 
 integer, allocatable :: necp_term(:,:),necp_power(:,:,:)
 double precision, allocatable :: ecp_coef(:,:,:),ecp_exponent(:,:,:)

end module gauss_ecp_mod
