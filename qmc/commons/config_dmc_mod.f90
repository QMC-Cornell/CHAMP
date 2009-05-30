module config_dmc_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: xoldw(:,:,:,:),voldw(:,:,:,:)
 double precision, allocatable :: psidow(:,:),psijow(:,:),peow(:,:),peiow(:,:),d2ow(:,:)

end module config_dmc_mod
