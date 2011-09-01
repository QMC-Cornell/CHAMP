module config_dmc_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: xoldw(:,:,:,:),voldw(:,:,:,:)
 double precision, allocatable :: psidow(:,:),psijow(:,:),peow(:,:),peiow(:,:),d2ow(:,:)
 double precision, allocatable :: pot_ee_oldw(:,:,:)

end module config_dmc_mod
