module distance_mod

 use constants_mod
 implicit none
 save

 integer                       :: iring_coulomb
 double precision, allocatable :: rshift(:,:,:)!,rvec_en(:,:,:),r_en(:,:),rvec_ee(:,:),r_ee(:)
 double precision, pointer :: rvec_en(:,:,:),r_en(:,:),rvec_ee(:,:),r_ee(:)
 double precision, allocatable :: pot_ee(:)

end module distance_mod
