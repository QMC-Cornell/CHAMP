module pseudo_mod

 use constants_mod
 implicit none
 save

 integer nloc
 double precision, allocatable ::  vps(:,:,:)
 integer, allocatable ::  npotd(:), lpotp1(:)
     
end module pseudo_mod
