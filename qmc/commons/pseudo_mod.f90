module pseudo_mod

 use constants_mod
 implicit none
 save

 integer nloc
 double precision, allocatable ::  vps(:,:,:)
 integer, allocatable ::  npotd(:), lpotp1(:)
 integer :: MPS_L=4
 integer :: MPS_GRID=3501
     
end module pseudo_mod
