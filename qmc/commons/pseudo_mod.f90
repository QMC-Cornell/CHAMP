module pseudo_mod

 use constants_mod
 implicit none
 save

 integer nloc
 double precision, allocatable ::  vps(:,:,:)
! double precision vpso(MELEC,MCENT,MPS_L,MFORCE)
 integer, allocatable ::  npotd(:), lpotp1(:)
     
end module pseudo_mod
