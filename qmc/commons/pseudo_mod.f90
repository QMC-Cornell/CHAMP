module pseudo_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable ::  vps(:,:,:)
! double precision vpso(MELEC,MCENT,MPS_L,MFORCE)
 integer npotd(MCTYPE),lpotp1(MCTYPE),nloc
     
end module pseudo_mod
