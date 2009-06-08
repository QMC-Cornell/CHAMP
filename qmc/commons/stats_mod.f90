module stats_mod

 use constants_mod
 implicit none
 save

 double precision dfus2ac,dfus2un,dr2ac,dr2un,acc,acc_int,try_int
 integer nbrnch,nodecr
 double precision, allocatable :: dfus2unf(:)

end module stats_mod
