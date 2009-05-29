module jaspar4_mod

 use constants_mod
 implicit none
 save

 integer norda,nordb,nordc
 double precision, allocatable :: a4(:,:,:)

 double precision, allocatable :: a4_sav(:,:)
 double precision, allocatable :: a4_best(:,:)

end module jaspar4_mod
