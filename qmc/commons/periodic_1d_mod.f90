module periodic_1d_mod

 use constants_mod
 implicit none
 save

 integer ngvecs_1d
 
 double precision, allocatable :: gvec_1d(:), gamma_gvec(:)

 double precision alattice, ewald_1d_cutoff

end module periodic_1d_mod
