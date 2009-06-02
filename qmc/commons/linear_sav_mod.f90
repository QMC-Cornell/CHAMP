module linear_sav_mod
 implicit none
 save

 double precision, allocatable :: ham_sav(:,:),ovlp_sav(:,:),renorm_ovlp(:)

end module linear_sav_mod
