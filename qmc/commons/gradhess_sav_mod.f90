module gradhess_sav_mod
 implicit none
 save

 double precision, allocatable :: hess_sav(:,:),grad_sav(:)

end module gradhess_sav_mod
