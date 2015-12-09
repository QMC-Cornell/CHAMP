module pjas_mod

 use constants_mod
 implicit none
 save

 integer                                :: param_pjas_nb
 real(dp), allocatable                  :: pjas_parms(:,:), pjas_parms_best(:)

end module pjas_mod
