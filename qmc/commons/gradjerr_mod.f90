module gradjerr_mod

 implicit none
 save
 
 double precision              :: e_bsum
 double precision, allocatable :: grad_jas_bcum(:),grad_jas_bcm2(:),dj_e_bsum(:),dj_bsum(:),dj_e_save(:),dj_save(:)

end module gradjerr_mod
