module gradhessdero_mod

 implicit none
 save

 double precision, allocatable :: deti_det_old(:),gvalue_old(:),denergy_old(:)
 double precision, allocatable :: d1d2a_old(:),d2d2a_old(:),d1d2b_old(:),d2d2b_old(:),didk_old(:)
 double precision, allocatable :: detij_det_old(:,:)

end module gradhessdero_mod
