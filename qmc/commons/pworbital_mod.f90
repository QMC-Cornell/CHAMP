module pworbital_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: c_rp(:,:),c_rm(:,:),c_ip(:,:),c_im(:,:)
! integer, allocatable :: ngorb(MORB),isortg(NGVECX,MORB),isortk(MKPTS) ! JT not used
 integer              :: icmplx

end module pworbital_mod
