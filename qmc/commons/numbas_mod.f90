module numbas_mod

 use constants_mod
 implicit none
 save

 integer numr
 double precision, allocatable :: exp_h_bas(:),r0_bas(:),rwf(:,:,:,:),d2rwf(:,:,:,:)
 integer, allocatable :: nrbas(:),igrid(:),nr(:)
 integer, allocatable :: iwrwf(:,:)
 integer :: MRWF_PTS
 integer :: MRWF

end module numbas_mod
