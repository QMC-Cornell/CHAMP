module numbas_mod

 use constants_mod
 implicit none
 save

 integer numr
 double precision, allocatable :: exp_h_bas(:),r0_bas(:),rwf(:,:,:,:),d2rwf(:,:,:,:)
 integer, allocatable :: nrbas(:),nrbas_analytical(:),nrbas_numerical(:),igrid(:),nr(:)
 integer, allocatable :: iwrwf(:,:)
 integer :: MRWF_PTS
 integer :: MRWF = 0 ! must be initialized to zero, actual value is calculated at execution

end module numbas_mod
