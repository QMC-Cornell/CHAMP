module confg_mod

 implicit none
 save

 integer :: ndata
 double precision :: eguess,wghtsm,cuspwt
 double precision, allocatable :: x(:,:,:),psid(:),psij(:)
 double precision, allocatable :: psio(:),eold(:),uwdiff(:),wght(:),dvpdv(:)

end module confg_mod
