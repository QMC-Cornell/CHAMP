module vardep_mod

 implicit none
 save

 integer, allocatable :: nvdepend(:,:),iwdepend(:,:,:)
 double precision, allocatable :: cdep(:,:,:)

end module vardep_mod
