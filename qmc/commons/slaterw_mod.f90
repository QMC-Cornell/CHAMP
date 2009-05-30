module slaterw_mod

 implicit none
 save

 double precision, allocatable :: slmuiw(:,:,:),slmdiw(:,:,:),fpuw(:,:,:,:),fpdw(:,:,:,:)
 double precision, allocatable :: detuw(:,:),detdw(:,:),ddeti_detiw(:,:,:,:)

end module slaterw_mod
