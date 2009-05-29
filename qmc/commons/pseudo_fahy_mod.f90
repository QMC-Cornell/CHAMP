module pseudo_fahy_mod

 implicit none
 save

 double precision, allocatable :: potl(:,:),ptnlc(:,:,:),dradl(:),drad(:),rcmax(:)
 integer, allocatable :: npotl(:),nlrad(:)

end module pseudo_fahy_mod
