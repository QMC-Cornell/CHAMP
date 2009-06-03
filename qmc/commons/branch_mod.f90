module branch_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: wtgen(:),ff(:),eoldw(:,:)
 double precision, allocatable :: pwt(:,:),wthist(:,:,:), wt(:)
 double precision eigv,eest,wdsumo,wgdsumo,fprod
 integer nwalk

end module branch_mod
