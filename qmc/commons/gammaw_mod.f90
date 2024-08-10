module gammaw_mod
 implicit none
 save

 double precision, allocatable, target :: orbw(:,:,:), dorbw(:,:,:,:), ddorbw(:,:,:)
 double precision, allocatable, target :: aiupw(:,:,:), aidnw(:,:,:)
 double precision, allocatable, target :: deta_upw(:), deta_dnw(:)
 double precision, allocatable, target :: tupw(:,:,:), tdnw(:,:,:)
 double precision, allocatable, target :: yupw(:,:,:), ydnw(:,:,:)
 double precision, allocatable, target :: chiw(:)
 double precision, allocatable, target :: detupw(:,:), detdnw(:,:)
 double precision, allocatable, target :: invupw(:,:,:), invdnw(:,:,:)

end module gammaw_mod
