module optim_mod

 use constants_mod

 integer, allocatable  :: lo(:),npoint(:)
 integer, allocatable  :: iwjasa(:,:),iwjasb(:,:),iwjasc(:,:)
 integer, allocatable  :: iwjasf(:,:),iwbase(:),iwbasi(:),iworb(:)
 integer, allocatable  :: iwcsf(:),iebase(:,:),iebasi(:,:),ieorb(:,:)
 integer, allocatable  :: imnbas(:)
 integer               :: nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj
 integer, allocatable  :: nparma(:),nparmb(:),nparmc(:),nparmf(:)
 integer               :: necn,nebase
 integer               :: nparmd !JT add nparmd to replace MPARMD
 integer               :: nparmjs !JT

end module optim_mod
