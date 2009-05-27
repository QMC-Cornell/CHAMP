module optim_mod

 use constants_mod

 integer               :: lo(MORB),npoint(MORB)
 integer               :: iwjasa(MPARMJ,NCTYP3X),iwjasb(MPARMJ,3),iwjasc(MPARMJ,MCTYPE)
 integer, allocatable  :: iwjasf(:,:),iwbase(:),iwbasi(:),iworb(:)
 integer, allocatable  :: iwcsf(:),iebase(:,:),iebasi(:,:),ieorb(:,:)
 integer, allocatable  :: imnbas(:)
 integer               :: nparml,nparme,nparmcsf,nparms,nparmg,nparm_read,nparmj
 integer               :: nparma(NCTYP3X),nparmb(3),nparmc(MCTYPE),nparmf(MCTYPE)
 integer               :: necn,nebase
 integer               :: nparmd !JT add nparmd to replace MPARMD

end module optim_mod
