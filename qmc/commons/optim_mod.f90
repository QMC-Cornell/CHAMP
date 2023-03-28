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

 integer, allocatable  :: ipivot(:)
 integer               :: noutput, nstep_fit, ibold
 real(dp)              :: pmarquardt, tau_fit
 logical               :: analytic_jacobian, cholesky
!      read(5,*) (ipivot(j),j=1,norb)
!      read(5,*) eguess
!      read(5,*) pmarquardt,tau,noutput,nstep,ibold
!      read(5,'(2l2)') analytic,cholesky

end module optim_mod
