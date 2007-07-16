module cslater_mod

  complex*16, allocatable ::  cslmui(:,:),cslmdi(:,:)
  complex*16, allocatable ::  cfpu(:,:,:),cfpd(:,:,:)
  complex*16, allocatable ::  cfppu(:,:),cfppd(:,:)
  complex*16, allocatable ::  cdetu(:),cdetd(:)
  complex*16, allocatable ::  cddeti_deti(:,:,:),cd2edeti_deti(:,:)
  complex*16, allocatable ::  cdeti_det(:),cddeti_det(:,:,:),cd2deti_det(:)
  complex*16                  cd2det_det

end module cslater_mod
