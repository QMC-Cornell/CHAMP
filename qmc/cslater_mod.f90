module cslater_mod

  use types_mod
  complex(dpc), allocatable ::  cslmui(:,:),cslmdi(:,:)
  complex(dpc), allocatable ::  cfpu(:,:,:),cfpd(:,:,:)
  complex(dpc), allocatable ::  cfppu(:,:),cfppd(:,:)
  complex(dpc), allocatable ::  cdetu(:),cdetd(:)
  complex(dpc), allocatable ::  cddeti_deti(:,:,:),cd2edeti_deti(:,:)
  complex(dpc), allocatable ::  cdeti_det(:),cddeti_det(:,:,:),cd2deti_det(:)
  complex(dpc)                  cd2det_det

end module cslater_mod
