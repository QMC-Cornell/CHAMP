module cslater_cf_mod

  use types_mod

  complex(dpc), allocatable ::  cslmui(:,:),cslmdi(:,:)
  complex(dpc), allocatable ::  cfpu(:,:,:,:),cfpd(:,:,:,:)
  complex(dpc), allocatable ::  cdetu(:),cdetd(:)
  complex(dpc), allocatable ::  cddeti_deti(:,:,:),cd2edeti_deti(:,:)

end module cslater_cf_mod
