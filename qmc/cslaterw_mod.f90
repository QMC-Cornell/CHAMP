module cslaterw_mod

  use types_mod

  complex(dpc), allocatable ::  cslmuiw(:,:,:),cslmdiw(:,:,:)
  complex(dpc), allocatable ::  cfpuw(:,:,:,:),cfpdw(:,:,:,:)
  complex(dpc), allocatable ::  cdetuw(:,:),cdetdw(:,:)
  complex(dpc), allocatable ::  cddeti_detiw(:,:,:,:)

end module cslaterw_mod
