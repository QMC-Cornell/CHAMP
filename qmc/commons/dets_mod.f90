module dets_mod

 use constants_mod

 integer               :: nup, ndn, nupdn
 integer               :: nup_square, ndn_square, nupdn_square
 integer               :: ncsf,ndet
 real(dp), allocatable :: csf_coef(:,:),cdet_in_csf(:,:)
 integer,  allocatable :: ndet_in_csf(:),iwdet_in_csf(:,:)

end module dets_mod
