module dets_mod

 use constants_mod

 integer               :: nup, ndn, nupdn
 integer               :: iantiferromagnetic
 integer               :: nup_square, ndn_square, nupdn_square
 integer               :: ncsf,ndet
 double precision, allocatable :: csf_coef(:,:),cdet_in_csf(:,:)
 integer,  allocatable :: ndet_in_csf(:),iwdet_in_csf(:,:)

 double precision, allocatable :: csf_coef_sav(:)
 double precision, allocatable :: csf_coef_best(:)

end module dets_mod
