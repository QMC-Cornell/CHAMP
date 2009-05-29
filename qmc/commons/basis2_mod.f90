module basis2_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: zex2(:,:,:)
 integer, allocatable :: n_bas(:), l_bas(:), m_bas(:)
 integer, allocatable :: icenter_basis(:), ictype_basis(:)
 integer, allocatable :: nbasis_ctype(:), n_bas2(:,:), iwrwf2(:)
 integer              :: mbasis_ctype !JT

end module basis2_mod
