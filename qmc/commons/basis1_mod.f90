module basis1_mod

 use constants_mod

 real(dp)              :: betaq
 real(dp), allocatable :: zex(:,:)
 integer, allocatable  :: n1s(:),n2s(:),n2p(:,:)
 integer, allocatable  :: n3s(:),n3p(:,:),n3d(:,:)
 integer, allocatable  :: n4s(:),n4p(:,:),n4d(:,:),n4f(:,:)
 integer, allocatable  :: n5s(:),n5p(:,:),n5d(:,:),n5f(:,:),n5g(:,:)
 integer, allocatable  :: n6d(:,:),n6f(:,:),n6g(:,:),n6h(:,:)
 integer, allocatable  :: n7g(:,:),n7h(:,:),n7i(:,:)
 integer, allocatable  :: n8i(:,:),n8j(:,:)
 integer, allocatable  :: n9k(:,:)
 integer, allocatable  :: n10l(:,:)
 integer, allocatable  :: n11m(:,:)
 integer, allocatable  :: n12n(:,:)
 integer, allocatable  :: n13o(:,:)
 integer, allocatable  :: nsa(:),npa(:,:),nda(:,:)

end module basis1_mod
