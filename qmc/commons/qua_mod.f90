module qua_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: xq0(:),yq0(:),zq0(:)          
 double precision, allocatable :: xq(:),yq(:),zq(:),wq(:)
 double precision, allocatable :: quadr(:,:,:), quadx(:,:,:,:)
 integer :: nquad
 integer :: MPS_QUAD = 86
 logical :: l_do_tmoves=.false.
 integer :: iaccept_tmove

end module qua_mod
