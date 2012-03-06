module qua_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: xq0(:),yq0(:),zq0(:)          
 double precision, allocatable :: xq(:),yq(:),zq(:),wq(:)
 integer :: nquad
 integer :: MPS_QUAD = 86
 logical :: l_do_tmoves=.false.

end module qua_mod
