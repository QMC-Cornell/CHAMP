module dorb_mod

 use constants_mod
 implicit none
 save

 integer, allocatable :: iworbd(:,:),iworbdup(:,:),iworbddn(:,:)
 integer, allocatable :: iwdetup(:),iwdetdn(:)
 integer ndetup, ndetdn, ndetupdn

end module dorb_mod
