module jacobsave_mod

 use constants_mod
 implicit none
 save

 double precision              :: ajacob
 double precision, allocatable :: ajacold(:,:)

end module jacobsave_mod
