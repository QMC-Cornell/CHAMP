module jasn_mod

! These are for saving pieces of the Jastrow to avoid recalculating unnecessary terms when 1 electron is moved
! fsumn
! d2n
! fsn(nelec,nelec)
! fijn(3,nelec,nelec)
! fjn(3,nelec)
! d2ijn(nelec,nelec)

 use constants_mod
 implicit none
 save

! double precision :: fsumn,d2n
! double precision, allocatable :: fsn(:,:),fijn(:,:,:),fjn(:,:),d2ijn(:,:)
! double precision, allocatable :: lapjijn(:,:),lapjn(:)
 double precision, pointer :: fsumn,d2n
 double precision, pointer :: fsn(:,:),fijn(:,:,:),fjn(:,:),d2ijn(:,:)
 double precision, pointer :: lapjijn(:,:),lapjn(:)

end module jasn_mod
