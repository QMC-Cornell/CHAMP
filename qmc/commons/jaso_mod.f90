module jaso_mod

! These are for saving pieces of the Jastrow to avoid recalculating unnecessary terms when 1 electron is moved
! fsumo
! d2o
! fso(nelec,nelec)
! fijo(3,nelec,nelec)
! fjo(3,nelec)
! d2ijo(nelec,nelec)

 use constants_mod
 implicit none
 save

! double precision fsumo,d2o
! double precision, allocatable :: fso(:,:),fijo(:,:,:),fjo(:,:),d2ijo(:,:)
! double precision, allocatable :: lapjijo(:,:),lapjo(:)
 double precision, pointer :: fsumo,d2o
 double precision, pointer :: fso(:,:),fijo(:,:,:),fjo(:,:),d2ijo(:,:)
 double precision, pointer :: lapjijo(:,:),lapjo(:)

end module jaso_mod
