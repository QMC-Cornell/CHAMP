module jasow_mod

! These arrays are for saving pieces of the Jastrow to avoid recalculating unnecessary terms when 1 electron is moved
! fsumow(MWALK)
! fsow(nelec,nelec,MWALK)
! fijow(3,nelec,nelec,MWALK)
! fjow(3,nelec,MWALK)
  
 implicit none
 save
 
 double precision, allocatable :: fsumow(:),fsow(:,:,:),fijow(:,:,:,:),fjow(:,:,:)

end module jasow_mod
