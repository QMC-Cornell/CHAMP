module gamma_mod
 implicit none
 save

 double precision, allocatable, target :: gup(:,:), gdn(:,:)
 double precision, allocatable, target :: psp_nonloc_orb(:,:)
 double precision, allocatable, target :: psp_nonloc_pot
 integer, allocatable, target ::  occup(:),  occdn(:)
 integer, allocatable, target :: ioccup(:), ioccdn(:)
 integer                      :: noccup, noccdn

end module gamma_mod
