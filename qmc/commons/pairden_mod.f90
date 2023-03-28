module pairden_mod

 implicit none
 save
 
 double precision, allocatable :: xx0probut(:,:,:),xx0probuu(:,:,:),xx0probud(:,:,:),xx0probdt(:,:,:)
 double precision, allocatable :: xx0probdu(:,:,:),xx0probdd(:,:,:),den2d_t(:,:),den2d_d(:,:),den2d_u(:,:)
 double precision, allocatable :: pot_ee2d_t(:,:), pot_ee2d_u(:,:), pot_ee2d_d(:,:)
 double precision, allocatable :: dos(:)
 double precision              :: dos_dele
 double precision              :: delxi(3),xmax,xfix(3)
 integer ifixe
 integer, parameter :: NAX = 50

end module pairden_mod
