module jaspar3_mod

 use constants_mod
 implicit none
 save

 integer nord
 double precision, allocatable :: a(:,:),b(:,:,:),c(:,:,:),fck(:,:,:),scalek(:)
 double precision, allocatable :: b_sav(:,:),c_sav(:,:)
 double precision, allocatable :: b_best(:,:),c_best(:,:)
 double precision              :: scalek_best
           
end module jaspar3_mod
