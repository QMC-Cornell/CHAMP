module zigzag_mod

 implicit none
 save
 
 double precision, allocatable :: zzpairden_t(:,:), zzpariden_u(:,:), zzpairden_d(:,:)
 double precision, allocatable :: zzdenij_t(:,:), zzdenij_u(:,:), zzdenij_d(:,:)
 double precision              :: zigzag_avg, zigzag_abs, zigzag2_avg
 integer                       :: izigzag

end module zigzag_mod
