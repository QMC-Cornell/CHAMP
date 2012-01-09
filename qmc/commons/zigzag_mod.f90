module zigzag_mod

 implicit none
 save
 
 double precision, allocatable :: zzpairden_t(:,:), zzpariden_u(:,:), zzpairden_d(:,:)
 double precision, allocatable :: zzdenij_t(:,:), zzdenij_u(:,:), zzdenij_d(:,:)
 double precision, allocatable :: xold_sav(:,:), xnew_sav(:,:), zzposold(:,:), zzposnew(:,:)
 integer, allocatable          :: iold_indices(:), inew_indices(:)
 double precision              :: zzcum, zzsum, zzcm2, zz2cum, zz2sum, zz2cm2
 integer                       :: izigzag

end module zigzag_mod
