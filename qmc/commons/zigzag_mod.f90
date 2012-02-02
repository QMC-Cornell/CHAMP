module zigzag_mod

 implicit none
 save
 
 double precision, allocatable :: zzpairden_t(:,:)
 double precision, allocatable :: zzdenij_t(:,:)
 double precision, allocatable :: zzcorrij(:), zzcorr(:)
 double precision, allocatable :: xold_sav(:,:), xnew_sav(:,:), zzposold(:,:), zzposnew(:,:)
 integer, allocatable          :: iold_indices(:), inew_indices(:)
 double precision, allocatable :: zzcum(:), zzsum(:), zzcm2(:)
 integer                       :: izigzag
 integer, parameter            :: nzzvars = 3
end module zigzag_mod
