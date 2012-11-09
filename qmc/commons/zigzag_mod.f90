module zigzag_mod

 implicit none
 save
 
 double precision, allocatable :: zzpairden_t(:,:)
 double precision, allocatable :: zzpairdenij_t(:,:)
 double precision, allocatable :: zzcorrij(:), zzcorr(:)
 double precision, allocatable :: zzcorr1(:), zzcorr2(:)
 double precision, allocatable :: yycorrij(:), yycorr(:)
 double precision, allocatable :: yycorr1(:), yycorr2(:)
 double precision, allocatable :: znncorr(:), zn2ncorr(:)
 double precision, allocatable :: xold_sav(:,:), xnew_sav(:,:), zzposold(:,:), zzposnew(:,:)
 integer, allocatable          :: iold_indices(:), inew_indices(:)
 double precision, allocatable :: zzcum(:), zzsum(:), zzcm2(:)
 double precision              :: zzdelyr
 integer                       :: izigzag
 integer, parameter            :: nzzvars = 24
end module zigzag_mod
