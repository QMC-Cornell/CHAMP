module jaspar6_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: asymp_jasa(:,:),asymp_jasb(:,:),dasymp_jasa(:,:),dasymp_jasb(:,:)
 double precision, allocatable :: d2asymp_jasa(:,:),d2asymp_jasb(:,:),asymp_r_en(:),dasymp_r_en(:),d2asymp_r_en(:)
 double precision, allocatable :: asymp_r_ee(:),dasymp_r_ee(:),d2asymp_r_ee(:)
 double precision              :: cutjas_en,cutjasi_en
 double precision, allocatable :: c1_jas6_en(:),c2_jas6_en(:)
 double precision              :: cutjas_ee,cutjasi_ee
 double precision, allocatable :: c1_jas6_ee(:),c2_jas6_ee(:)

end module jaspar6_mod
