module jaspar6_mod

 use constants_mod
 implicit none
 save

 double precision asymp_jasa(MCTYPE,MWF),asymp_jasb(2,MWF),dasymp_jasa(MCTYPE,MWF),dasymp_jasb(2,MWF),    &
 d2asymp_jasa(MCTYPE,MWF),d2asymp_jasb(2,MWF),asymp_r_en(MWF),dasymp_r_en(MWF),d2asymp_r_en(MWF),         &
 asymp_r_ee(MWF),dasymp_r_ee(MWF),d2asymp_r_ee(MWF),cutjas_en,cutjasi_en,c1_jas6_en(MWF),c2_jas6_en(MWF), &
 cutjas_ee,cutjasi_ee,c1_jas6_ee(MWF),c2_jas6_ee(MWF)

end module jaspar6_mod
