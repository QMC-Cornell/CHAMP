module branch_mod

 use constants_mod
 implicit none
 save

 double precision wtgen(0:MFPRD1),ff(0:MFPRD1),eoldw(MWALK,MFORCE)
 double precision pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE)
 double precision wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod
 integer nwalk

end module branch_mod
