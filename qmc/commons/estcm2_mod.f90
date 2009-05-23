module estcm2_mod

 use constants_mod
 implicit none
 save

 double precision wcm2,wfcm2,wgcm2(MFORCE),wdcm2,wgdcm2, wcm21
 double precision wfcm21,wgcm21(MFORCE),wdcm21, ecm2,efcm2,egcm2(MFORCE), ecm21
 double precision efcm21,egcm21(MFORCE),ei1cm2,ei2cm2,ei3cm2, pecm2(MFORCE),peicm2(MFORCE),tpbcm2(MFORCE)
 double precision tjfcm2(MFORCE),r2cm2,ricm2

end module estcm2_mod
