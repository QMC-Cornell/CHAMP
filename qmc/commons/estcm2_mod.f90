module estcm2_mod

 use constants_mod
 implicit none
 save

 double precision, allocatable :: wgcm2(:),wgcm21(:),egcm2(:),egcm21(:),pecm2(:),peicm2(:),tpbcm2(:),tjfcm2(:)
 double precision wcm2,wfcm2,wdcm2,wgdcm2, wcm21
 double precision wfcm21,wdcm21,ecm2,efcm2,ecm21
 double precision efcm21,ei1cm2,ei2cm2,ei3cm2
 double precision r1cm2,r2cm2,r3cm2,r4cm2,ricm2

end module estcm2_mod
