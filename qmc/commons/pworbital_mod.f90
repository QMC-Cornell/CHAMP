module pworbital_mod

 use constants_mod
 implicit none
 save

 double precision c_rp(NGVECX,MORB_OCC),c_rm(NGVECX,MORB_OCC),c_ip(NGVECX,MORB_OCC),c_im(NGVECX,MORB_OCC)
 integer ngorb(MORB),isortg(NGVECX,MORB),isortk(MKPTS),icmplx

end module pworbital_mod
