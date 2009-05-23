module distance_mod

 use constants_mod
 implicit none
 save

 double precision rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

end module distance_mod
