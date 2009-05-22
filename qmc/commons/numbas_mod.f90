module numbas_mod

 use constants_mod
 implicit none
 save

 double precision exp_h_bas(MCTYPE),r0_bas(MCTYPE),rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
 double precision d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
 integer numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE)
 integer iwrwf(MBASIS_CTYPE,MCTYPE)

end module numbas_mod
