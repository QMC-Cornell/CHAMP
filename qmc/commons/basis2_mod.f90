module basis2_mod

 use constants_mod
 implicit none
 save

 double precision zex2(MRWF,MCTYPE,MWF)
 integer n_bas(MBASIS),l_bas(MBASIS),m_bas(MBASIS)
 integer icenter_basis(MBASIS),ictype_basis(MBASIS)
 integer nbasis_ctype(MCTYPE),n_bas2(MRWF,MCTYPE),iwrwf2(MBASIS)

end module basis2_mod
