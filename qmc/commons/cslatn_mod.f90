module cslatn_mod       

 use constants_mod

 real(dp)              :: cslmin(MMAT_DIM,MDETUD),cdetn(MDETUD)
 real(dp)              :: cddeti_detin(3,MELEC,MDETUD),cd2edeti_detin(MELEC,MDETUD)
 real(dp), allocatable :: cdorb(:,:),cddorb(:)

end module cslatn_mod       
