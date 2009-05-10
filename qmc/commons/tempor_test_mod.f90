module tempor_test_mod       

 use constants_mod

 integer               :: igvec_dft(3,NGVEC_BIGX),iwgvec(NGVEC2X)
 real(dp)              :: c_real(NGVEC2X),c_imag(NGVEC2X)
 integer               :: map_gvecdft_gvec(NGVEC2X),isign_gvecdft_gvec(NGVEC2X)
 real(dp), allocatable :: orb(:), dorb(:,:), ddorb(:)
 real(dp)              :: orb_si(MORB),dorb_si(3,MORB),ddorb_si(MORB)
 integer               :: iflag(MORB)
 real(dp)              :: rnorm,r(3),rkvec_tmp(3),rkvec_tmp2(3)
 integer               :: ngg(MKPTS),ngvec_dft



end module tempor_test_mod       
