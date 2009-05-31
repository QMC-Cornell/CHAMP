module tempor_test_mod       

 use constants_mod
 implicit none
 save

! igvec_dft needs to be dimensioned with NGVEC_BIGX since in file orbitals_pw_pwscf the
! list of g vectors at the top is longer than what is actually used.
! The other arrays are dimensioned NGVEC2X rather than NGVECX because planewave code does not
! combine coefs. of G and -G, whereas QMC code does.
 integer               :: igvec_dft(3,NGVEC_BIGX),iwgvec(NGVEC2X)
 real(dp)              :: c_real(NGVEC2X),c_imag(NGVEC2X)
 integer               :: map_gvecdft_gvec(NGVEC2X),isign_gvecdft_gvec(NGVEC2X)
 real(dp), allocatable :: orb(:), dorb(:,:), ddorb(:)
 real(dp), allocatable :: orb_si(:),dorb_si(:,:),ddorb_si(:)
 integer, allocatable  :: iflag(:)
 real(dp)              :: rnorm,r(3),rkvec_tmp(3),rkvec_tmp2(3)
 integer               :: ngg(MKPTS),ngvec_dft



end module tempor_test_mod       
