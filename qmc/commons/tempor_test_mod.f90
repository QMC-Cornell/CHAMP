module tempor_test_mod       

 use constants_mod
 implicit none
 save

! igvec_dft needs to be dimensioned with NGVEC_BIGX since in file orbitals_pw_pwscf the
! list of g vectors at the top is longer than what is actually used.
! The other arrays are dimensioned NGVEC2X rather than NGVECX because planewave code does not
! combine coefs. of G and -G, whereas QMC code does.
 integer, allocatable  :: igvec_dft(:,:),iwgvec(:)
 real(dp), allocatable :: c_real(:),c_imag(:)
 integer, allocatable  :: map_gvecdft_gvec(:),isign_gvecdft_gvec(:)
 real(dp), allocatable :: orb(:), dorb(:,:), ddorb(:)
 real(dp), allocatable :: orb_si(:),dorb_si(:,:),ddorb_si(:)
 integer, allocatable  :: iflag(:)
 real(dp)              :: rnorm,r(3),rkvec_tmp(3),rkvec_tmp2(3)
 integer               :: ngvec_dft
 integer, allocatable  :: ngg(:)



end module tempor_test_mod       
