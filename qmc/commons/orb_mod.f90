module orb_mod       

 use constants_mod

! real(dp), allocatable :: orb(:,:), dorb(:,:,:), ddorb(:,:)
 real(dp), pointer :: orb(:,:), dorb(:,:,:), ddorb(:,:)

end module orb_mod       
