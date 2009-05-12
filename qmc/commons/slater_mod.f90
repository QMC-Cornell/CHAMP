module slater_mod      

 use constants_mod

 real(dp), allocatable :: slmui(:,:), slmdi(:,:)
 real(dp), allocatable :: fpu(:,:,:), fpd(:,:,:)
 real(dp), allocatable :: fppu(:,:), fppd(:,:)
 real(dp), allocatable :: detu(:), detd(:)
 real(dp), allocatable :: ddeti_deti(:,:,:), d2edeti_deti(:,:)
 real(dp), allocatable :: deti_det(:), ddeti_det(:,:,:)
 real(dp), allocatable :: d2deti_det(:), detij_det(:,:)
 real(dp)              :: d2det_det

end module slater_mod      
