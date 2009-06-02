module stepv_mod

 use constants_mod
 implicit none
 save

 integer, parameter :: NRAD = 1001
 double precision, parameter :: RADMAX = 8.d0
 double precision :: delri =(NRAD-1)/RADMAX
 double precision, allocatable :: try(:),suc(:),trunfb(:),rprob(:),ekin(:),ekin2(:)

end module stepv_mod
