module config_mod

 use constants_mod

 real(dp), allocatable :: xold(:,:), xnew(:,:)
 real(dp), allocatable :: vold(:,:), vnew(:,:)
 real(dp), allocatable :: psi2o(:), psi2n(:)
 real(dp), allocatable :: eold(:), enew(:)
 real(dp)              :: peo,pen,peio,pein,tjfn,tjfo,psido,psijo
 real(dp), allocatable :: rmino(:), rminn(:)
 real(dp), allocatable :: rvmino(:,:), rvminn(:,:)
 real(dp), allocatable :: rminon(:), rminno(:)
 real(dp), allocatable :: rvminon(:,:), rvminno(:,:)
 real(dp), allocatable :: delttn(:)
 integer,  allocatable :: nearesto(:), nearestn(:)

end module config_mod
