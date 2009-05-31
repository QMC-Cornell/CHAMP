module orbital_grid_mod

  use types_mod
! use all_tools_mod

  save

! Declaration of shared variables and default values
  real(dp), allocatable ::  orb_num(:,:,:,:),dorb_num(:,:,:,:,:),ddorb_num(:,:,:,:)
  integer ::                ngrid_orbx,ngrid_orby,ngrid_orbz
  integer ::                igrad_lap

end module orbital_grid_mod

!     common /orbital_per_num/ orb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
!    &,dorb_num(3,MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
!    &,ddorb_num(MORB_OCC,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1,0:MGRID_ORB_PER-1)
!    &,ngrid_orbx,ngrid_orby,ngrid_orbz
