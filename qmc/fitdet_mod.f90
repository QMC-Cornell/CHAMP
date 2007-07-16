module fitdet_mod

  complex*16,       allocatable ::  cvd_sav(:,:,:)
  double precision, allocatable ::  vd_sav(:,:,:)
  double precision, allocatable ::  psid_sav(:)
  double precision, allocatable ::  d2d_sav(:)
  double precision, allocatable ::  div_vd_sav(:,:)
  complex*16,       allocatable ::  cvk_sav(:,:,:)
  double precision, allocatable ::  psik_sav(:)
  double precision, allocatable ::  div_vk_sav(:,:)
  double precision, allocatable ::  d2k_sav(:)
  integer                           iconfg
  integer                           isaved

end module fitdet_mod
