module fragments_mod
  use types_mod

  integer :: nfrag
  logical :: l_fragments
  integer, allocatable :: iwfragnucl(:), iwfragelec(:), nelecfrag(:)
  real(dp), allocatable :: enefrag(:), pecent_frag(:), etrialfrag(:)

end module fragments_mod
