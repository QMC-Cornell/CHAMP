module stats_index_mod
  implicit none

  integer, parameter :: elocal0_   =1
  integer, parameter :: elocalf_   =2
  integer, parameter :: elocalg_   =3
  integer, parameter :: ekinpb_    =4
  integer, parameter :: ekinjf_    =5
  integer, parameter :: epot_      =6
  integer, parameter :: epot_ee_   =7
  integer, parameter :: taueff_    =8
  integer, parameter :: R1_        =9
  integer, parameter :: R2_        =10
  integer, parameter :: R3_        =11
  integer, parameter :: R4_        =12
  integer, parameter :: RI_        =13
  integer, parameter :: ws0_       =14
  integer, parameter :: wsf_       =15
  integer, parameter :: wsg_       =16
  integer, parameter :: wtsq_      =17
  integer, parameter :: nest_static=17
  integer            :: offset_enefrag_
  integer            :: offset_taufrag_
end module stats_index_mod
