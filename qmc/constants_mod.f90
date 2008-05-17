module constants_mod

  use types_mod

! string lengths
  integer, parameter    ::  max_string_len      = 100
  integer, parameter    ::  max_string_len_obj  = 50  ! important for speed!
  integer, parameter    ::  max_string_len_type = 15
  integer, parameter    ::  max_string_len_rout = 60
  integer, parameter    ::  max_string_len_file = 200

! max length of string arrays
  integer, parameter    ::  max_string_array_len = 1000 ! added for pathscale compiler
!added by WAS for pathscale compiler
  integer, parameter    ::  max_int_array_len = 1000 ! added for pathscale compiler
  integer, parameter    ::  max_double_array_len = 1000 ! added for pathscale compiler
!

! numbers
  real(dp), parameter   :: zero    = 0.d0
  real(dp), parameter   :: third   = 1.d0/3.d0
  real(dp), parameter   :: half    = 0.5d0
  real(dp), parameter   :: one     = 1.d0
  real(dp), parameter   :: two     = 2.d0
  real(dp), parameter   :: three   = 3.d0
  real(dp), parameter   :: four    = 4.d0
  real(dp), parameter   :: ten     = 10.d0
  real(dp), parameter   :: tenth   = 0.1d0

! mathematical constants
  real(dp), parameter   :: pi1        = 3.141592653589793d0
  real(dp), parameter   :: sqrt_pi    = 1.772453850905516d0
  real(dp), parameter   :: oneover2pi = 0.15915494309189535d0
  real(dp), parameter   :: oneover4pi = 0.07957747154594767d0

! conversion of units
  real(dp), parameter   :: angstrom_to_bohr = 1.8897261249935897d0

! maximun numbers of nodes
  integer, parameter    ::  max_nodes_nb = 500

! maximun numbers of routines
  integer, parameter    ::  max_routines_nb = 500

! maximun numbers of objects
  integer, parameter    ::  max_objects_nb = 1020

! dimension parameters
  include 'vmc/vmc.h'
  include 'vmc/pseudo.h'
  include 'vmc/ewald.h'
  include 'vmc/force.h'
  include 'fit/fit.h'
  include 'dmc/dmc.h'
  include 'vmc/numbas.h'
  include 'vmc/numorb.h'
  include 'vmc/MPI/mpi_qmc.h'

# if defined (MPI)
   include 'mpif.h'
# endif

end module constants_mod
