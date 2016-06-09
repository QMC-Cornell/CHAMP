module variables_mod

  use constants_mod

! mode
  integer, parameter              :: modes_nb = 11
! Variable string lengths work with Pathscale 3 but not Pathscale 2.5 in Lyngby, so temporarily make all string lengths equal
! character(len=max_string_len)   :: modes (modes_nb) = (/'fit', 'fit_mpi',         &
!                                    'vmc', 'vmc_mov1', 'vmc_mpi', 'vmc_mov1_mpi',  &
!                                    'dmc', 'dmc_mov1', 'dmc_mov1_mpi1', 'dmc_mov1_mpi2', 'dmc_mov1_mpi3'/)
  character(len=max_string_len)   :: modes (modes_nb) = (/'fit          ', 'fit_mpi      ',              &
                                     'vmc          ', 'vmc_mov1     ', 'vmc_mpi      ', 'vmc_mov1_mpi ', &
                                     'dmc          ', 'dmc_mov1     ', 'dmc_mov1_mpi1', 'dmc_mov1_mpi2', 'dmc_mov1_mpi3'/)
  logical                         :: l_mode_mpi            = .false.
  logical                         :: l_mode_fit            = .false.
  logical                         :: l_mode_fit_mpi        = .false.
  logical                         :: l_mode_vmc            = .false.
  logical                         :: l_mode_vmc_mov1       = .false.
  logical                         :: l_mode_vmc_mov1_mpi   = .false.
  logical                         :: l_mode_dmc            = .false.
  logical                         :: l_mode_dmc_mov1       = .false.
  logical                         :: l_mode_dmc_mov1_mpi1  = .false.
  logical                         :: l_mode_dmc_mov1_mpi2  = .false.
  logical                         :: l_mode_dmc_mov1_mpi3  = .false.

! general common variables
  logical                         :: use_parser = .false.
  integer                         :: nstep_total
  integer                         :: walk_nb = 1
  logical                         :: equilibration  = .false.
  logical                         :: l_end_of_block = .false.
  integer                         :: spin_nb = 2
  logical                         :: debug = .false.
  logical                         :: l_warning = .false.

! label of last occupied orbital
  integer                         :: orb_occ_last_in_wf_lab
  logical                         :: l_print_orbital_excitations = .false.

! CPU time
  real(dp)                :: cpu_time_start
  real(dp)                :: cpu_time_last

! memory used
  real(dp)                :: memory = 0.d0

! current routine name
  character(len=max_string_len)   :: here       = 'undefined'

! Has a MC run already been done?
  logical                         :: run_done   = .false.

! number of parameters to optimized
  integer                         :: param_nb = 0

! diagonal stabilizer for optimization
  real(dp)               :: diag_stab = 0.d0

! wave function optimization
  logical                         :: l_opt      = .false.

  integer                         :: iter_global_optimization
  logical                         :: l_lastrun_global_optimization
  logical                         :: l_opt_nwt = .false.
  logical                         :: l_opt_lin  = .false.
  logical                         :: l_opt_ptb = .false.
  logical                         :: l_opt_reb = .false.
  logical                         :: l_opt_ovlp_fn = .false.
  logical                         :: l_opt_grad = .false.
  logical                         :: l_diagonal_overlap = .false.

  logical                         :: l_opt_orb  = .false.
  logical                         :: l_opt_csf  = .false.
  logical                         :: l_opt_jas  = .false.
  logical                         :: l_opt_exp  = .false.
  logical                         :: l_opt_geo  = .false.
  real(dp)                        :: increase_blocks_limit 

!! pjas WAS
  logical                         :: l_opt_pjasen  = .false.
  logical                         :: l_opt_pjasee  = .false.
  logical                         :: l_opt_pjas  = .false.

  logical                         :: l_opt_jas_2nd_deriv = .false.
  logical                         :: l_opt_orb_eig  = .false.
  logical                         :: l_opt_orb_energy  = .false.
  logical                         :: l_stab  = .true.
  character(len=max_string_len)   :: stabilization = 'identity'
  logical                         :: l_casscf = .false.
  logical                         :: l_check_redundant_orbital_derivative = .true.
  logical                         :: l_active_orb_ortho_constraint = .true.

  logical                         :: l_deriv2nd = .true.
  logical                         :: l_deriv_bound = .false.
  real(dp)                        :: deriv_bound_value = 10.d0
  logical                         :: l_eloc_bound = .false.
  real(dp)                        :: eloc_bound_value = 10.d0

  real(dp)                        :: lambda = 0.3d0

! total number of configurations (walkers)
  integer   :: nconf_total = 1

! weight of walkers
  real(dp), allocatable   :: walker_weights (:)

! probability of acceptance and rejection
  real(dp)    :: prob_acc
  real(dp)    :: prob_rej

! compute only occupied orbitals?
  integer :: iorb_used = 1

! print out orbitals for periodic systems?
  logical :: l_print_orbitals_pw = .false.

  integer :: iter_opt_max_nb = 100

! number of Jastrow and CSF parameters given in input
!JT  integer :: nparmj_input = 0
!JT  integer :: nparmcsf_input = 0

! forces
  integer                                     :: forces_nb = 0
  integer, allocatable                        :: forces_cent (:)
  integer, allocatable                        :: forces_direct (:)

! function
  real(dp), external     :: erf_spline
  real(dp), external     :: ei2

  logical                          :: print_radial_probability = .false.

  logical                          :: l_equilibration = .false.

!  real(dp)                       :: eloc_av_err ! temporary

end module variables_mod
