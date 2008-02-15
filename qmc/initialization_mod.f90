module initialization_mod

  use all_tools_mod
  use orbitals_mod

  contains

! ===================================================================================
  subroutine initialization
! -----------------------------------------------------------------------------------
! Description   : initialization of some global variables after reading Cyrus's input
!
! Created       : J. Toulouse, 03 Ocr 2006
! -----------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'initialization'
  integer bas_i

  write(6,*)
  write(6,'(a)') 'Beginning of global initialization -----------------------------------------------------------------------'

! mode
  if (index(mode, 'mpi') /= 0)           l_mode_mpi = .true.
  if (index(mode, 'fit') /= 0)           l_mode_fit = .true.
  if (index(mode, 'fit_mpi') /= 0)       l_mode_fit_mpi = .true.
  if (index(mode, 'vmc') /= 0)           l_mode_vmc = .true.
  if (index(mode, 'vmc_mov1') /= 0)      l_mode_vmc_mov1 = .true.
  if (index(mode, 'vmc_mov1_mpi') /= 0)  l_mode_vmc_mov1_mpi = .true.
  if (index(mode, 'dmc') /= 0)           l_mode_dmc = .true.
  if (index(mode, 'dmc_mov1') /= 0)      l_mode_dmc_mov1 = .true.
  if (index(mode, 'dmc_mov1_mpi1') /= 0) l_mode_dmc_mov1_mpi1 = .true.
  if (index(mode, 'dmc_mov1_mpi2') /= 0) l_mode_dmc_mov1_mpi2 = .true.

! set nstep_total
  nstep_input = nstep
!  nstep = max(1,int(nstep/nproc))
  nstep_total = nstep * nproc

! total number of configurations
  nconf_total = nconf * nproc

! total number of orbitals
  call object_provide ('norb')
  orb_tot_nb = norb
  call object_modified ('orb_tot_nb')

! number of orbitals to compute
  call object_provide ('orb_occ_last_in_wf_lab')
  norb = orb_occ_last_in_wf_lab
  write(6,'(a,i)') ' Number of computed orbitals initialized to ', norb
  call object_modified ('norb')

! orbital coefficients on normalized and orthonormalized basis functions
  call coef_orb_on_norm_basis_from_coef (1)
!  if (trim(basis_functions_varied) == 'orthonormalized') then
!   call coef_orb_on_ortho_basis_from_coef
!  endif

! number of walkers
  if (index(mode,'vmc') /= 0) then
   nwalk = 1
  elseif (index(mode,'dmc') /= 0) then
   nwalk = nconf
  endif
  call object_modified ('nwalk')
  write(6,*) 
  write(6,'(a,i)') ' Number of walkers initialized to ', nwalk

! maximal number of optimization iterations
  iter_opt_max_nb = nopt_iter

! save number of Jastrow and CSF parameters in the input
  nparmj_input = nparmj
  nparmcsf_input = nparmcsf

  write(6,'(a)') 'End of global initialization -----------------------------------------------------------------------------'

 end subroutine initialization

end module initialization_mod
