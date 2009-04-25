module initialization_mod

  use all_tools_mod
  use orbitals_mod

  contains

! ===================================================================================
  subroutine initialization_before_parser
! -----------------------------------------------------------------------------------
! Description   : initialization of some global variables before parser
!
! Created       : J. Toulouse, 11 Mar 2009
! -----------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

!  write(6,*)
!  write(6,'(a)') 'Beginning of global initialization -----------------------------------------------------------------------'

! constant
  pi=four*datan(one)

! default values for input variables
  iperiodic=0
  ibasis=1
  hb = 0.5d0
  etrial = 0.d0
  nstep=1000
  call object_modified ('nstep')
  nblk=10
  nblkeq=1
  nconf=100
  call object_modified ('nconf')
  nstep_total = nstep * nproc
  nconf_total = nconf * nproc
  nconf_new=0
  idump=0
  irstar=0
  isite=1
  ipr=-2
  imetro=6
  delta=1.
  deltar=5.
  deltat=1.
  fbias=1.
  idmc=2
  ipq=1
  itau_eff=1
  iacc_rej=1
  icross=1
  icuspg=0
  idiv_v=0
  icut_br=0
  icut_e=0
  nfprod=50
  tau=0.01
  nloc=0
  call object_modified ('nloc')
  numr=-3
  call object_modified ('numr')
  nforce=1
  call object_modified ('nforce')
  nefp=0
  nquad=6
  ndim=3
  call object_modified ('ndim')
  inum_orb=0
  iorb_used=0
  iorb_format='unused'

! csfs
  ncsf=1
  csf_coef(1:ncsf,1) = 1.d0
  ndet_in_csf (1:ncsf) = 1
  iwdet_in_csf(1,1:ncsf) = 1
  cdet_in_csf(1,1:ncsf) = 1.d0
  call object_modified ('ncsf')
  call object_modified ('csf_coef')
  call object_modified ('ndet_in_csf')
  call object_modified ('iwdet_in_csf')
  call object_modified ('cdet_in_csf')

! Jastrow variables
  ianalyt_lap=1
  ijas=4
  isc=2
  call object_modified ('isc')
  nspin1=1
  nspin2=1
  nord=5
  ifock=0
  scalek(1)=0.5d0
  call object_modified ('scalek')
  a21=0.d0
  norda=5
  nordb=5
  nordc=5
  call object_modified ('norda')
  call object_modified ('nordb')
  call object_modified ('nordc')

! correlated sampling index
  iwf = 1
  nwftype=1
  iwftype(1)=1
  call object_modified ('iwf')
 
! optimization
  nparm=0
  nblk_max=10000
  increase_blocks_limit=nblk_max
  add_diag(1) = 1.d-8
  diag_stab = add_diag(1)
  p_var=0.d0
  call object_modified ('p_var')
!  ndata=1000
  icusp=1
  icusp2=1
  nsig=5
  ncalls=1000
  iopt=21101
  ipr_opt=1
  i3body=0
  irewgt=0
  iaver=0
  istrech=0
  ipos=0
!  idcds=0
!  idcdr=0
!  idcdt=0
!  id2cds=0
!  id2cdr=0
!  id2cdt=0
!  idbds=0
!  idbdr=0
!  idbdt=0
  nparml=0
  nparmb=5
  nparmc=15
  nparmf=0
  nparmcsf=0
  nparms=0
  nparmg=0

! initialize saved configuration indice iconfg (necessary for some compilers)
!  isaved=0
!  iconfg=1

!  write(6,'(a)') 'End of global initialization -----------------------------------------------------------------------------'

 end subroutine initialization_before_parser

! ===================================================================================
  subroutine initialization
! -----------------------------------------------------------------------------------
! Description   : initialization of some global variables after reading Cyrus's input
!
! Created       : J. Toulouse, 03 Oct 2006
! -----------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

  write(6,*)
  write(6,'(a)') 'Beginning of global initialization -----------------------------------------------------------------------'

! set nstep_total
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
  write(6,'(a,i8)') ' Number of computed orbitals initialized to ', norb
  call object_modified ('norb')

! orbital coefficients on normalized and orthonormalized basis functions
!  if(.not.(ibasis.ge.3 .and. ibasis.le.6)) then ! do not do it for quantum dots, rings, etc...
   if (inum_orb == 0) then ! only if non-numerical orbitals
     call coef_orb_on_norm_basis_from_coef (1)
   endif
!  if (trim(basis_functions_varied) == 'orthonormalized') then
!   call coef_orb_on_ortho_basis_from_coef
!  endif
!  endif

! number of walkers
  if (index(mode,'vmc') /= 0) then
   nwalk = 1
  elseif (index(mode,'dmc') /= 0) then
   nwalk = nconf
  endif
  call object_modified ('nwalk')
  write(6,*)
  write(6,'(a,i8)') ' Number of walkers initialized to ', nwalk

! maximal number of optimization iterations
  iter_opt_max_nb = nopt_iter

! save number of Jastrow and CSF parameters in the input
  nparmj_input = nparmj
  nparmcsf_input = nparmcsf

! correlated sampling
  iwf = 1
  nwftype=1
  iwftype(1)=1

  write(6,'(a)') 'End of global initialization -----------------------------------------------------------------------------'

 end subroutine initialization

! ===================================================================================
  subroutine set_mode
! -----------------------------------------------------------------------------------
! Description   : set logical variables for mode
!
! Created       : J. Toulouse, 13 Jun 2008
! -----------------------------------------------------------------------------------
  implicit none
  include 'commons.h'


  l_mode_mpi = .false.
  l_mode_fit = .false.
  l_mode_fit_mpi = .false.
  l_mode_vmc = .false.
  l_mode_vmc_mov1 = .false.
  l_mode_vmc_mov1_mpi = .false.
  l_mode_dmc = .false.
  l_mode_dmc_mov1 = .false.
  l_mode_dmc_mov1_mpi1 = .false.
  l_mode_dmc_mov1_mpi2 = .false.
  l_mode_dmc_mov1_mpi3 = .false.

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
  if (index(mode, 'dmc_mov1_mpi3') /= 0) l_mode_dmc_mov1_mpi3 = .true.

 end subroutine set_mode

end module initialization_mod
