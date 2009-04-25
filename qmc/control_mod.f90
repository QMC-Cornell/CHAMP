module control_mod

  use all_tools_mod
  use walkers_mod

! threshold on statistical error on energy
  character(len=max_string_len)  :: title = 'no title'
  character(len=max_string_len)  :: seed = '1837465927472523'

  real(dp)                :: error_threshold = 1.d30

  logical                 :: l_nstep_all_cpus = .true.
  logical                 :: l_hopping_moves  = .false.
  real(dp)                :: proba_hopping_moves = 0.d0

  contains

!===========================================================================
  subroutine control_menu
!---------------------------------------------------------------------------
! Description : menu for Monte Carlo parameters
!
! Created     : J. Toulouse, 24 Oct 2005
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere= 'control_menu'
  integer                        :: irn(4)

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of control menu --------------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for control menu:'
   write(6,'(a)') 'control'
   write(6,'(a)') ' title = [string]: title for the job'
   write(6,'(a)') ' seed = [integer]'
   write(6,'(a)') ' iperiodic = [integer]: 0 for finite system, >0 for periodic system (default: 0)'
   write(6,'(a)') ' etrial = [real]: trial energy'
   write(6,'(a)') ' nstep = [integer]: number of steps per block'
   write(6,'(a)') ' nblk = [integer]: number of blocks'
   write(6,'(a)') ' nblkeq = [integer]: number of equilibration blocks'
   write(6,'(a)') ' nconf = [integer]: target walker population'
   write(6,'(a)') ' nconf_new = [integer]: number of walkers to write out'
   write(6,'(a)') ' idump = [integer]: write dump file for restart (default: 0)'
   write(6,'(a)') ' irstar = [integer]: restart from dump file (default: 0)'
   write(6,'(a)') ' isite = [integer]: start from random configuration (default: 1)'
   write(6,'(a)') ' ipr = [integer]: print level (default: -2)'
   write(6,'(a)') ' imetro = [integer] type of Metropolis algorithm (default=6)'
   write(6,'(a)') ' delta = [real]'
   write(6,'(a)') ' deltar = [real]'
   write(6,'(a)') ' deltat = [real]'
   write(6,'(a)') ' fbias = [real]'
   write(6,'(a)') ' idmc = [integer] type of DMC algorithm (default=2)'
   write(6,'(a)') ' ipq = [integer]'
   write(6,'(a)') ' itau_eff = [integer]'
   write(6,'(a)') ' iacc_rej = [integer]'
   write(6,'(a)') ' icross = [integer]'
   write(6,'(a)') ' icuspg = [integer]'
   write(6,'(a)') ' idiv_v = [integer]'
   write(6,'(a)') ' icut_br = [integer]'
   write(6,'(a)') ' icut_e = [integer]'
   write(6,'(a)') ' nfprod = [integer] number of products undone for population-control bias correction'
   write(6,'(a)') ' tau = [real] time step for DMC (default: 0.01)'
   write(6,'(a)') ' nefp = [integer] (default: 0)'
   write(6,'(a)') ' vmc ... end : control menu for vmc'
   write(6,'(a)') ' error_threshold = [real] : montecarlo run until statistical error on energy reaches error_threshold'
   write(6,'(a)') ' nstep_total = [real]: For MPI, total number of steps per block for all CPUs'
   write(6,'(a)') 'end'
   write(6,*)

  case ('title')
   call get_next_value (title)
   write (6,'(2a)') ' job title = ',trim(title)

  case ('seed')
   call get_next_value (seed)

  case ('iperiodic'); call get_next_value (iperiodic)
  case ('etrial'); call get_next_value (etrial)
  case ('nstep');  call get_next_value (nstep); call object_modified ('nstep')
  case ('nblk');   call get_next_value (nblk)
  case ('nblkeq'); call get_next_value (nblkeq)
  case ('nconf');  call get_next_value (nconf); call object_modified ('nconf')
  case ('nconf_new'); call get_next_value (nconf_new)
  case ('idump');  call get_next_value (idump)
  case ('irstar'); call get_next_value (irstar)
  case ('isite');  call get_next_value (isite)
  case ('ipr');    call get_next_value (ipr)
  case ('imetro'); call get_next_value (imetro)
  case ('delta');  call get_next_value (delta)
  case ('deltar'); call get_next_value (deltar)
  case ('deltat'); call get_next_value (deltat)
  case ('fbias');  call get_next_value (fbias)
  case ('idmc');  call get_next_value (idmc)
  case ('ipq');  call get_next_value (ipq)
  case ('itau_eff');  call get_next_value (itau_eff)
  case ('iacc_rej');  call get_next_value (iacc_rej)
  case ('icross');  call get_next_value (icross)
  case ('icuspg');  call get_next_value (icuspg)
  case ('idiv_v');  call get_next_value (idiv_v)
  case ('icut_br');  call get_next_value (icut_br)
  case ('icut_e');  call get_next_value (icut_e)
  case ('nfprod');  call get_next_value (nfprod)
  case ('tau');  call get_next_value (tau)
  case ('nefp');  call get_next_value (nefp)

  case ('vmc')
   call vmc_menu

  case ('error_threshold')
   call get_next_value (error_threshold)
   write (6,'(a,es15.8)') ' setting target statistical error on the energy to ', error_threshold

  case ('nstep_total')
   call get_next_value (nstep_total)
   nstep = max(1,int(nstep_total/nproc))
   nstep_total = nstep * nproc
   write(6,'(a,i8)') ' nstep total =', nstep_total
   write(6,'(a,i8)') ' nstep per CPU =', nstep

  case ('print_orbitals_pw')
   call get_next_value (l_print_orbitals_pw)
   if (l_print_orbitals_pw) then
      select case (inum_orb)
      case(0)
          write(6,'(a)') ' Printing out plane wave orbitals'
      case(4)
          write(6,'(a)') ' Printing out Lagrange polynomial orbitals'
      case(5)
          write(6,'(a)') ' Printing out pp-spline orbitals'
      case(6)
          write(6,'(a)') ' Printing out B-spline orbitals'
      end select
      call print_orbitals_pw
   endif

!  case ('nstep_all_cpus')
!   call get_next_value (l_nstep_all_cpus)
!
!   if (.not. l_nstep_all_cpus) then
!    nstep = nstep_input
!    nstep_total = nstep * nproc
!    write(6,'(2a,i)') trim(lhere),': nstep total =', nstep_total
!    write(6,'(2a,i)') trim(lhere),': nstep per CPU =', nstep
!   endif

  case ('end')
   exit
  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines

! printing
  if (use_parser) then

  if (iperiodic > 0) then
     write(6,'(a)') ' type of system: periodic system' 
  else
     write(6,'(a)') ' type of system: finite system' 
  endif

  write(6,'(a,f12.6)') ' trial energy = ',etrial

! Random seed
  write(6,*)
  read(seed,'(4i4)') irn
  if (index(mode,'vmc').ne.0 .or. index(mode,'dmc').ne.0) then
    write (6,'(a,4i4)') ' random number seed = ',irn
    call setrn(irn)
  endif

! Basic QMC parameters
  if (index(mode,'vmc').ne.0 .or. index(mode,'dmc').ne.0) then
   if(irstar.eq.1) nblkeq=0
   write (6,'(a,i10)') ' number of steps per block              = ',nstep
   write (6,'(a,i10)') ' number of blocks after equilibration   = ',nblk
   write (6,'(a,i10)') ' number of blocks before equilibration  = ',nblkeq
   if (index(mode,'dmc').ne.0) then
     write(6,'(a,i5)') ' target walker population per processor = ',nconf
     call require (lhere, 'nconf >0', nconf >0)
   endif
   write(6,'(a,i5)')   ' number of configurations saved         = ',nconf_new
   if(irstar.eq.1) write(6,'(a)') ' job is starting from restart file'
   if(idump.eq.1) write(6,'(a)') ' job will write restart file at end of run'
  endif

  if (mode.eq.'dmc' .or. mode.eq.'dmc_mov1' .or. mode.eq.'dmc_mov1_mpi1') then
    nconf_global=nconf
   elseif (mode.eq.'dmc_mov1_mpi2' .or. mode.eq.'dmc_mov1_mpi3') then
    nconf_global=nconf*nproc
  endif

! set nstep_total
  nstep_total = nstep * nproc

! total number of configurations
  nconf_total = nconf * nproc

! number of walkers
  if (index(mode,'vmc') /= 0) then
   nwalk = 1
  elseif (index(mode,'dmc') /= 0) then
   nwalk = nconf
  endif
  call object_modified ('nwalk')
  write(6,'(a,i8)') ' number of walkers initialized to ', nwalk

! Metropolis parameters
  if (index(mode,'vmc').ne.0) then
    write(6,*)
    deltai=one/delta
    if(deltar.lt.one) then
       write(6,*) '**Warning value of deltar reset to 2.'
       deltar=two
    endif
    if(deltat.lt.zero .or. deltat.gt.two) then
       write(6,*) '**Warning value of deltat reset to 2.'
       deltat=two
    endif
! Truncate fbias so that it is never negative, and the quantity sampled is never negative
    fbias=dmin1(two,dmax1(zero,fbias))
    write(6,'(a,i10)')   ' version of Metropolis  = ', imetro
    write(6,'(a,f10.5)') ' step size              = ', delta
    write(6,'(a,f10.5)') ' radial step multiplier = ', deltar
    write(6,'(a,f10.5)') ' cos(theta) step size   = ', deltat
    write(6,'(a,f10.5)') ' force bias             = ', fbias
    if(imetro.ne.1 .and. imetro.ne.6) stop 'imetro must be 1 or 6 (accel. Metropolis)'
    if(imetro.ne.1 .and. iperiodic.gt.0) stop 'In order to do VMC calculation for periodic system run dmc or dmc.mov1 with idmc < 0 or run vmc with imetro=1'
  endif

! DMC parameters
  if (index(mode,'dmc').ne.0) then
    write(6,*)
    write(6,'(a)')   ' DMC algorithm parameters:'
    write(6,'(a,9i4)') ' idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e =', idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
    if(idmc.lt.0) write(6,'(a)') ' running DMC program in VMC mode'
    if(iabs(idmc).ne.1 .and. iabs(idmc).ne.2) stop 'iabs(idmc) must be 1 or 2'
    rttau=dsqrt(tau)
    write(6,'(a,i5)')    ' nfprod = ', nfprod
    write(6,'(a,f10.5)') ' time-step tau = ', tau
  endif

  endif ! if use_parser


  write(6,'(a)') 'End of control menu --------------------------------------------------------------------------------------'

  end subroutine control_menu

!===========================================================================
  subroutine vmc_menu
!---------------------------------------------------------------------------
! Description : menu for vmc
!
! Created     : J. Toulouse, 12 Apr 2007
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere= 'vmc_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))
  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for vmc menu:'
   write(6,'(a)') 'vmc'
   write(6,'(a)') ' proba_hopping_moves = [real] : probability of hopping electron moves between the atoms (default = 0)'
   write(6,'(a)') 'end'
   write(6,*)

  case ('proba_hopping_moves')
   call get_next_value (proba_hopping_moves)
   call require (lhere, 'proba_hopping_moves >= 0', proba_hopping_moves >= 0) !fp
   call require (lhere, 'proba_hopping_moves <= 1', proba_hopping_moves <= 1) !fp
   if (proba_hopping_moves > 0) then
     l_hopping_moves = .true.
     write(6,'(a,es15.8)') 'Hopping moves between atoms will be used in VMC with probability = ',proba_hopping_moves
   endif

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines

  end subroutine vmc_menu

end module control_mod
