module control_mod

  use all_tools_mod
!  use walkers_mod

! threshold on statistical error on energy
  character(len=max_string_len)  :: title = 'no title'
  character(len=max_string_len)  :: seed = '1837465927472523'
! character(len=max_string_len) :: drift_type='unr93', ene_int='no_ene_int', rewt_type='sym', tmoves_type='det_balance1'
! character(len=max_string_len) :: drift_type='unr93', ene_int='unr93', rewt_type='sym', tmoves_type='det_balance1'
  character(len=max_string_len) :: drift_type='unr93', ene_int='ene_int_v9', rewt_type='sym', tmoves_type='det_balance1'

  real(dp)                :: error_threshold = 1.d30
  real(dp)                :: limit_rewt_dmc=10.d0, c_rewt=3.5d0, p_rewt=1.d0, adrift=0.5d0

  logical                 :: l_nstep_all_cpus = .true.
  logical                 :: l_hopping_moves  = .false.
  real(dp)                :: proba_hopping_moves = 0.d0

! branching in DMC?
  logical                 :: l_branching = .true.
  real(dp)                :: wt_lambda = 1.d0
  logical                 :: l_population_control = .true.
  logical                 :: l_reset_etrial = .true.

  contains

!===========================================================================
  subroutine control_menu
!---------------------------------------------------------------------------
! Description : menu for Monte Carlo parameters
!
! Created     : J. Toulouse, 24 Oct 2005
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere= 'control_menu'
  integer                        :: irn(4)

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of control menu --------------------------------------------------------------------------------'

! initialization
  MWALK=-1

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
   write(6,'(a)') ' nconf = [integer]: target walker population (default=100)'
   write(6,'(a)') ' mwalk = [integer]: max number walkers (default=nconf+50)'
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
   write(6,'(a)') ' branching = [logical] do branching in DMC? (default=true)'
   write(6,'(a)') ' wt_lambda = [real] Power to which DMC wts are raised unit physical time later. (default=1.d0)'
   write(6,'(a)') ' population_control = [logical] do population control in DMC? (default=true)'
   write(6,'(a)') ' ipq = [integer]'
   write(6,'(a)') ' itau_eff = [integer]'
   write(6,'(a)') ' itau_integ = [integer]'
   write(6,'(a)') ' iacc_rej = [integer]'
   write(6,'(a)') ' icross = [integer]'
   write(6,'(a)') ' icuspg = [integer]'
   write(6,'(a)') ' idiv_v = [integer]'
   write(6,'(a)') ' icut_br = [integer]'
   write(6,'(a)') ' icut_e = [integer]'
   write(6,'(a)') ' nfprod = [integer] number of products undone for population-control bias correction'
   write(6,'(a)') ' tau = [real] time step for DMC (default: 0.01)'
   write(6,'(a)') ' adrift = [real] parameter for calculating average velocity over time step for DMC (default: 0.5)'
   write(6,'(a)') ' drift_type = [string] Either unr93 or quadratic (default: unr93)'
   write(6,'(a)') ' ene_int = [string] : Type of energy integration (if any) used for local energy over time step tau'
   write(6,'(a)') ' limit_rewt_dmc = [real] : The reweight factor is not allowed to be larger than 1+limit_rewt_dmc*energy_sigma*tau'
   write(6,'(a)') ' c_rewt = [real] : The reweight exponent is divided by a factor that depends on c_rewt'
   write(6,'(a)') ' p_rewt = [real] : The reweight exponent is divided by a factor that depends on p_rewt'
   write(6,'(a)') ' rewt_type = [string] : Either "sym" or "pq"'
   write(6,'(a)') ' tmoves = [logical] do version 1 of size consistent tmoves of Casula, Moroni, Sorella, Filippi, JCP10'
   write(6,'(a)') ' tmoves_type [string] : Either casula1 or det_balance1'
   write(6,'(a)') ' vmc ... end : control menu for vmc'
   write(6,'(a)') ' error_threshold = [real] : montecarlo run until statistical error on energy reaches error_threshold'
   write(6,'(a)') ' nstep_total = [real]: For MPI, total number of steps per block for all CPUs'
   write(6,'(a)') ' reset_etrial = [logical] reset etrial after each set of equilibration blocks'
   write(6,'(a)') ' improved_gf = [logical] For all electron calculations, use improved Greens function which approximates the pair Greens function.'
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
  case ('mwalk');  call get_next_value (MWALK)
  case ('idump');  call get_next_value (idump)
  case ('irstar'); call get_next_value (irstar)
  case ('isite');  call get_next_value (isite)
  case ('ipr');    call get_next_value (ipr)
  case ('imetro'); call get_next_value (imetro)
  case ('delta');  call get_next_value (delta)
  case ('deltar'); call get_next_value (deltar)
  case ('deltat'); call get_next_value (deltat)
  case ('fbias');  call get_next_value (fbias)
  case ('idmc');   call get_next_value (idmc)
  case ('branching'); call get_next_value (l_branching)
  case ('wt_lambda'); call get_next_value (wt_lambda)
  case ('population_control'); call get_next_value (l_population_control)
  case ('reset_etrial'); call get_next_value (l_reset_etrial)
  case ('improved_gf'); call get_next_value (l_improved_gf)
  case ('ipq');    call get_next_value (ipq)
  case ('itau_eff');  call get_next_value (itau_eff)
  case ('itau_integ');  call get_next_value (itau_integ)
  case ('iacc_rej');  call get_next_value (iacc_rej)
  case ('icross');  call get_next_value (icross)
  case ('icuspg');  call get_next_value (icuspg)
  case ('idiv_v');  call get_next_value (idiv_v)
  case ('icut_br'); call get_next_value (icut_br)
  case ('icut_e');  call get_next_value (icut_e)
  case ('nfprod');  call get_next_value (nfprod); call object_modified ('nfprod')
  case ('tau');     call get_next_value (tau)
  case ('adrift');     call get_next_value (adrift)
   call require (lhere, 'adrift > 0', adrift > 0)
  case ('drift_type');     call get_next_value (drift_type)
   call require (lhere, 'drift_type=unr93 or quadratic', drift_type=='unr93' .or. drift_type=='quadratic')
  case ('ene_int') ; call get_next_value (ene_int)
!  call require (lhere, 'ene_int=unr93 or new_ene_int new_ene_int2 or no_ene_int or alfe', ene_int=='unr93' .or. ene_int=='new_ene_int' .or. ene_int=='new_ene_int2' .or. ene_int=='new_ene_int3' .or. ene_int=='new_ene_int4' .or. ene_int=='new_ene_int5' .or. ene_int=='new_ene_int7' .or. ene_int=='new_ene_int8' .or. ene_int=='new_ene_int9' .or. ene_int=='new_ene_int10' .or. ene_int=='no_ene_int' .or. ene_int=='alfe')
  case ('limit_rewt_dmc') ; call get_next_value (limit_rewt_dmc)
   call require (lhere, 'limit_rewt_dmc > 0', limit_rewt_dmc > 0)
  case ('c_rewt') ; call get_next_value (c_rewt)
   call require (lhere, 'c_rewt > 0', c_rewt > 0)
  case ('p_rewt') ; call get_next_value (p_rewt)
   call require (lhere, 'p_rewt > 0', p_rewt > 0)
  case ('rewt_type') ; call get_next_value (rewt_type)
   call require (lhere, 'rewt_type = sym or pq', rewt_type=='sym' .or. rewt_type=='pq')
  case ('tmoves');  call get_next_value (tmoves)
  case ('tmoves_type') ; call get_next_value (tmoves_type)
   call require (lhere, 'tmoves_type = casula1 or det_balance1', tmoves_type=='casula1' .or. tmoves_type=='det_balance1')

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
          write(6,'(a)') ' Printing out smoothing B-spline orbitals'
      case(8)
          write(6,'(a)') ' Printing out interpolating B-spline orbitals'
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

! Duplicated these here because there is a "if (use_parser)" before the later print statements
  if (index(mode,'dmc').ne.0) then
    write(6,'(/,a)') ' DMC algorithm parameters:'
    write(6,'(''adrift= '',f6.2)') adrift
    write(6,'(''drift_type= '',a)') trim(drift_type)
    write(6,'(''ene_int= '',a)') trim(ene_int)
    if(index(ene_int,'ene_int_v').ne.0 .or. index(ene_int,'ene_int_e').ne.0) then
      write(6,'(''c_rewt='',f6.3)') c_rewt
    endif
    if(ene_int=='ene_int_v' .or. ene_int=='ene_int_e' .or. ene_int=='ene_int_v8' .or. ene_int=='ene_int_e8') then
      write(6,'(''p_rewt='',f6.3)') p_rewt
    endif
    write(6,'(''itau_integ='',i3)') itau_integ
    write(6,'(''limit_rewt_dmc='',es9.1)') limit_rewt_dmc
    write(6,'(''rewt_type='',a,/)') trim(rewt_type)
    write(6,'(''tmoves='',l)') tmoves
    write(6,'(''tmoves_type='',a,/)') trim(tmoves_type)
  endif

! default value of MWALK
  if (MWALK == -1) then
    MWALK = nconf + 50
  endif
  if (MWALK < nconf) then
    call die (lhere, 'MWALK='+MWALK+' < nconf='+nconf)
  endif

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
! JT: Metropolis parameters need to be set also for DMC
! in the case where the starting walkers are generated by a initial VMC run
!  if (index(mode,'vmc').ne.0) then
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
!  endif

! DMC parameters
  if (index(mode,'dmc').ne.0) then
    write(6,*)
    write(6,'(a)')   ' DMC algorithm parameters:'
    write(6,'(a,9i4)') ' idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e =', idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
    write(6,'(a,9i4)') ' itau_integ =', itau_integ
    if(idmc.lt.0) write(6,'(a)') ' running DMC program in VMC mode'
    if(iabs(idmc).lt.1 .or. iabs(idmc).gt.3) stop 'iabs(idmc) must be 1 or 2 or 2'
    rttau=dsqrt(tau)
    rtrttau=dsqrt(rttau)
    write(6,'(a,i5)')    ' nfprod = ', nfprod
    write(6,'(a,f10.5)') ' time-step tau = ', tau
    if (.not. l_branching) then
      write(6,'(a)') ' Warning: branching turned off in DMC'
    endif
    if (wt_lambda /= 1.d0) then
      write(6,'(a,f10.5)') ' Warning: the weights of each previous generation will be raised to the power wt_lambda=',wt_lambda
    endif
    if (.not. l_population_control) then
      write(6,'(a)') ' Warning: population control turned off in DMC'
    endif
    write(6,'(''DMC ene_int= '',a)') ene_int
    if(index(ene_int,'ene_int_v').ne.0 .or. index(ene_int,'ene_int_e').ne.0) then
      write(6,'(''c_rewt='',f6.3)') c_rewt
    endif
    if(ene_int=='ene_int_v' .or. ene_int=='ene_int_e' .or. ene_int=='ene_int_v8' .or. ene_int=='ene_int_e8') then
      write(6,'(''p_rewt='',f6.3)') p_rewt
    endif
    write(6,'(''DMC limit_rewt_dmc='',es9.1)') limit_rewt_dmc
    write(6,'(''tmoves='',l)') tmoves
    write(6,'(''tmoves_type='',a)') tmoves_type
    write(6,'(''reset_etrial='',l)') l_reset_etrial
    write(6,'(''improved_gf='',l)') l_improved_gf
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
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere= 'vmc_menu'

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
