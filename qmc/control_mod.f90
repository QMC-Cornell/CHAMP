module control_mod

  use all_tools_mod
  use walkers_mod

! threshold on statistical error on energy
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
   write(6,'(a)') ' vmc ... end : control menu for vmc'
   write(6,'(a)') ' error_threshold = [real] : montecarlo run until statistical error on energy reaches error_threshold'
   write(6,'(a)') ' nstep_total = [real]: For MPI, total number of steps per block for all CPUs'
   write(6,'(a)') 'end'
   write(6,*)

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
