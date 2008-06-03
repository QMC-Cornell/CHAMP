program champ

! modules use
  use all_tools_mod
  use f2kcli
  use main_menu_mod
  use catalog_routines_mod
  use initialization_mod
  use mpi_mod
  use dmc_mod

  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'champ'
  character (len=max_string_len) :: date = ''
  character*16                             :: mode
  integer            :: narg, iarg
  character(len=256) :: command_line_arguments
  character(len=256) :: executable_name
  character(len=256) :: argument
  character(len=max_string_len_file)       :: input_file_name = ''
  integer iostat, mode_i

  integer ierr
  common /contr3/ mode

! initialization
  mode = 'vmc_mov1'

# if defined (MPI)
  call mpi_init (ierr)
  if (ierr /= 0) then
   call die (lhere, 'error in mpi_init.')
  endif

  call mpi_comm_rank (MPI_COMM_WORLD,idtask,ierr)
  if (ierr /= 0) then
   call die (lhere, 'error in mpi_comm_rank.')
  endif

  call mpi_comm_size (MPI_COMM_WORLD,nproc,ierr)
  if (ierr /= 0) then
   call die (lhere, 'error in mpi_comm_size.')
  endif

  if(nproc > MPROC) then
   call die (lhere, 'nproc > MPROC.')
  endif

! Close standard output if not the master
  if (idtask.ne.0) then
   close(6)
! Warning temporarily commented out
  open(6,file='/dev/null')
!   open(6,file='slave.out')
  endif
# endif

  call get_date (date)
  write(6,'(2a)') 'PROGRAM CHAMP run on ',date
  write(6,'(2a)') 'SVN VERSION $Rev$ on $Date$'

# if defined (MPI)
  write(6,'(a,i4,a)') 'MPI version running on ', nproc, ' processors.'
# endif

  call object_modified ('nproc')

! begin
  call cpu_time (cpu_time_start)
  cpu_time_last = cpu_time_start

! number of arguments on command line
  narg = command_argument_count()

! get command line arguments
  call get_command (command_line_arguments)

! get executable name
  call get_command_argument (0, executable_name)

  write(6,'(2a)') 'Executable: ',trim(executable_name)
  write(6,'(2a)') 'Command line arguments: ',trim(command_line_arguments)

! loop over arguments on command line
  iarg = 0
  do
   iarg = iarg + 1
   if (iarg > narg) exit
   call get_command_argument (iarg, argument)

!  use parser for input
   if (trim(argument) == '-p' .or. trim(argument) == '-parser') then
    use_parser = .true.
   endif

!  specified mode of champ
   if (trim(argument) == '-m' .or. trim(argument) == '-mode') then
    iarg = iarg + 1
    call get_command_argument (iarg, argument)
    if (iarg > narg) then
     call die (lhere, 'value for option "-m{ode}" missing.')
    endif
    mode = argument
   endif

!  input file (necessary for MPI version)
   if (trim(argument) == '-i' .or. trim(argument) == '-input') then
    iarg = iarg + 1
    call get_command_argument (iarg, argument)
    if (iarg > narg) then
     call die (lhere, 'value for option "-i{input}" missing.')
    endif
    input_file_name  = argument
   endif

  enddo ! end loop over arguments

! check mode
  write(6,'(3a)') 'The mode is >',trim(mode),'<.'
  if (.not. elt_in_array (modes, mode)) then
   write(6,'(3a)') 'This mode is unknown.'
   write(6,'(20a)') 'The available modes are:'
   do mode_i = 1, modes_nb
     write(6,'(2x,a)') trim(modes (mode_i))
   enddo
   call die (lhere)
  endif

# if defined (MPI)
!  check mode
   if (index(mode, 'mpi') == 0) then
     write(6,'(4a)') trim(lhere), ': mode = ',trim(mode), ' is incorrect for the MPI version.'
     write(6,'(2a)') trim(lhere), ': the available modes for the MPI version are fit_mpi, vmc_mpi, vmc_mov1_mpi, dmc_mov1_mpi1, dmc_mov1_mpi2 and dmc_mov1_mpi3.'
     call die (lhere)
   endif

! check input file
   if (trim(input_file_name) == '') then
     call die (lhere, 'input file name must be given in the command line by "-i{nput} filename"')
   endif
# endif

! Open input file if given by '-i ...'
  write(6,*)
  if (trim(input_file_name) /= '') then
     write(6,'(3a)') 'Opening input file >',trim(input_file_name),'<.'
     open(5, file=trim(input_file_name), status='old', iostat=iostat)
     if (iostat /= 0) then
       call die (lhere, 'error on opening file >'+trim(input_file_name)+'<')
     endif
  endif

! Build production tree
  call build_tree

! pre-catalog some objects for speed
  call catalog_objects

! Define all averages and errors
  call define_averages_and_errors

! Read input
  write(6,'(2a)') 'Reading input file...'
  write(6,*)

  if (trim(mode) == 'parser') then
    call main_menu
  else
    call read_input
    if (index(mode,'fit') /= 0) then
      call fit
      run_done = .true.
    else
      call read_up_to_end

!     initialization of some global variables
      call initialization
      call main_menu
    endif

!   do default run is not already done
    if(.not. run_done) then
      call run_default
    endif

  endif

   if (l_warning) then
    write(6,*)
    write(6,'(a)') 'Some warnings were encountered.'
   endif

   write(6,*)
   call print_cpu_time
   write(6,'(a)') 'The program ended normally.'

# if defined (MPI)
   call mpi_finalize (ierr)
   if (ierr /= 0) then
    call die (lhere, 'error in mpi_finalize.')
   endif
# endif

end program champ

!---------------------------------------------------------------------------
  subroutine run
!---------------------------------------------------------------------------
! Description : montecarlo run
!
! Created     : J. Toulouse, 01 Apr 2007
!---------------------------------------------------------------------------
  use main_menu_mod

  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save  :: lhere = 'run'

! begin
  if (l_mode_fit) then
   call fit

  elseif (l_mode_vmc) then
   call vmc_run

  elseif (l_mode_dmc .and. .not. l_mode_mpi) then
   call maindmc
!   call dmc_mov1_clean
!   call dmc_mov1

  elseif (l_mode_dmc_mov1_mpi1 .or. l_mode_dmc_mov1_mpi2 .or. l_mode_dmc_mov1_mpi3) then
   call maindmc_mov1_mpi

  else
   call die (lhere, 'unknown mode >'+trim(mode)+'<')
  endif

  run_done = .true.

  end subroutine run

!---------------------------------------------------------------------------
  subroutine run_default
!---------------------------------------------------------------------------
! Description : default montecarlo run
!
! Created     : J. Toulouse, 12 Feb 2008
!---------------------------------------------------------------------------
  use main_menu_mod

  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save  :: lhere = 'run'

! begin
  if (l_mode_fit) then
   call fit

  elseif (l_mode_vmc) then
   call opt_wf

  elseif (l_mode_dmc .and. .not. l_mode_mpi) then
   call maindmc
!  call dmc_mov1_clean
!  call dmc_mov1

  elseif (l_mode_dmc_mov1_mpi1 .or. l_mode_dmc_mov1_mpi2 .or. l_mode_dmc_mov1_mpi3) then
   call maindmc_mov1_mpi

  else
   call die (lhere, 'unknown mode >'+trim(mode)+'<')
  endif

  run_done = .true.

  end subroutine run_default
