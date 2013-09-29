module main_menu_mod

  use all_tools_mod
  use control_mod
  use basis_mod
  use orbitals_mod
  use jastrow_mod
  use psi_mod
  use periodic_jastrow_mod 
  use optimization_mod
  use grid_mod
  use dipole_moment_mod
  use density_mod
  use intracule_mod
  use extracule_mod
  use forces_mod
  use print_mod
  use testing_mod
  use dmc_mod
  use debug_mod
  use backflow_mod, only: backflow_menu !fp
  use correlated_sampling_mod
  use cusp_mod

  character(len=max_string_len_file)  :: include_file = ''
  logical                             :: in_include_file = .false.
  character(len=max_string_len)       :: current_line_save
  integer                             :: position_in_current_line_save

  contains

!===========================================================================
  subroutine main_menu
!---------------------------------------------------------------------------
! Description : read main menu
! Created     : J. Toulouse, 5 Oct 2005
! Revised     : J. Toulouse, 13 Oct 2005: improve parser
!---------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save  :: lhere = 'main_menu'
  character(len=max_string_len)  :: command

! begin

! main menu
  do

  call get_next_command (command)

  select case(trim(command))

  case ('help')
   write(6,*)
   write(6,'(a)') 'HELP for main menu:'
   write(6,'(a)') ' mode ... end: menu to control the mode of run'
   write(6,'(a)') ' include ... end: menu to include external files'
   write(6,'(a)') ' control ... end: menu for control parameters'
   write(6,'(a)') ' nuclei... end: menu for nuclei information'
   write(6,'(a)') ' basis ... end: menu to input information for basis'
   write(6,'(a)') ' orbitals ... end: menu to input information for orbitals'
   write(6,'(a)') ' csfs ... end: menu to input csfs information'
   write(6,'(a)') ' jastrow ... end: menu to input jastrow information'
   write(6,'(a)') ' wavefunction ... end: menu related to wave function '
   write(6,'(a)') ' walkers ... end: menu for walkers'
   write(6,'(a)') ' optimization ... end: menu for wave function optimization'
   write(6,'(a)') ' print ... end: menu for printing objects'
   write(6,'(a)') ' average ... end: menu for calculating MC averages'
   write(6,'(a)') ' grid ... end: menu for settinp up spatial grid'
   write(6,'(a)') ' dipole_moment ... end: menu for calculation of dipole moment'
   write(6,'(a)') ' radial_probability ... end: menu for histogram calculation of radial density of spherically symmetric systems'
   write(6,'(a)') ' density_xy_z ... end: menu for histogram calculation of density of axially symmetric systems'
   write(6,'(a)') ' density ... end: menu for ZVZB  calculation of radial density of spherically symmetric systems'
   write(6,'(a)') ' density_3d ... end: menu for ZVZB calculation of 3d density'
   write(6,'(a)') ' density_fourier ... end: menu for density in fourier space for periodic calculations'
   write(6,'(a)') ' intracule ... end: menu for calculation of 1D intracules'
   write(6,'(a)') ' intracule3d ... end: menu for calculation of 3D intracules'
   write(6,'(a)') ' extracule ... end: menu for calculation of 1D extracules'
   write(6,'(a)') ' forces    ... end: menu for calculation of forces'
   write(6,'(a)') ' cusp      ... end: menu for e-N cusp options'
   write(6,'(a)') ' debug ... end: menu for debugging'
   write(6,'(a)') ' statistics ... end: menu for printing timing statistics'
   write(6,'(a)') ' backflow ... end: menu for backflow transformation' !fp
   write(6,'(a)') ' correlated_sampling ... end: menu for correlated sampling (forces)'
   write(6,*)

  case ('test')              ; call testing
  case ('include')           ; call include_menu
  case ('mode')              ; call mode_menu
  case ('trace')             ; l_trace = .true.
  case ('control')           ; call control_menu
  case ('nuclei')            ; call nuclei_menu
  case ('basis')             ; call basis_menu
  case ('orbitals')          ; call orbitals_menu
  case ('csfs')              ; call csfs_menu
  case ('jastrow')           ; call jastrow_menu
  case ('wavefunction')      ; call wavefunction_menu
  case ('walkers')           ; call walkers_menu
  case ('periodic_jastrow')  ; call periodic_jastrow_menu
  case ('optimization')      ; call optimization_menu
  case ('backflow')          ; call backflow_menu      !fp
  case ('correlated_sampling'); call correlated_sampling_menu
  case ('print')             ; call print_menu
  case ('average')           ; call average_menu
  case ('dipole_moment')     ; call dipole_moment_menu
  case ('radial_probability'); call radial_probability_menu
  case ('density')           ; call dens_menu
  case ('density_xy_z')      ; call dens_xy_z_menu
  case ('density_3d')        ; call dens_3d_menu
  case ('density_fourier')   ; call dens_fourier_menu
  case ('intracule')         ; call intra_menu
  case ('intracule_3d')      ; call intra_3d_menu
  case ('extracule')         ; call extracule_menu
  case ('forces')            ; call forces_menu
  case ('cusp')              ; call orb_cusp_menu
  case ('grid')              ; call grid_menu
  case ('node')              ; call node_menu
  case ('debug')             ; call debug_menu
  case ('statistics')        ; call statistics_menu
  case ('norun')             ; run_done = .true.
  case ('run')               ; call run
  case ('exit')              ; write(6,*); write(6,'(a)') 'Exit of menu.';  write(6,*); exit
  case default
   call die (lhere,'unknown command >'+trim(command)+'<')
  end select

  enddo

! if we were in include file, close it and go back to main menu
  if (in_include_file) then
   call close_include
  endif

  end subroutine main_menu

!===========================================================================
  subroutine mode_menu
!---------------------------------------------------------------------------
! Description : menu for mode or 'mode'
!
! Created     : J. Toulouse, 06 Oct 2005
!---------------------------------------------------------------------------
  use contr3_mod
  implicit none

  character (len=max_string_len_rout), save :: lhere = 'mode_menu'

! begin
! loop over menu lines
  do
  call get_next_word (word)

  if(trim(word) == 'help') then
   write(6,*) 'mode-h: menu for mode'
   write(6,*) 'mode-h: mode'
   write(6,*) 'mode-h:  vmc_mov1:dmc_mov1'
   write(6,*) 'mode-h: end'

  elseif(trim(word) == 'vmc_mov1') then
   mode = trim(word)

  elseif(trim(word) == 'dmc_mov1') then
   mode = trim(word)

  elseif(trim(word) == 'end') then
   exit

  else

   write(6,*) trim(lhere),': unknown word = ',trim(word)
   call die (lhere)

  endif


  enddo ! end loop over menu lines


  end subroutine mode_menu

!===========================================================================
  subroutine include_menu
!---------------------------------------------------------------------------
! Description : menu for include
! Description : only one level of include file is allowed
!
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none


! begin
! loop over menu lines
  do
  call get_next_word (word)

  if(trim(word) == 'help') then
   write(6,*) 'include-h: menu for include'
   write(6,*) 'include-h: include'
   write(6,*) 'include-h:  [string] file to  be included'
   write(6,*) 'include-h: end'

  elseif(trim(word) == 'end') then
   exit

  else

   include_file = trim(word)
   call open_include

  endif

  enddo ! end loop over menu lines

  end subroutine include_menu

!===========================================================================
  subroutine open_include
!---------------------------------------------------------------------------
! Description : open include file
!
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none

  character (len=max_string_len_rout), save :: lhere = 'open_include'
  integer iostat

! begin
  in_include_file = .true.
  write(6,*) trim(lhere),': include file >',trim(include_file),'<'

  current_line_save = trim(current_line)
  current_line = ''
  position_in_current_line_save = position_in_current_line
  position_in_current_line = 0

  unit_input = 50
  open(file=trim(include_file), unit=unit_input, iostat=iostat)

  if (iostat  /= 0) then
   call die (lhere, 'error in openning file >'+trim(include_file)+'<.')
  endif

  call main_menu

  end subroutine open_include

!===========================================================================
  subroutine close_include
!---------------------------------------------------------------------------
! Description : close include file
!
! Created     : J. Toulouse, 13 Oct 2005
!---------------------------------------------------------------------------
  implicit none


! begin
  current_line = current_line_save
  position_in_current_line = position_in_current_line_save

  if(unit_input /= 5) then
    close(unit_input)
    unit_input = 5
  endif

  in_include_file = .false.

  end subroutine close_include

!===========================================================================
  subroutine radial_probability_menu
!---------------------------------------------------------------------------
! Description : menu for radial probability for atoms
!
! Created     : J. Toulouse, 22 Oct 2005
!---------------------------------------------------------------------------
  use stepv_mod
  implicit none

  character (len=max_string_len_rout), save :: lhere = 'radial_probability_menu'

! begin
! loop over menu lines
  do
  call get_next_word (word)

  if(trim(word) == 'help') then
   write(6,'(a)') 'menu for radial probability:'
   write(6,'(a)') ' radial_probability'
   write(6,'(a)') ' print = [boolean] print radial probability?'
   write(6,'(a)') 'end'

  elseif(trim(word) == 'print') then
   call get_next_value (print_radial_probability)

  elseif(trim(word) == 'end') then
   exit

  else
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  endif

  enddo ! end loop over menu lines

  end subroutine radial_probability_menu


!===========================================================================
  subroutine node_menu
!---------------------------------------------------------------------------
! Description : write info for nodes
!
! Created     : J. Toulouse, 16 Dec 2005
!---------------------------------------------------------------------------
  implicit none

  character (len=max_string_len_rout), save :: lhere = 'node_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  if(trim(word) == 'help') then
   write(6,*) 'node-h: menu for printing info on nodes:'
   write(6,*) 'node-h: node'
   write(6,*) 'node-h:  node_name'
   write(6,*) 'node-h: end'

  elseif(.not. trim(word) == 'end') then

   call node_info (word)

  elseif(trim(word) == 'end') then
   exit

  else

   write(6,*) trim(lhere),': unknown word = ',trim(word)
   call die (lhere)

  endif


  enddo ! end loop over menu lines

  end subroutine node_menu

!===========================================================================
  subroutine statistics_menu
!---------------------------------------------------------------------------
! Description : menu for printing statistics
!
! Created     : J. Toulouse, 18 Dec 2005
!---------------------------------------------------------------------------
  implicit none

  character*(max_string_len_rout), save :: lhere = 'statistics_menu'
  character*(max_string_len) word

! begin

! loop over menu lines
  do
  call get_next_word (word)

  if(trim(word) == 'help') then
   write(6,'(2a)') trim(lhere), ': help for menu statistics:'
   write(6,'(2a)') trim(lhere), ': statistics'
   write(6,'(2a)') trim(lhere), ':   nodes: print statistics for each executed node routine'
   write(6,'(2a)') trim(lhere), ':   routines: print statistics for each executed routine'
   write(6,'(2a)') trim(lhere), ': end'

  elseif(trim(word) == 'nodes') then
   call nodes_statistics

  elseif(trim(word) == 'routines') then
   call routines_statistics

  elseif(trim(word) == 'end') then
   exit

  else

   write(6,'(3a)') trim(lhere),': unknown word = ',trim(word)
   call die (lhere)

  endif


  enddo ! end loop over menu lines

  end subroutine statistics_menu


!===============================================================================
  subroutine node_info (node_name)
! ------------------------------------------------------------------------------
! Description   :
!
! Created       : J. Toulouse, 14 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'node_info'
  character(len=max_string_len) :: node_name
  integer node_ind, obj_i

! begin
! index of node
  node_ind = node_index (node_name)

! if node not found, die
   if (node_ind == 0) then
      write(6,*) trim(lhere),': entering node ', trim(node_name),' that not catalogued'
      call die (lhere)
   endif

   write(6,'(a,a,i3,a,a)') trim(lhere),': node # ', node_ind,': ', trim(node_name)

   do obj_i = 1, size(nodes(node_ind)%objects_create_index (:))
    write(6,'(3a)') trim(lhere),': created object ', objects (nodes(node_ind)%objects_create_index (obj_i))%name
   enddo

   write(6,*)

   do obj_i = 1, size(nodes(node_ind)%objects_needed_index (:))
    write(6,'(3a)') trim(lhere),': needed object ', objects (nodes(node_ind)%objects_needed_index (obj_i))%name
   enddo

   write(6,*)

   write(6,'(a,a,i3)') trim(lhere),': number of calls =',nodes(node_ind)%calls_nb


   write(6,'(a,a,es15.8)') trim(lhere),': cpu duration =',nodes(node_ind)%cpu_duration

   write(6,*) trim(lhere),': cpu duration =', trim(cpu_to_string( nodes(node_ind)%cpu_duration))


 end subroutine node_info

end module main_menu_mod
