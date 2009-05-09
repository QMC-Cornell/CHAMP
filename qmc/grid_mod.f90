module grid_mod

  use all_tools_mod
  use orbitals_mod

! Declaration of global variables and default values

! grid r
  real(dp)                  :: grid_r_step
  real(dp)                  :: grid_r_min = 0.d0
  real(dp)                  :: grid_r_max
  integer                   :: grid_r_nb = 0.d0
  real(dp), allocatable     :: grid_r (:)

! grid xyz
  real(dp)                  :: grid_x_step
  real(dp)                  :: grid_x_max
  real(dp)                  :: grid_x_min
  integer                   :: grid_x_nb = 0.d0
  real(dp)                  :: grid_y_step
  real(dp)                  :: grid_y_max
  real(dp)                  :: grid_y_min
  integer                   :: grid_y_nb = 0.d0
  real(dp)                  :: grid_z_step
  real(dp)                  :: grid_z_max
  real(dp)                  :: grid_z_min
  integer                   :: grid_z_nb = 0.d0
  integer                   :: grid_xyz_nb = 0.d0
  real(dp), allocatable     :: grid_xyz (:,:)
  real(dp), allocatable     :: grid_xyz_index (:,:,:)

! old grid for orbitals
  integer                   :: grid_on_x_nb
  real(dp)                  :: grid_max
  real(dp)                  :: grid_step
  real(dp), allocatable     :: grid_on_x (:)
  character*(max_string_len_file)   :: orb_grid_file = ''
  real(dp), allocatable     :: orb_on_x (:,:)

  contains

!===========================================================================
  subroutine grid_menu
!---------------------------------------------------------------------------
! Description : menu for grids
!
! Created     : J. Toulouse, 04 Mar 2006
!---------------------------------------------------------------------------
  implicit none

  character(len=max_string_len_rout), save :: lhere = 'grid_menu'

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of grid menu -----------------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)') 'HELP for grid menu:'
   write(6,'(a)') ' grid'
   write(6,'(a)') '  grid_r ... end   : submenu for grid over spatial radial coordinate r'
   write(6,'(a)') '  grid_xyz ... end : submenu for grid over spatial cartesian coordinates x y z'
   write(6,'(a)') ' end'

  case ('grid_r')
   call grid_r_menu

  case ('grid_xyz')
   call grid_xyz_menu

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines

  write(6,'(a)') 'End of grid menu -----------------------------------------------------------------------------------------'

  end subroutine grid_menu

!===========================================================================
  subroutine grid_r_menu
!---------------------------------------------------------------------------
! Description : menu for radial grid
!
! Created     : J. Toulouse, 04 Mar 2006
!---------------------------------------------------------------------------
  implicit none

  character(len=max_string_len_rout), save :: lhere = 'grid_r_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)') 'HELP for menu grid_r:'
   write(6,'(a)') ' grid_r'
   write(6,'(a)') '  grid_r_step   = [real] : grid spacing'
   write(6,'(a)') '  grid_r_min    = [real] : grid min value (default = 0.)'
   write(6,'(a)') '  grid_r_max    = [real] : grid max value'
   write(6,'(a)') ' end'

  case ('grid_r_step')
   call get_next_value (grid_r_step)

  case ('grid_r_min')
   call get_next_value (grid_r_min)

  case ('grid_r_max')
   call get_next_value (grid_r_max)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines

! Parameters for grid over r
  call require (lhere, 'grid_r_step > 0', grid_r_step > 0.d0 )
  call require (lhere, 'grid_r_max > 0', grid_r_max > 0.d0 )
  call require (lhere, 'grid_r_min >= 0', grid_r_min >= 0.d0 )

  grid_r_nb = int((grid_r_max-grid_r_min)/grid_r_step) + 1
  call require (lhere, 'grid_r_nb > 0', grid_r_nb > 0 )

  write (6,'(a)') ' Parameters for radial grid:'
  write (6,'(a,es15.8)') ' grid_r_step = ', grid_r_step
  write (6,'(a,es15.8)') ' grid_r_min = ', grid_r_min
  write (6,'(a,es15.8)') ' grid_r_max = ', grid_r_max
  write (6,'(a,i8)') ' grid_r_nb=', grid_r_nb

  call object_modified ('grid_r_step')
  call object_modified ('grid_r_min')
  call object_modified ('grid_r_max')
  call object_modified ('grid_r_nb')

  end subroutine grid_r_menu

!===========================================================================
  subroutine grid_xyz_menu
!---------------------------------------------------------------------------
! Description : menu for grid over x, y, z
!
! Created     : J. Toulouse, 04 Mar 2006
!---------------------------------------------------------------------------
  implicit none

  character(len=max_string_len_rout), save :: lhere = 'grid_xyz_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)') 'HELP for menu grid_xyz:'
   write(6,'(a)') ' grid_xyz'
   write(6,'(a)') '  grid_x_step   = [real] grid x spacing'
   write(6,'(a)') '  grid_x_max    = [real] grid x max value'
   write(6,'(a)') '  grid_x_min    = [real] grid x min value'
   write(6,'(a)') '  grid_y_step   = [real] grid y spacing'
   write(6,'(a)') '  grid_y_max    = [real] grid y max value'
   write(6,'(a)') '  grid_y_min    = [real] grid y min value'
   write(6,'(a)') '  grid_z_step   = [real] grid z spacing'
   write(6,'(a)') '  grid_z_max    = [real] grid z max value'
   write(6,'(a)') '  grid_z_min    = [real] grid z min value'
   write(6,'(a)') ' end'

  case ('grid_x_step')
   call get_next_value (grid_x_step)
  case ('grid_x_max')
   call get_next_value (grid_x_max)
  case ('grid_x_min')
   call get_next_value (grid_x_min)

  case ('grid_y_step')
   call get_next_value (grid_y_step)
  case ('grid_y_max')
   call get_next_value (grid_y_max)
  case ('grid_y_min')
   call get_next_value (grid_y_min)

  case ('grid_z_step')
   call get_next_value (grid_z_step)
  case ('grid_z_max')
   call get_next_value (grid_z_max)
  case ('grid_z_min')
   call get_next_value (grid_z_min)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines

! Parameters for grid xyz
  call require (lhere, 'grid_x_step > 0', grid_x_step > 0.d0 ) !fp 
  call require (lhere, 'grid_x_min <= grid_x_max ', grid_x_min <= grid_x_max) !fp
  call require (lhere, 'grid_y_step > 0', grid_y_step > 0.d0 ) !fp
  call require (lhere, 'grid_y_min <= grid_y_max ', grid_y_min <= grid_y_max) !fp
  call require (lhere, 'grid_z_step > 0', grid_z_step > 0.d0 ) !fp
  call require (lhere, 'grid_z_min <= grid_z_max ', grid_z_min <= grid_z_max) !fp

  grid_x_nb = int((grid_x_max - grid_x_min)/grid_x_step) + 1
  grid_y_nb = int((grid_y_max - grid_y_min)/grid_y_step) + 1
  grid_z_nb = int((grid_z_max - grid_z_min)/grid_z_step) + 1
  grid_xyz_nb = grid_x_nb * grid_y_nb * grid_z_nb
  call require (lhere, 'grid_x_nb > 0', grid_x_nb > 0 )
  call require (lhere, 'grid_y_nb > 0', grid_y_nb > 0 )
  call require (lhere, 'grid_z_nb > 0', grid_z_nb > 0 )
  call require (lhere, 'grid_xyz_nb > 0', grid_xyz_nb > 0 )

  write (6,'(a)')  ' Parameters for grid xyz:'
  write (6,'(a,es15.8)')  ' grid_x_step = ', grid_x_step
  write (6,'(a,es15.8)')  ' grid_x_min  = ', grid_x_min
  write (6,'(a,es15.8)')  ' grid_x_max  = ', grid_x_max
  write (6,'(a,es15.8)')  ' grid_x_nb   = ', grid_x_nb
  write (6,'(a,es15.8)')  ' grid_y_step = ', grid_y_step
  write (6,'(a,es15.8)')  ' grid_y_min  = ', grid_y_min
  write (6,'(a,es15.8)')  ' grid_y_max  = ', grid_y_max
  write (6,'(a,es15.8)')  ' grid_y_nb   = ', grid_y_nb
  write (6,'(a,es15.8)')  ' grid_z_step = ', grid_z_step
  write (6,'(a,es15.8)')  ' grid_z_min  = ', grid_z_min
  write (6,'(a,es15.8)')  ' grid_z_max  = ', grid_z_max
  write (6,'(a,es15.8)')  ' grid_z_nb   = ', grid_z_nb
  write (6,'(a,es15.8)')  ' grid_xyz_nb = ', grid_xyz_nb

  call object_modified ('grid_x_step')
  call object_modified ('grid_x_max')
  call object_modified ('grid_x_min')
  call object_modified ('grid_x_nb')
  call object_modified ('grid_y_step')
  call object_modified ('grid_y_max')
  call object_modified ('grid_y_min')
  call object_modified ('grid_y_nb')
  call object_modified ('grid_z_step')
  call object_modified ('grid_z_max')
  call object_modified ('grid_z_min')
  call object_modified ('grid_z_nb')
  call object_modified ('grid_xyz_nb')

  end subroutine grid_xyz_menu

! ==============================================================================
  subroutine grid_r_bld
! ------------------------------------------------------------------------------
! Descriptions   : build grid over r
!
! Created        : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  integer grid_i
  character(len=max_string_len_rout), save :: lhere = 'grid_r_bld'

! header
  if (header_exe) then

   call object_create ('grid_r')

   call object_needed ('grid_r_nb')
   call object_needed ('grid_r_step')
   call object_needed ('grid_r_min')

   return

  endif

! begin
  call require (lhere, 'grid_r_nb > 0', grid_r_nb > 0)
  call require (lhere, 'grid_r_step > 0', grid_r_step > 0)

! allocation
  call object_alloc ('grid_r', grid_r, grid_r_nb)

  do grid_i = 1, grid_r_nb
   grid_r (grid_i) = grid_r_min + (grid_i - 1) * grid_r_step
  enddo

  end subroutine grid_r_bld

! ==============================================================================
  subroutine grid_xyz_bld
! ------------------------------------------------------------------------------
!
! Created        : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'grid_xyz_bld'
  integer grid_i, grid_x_i, grid_y_i, grid_z_i

! header
  if (header_exe) then

   call object_create ('grid_xyz')
   call object_create ('grid_xyz_index')

   call object_needed ('grid_x_nb')
   call object_needed ('grid_x_min')
   call object_needed ('grid_x_step')
   call object_needed ('grid_y_nb')
   call object_needed ('grid_y_min')
   call object_needed ('grid_y_step')
   call object_needed ('grid_z_nb')
   call object_needed ('grid_z_min')
   call object_needed ('grid_z_step')
   call object_needed ('grid_xyz_nb')

   return

  endif

! begin

! allocation
  call object_alloc ('grid_xyz', grid_xyz, 3, grid_xyz_nb)
  call object_alloc ('grid_xyz_index', grid_xyz_index, grid_x_nb, grid_y_nb, grid_z_nb)

  grid_i = 0

  do grid_x_i = 1, grid_x_nb
   do grid_y_i = 1, grid_y_nb
    do grid_z_i = 1, grid_z_nb
     grid_i = grid_i + 1
     grid_xyz_index (grid_x_i, grid_y_i, grid_z_i) = grid_i
     grid_xyz (1, grid_i) = grid_x_min + grid_x_step * (grid_x_i - 1)
     grid_xyz (2, grid_i) = grid_y_min + grid_y_step * (grid_y_i - 1)
     grid_xyz (3, grid_i) = grid_z_min + grid_z_step * (grid_z_i - 1)
    enddo
   enddo
  enddo

  call require (lhere, 'grid_i == grid_xyz_nb', grid_i == grid_xyz_nb)

  end subroutine grid_xyz_bld

!===========================================================================
  subroutine grid_orb_menu
!---------------------------------------------------------------------------
! Description : menu for grid
!
! Created     : J. Toulouse, 05 Jan 2005
!---------------------------------------------------------------------------
  implicit none

  character*(max_string_len_rout) here
  character*(max_string_len) word

! begin
  here = 'grid_orb_menu'

! loop over menu lines
  do
  call get_next_word (word)

  if(trim(word) == 'help') then
   write(6,*) 'grid-h: menu for grid'
   write(6,*) 'grid-h: grid'
   write(6,*) 'grid-h:  max = [real] maximal value of the grid'
   write(6,*) 'grid-h:  step = [real] step of the grid'
   write(6,*) 'grid-h:  orbitals = [string] file to write orbitals'
   write(6,*) 'grid-h: end'

  elseif(trim(word) == 'max') then
   call get_next_value (grid_max)
   call object_modified ('grid_max')

  elseif(trim(word) == 'step') then
   call get_next_value (grid_step)
   call object_modified ('grid_step')

  elseif(trim(word) == 'orbitals') then
   call get_next_value (orb_grid_file)
   call orb_on_x_wrt

  elseif(trim(word) == 'end') then
   exit

  else

   write(6,*) trim(here),'unknown word = ',trim(word)
   call die(here)

  endif

  enddo ! end loop over menu lines

  end subroutine grid_orb_menu

! ==============================================================================
  subroutine grid_on_x_bld
! ------------------------------------------------------------------------------
! Description   : build grid on x coordinate
!
! Created       : J. Toulouse, 05 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: here = 'grid_on_x_bld'
  integer i

! header
  if (header_exe) then

   call object_create ('grid_on_x')
   call object_create ('grid_on_x_nb')

   call object_needed ('grid_max')
   call object_needed ('grid_step')

   return

  endif

! begin

! requirements
  if (grid_max < 0) then
   write(6,*) trim(here),': grid_max',grid_max,' <= 0'
   call die (here)
  endif
  if (grid_step < 0) then
   write(6,*) trim(here),': grid_step',grid_step,' <= 0'
   call die (here)
  endif

! number of grid points
  grid_on_x_nb = int(grid_max/grid_step) + 1

! allocations
  call object_alloc ('grid_on_x', grid_on_x, grid_on_x_nb)

  do i = 1, grid_on_x_nb
     grid_on_x (i) = grid_step * (i-1)
  enddo

  end subroutine grid_on_x_bld

! ==============================================================================
  subroutine orb_on_x_bld
! ------------------------------------------------------------------------------
! Description   : build orbitals on a x-grid
!
! Created       : J. Toulouse, 05 Jan 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: here = 'orb_on_x_bld'
  integer orb_i, grd_i
  integer iel
  real(dp) orbl(MORB)

! header
  if (header_exe) then

   call object_create ('orb_on_x')

   call object_needed ('orb_tot_nb')
   call object_needed ('grid_on_x_nb')
   call object_needed ('grid_on_x')

   return

  endif

! begin
!  write(6,*) trim(here),': entering'

! requirements
  if (nbasis <= 0) then
   write(6,*) trim(here),': nbasis=',nbasis,' <= 0'
   call die (here)
  endif

  if (ncent /= 1 ) then
   write(6,*) trim(here),': ncent=',ncent,' /= 1'
   write(6,*) trim(here),': implemented only for one center'
   call die (here)
  endif

   if (iperiodic /= 0) then
   write(6,*) trim(here),': iperiodic=',iperiodic,' /= 0'
   write(6,*) trim(here),': implemented only for no periodic case'
   call die (here)
   endif

   if (inum_orb /= 0) then
   write(6,*) trim(here),': inum_orb=',inum_orb,' /= 0'
   write(6,*) trim(here),': implemented only for no numerical orbitals'
   call die (here)
   endif
! end requirements

! allocations
  call alloc ('orb_on_x', orb_on_x, orb_tot_nb, grid_on_x_nb)

  iel = 1
  rvec_en = 0.d0
  r_en = 0.d0
  iwf = 1

  do grd_i = 1, grid_on_x_nb
    rvec_en (1,1,1) = grid_on_x (grd_i)
    r_en (1,1) = grid_on_x (grd_i)
    call orbitals_loc_anae(iel,rvec_en,r_en,orbl)
    write(6,*) orbl(1)
    do orb_i = 1, orb_tot_nb
      orb_on_x (orb_i, grd_i) = orbl (orb_i)
    enddo
  enddo

 end subroutine orb_on_x_bld

! ==============================================================================
  subroutine orb_on_x_wrt
! ------------------------------------------------------------------------------
! Description   : build orbitals on a x-grid
!
! Created       : J. Toulouse, 05 Jan 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_file) file
  integer orb_i, grd_i
  integer unit

! begin
  call object_provide ('orb_tot_nb')
  call object_provide ('grid_on_x_nb')
  call object_provide ('grid_on_x')
  call object_provide ('orb_on_x')

  file = orb_grid_file
  unit = 20
  open(file=trim(file),unit=unit)

  do grd_i = 1, grid_on_x_nb
   write(unit,'(es15.8,100es15.8)') grid_on_x (grd_i),  (orb_on_x (orb_i, grd_i), orb_i=1,orb_tot_nb)
  enddo

  close (unit)

  end subroutine orb_on_x_wrt

end module grid_mod
