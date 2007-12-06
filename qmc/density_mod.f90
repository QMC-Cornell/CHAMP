module density_mod

  use all_tools_mod
  use grid_mod
  use electrons_mod
  use psi_mod

! Declaration of global variables and default values
  character(len=max_string_len_file):: dens_file_out  = ''
  character(len=max_string_len)     :: dens_estimator  = 'zv1'
  real(dp), allocatable     :: dens_zv1 (:)
  real(dp), allocatable     :: dens_zv1_av (:)
  real(dp), allocatable     :: dens_zv1_av_err (:)
  real(dp), allocatable     :: dens (:)
  real(dp), allocatable     :: dens_err (:)

  character(len=max_string_len_file):: dens_3d_file_out  = ''
  character(len=max_string_len)     :: dens_3d_estimator  = 'zv2'

  real(dp), allocatable     :: dens_3d_histo (:)
  real(dp), allocatable     :: dens_3d_histo_av (:)
  real(dp), allocatable     :: dens_3d_histo_av_err (:)
  real(dp), allocatable     :: dens_3d_zv1 (:)
  real(dp), allocatable     :: dens_3d_zv1_av (:)
  real(dp), allocatable     :: dens_3d_zv1_av_err (:)
  real(dp), allocatable     :: dens_3d_zv2 (:)
  real(dp), allocatable     :: dens_3d_zv2_av (:)
  real(dp), allocatable     :: dens_3d_zv2_av_err (:)
  real(dp), allocatable     :: dens_3d (:)
  real(dp), allocatable     :: dens_3d_err (:)

  contains

!===========================================================================
  subroutine dens_menu
!---------------------------------------------------------------------------
! Description : menu for spherically averaged density calculation
!
! Created     : J. Toulouse, 07 Nov 2006
!---------------------------------------------------------------------------
  implicit none

  character(len=max_string_len_rout), save :: lhere = 'dens_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(2a)') trim(lhere),': menu for density'
   write(6,'(2a)') trim(lhere),': density'
   write(6,'(2a)') trim(lhere),':  estimator = [string] {histogram|zv1|zv2}'
   write(6,'(2a)') trim(lhere),':  file    = [string] file in which density will be written'
   write(6,'(2a)') trim(lhere),': end'

  case ('file')
   call get_next_value (dens_file_out)

  case ('estimator')
   call get_next_value (dens_estimator)

  case ('end')
   exit

  case default
   write(6,'(3a)') trim(lhere),': unknown keyword = ',trim(word)
   call die (lhere)
  end select

  enddo ! end loop over menu lines


! File
  if (trim(dens_file_out) /= '') then
   write (6,'(3a)') trim(lhere),': density will be written on file: ',trim(dens_file_out)
  else
   write (6,'(3a)') trim(lhere),': file for writing density not specified'
   call die (lhere)
  endif

  select case(trim(dens_estimator))
   case ('histogram')
   call object_average_request ('dens_histo_av')
   call object_error_request ('dens_histo_av_err')

   case ('zv1')
   call object_average_request ('dens_zv1_av')
   call object_error_request ('dens_zv1_av_err')

   case ('zv2')
   call object_average_request ('dens_zv2_av')
   call object_error_request ('dens_zv2_av_err')

   case default
   write(6,'(4a)') trim(lhere), ': unknown estimator >',trim(dens_estimator),'<'
   call die (lhere)
  end select

  call routine_write  ('dens_wrt')

  end subroutine dens_menu

!===========================================================================
  subroutine dens_3d_menu
!---------------------------------------------------------------------------
! Description : menu for 3D density calculation
!
! Created     : J. Toulouse, 07 Mar 2006
!---------------------------------------------------------------------------
  implicit none

  character(len=max_string_len_rout), save :: lhere = 'dens_3d_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(2a)') trim(lhere),': menu for dens_3d'
   write(6,'(2a)') trim(lhere),': density_3d'
   write(6,'(2a)') trim(lhere),':  estimator = [string] {histogram|zv1|zv2}'
   write(6,'(2a)') trim(lhere),':  file    = [string] file in which density will be written'
   write(6,'(2a)') trim(lhere),': end'

  case ('file')
   call get_next_value (dens_3d_file_out)

  case ('estimator')
   call get_next_value (dens_3d_estimator)

  case ('end')
   exit

  case default
   write(6,'(3a)') trim(lhere),': unknown keyword = ',trim(word)
   call die (lhere)
  end select

  enddo ! end loop over menu lines


! File
  if (trim(dens_3d_file_out) /= '') then
   write (6,'(3a)') trim(lhere),': density will be written on file: ',trim(dens_3d_file_out)
  else
   write (6,'(3a)') trim(lhere),': file for writing 3D density not specified'
   call die (lhere)
  endif

  select case(trim(dens_3d_estimator))
   case ('histogram')
   call object_average_request ('dens_3d_histo_av')
   call object_error_request ('dens_3d_histo_av_err')

   case ('zv1')
   call object_average_request ('dens_3d_zv1_av')
   call object_error_request ('dens_3d_zv1_av_err')

   case ('zv2')
   call object_average_request ('dens_3d_zv2_av')
   call object_error_request ('dens_3d_zv2_av_err')

   case default
   write(6,'(4a)') trim(lhere), ': unknown estimator >',trim(dens_3d_estimator),'<'
   call die (lhere)
  end select

  call routine_write  ('dens_3d_wrt')

  end subroutine dens_3d_menu

! ==============================================================================
  subroutine dens_zv1_bld
! ------------------------------------------------------------------------------
! Description   : first-order renormalized improved estimator of spherically averaged density
!
! Created       : J. Toulouse, 07 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, dim_i
  real(dp)              :: di, r
  real(dp)              :: dotproduct, dens_temp

! begin

! header
  if (header_exe) then

   call object_create ('dens_zv1')

   call object_needed ('grid_r_nb')
   call object_needed ('grid_r')
   call object_needed ('nelec')
   call object_needed ('xold')
   call object_needed ('dist_e')
   call object_needed ('grd_psi_over_psi_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_zv1', dens_zv1, grid_r_nb)
  call object_alloc ('dens_zv1_av', dens_zv1_av, grid_r_nb)
  call object_alloc ('dens_zv1_av_err', dens_zv1_av_err, grid_r_nb)

  dens_zv1 (:) = 0.d0

  do elec_i = 1, nelec

!     distance |r_i|
      di = dist_e (elec_i)

!     dot product: drift_i . r_i
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_i, 1) * xold (dim_i, elec_i)
      enddo

      dens_temp = - oneover2pi * dotproduct / di**3

       do grid_i = 1, grid_r_nb
        r = grid_r (grid_i)

        if ( di >= r ) then
           dens_zv1 (grid_i) = dens_zv1 (grid_i) + dens_temp
        endif

     enddo !grid_i

  enddo !elec_i

  end subroutine dens_zv1_bld

! ==============================================================================
  subroutine dens_bld
! ------------------------------------------------------------------------------
! Description   : spherically averaged density density  n(r)
!
! Created       : J. Toulouse, 07 Nov 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'dens_bld'

! header
  if (header_exe) then

   call object_create ('dens')
   call object_create ('dens_err')

   call object_needed ('grid_r_nb')

   return

  endif

! begin

! allocations
  call object_alloc ('dens', dens, grid_r_nb)
  call object_alloc ('dens_err', dens_err, grid_r_nb)

  select case(trim(dens_estimator))
   case ('zv1')
   call object_provide (lhere,'dens_zv1_av')
   call object_provide (lhere,'dens_zv1_av_err')
   dens (:)     = dens_zv1_av (:)
   dens_err (:) = dens_zv1_av_err (:)

   case default
   write(6,'(4a)') trim(lhere), ': unknown estimator >',trim(dens_estimator),'<'
   call die (lhere)
  end select

  end subroutine dens_bld

! ==============================================================================
  subroutine dens_3d_histo_bld
! ------------------------------------------------------------------------------
! Description   : histo estimator of 3D density
!
! Created       : J. Toulouse, 05 Mar 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: grid_x_i, grid_y_i, grid_z_i
  integer                       :: elec_i, dim_i
  real(dp)              :: xi, yi, zi

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_histo')

   call object_needed ('grid_xyz_nb')
   call object_needed ('grid_xyz')
   call object_needed ('grid_xyz_index')
   call object_needed ('grid_x_nb')
   call object_needed ('grid_y_nb')
   call object_needed ('grid_z_nb')
   call object_needed ('grid_x_step')
   call object_needed ('grid_y_step')
   call object_needed ('grid_z_step')
   call object_needed ('grid_x_min')
   call object_needed ('grid_y_min')
   call object_needed ('grid_z_min')
   call object_needed ('nelec')
   call object_needed ('xold')

   return

  endif

! allocations
  call object_alloc ('dens_3d_histo', dens_3d_histo, grid_xyz_nb)
  call object_alloc ('dens_3d_histo_av', dens_3d_histo_av, grid_xyz_nb)
  call object_alloc ('dens_3d_histo_av_err', dens_3d_histo_av_err, grid_xyz_nb)

  dens_3d_histo (:) = 0.d0

  do elec_i = 1, nelec

      xi =  xold(1,elec_i)
      yi =  xold(2,elec_i)
      zi =  xold(3,elec_i)

      grid_x_i = floor((xi - grid_x_min)/grid_x_step - 0.5d0)  + 2
      grid_y_i = floor((yi - grid_y_min)/grid_y_step - 0.5d0)  + 2
      grid_z_i = floor((zi - grid_z_min)/grid_z_step - 0.5d0)  + 2

      if (grid_x_i < 1 .or. grid_x_i > grid_x_nb) cycle
      if (grid_y_i < 1 .or. grid_y_i > grid_y_nb) cycle
      if (grid_z_i < 1 .or. grid_z_i > grid_z_nb) cycle

      grid_i = grid_xyz_index (grid_x_i, grid_y_i, grid_z_i)

      dens_3d_histo (grid_i) = dens_3d_histo (grid_i) + 1.d0/(grid_x_step*grid_y_step*grid_z_step)

  enddo !elec_i

  end subroutine dens_3d_histo_bld

! ==============================================================================
  subroutine dens_3d_zv1_bld
! ------------------------------------------------------------------------------
! Description   : first-order renormalized improved estimator of 3D density
!
! Created       : J. Toulouse, 07 Mar 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, dim_i
  real(dp)              :: di, di2, r
  real(dp)              :: dotproduct, dens_temp

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_zv1')

   call object_needed ('grid_xyz_nb')
   call object_needed ('grid_xyz')
   call object_needed ('nelec')
   call object_needed ('xold')
   call object_needed ('grd_psi_over_psi_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_3d_zv1', dens_3d_zv1, grid_xyz_nb)
  call object_alloc ('dens_3d_zv1_av', dens_3d_zv1_av, grid_xyz_nb)
  call object_alloc ('dens_3d_zv1_av_err', dens_3d_zv1_av_err, grid_xyz_nb)

  dens_3d_zv1 (:) = 0.d0

  do elec_i = 1, nelec

      do grid_i = 1, grid_xyz_nb

!     (ri - r)
      di2 = 0.d0
      do dim_i = 1,ndim
         di2 = di2 + (xold(dim_i,elec_i) - grid_xyz (dim_i, grid_i))**2
      enddo ! dim_i
      di = dsqrt(di2)

!     dot product: drift_i . (ri -r)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_i, 1) * (xold(dim_i,elec_i) - grid_xyz (dim_i, grid_i))
      enddo

      dens_3d_zv1 (grid_i) = dens_3d_zv1 (grid_i) - oneover2pi * dotproduct / di**3

     enddo !grid_i

  enddo !elec_i

  end subroutine dens_3d_zv1_bld

! ==============================================================================
  subroutine dens_3d_zv2_bld
! ------------------------------------------------------------------------------
! Description   : second-order renormalized improved estimator of 3D density
!
! Created       : J. Toulouse, 07 Mar 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, dim_i
  real(dp)              :: di, di2, r
  real(dp)              :: dotproduct, dens_temp

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_zv2')

   call object_needed ('grid_xyz_nb')
   call object_needed ('grid_xyz')
   call object_needed ('nelec')
   call object_needed ('lap_psi_over_psi_wlk')
   call object_needed ('grd_psi_over_psi_sq_wlk')
   call object_needed ('xold')

   return

  endif

! allocations
  call object_alloc ('dens_3d_zv2', dens_3d_zv2, grid_xyz_nb)
  call object_alloc ('dens_3d_zv2_av', dens_3d_zv2_av, grid_xyz_nb)
  call object_alloc ('dens_3d_zv2_av_err', dens_3d_zv2_av_err, grid_xyz_nb)

  dens_3d_zv2 (:) = 0.d0

  do elec_i = 1, nelec

      do grid_i = 1, grid_xyz_nb

!     (rij - r)
      di2 = 0.d0
      do dim_i = 1,ndim
         di2 = di2 + (xold (dim_i, elec_i)  - grid_xyz (dim_i, grid_i))**2
      enddo ! dim_i
      di = dsqrt(di2)

      dens_3d_zv2 (grid_i) = dens_3d_zv2 (grid_i)  - oneover2pi * (lap_psi_over_psi_wlk (elec_i, 1) + grd_psi_over_psi_sq_wlk (elec_i, 1)) / di

     enddo !grid_i

  enddo !elec_i

  end subroutine dens_3d_zv2_bld

! ==============================================================================
  subroutine dens_3d_bld
! ------------------------------------------------------------------------------
! Description   : 3D density density  n(r)
!
! Created       : J. Toulouse, 07 Mar 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'dens_3d_bld'

! header
  if (header_exe) then

   call object_create ('dens_3d')
   call object_create ('dens_3d_err')

   call object_needed ('grid_xyz_nb')

   return

  endif

! begin

! allocations
  call object_alloc ('dens_3d', dens_3d, grid_xyz_nb)
  call object_alloc ('dens_3d_err', dens_3d_err, grid_xyz_nb)

  select case(trim(dens_3d_estimator))
   case ('histogram')
   call object_provide (lhere,'dens_3d_histo_av')
   call object_provide (lhere,'dens_3d_histo_av_err')
   dens_3d (:)     = dens_3d_histo_av (:)
   dens_3d_err (:) = dens_3d_histo_av_err (:)

   case ('zv1')
   call object_provide (lhere,'dens_3d_zv1_av')
   call object_provide (lhere,'dens_3d_zv1_av_err')
   dens_3d (:)     = dens_3d_zv1_av (:)
   dens_3d_err (:) = dens_3d_zv1_av_err (:)

   case ('zv2')
   call object_provide (lhere,'dens_3d_zv2_av')
   call object_provide (lhere,'dens_3d_zv2_av_err')
   dens_3d (:)     = dens_3d_zv2_av (:)
   dens_3d_err (:) = dens_3d_zv2_av_err (:)

   case default
   write(6,'(4a)') trim(lhere), ': unknown estimator >',trim(dens_3d_estimator),'<'
   call die (lhere)
  end select

  end subroutine dens_3d_bld

! ========================================================================
  subroutine dens_wrt
! ------------------------------------------------------------------------
! Description    : write spherically averaged density density on file
!
! Created        : J. Toulouse, 07 Nov 2006
! ------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save   :: lhere = 'dens_wrt'
  integer                                    :: unit
  character(len=max_string_len_file)         :: file
  integer                                    :: grid_i

! begin

! return if not main node
# if defined (MPI)
   if (idtask /= 0) return
# endif

  call object_provide ('grid_r_nb')
  call object_provide ('grid_r')
  call object_provide ('dens')


! open file
  file = dens_file_out
  unit = 21
  open(file=trim(file),unit=unit)

  write(unit,'(a)')         'spherically averaged density calculated with CHAMP'
  write(unit,'(a,i5)')      'number of electrons       =',nelec
  write(unit,'(a,i20)')     'number of steps per block =',nstep_total
  write(unit,'(a,i20)')     'number of blocks          =',block_iterations_nb
  write(unit,'(a,3e25.15)') 'grid_r_step               =',grid_r_step
  write(unit,'(a,3e25.15)') 'grid_r_max                =',grid_r_max
  write(unit,'(a,i25)')     'grid_r_nb                 =',grid_r_nb
  write(unit,'(a,f)')       'dist_e_min                =',dist_e_min
  write(unit,'(a,f)')       'dist_e_max                =',dist_e_max


  write(unit,*) ''
  write(unit,'(a)') '              r                        n(r)                error on n(r)'

  do grid_i = 1, grid_r_nb

       write(unit,'(3e25.15)') grid_r (grid_i), dens (grid_i), dens_err (grid_i)
  enddo ! grid_i

  close(unit)

  end subroutine dens_wrt

! ========================================================================
  subroutine dens_3d_wrt
! ------------------------------------------------------------------------
! Description    : write 3D density density on file
!
! Created        : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save   :: lhere = 'dens_3d_wrt'
  integer                                    :: unit
  character(len=max_string_len_file)         :: file
  integer                                    :: grid_i, grid_x_i, grid_y_i, grid_z_i

! begin
  call object_provide ('grid_xyz_nb')
  call object_provide ('grid_xyz')
  call object_provide ('dens_3d')


! open file
  file = dens_3d_file_out
  unit = 21
  open(file=trim(file),unit=unit)

  write(unit,'(a)')         '3D density generated by CHAMP'
  write(unit,'(a,i5)')      'number of electrons       =',nelec
  write(unit,'(a,i20)')     'number of steps per block =',nstep_total
  write(unit,'(a,i20)')     'number of blocks          =',block_iterations_nb
  write(unit,'(a,3e25.15)') 'grid_x_max                =',grid_x_max
  write(unit,'(a,3e25.15)') 'grid_x_min                =',grid_x_min
  write(unit,'(a,3e25.15)') 'grid_x_step               =',grid_x_step
  write(unit,'(a,3e25.15)') 'grid_y_max                =',grid_y_max
  write(unit,'(a,3e25.15)') 'grid_y_min                =',grid_y_min
  write(unit,'(a,3e25.15)') 'grid_x_step               =',grid_y_step
  write(unit,'(a,3e25.15)') 'grid_z_max                =',grid_z_max
  write(unit,'(a,3e25.15)') 'grid_z_min                =',grid_z_min
  write(unit,'(a,3e25.15)') 'grid_z_step               =',grid_z_step
  write(unit,'(a,i25)')     'grid_xyz_nb               =',grid_xyz_nb

  write(unit,*) ''
  write(unit,'(a)') '             x                        y                        z                       n(r)                   error n(r)'

  do grid_x_i = 1, grid_x_nb
   do grid_y_i = 1, grid_y_nb
    do grid_z_i = 1, grid_z_nb

       grid_i = grid_xyz_index (grid_x_i, grid_y_i, grid_z_i)

       write(unit,'(5e25.15)') grid_xyz(1,grid_i), grid_xyz(2,grid_i), grid_xyz(3,grid_i), dens_3d (grid_i), dens_3d_err(grid_i)
    enddo ! grid_z_i
   enddo  ! grid_y_i
       write(unit,'(a)') ''
  enddo ! grid_z_i

  close(unit)

  end subroutine dens_3d_wrt

! ========================================================================
  subroutine dens_3d_wrt_new
! ------------------------------------------------------------------------
! Description    : write 3D density density on file
!
! Created        : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save   :: lhere = 'dens_3d_wrt'
  integer                                    :: unit
  character(len=max_string_len_file)         :: file
  integer                                    :: grid_i, grid_x_i, grid_y_i, grid_z_i

! begin
  call object_provide ('grid_xyz_nb')
  call object_provide ('grid_xyz')
  call object_provide ('dens_3d')


! open file
  file = dens_3d_file_out
  unit = 21
  open(file=trim(file),unit=unit)

  do grid_x_i = 1, grid_x_nb
   do grid_y_i = 1, grid_y_nb
    do grid_z_i = 1, grid_z_nb

       grid_i = grid_xyz_index (grid_x_i, grid_y_i, grid_z_i)

!       write(unit,'(5e25.15)') grid_xyz(1,grid_i), grid_xyz(2,grid_i), grid_xyz(3,grid_i), dens_3d (grid_i), dens_3d_err(grid_i)
       write(unit,'(5e25.15)') grid_xyz(1,grid_i), grid_xyz(3,grid_i), dens_3d (grid_i)
    enddo ! grid_z_i
   enddo  ! grid_y_i
       write(unit,'(a)') ''
  enddo ! grid_z_i

  close(unit)

  end subroutine dens_3d_wrt_new

end module density_mod
