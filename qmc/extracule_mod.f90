module extracule_mod

  use all_tools_mod

! Declaration of global variables and default values
  logical                           :: do_extracule = .false.
  character*(max_string_len_file)   :: file_extracule_out  = ''
  real(dp)                  :: dist_step = 0.d0
  real(dp)                  :: dist_max  = 0.d0
  integer                           :: grid_extra_nb = 0

  integer                           :: extracule_block_nb = 0

  real(dp), allocatable     :: grid_extra(:)
  real(dp), allocatable     :: sum_extracule_improved_num(:)
  real(dp), allocatable     :: extracule_improved_num(:)
  real(dp), allocatable     :: extracule_improved_norm(:)
  real(dp), allocatable     :: sum_extracule_improved_norm(:)
  real(dp), allocatable     :: sum_extracule_improved_norm_square(:)
  real(dp), allocatable     :: extracule_improved_f_av(:)
  real(dp), allocatable     :: extracule_improved_f_var(:)
  real(dp), allocatable     :: extracule_improved_f_av_err(:)
  real(dp), allocatable     :: extracule_improved_f4pir2_av(:)
  real(dp), allocatable     :: extracule_improved_f4pir2_av_err(:)
  integer                           :: extracule_improved_nb
  logical                           :: init_extracule_improved_step = .true.
  logical                           :: init_extracule_improved_block = .true.

  real(dp), allocatable     :: extracule (:)
  real(dp), allocatable     :: extracule_err (:)
  real(dp), allocatable     :: extracule_4pir2 (:)
  real(dp), allocatable     :: extracule_4pir2_err (:)
  real(dp), allocatable     :: extracule_zv1 (:)
  real(dp), allocatable     :: extracule_zv1_av (:)
  real(dp), allocatable     :: extracule_zv1_av_err (:)

  contains
  
!===========================================================================
  subroutine extracule_menu
!---------------------------------------------------------------------------
! Description : menu for extracule calculation
!
! Created     : J. Toulouse, 04 Mar 2006
!---------------------------------------------------------------------------
  implicit none

  character*(max_string_len_rout), save :: lhere = 'extracule_menu'
  character*(max_string_len) line, word

  integer iostat

! begin

! loop over menu lines
  do 
  call get_next_word (word)

  if(trim(word) == 'help') then
   write(6,*) 'extracule-h: menu for extracule'
   write(6,*) 'extracule-h: extracule'
   write(6,*) 'extracule-h:  file    = [string] file in which extracule will be written'
   write(6,*) 'extracule-h:  step    = [real] step on e-e distance'
   write(6,*) 'extracule-h:  distmax = [real] maximum e-e distance'
   write(6,*) 'extracule-h: end'

  elseif(trim(word) == 'file') then
   call get_next_value (file_extracule_out)

  elseif(trim(word) == 'step') then
   call get_next_value (dist_step)

  elseif(trim(word) == 'distmax') then
   call get_next_value (dist_max)

  elseif(trim(word) == 'end') then
   exit

  else

   write(6,*) trim(lhere),'unknown word = ',trim(word)
   call die(here)

  endif
  
  enddo ! end loop over menu lines


! File
  if (trim(file_extracule_out) /= '') then
   write (6,*) trim(lhere),': extracule will be written on file: ',trim(file_extracule_out)
  else
   write (6,*) trim(lhere),': file for writing extracule not specified'
   call die(here)
  endif

! Parameters for grid over (ri+rj)/2
  if(dist_step <= 0.d0) then
   write (6,*) trim(lhere), ': dist_step=',dist_step,' <= 0'
   call die(here)
  endif
  call object_modified ('dist_step')

  if(dist_max <= 0.d0) then
   write (6,*) trim(lhere), ': dist_max=',dist_max,' <= 0'
   call die(here)
  endif

  grid_extra_nb = int(dist_max/dist_step) + 1
  call object_modified ('grid_extra_nb')

  if(grid_extra_nb <= 0) then
    write (6,*) trim(lhere), ': grid_extra_nb=',grid_extra_nb,' <= 0'
    call die(here)
  endif

  write (6,*) trim(lhere), ': Parameters for extracule grid:'
  write (6,*) trim(lhere), ': step = ', dist_step
  write (6,*) trim(lhere), ': dist_max = ', dist_max
  write (6,*) trim(lhere), ': grid_extra_nb=',grid_extra_nb

  call object_average_request ('extracule_zv1_av')
  call object_error_request ('extracule_zv1_av_err')
  call routine_write  ('extracule_wrt')

  end subroutine extracule_menu

! ==============================================================================
  subroutine grid_extra_bld
! ------------------------------------------------------------------------------
! Descriptions   : Grid over rij for calculate of extracule with improved estimator
!
! Created        : J. Toulouse, 06 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  integer                       :: grid_i

! header
  if (header_exe) then

   call object_create ('grid_extra')

   call object_needed ('grid_extra_nb')
   call object_needed ('dist_step')

   return

  endif

! begin
  if(grid_extra_nb <= 0) then
    write(6,*) trim(here),': grid_extra_nb=',grid_extra_nb,' <= 0'
   call die(here)
  endif

  if(dist_step <= 0) then
    write(6,*) trim(here),': dist_step=',dist_step,' <=0'
    call die(here)
  endif

! allocation
  call object_alloc ('grid_extra', grid_extra, grid_extra_nb)

  do grid_i = 1, grid_extra_nb
   grid_extra (grid_i) = (grid_i - 1) * dist_step
  enddo

  end subroutine grid_extra_bld

! ==============================================================================
  subroutine extracule_zv1_bld
! ------------------------------------------------------------------------------
! Description   : first-order renormalized improved estimator of extracule
!
! Created       : J. Toulouse, 05 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: rijcm, rijcm2, r
  real(dp)              :: dotproduct, extracule_temp
  real(dp)              :: oneover2pi

! begin

! header
  if (header_exe) then

   call object_create ('extracule_zv1')

   call object_needed ('grid_extra')
   call object_needed ('grid_extra_nb')
   call object_needed ('nelec')
   call object_needed ('xold')
   call object_needed ('vold')

   return

  endif

! allocations
  call object_alloc ('extracule_zv1', extracule_zv1, grid_extra_nb)
  call object_alloc ('extracule_zv1_av', extracule_zv1_av, grid_extra_nb)
  call object_alloc ('extracule_zv1_av_err', extracule_zv1_av_err, grid_extra_nb)

  oneover2pi = 1.d0/(2.d0*pi)

  extracule_zv1 (:) = 0.d0

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     norm of (ri+rj)/2
      rijcm2 = 0.d0
      do dim_i = 1,ndim
         rijcm2 = rijcm2 + ((xold(dim_i,elec_i) + xold(dim_i,elec_j))/2.d0)**2
      enddo ! dim_i
      rijcm = dsqrt(rijcm2)

!     dot product: drift_j . (rj + ri)/2
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  vold(dim_i,elec_j) * (xold(dim_i,elec_j) + xold(dim_i,elec_i))/2.d0
      end do

      extracule_temp = - oneover2pi * dotproduct / rijcm**3

      do grid_i = 1, grid_extra_nb
        r = grid_extra(grid_i)
        if ( rijcm >= r ) then
           extracule_zv1 (grid_i) = extracule_zv1 (grid_i) + extracule_temp
        endif
     end do !grid_i

    end do !elec_j
  end do !elec_i

  end subroutine extracule_zv1_bld

! ==============================================================================
  subroutine extracule_bld
! ------------------------------------------------------------------------------
! Description   : extracule density  E(r)
!
! Created       : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------------
  implicit none

! header
  if (header_exe) then

   call object_create ('extracule')
   call object_create ('extracule_err')

   call object_needed ('grid_extra_nb')
   call object_needed ('extracule_zv1_av')
   call object_needed ('extracule_zv1_av_err')

   return

  endif

! begin

! allocations
  call object_alloc ('extracule', extracule, grid_extra_nb)
  call object_alloc ('extracule_err', extracule_err, grid_extra_nb)

  extracule (:)     = extracule_zv1_av (:)
  extracule_err (:) = extracule_zv1_av_err (:)

  end subroutine extracule_bld

! ==============================================================================
  subroutine extracule_4pir2_bld
! ------------------------------------------------------------------------------
! Description   : 4 pi r^2 E(r)
!
! Created       : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  integer                       :: grid_i
  real(dp)              :: volume_elt

  
! begin

! header
  if (header_exe) then

   call object_create ('extracule_4pir2')

   call object_needed ('grid_extra_nb')
   call object_needed ('grid_extra')
   call object_needed ('extracule')

   return

  endif

! allocations
  call object_alloc ('extracule_4pir2', extracule_4pir2, grid_extra_nb)
  call object_alloc ('extracule_4pir2_err', extracule_4pir2_err, grid_extra_nb)

! average and statistical error of extracule f(r) 4 pi r^2
  do grid_i = 1, grid_extra_nb
   volume_elt = 4.d0*pi1*grid_extra (grid_i)*grid_extra(grid_i)
   extracule_4pir2 (grid_i)   = extracule (grid_i)*volume_elt
   extracule_4pir2_err (grid_i)  = extracule_err (grid_i)*volume_elt
  enddo

  end subroutine extracule_4pir2_bld

! ========================================================================
  subroutine extracule_wrt
! ------------------------------------------------------------------------
! Description    : write extracule density on file
! Description    : E(r) 4 pi r^2 : extracule with volume element
! Description    : E(r)          : extracule without volume element
!
! Created        : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout),save    :: lhere = 'extracule_wrt'
  integer                                    :: unit
  character(len=max_string_len_file)         :: file
  integer                                    :: grid_i

! begin
  call object_provide ('grid_extra')
  call object_provide ('extracule')
  call object_provide ('extracule_4pir2')


! open file
  file = file_extracule_out
  unit = 20
  open(file=trim(file),unit=unit)

  write(unit,'(a)')         'extracule generated by CHAMP'
  write(unit,'(a,i5)')      'number of electrons       =',nelec
  write(unit,'(a,i20)')      'number of steps per block =',nstep
  write(unit,'(a,i20)')      'number of blocks          =',block_iterations_nb
  write(unit,'(a,3e25.15)')  'distance step             =',dist_step
  write(unit,'(a,3e25.15)')  'max distance              =',dist_max  
  write(unit,'(a,i25)')     'number of grid points     =',grid_extra_nb

  write(unit,*) ''
  write(unit,'(a)') '             r                   4 pi r^2 E(r)          error 4pi r^2 E(r)             E(r)                  error E(r)'

  do grid_i = 1, grid_extra_nb

  write(unit,'(3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,3e25.15)') grid_extra(grid_i), &
        extracule_4pir2(grid_i),extracule_4pir2_err(grid_i),extracule(grid_i), extracule_err(grid_i)
  enddo


  close(unit)

  end subroutine extracule_wrt

end module extracule_mod
