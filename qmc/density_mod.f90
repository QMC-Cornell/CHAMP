module density_mod

  use all_tools_mod
  use grid_mod
  use electrons_mod
  use psi_mod
  use eloc_mod
  use deriv_mod

! Declaration of global variables and default values
! For spherically symmetric densities.  Estimators, histogram, zv1, zv5, zvzb1, zvzb5 are implemented.
  character(len=max_string_len_file):: dens_file_out  = ''
  character(len=max_string_len)     :: dens_estimator  = 'zv1'
  real(dp)                          :: dens_exp = 1.d0
  real(dp), allocatable             :: dens_histo (:,:)
  real(dp), allocatable             :: dens_histo_av (:)
  real(dp), allocatable             :: dens_histo_av_err (:)
  real(dp), allocatable             :: dens_zv1 (:,:)
  real(dp), allocatable             :: dens_zv1_av (:)
  real(dp), allocatable             :: dens_zv1_av_err (:)
  real(dp), allocatable             :: dens_zv5 (:,:)
  real(dp), allocatable             :: dens_zv5_av (:)
  real(dp), allocatable             :: dens_zv5_av_err (:)
  real(dp), allocatable             :: dens_q1 (:,:)
  real(dp), allocatable             :: dens_q1_av (:)
  real(dp), allocatable             :: dens_q1_eloc (:,:)
  real(dp), allocatable             :: dens_q1_eloc_av (:)
  real(dp), allocatable             :: dens_zb1_av (:)
  real(dp), allocatable             :: dens_zvzb1_av (:)
  real(dp), allocatable             :: dens_zvzb1_av_err (:)
  real(dp), allocatable             :: dens_q5 (:,:)
  real(dp), allocatable             :: dens_q5_av (:)
  real(dp), allocatable             :: dens_q5_eloc (:,:)
  real(dp), allocatable             :: dens_q5_eloc_av (:)
  real(dp), allocatable             :: dens_zb5_av (:)
  real(dp), allocatable             :: dens_zvzb5_av (:)
  real(dp), allocatable             :: dens_zvzb5_av_err (:)
  real(dp), allocatable             :: dens (:)
  real(dp), allocatable             :: dens_err (:)

! For axially symmetric densities.  Estimators: at present only histogram is implemented.
  character(len=max_string_len_file):: dens_xy_z_file_out  = ''
  character(len=max_string_len)     :: dens_xy_z_estimator  = 'histogram'
  character(len=max_string_len)     :: point_group  = 'c_inf_v'
  real(dp), allocatable             :: dens_xy_z_histo (:,:)
  real(dp), allocatable             :: dens_xy_z_histo_av (:)
  real(dp), allocatable             :: dens_xy_z_histo_av_err (:)
  real(dp), allocatable             :: dens_xy_z (:)
  real(dp), allocatable             :: dens_xy_z_err (:)

! Does not assume any symmetry for the density.  Estimators, histogram, zv1, zv2 are implemented.
  character(len=max_string_len_file):: dens_3d_file_out  = ''
  character(len=max_string_len)     :: dens_3d_estimator  = 'zv2'
  real(dp), allocatable             :: dens_3d_histo (:,:)
  real(dp), allocatable             :: dens_3d_histo_av (:)
  real(dp), allocatable             :: dens_3d_histo_av_err (:)
  real(dp), allocatable             :: dens_3d_zv1 (:,:)
  real(dp), allocatable             :: dens_3d_zv1_av (:)
  real(dp), allocatable             :: dens_3d_zv1_av_err (:)
  real(dp), allocatable             :: dens_3d_zv2 (:,:)
  real(dp), allocatable             :: dens_3d_zv2_av (:)
  real(dp), allocatable             :: dens_3d_zv2_av_var (:)
  real(dp), allocatable             :: dens_3d_zv2_av_err (:)
  real(dp), allocatable             :: dens_3d_zv2_linear_av (:)
  real(dp), allocatable             :: dens_3d_zv2_linear_av_var (:)
  real(dp), allocatable             :: dens_3d_zv2_linear_av_err (:)
  real(dp), allocatable             :: dens_3d_zv2_linear_coef (:,:)
  real(dp), allocatable             :: dens_3d_zv2_deloc_covar (:,:)
  real(dp), allocatable             :: dens_3d_zv2_deloc (:,:)
  real(dp), allocatable             :: dens_3d_zv2_deloc_av (:,:)
  real(dp), allocatable             :: dens_3d_zv2_av_deloc_av_covar (:,:)
  real(dp), allocatable             :: dens_3d (:)
  real(dp), allocatable             :: dens_3d_err (:)

! For periodic densities in fourier space
  character(len=max_string_len_file):: dens_fourier_file_out  = ''
  complex(dpc), allocatable         :: dens_fourier (:)
  complex(dpc), allocatable         :: dens_fourier_av (:)
  complex(dpc), allocatable         :: dens_fourier_av_err (:)

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
  write(6,*)
  write(6,'(a)') 'Beginning of density menu --------------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)') 'HELP for menu density'
   write(6,'(a)') 'density'
   write(6,'(a)') ' estimator = [string] : choice of estimator {histogram|zv1|zv5|zvzb1|zvzb5} (default=zv1)'
   write(6,'(a)') ' exponent  = [real]   : exponent for zv5 or zvzb5 estimator (should be 2*sqrt(2*I)) (default = 1.0)'
   write(6,'(a)') ' file      = [string] : file in which density will be written'
   write(6,'(a)') 'end'

  case ('file')
   call get_next_value (dens_file_out)

  case ('estimator')
   call get_next_value (dens_estimator)

  case ('exponent')
   call get_next_value (dens_exp)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines


! File
  if (trim(dens_file_out) /= '') then
   write (6,'(3a)') ' density will be written on file >',trim(dens_file_out),'<.'
  else
   call die (lhere, 'file for writing density not specified.')
  endif

  select case(trim(dens_estimator))
   case ('histogram')
   write (6,'(a)') ' density will be calculated with histogram estimator'
   call object_average_request ('dens_histo_av')
   call object_error_request ('dens_histo_av_err')

   case ('zv1')
   write (6,'(a)') ' density will be calculated with ZV1 estimator (simplest improved estimator)'
   call object_average_request ('dens_zv1_av')
   call object_error_request ('dens_zv1_av_err')

   case ('zv2')
   call die (lhere, 'estimator "zv2" not yet implemented.')
   call object_average_request ('dens_zv2_av')
   call object_error_request ('dens_zv2_av_err')

   case ('zv5')
   write (6,'(a,es15.8)') ' density will be calculated with ZV5 estimator, i.e. improved estimator with exponential decay with exponent=',dens_exp
   call object_average_request ('dens_zv5_av')
   call object_error_request ('dens_zv5_av_err')

   case ('zvzb1')
   write (6,'(a)') ' density will be calculated with ZVZB1 estimator'
   call object_average_request ('dens_zv1_av')
   call object_average_request ('dens_q1_av')
   call object_average_request ('dens_q1_eloc_av')
   call object_error_request ('dens_zvzb1_av_err')

   case ('zvzb5')
   write (6,'(a)') ' density will be calculated with ZVZB5 estimator'
   call object_average_request ('dens_zv5_av')
   call object_average_request ('dens_q5_av')
   call object_average_request ('dens_q5_eloc_av')
   call object_error_request ('dens_zvzb5_av_err')

   case default
   call die (lhere, 'unknown estimator >'+trim(dens_estimator)+'<')
  end select

  call routine_write_block_request  ('dens_wrt')

  write(6,'(a)') 'End of density menu --------------------------------------------------------------------------------------'

  end subroutine dens_menu

!===========================================================================
  subroutine dens_xy_z_menu
!---------------------------------------------------------------------------
! Description : menu for axially averaged density calculation
!
! Created     : C. Umrigar, 07 Jul 2012, by imitating dens_menu of J. Toulouse
!---------------------------------------------------------------------------
  implicit none

  character(len=max_string_len_rout), save :: lhere = 'dens_xy_z_menu'

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of dens_xy_z menu --------------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)') 'HELP for menu density'
   write(6,'(a)') 'density'
   write(6,'(a)') ' exponent  = [real] : exponent for zv5 or zvzb5 estimator (should be 2*sqrt(2*I)) (default = 1.0)'
   write(6,'(a)') ' file      = [string] : file in which density will be written'
   write(6,'(a)') 'end'

  case ('file')
   call get_next_value (dens_xy_z_file_out)

  case ('estimator')
   call get_next_value (dens_xy_z_estimator)

  case ('point_group')
   call get_next_value (point_group)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines


! File
  if (trim(dens_xy_z_file_out) /= '') then
   write (6,'(3a)') ' density will be written on file >',trim(dens_xy_z_file_out),'<.'
  else
   call die (lhere, 'file for writing density not specified.')
  endif

  select case(trim(dens_xy_z_estimator))
   case ('histogram')
   write (6,'(a)') ' density will be calculated with histogram estimator'
   call object_average_request ('dens_xy_z_histo_av')
   call object_error_request ('dens_xy_z_histo_av_err')

   case default
   call die (lhere, 'unknown estimator >'+trim(dens_xy_z_estimator)+'<')
  end select

  select case(trim(point_group))
   case ('c_inf_v')
   write (6,'(a)') ' density has c_inf_v symmetry'

   case ('d_inf_h')
   write (6,'(a)') ' density has d_inf_h symmetry'

   case default
   call die (lhere, 'unknown point_group >'+trim(point_group)+'<')
  end select


  call routine_write_block_request  ('dens_xy_z_wrt')

  write(6,'(a)') 'End of density menu --------------------------------------------------------------------------------------'

  end subroutine dens_xy_z_menu

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
  write(6,*)
  write(6,'(a)') 'Beginning of density_3d menu -----------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)') 'HELP for density_3d menu'
   write(6,'(a)') 'density_3d'
   write(6,'(a)') '  estimator = [string] : can be "histogram", "zv1",  "zv2" or "zv2_linear" (default=zv2)'
   write(6,'(a)') '  file    = [string]   : file in which density will be written'
   write(6,'(a)') 'end'

  case ('file')
   call get_next_value (dens_3d_file_out)

  case ('estimator')
   call get_next_value (dens_3d_estimator)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines


! File
  if (trim(dens_3d_file_out) /= '') then
   write (6,'(3a)') ' density will be written on file >',trim(dens_3d_file_out),'<.'
  else
   call die (lhere, 'file for writing 3D density not specified.')
  endif

  select case(trim(dens_3d_estimator))
   case ('histogram')
   write(6,'(a)') ' density will be calculated with histogram estimator.'
   call object_average_request ('dens_3d_histo_av')
   call object_error_request ('dens_3d_histo_av_err')

   case ('zv1')
   write(6,'(a)') ' density will be calculated with the zero-variance estimator using the drift velocity.'
   call object_average_request ('dens_3d_zv1_av')
   call object_error_request ('dens_3d_zv1_av_err')

   case ('zv2')
   write(6,'(a)') ' density will be calculated with the zero-variance estimator using the drift velocity and the laplacian of the wave function.'
   call object_average_request ('dens_3d_zv2_av')
   call object_error_request ('dens_3d_zv2_av_err')

   case ('zv2_linear')
   write(6,'(a)') ' density will be calculated with the zero-variance estimator using the drift velocity, the laplacian of the wave function and the wave function derivatives wrt parameters.'
   call object_average_request ('dens_3d_zv2_av')
   call object_error_request ('dens_3d_zv2_av_err')
   call object_average_request ('deloc_av')
   call object_average_request ('deloc_deloc_av')
   call object_average_request ('dens_3d_zv2_deloc_av')
   call object_covariance_request ('deloc_av_deloc_av_covar')
   call object_covariance_request ('dens_3d_zv2_av_deloc_av_covar')
   call object_error_request ('dens_3d_zv2_linear_av_err')

   case default
   call die (lhere, 'unknown estimator >'+trim(dens_3d_estimator)+'<.')
  end select

  call routine_write_block_request ('dens_3d_wrt')

  write(6,'(a)') 'End of density_3d menu -----------------------------------------------------------------------------------'

  end subroutine dens_3d_menu

!===========================================================================
  subroutine dens_fourier_menu
!---------------------------------------------------------------------------
! Description : menu for density calculation for periodic calculations (in Fourier space)
!
! Created     : J. Toulouse, 28 Sep 2013
!---------------------------------------------------------------------------
  implicit none

  character(len=max_string_len_rout), save :: lhere = 'dens_fourier_menu'

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of density_fourier menu -----------------------------------------------------------------------------'

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)') 'HELP for density_fourier menu'
   write(6,'(a)') 'density_fourier'
   write(6,'(a)') '  file    = [string]   : file in which density will be written'
   write(6,'(a)') 'end'

  case ('file')
   call get_next_value (dens_fourier_file_out)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines


! File
  if (trim(dens_fourier_file_out) /= '') then
   write (6,'(3a)') ' density will be written on file >',trim(dens_fourier_file_out),'<.'
  else
   call die (lhere, 'file for writing periodic density not specified.')
  endif

  call object_average_request ('dens_fourier_av')
  call object_error_request ('dens_fourier_av_err')
  call routine_write_block_request ('dens_fourier_wrt')

  write(6,'(a)') 'End of density_fourier menu -----------------------------------------------------------------------------------'

  end subroutine dens_fourier_menu

! ==============================================================================
  subroutine dens_histo_bld
! ------------------------------------------------------------------------------
! Description   : histogram estimator of radial density
!
! Created       : J. Toulouse, 02 Oct 2008
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer               :: grid_i
  integer               :: elec_i, walk_i
! real(dp)              :: di, d4pir2
  real(dp)              :: di, factor

! begin

! header
  if (header_exe) then

   call object_create ('dens_histo')
   call object_average_walk_define ('dens_histo', 'dens_histo_av')
   call object_error_define ('dens_histo_av', 'dens_histo_av_err')

   call object_needed ('grid_r')
   call object_needed ('grid_r_step')
   call object_needed ('grid_r_nb')
   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('dist_e_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_histo', dens_histo, grid_r_nb, nwalk)
  call object_alloc ('dens_histo_av', dens_histo_av, grid_r_nb)
  call object_alloc ('dens_histo_av_err', dens_histo_av_err, grid_r_nb)

  dens_histo (:,:) = 0.d0

  do walk_i = 1, nwalk
    do elec_i = 1, nelec

!     distance e-origin
      di = dist_e_wlk (elec_i, walk_i)

      grid_i = floor(di/grid_r_step - 0.5d0) + 2

      if (grid_i <= grid_r_nb) then
!      d4pir2 = 4.d0*pi*grid_r(grid_i)**2
       if(grid_r(grid_i).ne.0.d0) then
         factor = 4.d0*pi*grid_r(grid_i)**2 * grid_r_step     ! 4 pi r^2 deltar
       else
         factor = 1.333333333333d0*pi*(0.5d0*grid_r_step)**3  ! (4/3) pi (deltar/2)^3
       endif
       dens_histo (grid_i, walk_i) = dens_histo (grid_i, walk_i) + 1.d0/factor
      endif

    enddo ! elec_i
  enddo ! walk_i

  end subroutine dens_histo_bld

! ==============================================================================
  subroutine dens_zv1_bld
! ------------------------------------------------------------------------------
! Description   : first-order renormalized improved estimator of spherically averaged density
!
! Created       : J. Toulouse, 07 Nov 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer grid_i, elec_i, dim_i, walk_i
  real(dp) di, r, dotproduct, dens_temp

! begin

! header
  if (header_exe) then

   call object_create ('dens_zv1')
   call object_average_walk_define ('dens_zv1', 'dens_zv1_av')
   call object_error_define ('dens_zv1_av', 'dens_zv1_av_err')

   call object_needed ('nwalk')
   call object_needed ('grid_r_nb')
   call object_needed ('grid_r')
   call object_needed ('nelec')
   call object_needed ('coord_elec_wlk')
   call object_needed ('dist_e_wlk')
   call object_needed ('grd_psi_over_psi_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_zv1', dens_zv1, grid_r_nb, nwalk)
  call object_alloc ('dens_zv1_av', dens_zv1_av, grid_r_nb)
  call object_alloc ('dens_zv1_av_err', dens_zv1_av_err, grid_r_nb)

  dens_zv1 (:,:) = 0.d0

  do walk_i = 1, nwalk
    do elec_i = 1, nelec

!     distance |r_i|
      di = dist_e_wlk (elec_i, walk_i)

!     dot product: drift_i . r_i
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_i, walk_i) * coord_elec_wlk (dim_i, elec_i, walk_i)
      enddo

      dens_temp = - oneover2pi * dotproduct / di**3

       do grid_i = 1, grid_r_nb
        r = grid_r (grid_i)

        if ( di >= r ) then
           dens_zv1 (grid_i, walk_i) = dens_zv1 (grid_i, walk_i) + dens_temp
        endif

     enddo ! grid_i

    enddo ! elec_i
  enddo ! walk_i

  end subroutine dens_zv1_bld

! ==============================================================================
  subroutine dens_zv5_bld
! ------------------------------------------------------------------------------
! Description   : improved zero-variance estimator of radial density
! Description   : Q(R) = (-1/4pi) sum_i dOmega_u/(4pi) 1/|ri -u| f(ri,u)
! Description   : with f(ri,u) =exp [-c |ri-u|]
!
! Created       : J. Toulouse, 02 Oct 2008
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer               :: grid_i
  integer               :: elec_i, dim_i, walk_i
  real(dp)              :: cc, cc2
  real(dp)              :: r, r2, u, rpu, ru, rmu, abs_rmu
  real(dp)              :: mcabsrmu, mcrpu, exp_mcabsrmu, exp_mcrpu
  real(dp)              :: dotproduct
  real(dp)              :: dens_temp1, dens_temp2

! begin

! header
  if (header_exe) then

   call object_create ('dens_zv5')
   call object_average_walk_define ('dens_zv5', 'dens_zv5_av')
   call object_error_define ('dens_zv5_av', 'dens_zv5_av_err')

   call object_needed ('nwalk')
   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_e_wlk')
   call object_needed ('coord_elec_wlk')
   call object_needed ('grd_psi_over_psi_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_zv5', dens_zv5, grid_r_nb, nwalk)
  call object_alloc ('dens_zv5_av', dens_zv5_av, grid_r_nb)
  call object_alloc ('dens_zv5_av_err', dens_zv5_av_err, grid_r_nb)

  dens_zv5 (:,:) = 0.d0

  cc = dens_exp
  cc2 = cc*cc

  do walk_i = 1, nwalk
    do elec_i = 1, nelec

!     distance e-origin
      r = dist_e_wlk (elec_i, walk_i)
      r2 = r*r
!
!     dot product: drift_i . ri
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk  (dim_i, elec_i, walk_i) * coord_elec_wlk (dim_i, elec_i, walk_i)
      enddo

!     loop over grid points
      do grid_i = 1, grid_r_nb

        u = grid_r (grid_i)

        rpu = r + u
        rmu   = r - u
        abs_rmu = dabs(rmu)
        ru = r * u

        mcabsrmu = -cc*abs_rmu
        mcrpu = -cc*rpu

        exp_mcabsrmu = dexp(mcabsrmu)
        exp_mcrpu = dexp(mcrpu)

        dens_temp1 = (exp_mcabsrmu - exp_mcrpu)/(2.d0*cc*ru)
        dens_temp2 = 0.5d0 * (sign(1.d0,rmu) * exp_mcabsrmu - exp_mcrpu)

        dens_zv5 (grid_i, walk_i) = dens_zv5 (grid_i, walk_i) -          &
            oneover2pi*((dotproduct/r2)*(dens_temp1 + dens_temp2/u)  &
            - 0.5d0*cc2*dens_temp1)

     enddo !grid_i

   enddo !elec_i
  enddo !walk_i

  end subroutine dens_zv5_bld

! ==============================================================================
  subroutine dens_q1_bld
! ------------------------------------------------------------------------------
! Description   : term for zero-bias correction for density
!
! Created       : J. Toulouse, 04 Jan 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer grid_i, elec_i, walk_i
  real(dp) di, r

! begin

! header
  if (header_exe) then

   call object_create ('dens_q1')
   call object_average_walk_define ('dens_q1', 'dens_q1_av')

   call object_needed ('nwalk')
   call object_needed ('grid_r_nb')
   call object_needed ('grid_r')
   call object_needed ('nelec')
   call object_needed ('dist_e_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_q1', dens_q1, grid_r_nb, nwalk)
  call object_alloc ('dens_q1_av', dens_q1_av, grid_r_nb)

  dens_q1 (:,:) = 0.d0

  do walk_i = 1, nwalk
    do elec_i = 1, nelec

!     distance |r_i|
      di = dist_e_wlk (elec_i, walk_i)

       do grid_i = 1, grid_r_nb
        r = grid_r (grid_i)

        if ( di >= r ) then
           dens_q1 (grid_i, walk_i) = dens_q1 (grid_i, walk_i) + (1.d0/di)
        else
           dens_q1 (grid_i, walk_i) = dens_q1 (grid_i, walk_i) + (1.d0/r)
        endif

     enddo ! grid_i

    enddo ! elec_i
  enddo ! walk_i

  end subroutine dens_q1_bld

! ==============================================================================
  subroutine dens_q1_eloc_bld
! ------------------------------------------------------------------------------
! Description   : term for zero-bias correction for density
!
! Created       : J. Toulouse, 04 Jan 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer walk_i

! begin

! header
  if (header_exe) then

   call object_create ('dens_q1_eloc')
   call object_average_walk_define ('dens_q1_eloc', 'dens_q1_eloc_av')

   call object_needed ('nwalk')
   call object_needed ('grid_r_nb')
   call object_needed ('dens_q1')
   call object_needed ('eloc_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_q1_eloc', dens_q1_eloc, grid_r_nb, nwalk)
  call object_alloc ('dens_q1_eloc_av', dens_q1_eloc_av, grid_r_nb)

  do walk_i = 1, nwalk
    dens_q1_eloc (:, walk_i) = dens_q1 (:, walk_i) * eloc_wlk (walk_i)
  enddo

  end subroutine dens_q1_eloc_bld

! ==============================================================================
  subroutine dens_zb1_av_bld
! ------------------------------------------------------------------------------
! Description   : zero-bias correction for density
!
! Created       : J. Toulouse, 04 Jan 2010
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('dens_zb1_av')

   call object_needed ('grid_r_nb')
   call object_needed ('dens_q1_eloc_av')
   call object_needed ('dens_q1_av')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('dens_zb1_av', dens_zb1_av, grid_r_nb)

  dens_zb1_av (:) = -oneover2pi * (dens_q1_eloc_av (:) - dens_q1_av (:) * eloc_av)

  end subroutine dens_zb1_av_bld

! ==============================================================================
  subroutine dens_zvzb1_av_bld
! ------------------------------------------------------------------------------
! Description   : zero-variance zero-bias estimator for density
!
! Created       : J. Toulouse, 04 Jan 2010
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('dens_zvzb1_av')
   call object_error_define ('dens_zvzb1_av', 'dens_zvzb1_av_err')

   call object_needed ('grid_r_nb')
   call object_needed ('dens_zv1_av')
   call object_needed ('dens_zb1_av')

   return

  endif

! allocations
  call object_alloc ('dens_zvzb1_av', dens_zvzb1_av, grid_r_nb)
  call object_alloc ('dens_zvzb1_av_err', dens_zvzb1_av_err, grid_r_nb)

  dens_zvzb1_av (:) = dens_zv1_av (:) + dens_zb1_av (:)

  end subroutine dens_zvzb1_av_bld

! ==============================================================================
  subroutine dens_q5_bld
! ------------------------------------------------------------------------------
! Description   : term for zero-bias correction for density
!
! Created       : J. Toulouse, 07 Jan 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer               :: walk_i, elec_i, grid_i
  real(dp)              :: r
  real(dp)              :: u, rpu, ru, rmu, abs_rmu
  real(dp)              :: mcabsrmu, mcrpu, exp_mcabsrmu, exp_mcrpu

! begin

! header
  if (header_exe) then

   call object_create ('dens_q5')
   call object_average_walk_define ('dens_q5', 'dens_q5_av')

   call object_needed ('nwalk')
   call object_needed ('grid_r_nb')
   call object_needed ('grid_r')
   call object_needed ('nelec')
   call object_needed ('dist_e_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_q5', dens_q5, grid_r_nb, nwalk)
  call object_alloc ('dens_q5_av', dens_q5_av, grid_r_nb)

  dens_q5 (:,:) = 0.d0

  do walk_i = 1, nwalk
    do elec_i = 1, nelec

!     distance |r_i|
      r = dist_e_wlk (elec_i, walk_i)

!     loop over grid points
      do grid_i = 1, grid_r_nb

        u = grid_r (grid_i)

        rpu = r + u
        rmu   = r - u
        abs_rmu = dabs(rmu)
        ru = r * u

        mcabsrmu = -dens_exp*abs_rmu
        mcrpu = -dens_exp*rpu

        exp_mcabsrmu = dexp(mcabsrmu)
        exp_mcrpu = dexp(mcrpu)

        dens_q5 (grid_i, walk_i) = dens_q5 (grid_i, walk_i) + (exp_mcabsrmu - exp_mcrpu)/(2.d0*dens_exp*ru)

     enddo ! grid_i

    enddo ! elec_i
  enddo ! walk_i

  end subroutine dens_q5_bld

! ==============================================================================
  subroutine dens_q5_eloc_bld
! ------------------------------------------------------------------------------
! Description   : term for zero-bias correction for density
!
! Created       : J. Toulouse, 07 Jan 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer walk_i

! begin

! header
  if (header_exe) then

   call object_create ('dens_q5_eloc')
   call object_average_walk_define ('dens_q5_eloc', 'dens_q5_eloc_av')

   call object_needed ('nwalk')
   call object_needed ('grid_r_nb')
   call object_needed ('dens_q5')
   call object_needed ('eloc_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_q5_eloc', dens_q5_eloc, grid_r_nb, nwalk)
  call object_alloc ('dens_q5_eloc_av', dens_q5_eloc_av, grid_r_nb)

  do walk_i = 1, nwalk
    dens_q5_eloc (:, walk_i) = dens_q5 (:, walk_i) * eloc_wlk (walk_i)
  enddo

  end subroutine dens_q5_eloc_bld

! ==============================================================================
  subroutine dens_zb5_av_bld
! ------------------------------------------------------------------------------
! Description   : zero-bias correction for density
!
! Created       : J. Toulouse, 07 Jan 2010
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('dens_zb5_av')

   call object_needed ('grid_r_nb')
   call object_needed ('dens_q5_eloc_av')
   call object_needed ('dens_q5_av')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('dens_zb5_av', dens_zb5_av, grid_r_nb)

  dens_zb5_av (:) = -oneover2pi * (dens_q5_eloc_av (:) - dens_q5_av (:) * eloc_av)

  end subroutine dens_zb5_av_bld

! ==============================================================================
  subroutine dens_zvzb5_av_bld
! ------------------------------------------------------------------------------
! Description   : zero-variance zero-bias estimator for density
!
! Created       : J. Toulouse, 07 Jan 2010
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('dens_zvzb5_av')
   call object_error_define ('dens_zvzb5_av', 'dens_zvzb5_av_err')

   call object_needed ('grid_r_nb')
   call object_needed ('dens_zv5_av')
   call object_needed ('dens_zb5_av')

   return

  endif

! allocations
  call object_alloc ('dens_zvzb5_av', dens_zvzb5_av, grid_r_nb)
  call object_alloc ('dens_zvzb5_av_err', dens_zvzb5_av_err, grid_r_nb)

  dens_zvzb5_av (:) = dens_zv5_av (:) + dens_zb5_av (:)

  end subroutine dens_zvzb5_av_bld

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
   case ('histogram')
   call object_provide (lhere,'dens_histo_av')
   call object_provide (lhere,'dens_histo_av_err')
   dens (:)     = dens_histo_av (:)
   dens_err (:) = dens_histo_av_err (:)

   case ('zv1')
   call object_provide (lhere,'dens_zv1_av')
   call object_provide (lhere,'dens_zv1_av_err')
   dens (:)     = dens_zv1_av (:)
   dens_err (:) = dens_zv1_av_err (:)

   case ('zv5')
   call object_provide (lhere,'dens_zv5_av')
   call object_provide (lhere,'dens_zv5_av_err')
   dens (:)     = dens_zv5_av (:)
   dens_err (:) = dens_zv5_av_err (:)

   case ('zvzb1')
   call object_provide (lhere,'dens_zvzb1_av')
   call object_provide (lhere,'dens_zvzb1_av_err')
   dens (:)     = dens_zvzb1_av (:)
   dens_err (:) = dens_zvzb1_av_err (:)

   case ('zvzb5')
   call object_provide (lhere,'dens_zvzb5_av')
   call object_provide (lhere,'dens_zvzb5_av_err')
   dens (:)     = dens_zvzb5_av (:)
   dens_err (:) = dens_zvzb5_av_err (:)

   case default
   call die (lhere, 'unknown estimator >'+trim(dens_estimator)+'<.')
  end select

  end subroutine dens_bld

! ==============================================================================
  subroutine dens_xy_z_histo_bld
! ------------------------------------------------------------------------------
! Description   : histo estimator of axial density
!
! Created       : C. Umrigar, 07 Jul 2012, by imitating dens_histo_blk of J. Toulouse
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer grid_i, grid_xy_i, grid_z_i
  integer elec_i, walk_i
! real(dp) xyi, zi, d2pixy
  real(dp) xyi, zi, factor

! begin

! header
  if (header_exe) then

   call object_create ('dens_xy_z_histo')
   call object_average_walk_define ('dens_xy_z_histo', 'dens_xy_z_histo_av')
   call object_error_define ('dens_xy_z_histo_av', 'dens_xy_z_histo_av_err')

   call object_needed ('nwalk')
   call object_needed ('grid_xy_z_nb')
   call object_needed ('grid_xy_z')
   call object_needed ('grid_xy_z_index')
   call object_needed ('grid_xy_nb')
   call object_needed ('grid_z_nb')
   call object_needed ('grid_xy_step')
   call object_needed ('grid_z_step')
   call object_needed ('grid_xy_min')
   call object_needed ('grid_z_min')
   call object_needed ('nelec')
   call object_needed ('coord_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_xy_z_histo', dens_xy_z_histo, grid_xy_z_nb, nwalk)
  call object_alloc ('dens_xy_z_histo_av', dens_xy_z_histo_av, grid_xy_z_nb)
  call object_alloc ('dens_xy_z_histo_av_err', dens_xy_z_histo_av_err, grid_xy_z_nb)

  dens_xy_z_histo (:,:) = 0.d0

  do walk_i = 1, nwalk
    do elec_i = 1, nelec

      xyi = sqrt(coord_elec_wlk(1,elec_i,walk_i)**2+coord_elec_wlk(2,elec_i,walk_i)**2)
      zi =  coord_elec_wlk(3,elec_i,walk_i)

!     grid_xy_i = floor((xyi - grid_xy_min)/grid_xy_step - 0.5d0)  + 2
      grid_xy_i = nint((xyi - grid_xy_min)/grid_xy_step) + 1
      if(point_group.eq.'c_inf_v') then
!       grid_z_i = floor((zi - grid_z_min)/grid_z_step - 0.5d0)  + 2
        grid_z_i = nint((zi - grid_z_min)/grid_z_step) + 1
      elseif(point_group.eq.'d_inf_h') then
!       grid_z_i = floor((abs(zi) - grid_z_min)/grid_z_step - 0.5d0)  + 2
        grid_z_i = nint((abs(zi) - grid_z_min)/grid_z_step) + 1
      endif

      if (grid_xy_i < 1 .or. grid_xy_i > grid_xy_nb) cycle
      if (grid_z_i < 1 .or. grid_z_i > grid_z_nb) cycle

      grid_i = grid_xy_z_index (grid_xy_i, grid_z_i)
!     d2pixy = 2.d0*pi*grid_xy_z(1,grid_i)

!     if(grid_xy_z(1,grid_i).ne.0.d0) then
!       if(point_group.eq.'c_inf_v') then
!         factor = 2.d0*pi*grid_xy_z(1,grid_i) * grid_xy_step * grid_z_step  ! 2 pi r deltaxy deltaz
!       elseif(point_group.eq.'d_inf_h') then
!         factor = 4.d0*pi*grid_xy_z(1,grid_i) * grid_xy_step * grid_z_step  ! 2 pi r deltaxy deltaz
!       endif
!     else
!       factor = pi*(0.5d0*grid_xy_step)**2 * grid_z_step                  ! pi (deltaxy/2)^2 deltaz
!     endif

      if(grid_xy_z(1,grid_i).ne.0.d0) then
        factor = 2.d0*pi*grid_xy_z(1,grid_i) * grid_xy_step * grid_z_step  ! 2 pi r deltaxy deltaz
      else
        factor = pi*(0.5d0*grid_xy_step)**2 * grid_z_step                  ! pi (deltaxy/2)^2 deltaz
      endif

      if(grid_xy_z(2,grid_i).ne.0.d0 .and. point_group.eq.'d_inf_h') then  ! Since we are adding density at z and -z
        factor=2*factor
      endif

!     dens_xy_z_histo (grid_i, walk_i) = dens_xy_z_histo (grid_i, walk_i) + 1.d0/(d2pixy * grid_xy_step * grid_z_step)
      dens_xy_z_histo (grid_i, walk_i) = dens_xy_z_histo (grid_i, walk_i) + 1.d0/factor

    enddo ! elec_i
  enddo ! walk_i

  end subroutine dens_xy_z_histo_bld

! ==============================================================================
  subroutine dens_xy_z_bld
! ------------------------------------------------------------------------------
! Description   : axially averaged density density  n(r)
!
! Created       : C. Umrigar, 07 Jul 2012, by imitating dens_bld of J. Toulouse
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'dens_xy_z_bld'

! header
  if (header_exe) then

   call object_create ('dens_xy_z')
   call object_create ('dens_xy_z_err')

   call object_needed ('grid_xy_z_nb')

   return

  endif

! begin

! allocations
  call object_alloc ('dens_xy_z', dens_xy_z, grid_xy_z_nb)
  call object_alloc ('dens_xy_z_err', dens_xy_z_err, grid_xy_z_nb)

  select case(trim(dens_xy_z_estimator))
   case ('histogram')
   call object_provide (lhere,'dens_xy_z_histo_av')
   call object_provide (lhere,'dens_xy_z_histo_av_err')
   dens_xy_z (:)     = dens_xy_z_histo_av (:)
   dens_xy_z_err (:) = dens_xy_z_histo_av_err (:)

!  case ('zv1')
!  call object_provide (lhere,'dens_xy_z_zv1_av')
!  call object_provide (lhere,'dens_xy_z_zv1_av_err')
!  dens (:)     = dens_xy_z_zv1_av (:)
!  dens_xy_z_err (:) = dens_xy_z_zv1_av_err (:)

!  case ('zv5')
!  call object_provide (lhere,'dens_xy_z_zv5_av')
!  call object_provide (lhere,'dens_xy_z_zv5_av_err')
!  dens (:)     = dens_xy_z_zv5_av (:)
!  dens_xy_z_err (:) = dens_xy_z_zv5_av_err (:)

!  case ('zvzb1')
!  call object_provide (lhere,'dens_xy_z_zvzb1_av')
!  call object_provide (lhere,'dens_xy_z_zvzb1_av_err')
!  dens (:)     = dens_xy_z_zvzb1_av (:)
!  dens_xy_z_err (:) = dens_xy_z_zvzb1_av_err (:)

!  case ('zvzb5')
!  call object_provide (lhere,'dens_xy_z_zvzb5_av')
!  call object_provide (lhere,'dens_xy_z_zvzb5_av_err')
!  dens (:)     = dens_xy_z_zvzb5_av (:)
!  dens_xy_z_err (:) = dens_xy_z_zvzb5_av_err (:)

   case default
   call die (lhere, 'unknown estimator >'+trim(dens_xy_z_estimator)+'<.')
  end select

  end subroutine dens_xy_z_bld

! ==============================================================================
  subroutine dens_3d_histo_bld
! ------------------------------------------------------------------------------
! Description   : histo estimator of 3D density
!
! Created       : J. Toulouse, 05 Mar 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer grid_i, grid_x_i, grid_y_i, grid_z_i
  integer elec_i, walk_i
  real(dp) xi, yi, zi

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_histo')
   call object_average_walk_define ('dens_3d_histo', 'dens_3d_histo_av')
   call object_error_define ('dens_3d_histo_av', 'dens_3d_histo_av_err')

   call object_needed ('nwalk')
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
   call object_needed ('coord_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_3d_histo', dens_3d_histo, grid_xyz_nb, nwalk)
  call object_alloc ('dens_3d_histo_av', dens_3d_histo_av, grid_xyz_nb)
  call object_alloc ('dens_3d_histo_av_err', dens_3d_histo_av_err, grid_xyz_nb)

  dens_3d_histo (:,:) = 0.d0

  do walk_i = 1, nwalk
    do elec_i = 1, nelec

      xi =  coord_elec_wlk(1,elec_i,walk_i)
      yi =  coord_elec_wlk(2,elec_i,walk_i)
      zi =  coord_elec_wlk(3,elec_i,walk_i)

      grid_x_i = floor((xi - grid_x_min)/grid_x_step - 0.5d0)  + 2
      grid_y_i = floor((yi - grid_y_min)/grid_y_step - 0.5d0)  + 2
      grid_z_i = floor((zi - grid_z_min)/grid_z_step - 0.5d0)  + 2

      if (grid_x_i < 1 .or. grid_x_i > grid_x_nb) cycle
      if (grid_y_i < 1 .or. grid_y_i > grid_y_nb) cycle
      if (grid_z_i < 1 .or. grid_z_i > grid_z_nb) cycle

      grid_i = grid_xyz_index (grid_x_i, grid_y_i, grid_z_i)

      dens_3d_histo (grid_i, walk_i) = dens_3d_histo (grid_i, walk_i) + 1.d0/(grid_x_step*grid_y_step*grid_z_step)

    enddo ! elec_i
  enddo ! walk_i

  end subroutine dens_3d_histo_bld

! ==============================================================================
  subroutine dens_3d_zv1_bld
! ------------------------------------------------------------------------------
! Description   : first-order renormalized improved estimator of 3D density
!
! Created       : J. Toulouse, 07 Mar 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer grid_i, elec_i, dim_i, walk_i
  real(dp) di, di2, dotproduct

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_zv1')
   call object_average_walk_define ('dens_3d_zv1', 'dens_3d_zv1_av')
   call object_error_define ('dens_3d_zv1_av', 'dens_3d_zv1_av_err')

   call object_needed ('nwalk')
   call object_needed ('grid_xyz_nb')
   call object_needed ('grid_xyz')
   call object_needed ('nelec')
   call object_needed ('coord_elec_wlk')
   call object_needed ('grd_psi_over_psi_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_3d_zv1', dens_3d_zv1, grid_xyz_nb, nwalk)
  call object_alloc ('dens_3d_zv1_av', dens_3d_zv1_av, grid_xyz_nb)
  call object_alloc ('dens_3d_zv1_av_err', dens_3d_zv1_av_err, grid_xyz_nb)

  dens_3d_zv1 (:,:) = 0.d0

  do walk_i = 1, nwalk
    do elec_i = 1, nelec

      do grid_i = 1, grid_xyz_nb

!     (ri - r)
      di2 = 0.d0
      do dim_i = 1,ndim
         di2 = di2 + (coord_elec_wlk (dim_i, elec_i, walk_i) - grid_xyz (dim_i, grid_i))**2
      enddo ! dim_i
      di = dsqrt(di2)

!     dot product: drift_i . (ri -r)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct + grd_psi_over_psi_wlk (dim_i, elec_i, walk_i) * (coord_elec_wlk (dim_i, elec_i, walk_i) - grid_xyz (dim_i, grid_i))
      enddo

      dens_3d_zv1 (grid_i, walk_i) = dens_3d_zv1 (grid_i, walk_i) - oneover2pi * dotproduct / di**3

     enddo ! grid_i

    enddo ! elec_i
  enddo ! walk_i

  end subroutine dens_3d_zv1_bld

! ==============================================================================
  subroutine dens_3d_zv2_bld
! ------------------------------------------------------------------------------
! Description   : second-order renormalized improved estimator of 3D density
!
! Created       : J. Toulouse, 07 Mar 2006
! Modified      : J. Toulouse, 13 Aug 2008: walkers
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer grid_i, elec_i, dim_i, walk_i
  real(dp) di, di2

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_zv2')
   call object_average_walk_define ('dens_3d_zv2', 'dens_3d_zv2_av')
   call object_variance_define ('dens_3d_zv2_av', 'dens_3d_zv2_av_var')
   call object_error_define ('dens_3d_zv2_av', 'dens_3d_zv2_av_err')

   call object_needed ('nwalk')
   call object_needed ('grid_xyz_nb')
   call object_needed ('grid_xyz')
   call object_needed ('nelec')
   call object_needed ('lap_psi_over_psi_wlk')
   call object_needed ('grd_psi_over_psi_sq_wlk')
   call object_needed ('coord_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('dens_3d_zv2', dens_3d_zv2, grid_xyz_nb, nwalk)
  call object_alloc ('dens_3d_zv2_av', dens_3d_zv2_av, grid_xyz_nb)
  call object_alloc ('dens_3d_zv2_av_var', dens_3d_zv2_av_var, grid_xyz_nb)
  call object_alloc ('dens_3d_zv2_av_err', dens_3d_zv2_av_err, grid_xyz_nb)

  dens_3d_zv2 (:,:) = 0.d0

  do walk_i = 1, nwalk
    do elec_i = 1, nelec

      do grid_i = 1, grid_xyz_nb

!     (rij - r)
      di2 = 0.d0
      do dim_i = 1,ndim
         di2 = di2 + (coord_elec_wlk (dim_i, elec_i, walk_i)  - grid_xyz (dim_i, grid_i))**2
      enddo ! dim_i
      di = dsqrt(di2)

      dens_3d_zv2 (grid_i, walk_i) = dens_3d_zv2 (grid_i, walk_i) - oneover2pi * (lap_psi_over_psi_wlk (elec_i, walk_i) + grd_psi_over_psi_sq_wlk (elec_i, walk_i)) / di

     enddo ! grid_i

    enddo ! elec_i
  enddo ! walk_i

  end subroutine dens_3d_zv2_bld

! ==============================================================================
  subroutine dens_3d_zv2_deloc_bld
! ------------------------------------------------------------------------------
! Description   : dens_3d_zv2 * deloc
!
! Created       : J. Toulouse, 13 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i, grid_i

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_zv2_deloc')
   call object_average_define ('dens_3d_zv2_deloc', 'dens_3d_zv2_deloc_av')
   call object_covariance_define ('dens_3d_zv2_av', 'deloc_av', 'dens_3d_zv2_av_deloc_av_covar')

   call object_needed ('grid_xyz_nb')
   call object_needed ('param_nb')
   call object_needed ('dens_3d_zv2')
   call object_needed ('deloc')

   return

  endif

! allocations
  call object_alloc ('dens_3d_zv2_deloc', dens_3d_zv2_deloc, grid_xyz_nb, param_nb)
  call object_alloc ('dens_3d_zv2_deloc_av', dens_3d_zv2_deloc_av, grid_xyz_nb, param_nb)
  call object_alloc ('dens_3d_zv2_av_deloc_av_covar', dens_3d_zv2_av_deloc_av_covar, grid_xyz_nb, param_nb)

  do grid_i = 1, grid_xyz_nb
   do param_i = 1, param_nb
    dens_3d_zv2_deloc (grid_i, param_i) = dens_3d_zv2 (grid_i, 1) * deloc (param_i)
   enddo ! param_i
  enddo ! grid_i

  end subroutine dens_3d_zv2_deloc_bld

! ==============================================================================
  subroutine dens_3d_zv2_deloc_covar_bld
! ------------------------------------------------------------------------------
! Description   : covariance : < dens_3d_zv2 * deloc > - < dens_3d_zv2 > * < deloc >
!
! Created       : J. Toulouse, 13 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i, grid_i

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_zv2_deloc_covar')

   call object_needed ('grid_xyz_nb')
   call object_needed ('param_nb')
   call object_needed ('dens_3d_zv2_deloc_av')
   call object_needed ('dens_3d_zv2_av')
   call object_needed ('deloc_av')

   return

  endif

! allocations
  call object_alloc ('dens_3d_zv2_deloc_covar', dens_3d_zv2_deloc_covar, grid_xyz_nb, param_nb)

  do grid_i = 1, grid_xyz_nb
   do param_i = 1, param_nb
    dens_3d_zv2_deloc_covar (grid_i, param_i) = dens_3d_zv2_deloc_av (grid_i, param_i) - dens_3d_zv2_av (grid_i) * deloc_av (param_i)
   enddo ! param_i
  enddo ! grid_i

  end subroutine dens_3d_zv2_deloc_covar_bld

! ==============================================================================
  subroutine dens_3d_zv2_linear_coef_bld
! ------------------------------------------------------------------------------
! Description   : coefficient minimizing the variance of estimator dens_3d_zv2_linear
!
! Created       : J. Toulouse, 13 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer grid_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_zv2_linear_coef')

   call object_needed ('grid_xyz_nb')
   call object_needed ('param_nb')
   call object_needed ('dens_3d_zv2_deloc_covar')
   call object_needed ('deloc_deloc_covar_inv')

   return

  endif

! allocations
  call object_alloc ('dens_3d_zv2_linear_coef', dens_3d_zv2_linear_coef, grid_xyz_nb, param_nb)

  do grid_i = 1, grid_xyz_nb
   do param_i = 1, param_nb
    dens_3d_zv2_linear_coef (grid_i, param_i) = 0.d0
    do param_j = 1, param_nb
     dens_3d_zv2_linear_coef (grid_i, param_i) = dens_3d_zv2_linear_coef (grid_i, param_i) - deloc_deloc_covar_inv (param_i, param_j) * dens_3d_zv2_deloc_covar (grid_i, param_j)
    enddo ! param_j
   enddo ! param_i
  enddo ! grid_i

  end subroutine dens_3d_zv2_linear_coef_bld

! ==============================================================================
  subroutine dens_3d_zv2_linear_av_bld
! ------------------------------------------------------------------------------
! Description   : average of second-order renormalized improved estimator of 3D density
! Description   : including wave function derivatives wrt to parameters
!
! Created       : J. Toulouse, 13 Aug 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer grid_i, param_i

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_zv2_linear_av')

   call object_needed ('grid_xyz_nb')
   call object_needed ('param_nb')
   call object_needed ('dens_3d_zv2_av')
   call object_needed ('dens_3d_zv2_linear_coef')
   call object_needed ('deloc_av')

   return

  endif

! allocations
  call object_alloc ('dens_3d_zv2_linear_av', dens_3d_zv2_linear_av, grid_xyz_nb)

  do grid_i = 1, grid_xyz_nb
   dens_3d_zv2_linear_av (grid_i) = dens_3d_zv2_av (grid_i)
   do param_i = 1, param_nb
    dens_3d_zv2_linear_av (grid_i) = dens_3d_zv2_linear_av (grid_i) + dens_3d_zv2_linear_coef (grid_i, param_i) * deloc_av (param_i)
   enddo ! param_i
  enddo ! grid_i

  end subroutine dens_3d_zv2_linear_av_bld

! ==============================================================================
  subroutine dens_3d_zv2_linear_av_var_bld
! ------------------------------------------------------------------------------
! Description   : variance of average of zero-variance estimator of force
!
! Created       : J. Toulouse, 24 Jul 2008
! ------------------------------------------------------------------------------
  implicit none

! local
  integer grid_i, param_i, param_j

! begin

! header
  if (header_exe) then

   call object_create ('dens_3d_zv2_linear_av_var')
   call object_error_define_from_variance ('dens_3d_zv2_linear_av_var', 'dens_3d_zv2_linear_av_err')

   call object_needed ('grid_xyz_nb')
   call object_needed ('param_nb')
   call object_needed ('dens_3d_zv2_av_var')
   call object_needed ('dens_3d_zv2_linear_coef')
   call object_needed ('dens_3d_zv2_av_deloc_av_covar')
   call object_needed ('deloc_av_deloc_av_covar')

   return

  endif

! allocations
  call object_alloc ('dens_3d_zv2_linear_av_var', dens_3d_zv2_linear_av_var, grid_xyz_nb)
  call object_alloc ('dens_3d_zv2_linear_av_err', dens_3d_zv2_linear_av_err, grid_xyz_nb)

  do grid_i = 1, grid_xyz_nb
   dens_3d_zv2_linear_av_var (grid_i) = dens_3d_zv2_av_var (grid_i)
   do param_i = 1, param_nb
    dens_3d_zv2_linear_av_var (grid_i) = dens_3d_zv2_linear_av_var (grid_i) + 2.d0 * dens_3d_zv2_linear_coef (grid_i, param_i) * dens_3d_zv2_av_deloc_av_covar (grid_i, param_i)
     do param_j = 1, param_nb
      dens_3d_zv2_linear_av_var (grid_i) = dens_3d_zv2_linear_av_var (grid_i) + dens_3d_zv2_linear_coef (grid_i, param_i) * dens_3d_zv2_linear_coef (grid_i, param_j) * deloc_av_deloc_av_covar (param_i, param_j)
     enddo ! param_j
   enddo ! param_i
  enddo ! grid_i

  end subroutine dens_3d_zv2_linear_av_var_bld

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

   case ('zv2_linear')
   call object_provide (lhere,'dens_3d_zv2_linear_av')
   call object_provide (lhere,'dens_3d_zv2_linear_av_err')
   dens_3d (:)     = dens_3d_zv2_linear_av (:)
   dens_3d_err (:) = dens_3d_zv2_linear_av_err (:)

   case default
   call die (lhere, 'unknown estimator >'+trim(dens_3d_estimator)+'<.')
  end select

  end subroutine dens_3d_bld

! ==============================================================================
  subroutine dens_fourier_bld
! ------------------------------------------------------------------------------
! Description   : estimator of periodic density (Fourier transform)
!
! Created       : J. Toulouse, 28 Sep 2013
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer gvec_i, elec_i, dim_i
  real(dp) dotproduct

! begin

! header
  if (header_exe) then

   call object_create ('dens_fourier')
   call object_average_define ('dens_fourier', 'dens_fourier_av')
   call object_error_define ('dens_fourier_av', 'dens_fourier_av_err')

   call object_needed ('ndim')
   call object_needed ('gvec')
   call object_needed ('ngvec_big')
   call object_needed ('nelec')
   call object_needed ('coord_elec')
   call object_needed ('vcell')

   return

  endif

! allocations
  call object_alloc ('dens_fourier', dens_fourier, ngvec_big)
  call object_alloc ('dens_fourier_av', dens_fourier_av, ngvec_big)
  call object_alloc ('dens_fourier_av_err', dens_fourier_av_err, ngvec_big)

  do gvec_i = 1, ngvec_big

    dens_fourier(gvec_i) = 0.d0

    do elec_i = 1, nelec
 
      dotproduct = 0.d0
      do dim_i = 1, ndim
       dotproduct = dotproduct + gvec(dim_i,gvec_i)*coord_elec(dim_i,elec_i)
      enddo ! dim_i

    dens_fourier(gvec_i) = dens_fourier(gvec_i) + cdexp(-(0.0d0,1.0d0)*dotproduct)

    enddo ! elec_i

    dens_fourier(gvec_i) = dens_fourier(gvec_i) / vcell

  enddo ! gvec_i

  end subroutine dens_fourier_bld

! ========================================================================
  subroutine dens_wrt
! ------------------------------------------------------------------------
! Description    : write spherically averaged density density on file
!
! Created        : J. Toulouse, 07 Nov 2006
! ------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer unit, grid_i

! begin

! return if not main node
# if defined (MPI)
   if (idtask /= 0) return
# endif

! provide necessary objects
  call object_provide ('grid_r_nb')
  call object_provide ('grid_r')
  call object_provide ('dens')
  call object_provide ('dens_err')

! open file
  unit = 0
  call open_file_or_die (dens_file_out, unit)

  write(unit,'(9a)')        'spherically averaged density generated by CHAMP using estimator "',trim(dens_estimator),'" and mode "',trim(mode)
  write(unit,'(a,i5)')      'number of electrons       =',nelec
  write(unit,'(a,i15)')     'number of steps per block =',nstep_total
  write(unit,'(a,i15)')     'number of blocks          =',block_iterations_nb
  write(unit,'(a,3e15.8)')  'grid_r_step               =',grid_r_step
  write(unit,'(a,3e15.8)')  'grid_r_min                =',grid_r_min
  write(unit,'(a,3e15.8)')  'grid_r_max                =',grid_r_max
  write(unit,'(a,i15)')     'grid_r_nb                 =',grid_r_nb
  write(unit,'(a,es15.8)')  'dist_e_min                =',dist_e_min
  write(unit,'(a,es15.8)')  'dist_e_max                =',dist_e_max

  write(unit,*) ''
  write(unit,'(a)') '    r                 n(r)            error of n(r)'

  do grid_i = 1, grid_r_nb
       write(unit,'(es13.6,2es21.14)') grid_r (grid_i), dens (grid_i), dens_err (grid_i)
  enddo ! grid_i

  close(unit)

  end subroutine dens_wrt

! ========================================================================
  subroutine dens_xy_z_wrt
! ------------------------------------------------------------------------
! Description    : write axial density density on file
!
! Created       : C. Umrigar, 07 Jul 2012, by imitating dens_wrt of J. Toulouse
! ------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer                                    :: unit
  integer                                    :: grid_i, grid_xy_i, grid_z_i

! begin

! return if not main node
# if defined (MPI)
   if (idtask /= 0) return
# endif

! provide necessary objects
  call object_provide ('grid_xy_z_nb')
  call object_provide ('grid_xy_z')
  call object_provide ('dens_xy_z')
  call object_provide ('dens_xy_z_err')

! open file
  unit = 0
  call open_file_or_die (dens_xy_z_file_out, unit)

  write(unit,'(9a)')        'Axially symmetric density generated by CHAMP using estimator "',trim(dens_xy_z_estimator),'" and mode "',trim(mode)
  write(unit,'(a,i5)')      'number of electrons       =',nelec
  write(unit,'(a,i20)')     'number of walkers         =',nwalk
  write(unit,'(a,i20)')     'number of steps per block =',nstep_total
  write(unit,'(a,i20)')     'number of blocks          =',block_iterations_nb
  write(unit,'(a,3es12.4)') 'grid_xy_max               =',grid_xy_max
  write(unit,'(a,3es12.4)') 'grid_xy_min               =',grid_xy_min
  write(unit,'(a,3es12.4)') 'grid_xy_step              =',grid_xy_step
  write(unit,'(a,3es12.4)') 'grid_z_max                =',grid_z_max
  write(unit,'(a,3es12.4)') 'grid_z_min                =',grid_z_min
  write(unit,'(a,3es12.4)') 'grid_z_step               =',grid_z_step
  write(unit,'(a,i12)')     'grid_xy_z_nb              =',grid_xy_z_nb

  write(unit,*) ''
  write(unit,'(a)') '   xy          z           n(r)        error_n(r)'

  do grid_xy_i = 1, grid_xy_nb
    do grid_z_i = 1, grid_z_nb
      grid_i = grid_xy_z_index (grid_xy_i, grid_z_i)
!     write(unit,'(5e25.15)') grid_xy_z(1,grid_i), grid_xy_z(2,grid_i), grid_xy_z(3,grid_i), dens_xy_z (grid_i), dens_xy_z_err(grid_i)
      write(unit,'(2f10.6,es17.8,es10.2)') grid_xy_z(1,grid_i), grid_xy_z(2,grid_i), dens_xy_z (grid_i), dens_xy_z_err(grid_i)
    enddo ! grid_z_i
    write(unit,'(a)') ''
  enddo ! grid_xy_i

  close(unit)

  end subroutine dens_xy_z_wrt

! ========================================================================
  subroutine dens_3d_wrt
! ------------------------------------------------------------------------
! Description    : write 3D density on file
!
! Created        : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer                                    :: unit
  integer                                    :: grid_i, grid_x_i, grid_y_i, grid_z_i

! begin

! return if not main node
# if defined (MPI)
   if (idtask /= 0) return
# endif

! provide necessary objects
  call object_provide ('grid_xyz_nb')
  call object_provide ('grid_xyz')
  call object_provide ('dens_3d')
  call object_provide ('dens_3d_err')

! open file
  unit = 0
  call open_file_or_die (dens_3d_file_out, unit)

  write(unit,'(9a)')        '#D density generated by CHAMP using estimator "',trim(dens_3d_estimator),'" and mode "',trim(mode)
  write(unit,'(a,i5)')      'number of electrons       =',nelec
  write(unit,'(a,i20)')     'number of walkers         =',nwalk
  write(unit,'(a,i20)')     'number of steps per block =',nstep_total
  write(unit,'(a,i20)')     'number of blocks          =',block_iterations_nb
  write(unit,'(a,3es12.4)') 'grid_x_max                =',grid_x_max
  write(unit,'(a,3es12.4)') 'grid_x_min                =',grid_x_min
  write(unit,'(a,3es12.4)') 'grid_x_step               =',grid_x_step
  write(unit,'(a,3es12.4)') 'grid_y_max                =',grid_y_max
  write(unit,'(a,3es12.4)') 'grid_y_min                =',grid_y_min
  write(unit,'(a,3es12.4)') 'grid_y_step               =',grid_y_step
  write(unit,'(a,3es12.4)') 'grid_z_max                =',grid_z_max
  write(unit,'(a,3es12.4)') 'grid_z_min                =',grid_z_min
  write(unit,'(a,3es12.4)') 'grid_z_step               =',grid_z_step
  write(unit,'(a,i12)')     'grid_xyz_nb               =',grid_xyz_nb

  write(unit,*) ''
  write(unit,'(a)') '   x         y         z           n(r)         error_n(r)'

  do grid_x_i = 1, grid_x_nb
   do grid_y_i = 1, grid_y_nb
    do grid_z_i = 1, grid_z_nb

       grid_i = grid_xyz_index (grid_x_i, grid_y_i, grid_z_i)
       write(unit,'(3es13.6,es18.8,es10.2)') grid_xyz(1,grid_i), grid_xyz(2,grid_i), grid_xyz(3,grid_i), dens_3d (grid_i), dens_3d_err(grid_i)
    enddo ! grid_z_i
    write(unit,'(a)') ''
   enddo  ! grid_y_i
   write(unit,'(a)') ''
  enddo ! grid_x_i

  close(unit)

  end subroutine dens_3d_wrt

! ========================================================================
  subroutine dens_fourier_wrt
! ------------------------------------------------------------------------
! Description    : write density in fourier space on file
!
! Created        : J. Toulouse, 28 Sep 2013
! ------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer                                    :: unit
  integer                                    :: gvec_i

! begin

! return if not main node
# if defined (MPI)
   if (idtask /= 0) return
# endif

! provide necessary objects
  call object_provide ('ngvec_big')
  call object_provide ('gvec')
  call object_provide ('dens_fourier_av')
  call object_provide ('dens_fourier_av_err')

! open file
  unit = 0
  call open_file_or_die (dens_fourier_file_out, unit)

  write(unit,'(9a)')        '#D fourier density generated by CHAMP using mode ',trim(mode)
  write(unit,'(a,i5)')      'number of electrons       =',nelec
  write(unit,'(a,i20)')     'number of walkers         =',nwalk
  write(unit,'(a,i20)')     'number of steps per block =',nstep_total
  write(unit,'(a,i20)')     'number of blocks          =',block_iterations_nb
  write(unit,'(a,i12)')     'ngvec_big                 =',ngvec_big

  write(unit,*) ''
  write(unit,'(a)') '                        gvec                     density Fourier components       statistical error'

  do gvec_i = 1, ngvec_big
    write(unit,'(i6,3f12.6,es17.8,a,es16.8,es10.2,a,es9.2)') gvec_i, gvec(1,gvec_i), gvec(2,gvec_i), gvec(3,gvec_i), real(dens_fourier_av (gvec_i)), ' +i*',aimag(dens_fourier_av (gvec_i)), real(dens_fourier_av_err(gvec_i)), ' +i*',aimag(dens_fourier_av_err(gvec_i))
  enddo

  close(unit)

  end subroutine dens_fourier_wrt

end module density_mod
