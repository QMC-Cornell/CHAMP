 module intracule_mod


  use all_tools_mod
  use grid_mod
  use electrons_mod
  use psi_mod

! Declaration of global variables and default values
  character(len=max_string_len_file)  :: intra_file_out  = ''
  character(len=max_string_len_file)  :: intra_3d_file_out  = ''

  character(len=max_string_len)     :: intra_estimator  = 'zv1'
  character(len=max_string_len)     :: intra_3d_estimator  = 'zv2'

  integer, allocatable              :: spin_two_elec_wlk (:,:,:)

  integer                           :: intra_block_nb = 0

  real(dp), allocatable             :: sum_intra_imp_num(:)
  real(dp), allocatable             :: intra_imp_num(:)
  real(dp), allocatable             :: intra_imp_norm(:)
  real(dp), allocatable             :: sum_intra_imp_norm(:)
  real(dp), allocatable             :: sum_intra_imp_norm_square(:)
  real(dp), allocatable             :: intra_imp_f_av(:)
  real(dp), allocatable             :: intra_imp_f_var(:)
  real(dp), allocatable             :: intra_imp_f_av_err(:)
  real(dp), allocatable             :: intra_imp_f4pir2_av(:)
  real(dp), allocatable             :: intra_imp_f4pir2_av_err(:)
  integer                           :: intra_imp_nb
  logical                           :: init_intra_imp_step = .true.
  logical                           :: init_intra_imp_block = .true.

! intracules
  real(dp)                  :: intra_exp = 1.d0
  real(dp), allocatable     :: intra (:)
  real(dp), allocatable     :: intra_err (:)
  real(dp), allocatable     :: intra_4pir2 (:)
  real(dp), allocatable     :: intra_4pir2_err (:)
  real(dp), allocatable     :: intra_histo (:)
  real(dp), allocatable     :: intra_histo_av (:)
  real(dp), allocatable     :: intra_histo_av_err (:)
  real(dp), allocatable     :: intra_zv1 (:)
  real(dp), allocatable     :: intra_zv1_av (:)
  real(dp), allocatable     :: intra_zv1_av_err (:)
  real(dp), allocatable     :: intra_zv2 (:)
  real(dp), allocatable     :: intra_zv2_av (:)
  real(dp), allocatable     :: intra_zv2_av_err (:)

! spin-resolved intracules
  real(dp), allocatable     :: intra_sp (:,:)
  real(dp), allocatable     :: intra_sp_err (:,:)
  real(dp), allocatable     :: intra_sp_4pir2 (:,:)
  real(dp), allocatable     :: intra_sp_4pir2_err (:,:)
  real(dp), allocatable     :: intra_sp_histo (:,:)
  real(dp), allocatable     :: intra_sp_histo_av (:,:)
  real(dp), allocatable     :: intra_sp_histo_av_err (:,:)
  real(dp), allocatable     :: intra_sp_zv1 (:,:)
  real(dp), allocatable     :: intra_sp_zv1_av (:,:)
  real(dp), allocatable     :: intra_sp_zv1_av_err (:,:)
  real(dp), allocatable     :: intra_sp_zv2 (:,:)
  real(dp), allocatable     :: intra_sp_zv2_av (:,:)
  real(dp), allocatable     :: intra_sp_zv2_av_err (:,:)
  real(dp), allocatable     :: intra_sp_zv3 (:,:)
  real(dp), allocatable     :: intra_sp_zv3_av (:,:)
  real(dp), allocatable     :: intra_sp_zv3_av_err (:,:)
  real(dp), allocatable     :: intra_sp_zv4 (:,:)
  real(dp), allocatable     :: intra_sp_zv4_av (:,:)
  real(dp), allocatable     :: intra_sp_zv4_av_err (:,:)
  real(dp), allocatable     :: intra_sp_zv5 (:,:,:)
  real(dp), allocatable     :: intra_sp_zv5_av (:,:)
  real(dp), allocatable     :: intra_sp_zv5_av_err (:,:)
  real(dp), allocatable     :: intra_sp_q1 (:,:)
  real(dp), allocatable     :: intra_sp_q1_av (:,:)
  real(dp), allocatable     :: intra_sp_q1_eloc (:,:)
  real(dp), allocatable     :: intra_sp_q1_eloc_av (:,:)
  real(dp), allocatable     :: intra_sp_zb1_av (:,:)
  real(dp), allocatable     :: intra_sp_zvzb1 (:,:)
  real(dp), allocatable     :: intra_sp_zvzb1_av (:,:)
  real(dp), allocatable     :: intra_sp_zvzb1_av_err (:,:)
  real(dp), allocatable     :: intra_sp_zvzb3 (:,:)
  real(dp), allocatable     :: intra_sp_zvzb3_av (:,:)
  real(dp), allocatable     :: intra_sp_zvzb3_av_err (:,:)
  real(dp), allocatable     :: intra_sp_zvzb4 (:,:)
  real(dp), allocatable     :: intra_sp_zvzb4_av (:,:)
  real(dp), allocatable     :: intra_sp_zvzb4_av_err (:,:)
  real(dp), allocatable     :: intra_sp_zv1zb3 (:,:)
  real(dp), allocatable     :: intra_sp_zv1zb3_av (:,:)
  real(dp), allocatable     :: intra_sp_zv1zb3_av_err (:,:)
  real(dp), allocatable     :: intra_sp_q5_wlk (:,:,:)
  real(dp), allocatable     :: intra_sp_q5_av (:,:)
  real(dp), allocatable     :: intra_sp_q5_eloc_wlk (:,:,:)
  real(dp), allocatable     :: intra_sp_q5_eloc_av (:,:)
  real(dp), allocatable     :: intra_sp_zb5_av (:,:)
  real(dp), allocatable     :: intra_sp_zvzb5    (:,:,:)
  real(dp), allocatable     :: intra_sp_zvzb5_av (:,:)
  real(dp), allocatable     :: intra_sp_zvzb5_av_err (:,:)

! 3D intracules
  real(dp), allocatable     :: intra_3d_histo (:)
  real(dp), allocatable     :: intra_3d_histo_av (:)
  real(dp), allocatable     :: intra_3d_histo_av_err (:)
  real(dp), allocatable     :: intra_3d_zv1 (:)
  real(dp), allocatable     :: intra_3d_zv1_av (:)
  real(dp), allocatable     :: intra_3d_zv1_av_err (:)
  real(dp), allocatable     :: intra_3d_zv2 (:)
  real(dp), allocatable     :: intra_3d_zv2_av (:)
  real(dp), allocatable     :: intra_3d_zv2_av_err (:)
  real(dp), allocatable     :: intra_3d (:)
  real(dp), allocatable     :: intra_3d_err (:)

  contains

!===========================================================================
  subroutine intra_menu
!---------------------------------------------------------------------------
! Description : menu for intracule calculation
!
! Created     : J. Toulouse, 5 Oct 2005
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

  character(len=max_string_len_rout), save :: lhere = 'intra_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))
   case ('help')
   write(6,'(2a)') trim(lhere),': help for intracule menu'
   write(6,'(2a)') trim(lhere),': intracule'
   write(6,'(2a)') trim(lhere),':  estimator = {histogram|zv1|zv2|zv3|zv4|zv5|zvzb1|zvzb3|zvzb4|zv1zb3|zvzb5} (default: zv1)'
   write(6,'(2a)') trim(lhere),':  exponent  = [real] (default: 1) exponent for zv3, zv4, zvzb3, zvzb4, zv1zb3'
   write(6,'(2a)') trim(lhere),':  file      = [string] : file in which intracule will be written'
   write(6,'(2a)') trim(lhere),': end'

  case ('file')
   call get_next_value (intra_file_out)

  case ('estimator')
   call get_next_value (intra_estimator)

  case ('exponent')
   call get_next_value (intra_exp)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines


! File
  if (trim(intra_file_out) /= '') then
   write (6,'(3a)') trim(lhere),': intracule will be written on file: ',trim(intra_file_out)
  else
   write (6,'(2a)') trim(lhere),': file for writing intracule not specified'
   call die (lhere)
  endif

  select case(trim(intra_estimator))
   case ('histogram')
   write (6,*) trim(lhere),': intracule will be calculated with estimator histogram'
   call object_average_request ('intra_sp_histo_av')
   call object_error_request ('intra_sp_histo_av_err')

   case ('zv1')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zv1'
   call object_average_request ('intra_sp_zv1_av')
   call object_error_request ('intra_sp_zv1_av_err')

   case ('zv2')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zv2'
   call object_average_request ('intra_sp_zv2_av')
   call object_error_request ('intra_sp_zv2_av_err')

   case ('zv3')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zv3'
   write (6,*) trim(lhere),': exponent=',intra_exp
   call object_average_request ('intra_sp_zv3_av')
   call object_error_request ('intra_sp_zv3_av_err')

   case ('zv4')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zv4'
   write (6,*) trim(lhere),': exponent=',intra_exp
   call object_average_request ('intra_sp_zv4_av')
   call object_error_request ('intra_sp_zv4_av_err')

   case ('zv5')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zv5'
   write (6,*) trim(lhere),': exponent=',intra_exp
   call object_average_request ('intra_sp_zv5_av')
   call object_error_request ('intra_sp_zv5_av_err')

   case ('zvzb1')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zvzb1'
!   call object_average_request ('intra_sp_zv1_av')
!   call object_average_request ('intra_sp_q1_av')
!   call object_average_request ('intra_sp_q1_eloc_av')
   call object_average_request ('intra_sp_zvzb1_av')
   call object_error_request ('intra_sp_zvzb1_av_err')

   case ('zvzb3')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zvzb3'
   write (6,*) trim(lhere),': exponent=',intra_exp
   call object_average_request ('intra_sp_zvzb3_av')
   call object_error_request ('intra_sp_zvzb3_av_err')

   case ('zvzb4')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zvzb4'
   write (6,*) trim(lhere),': exponent=',intra_exp
   call object_average_request ('intra_sp_zvzb4_av')
   call object_error_request ('intra_sp_zvzb4_av_err')

   case ('zv1zb3')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zv1zb3'
   write (6,*) trim(lhere),': exponent=',intra_exp
   call object_average_request ('intra_sp_zv1zb3_av')
   call object_error_request ('intra_sp_zv1zb3_av_err')

   case ('zvzb5')
   write (6,*) trim(lhere),': intracule will be calculated with estimator zvzb5'
   write (6,*) trim(lhere),': exponent=',intra_exp
!   call object_average_request ('intra_sp_zv5_av')
!   call object_average_request ('intra_sp_q5_av')
!   call object_average_request ('intra_sp_q5_eloc_av')
   call object_average_request ('intra_sp_zvzb5_av')
   call object_error_request ('intra_sp_zvzb5_av_err')

   case default
   write(6,*) trim(lhere), ': unknown estimator >',trim(intra_estimator),'<'
   call die (lhere)
  end select

  call routine_write_block_request  ('intra_wrt')

  end subroutine intra_menu

!===========================================================================
  subroutine intra_3d_menu
!---------------------------------------------------------------------------
! Description : menu for 3D intracule calculation
!
! Created     : J. Toulouse, 04 Mar 2006
!---------------------------------------------------------------------------
  implicit none

  character(len=max_string_len_rout), save :: lhere = 'intra_3d_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  if(trim(word) == 'help') then
   write(6,*) 'intra_3d-h: menu for intra_3d'
   write(6,*) 'intra_3d-h: intra_3d'
   write(6,*) 'intra_3d-h:  estimator = [string] {histogram|zv1|zv2|}'
   write(6,*) 'intra_3d-h:  file    = [string] file in which intracule will be written'
   write(6,*) 'intra_3d-h: end'

  elseif(trim(word) == 'file') then
   call get_next_value (intra_3d_file_out)

  elseif(trim(word) == 'estimator') then
   call get_next_value (intra_3d_estimator)

  elseif(trim(word) == 'end') then
   exit

  else

   write(6,*) trim(lhere),': unknown word = ',trim(word)
   call die(here)

  endif

  enddo ! end loop over menu lines


! File
  if (trim(intra_3d_file_out) /= '') then
   write (6,*) trim(lhere),': intracule will be written on file: ',trim(intra_3d_file_out)
  else
   write (6,*) trim(lhere),': file for writing 3D intracule not specified'
   call die(here)
  endif

  select case(trim(intra_3d_estimator))
   case ('histogram')
   call object_average_request ('intra_3d_histo_av')
   call object_error_request ('intra_3d_histo_av_err')

   case ('zv1')
   call object_average_request ('intra_3d_zv1_av')
   call object_error_request ('intra_3d_zv1_av_err')

   case ('zv2')
   call object_average_request ('intra_3d_zv2_av')
   call object_error_request ('intra_3d_zv2_av_err')

   case default
   write(6,*) trim(lhere), ': unknown estimator >',trim(intra_3d_estimator),'<'
   call die (lhere)
  end select
  call routine_write_block_request  ('intra_3d_wrt')

  end subroutine intra_3d_menu

! ==============================================================================
  subroutine spin_two_elec_wlk_bld
! ------------------------------------------------------------------------------
! Description   : spin_two_elec_wlk = 1 for two electrons of same spins
! Description   :                   = 2 for two electrons of opposite spins
!
! Created       : J. Toulouse, 16 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: elec_i, elec_j

! begin

! header
  if (header_exe) then

   call object_create ('spin_two_elec_wlk')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('nup')
   call object_needed ('ndn')

   return

  endif

! allocations
  call object_alloc ('spin_two_elec_wlk', spin_two_elec_wlk, nelec, nelec, nwalk)

  do elec_i = 1, nelec
    do elec_j = 1, nelec

!   the two electrons are of same spin
    if ( (elec_i <= nup .and. elec_j <= nup) .or. (elec_i > nup .and. elec_j > nup) ) then

      spin_two_elec_wlk (elec_i, elec_j, :) = 1

!   the two electrons are of opposite spin
    else

      spin_two_elec_wlk (elec_i, elec_j, :) = 2

    endif

    enddo ! elec_j
  enddo ! elec_i

  end subroutine spin_two_elec_wlk_bld

! ==============================================================================
  subroutine intra_sp_histo_bld
! ------------------------------------------------------------------------------
! Description   : histogram estimator of spin-resolved intracule
!
! Created       : J. Toulouse, 05 Mar 2006
! Revised       : J. Toulouse, 20 Mar 2006, spin-resolved
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: dij, dij2
  real(dp)              :: d4pir2
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_histo')

   call object_needed ('grid_r')
   call object_needed ('grid_r_step')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('spin_two_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_sp_histo', intra_sp_histo, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_histo_av', intra_sp_histo_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_histo_av_err', intra_sp_histo_av_err, spin_nb, grid_r_nb)

  intra_sp_histo (:,:) = 0.d0

  do elec_j = 1, nelec
    do elec_i = elec_j+1, nelec

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      dij = dist_ee_wlk (elec_i, elec_j, 1)

      grid_i = floor(dij/grid_r_step - 0.5d0) + 2

      if (grid_i <= grid_r_nb) then
       d4pir2 = 4.d0*pi*grid_r(grid_i)**2
       intra_sp_histo (spin, grid_i) = intra_sp_histo (spin, grid_i) + 1.d0/(d4pir2 * grid_r_step)
      endif

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_histo_bld

! ==============================================================================
  subroutine intra_sp_zv1_bld
! ------------------------------------------------------------------------------
! Description   : zero-variance improved estimator of spin_resolved intracule
!
! Created       : J. Toulouse, 06 Mar 2006
! Revised       : J. Toulouse, 20 Mar 2006, spin-resolved
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: dij, dij2, r
  real(dp)              :: dotproduct, intra_temp
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zv1')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zv1', intra_sp_zv1, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv1_av', intra_sp_zv1_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv1_av_err', intra_sp_zv1_av_err, spin_nb, grid_r_nb)

  intra_sp_zv1 (:,:) = 0.d0

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      dij = dist_ee_wlk (elec_i, elec_j, 1)

!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_j, 1) * vec_ee_xyz_wlk (dim_i, elec_j, elec_i, 1)
      enddo

      intra_temp = - oneover4pi * dotproduct / dij**3


      do grid_i = 1, grid_r_nb
        r = grid_r(grid_i)

        if ( dij >= r ) then
           intra_sp_zv1 (spin, grid_i) = intra_sp_zv1 (spin, grid_i) + intra_temp
        endif
     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_zv1_bld

! ==============================================================================
  subroutine intra_sp_zv2_bld
! ------------------------------------------------------------------------------
! Description   : second-order renormalized improved estimator of intracule
!
! Created       : J. Toulouse, 06 Mar 2006
! Revised       : J. Toulouse, 20 Mar 2006, spin-resolved
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: dij, dij2, r
  real(dp)              :: dotproduct, intra_temp
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zv2')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('grd_psi_over_psi_sq_wlk')
   call object_needed ('lap_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zv2', intra_sp_zv2, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv2_av', intra_sp_zv2_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv2_av_err', intra_sp_zv2_av_err, spin_nb, grid_r_nb)

  intra_sp_zv2 (:,:) = 0.d0

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      dij = dist_ee_wlk (elec_i, elec_j, 1)

      intra_temp = - oneover4pi * (lap_psi_over_psi_wlk (elec_j, 1) + grd_psi_over_psi_sq_wlk (elec_j, 1))

      do grid_i = 1, grid_r_nb
        r = grid_r (grid_i)
        if ( dij >= r ) then
           intra_sp_zv2 (spin, grid_i) = intra_sp_zv2 (spin, grid_i) + intra_temp / dij
        else
           intra_sp_zv2 (spin, grid_i) = intra_sp_zv2 (spin, grid_i) + intra_temp / r
        endif
     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_zv2_bld

! ==============================================================================
  subroutine intra_sp_zv3_bld
! ------------------------------------------------------------------------------
! Description   : improved zero-variance estimator of spin_resolved intracule
! Description   : Q(R) = (-1/4pi) sum_i sim_{j<i} 1/|rij -u| f(rij,u)
! Description   : with f(rij,u) =exp [-c (rij-u)^2]
!
! Created       : J. Toulouse, 03 Aug 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: cc, sqrtc, sqrtpi_c
  real(dp)              :: r, u, rpu, rpu2, ru, rmu, rmu2, abs_rmu
  real(dp)              :: exp_mcrpu2, exp_mcrmu2
  real(dp)              :: dotproduct
  real(dp)              :: intra_temp1, intra_temp2, intra_temp3
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zv3')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zv3', intra_sp_zv3, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv3_av', intra_sp_zv3_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv3_av_err', intra_sp_zv3_av_err, spin_nb, grid_r_nb)

  intra_sp_zv3 (:,:) = 0.d0

  cc= intra_exp

  sqrtc  = dsqrt(cc)
  sqrtpi_c = sqrt_pi/sqrtc

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      r = dist_ee_wlk (elec_i, elec_j, 1)
!
!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_j, 1) * vec_ee_xyz_wlk (dim_i, elec_j, elec_i, 1)
      enddo


!     loop over grid points
      do grid_i = 1, grid_r_nb

        u = grid_r (grid_i)

        rpu = r + u
        rpu2 = rpu * rpu

        rmu   = r - u
        rmu2   = rmu * rmu
        abs_rmu = dabs(rmu)

        ru = r * u

        exp_mcrpu2 = dexp(-cc*rpu2)
        exp_mcrmu2 = dexp(-cc*rmu2)

        intra_temp1 = 0.5d0 * (exp_mcrpu2 - sign(1.d0,rmu) * exp_mcrmu2)
        intra_temp2 = 0.5d0 * (exp_mcrpu2 - (abs_rmu/rpu)* exp_mcrmu2)
        intra_temp3 = sqrtpi_c * (erf_spline(sqrtc*rpu) - erf_spline(sqrtc*abs_rmu))/(4.d0*ru)

        intra_sp_zv3 (spin, grid_i) = intra_sp_zv3 (spin, grid_i) +          &
            oneover4pi*((dotproduct/r)*(intra_temp1 - u*intra_temp3) - cc*rpu*intra_temp2)/ru

     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_zv3_bld

! ==============================================================================
  subroutine intra_sp_zv4_bld
! ------------------------------------------------------------------------------
! Description   : improved estimator of spin_resolved intracule determined from 1D
! Description   : not good
!
! Created       : J. Toulouse, 17 Aug 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: r, u, rmu, rmu2
  real(dp)              :: dotproduct
  real(dp)              :: cc
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zv4')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zv4', intra_sp_zv4, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv4_av', intra_sp_zv4_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv4_av_err', intra_sp_zv4_av_err, spin_nb, grid_r_nb)

  intra_sp_zv4 (:,:) = 0.d0

  cc= intra_exp

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      r = dist_ee_wlk (elec_i, elec_j, 1)

!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_j, 1) * vec_ee_xyz_wlk (dim_i, elec_j, elec_i, 1)
      enddo

      do grid_i = 1, grid_r_nb
        u = grid_r(grid_i)

        if ( r >= u ) then
           intra_sp_zv4 (spin, grid_i) = intra_sp_zv4 (spin, grid_i) - oneover4pi * dotproduct/r**3
        else
           rmu = r-u
           rmu2 = rmu * rmu
           intra_sp_zv4 (spin, grid_i) = intra_sp_zv4 (spin, grid_i) - oneover4pi * cc*                      &
                                (1.d0 + 2.d0*rmu/r - 2.d0*cc*rmu2 + 2.d0*dotproduct*rmu/r) * dexp(-cc*rmu2)
        endif
     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_zv4_bld

! ==============================================================================
  subroutine intra_sp_zv5_bld
! ------------------------------------------------------------------------------
! Description   : improved zero-variance estimator of spin_resolved intracule
! Description   : Q(R) = (-1/4pi) sum_i sim_{j<i} 1/|rij -u| f(rij,u)
! Description   : with f(rij,u) =exp [-c |rij-u|]
!
! Created       : J. Toulouse, 08 Nov 2006
! Modified      : J. Toulouse, 16 Nov 2006, loop over walkers
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer               :: grid_i
  integer               :: elec_i, elec_j, dim_i, walk_i
  real(dp)              :: cc, cc2
  real(dp)              :: r, r2, u, rpu, ru, rmu, abs_rmu
  real(dp)              :: mcabsrmu, mcrpu, exp_mcabsrmu, exp_mcrpu
  real(dp)              :: dotproduct
  real(dp)              :: intra_temp1, intra_temp2
  integer               :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zv5')
   call object_average_walk_define ('intra_sp_zv5', 'intra_sp_zv5_av')
   call object_error_define ('intra_sp_zv5_av', 'intra_sp_zv5_av_err')

   call object_needed ('nwalk')
   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zv5', intra_sp_zv5, spin_nb, grid_r_nb, nwalk)
  call object_alloc ('intra_sp_zv5_av', intra_sp_zv5_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv5_av_err', intra_sp_zv5_av_err, spin_nb, grid_r_nb)

  intra_sp_zv5 (:,:,:) = 0.d0

  cc = intra_exp
  cc2 = cc*cc

  do walk_i = 1, nwalk
   do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, walk_i)

!     distance e-e
      r = dist_ee_wlk (elec_i, elec_j, walk_i)
      r2 = r*r
!
!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk  (dim_i, elec_j, walk_i) * vec_ee_xyz_wlk (dim_i, elec_j, elec_i, walk_i)
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

        intra_temp1 = (exp_mcabsrmu - exp_mcrpu)/(2.d0*cc*ru)
        intra_temp2 = 0.5d0 * (sign(1.d0,rmu) * exp_mcabsrmu - exp_mcrpu)

        intra_sp_zv5 (spin, grid_i, walk_i) = intra_sp_zv5 (spin, grid_i, walk_i) -          &
            oneover4pi*((dotproduct/r2)*(intra_temp1 + intra_temp2/u)  &
            - 0.5d0*cc2*intra_temp1)

     enddo !grid_i

    enddo !elec_j
   enddo !elec_i
  enddo !walk_i

  end subroutine intra_sp_zv5_bld

! ==============================================================================
  subroutine intra_sp_q1_bld
! ------------------------------------------------------------------------------
! Description   : term for zero-bias correction for spin-resolved intracule
!
! Created       : J. Toulouse, 20 Jul 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j
  real(dp)              :: dij,  r
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_q1')
   call object_average_define ('intra_sp_q1', 'intra_sp_q1_av')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('spin_two_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_sp_q1', intra_sp_q1, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_q1_av', intra_sp_q1_av, spin_nb, grid_r_nb)

  intra_sp_q1 (:,:) = 0.d0

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      dij = dist_ee_wlk (elec_i, elec_j, 1)

      do grid_i = 1, grid_r_nb
        r = grid_r(grid_i)

        if ( dij >= r ) then
           intra_sp_q1 (spin, grid_i) = intra_sp_q1 (spin, grid_i) + (1.d0/dij)
         else
           intra_sp_q1 (spin, grid_i) = intra_sp_q1 (spin, grid_i) + (1.d0/r)
        endif
     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_q1_bld

! ==============================================================================
  subroutine intra_sp_q1_eloc_bld
! ------------------------------------------------------------------------------
! Description   : term for zero-bias correction for spin-resolved intracule
!
! Created       : J. Toulouse, 20 Jul 2006
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_q1_eloc')
   call object_average_define ('intra_sp_q1_eloc', 'intra_sp_q1_eloc_av')

   call object_needed ('grid_r_nb')
   call object_needed ('intra_sp_q1')
   call object_needed ('eloc')

   return

  endif

! allocations
  call object_alloc ('intra_sp_q1_eloc', intra_sp_q1_eloc, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_q1_eloc_av', intra_sp_q1_eloc_av, spin_nb, grid_r_nb)

  intra_sp_q1_eloc (:,:) = intra_sp_q1 (:,:) * eloc

  end subroutine intra_sp_q1_eloc_bld

! ==============================================================================
  subroutine intra_sp_zb1_av_bld
! ------------------------------------------------------------------------------
! Description   : zero-bias correction for spin-resolved intracule
!
! Created       : J. Toulouse, 20 Jul 2006
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zb1_av')

   call object_needed ('grid_r_nb')
   call object_needed ('intra_sp_q1_eloc_av')
   call object_needed ('intra_sp_q1_av')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zb1_av', intra_sp_zb1_av, spin_nb, grid_r_nb)

  intra_sp_zb1_av (:,:) = -oneover4pi * (intra_sp_q1_eloc_av (:,:) - intra_sp_q1_av (:,:) * eloc_av )

  end subroutine intra_sp_zb1_av_bld

! ==============================================================================
  subroutine intra_sp_zvzb1_av_bld
! ------------------------------------------------------------------------------
! Description   : zero-variance zero-bias correction for spin-resolved intracule
!
! Created       : J. Toulouse, 20 Jul 2006
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zvzb1_av')
   call object_error_define ('intra_sp_zvzb1_av', 'intra_sp_zvzb1_av_err')

   call object_needed ('grid_r_nb')
   call object_needed ('intra_sp_zv1_av')
   call object_needed ('intra_sp_zb1_av')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zvzb1_av', intra_sp_zvzb1_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zvzb1_av_err', intra_sp_zvzb1_av_err, spin_nb, grid_r_nb)

  intra_sp_zvzb1_av (:,:) = intra_sp_zv1_av (:,:) + intra_sp_zb1_av (:,:)

  end subroutine intra_sp_zvzb1_av_bld

! ==============================================================================
  subroutine intra_sp_q5_wlk_bld
! ------------------------------------------------------------------------------
! Description   : term for zero-bias correction for spin-resolved intracule
! Description   : Q = sum_{i,j} \int dOmega_u/(4pi) exp [-zeta |rij-u|] / |rij-u|
!
! Created       : J. Toulouse, 21 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer               :: walk_i, elec_i, elec_j, grid_i
  real(dp)              :: dij,  r
  real(dp)              :: u, rpu, ru, rmu, abs_rmu
  real(dp)              :: mcabsrmu, mcrpu, exp_mcabsrmu, exp_mcrpu
  integer               :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_q5_wlk')
   call object_average_walk_define ('intra_sp_q5_wlk', 'intra_sp_q5_av')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('spin_two_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_sp_q5_wlk', intra_sp_q5_wlk, spin_nb, grid_r_nb, nwalk)
  call object_alloc ('intra_sp_q5_av', intra_sp_q5_av, spin_nb, grid_r_nb)

  intra_sp_q5_wlk (:,:,:) = 0.d0

  do walk_i = 1, nwalk

   do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, walk_i)

!     distance e-e
      r = dist_ee_wlk (elec_i, elec_j, walk_i)

!     loop over grid points
      do grid_i = 1, grid_r_nb

        u = grid_r (grid_i)

        rpu = r + u
        rmu   = r - u
        abs_rmu = dabs(rmu)
        ru = r * u

        mcabsrmu = -intra_exp*abs_rmu
        mcrpu = -intra_exp*rpu

        exp_mcabsrmu = dexp(mcabsrmu)
        exp_mcrpu = dexp(mcrpu)

        intra_sp_q5_wlk (spin, grid_i, walk_i) = intra_sp_q5_wlk (spin, grid_i, walk_i) + (exp_mcabsrmu - exp_mcrpu)/(2.d0*intra_exp*ru)

     enddo !grid_i

    enddo !elec_j
   enddo !elec_i
  enddo !walk_i

  end subroutine intra_sp_q5_wlk_bld

! ==============================================================================
  subroutine intra_sp_q5_eloc_wlk_bld
! ------------------------------------------------------------------------------
! Description   : term for zero-bias correction for spin-resolved intracule
! Description   : Q * Eloc
!
! Created       : J. Toulouse, 21 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer walk_i

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_q5_eloc_wlk')
   call object_average_walk_define ('intra_sp_q5_eloc_wlk', 'intra_sp_q5_eloc_av')

   call object_needed ('grid_r_nb')
   call object_needed ('nwalk')
   call object_needed ('intra_sp_q5_wlk')
   call object_needed ('eloc_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_sp_q5_eloc_wlk', intra_sp_q5_eloc_wlk, spin_nb, grid_r_nb, nwalk)
  call object_alloc ('intra_sp_q5_eloc_av', intra_sp_q5_eloc_av, spin_nb, grid_r_nb)

  do walk_i = 1, nwalk
   intra_sp_q5_eloc_wlk (:,:,walk_i) = intra_sp_q5_wlk (:,:,walk_i) * eloc_wlk (walk_i)
  enddo

  end subroutine intra_sp_q5_eloc_wlk_bld

! ==============================================================================
  subroutine intra_sp_zb5_av_bld
! ------------------------------------------------------------------------------
! Description   : zero-bias correction for spin-resolved intracule
! Description   : (-1/4pi) [ < Q * Eloc > - < Q > * < Eloc> ]
!
! Created       : J. Toulouse, 21 Nov 2006
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zb5_av')

   call object_needed ('grid_r_nb')
   call object_needed ('intra_sp_q5_eloc_av')
   call object_needed ('intra_sp_q5_av')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zb5_av', intra_sp_zb5_av, spin_nb, grid_r_nb)

  intra_sp_zb5_av (:,:) = -oneover4pi * (intra_sp_q5_eloc_av (:,:) - intra_sp_q5_av (:,:) * eloc_av )

  end subroutine intra_sp_zb5_av_bld

! ==============================================================================
  subroutine intra_sp_zvzb5_av_bld
! ------------------------------------------------------------------------------
! Description   : zero-variance zero-bias correction for spin-resolved intracule
!
! Created       : J. Toulouse, 20 Jul 2006
! ------------------------------------------------------------------------------
  implicit none

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zvzb5_av')
   call object_error_define ('intra_sp_zvzb5_av', 'intra_sp_zvzb5_av_err')

   call object_needed ('grid_r_nb')
   call object_needed ('intra_sp_zv5_av')
   call object_needed ('intra_sp_zb5_av')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zvzb5_av', intra_sp_zvzb5_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zvzb5_av_err', intra_sp_zvzb5_av_err, spin_nb, grid_r_nb)

  intra_sp_zvzb5_av (:,:) = intra_sp_zv5_av (:,:) + intra_sp_zb5_av (:,:)

  end subroutine intra_sp_zvzb5_av_bld

! ==============================================================================
  subroutine intra_sp_zvzb1_bld
! ------------------------------------------------------------------------------
! Description   : first Zero-Variance Zero-Bias improved estimator of spin_resolved intracule
!
! Created       : J. Toulouse, 20 Jul 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: dij, dij2, r
  real(dp)              :: dotproduct, intra_temp
  real(dp)              :: zbfac
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zvzb1')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')

   call object_needed ('eloc')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zvzb1', intra_sp_zvzb1, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zvzb1_av', intra_sp_zvzb1_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zvzb1_av_err', intra_sp_zvzb1_av_err, spin_nb, grid_r_nb)

  intra_sp_zvzb1 (:,:) = 0.d0

  zbfac = -(eloc - eloc_av)*oneover4pi

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      dij = dist_ee_wlk (elec_i, elec_j, 1)

!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_j, 1) * vec_ee_xyz_wlk (dim_i, elec_j, elec_i, 1)
      enddo

      intra_temp = - oneover4pi * dotproduct / dij**3

      do grid_i = 1, grid_r_nb
        r = grid_r(grid_i)

        if ( dij >= r ) then
           intra_sp_zvzb1 (spin, grid_i) = intra_sp_zvzb1 (spin, grid_i) + intra_temp + zbfac*(1.d0/dij)
         else
           intra_sp_zvzb1 (spin, grid_i) = intra_sp_zvzb1 (spin, grid_i) + zbfac*(1.d0/r)
        endif
     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_zvzb1_bld

! ==============================================================================
  subroutine intra_sp_zvzb3_bld
! ------------------------------------------------------------------------------
! Description   : improved zero-variance zer-bias estimator of spin_resolved intracule
! Description   : Q(R) = (-1/4pi) sum_i sim_{j<i} 1/|rij -u| f(rij,u)
! Description   : with f(rij,u) =exp [-c (rij-u)^2]
!
! Created       : J. Toulouse, 03 Aug 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: cc, sqrtc, sqrtpi_c
  real(dp)              :: r, u, rpu, rpu2, ru, rmu, rmu2, abs_rmu
  real(dp)              :: exp_mcrpu2, exp_mcrmu2
  real(dp)              :: dotproduct
  real(dp)              :: intra_temp1, intra_temp2, intra_temp3
  real(dp)              :: zbfac
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zvzb3')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')

   call object_needed ('eloc')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zvzb3', intra_sp_zvzb3, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zvzb3_av', intra_sp_zvzb3_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zvzb3_av_err', intra_sp_zvzb3_av_err, spin_nb, grid_r_nb)

  intra_sp_zvzb3 (:,:) = 0.d0

  cc= intra_exp

  sqrtc  = dsqrt(cc)
  sqrtpi_c = sqrt_pi/sqrtc

  zbfac = -(eloc - eloc_av)

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      r = dist_ee_wlk (elec_i, elec_j, 1)
!
!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_j, 1) * vec_ee_xyz_wlk (dim_i, elec_j, elec_i, 1)
      enddo


!     loop over grid points
      do grid_i = 1, grid_r_nb

        u = grid_r (grid_i)

        rpu = r + u
        rpu2 = rpu * rpu

        rmu   = r - u
        rmu2   = rmu * rmu
        abs_rmu = dabs(rmu)

        ru = r * u

        exp_mcrpu2 = dexp(-cc*rpu2)
        exp_mcrmu2 = dexp(-cc*rmu2)

        intra_temp1 = 0.5d0 * (exp_mcrpu2 - sign(1.d0,rmu) * exp_mcrmu2)
        intra_temp2 = 0.5d0 * (exp_mcrpu2 - (abs_rmu/rpu)* exp_mcrmu2)
        intra_temp3 = sqrtpi_c * (erf_spline(sqrtc*rpu) - erf_spline(sqrtc*abs_rmu))/(4.d0*ru)

        intra_sp_zvzb3 (spin, grid_i) = intra_sp_zvzb3 (spin, grid_i) +          &
            oneover4pi*(((dotproduct/r)*(intra_temp1 - u*intra_temp3) - cc*rpu*intra_temp2)/ru + zbfac*intra_temp3)

     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_zvzb3_bld

! ==============================================================================
  subroutine intra_sp_zvzb4_bld
! ------------------------------------------------------------------------------
! Description   : improved estimator of spin_resolved intracule determined from 1D
!
! Created       : J. Toulouse, 17 Aug 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: r, u, rmu, rmu2
  real(dp)              :: dotproduct
  real(dp)              :: cc
  real(dp)              :: deloc
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zvzb4')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')

   call object_needed ('eloc')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zvzb4', intra_sp_zvzb4, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zvzb4_av', intra_sp_zvzb4_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zvzb4_av_err', intra_sp_zvzb4_av_err, spin_nb, grid_r_nb)

  intra_sp_zvzb4 (:,:) = 0.d0

  cc= intra_exp

  deloc = eloc - eloc_av

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      r = dist_ee_wlk (elec_i, elec_j, 1)

!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_j, 1) * vec_ee_xyz_wlk (dim_i, elec_j, elec_i, 1)
      enddo

      do grid_i = 1, grid_r_nb
        u = grid_r(grid_i)

        if ( r >= u ) then
           intra_sp_zvzb4 (spin, grid_i) = intra_sp_zvzb4 (spin, grid_i) - oneover4pi * (dotproduct/r**3 +deloc/r)
        else
           rmu = r-u
           rmu2 = rmu * rmu
           intra_sp_zvzb4 (spin, grid_i) = intra_sp_zvzb4 (spin, grid_i) - oneover4pi * (cc*                      &
                                (1.d0 + 2.d0*rmu/r - 2.d0*cc*rmu2 + 2.d0*dotproduct*rmu/r) + deloc/u) * dexp(-cc*rmu2)
        endif
     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_zvzb4_bld

! ==============================================================================
  subroutine intra_sp_zv1zb3_bld
! ------------------------------------------------------------------------------
! Description   : Mixed Zero-Variance Zero-Bias improved estimator of spin_resolved intracule
!
! Created       : J. Toulouse, 08 Oct 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: r, u
  real(dp)              :: dotproduct, intra_temp
  real(dp)              :: zbfac
  integer                       :: spin
  real(dp)              :: sqrtpi_c, sqrtc, delta_erf

#ifdef PGI
  real(dp), external    :: derf
#endif

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zv1zb3')

   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')

   call object_needed ('eloc')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zv1zb3', intra_sp_zv1zb3, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv1zb3_av', intra_sp_zv1zb3_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zv1zb3_av_err', intra_sp_zv1zb3_av_err, spin_nb, grid_r_nb)

  sqrtc  = dsqrt(intra_exp)
  sqrtpi_c = sqrt_pi/sqrtc

  zbfac = -(eloc - eloc_av)*oneover4pi

  intra_sp_zv1zb3 (:,:) = 0.d0

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, 1)

!     distance e-e
      r = dist_ee_wlk (elec_i, elec_j, 1)

!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_j, 1) * vec_ee_xyz_wlk (dim_i, elec_j, elec_i, 1)
      enddo

      intra_temp = - oneover4pi * dotproduct / r**3

      do grid_i = 1, grid_r_nb
        u = grid_r(grid_i)

! (1/4pi) int dOmega_u 1/|rij -u| exp[-c(rij-u)^2]
        delta_erf = sqrtpi_c*(derf(sqrtc*(r+u)) - derf(sqrtc*dabs(r-u)))/(4.d0*r*u)

        if ( r >= u ) then
           intra_sp_zv1zb3 (spin, grid_i) = intra_sp_zv1zb3 (spin, grid_i) + intra_temp + zbfac*delta_erf
         else
           intra_sp_zv1zb3 (spin, grid_i) = intra_sp_zv1zb3 (spin, grid_i) + zbfac*delta_erf
        endif
     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_sp_zv1zb3_bld

! ==============================================================================
  subroutine intra_sp_zvzb5_bld
! ------------------------------------------------------------------------------
! Description   : improved zero-variance zero-bias estimator of spin_resolved intracule
! Description   : Q(R) = (-1/4pi) sum_i sim_{j<i} 1/|rij -u| f(rij,u)
! Description   : with f(rij,u) =exp [-c |rij-u|]
! Description   : eloc_av comes from previous iterations
!
! Created       : J. Toulouse, 10 Nov 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i, walk_i
  real(dp)              :: cc, cc2
  real(dp)              :: r, r2, u, rpu, ru, rmu, abs_rmu
  real(dp)              :: mcabsrmu, mcrpu, exp_mcabsrmu, exp_mcrpu
  real(dp)              :: dotproduct
  real(dp)              :: intra_temp1, intra_temp2
  real(dp)              :: zbfac
  integer                       :: spin

! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_zvzb5')
   call object_average_walk_define ('intra_sp_zvzb5', 'intra_sp_zvzb5_av')
   call object_error_define ('intra_sp_zvzb5_av', 'intra_sp_zvzb5_av_err')

   call object_needed ('nwalk')
   call object_needed ('grid_r')
   call object_needed ('grid_r_nb')
   call object_needed ('nelec')
   call object_needed ('dist_ee_wlk')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('grd_psi_over_psi_wlk')
   call object_needed ('spin_two_elec_wlk')
   call object_needed ('eloc_wlk')
   call object_needed ('eloc_av')

   return

  endif

! allocations
  call object_alloc ('intra_sp_zvzb5', intra_sp_zvzb5, spin_nb, grid_r_nb, nwalk)
  call object_alloc ('intra_sp_zvzb5_av', intra_sp_zvzb5_av, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_zvzb5_av_err', intra_sp_zvzb5_av_err, spin_nb, grid_r_nb)

  intra_sp_zvzb5 (:,:,:) = 0.d0

  cc = intra_exp
  cc2 = cc*cc


  do walk_i = 1, nwalk

   zbfac = eloc_wlk (walk_i) - eloc_av

   do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

!     same or opposite spin?
      spin = spin_two_elec_wlk (elec_i, elec_j, walk_i)

!     distance e-e
      r = dist_ee_wlk (elec_i, elec_j, walk_i)
      r2 = r*r
!
!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  grd_psi_over_psi_wlk (dim_i, elec_j, walk_i) * vec_ee_xyz_wlk (dim_i, elec_j, elec_i, walk_i)
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

        intra_temp1 = (exp_mcabsrmu - exp_mcrpu)/(2.d0*cc*ru)
        intra_temp2 = 0.5d0 * (sign(1.d0,rmu) * exp_mcabsrmu - exp_mcrpu)

        intra_sp_zvzb5 (spin, grid_i, walk_i) = intra_sp_zvzb5 (spin, grid_i, walk_i) -          &
            oneover4pi*((dotproduct/r2)*(intra_temp1 + intra_temp2/u)  &
            - 0.5d0*cc2*intra_temp1 + zbfac*intra_temp1)

     enddo !grid_i

    enddo !elec_j
   enddo !elec_i
  enddo !walk_i

  end subroutine intra_sp_zvzb5_bld

! ==============================================================================
  subroutine intra_sp_bld
! ------------------------------------------------------------------------------
! Description   : spin-resolved intracule densities  I(r)
!
! Created       : J. Toulouse, 05 Oct 2005
! Revised       : J. Toulouse, 20 Mar 2006, spin-resolved
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'intra_sp_bld'

! header
  if (header_exe) then

   call object_create ('intra_sp')
   call object_create ('intra_sp_err')

   call object_needed ('grid_r_nb')

   return

  endif

! begin

! allocations
  call object_alloc ('intra_sp', intra_sp, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_err', intra_sp_err, spin_nb, grid_r_nb)

  select case(trim(intra_estimator))
   case ('histogram')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_histo_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_histo_av_err_index)
   intra_sp (:,:)     = intra_sp_histo_av (:,:)
   intra_sp_err (:,:) = intra_sp_histo_av_err (:,:)

   case ('zv1')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv1_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv1_av_err_index)
   intra_sp (:,:)     = intra_sp_zv1_av (:,:)
   intra_sp_err (:,:) = intra_sp_zv1_av_err (:,:)

   case ('zv2')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv2_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv2_av_err_index)
   intra_sp (:,:)     = intra_sp_zv2_av (:,:)
   intra_sp_err (:,:) = intra_sp_zv2_av_err (:,:)

   case ('zv3')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv3_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv3_av_err_index)
   intra_sp (:,:)     = intra_sp_zv3_av (:,:)
   intra_sp_err (:,:) = intra_sp_zv3_av_err (:,:)

   case ('zv4')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv4_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv4_av_err_index)
   intra_sp (:,:)     = intra_sp_zv4_av (:,:)
   intra_sp_err (:,:) = intra_sp_zv4_av_err (:,:)

   case ('zv5')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv5_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv5_av_err_index)
   intra_sp (:,:)     = intra_sp_zv5_av (:,:)
   intra_sp_err (:,:) = intra_sp_zv5_av_err (:,:)

   case ('zvzb1')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zvzb1_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zvzb1_av_err_index)
   intra_sp (:,:)     = intra_sp_zvzb1_av (:,:)
   intra_sp_err (:,:) = intra_sp_zvzb1_av_err (:,:)

   case ('zvzb3')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zvzb3_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zvzb3_av_err_index)
   intra_sp (:,:)     = intra_sp_zvzb3_av (:,:)
   intra_sp_err (:,:) = intra_sp_zvzb3_av_err (:,:)

   case ('zvzb4')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zvzb4_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zvzb4_av_err_index)
   intra_sp (:,:)     = intra_sp_zvzb4_av (:,:)
   intra_sp_err (:,:) = intra_sp_zvzb4_av_err (:,:)

   case ('zv1zb3')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv1zb3_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zv1zb3_av_err_index)
   intra_sp (:,:)     = intra_sp_zv1zb3_av (:,:)
   intra_sp_err (:,:) = intra_sp_zv1zb3_av_err (:,:)

   case ('zvzb5')
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zvzb5_av_index)
   call object_provide_by_index (intra_sp_bld_index, intra_sp_zvzb5_av_err_index)
   intra_sp (:,:)     = intra_sp_zvzb5_av (:,:)
   intra_sp_err (:,:) = intra_sp_zvzb5_av_err (:,:)

   case default
   write(6,*) trim(lhere), ': unknown estimator >',trim(intra_estimator),'<'
   call die (lhere)
  end select

  end subroutine intra_sp_bld

! ==============================================================================
  subroutine intra_sp_4pir2_bld
! ------------------------------------------------------------------------------
! Description   : 4 pi r^2 I(r), spin-resolved
!
! Created       : J. Toulouse, 05 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  integer                       :: grid_i
  real(dp)              :: volume_elt


! begin

! header
  if (header_exe) then

   call object_create ('intra_sp_4pir2')
   call object_create ('intra_sp_4pir2_err')

   call object_needed ('grid_r_nb')
   call object_needed ('grid_r')
   call object_needed ('intra_sp')

   return

  endif

! allocations
  call object_alloc ('intra_sp_4pir2', intra_sp_4pir2, spin_nb, grid_r_nb)
  call object_alloc ('intra_sp_4pir2_err', intra_sp_4pir2_err, spin_nb, grid_r_nb)

! average and statistical error of intracule f(r) 4 pi r^2
  do grid_i = 1, grid_r_nb
   volume_elt = 4.d0*pi1*grid_r (grid_i)*grid_r(grid_i)
   intra_sp_4pir2 (:,grid_i)   = intra_sp (:,grid_i)*volume_elt
   intra_sp_4pir2_err (:,grid_i)  = intra_sp_err (:,grid_i)*volume_elt
  enddo

  end subroutine intra_sp_4pir2_bld

! ==============================================================================
  subroutine intra_bld
! ------------------------------------------------------------------------------
! Description   : intracule density  I(r) = I_{same spins} (r) +I_{opposite spins} (r)
!
! Created       : J. Toulouse, 05 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'intra_bld'
  integer :: spin

! header
  if (header_exe) then

   call object_create ('intra')
   call object_create ('intra_err')

   call object_needed ('intra_sp')
   call object_needed ('intra_sp_err')
   call object_needed ('grid_r_nb')

   return

  endif

! begin

! allocations
  call object_alloc ('intra', intra, grid_r_nb)
  call object_alloc ('intra_err', intra_err, grid_r_nb)

  intra (:) = 0.d0
  intra_err (:) = 0.d0

!  intra (:) = sum(intra_sp (spin,:), spin = 1, spin_nb)
!  intra_err (:) = sum(intra_sp_err (spin,:), spin = 1, spin_nb)

  do spin = 1, spin_nb
   intra (:) = intra (:) + intra_sp (spin,:)
   intra_err (:) = intra_err (:) + intra_sp_err (spin,:)
  enddo

  end subroutine intra_bld

! ==============================================================================
  subroutine intra_4pir2_bld
! ------------------------------------------------------------------------------
! Description   : 4 pi r^2 I(r)
!
! Created       : J. Toulouse, 05 Oct 2005
! ------------------------------------------------------------------------------
  implicit none

! local
  integer                       :: grid_i
  real(dp)              :: volume_elt


! begin

! header
  if (header_exe) then

   call object_create ('intra_4pir2')
   call object_create ('intra_4pir2_err')

   call object_needed ('grid_r_nb')
   call object_needed ('grid_r')
   call object_needed ('intra')

   return

  endif

! allocations
  call object_alloc ('intra_4pir2', intra_4pir2, grid_r_nb)
  call object_alloc ('intra_4pir2_err', intra_4pir2_err, grid_r_nb)

! average and statistical error of intracule f(r) 4 pi r^2
  do grid_i = 1, grid_r_nb
   volume_elt = 4.d0*pi1*grid_r (grid_i)*grid_r(grid_i)
   intra_4pir2 (grid_i)   = intra (grid_i)*volume_elt
   intra_4pir2_err (grid_i)  = intra_err (grid_i)*volume_elt
  enddo

  end subroutine intra_4pir2_bld

! ==============================================================================
  subroutine intra_3d_histo_bld
! ------------------------------------------------------------------------------
! Description   : histo estimator of 3D intracule
!
! Created       : J. Toulouse, 05 Mar 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: grid_x_i, grid_y_i, grid_z_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: xij, yij, zij

! begin

! header
  if (header_exe) then

   call object_create ('intra_3d_histo')

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
  call object_alloc ('intra_3d_histo', intra_3d_histo, grid_xyz_nb)
  call object_alloc ('intra_3d_histo_av', intra_3d_histo_av, grid_xyz_nb)
  call object_alloc ('intra_3d_histo_av_err', intra_3d_histo_av_err, grid_xyz_nb)

  intra_3d_histo (:) = 0.d0

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

      xij =  xold(1,elec_i) - xold(1,elec_j)
      yij =  xold(2,elec_i) - xold(2,elec_j)
      zij =  xold(3,elec_i) - xold(3,elec_j)

      grid_x_i = floor((xij - grid_x_min)/grid_x_step - 0.5d0)  + 2
      grid_y_i = floor((yij - grid_y_min)/grid_y_step - 0.5d0)  + 2
      grid_z_i = floor((zij - grid_z_min)/grid_z_step - 0.5d0)  + 2

      if (grid_x_i < 1 .or. grid_x_i > grid_x_nb) cycle
      if (grid_y_i < 1 .or. grid_y_i > grid_y_nb) cycle
      if (grid_z_i < 1 .or. grid_z_i > grid_z_nb) cycle


      grid_i = grid_xyz_index (grid_x_i, grid_y_i, grid_z_i)

      intra_3d_histo (grid_i) = intra_3d_histo (grid_i) + 0.5d0/(grid_x_step*grid_y_step*grid_z_step)

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_3d_histo_bld

! ==============================================================================
  subroutine intra_3d_zv1_bld
! ------------------------------------------------------------------------------
! Description   : first-order renormalized imp estimator of 3D intracule
!
! Created       : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: dij, dij2, r
  real(dp)              :: dotproduct, intra_temp

! begin

! header
  if (header_exe) then

   call object_create ('intra_3d_zv1')

   call object_needed ('grid_xyz_nb')
   call object_needed ('grid_xyz')
   call object_needed ('nelec')
   call object_needed ('xold')
   call object_needed ('vold')

   return

  endif

! allocations
  call object_alloc ('intra_3d_zv1', intra_3d_zv1, grid_xyz_nb)
  call object_alloc ('intra_3d_zv1_av', intra_3d_zv1_av, grid_xyz_nb)
  call object_alloc ('intra_3d_zv1_av_err', intra_3d_zv1_av_err, grid_xyz_nb)

  intra_3d_zv1 (:) = 0.d0

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

      do grid_i = 1, grid_xyz_nb

!     (rij - r)
      dij2 = 0.d0
      do dim_i = 1,ndim
         dij2 = dij2 + (xold(dim_i,elec_i) - xold(dim_i,elec_j)  - grid_xyz (dim_i, grid_i))**2
      enddo ! dim_i
      dij = dsqrt(dij2)

!     dot product: drift_j . (rij -r)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  vold(dim_i,elec_j) * (xold(dim_i,elec_j) - xold(dim_i,elec_i) - grid_xyz (dim_i, grid_i))
      enddo

      intra_temp = - oneover4pi * dotproduct / dij**3

      intra_3d_zv1 (grid_i) = intra_3d_zv1 (grid_i) + intra_temp

     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_3d_zv1_bld

! ==============================================================================
  subroutine intra_3d_zv2_bld
! ------------------------------------------------------------------------------
! Description   : second-order renormalized improved estimator of 3D intracule
!
! Created       : J. Toulouse, 07 Mar 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: dij, dij2, r
  real(dp)              :: dotproduct, intra_temp

! begin

! header
  if (header_exe) then

   call object_create ('intra_3d_zv2')

   call object_needed ('grid_xyz_nb')
   call object_needed ('grid_xyz')
   call object_needed ('nelec')
   call object_needed ('vec_ee_xyz_wlk')
   call object_needed ('lap_psi_over_psi_wlk')
   call object_needed ('grd_psi_over_psi_sq_wlk')

   return

  endif

! allocations
  call object_alloc ('intra_3d_zv2', intra_3d_zv2, grid_xyz_nb)
  call object_alloc ('intra_3d_zv2_av', intra_3d_zv2_av, grid_xyz_nb)
  call object_alloc ('intra_3d_zv2_av_err', intra_3d_zv2_av_err, grid_xyz_nb)

  intra_3d_zv2 (:) = 0.d0

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle

      do grid_i = 1, grid_xyz_nb

!     (rij - r)
      dij2 = 0.d0
      do dim_i = 1,ndim
         dij2 = dij2 + (vec_ee_xyz_wlk (dim_i, elec_i, elec_j, 1)  - grid_xyz (dim_i, grid_i))**2
      enddo ! dim_i
      dij = dsqrt(dij2)

      intra_temp = - oneover4pi * (lap_psi_over_psi_wlk (elec_j, 1) + grd_psi_over_psi_sq_wlk (elec_j, 1)) / dij

      intra_3d_zv2 (grid_i) = intra_3d_zv2 (grid_i) + intra_temp

     enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_3d_zv2_bld

! ==============================================================================
  subroutine intra_3d_bld
! ------------------------------------------------------------------------------
! Description   : 3D intracule density  I(r)
!
! Created       : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'intra_3d_bld'

! header
  if (header_exe) then

   call object_create ('intra_3d')
   call object_create ('intra_3d_err')

   call object_needed ('grid_xyz_nb')

   return

  endif

! begin

! allocations
  call object_alloc ('intra_3d', intra_3d, grid_xyz_nb)
  call object_alloc ('intra_3d_err', intra_3d_err, grid_xyz_nb)

  select case(trim(intra_3d_estimator))
   case ('histogram')
   call object_provide (lhere,'intra_3d_histo_av')
   call object_provide (lhere,'intra_3d_histo_av_err')
   intra_3d (:)     = intra_3d_histo_av (:)
   intra_3d_err (:) = intra_3d_histo_av_err (:)

   case ('zv1')
   call object_provide (lhere,'intra_3d_zv1_av')
   call object_provide (lhere,'intra_3d_zv1_av_err')
   intra_3d (:)     = intra_3d_zv1_av (:)
   intra_3d_err (:) = intra_3d_zv1_av_err (:)

   case ('zv2')
   call object_provide (lhere,'intra_3d_zv2_av')
   call object_provide (lhere,'intra_3d_zv2_av_err')
   intra_3d (:)     = intra_3d_zv2_av (:)
   intra_3d_err (:) = intra_3d_zv2_av_err (:)

   case default
   write(6,*) trim(lhere), ': unknown estimator >',trim(intra_3d_estimator),'<'
   call die (lhere)
  end select

  end subroutine intra_3d_bld

! ==============================================================================
  subroutine intra_zv1_step_vmc
! ------------------------------------------------------------------------------
! Description   : Calculate first-order renormalized imp estimator of intracule
! Description   : at each step of the VMC dynamics
!
! Created       : J. Toulouse, 05 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout) :: here
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i
  real(dp)              :: dij, dij2, r
  real(dp)              :: dotproduct, intra_temp


! begin
  here = 'intra_zv1_step_vmc'
!  write(6,*) trim(here),': entering'

! objects needed
  if(.not.allocated(grid_r)) then
    call grid_r_bld
  endif

! allocations
  call alloc ('intra_imp_num', intra_imp_num, grid_r_nb)

! initializations
  if(init_intra_imp_step) then
   intra_imp_num = 0.d0
   intra_imp_nb = 0
   init_intra_imp_step = .false.
  endif

  intra_imp_nb = intra_imp_nb + 1

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle


!     distance e-e
      dij2 = 0.d0
      do dim_i = 1,ndim
         dij2 = dij2 + (xold(dim_i,elec_i) - xold(dim_i,elec_j))**2
      enddo ! dim_i
      dij = dsqrt(dij2)

!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +  vold(dim_i,elec_j) * (xold(dim_i,elec_j) - xold(dim_i,elec_i))
      enddo

      intra_temp = oneover4pi * dotproduct / dij**3

      do grid_i = 1, grid_r_nb
         r = grid_r(grid_i)
         if ( dij >= r ) then
           intra_imp_num(grid_i) = intra_imp_num(grid_i) - intra_temp
         endif
      enddo !grid_i

    enddo !elec_j
  enddo !elec_i

  end subroutine intra_zv1_step_vmc

! ==============================================================================
  subroutine intra_zv1_step_dmc
! ------------------------------------------------------------------------------
! Description   : Calculate first-order renormalized imp estimator of intracule
! Description   : at each step of the DMC dynamics
!
! Created       : J. Toulouse, 06 Oct 2005
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout) :: here
  integer                       :: grid_i
  integer                       :: elec_i, elec_j, dim_i, walk_i
  real(dp)              :: dij, dij2, r
  real(dp)              :: dotproduct, intra_temp
  real(dp)              :: weight


! begin
  here = 'intra_zv1_step_dmc'
!  write(6,*) trim(here),': entering'

! objects needed
  if(.not.allocated(grid_r)) then
    call grid_r_bld
  endif

! allocations
  call alloc ('intra_imp_num', intra_imp_num, grid_r_nb)

! initializations
  if(init_intra_imp_step) then
   intra_imp_num = 0.d0
   intra_imp_nb = 0
   init_intra_imp_step = .false.
  endif

  intra_imp_nb = intra_imp_nb + 1

  do walk_i =1, nwalk

   weight = wt(walk_i)*fprod

!  write(6,*) trim(here), 'walk_i=',walk_i, 'wt =',wt(walk_i)

  do elec_j = 1, nelec
    do elec_i = 1, nelec

      if ( elec_i == elec_j ) cycle


!     distance e-e
      dij2 = 0.d0
      do dim_i = 1,ndim
         dij2 = dij2 + (xoldw(dim_i,elec_i,walk_i,1) - xoldw(dim_i,elec_j,walk_i,1))**2
      enddo ! dim_i
      dij = dsqrt(dij2)

!     dot product: drift_j . (rj - ri)
      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct +   &
            voldw(dim_i,elec_j,walk_i,1) * (xoldw(dim_i,elec_j,walk_i,1) - xoldw(dim_i,elec_i,walk_i,1))
      enddo

      intra_temp = oneover4pi * dotproduct / dij**3

      do grid_i = 1, grid_r_nb
         r = grid_r(grid_i)
         if ( dij >= r ) then
           intra_imp_num(grid_i) = intra_imp_num(grid_i) - intra_temp * weight
         endif
      enddo !grid_i

    enddo !elec_j
  enddo !elec_i

 enddo ! walk_i

  end subroutine intra_zv1_step_dmc

! ========================================================================
  subroutine intra_imp_block_vmc
! ------------------------------------------------------------------------
! Description    : compute intracule at each block of MC dynamics
! Description    : f(r) 4 pi r^2 : intracule with volume element
! Description    : f(r)          : intracule without volume element
! Description    : using imp estimators
!
! Created        : J. Toulouse, 06 Oct 2005
! ------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout)         :: here
  integer                                    :: grid_i
  real(dp)                           :: volume_elt

! begin
  here = "intra_imp_block_vmc"


! allocations
  call alloc('sum_intra_imp_num', sum_intra_imp_num, grid_r_nb)
  call alloc('intra_imp_norm', intra_imp_norm, grid_r_nb)
  call alloc('sum_intra_imp_norm', sum_intra_imp_norm, grid_r_nb)
  call alloc('sum_intra_imp_norm_square', sum_intra_imp_norm_square, grid_r_nb)
  call alloc('intra_imp_f_av', intra_imp_f_av, grid_r_nb)
  call alloc('intra_imp_f_var', intra_imp_f_var, grid_r_nb)
  call alloc('intra_imp_f_av_err', intra_imp_f_av_err, grid_r_nb)
  call alloc('intra_imp_f4pir2_av', intra_imp_f4pir2_av, grid_r_nb)
  call alloc('intra_imp_f4pir2_av_err', intra_imp_f4pir2_av_err, grid_r_nb)

! initializations
  if(init_intra_imp_block) then
   intra_block_nb = 0
   sum_intra_imp_num = 0.d0
   sum_intra_imp_norm = 0.d0
   sum_intra_imp_norm_square = 0.d0
   init_intra_imp_block = .false.
  endif

! block number
  intra_block_nb = intra_block_nb + 1

! intracule normalized to N(N-1)/2
  intra_imp_norm = intra_imp_num/intra_imp_nb

! sums over the blocks of intra_imp_norm and square of intra_imp_norm
  sum_intra_imp_norm = sum_intra_imp_norm + intra_imp_norm
  sum_intra_imp_norm_square = sum_intra_imp_norm_square + intra_imp_norm**2

! average of intra_imp_norm = average of intracule f(r)
   intra_imp_f_av =  sum_intra_imp_norm/intra_block_nb

! variance of the average of intracule f(r)
  intra_imp_f_var = sum_intra_imp_norm_square/intra_block_nb - intra_imp_f_av**2

! statistical error on the average of intracule f(r)
   if (intra_block_nb == 1 ) then
    intra_imp_f_av_err = 0.d0
   else
    intra_imp_f_av_err = dsqrt(intra_imp_f_var/(intra_block_nb-1))
   endif

! average and statistical error of intracule f(r) 4 pi r^2
  do grid_i = 1, grid_r_nb
   volume_elt = 4.d0*pi1*grid_r(grid_i)*grid_r(grid_i)
   intra_imp_f4pir2_av(grid_i)   = intra_imp_f_av(grid_i)*volume_elt
   intra_imp_f4pir2_av_err(grid_i)  = intra_imp_f_av_err(grid_i)*volume_elt
  enddo

! Reinitialization of intra_step
  init_intra_imp_step = .true.

! Write intracule
  call intra_imp_write

  end subroutine intra_imp_block_vmc

! ========================================================================
  subroutine intra_imp_block_dmc
! ------------------------------------------------------------------------
! Description    : compute intracule at each block of MC dynamics
! Description    : f(r) 4 pi r^2 : intracule with volume element
! Description    : f(r)          : intracule without volume element
! Description    : using imp estimators
!
! Created        : J. Toulouse, 06 Oct 2005
! ------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout)         :: here
  integer                                    :: grid_i
  real(dp)                           :: volume_elt

! begin
  here = "intra_imp_block_dmc"


! allocations
  call alloc('sum_intra_imp_num', sum_intra_imp_num, grid_r_nb)
  call alloc('intra_imp_norm', intra_imp_norm, grid_r_nb)
  call alloc('sum_intra_imp_norm', sum_intra_imp_norm, grid_r_nb)
  call alloc('sum_intra_imp_norm_square', sum_intra_imp_norm_square, grid_r_nb)
  call alloc('intra_imp_f_av', intra_imp_f_av, grid_r_nb)
  call alloc('intra_imp_f_var', intra_imp_f_var, grid_r_nb)
  call alloc('intra_imp_f_av_err', intra_imp_f_av_err, grid_r_nb)
  call alloc('intra_imp_f4pir2_av', intra_imp_f4pir2_av, grid_r_nb)
  call alloc('intra_imp_f4pir2_av_err', intra_imp_f4pir2_av_err, grid_r_nb)

! initializations
  if(init_intra_imp_block) then
   intra_block_nb = 0
   sum_intra_imp_num = 0.d0
   sum_intra_imp_norm = 0.d0
   sum_intra_imp_norm_square = 0.d0
   init_intra_imp_block = .false.
  endif

! block number
  intra_block_nb = intra_block_nb + 1

! intracule normalized to N(N-1)/2
! divide by the sum of weigths * number of steps
!   write(6,*) trim(here), ': wgsum=',wgsum(1)
   intra_imp_norm = intra_imp_num/wgsum(1)

! sums over the blocks of intra_imp_norm and square of intra_imp_norm
  sum_intra_imp_norm = sum_intra_imp_norm + intra_imp_norm
  sum_intra_imp_norm_square = sum_intra_imp_norm_square + intra_imp_norm**2

! average of intra_imp_norm = average of intracule f(r)
   sum_intra_imp_num = sum_intra_imp_num + intra_imp_num
   intra_imp_f_av =  sum_intra_imp_num/wgcum(1)

! variance of the average of intracule f(r)
  intra_imp_f_var = sum_intra_imp_norm_square/intra_block_nb - (sum_intra_imp_norm/intra_block_nb)**2

! statistical error on the average of intracule f(r)
   if (intra_block_nb == 1 ) then
    intra_imp_f_av_err = 0.d0
   else
    intra_imp_f_av_err = dsqrt(intra_imp_f_var/(intra_block_nb-1))
   endif

! average and statistical error of intracule f(r) 4 pi r^2
  do grid_i = 1, grid_r_nb
   volume_elt = 4.d0*pi*grid_r(grid_i)*grid_r(grid_i)
   intra_imp_f4pir2_av(grid_i)   = intra_imp_f_av(grid_i)*volume_elt
   intra_imp_f4pir2_av_err(grid_i)  = intra_imp_f_av_err(grid_i)*volume_elt
  enddo


! Reinitialization of intra_step
  init_intra_imp_step = .true.

! Write intracule
  call intra_imp_write

  end subroutine intra_imp_block_dmc

! ========================================================================
  subroutine intra_imp_write
! ------------------------------------------------------------------------
! Description    : write intracule density on file
! Description    : f(r) 4 pi r^2 : intracule with volume element
! Description    : f(r)          : intracule without volume element
! Description    : using imp estimator
!
! Created        : J. Toulouse, 06 Oct 2005
! ------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save   :: lhere
  integer                                    :: unit
  character(len=max_string_len_file)         :: file
  integer                                    :: grid_i

! begin
  here ='intra_imp_write'

! open file
!  write(6,*) trim(lhere),': entering'
  file = intra_file_out
  unit = 20
  open(file=trim(file),unit=unit)

  write(unit,'(a)')         'intracule generated by CHAMP'
  write(unit,'(a,i5)')      'number of electrons       =',nelec
  write(unit,'(a,i20)')      'number of steps per block =',nstep
  write(unit,'(a,i20)')      'number of blocks          =',intra_block_nb
  write(unit,'(a,3e25.15)') 'ee distance step          =',grid_r_step
  write(unit,'(a,3e25.15)') 'ee max distance           =',grid_r_max
  write(unit,'(a,i25)')     'number of grid points     =',grid_r_nb

  write(unit,*) ''
  write(unit,'(a)') '             r                     f4pir2_av              f4pir2_av_err               f_av                    f_av_err            sum_intra_norm    sum_intra_norm_square'

  do grid_i = 1, grid_r_nb

  write(unit,'(3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,3e25.15)') grid_r(grid_i), &
        intra_imp_f4pir2_av(grid_i),intra_imp_f4pir2_av_err(grid_i),intra_imp_f_av(grid_i), &
        intra_imp_f_av_err(grid_i),sum_intra_imp_norm(grid_i),     &
        sum_intra_imp_norm_square(grid_i)


  enddo


  close(unit)

!  write(6,*) trim(lhere),': intracule has been written on file: >',trim(file),'<'

  end subroutine intra_imp_write


! ========================================================================
  subroutine intra_wrt
! ------------------------------------------------------------------------
! Description    : write intracule density on file
! Description    : f(r) 4 pi r^2 : intracule with volume element
! Description    : f(r)          : intracule without volume element
!
! Created        : J. Toulouse, 06 Oct 2005
! ------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save   :: lhere = 'intra_wrt'
  integer                                    :: unit
  character(len=max_string_len_file)         :: file
  integer                                    :: grid_i

! begin
  call object_provide ('grid_r')
  call object_provide ('intra')
  call object_provide ('intra_err')
  call object_provide ('intra_4pir2')
  call object_provide ('intra_4pir2_err')
  call object_provide ('intra_sp')
  call object_provide ('intra_sp_err')
  call object_provide ('intra_sp_4pir2')
  call object_provide ('intra_sp_4pir2_err')
  call object_provide ('dist_ee_min')
  call object_provide ('dist_ee_max')


! open file
  file = intra_file_out
  unit = 20
  open(file=trim(file),unit=unit)

  write(unit,'(a)')         'intracule generated by CHAMP'
  write(unit,'(a,i25)')     'number of electrons           =',nelec
  write(unit,'(a,i25)')     'number of spin-up electrons   =',nup
  write(unit,'(a,i25)')     'number of spin-down electrons =',ndn
  write(unit,'(a,i25)')     'number of steps per block     =',nstep_total
  write(unit,'(a,i25)')     'number of blocks              =',block_iterations_nb
  write(unit,'(a,3e25.15)') 'grid_r_step                   =',grid_r_step
  write(unit,'(a,3e25.15)') 'grid_r_min                    =',grid_r_min
  write(unit,'(a,3e25.15)') 'grid_r_max                    =',grid_r_max
  write(unit,'(a,i25)')     'grid_r_nb                     =',grid_r_nb
  write(unit,'(a,f)')       'dist_ee_min                   =',dist_ee_min
  write(unit,'(a,f)')       'dist_ee_max                   =',dist_ee_max

  write(unit,*) ''
  write(unit,'(a)') '             r                   4 pi r2 I(r)               error                    I(r)                     error           4 pi r2 I(r) same spins&
&          error               I(r) same spins               error          4 pi r2 I(r) opposite spins       error              I(r) opposite spins            error'

  do grid_i = 1, grid_r_nb

  write(unit,'(3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,3e25.15,10e25.15)') grid_r(grid_i), &
        intra_4pir2(grid_i),intra_4pir2_err(grid_i),intra(grid_i), intra_err(grid_i),             &
        intra_sp_4pir2(1,grid_i),intra_sp_4pir2_err(1,grid_i),intra_sp(1,grid_i), intra_sp_err(1,grid_i),      &
        intra_sp_4pir2(2,grid_i),intra_sp_4pir2_err(2,grid_i),intra_sp(2,grid_i), intra_sp_err(2,grid_i)
  enddo


  close(unit)

  end subroutine intra_wrt

! ========================================================================
  subroutine intra_3d_wrt
! ------------------------------------------------------------------------
! Description    : write 3D intracule density on file
!
! Created        : J. Toulouse, 04 Mar 2006
! ------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save   :: lhere = 'intra_3d_wrt'
  integer                                    :: unit
  character(len=max_string_len_file)         :: file
  integer                                    :: grid_i, grid_x_i, grid_y_i, grid_z_i

! begin
  call object_provide ('grid_xyz_nb')
  call object_provide ('grid_xyz')
  call object_provide ('intra_3d')


! open file
  file = intra_3d_file_out
  unit = 21
  open(file=trim(file),unit=unit)

  write(unit,'(a)')         '3D intracule generated by CHAMP'
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
  write(unit,'(a)') '             x                        y                        z                       I(r)                   error I(r)'

  do grid_x_i = 1, grid_x_nb
   do grid_y_i = 1, grid_y_nb
    do grid_z_i = 1, grid_z_nb

       grid_i = grid_xyz_index (grid_x_i, grid_y_i, grid_z_i)

       write(unit,'(5e25.15)') grid_xyz(1,grid_i), grid_xyz(2,grid_i), grid_xyz(3,grid_i), intra_3d (grid_i), intra_3d_err(grid_i)
    enddo ! grid_z_i
   enddo  ! grid_y_i
       write(unit,'(a)') ''
  enddo ! grid_z_i

  close(unit)

  end subroutine intra_3d_wrt

end module intracule_mod
