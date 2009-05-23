module electrons_mod

  use all_tools_mod
  use walkers_mod

! Declaration of global variables and default values
  real(dp), allocatable     :: coord_elec_wlk (:,:,:)

  real(dp), allocatable     :: dist_e_wlk (:,:)
  real(dp)                  :: dist_e_min = 1000.d0
  real(dp)                  :: dist_e_max = 0.d0

  real(dp), allocatable     :: vec_ee_xyz_wlk (:,:,:,:)
  real(dp), allocatable     :: dist_ee_wlk (:,:,:)
  real(dp)                  :: dist_ee_min = 1000.d0
  real(dp)                  :: dist_ee_max = 0.d0

  real(dp), allocatable     :: vec_en_xyz_wlk (:,:,:,:)
  real(dp), allocatable     :: dist_en_wlk (:,:,:)
  real(dp), allocatable     :: grd_dist_en (:,:,:)
  real(dp), allocatable     :: lap_dist_en (:,:)

  real(dp), allocatable     :: elec_nb_close_atom_wlk (:,:)
  real(dp), allocatable     :: elec_spin_nb_close_atom_wlk (:,:,:)
  real(dp), allocatable     :: elec_up_nb_close_atom_wlk (:,:)
  real(dp), allocatable     :: elec_dn_nb_close_atom_wlk (:,:)
  real(dp), allocatable     :: elec_nb_close_atom_av (:)
  real(dp), allocatable     :: elec_spin_nb_close_atom_av (:,:)
  real(dp), allocatable     :: elec_up_nb_close_atom_av (:)
  real(dp), allocatable     :: elec_dn_nb_close_atom_av (:)
  real(dp), allocatable     :: frac_elec_nb_close_atom_input_wlk (:)
  real(dp), allocatable     :: frac_elec_spin_nb_close_atom_input_wlk (:)
  real(dp)                  :: frac_elec_nb_close_atom_input_av
  real(dp)                  :: frac_elec_spin_nb_close_atom_input_av

  integer                   :: elec_nb_close_atom_encount_nb = 0
  real(dp), allocatable     :: elec_nb_close_atom_encount (:,:)
  real(dp), allocatable     :: frac_elec_nb_close_atom_encount_wlk (:,:)
  real(dp), allocatable     :: frac_elec_nb_close_atom_encount_av (:)
  integer                   :: elec_spin_nb_close_atom_encount_nb = 0
  real(dp), allocatable     :: elec_spin_nb_close_atom_encount (:,:,:)
  real(dp), allocatable     :: frac_elec_spin_nb_close_atom_encount_wlk (:,:)
  real(dp), allocatable     :: frac_elec_spin_nb_close_atom_encount_av (:)

  contains

! ==============================================================================
  subroutine coord_elec_wlk_bld
! ------------------------------------------------------------------------------
! Description   : Cartesian coordinates of electrons
!
! Created       : J. Toulouse, 16 Nov 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! header
  if (header_exe) then

   call object_create ('coord_elec_wlk')

   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('nwalk')

   return

  endif

! begin

! allocations
  call object_alloc ('coord_elec_wlk', coord_elec_wlk, ndim, nelec, nwalk)

  if (index(mode, 'vmc') /= 0) then
   call object_provide_by_index (coord_elec_wlk_bld_index, xold_index)
   coord_elec_wlk (:,:,1) = xold (1:ndim, 1:nelec)

  elseif (index(mode, 'dmc') /= 0) then
   call object_provide_by_index (coord_elec_wlk_bld_index, xoldw_index)
   coord_elec_wlk (:,:,:) = xoldw (1:ndim, 1:nelec, 1:nwalk, 1)

  else
   call die (here,'mode='+trim(mode)+' should contain either vmc or dmc.')
  endif

  end subroutine coord_elec_wlk_bld

! ==============================================================================
  subroutine dist_e_wlk_bld
! ------------------------------------------------------------------------------
! Description   : distance |ri - 0| between one electron and origin
!
! Created       : J. Toulouse, 07 Nov 2006
! Modified      : J. Toulouse, 02 Oct 2008: add walkers for DMC
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer                       :: walk_i, elec_i, dim_i

! begin

! header
  if (header_exe) then

   call object_create ('dist_e_wlk')
   call object_create ('dist_e_min')
   call object_create ('dist_e_max')

   call object_needed ('nwalk')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('coord_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('dist_e_wlk', dist_e_wlk, nelec, nwalk)
  call object_associate ('dist_e_min', dist_e_min)
  call object_associate ('dist_e_max', dist_e_max)

  dist_e_wlk (:,:) = 0.d0

  do walk_i = 1, nwalk
   do elec_i = 1, nelec

      do dim_i = 1, ndim
        dist_e_wlk (elec_i, walk_i) = dist_e_wlk (elec_i, walk_i) + coord_elec_wlk (dim_i, elec_i, walk_i)**2
      enddo ! dim_i

      dist_e_wlk (elec_i, walk_i) = dsqrt (dist_e_wlk (elec_i, walk_i))

!     minimal e distance in the simulation
      if (dist_e_wlk (elec_i, walk_i) < dist_e_min) then
        dist_e_min = dist_e_wlk (elec_i, walk_i)
      endif

!     maximal e distance in the simulation
      if (dist_e_wlk (elec_i, walk_i) > dist_e_max) then
        dist_e_max = dist_e_wlk (elec_i, walk_i)
      endif

   enddo ! elec_i
  enddo ! walk_i

  end subroutine dist_e_wlk_bld

! ==============================================================================
  subroutine vec_ee_xyz_wlk_bld
! ------------------------------------------------------------------------------
! Description   : cartesian vector (ri - rj) for two electrons
!
! Created       : J. Toulouse, 06 Mar 2006
! Modified      : J. Toulouse, 16 Nov 2006, walkers
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer                       :: elec_i, elec_j, dim_i, walk_i

! begin

! header
  if (header_exe) then

   call object_create ('vec_ee_xyz_wlk')

   call object_needed ('nwalk')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('coord_elec_wlk')

   return

  endif

! allocations
  call object_alloc ('vec_ee_xyz_wlk', vec_ee_xyz_wlk, ndim, nelec, nelec, nwalk)

  do walk_i = 1, nwalk
  do elec_j = 1, nelec
    do elec_i = elec_j+1, nelec
      do dim_i = 1, ndim
        vec_ee_xyz_wlk (dim_i, elec_i, elec_j, walk_i) = coord_elec_wlk (dim_i, elec_i, walk_i) - coord_elec_wlk (dim_i, elec_j, walk_i)
        vec_ee_xyz_wlk (dim_i, elec_j, elec_i, walk_i) = - vec_ee_xyz_wlk (dim_i, elec_i, elec_j, walk_i)
      enddo ! dim_i
    enddo !elec_j
  enddo !elec_i
  enddo !walk_i

  end subroutine vec_ee_xyz_wlk_bld

! ==============================================================================
  subroutine dist_ee_wlk_bld
! ------------------------------------------------------------------------------
! Description   : distance |ri - rj| between two electrons
!
! Created       : J. Toulouse, 06 Mar 2006
! Modified      : J. Toulouse, 16 Mar 2006, walkers
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer                       :: elec_i, elec_j, dim_i, walk_i

! begin

! header
  if (header_exe) then

   call object_create ('dist_ee_wlk')

   call object_needed ('nwalk')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('vec_ee_xyz_wlk')

   return

  endif

! allocations
  call object_alloc ('dist_ee_wlk', dist_ee_wlk, nelec, nelec, nwalk)
  call object_associate ('dist_ee_min', dist_ee_min)
  call object_associate ('dist_ee_max', dist_ee_max)

  dist_ee_wlk (:,:,:) = 0.d0

  do walk_i = 1, nwalk

  do elec_j = 1, nelec
    do elec_i = elec_j+1, nelec

      do dim_i = 1, ndim
        dist_ee_wlk (elec_i, elec_j, walk_i) = dist_ee_wlk (elec_i, elec_j, walk_i) + vec_ee_xyz_wlk (dim_i, elec_i, elec_j, walk_i)**2
      enddo ! dim_i

      dist_ee_wlk (elec_i, elec_j, walk_i) = dsqrt (dist_ee_wlk (elec_i, elec_j,walk_i))

      dist_ee_wlk (elec_j, elec_i, walk_i) = dist_ee_wlk (elec_i, elec_j, walk_i)

!     minimal e-e distance in the simulation
      if (dist_ee_wlk (elec_i, elec_j, walk_i) < dist_ee_min) then
        dist_ee_min = dist_ee_wlk (elec_i, elec_j, walk_i)
      endif

!     maximal e-e distance in the simulation
      if (dist_ee_wlk (elec_i, elec_j, walk_i) > dist_ee_max) then
        dist_ee_max = dist_ee_wlk (elec_i, elec_j, walk_i)
      endif

    enddo !elec_j
  enddo !elec_i

  enddo !walk_i

  call object_modified_by_index (dist_ee_min_index)
  call object_modified_by_index (dist_ee_max_index)

  end subroutine dist_ee_wlk_bld

! ==============================================================================
  subroutine vec_en_xyz_wlk_bld
! ------------------------------------------------------------------------------
! Description   : cartesian vector between electrons and nuclei
!
! Created       : J. Toulouse, 30 Mar 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer                       :: elec_i, cent_i, dim_i, walk_i

! begin

! header
  if (header_exe) then

   call object_create ('vec_en_xyz_wlk')

   call object_needed ('nwalk')
   call object_needed ('ncent')
   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('coord_elec_wlk')
   call object_needed ('cent')

   return

  endif

! allocations
  call object_alloc ('vec_en_xyz_wlk', vec_en_xyz_wlk, ndim, nelec, ncent, nwalk)

  do walk_i = 1, nwalk
   do cent_i = 1, ncent
    do elec_i = 1, nelec
      do dim_i = 1, ndim
        vec_en_xyz_wlk (dim_i, elec_i, cent_i, walk_i) = coord_elec_wlk (dim_i, elec_i, walk_i) - cent (dim_i, cent_i)
      enddo ! dim_i
    enddo ! elec_i
   enddo ! cent_i
  enddo ! walk_i

  end subroutine vec_en_xyz_wlk_bld

! ==============================================================================
  subroutine dist_en_wlk_bld
! ------------------------------------------------------------------------------
! Description   : electron-nuclei distances
!
! Created       : J. Toulouse, 30 Mar 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer  :: elec_i, cent_i, dim_i, walk_i

! header
  if (header_exe) then

   call object_create ('dist_en_wlk')

   call object_needed ('nwalk')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('ncent')
   call object_needed ('vec_en_xyz_wlk')

   return

  endif

! begin

! allocations
  call object_alloc ('dist_en_wlk', dist_en_wlk, nelec, ncent, nwalk)

  dist_en_wlk (:,:,:) = 0.d0

  do walk_i = 1, nwalk
    do cent_i = 1, ncent
      do elec_i = 1, nelec
        do dim_i = 1, ndim
          dist_en_wlk (elec_i, cent_i, walk_i) = dist_en_wlk (elec_i, cent_i, walk_i) + vec_en_xyz_wlk (dim_i, elec_i, cent_i, walk_i)**2
        enddo ! dim_i
        dist_en_wlk (elec_i, cent_i, walk_i) = dsqrt (dist_en_wlk (elec_i, cent_i, walk_i))
       enddo ! elec_i
    enddo ! cent_i
  enddo ! walk_i

  end subroutine dist_en_wlk_bld

! ==============================================================================
  subroutine grd_dist_en_bld
! ------------------------------------------------------------------------------
! Description   : gradients of electon-nuclei distances
!
! Created       : J. Toulouse, 27 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer                       :: cent_i, elec_i, dim_i

! header
  if (header_exe) then

   call object_create ('grd_dist_en')

   call object_needed ('ncent')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('rvec_en')
   call object_needed ('r_en')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_dist_en', grd_dist_en, ndim, nelec, ncent)

  do cent_i = 1, ncent
   do elec_i = 1, nelec
      do dim_i = 1, ndim
        grd_dist_en (dim_i, elec_i, cent_i) = rvec_en (dim_i, elec_i, cent_i) / r_en (elec_i, cent_i)
      enddo ! dim_i
   enddo ! elec_i
  enddo ! cent_i

  end subroutine grd_dist_en_bld

! ==============================================================================
  subroutine lap_dist_en_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of electon-nuclei distances
!
! Created       : J. Toulouse, 27 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer                       :: cent_i, elec_i, dim_i

! header
  if (header_exe) then

   call object_create ('lap_dist_en')

   call object_needed ('ncent')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('r_en')
   call object_needed ('grd_dist_en')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_dist_en', lap_dist_en, nelec, ncent)

  do cent_i = 1, ncent
   do elec_i = 1, nelec
      lap_dist_en (elec_i, cent_i) = 0.d0
      do dim_i = 1, ndim
        lap_dist_en (elec_i, cent_i) = lap_dist_en (elec_i, cent_i) + (1.d0 - grd_dist_en (dim_i, elec_i, cent_i)**2) / r_en (elec_i, cent_i)
      enddo ! dim_i
   enddo ! elec_i
  enddo ! cent_i

!  write(6,*) trim(here),': lap_dist_en=',lap_dist_en

  end subroutine lap_dist_en_bld

! ==============================================================================
  subroutine elec_nb_close_atom_wlk_bld
! ------------------------------------------------------------------------------
! Description   : number of electons closest to each atom
!
! Created       : J. Toulouse, 29 Mar 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer walk_i, elec_i, cent_i, cent_closest_i

! header
  if (header_exe) then

   call object_create ('elec_nb_close_atom_wlk')
   call object_create ('elec_up_nb_close_atom_wlk')
   call object_create ('elec_dn_nb_close_atom_wlk')
   call object_average_walk_define ('elec_nb_close_atom_wlk', 'elec_nb_close_atom_av')
   call object_average_walk_define ('elec_up_nb_close_atom_wlk', 'elec_up_nb_close_atom_av')
   call object_average_walk_define ('elec_dn_nb_close_atom_wlk', 'elec_dn_nb_close_atom_av')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ncent')
   call object_needed ('dist_en_wlk')

   return

  endif

! begin

! allocations
  call object_alloc ('elec_nb_close_atom_wlk', elec_nb_close_atom_wlk, ncent, nwalk)
  call object_alloc ('elec_up_nb_close_atom_wlk', elec_up_nb_close_atom_wlk, ncent, nwalk)
  call object_alloc ('elec_dn_nb_close_atom_wlk', elec_dn_nb_close_atom_wlk, ncent, nwalk)
  call object_alloc ('elec_nb_close_atom_av', elec_nb_close_atom_av, ncent)
  call object_alloc ('elec_up_nb_close_atom_av', elec_up_nb_close_atom_av, ncent)
  call object_alloc ('elec_dn_nb_close_atom_av', elec_dn_nb_close_atom_av, ncent)

  elec_nb_close_atom_wlk (:,:) = 0.d0
  elec_up_nb_close_atom_wlk (:,:) = 0.d0
  elec_dn_nb_close_atom_wlk (:,:) = 0.d0

  do walk_i = 1, nwalk

! spin up
    do elec_i = 1, nup
      cent_closest_i = 1
      do cent_i = 1, ncent
        if (dist_en_wlk (elec_i, cent_i, walk_i) < dist_en_wlk (elec_i, cent_closest_i, walk_i)) then
          cent_closest_i = cent_i
        endif
      enddo ! cent_i
      elec_up_nb_close_atom_wlk (cent_closest_i, walk_i) =  elec_up_nb_close_atom_wlk (cent_closest_i, walk_i) + 1.d0
    enddo ! elec_i

! spin dn
    do elec_i = 1, ndn
      cent_closest_i = 1
      do cent_i = 1, ncent
        if (dist_en_wlk (nup+elec_i, cent_i, walk_i) < dist_en_wlk (nup+elec_i, cent_closest_i, walk_i)) then
          cent_closest_i = cent_i
        endif
      enddo ! cent_i
      elec_dn_nb_close_atom_wlk (cent_closest_i, walk_i) =  elec_dn_nb_close_atom_wlk (cent_closest_i, walk_i) + 1.d0
    enddo ! elec_i

  enddo ! walk_i

  elec_nb_close_atom_wlk (:,:) = elec_up_nb_close_atom_wlk (:,:) + elec_dn_nb_close_atom_wlk (:,:)

!  do walk_i = 1, nwalk
!    do cent_i = 1, ncent
!      write(6,'(2a,i5,a,i3,a,f)') trim(here),': walker #', walk_i,', atom #',cent_i,', number of closest electrons =',elec_nb_close_atom_wlk (cent_i, walk_i)
!    enddo
!  enddo

  end subroutine elec_nb_close_atom_wlk_bld

! ==============================================================================
  subroutine elec_spin_nb_close_atom_wlk_bld
! ------------------------------------------------------------------------------
! Description   : number of electons closest to each atom (including spins)
!
! Created       : J. Toulouse, 09 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer walk_i, elec_i, cent_i, cent_closest_i

! header
  if (header_exe) then

   call object_create ('elec_spin_nb_close_atom_wlk')
   call object_average_walk_define ('elec_spin_nb_close_atom_wlk', 'elec_spin_nb_close_atom_av')

   call object_needed ('nwalk')
   call object_needed ('nelec')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ncent')
   call object_needed ('dist_en_wlk')

   return

  endif

! begin

! allocations
  call object_alloc ('elec_spin_nb_close_atom_wlk', elec_spin_nb_close_atom_wlk, spin_nb, ncent, nwalk)
  call object_alloc ('elec_spin_nb_close_atom_av', elec_spin_nb_close_atom_av, spin_nb, ncent)

  elec_spin_nb_close_atom_wlk (:,:,:) = 0.d0

  do walk_i = 1, nwalk

! spin up
    do elec_i = 1, nup
      cent_closest_i = 1
      do cent_i = 1, ncent
        if (dist_en_wlk (elec_i, cent_i, walk_i) < dist_en_wlk (elec_i, cent_closest_i, walk_i)) then
          cent_closest_i = cent_i
        endif
      enddo ! cent_i
      elec_spin_nb_close_atom_wlk (1, cent_closest_i, walk_i) =  elec_spin_nb_close_atom_wlk (1, cent_closest_i, walk_i) + 1.d0
    enddo ! elec_i

! spin dn
    do elec_i = 1, ndn
      cent_closest_i = 1
      do cent_i = 1, ncent
        if (dist_en_wlk (nup+elec_i, cent_i, walk_i) < dist_en_wlk (nup+elec_i, cent_closest_i, walk_i)) then
          cent_closest_i = cent_i
        endif
      enddo ! cent_i
      elec_spin_nb_close_atom_wlk (2, cent_closest_i, walk_i) =  elec_spin_nb_close_atom_wlk (2, cent_closest_i, walk_i) + 1.d0
    enddo ! elec_i

  enddo ! walk_i

  end subroutine elec_spin_nb_close_atom_wlk_bld

! ==============================================================================
  subroutine frac_elec_nb_close_atom_input_wlk_bld
! ------------------------------------------------------------------------------
! Description   : fraction of configurations with the number of electons closest to each atom
! Description   : identical to the one given in the input, i.e. elec_nb_close_atom_input
!
! Created       : J. Toulouse, 02 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer walk_i

! header
  if (header_exe) then

   call object_create ('frac_elec_nb_close_atom_input_wlk')
   call object_average_walk_define ('frac_elec_nb_close_atom_input_wlk', 'frac_elec_nb_close_atom_input_av')

   call object_needed ('nwalk')
   call object_needed ('elec_nb_close_atom_wlk')
   call object_needed ('elec_nb_close_atom_input')

   return

  endif

! begin

! allocations
  call object_alloc ('frac_elec_nb_close_atom_input_wlk', frac_elec_nb_close_atom_input_wlk, nwalk)
  call object_associate ('frac_elec_nb_close_atom_input_av', frac_elec_nb_close_atom_input_av)

  do walk_i = 1, nwalk
    if (arrays_equal(elec_nb_close_atom_wlk (:, walk_i), elec_nb_close_atom_input (:))) then
     frac_elec_nb_close_atom_input_wlk (walk_i) = 1.d0
    else
     frac_elec_nb_close_atom_input_wlk (walk_i) = 0.d0
    endif
  enddo ! walk_i

  end subroutine frac_elec_nb_close_atom_input_wlk_bld

! ==============================================================================
  subroutine elec_nb_close_atom_encount_bld
! ------------------------------------------------------------------------------
! Description   : configurations with the number of electons closest to each atom
! Description   : encountered in the calculation
!
! Created       : J. Toulouse, 06 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer walk_i, elec_nb_close_atom_encount_i
  logical found

! header
  if (header_exe) then

   call object_create ('elec_nb_close_atom_encount_nb')
   call object_create ('elec_nb_close_atom_encount')
   call object_create ('frac_elec_nb_close_atom_encount_wlk')
   call object_average_walk_define ('frac_elec_nb_close_atom_encount_wlk', 'frac_elec_nb_close_atom_encount_av')

   call object_needed ('ncent')
   call object_needed ('nwalk')
   call object_needed ('elec_nb_close_atom_wlk')

   return

  endif

! begin
# if defined (MPI)
    call die (here, 'routine not implemented for MPI.')
# endif

! allocations
  call object_alloc ('elec_nb_close_atom_encount', elec_nb_close_atom_encount, ncent, elec_nb_close_atom_encount_nb)
  call object_alloc ('frac_elec_nb_close_atom_encount_wlk', frac_elec_nb_close_atom_encount_wlk, elec_nb_close_atom_encount_nb, nwalk)
  call object_alloc ('frac_elec_nb_close_atom_encount_av', frac_elec_nb_close_atom_encount_av, elec_nb_close_atom_encount_nb)

  frac_elec_nb_close_atom_encount_wlk (:,:) = 0.d0

  do walk_i = 1, nwalk
    found = .false.
    do elec_nb_close_atom_encount_i = 1, elec_nb_close_atom_encount_nb
      if (arrays_equal(elec_nb_close_atom_wlk (:, walk_i), elec_nb_close_atom_encount (:, elec_nb_close_atom_encount_i))) then
        frac_elec_nb_close_atom_encount_wlk (elec_nb_close_atom_encount_i, walk_i) = 1.d0
        found = .true.
        exit
      endif
    enddo ! elec_nb_close_atom_encount_i

!   add new configuration to the list of encountered configurations
    if (.not. found) then
      elec_nb_close_atom_encount_nb = elec_nb_close_atom_encount_nb + 1
      call object_alloc ('elec_nb_close_atom_encount', elec_nb_close_atom_encount, ncent, elec_nb_close_atom_encount_nb)
      call object_alloc ('frac_elec_nb_close_atom_encount_wlk', frac_elec_nb_close_atom_encount_wlk, elec_nb_close_atom_encount_nb, nwalk)
      call object_alloc ('frac_elec_nb_close_atom_encount_av', frac_elec_nb_close_atom_encount_av, elec_nb_close_atom_encount_nb)
      elec_nb_close_atom_encount (:, elec_nb_close_atom_encount_nb) = elec_nb_close_atom_wlk (:, walk_i)
      frac_elec_nb_close_atom_encount_wlk (elec_nb_close_atom_encount_nb, walk_i) = 1.d0
    endif
  enddo ! walk_i

  end subroutine elec_nb_close_atom_encount_bld

! ==============================================================================
  subroutine elec_spin_nb_close_atom_encount_bld
! ------------------------------------------------------------------------------
! Description   : configurations with the number of electons closest to each atom
! Description   : encountered in the calculation (including spins)
!
! Created       : J. Toulouse, 09 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer walk_i, elec_spin_nb_close_atom_encount_i
  logical found

! header
  if (header_exe) then

   call object_create ('elec_spin_nb_close_atom_encount_nb')
   call object_create ('elec_spin_nb_close_atom_encount')
   call object_create ('frac_elec_spin_nb_close_atom_encount_wlk')
   call object_average_walk_define ('frac_elec_spin_nb_close_atom_encount_wlk', 'frac_elec_spin_nb_close_atom_encount_av')

   call object_needed ('ncent')
   call object_needed ('nwalk')
   call object_needed ('elec_spin_nb_close_atom_wlk')

   return

  endif

! begin
# if defined (MPI)
    call die (here, 'routine not implemented for MPI.')
# endif

! allocations
  call object_alloc ('elec_spin_nb_close_atom_encount', elec_spin_nb_close_atom_encount, spin_nb, ncent, elec_spin_nb_close_atom_encount_nb)
  call object_alloc ('frac_elec_spin_nb_close_atom_encount_wlk', frac_elec_spin_nb_close_atom_encount_wlk, elec_spin_nb_close_atom_encount_nb, nwalk)
  call object_alloc ('frac_elec_spin_nb_close_atom_encount_av', frac_elec_spin_nb_close_atom_encount_av, elec_spin_nb_close_atom_encount_nb)

  frac_elec_spin_nb_close_atom_encount_wlk (:,:) = 0.d0

  do walk_i = 1, nwalk
    found = .false.
    do elec_spin_nb_close_atom_encount_i = 1, elec_spin_nb_close_atom_encount_nb
      if (arrays_equal(elec_spin_nb_close_atom_wlk (:,:, walk_i), elec_spin_nb_close_atom_encount (:,:, elec_spin_nb_close_atom_encount_i))) then
        frac_elec_spin_nb_close_atom_encount_wlk (elec_spin_nb_close_atom_encount_i, walk_i) = 1.d0
        found = .true.
        exit
      endif
    enddo ! elec_spin_nb_close_atom_encount_i

!   add new configuration to the list of encountered configurations
    if (.not. found) then
      elec_spin_nb_close_atom_encount_nb = elec_spin_nb_close_atom_encount_nb + 1
      call object_alloc ('elec_spin_nb_close_atom_encount', elec_spin_nb_close_atom_encount, spin_nb, ncent, elec_spin_nb_close_atom_encount_nb)
      call object_alloc ('frac_elec_spin_nb_close_atom_encount_wlk', frac_elec_spin_nb_close_atom_encount_wlk, elec_spin_nb_close_atom_encount_nb, nwalk)
      call object_alloc ('frac_elec_spin_nb_close_atom_encount_av', frac_elec_spin_nb_close_atom_encount_av, elec_spin_nb_close_atom_encount_nb)
      elec_spin_nb_close_atom_encount (:,:, elec_spin_nb_close_atom_encount_nb) = elec_spin_nb_close_atom_wlk (:,:, walk_i)
      frac_elec_spin_nb_close_atom_encount_wlk (elec_spin_nb_close_atom_encount_nb, walk_i) = 1.d0
    endif
  enddo ! walk_i

  end subroutine elec_spin_nb_close_atom_encount_bld

! ==============================================================================
  subroutine frac_elec_spin_nb_close_atom_input_wlk_bld
! ------------------------------------------------------------------------------
! Description   : fraction of configurations with the number of electons closest to each atom
! Description   : identical to the one given in the input, i.e. elec_up_nb_close_atom_input
! Description   : and elec_dn_nb_close_atom_input (including spin)
!
! Created       : J. Toulouse, 04 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer walk_i

! header
  if (header_exe) then

   call object_create ('frac_elec_spin_nb_close_atom_input_wlk')
   call object_average_walk_define ('frac_elec_spin_nb_close_atom_input_wlk', 'frac_elec_spin_nb_close_atom_input_av')

   call object_needed ('nwalk')
   call object_needed ('elec_up_nb_close_atom_wlk')
   call object_needed ('elec_dn_nb_close_atom_wlk')
   call object_needed ('elec_up_nb_close_atom_input')
   call object_needed ('elec_dn_nb_close_atom_input')

   return

  endif

! begin

! allocations
  call object_alloc ('frac_elec_spin_nb_close_atom_input_wlk', frac_elec_spin_nb_close_atom_input_wlk, nwalk)
  call object_associate ('frac_elec_spin_nb_close_atom_input_av', frac_elec_spin_nb_close_atom_input_av)

  do walk_i = 1, nwalk
    if (arrays_equal(elec_up_nb_close_atom_wlk (:, walk_i), elec_up_nb_close_atom_input (:)) .and. &
        arrays_equal(elec_dn_nb_close_atom_wlk (:, walk_i), elec_dn_nb_close_atom_input (:)) ) then
     frac_elec_spin_nb_close_atom_input_wlk (walk_i) = 1.d0
    else
     frac_elec_spin_nb_close_atom_input_wlk (walk_i) = 0.d0
    endif
  enddo ! walk_i

  end subroutine frac_elec_spin_nb_close_atom_input_wlk_bld

end module electrons_mod
