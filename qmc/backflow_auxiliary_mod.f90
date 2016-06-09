!see comprehensive notes of F. Petruzielo
module backflow_basis_mod
  use constants_mod    !give access to various common blocks
  use basic_tools_mod, only: die
  use objects_mod, only: object_create, object_needed, object_alloc, object_provide_in_node, object_modified
  use nodes_mod, only: header_exe
  !begin objects
  use electrons_mod, only: dist_en_wlk, vec_en_xyz_wlk, coord_elec_wlk
  use basis_mod, only: basis_fns_cent

  !Declaration of global variables
  real(dp), allocatable                          :: orbitals_bf(:, :)
  real(dp), allocatable                          :: orbital_coeff_bf(:, :)
  real(dp), allocatable                          :: d_orbitals_d_x_beta_bf(:, :, :)

  real(dp), allocatable                          :: normalized_bas_funcs_bf(:, :)
  real(dp), allocatable                          :: bas_funcs_bf(:, :)
  real(dp), allocatable                          :: bas_funcs_normalization_bf(:)
  real(dp), allocatable                          :: d_normalized_bas_funcs_d_x_beta_bf(:, :, :)
  real(dp), allocatable                          :: d_bas_funcs_d_x_beta_bf(:, :, :)
  real(dp), allocatable                          :: d2_bas_funcs_d_x_beta_d_x_gamma_bf(:, :, :, :)
  real(dp), allocatable                          :: d2_bas_funcs_d_x_beta_d_x_gamma_bf_num(:, :, :, :)

  real(dp), allocatable                          :: slat_bas_funcs_bf(:, :)
  real(dp), allocatable                          :: slat_bas_funcs_normalization_bf(:)
  real(dp), allocatable                          :: slat_bas_funcs_radial_bf(:, :)
  real(dp), allocatable                          :: d_slat_bas_funcs_d_x_beta_bf(:, :, :)
  real(dp), allocatable                          :: d_slat_bas_funcs_radial_d_x_beta_bf(:, :, :)
  real(dp), allocatable                          :: d_slat_bas_funcs_radial_d_x_bf(:, :)
  real(dp), allocatable                          :: d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf(:, :, :, :)
  real(dp), allocatable                          :: d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf(:, :, :, :)
  real(dp), allocatable                          :: d2_slat_bas_funcs_radial_d_x2_bf(:, :)

  real(dp), allocatable                          :: bas_funcs_3d_theta_bf(:, :)
  real(dp), allocatable                          :: bas_funcs_3d_phi_bf(:, :)
  real(dp), allocatable                          :: d_bas_funcs_3d_theta_d_x_beta_bf(:, :, :)
  real(dp), allocatable                          :: d_bas_funcs_3d_phi_d_x_beta_bf(:, :, :)
  real(dp), allocatable                          :: d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(:, :, :, :)
  real(dp), allocatable                          :: d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(:, :, :, :)

  !testing
  real(dp), allocatable                          :: orb_cyrus(:, :)
  real(dp), allocatable                          :: dorb_cyrus(:, :, :)
  real(dp), allocatable                          :: phin_cyrus(:, :)
  real(dp), allocatable                          :: dphin_cyrus(:, :, :)
  real(dp), allocatable                          :: d2phin_cyrus(:, :, :, :)

  !end objects

  !interfaces
  interface factorial
     module procedure factorial_array
     module procedure factorial_scalar
  end interface

contains

  ! ==============================================================================
  subroutine orbitals_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : Calculate orbitals evaluated at backflow transformed positions
    !                 by combining normalized basis function with lcao coefficients
    !                 This is a rank two array
    !                 The first index labels the electron
    !                 The second index labels the orbital
    !
    ! Created       : F. Petruzielo 13 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: orb_i, elec_i

    ! header
    if (header_exe) then

       call object_create('orbitals_bf')
       !testing
       call object_create('orb_cyrus')

       call object_needed('norb')
       call object_needed('orbital_coeff_bf')
       call object_needed('normalized_bas_funcs_bf')
       call object_needed('nelec')
       !testing
       call object_needed('orb')

       return

    endif

    ! allocation
    call object_alloc('orbitals_bf', orbitals_bf, nelec, norb)
    call object_alloc('orb_cyrus', orb_cyrus, nelec, norb)
    orbitals_bf(:, :) = 0.d0
    orb_cyrus(:, :) = 0.d0

    do elec_i = 1, nelec
       do orb_i = 1, norb
          orbitals_bf(elec_i, orb_i) = sum(normalized_bas_funcs_bf(:, elec_i) * orbital_coeff_bf(:, orb_i), 1)
       end do
    end do
    orb_cyrus(1:nelec, 1:norb) = orb(1:nelec, 1:norb)

  end subroutine orbitals_bf_bld

  ! ==============================================================================

  ! ==============================================================================
  subroutine d_orbitals_d_x_beta_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : Calculate derivative of orbitals as a function of backflow transformed coordinates
    !                 with respect to backflow transformed coordinates (so just like non backflow case, but r is replaced by x)
    !                 by combining normalized basis function with lcao coefficients
    !                 This is a rank three array
    !                 The first index labels the dimension
    !                 The second index labels the electron
    !                 The third index labels the orbital
    !
    ! Created       : F. Petruzielo 14 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: orb_i, elec_i, dim_i

    ! header
    if (header_exe) then

       call object_create('d_orbitals_d_x_beta_bf')
       !testing
       call object_create('dorb_cyrus')
       call object_needed('nelec')
       call object_needed('norb')
       call object_needed('ndim')
       call object_needed('orbital_coeff_bf')
       call object_needed('d_normalized_bas_funcs_d_x_beta_bf')

       !testing
       call object_needed('dorb')

       return

    endif

    ! allocation
    call object_alloc('d_orbitals_d_x_beta_bf', d_orbitals_d_x_beta_bf, ndim, nelec, norb)
    call object_alloc('dorb_cyrus', dorb_cyrus, ndim, nelec, norb)
    d_orbitals_d_x_beta_bf(:, :, :) = 0.d0
    dorb_cyrus(:, :, :) = 0.d0

    do elec_i = 1, nelec
       do orb_i = 1, norb
          do dim_i = 1, ndim
             d_orbitals_d_x_beta_bf(dim_i, elec_i, orb_i) = sum(d_normalized_bas_funcs_d_x_beta_bf(dim_i, :, elec_i) * orbital_coeff_bf(:, orb_i), 1)
          end do
       end do
    end do
    dorb_cyrus(1:ndim, 1:nelec, 1:norb) = dorb(1:ndim, 1:nelec, 1:norb)

  end subroutine d_orbitals_d_x_beta_bf_bld

  ! ==============================================================================

  ! ==============================================================================
  subroutine orbital_coeff_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : lcao coefficients of each normalized basis function for each orbital of
    !                 wavefunction specified by iwf
    !                 This is a rank two array
    !                 The first index labels the basis function
    !                 The second index labels the orbital
    !
    ! Created       : F. Petruzielo 13 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: orb_i

    ! header
    if (header_exe) then

       call object_create('orbital_coeff_bf')

       call object_needed('nbasis')
       call object_needed('norb')
       call object_needed('bas_funcs_normalization_bf')
       call object_needed('iwf')
       call object_needed('coef')

       return

    endif

    ! allocation
    call object_alloc('orbital_coeff_bf', orbital_coeff_bf, nbasis, norb)
    orbital_coeff_bf(:, :) = 0.d0

    !testing
    if (iwf .ne. 1) then
       stop "Are you sure of what's going on?"
    end if

    do orb_i = 1, norb
       orbital_coeff_bf(1:nbasis, orb_i) = coef(1:nbasis, orb_i, iwf) / bas_funcs_normalization_bf(1:nbasis)
    end do

  end subroutine orbital_coeff_bf_bld

  ! ==============================================================================

  ! ==============================================================================
  subroutine normalized_bas_funcs_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : Normalized basis functions evaluated at backflow transformed coordinates.
    !                 This is a rank two array
    !                 The first index labels the basis function
    !                 The second index labels the electron
    !                 Note that when one electron moves, all of the other electrons have new backflow transformed coordinates
    !                 so that all basis functions need to be re-evaluated to evaluate the orbitals and hence the wavefunction.
    !                 However, if only doing electron-nuclear backflow then can still do everything in terms of one-electron moves (use SM formula, etc.)
    !                 Should discuss this with Cyrus and Julien.
    !
    ! Created       : F. Petruzielo 10 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i

    ! header
    if (header_exe) then

       call object_create('normalized_bas_funcs_bf')

       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('bas_funcs_bf')
       call object_needed('bas_funcs_normalization_bf')

       return

    endif

    ! allocation
    call object_alloc('normalized_bas_funcs_bf', normalized_bas_funcs_bf, nbasis, nelec)
    normalized_bas_funcs_bf(:, :) = 0.d0

    do bas_i = 1, nbasis
       normalized_bas_funcs_bf(bas_i, :) = bas_funcs_normalization_bf(bas_i) * bas_funcs_bf(bas_i, :)
    enddo

  end subroutine normalized_bas_funcs_bf_bld

  ! ==============================================================================

  ! ==============================================================================
  subroutine d_normalized_bas_funcs_d_x_beta_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : derivative of normalized basis functions as a function of backflow transformed coordinates
    !                 with respect to backflow transformed coordinates. (so just like non backflow case, but r is replaced by x)
    !                 This is a rank three array
    !                 The first index labels the dimension (beta)
    !                 The second index labels the basis function
    !                 The third index labeles the electron
    !
    ! Created       : F. Petruzielo 14 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i

    ! header
    if (header_exe) then

       call object_create('d_normalized_bas_funcs_d_x_beta_bf')

       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('ndim')
       call object_needed('d_bas_funcs_d_x_beta_bf')
       call object_needed('bas_funcs_normalization_bf')

       return

    endif

    ! allocation
    call object_alloc('d_normalized_bas_funcs_d_x_beta_bf', d_normalized_bas_funcs_d_x_beta_bf, ndim, nbasis, nelec)
    d_normalized_bas_funcs_d_x_beta_bf(:, :, :) = 0.d0

    do bas_i = 1, nbasis
       d_normalized_bas_funcs_d_x_beta_bf(:, bas_i, :) = bas_funcs_normalization_bf(bas_i) * d_bas_funcs_d_x_beta_bf(:, bas_i, :)
    enddo

  end subroutine d_normalized_bas_funcs_d_x_beta_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine bas_funcs_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : Un-normalized basis functions evaluated at backflow transformed coordinates of each electron.
    !                 This is a rank two array
    !                 The first index labels the basis function
    !                 The second index labels the electron
    !
    ! Created       : F. Petruzielo 10 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    character(len=max_string_len_rout), save :: lhere = 'bas_funcs_bf_bld'

    ! header
    if (header_exe) then

       call object_create('bas_funcs_bf')
!testing
       call object_create('phin_cyrus')

       call object_needed('nbasis')
       call object_needed('ndim')
       call object_needed('numr')
       call object_needed('n_bas')
       call object_needed('nelec')
       !testing
       call object_needed('phin')

       return

    endif

    ! allocation
    call object_alloc('bas_funcs_bf', bas_funcs_bf, nbasis, nelec)
    bas_funcs_bf(:, :) = 0.d0

    call object_alloc('phin_cyrus', phin_cyrus, nbasis, nelec)
    phin_cyrus(:, :) = 0.d0


    phin_cyrus(1:nbasis, 1:nelec) = phin(1:nbasis, 1:nelec)

    if (numr .le. 0 .and. n_bas(1) > 0 .and. ndim==3) then
       !slater basis in 3d
       call object_provide_in_node(lhere, 'slat_bas_funcs_bf')
       bas_funcs_bf(:, :) = slat_bas_funcs_bf(:, :)
    else
       call die(lhere, "Only Slater basis in 3d implemented for backflow")
    end if

  end subroutine bas_funcs_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d_bas_funcs_d_x_beta_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : derivative of un-normalized basis functions as a function of backflow transformed coordinates
    !                 with respect to backflow transformed coordinates. (so just like non backflow case, but r is replaced by x)
    !                 This is a rank three array
    !                 The first index labels the dimension (beta)
    !                 The second index labels the basis function
    !                 The third index labeles the electron
    !
    ! Created       : F. Petruzielo 14 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    character(len=max_string_len_rout), save :: lhere = 'd_bas_funcs_d_x_beta_bf_bld'

    ! header
    if (header_exe) then

       call object_create('d_bas_funcs_d_x_beta_bf')
!testing
       call object_create('dphin_cyrus')

       call object_needed('nbasis')
       call object_needed('ndim')
       call object_needed('numr')
       call object_needed('n_bas')
       call object_needed('nelec')
       !testing
       call object_needed('dphin')

       return

    endif

    ! allocation
    call object_alloc('d_bas_funcs_d_x_beta_bf', d_bas_funcs_d_x_beta_bf, ndim, nbasis, nelec)
    d_bas_funcs_d_x_beta_bf(:, :, :) = 0.d0

    call object_alloc('dphin_cyrus', dphin_cyrus, ndim, nbasis, nelec)
    dphin_cyrus(:, :, :) = 0.d0

    dphin_cyrus(1:ndim, 1:nbasis, 1:nelec) = dphin(1:ndim, 1:nbasis, 1:nelec)

    if (numr .le. 0 .and. n_bas(1) > 0 .and. ndim==3) then
       !slater basis in 3d
       call object_provide_in_node(lhere, 'd_slat_bas_funcs_d_x_beta_bf')
       d_bas_funcs_d_x_beta_bf(:, :, :) = d_slat_bas_funcs_d_x_beta_bf(:, :, :)
    else
       call die(lhere, "Only Slater basis in 3d implemented for backflow")
    end if

  end subroutine d_bas_funcs_d_x_beta_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d2_bas_funcs_d_x_beta_d_x_gamma_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : second derivative of un-normalized basis functions as a function of backflow transformed coordinates
    !                 with respect to backflow transformed coordinates. (so just like non backflow case, but r is replaced by x)
    !                 This is a rank four array
    !                 The first index labels the dimension (gamma)
    !                 The second index labels the dimension (beta)
    !                 The third index labels the basis function
    !                 The fourth index labeles the electron
    !
    ! Created       : F. Petruzielo 15 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    character(len=max_string_len_rout), save :: lhere = 'd2_bas_funcs_d_x_beta_d_x_gamma_bf_bld'

    ! header
    if (header_exe) then

       call object_create('d2_bas_funcs_d_x_beta_d_x_gamma_bf')

       call object_needed('nbasis')
       call object_needed('ndim')
       call object_needed('numr')
       call object_needed('n_bas')
       call object_needed('nelec')

       return

    endif

    ! allocation
    call object_alloc('d2_bas_funcs_d_x_beta_d_x_gamma_bf', d2_bas_funcs_d_x_beta_d_x_gamma_bf, ndim, ndim, nbasis, nelec)
    d2_bas_funcs_d_x_beta_d_x_gamma_bf(:, :, :, :) = 0.d0

    if (numr .le. 0 .and. n_bas(1) > 0 .and. ndim==3) then
       !slater basis in 3d
       call object_provide_in_node(lhere, 'd2_slat_bas_funcs_d_x_beta_d_x_gamma_bf')
       d2_bas_funcs_d_x_beta_d_x_gamma_bf(:, :, :, :) = d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf(:, :, :, :)
    else
       call die(lhere, "Only Slater basis in 3d implemented for backflow")
    end if

  end subroutine d2_bas_funcs_d_x_beta_d_x_gamma_bf_bld

  !==============================================================================

  !==============================================================================================

  subroutine d2_bas_funcs_d_x_beta_d_x_gamma_bf_num_bld
    !---------------------------------------------------------------------------
    ! Description : build d2_bas_funcs_d_x_beta_d_x_gamma_bf_num. This is a rank four array.
    !               Numerical second derivative of bas_funcs_bf :
    !
    !               d^2 bas_funcs_bf
    !               ------------------
    !               d x_beta d x_gamma
    !
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               Note beta/gamma refer to cartesian coords. In 3-d, x_beta/x_gamma takes on values x,y,z.
    !               The first index is the over spatial dimension (gamma)
    !               The second index is the over spatial dimension (beta)
    !               The third index has dimensions equal to the number of basis functions.
    !               The fourth index has dimension equal to the number of electrons.
    !
    ! Created     : F. Petruzielo, 15 Apr 2009
    !---------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !need nwalk, ncent, nelec, ndim

    !local
    integer  :: dim_i, dim_j
    real(dp), allocatable :: d_bas_funcs_d_x_beta_bf_plus(:, :, :), d_bas_funcs_d_x_beta_bf_minus(:, :, :)
    real(dp), allocatable :: coord_elec_wlk_temp(:, :, :)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'd2_bas_funcs_d_x_beta_d_x_gamma_bf_num_bld'

    ! header
    if (header_exe) then

       call object_create('d2_bas_funcs_d_x_beta_d_x_gamma_bf_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('nbasis')
       call object_needed('d_bas_funcs_d_x_beta_bf')
       call object_needed('coord_elec_wlk')

       return
    endif

    call object_alloc('d2_bas_funcs_d_x_beta_d_x_gamma_bf_num', d2_bas_funcs_d_x_beta_d_x_gamma_bf_num, ndim, ndim, nbasis, nelec)
    d2_bas_funcs_d_x_beta_d_x_gamma_bf_num(:, :, :, :) = 0.d0

    allocate(d_bas_funcs_d_x_beta_bf_plus(ndim, nbasis, nelec))
    allocate(d_bas_funcs_d_x_beta_bf_minus(ndim, nbasis, nelec))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

    do dim_i = 1, ndim  !beta
       do dim_j = 1, ndim  !gamma
          coord_elec_wlk = coord_elec_wlk_temp
          !shift coordinates in a particular dimension, x,y,z, etc
          coord_elec_wlk(dim_j, :, 1) = coord_elec_wlk(dim_j, :, 1) + epsilon
          call object_modified('coord_elec_wlk')
          !get updated d_bas_funcs_d_x_beta_bf
          call object_provide_in_node(lhere, 'd_bas_funcs_d_x_beta_bf')
          d_bas_funcs_d_x_beta_bf_plus = d_bas_funcs_d_x_beta_bf
          coord_elec_wlk = coord_elec_wlk_temp
          !shift coordinates in a particular dimension, x,y,z, etc
          coord_elec_wlk(dim_j, :, 1) = coord_elec_wlk(dim_j, :, 1) - epsilon
          call object_modified('coord_elec_wlk')
          !get updated d_bas_funcs_d_x_beta_bf
          call object_provide_in_node(lhere, 'd_bas_funcs_d_x_beta_bf')
          d_bas_funcs_d_x_beta_bf_minus = d_bas_funcs_d_x_beta_bf
          !calculate derivative
          d2_bas_funcs_d_x_beta_d_x_gamma_bf_num(dim_j, dim_i, :, :) = (d_bas_funcs_d_x_beta_bf_plus(dim_i, :, :) - d_bas_funcs_d_x_beta_bf_minus(dim_i, :, :) ) / (2.d0*epsilon)
       end do
    end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'd_bas_funcs_d_x_beta_bf')

  end subroutine d2_bas_funcs_d_x_beta_d_x_gamma_bf_num_bld

  !==============================================================================================


  ! ==============================================================================

  subroutine bas_funcs_normalization_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : normalization of the basis functions
    !                 This is a rank one array
    !                 The first index labels the basis function
    !
    ! Created       : F. Petruzielo 13 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    character(len=max_string_len_rout), save :: lhere = 'bas_funcs_normalization_bf_bld'

    ! header
    if (header_exe) then

       call object_create('bas_funcs_normalization_bf')

       call object_needed('nbasis')
       call object_needed('ndim')
       call object_needed('numr')
       call object_needed('n_bas')

       return

    endif

    ! allocation
    call object_alloc('bas_funcs_normalization_bf', bas_funcs_normalization_bf, nbasis)
    bas_funcs_normalization_bf(:) = 0.d0

    if (numr .le. 0 .and. n_bas(1) > 0 .and. ndim==3) then
       !slater basis in 3d
       call object_provide_in_node(lhere, 'slat_bas_funcs_normalization_bf')
       bas_funcs_normalization_bf(:) = slat_bas_funcs_normalization_bf(:)
    else
       call die(lhere, "Only Slater basis in 3d implemented for backflow")
    end if

  end subroutine bas_funcs_normalization_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine slat_bas_funcs_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : Un-normalized slater basis functions evaluated at backflow transformed coordinates of each electron.
    !                 This is a rank two array
    !                 The first index labels the basis function
    !                 The second index labels the electron
    !
    ! Created       : F. Petruzielo 13 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local


    ! header
    if (header_exe) then

       call object_create('slat_bas_funcs_bf')

       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('slat_bas_funcs_radial_bf')
       call object_needed('bas_funcs_3d_theta_bf')
       call object_needed('bas_funcs_3d_phi_bf')

       return

    endif

    ! allocation
    call object_alloc('slat_bas_funcs_bf', slat_bas_funcs_bf, nbasis, nelec)
    slat_bas_funcs_bf(:, :) = 0.d0

    slat_bas_funcs_bf(:, :) = slat_bas_funcs_radial_bf(:, :) * bas_funcs_3d_theta_bf(:, :) * bas_funcs_3d_phi_bf(:, :)

  end subroutine slat_bas_funcs_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d_slat_bas_funcs_d_x_beta_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : deriviative of Un-normalized slater basis functions as a function of backflow transformed coordinate.
    !                 with respect to backflow transformed coordinate (so just like non backflow case, but r is replaced by x)
    !                 This is a rank three array
    !                 The first index labels the dimension
    !                 The second index labels the basis function
    !                 The third index labels the electron
    !
    ! Created       : F. Petruzielo 14 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: dim_i

    ! header
    if (header_exe) then

       call object_create('d_slat_bas_funcs_d_x_beta_bf')

       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('ndim')
       call object_needed('d_slat_bas_funcs_radial_d_x_beta_bf')
       call object_needed('d_bas_funcs_3d_theta_d_x_beta_bf')
       call object_needed('d_bas_funcs_3d_phi_d_x_beta_bf')
       call object_needed('slat_bas_funcs_radial_bf')
       call object_needed('bas_funcs_3d_theta_bf')
       call object_needed('bas_funcs_3d_phi_bf')

       return

    endif

    ! allocation
    call object_alloc('d_slat_bas_funcs_d_x_beta_bf', d_slat_bas_funcs_d_x_beta_bf, ndim, nbasis, nelec)
    d_slat_bas_funcs_d_x_beta_bf(:, :, :) = 0.d0

    do dim_i = 1, ndim
       d_slat_bas_funcs_d_x_beta_bf(dim_i, :, :) = d_slat_bas_funcs_radial_d_x_beta_bf(dim_i, :, :)  * bas_funcs_3d_phi_bf(:, :) * bas_funcs_3d_theta_bf(:, :) + slat_bas_funcs_radial_bf(:, :) * (d_bas_funcs_3d_phi_d_x_beta_bf(dim_i, :, :) * bas_funcs_3d_theta_bf(:, :) + d_bas_funcs_3d_theta_d_x_beta_bf(dim_i, :, :) *  bas_funcs_3d_phi_bf(:, :) )
    end do

  end subroutine d_slat_bas_funcs_d_x_beta_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : second deriviative of Un-normalized slater basis functions as a function of backflow transformed coordinate.
    !                 with respect to backflow transformed coordinate (so just like non backflow case, but r is replaced by x)
    !                 This is a rank four array
    !                 The first index labels the dimension (gamma)
    !                 The second index labels the dimension (beta)
    !                 The third index labels the basis function
    !                 The fourth index labels the electron
    !
    ! Created       : F. Petruzielo 14 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: dim_i, dim_j

    ! header
    if (header_exe) then

       call object_create('d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf')

       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('ndim')
       call object_needed('slat_bas_funcs_radial_bf')
       call object_needed('d_slat_bas_funcs_radial_d_x_beta_bf')
       call object_needed('d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf')
       call object_needed('bas_funcs_3d_theta_bf')
       call object_needed('bas_funcs_3d_phi_bf')
       call object_needed('d_bas_funcs_3d_theta_d_x_beta_bf')
       call object_needed('d_bas_funcs_3d_phi_d_x_beta_bf')
       call object_needed('d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf')
       call object_needed('d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf')

       return

    endif

    ! allocation
    call object_alloc('d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf', d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf, ndim, ndim, nbasis, nelec)
    d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf(:, :, :, :) = 0.d0

    do dim_i = 1, ndim !beta
       do dim_j = 1, ndim !gamma
          d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf(dim_j, dim_i, :, :) = d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf(dim_j, dim_i, :, :)  * bas_funcs_3d_phi_bf(:, :) * bas_funcs_3d_theta_bf(:, :) + d_slat_bas_funcs_radial_d_x_beta_bf(dim_i, :, :) * (d_bas_funcs_3d_phi_d_x_beta_bf(dim_j, :, :) * bas_funcs_3d_theta_bf(:, :) + d_bas_funcs_3d_theta_d_x_beta_bf(dim_j, :, :) *  bas_funcs_3d_phi_bf(:, :) ) +  d_slat_bas_funcs_radial_d_x_beta_bf(dim_j, :, :) * (d_bas_funcs_3d_phi_d_x_beta_bf(dim_i, :, :) * bas_funcs_3d_theta_bf(:, :) + d_bas_funcs_3d_theta_d_x_beta_bf(dim_i, :, :) *  bas_funcs_3d_phi_bf(:, :) ) + slat_bas_funcs_radial_bf(:, :) * ( d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(dim_j, dim_i, :, :) * bas_funcs_3d_theta_bf(:, :) + d_bas_funcs_3d_theta_d_x_beta_bf(dim_i, :, :) * d_bas_funcs_3d_phi_d_x_beta_bf(dim_j, :, :) + d_bas_funcs_3d_theta_d_x_beta_bf(dim_j, :, :) * d_bas_funcs_3d_phi_d_x_beta_bf(dim_i, :, :) + d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(dim_j, dim_i, :, :) * bas_funcs_3d_phi_bf(:, :) )
       end do
    end do

  end subroutine d2_slat_bas_funcs_d_x_beta_d_x_gamma_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine slat_bas_funcs_normalization_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : normalization of the 3-d slater basis functions
    !                 This is a rank one array
    !                 The first index labels the slater basis function
    !
    ! Created       : F. Petruzielo 13 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    ! header
    if (header_exe) then

       call object_create('slat_bas_funcs_normalization_bf')

       call object_needed('nbasis')
       call object_needed('n_bas')
       call object_needed('l_bas')
       call object_needed('zex')
       call object_needed('iwf')

       return

    endif

    ! allocation
    call object_alloc('slat_bas_funcs_normalization_bf', slat_bas_funcs_normalization_bf, nbasis)
    slat_bas_funcs_normalization_bf(:) = 0.d0

    slat_bas_funcs_normalization_bf(1:nbasis) = sqrt( (2 * l_bas(1:nbasis) + 1) / (4 * pi) * (2 * zex(1:nbasis, iwf))**(2 * n_bas(1:nbasis) + 1) / factorial(2 * n_bas(1:nbasis)))

  end subroutine slat_bas_funcs_normalization_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine slat_bas_funcs_radial_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : radial part of slater basis functions evaluated at backflow transformed coordinates of each electron.
    !                 This is a rank two array
    !                 The first index labels the basis function
    !                 The second index labels the electron
    !testing          Note we can just evaluate for unique n quantum number (implement this)
    !
    ! Created       : F. Petruzielo 13 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i

    ! header
    if (header_exe) then

       call object_create('slat_bas_funcs_radial_bf')

       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('n_bas')
       call object_needed('zex')
       call object_needed('iwf')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('dist_en_wlk')

       return

    endif

    ! allocation
    call object_alloc('slat_bas_funcs_radial_bf', slat_bas_funcs_radial_bf, nbasis, nelec)
    slat_bas_funcs_radial_bf(:, :) = 0.d0

    do bas_i = 1, nbasis
       slat_bas_funcs_radial_bf(bas_i, :) = dist_en_wlk(:, basis_fns_cent(bas_i), 1) ** (n_bas(bas_i) - 1) * exp(-zex(bas_i, iwf) * dist_en_wlk(:, basis_fns_cent(bas_i), 1))
    end do

  end subroutine slat_bas_funcs_radial_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d_slat_bas_funcs_radial_d_x_beta_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : derivative of radial part of slater basis functions as a function of backflow transformed distance 
    !                 with respect to backflow transformed coordinate (so use chain rule) (so just like non backflow case, but r is replaced by x).
    !                 This is a rank three array
    !                 The first index labels the dimension (beta)
    !                 The second index labels the basis function
    !                 The third index labels the electron
    !testing          Note we can just evaluate for unique n quantum number (implement this)
    !
    ! Created       : F. Petruzielo 15 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i, dim_i

    ! header
    if (header_exe) then

       call object_create('d_slat_bas_funcs_radial_d_x_beta_bf')

       call object_needed('ndim')
       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('d_slat_bas_funcs_radial_d_x_bf')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')

       return

    endif

    ! allocation
    call object_alloc('d_slat_bas_funcs_radial_d_x_beta_bf', d_slat_bas_funcs_radial_d_x_beta_bf, ndim, nbasis, nelec)
    d_slat_bas_funcs_radial_d_x_beta_bf(:, :, :) = 0.d0


    do bas_i = 1, nbasis
       do dim_i = 1, ndim
          d_slat_bas_funcs_radial_d_x_beta_bf(dim_i, bas_i, :) = vec_en_xyz_wlk(dim_i, :, basis_fns_cent(bas_i), 1) / dist_en_wlk(:, basis_fns_cent(bas_i), 1) * d_slat_bas_funcs_radial_d_x_bf(bas_i, :)
       end do
    end do

  end subroutine d_slat_bas_funcs_radial_d_x_beta_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : second derivative of radial part of slater basis functions as a function of backflow transformed distance 
    !                 with respect to backflow transformed coordinate (so use chain rule) (so just like non backflow case, but r is replaced by x).
    !                 This is a rank four array
    !                 The first index labels the dimension (gamma)
    !                 The second index labels the dimension (beta)
    !                 The third index labels the basis function
    !                 The fourth index labels the electron
    !testing          Note we can just evaluate for unique n quantum number (implement this)
    !
    ! Created       : F. Petruzielo 15 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i, dim_i, dim_j

    ! header
    if (header_exe) then

       call object_create('d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf')

       call object_needed('ndim')
       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('d2_slat_bas_funcs_radial_d_x2_bf')
       call object_needed('d_slat_bas_funcs_radial_d_x_bf')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')

       return

    endif

    ! allocation
    call object_alloc('d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf', d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf, ndim, ndim, nbasis, nelec)
    d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf(:, :, :, :) = 0.d0

    do bas_i = 1, nbasis
       do dim_i = 1, ndim
          !delta beta gamma
          d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf(dim_i, dim_i, bas_i, :) = d_slat_bas_funcs_radial_d_x_bf(bas_i, :) / dist_en_wlk(:, basis_fns_cent(bas_i), 1)
          do dim_j = 1, ndim
             d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf(dim_j, dim_i, bas_i, :) = d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf(dim_j, dim_i, bas_i, :) + vec_en_xyz_wlk(dim_i, :, basis_fns_cent(bas_i), 1) * vec_en_xyz_wlk(dim_j, :, basis_fns_cent(bas_i), 1) / dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2  * d2_slat_bas_funcs_radial_d_x2_bf(bas_i, :) -  vec_en_xyz_wlk(dim_i, :, basis_fns_cent(bas_i), 1) * vec_en_xyz_wlk(dim_j, :, basis_fns_cent(bas_i), 1) / dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3 * d_slat_bas_funcs_radial_d_x_bf(bas_i, :)
          end do
       end do
    end do

  end subroutine d2_slat_bas_funcs_radial_d_x_beta_d_x_gamma_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d_slat_bas_funcs_radial_d_x_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : derivative of radial part of slater basis functions as a function of backflow transformed distance
    !                 with respect to backflow transformed distance (so just like non backflow case, but r is replaced by x)
    !                 This is a rank two array
    !                 The first index labels the basis function
    !                 The second index labels the electron
    !testing          Note we can just evaluate for unique n quantum number (implement this)
    !
    ! Created       : F. Petruzielo 15 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i

    ! header
    if (header_exe) then

       call object_create('d_slat_bas_funcs_radial_d_x_bf')

       call object_needed('ndim')
       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('zex')
       call object_needed('n_bas')
       call object_needed('iwf')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')

       return

    endif

    ! allocation
    call object_alloc('d_slat_bas_funcs_radial_d_x_bf', d_slat_bas_funcs_radial_d_x_bf, nbasis, nelec)
    d_slat_bas_funcs_radial_d_x_bf(:, :) = 0.d0

    do bas_i = 1, nbasis
       d_slat_bas_funcs_radial_d_x_bf(bas_i, :) =  dist_en_wlk(:, basis_fns_cent(bas_i), 1) ** (n_bas(bas_i) - 2) * (n_bas(bas_i) - 1 - zex(bas_i, iwf) * dist_en_wlk(:, basis_fns_cent(bas_i), 1)) * exp(-zex(bas_i, iwf) * dist_en_wlk(:, basis_fns_cent(bas_i), 1))
    end do

  end subroutine d_slat_bas_funcs_radial_d_x_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d2_slat_bas_funcs_radial_d_x2_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : second derivative of radial part of slater basis functions as a function of backflow transformed distance
    !                 with respect to backflow transformed distance (so just like non backflow case, but r is replaced by x)
    !                 This is a rank two array
    !                 The first index labels the basis function
    !                 The second index labels the electron
    !testing          Note we can just evaluate for unique n quantum number (implement this)
    !
    ! Created       : F. Petruzielo 15 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i

    ! header
    if (header_exe) then

       call object_create('d2_slat_bas_funcs_radial_d_x2_bf')

       call object_needed('ndim')
       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('zex')
       call object_needed('n_bas')
       call object_needed('iwf')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')

       return

    endif

    ! allocation
    call object_alloc('d2_slat_bas_funcs_radial_d_x2_bf', d2_slat_bas_funcs_radial_d_x2_bf, nbasis, nelec)
    d2_slat_bas_funcs_radial_d_x2_bf(:, :) = 0.d0

    do bas_i = 1, nbasis
       d2_slat_bas_funcs_radial_d_x2_bf(bas_i, :) =  dist_en_wlk(:, basis_fns_cent(bas_i), 1) ** (n_bas(bas_i) - 3) * ( (n_bas(bas_i) - 1) * (n_bas(bas_i) - 2 - zex(bas_i, iwf) * dist_en_wlk(:, basis_fns_cent(bas_i), 1)) - (n_bas(bas_i) - 1 - zex(bas_i, iwf) * dist_en_wlk(:, basis_fns_cent(bas_i), 1)) * zex(bas_i, iwf) * dist_en_wlk(:, basis_fns_cent(bas_i), 1)) * exp(-zex(bas_i, iwf) * dist_en_wlk(:, basis_fns_cent(bas_i), 1))
    end do

  end subroutine d2_slat_bas_funcs_radial_d_x2_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine bas_funcs_3d_theta_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : 3d_theta part of analytic basis functions evaluated at backflow transformed coordinates of each electron.
    !                 This is a rank two array
    !                 The first index labels the basis function
    !                 The second index labels the electron
    !
    ! Created       : F. Petruzielo 13 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i
    character(len=max_string_len_rout), save :: lhere = 'bas_funcs_3d_theta_bf_bld'

    ! header
    if (header_exe) then

       call object_create('bas_funcs_3d_theta_bf')

       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('l_bas')
       call object_needed('m_bas')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('dist_en_wlk')

       return

    endif

    ! allocation
    call object_alloc('bas_funcs_3d_theta_bf', bas_funcs_3d_theta_bf, nbasis, nelec)
    bas_funcs_3d_theta_bf(:, :) = 0.d0


    if (maxval(l_bas) > 4) then
       call die(lhere, "Only implemented for l .le. 4")
    end if

    !calculate theta part of basis functions for the l,m values of each basis function
    do bas_i = 1, nbasis
       select case (l_bas(bas_i))

       case (0)
          select case(abs(m_bas(bas_i)))
          case(0)
             bas_funcs_3d_theta_bf(bas_i, :) = 1.d0
          end select

       case(1)
          select case(abs(m_bas(bas_i)))
          case(0)
             bas_funcs_3d_theta_bf(bas_i, :) = vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1) / dist_en_wlk(:, basis_fns_cent(bas_i), 1)
          case(1)
             bas_funcs_3d_theta_bf(bas_i, :) = 1 / dist_en_wlk(:, basis_fns_cent(bas_i), 1)
          end select

       case(2)
          select case(abs(m_bas(bas_i)))
          case(0)
             bas_funcs_3d_theta_bf(bas_i, :) = 1 - (3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2)
          case(1)
             bas_funcs_3d_theta_bf(bas_i, :) = (sqrt(3.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2
          case(2)
             bas_funcs_3d_theta_bf(bas_i, :) = sqrt(3.d0)/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2)
          end select

       case(3)
          select case(abs(m_bas(bas_i)))
          case(0)
             bas_funcs_3d_theta_bf(bas_i, :) = (-3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1) + 2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**3)/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3)
          case(1)
             bas_funcs_3d_theta_bf(bas_i, :) = -(sqrt(1.5d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3)
          case(2)
             bas_funcs_3d_theta_bf(bas_i, :) = (sqrt(15.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3)
          case(3)
             bas_funcs_3d_theta_bf(bas_i, :) = sqrt(2.5d0)/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3)
          end select

       case(4)
          select case(abs(m_bas(bas_i)))
          case(0)
             bas_funcs_3d_theta_bf(bas_i, :) = (3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 24*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 8*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4)/(8.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4)
          case(1)
             bas_funcs_3d_theta_bf(bas_i, :) = (sqrt(2.5d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(-3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) + 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4)
          case(2)
             bas_funcs_3d_theta_bf(bas_i, :) = -(sqrt(5.d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 6*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(4.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4)
          case(3)
             bas_funcs_3d_theta_bf(bas_i, :) = (sqrt(17.5d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4)
          case(4)
             bas_funcs_3d_theta_bf(bas_i, :) = sqrt(35.d0)/(8.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4)
          end select

       end select

    end do

  end subroutine bas_funcs_3d_theta_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d_bas_funcs_3d_theta_d_x_beta_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : derivative of 3d_theta part of analytic basis functions as a function of backflow transformed coordinate.
    !                 with respect to backflow transformed coordinate (so just like non backflow case, but r is replaced by x)
    !                 This is a rank three array
    !                 The first index labels the dimensions (beta)
    !                 The second index labels the basis function
    !                 The third index labels the electron
    !
    ! Created       : F. Petruzielo 14 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i
    character(len=max_string_len_rout), save :: lhere = 'd_bas_funcs_3d_theta_d_x_beta_bf_bld'

    ! header
    if (header_exe) then

       call object_create('d_bas_funcs_3d_theta_d_x_beta_bf')

       call object_needed('ndim')
       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('l_bas')
       call object_needed('m_bas')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('dist_en_wlk')

       return

    endif

    ! allocation
    call object_alloc('d_bas_funcs_3d_theta_d_x_beta_bf', d_bas_funcs_3d_theta_d_x_beta_bf, ndim, nbasis, nelec)
    d_bas_funcs_3d_theta_d_x_beta_bf(:, :, :) = 0.d0

    if (maxval(l_bas) > 4) then
       call die(lhere, "Only implemented for l .le. 4")
    end if

    !calculate derivative of theta part of basis functions for the l,m values of each basis function
    do bas_i = 1, nbasis
       select case (l_bas(bas_i))

       case (0)
          select case(abs(m_bas(bas_i)))
          case(0)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = 0.d0
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = 0.d0
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = 0.d0
          end select

       case(1)
          select case(abs(m_bas(bas_i)))
          case(0)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = -((vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = -((vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3
          case(1)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = -(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = -(vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = -(vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**3)
          end select

       case(2)
          select case(abs(m_bas(bas_i)))
          case(0)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = (-3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = (-3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4
          case(1)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = (-2*sqrt(3.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = (-2*sqrt(3.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (sqrt(3.d0)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4
          case(2)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = -((sqrt(3.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = -((sqrt(3.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = -((sqrt(3.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4)
          end select

       case(3)
          select case(abs(m_bas(bas_i)))
          case(0)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = (3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 5*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = (3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 5*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (-3*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4 - 6*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 5*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
          case(1)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = (sqrt(1.5d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 15*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = (sqrt(1.5d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 15*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (sqrt(1.5d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(11*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 15*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
          case(2)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = (-3*sqrt(15.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = (-3*sqrt(15.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (sqrt(15.d0)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
          case(3)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = (-3*sqrt(2.5d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = (-3*sqrt(2.5d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (-3*sqrt(2.5d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
          end select

       case(4)
          select case(abs(m_bas(bas_i)))
          case(0)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = (5*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2*(3*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 7*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = (5*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2*(3*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 7*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (-5*(3*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1) - 10*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**3 + 7*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**5))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
          case(1)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = (sqrt(2.5d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(3*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 14*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = (sqrt(2.5d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(3*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 14*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = -(sqrt(2.5d0)*(3*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**4 - 27*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 28*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
          case(2)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = (sqrt(5.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 14*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = (sqrt(5.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 14*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (sqrt(5.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(4*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 7*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
          case(3)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = -((sqrt(70.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = -((sqrt(70.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = (sqrt(17.5d0)*(dist_en_wlk(:, basis_fns_cent(bas_i), 1)**2 - 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
          case(4)
             d_bas_funcs_3d_theta_d_x_beta_bf(1, bas_i, :) = -(sqrt(35.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
             d_bas_funcs_3d_theta_d_x_beta_bf(2, bas_i, :) = -(sqrt(35.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
             d_bas_funcs_3d_theta_d_x_beta_bf(3, bas_i, :) = -(sqrt(35.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
          end select

       end select

    end do

  end subroutine d_bas_funcs_3d_theta_d_x_beta_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : 2nd derivative of 3d_theta part of analytic basis functions as a function of backflow transformed coordinate.
    !                 with respect to backflow transformed coordinate (so just like non backflow case, but r is replaced by x)
    !                 This is a rank four array
    !                 The first index labels the dimensions (gamma)
    !                 The second index labels the dimensions (beta)
    !                 The third index labels the basis function
    !                 The fourth index labels the electron
    !
    ! Created       : F. Petruzielo 15 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i
    character(len=max_string_len_rout), save :: lhere = 'd2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf_bld'

    ! header
    if (header_exe) then

       call object_create('d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf')

       call object_needed('ndim')
       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('l_bas')
       call object_needed('m_bas')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('dist_en_wlk')

       return

    endif

    ! allocation
    call object_alloc('d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf', d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf, ndim, ndim, nbasis, nelec)
    d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(:, :, :, :) = 0.d0

    if (maxval(l_bas) > 4) then
       call die(lhere, "Only implemented for l .le. 4")
    end if

    !calculate derivative of theta part of basis functions for the l,m values of each basis function
    do bas_i = 1, nbasis
       select case (l_bas(bas_i))

       case (0)
          select case(abs(m_bas(bas_i)))
          case(0)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 0
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = 0
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = 0
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = 0
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = 0
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = 0
          end select

       case(1)
          select case(abs(m_bas(bas_i)))
          case(0)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = -((vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(-2*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = -((vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -((vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - 2*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = -((vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (-3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5
          case(1)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = -((-2*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -((vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - 2*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = -((vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**5)
          end select

       case(2)
          select case(abs(m_bas(bas_i)))
          case(0)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = (-3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2*(-3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (12*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (-6*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = (-3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - 3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (-6*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
          case(1)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = (-2*sqrt(3.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(-3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (8*sqrt(3.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (-2*sqrt(3.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = (-2*sqrt(3.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - 3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (-2*sqrt(3.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (2*sqrt(3.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(-3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
          case(2)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = (sqrt(3.d0)*(3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (4*sqrt(3.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (4*sqrt(3.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -((sqrt(3.d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - 3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (4*sqrt(3.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = -((sqrt(3.d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**6)
          end select

       case(3)
          select case(abs(m_bas(bas_i)))
          case(0)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = (-3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(2*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**4 + 3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4 + vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 19*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (-3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) - 22*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*((vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 16*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 8*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = (3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 - 2*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**4 + 19*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 - 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4 - vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + 3*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*((vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 16*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 8*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(13*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) - 12*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
          case(1)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = -(sqrt(1.5d0)*(2*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 + vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 59*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2) - (vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 14*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)*(vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (-3*sqrt(1.5d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 24*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (-3*sqrt(1.5d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(11*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) - 14*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = (sqrt(1.5d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 - 2*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**4 + 59*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 - 14*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4 - vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + 13*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (-3*sqrt(1.5d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(11*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) - 14*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (sqrt(1.5d0)*(11*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 56*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 8*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
          case(2)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = (-3*sqrt(15.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(-4*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (15*sqrt(15.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (-3*sqrt(15.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = (-3*sqrt(15.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - 4*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (-3*sqrt(15.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (3*sqrt(15.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(-3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) + 2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
          case(3)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = (3*sqrt(2.5d0)*(4*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (15*sqrt(2.5d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (15*sqrt(2.5d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = (-3*sqrt(2.5d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - 4*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (15*sqrt(2.5d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (-3*sqrt(2.5d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**7)
          end select

       case(4)
          select case(abs(m_bas(bas_i)))
          case(0)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = (-5*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2*(9*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 - 3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**4 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4 + vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(6*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 29*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (15*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2*(-2*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) + 5*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (5*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 14*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = (-5*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2*(-3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 + 9*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**4 - 29*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4 + vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(6*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (5*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 14*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (15*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*(-(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 + 9*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 - 4*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
          case(1)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = -((sqrt(2.5d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(9*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 - 3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**4 + 8*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 11*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4 + vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(6*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 64*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (-6*sqrt(10.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 6*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (3*sqrt(2.5d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*((vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 16*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 11*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -((sqrt(2.5d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(-3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 + 9*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**4 - 64*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 11*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4 + vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(6*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + 8*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (3*sqrt(2.5d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*((vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 16*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 11*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (sqrt(10.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(15*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 25*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
          case(2)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = (sqrt(5.d0)*(-3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 - 2*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 34*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2) + (vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 13*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)*(vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (-2*sqrt(5.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 20*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (-2*sqrt(5.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(8*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) - 13*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = (sqrt(5.d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 - 3*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**4 + 68*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 - 13*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4 - 2*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*(vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + 6*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2)))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (-2*sqrt(5.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(8*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2) - 13*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (sqrt(5.d0)*(4*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)**2 - 29*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2 + 9*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**4))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
          case(3)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = -((sqrt(70.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(-5*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (6*sqrt(70.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = -((sqrt(70.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 5*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -((sqrt(70.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - 5*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = -((sqrt(70.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 5*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = (-3*sqrt(70.d0)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
          case(4)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = (sqrt(35.d0)*(5*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = (3*sqrt(35.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, bas_i, :) = (3*sqrt(35.d0)*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -(sqrt(35.d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - 5*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, bas_i, :) = (3*sqrt(35.d0)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1))/dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8
             d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 3, bas_i, :) = -(sqrt(35.d0)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 - 5*vec_en_xyz_wlk(3, :, basis_fns_cent(bas_i), 1)**2))/(2.*dist_en_wlk(:, basis_fns_cent(bas_i), 1)**8)
          end select

       end select

       !symmetry of partial derivatives
       d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 2, :, :) = d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 1, :, :)
       d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(1, 3, :, :) = d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 1, :, :)
       d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(2, 3, :, :) = d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf(3, 2, :, :)

    end do

  end subroutine d2_bas_funcs_3d_theta_d_x_beta_d_x_gamma_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine bas_funcs_3d_phi_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : 3d_phi part of analytic basis functions evaluated at backflow transformed coordinates of each electron.
    !                 This is a rank two array
    !                 The first index labels the basis function
    !                 The second index labels the electron
    !
    ! Created       : F. Petruzielo 13 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i
    character(len=max_string_len_rout), save :: lhere = 'bas_funcs_3d_phi_bf_bld'

    ! header
    if (header_exe) then

       call object_create('bas_funcs_3d_phi_bf')

       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('m_bas')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('vec_en_xyz_wlk')

       return

    endif

    ! allocation
    call object_alloc('bas_funcs_3d_phi_bf', bas_funcs_3d_phi_bf, nbasis, nelec)
    bas_funcs_3d_phi_bf(:, :) = 0.d0


    if (maxval(m_bas) > 4) then
       call die(lhere, "Only implemented for l .le. 4")
    end if

    do bas_i = 1, nbasis

       select case (m_bas(bas_i))

       case(-4)
          bas_funcs_3d_phi_bf(bas_i, :) = 4*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)
       case(-3)
          bas_funcs_3d_phi_bf(bas_i, :) = 3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1) - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**3
       case(-2)
          bas_funcs_3d_phi_bf(bas_i, :) = 2*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
       case(-1)
          bas_funcs_3d_phi_bf(bas_i, :) = vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
       case (0)
          bas_funcs_3d_phi_bf(bas_i, :) = 1
       case(1)
          bas_funcs_3d_phi_bf(bas_i, :) = vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)
       case(2)
          bas_funcs_3d_phi_bf(bas_i, :) = vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2
       case(3)
          bas_funcs_3d_phi_bf(bas_i, :) = vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**3 - 3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2
       case(4)
          bas_funcs_3d_phi_bf(bas_i, :) = vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**4 - 6*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**4
       end select

    end do

  end subroutine bas_funcs_3d_phi_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d_bas_funcs_3d_phi_d_x_beta_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : derivative of 3d_phi part of analytic basis functions as a function of backflow transformed coordinate.
    !                 with respect to backflow transformed coordinate (so just like non backflow case, but r is replaced by x)
    !                 This is a rank three array
    !                 The first index labels the dimension (beta)
    !                 The second index labels the basis function
    !                 The third index labels the electron
    !
    ! Created       : F. Petruzielo 14 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i
    character(len=max_string_len_rout), save :: lhere = 'd_bas_funcs_3d_phi_d_x_beta_bf_bld'

    ! header
    if (header_exe) then

       call object_create('d_bas_funcs_3d_phi_d_x_beta_bf')

       call object_needed('ndim')
       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('m_bas')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('dist_en_wlk')

       return

    endif

    ! allocation
    call object_alloc('d_bas_funcs_3d_phi_d_x_beta_bf', d_bas_funcs_3d_phi_d_x_beta_bf, ndim, nbasis, nelec)
    d_bas_funcs_3d_phi_d_x_beta_bf(:, :, :) = 0.d0


    if (maxval(m_bas) > 4) then
       call die(lhere, "Only implemented for l .le. 4")
    end if

    do bas_i = 1, nbasis

       select case (m_bas(bas_i))

       case(-4)
          d_bas_funcs_3d_phi_d_x_beta_bf(1, bas_i, :) = 12*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1) - 4*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**3
          d_bas_funcs_3d_phi_d_x_beta_bf(2, bas_i, :) = 4*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**3 - 3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)
       case(-3)
          d_bas_funcs_3d_phi_d_x_beta_bf(1, bas_i, :) = 6*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
          d_bas_funcs_3d_phi_d_x_beta_bf(2, bas_i, :) = 3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)
       case(-2)
          d_bas_funcs_3d_phi_d_x_beta_bf(1, bas_i, :) = 2*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
          d_bas_funcs_3d_phi_d_x_beta_bf(2, bas_i, :) = 2*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)
       case(-1)
          d_bas_funcs_3d_phi_d_x_beta_bf(1, bas_i, :) = 0
          d_bas_funcs_3d_phi_d_x_beta_bf(2, bas_i, :) = 1
       case (0)
          d_bas_funcs_3d_phi_d_x_beta_bf(1, bas_i, :) = 0
          d_bas_funcs_3d_phi_d_x_beta_bf(2, bas_i, :) = 0
       case(1)
          d_bas_funcs_3d_phi_d_x_beta_bf(1, bas_i, :) = 1
          d_bas_funcs_3d_phi_d_x_beta_bf(2, bas_i, :) = 0
       case(2)
          d_bas_funcs_3d_phi_d_x_beta_bf(1, bas_i, :) = 2*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)
          d_bas_funcs_3d_phi_d_x_beta_bf(2, bas_i, :) = -2*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
       case(3)
          d_bas_funcs_3d_phi_d_x_beta_bf(1, bas_i, :) = 3*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)
          d_bas_funcs_3d_phi_d_x_beta_bf(2, bas_i, :) = -6*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
       case(4)
          d_bas_funcs_3d_phi_d_x_beta_bf(1, bas_i, :) = 4*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**3 - 3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)
          d_bas_funcs_3d_phi_d_x_beta_bf(2, bas_i, :) = 4*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)*(-3*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)
       end select

    end do

  end subroutine d_bas_funcs_3d_phi_d_x_beta_bf_bld

  ! ==============================================================================

  ! ==============================================================================

  subroutine d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf_bld
    ! ------------------------------------------------------------------------------
    ! Description   : derivative of 3d_phi part of analytic basis functions as a function of backflow transformed coordinate.
    !                 with respect to backflow transformed coordinate (so just like non backflow case, but r is replaced by x)
    !                 This is a rank three array
    !                 The first index labels the dimension (beta)
    !                 The second index labels the basis function
    !                 The third index labels the electron
    !
    ! Created       : F. Petruzielo 14 Apr, 2009
    ! ------------------------------------------------------------------------------
    use all_modules_mod
    implicit none

    !local
    integer :: bas_i
    character(len=max_string_len_rout), save :: lhere = 'd2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf_bld'

    ! header
    if (header_exe) then

       call object_create('d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf')

       call object_needed('ndim')
       call object_needed('nbasis')
       call object_needed('nelec')
       call object_needed('m_bas')
       !testing
       call object_needed('basis_fns_cent')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('dist_en_wlk')

       return

    endif

    ! allocation
    call object_alloc('d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf', d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf, ndim, ndim, nbasis, nelec)
    d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(:, :, :, :) = 0.d0


    if (maxval(m_bas) > 4) then
       call die(lhere, "Only implemented for l .le. 4")
    end if

    do bas_i = 1, nbasis

       select case (m_bas(bas_i))

       case(-4)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 24*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = 12*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -24*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
       case(-3)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 6*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = 6*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -6*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
       case(-2)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 0
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = 2
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = 0
       case(-1)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 0
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = 0
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = 0
       case (0)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 0
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = 0
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = 0
       case(1)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 0
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = 0
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = 0
       case(2)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 2
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = 0
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -2
       case(3)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 6*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = -6*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = -6*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)
       case(4)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 1, bas_i, :) = 12*(vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 - vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, bas_i, :) = -24*vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)*vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)
          d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 2, bas_i, :) = 12*(-vec_en_xyz_wlk(1, :, basis_fns_cent(bas_i), 1)**2 + vec_en_xyz_wlk(2, :, basis_fns_cent(bas_i), 1)**2)
       end select

    end do

    !symmetry of partial derivatives
    d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(1, 2, :, :) = d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf(2, 1, :, :)

  end subroutine d2_bas_funcs_3d_phi_d_x_beta_d_x_gamma_bf_bld

  ! ==============================================================================

  ! =======================================================================
  function factorial_array(i)
    !-------------------------------------------------------------------------
    ! Description   : factorial function
    !                 Calculate the factorial of each element of a one dimensional array
    !
    ! Created       : F. Petruzielo 13 Apr, 2009
    !-------------------------------------------------------------------------
    implicit none

    ! input/output
    integer, dimension(:) :: i
    integer, dimension(size(i)) :: factorial_array

    !local
    integer :: j

    factorial_array(:) = 1

    do j = 1, size(i)
       factorial_array(j) = factorial_scalar(i(j))
    enddo

  end function factorial_array

  ! =======================================================================

  ! ========================================================================
  function factorial_scalar(i)
    !-------------------------------------------------------------------------
    ! Description   : factorial function
    !                 Calculate the factorial of a scalar
    !
    ! Created       : F. Petruzielo 13 Apr, 2009 (from J. Toulouse)
    !-------------------------------------------------------------------------
    implicit none

    ! input/output
    integer :: i
    integer :: factorial_scalar

    !local
    integer :: j
    character(len=max_string_len_rout), save :: lhere = 'factorial_scalar'

    if (i < 0) then
       call die (lhere, 'factorial with negative argument')
    endif

    factorial_scalar = 1

    do j = 2,i
       factorial_scalar = factorial_scalar * j
    enddo

  end function factorial_scalar
  ! ========================================================================



end module backflow_basis_mod
