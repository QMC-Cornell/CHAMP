!see comprehensive notes of F. Petruzielo
module backflow_mod
  use parser_tools_mod, only: word, get_next_word, get_next_value, get_next_value_list_object
  use constants_mod    !give access to various common blocks
  use basic_tools_mod, only: require, die
  use strings_tools_mod !require operator overloading of + to concatenate strings
  use objects_mod, only: object_modified, object_create, object_needed, object_alloc, object_provide, object_associate, object_provide_in_node
  use nodes_mod, only: header_exe
  use variables_mod, only: spin_nb
  !begin objects
  use electrons_mod, only: dist_en_wlk, vec_en_xyz_wlk, vec_ee_xyz_wlk, dist_ee_wlk, coord_elec_wlk

  !Declaration of global variables and default values

  integer                                        :: isc_bf
  real(dp)                                       :: scalek_bf
  logical                                        :: all_elec

  real(dp), allocatable                          :: xi(:, :, :)
  real(dp), allocatable                          :: d_xi_j_beta_d_r_i_alpha(:, :, :, :, :)
  real(dp), allocatable                          :: d_xi_j_beta_d_r_i_alpha_num(:, :, :, :, :)
  real(dp), allocatable                          :: lap_i_xi_j_beta(:,:,:,:)
  real(dp), allocatable                          :: lap_i_xi_j_beta_num(:,:,:,:)

  integer                                        :: threshold_l_nb
  real(dp), allocatable                          :: threshold_l(:)
  real(dp), allocatable                          :: smooth_cutoff_g(:,:,:)
  real(dp), allocatable                          :: d_smooth_cutoff_g_d_r(:,:,:)
  real(dp), allocatable                          :: d2_smooth_cutoff_g_d_r_2(:,:,:)

  logical                                        :: en_bf
  real(dp), allocatable                          :: xi_en(:,:,:)
  integer                                        :: order_mu_bf
  integer                                        :: read_d_param_mu_nb
  logical                                        :: read_d_param
  real(dp), allocatable                          :: read_d_param_mu(:)
  real(dp), allocatable                          :: mu(:,:,:)
  real(dp), allocatable                          :: d_param_mu(:,:,:)
  logical                                        :: mu_spin_dependent
  real(dp), allocatable                          :: asymp_mu_pseudo(:,:)
  real(dp), allocatable                          :: d_xi_en_j_beta_d_r_i_alpha(:,:,:,:)
  real(dp), allocatable                          :: d_xi_en_j_beta_d_r_i_alpha_num(:,:,:,:)
  real(dp), allocatable                          :: d_mu_d_r_j_inuc(:,:,:)
  real(dp), allocatable                          :: d_mu_d_r_j_alpha(:,:,:,:)
  real(dp), allocatable                          :: d2_mu_d_r_j_inuc_2(:,:,:)
  real(dp), allocatable                          :: lap_i_xi_en_j_beta(:,:,:)
  real(dp), allocatable                          :: lap_i_xi_en_j_beta_num(:,:,:)
  real(dp), allocatable                          :: lap_j_mu(:,:,:)

  logical                                        :: ee_bf
  real(dp), allocatable                          :: xi_ee(:,:,:)
  integer                                        :: order_eta_bf
  integer                                        :: read_c_param_eta_nb
  logical                                        :: read_c_param
  real(dp), allocatable                          :: read_c_param_eta(:)
  real(dp), allocatable                          :: eta(:,:,:)
  real(dp), allocatable                          :: c_param_eta(:,:)
  integer                                        :: eta_spin_dependence
  real(dp), allocatable                          :: asymp_eta(:)
  real(dp), allocatable                          :: d_xi_ee_j_beta_d_r_i_alpha(:,:,:,:,:)
  real(dp), allocatable                          :: d_xi_ee_j_beta_d_r_i_alpha_num(:,:,:,:,:)
  real(dp), allocatable                          :: d_eta_d_r_j_alpha_first(:,:,:,:)
  real(dp), allocatable                          :: d_eta_d_r_i_alpha_second(:,:,:,:)
  real(dp), allocatable                          :: d_eta_d_r_ji(:,:,:)
  real(dp), allocatable                          :: d2_eta_d_r_ji_2(:,:,:)
  real(dp), allocatable                          :: lap_i_xi_ee_j_beta(:,:,:,:)
  real(dp), allocatable                          :: lap_i_xi_ee_j_beta_num(:,:,:,:)
  real(dp), allocatable                          :: lap_j_eta_first(:,:,:)
  real(dp), allocatable                          :: lap_i_eta_second(:,:,:)

  logical                                        :: een_phi_bf
  real(dp), allocatable                          :: xi_een_phi(:,:,:)
  integer                                        :: order_phi_bf
  integer                                        :: read_a_param_phi_nb
  logical                                        :: read_a_param
  real(dp), allocatable                          :: read_a_param_phi(:)
  real(dp), allocatable                          :: phi(:,:,:,:)
  real(dp), allocatable                          :: a_param_phi(:,:,:)
  real(dp), allocatable                          :: a_param_phi_anti_parallel_cond(:,:)
  real(dp), allocatable                          :: a_param_phi_parallel_cond(:,:)
  integer                                        :: phi_spin_dependence
  integer                                        :: nb_a_param_phi
  integer, allocatable                           :: dep_a_param_phi_parallel(:)
  integer, allocatable                           :: dep_a_param_phi_anti_parallel(:)
  real(dp), allocatable                          :: d_xi_een_phi_j_beta_d_r_i_alpha(:,:,:,:,:)
  real(dp), allocatable                          :: d_xi_een_phi_j_beta_d_r_i_alpha_num(:,:,:,:,:)
  real(dp), allocatable                          :: d_phi_d_r_j_alpha_first(:,:,:,:,:)
  real(dp), allocatable                          :: d_phi_d_r_i_alpha_second(:,:,:,:,:)
  real(dp), allocatable                          :: d_phi_d_r_ji(:,:,:,:)
  real(dp), allocatable                          :: d_phi_d_r_j_inuc(:,:,:,:)
  real(dp), allocatable                          :: lap_i_xi_een_phi_j_beta(:,:,:,:)
  real(dp), allocatable                          :: lap_i_xi_een_phi_j_beta_num(:,:,:,:)
  real(dp), allocatable                          :: lap_j_phi_first(:,:,:,:)
  real(dp), allocatable                          :: lap_i_phi_second(:,:,:,:)
  real(dp), allocatable                          :: d2_phi_d_r_ji_2(:,:,:,:)
  real(dp), allocatable                          :: d2_phi_d_r_j_inuc_2(:,:,:,:)
  real(dp), allocatable                          :: d2_phi_d_r_j_inuc_d_r_ji(:,:,:,:)

  logical                                        :: een_theta_bf
  real(dp), allocatable                          :: xi_een_theta(:,:,:)
  integer                                        :: order_theta_bf
  integer                                        :: read_b_param_theta_nb
  logical                                        :: read_b_param
  real(dp), allocatable                          :: read_b_param_theta(:)
  real(dp), allocatable                          :: theta(:,:,:,:)
  real(dp), allocatable                          :: b_param_theta(:,:,:)
  real(dp), allocatable                          :: b_param_theta_cond(:,:)
  integer                                        :: theta_spin_dependence
  integer                                        :: nb_b_param_theta
  integer, allocatable                           :: dep_b_param_theta(:)
  real(dp), allocatable                          :: d_xi_een_theta_j_beta_d_r_i_alpha(:,:,:,:,:)
  real(dp), allocatable                          :: d_xi_een_theta_j_beta_d_r_i_alpha_num(:,:,:,:,:)
  real(dp), allocatable                          :: d_theta_d_r_j_alpha_first(:,:,:,:,:)
  real(dp), allocatable                          :: d_theta_d_r_i_alpha_second(:,:,:,:,:)
  real(dp), allocatable                          :: d_theta_d_r_ji(:,:,:,:)
  real(dp), allocatable                          :: d_theta_d_r_j_inuc(:,:,:,:)
  real(dp), allocatable                          :: lap_i_xi_een_theta_j_beta(:,:,:,:)
  real(dp), allocatable                          :: lap_i_xi_een_theta_j_beta_num(:,:,:,:)
  real(dp), allocatable                          :: lap_j_theta_first(:,:,:,:)
  real(dp), allocatable                          :: lap_i_theta_second(:,:,:,:)
  real(dp), allocatable                          :: d2_theta_d_r_ji_2(:,:,:,:)
  real(dp), allocatable                          :: d2_theta_d_r_j_inuc_2(:,:,:,:)
  real(dp), allocatable                          :: d2_theta_d_r_j_inuc_d_r_ji(:,:,:,:)

  real(dp), allocatable                          :: scaled_dist_een_ee_wlk(:,:,:)
  real(dp), allocatable                          :: scaled_dist_een_en_wlk(:,:,:)
  real(dp), allocatable                          :: scaled_dist_ee_wlk(:,:,:)
  real(dp), allocatable                          :: scaled_dist_en_wlk(:,:,:)
  real(dp)                                       :: asymp_scaled_dist_two_body
  real(dp), allocatable                          :: d_scaled_dist_en_wlk_d_r(:,:,:)
  real(dp), allocatable                          :: d_scaled_dist_ee_wlk_d_r(:,:,:)
  real(dp), allocatable                          :: d2_scaled_dist_ee_wlk_d_r_2(:,:,:)
  real(dp), allocatable                          :: d2_scaled_dist_en_wlk_d_r_2(:,:,:)
  real(dp), allocatable                          :: d_scaled_dist_een_ee_wlk_d_r(:,:,:)
  real(dp), allocatable                          :: d_scaled_dist_een_en_wlk_d_r(:,:,:)
  real(dp), allocatable                          :: d2_scaled_dist_een_ee_wlk_d_r_2(:,:,:)
  real(dp), allocatable                          :: d2_scaled_dist_een_en_wlk_d_r_2(:,:,:)

  !end objects

contains

  !===========================================================================
  subroutine backflow_menu
    !---------------------------------------------------------------------------
    ! Description : menu for backflow transformation on the wave function
    !
    ! Created     : F. Petruzielo, 17 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need scalek, isc, nctype
    include 'commons.h'

    ! local
    character(len=max_string_len_rout), save :: lhere = 'backflow_menu'

    ! begin
    write(6,*)
    write(6,'(a)') 'Beginning of backflow menu ---------------------------------------------------------------------------'

    !Add some defaults to the tree
    !object associate: add a pointer to this object so that it knows its name (called automatically when you call object_alloc for an array).
    !object_modified: when an object has been modified, validate it and invalidate all objects depending on it also catalog object if necessary
    en_bf = .true.
    call object_modified('en_bf')
    call object_associate('en_bf', en_bf)

    ee_bf = .true.
    call object_modified('ee_bf')
    call object_associate('ee_bf', ee_bf)

    een_phi_bf = .true.
    call object_modified('een_phi_bf')
    call object_associate('een_phi_bf', een_phi_bf)

    een_theta_bf = .true.
    call object_modified('een_theta_bf')
    call object_associate('een_theta_bf', een_theta_bf)

    mu_spin_dependent = .true.
    call object_modified('mu_spin_dependent')
    call object_associate('mu_spin_dependent', mu_spin_dependent)

    eta_spin_dependence = 3
    call object_modified('eta_spin_dependence')
    call object_associate('eta_spin_dependence', eta_spin_dependence)

    phi_spin_dependence = 3
    call object_modified('phi_spin_dependence')
    call object_associate('phi_spin_dependence', phi_spin_dependence)

    theta_spin_dependence = 3
    call object_modified('theta_spin_dependence')
    call object_associate('theta_spin_dependence', theta_spin_dependence)

    read_d_param=.false.
    call object_modified('read_d_param')
    call object_associate('read_d_param', read_d_param)

    read_d_param_mu_nb = 0
    call object_modified('read_d_param_mu_nb')
    call object_associate('read_d_param_mu_nb', read_d_param_mu_nb)

    read_c_param=.false.
    call object_modified('read_c_param')
    call object_associate('read_c_param', read_c_param)

    read_c_param_eta_nb = 0
    call object_modified('read_c_param_eta_nb')
    call object_associate('read_c_param_eta_nb', read_c_param_eta_nb)

    read_a_param=.false.
    call object_modified('read_a_param')
    call object_associate('read_a_param', read_a_param)

    read_a_param_phi_nb = 0
    call object_modified('read_a_param_phi_nb')
    call object_associate('read_a_param_phi_nb', read_a_param_phi_nb)

    read_b_param=.false.
    call object_modified('read_b_param')
    call object_associate('read_b_param', read_b_param)

    read_b_param_theta_nb = 0
    call object_modified('read_b_param_theta_nb')
    call object_associate('read_b_param_theta_nb', read_b_param_theta_nb)

    !use object provide because backflow_menu is not in the tree!
    call object_provide('scalek')
    scalek_bf = scalek(1)
    call object_modified('scalek_bf')
    call object_associate('scalek_bf',scalek_bf)

    call object_provide('isc')
    isc_bf = isc
    call object_modified('isc_bf')
    call object_associate('isc_bf',isc_bf)

    ! loop over menu lines
    do
       call get_next_word (word)

       select case (trim(word))
       case ('help')
          write(6,'(a)') 'HELP for backflow menu'
          write(6,'(a)') ': backflow'
          write(6,'(a)') ': en_bf = [logical] : backflow transformation includes electron-nuclear backflow (default=true)'
          write(6,'(a)') ': ee_bf = [logical] : backflow transformation includes electron-electron backflow (default=true)'
          write(6,'(a)') ': een_phi_bf = [logical] : backflow transformation includes electron-electron-nuclear phi backflow (default=true)'
          write(6,'(a)') ': een_theta_bf = [logical] : backflow transformation includes electron-electron-nuclear theta backflow (default=true)'
          write(6,'(a)') ': order_eta_bf = [integer] : order of expansion for eta function in the electron-electron backflow'
          write(6,'(a)') ': order_mu_bf = [integer] : order of expansion for mu function in the electron-nuclear backflow'
          write(6,'(a)') ': order_phi_bf = [integer] : order of expansion for phi function in the electron-electron-nuclear backflow'
          write(6,'(a)') ': order_theta_bf = [integer] : order of expansion for theta function in the electron-electron-nuclear backflow'
          write(6,'(a)') ': isc_bf = [integer] : form of the scaling functions used in backflow functions mu, eta, phi, and theta. (default=isc)'
          write(6,'(a)') ': scalek_bf = [real(dp)] : scale factor used in the scaling functions. (default=scalek)'
          write(6,'(a)') ': mu_spin_dependent = [logical] : mu is spin dependent. (default=true)'
          write(6,'(a)') ': eta_spin_dependence = [integer] :  1: up_up = down_down = up_down '
          write(6,'(a)') ':                                    2: up_up = down_down != up_down '
          write(6,'(a)') ':                                    3: up_up != down_down != up_down (default=3)'
          write(6,'(a)') ': phi_spin_dependence = [integer] :  1: up_up = down_down = up_down '
          write(6,'(a)') ':                                    2: up_up = down_down != up_down '
          write(6,'(a)') ':                                    3: up_up != down_down != up_down (default=3)'
          write(6,'(a)') ': theta_spin_dependence = [integer] :  1: up_up = down_down = up_down '
          write(6,'(a)') ':                                    2: up_up = down_down != up_down '
          write(6,'(a)') ':                                    3: up_up != down_down != up_down (default=3)'
          write(6,'(a)') ': all_elec = [logical] : this is an all electron calculation'
          write(6,'(a)') ': threshold_l = [real(dp), rank-1 array] : threshold_l threshold_l(1) ... threshold_l(nctype) end'
          write(6,'(a)') '                This specifies the distance from a nucleus beyond which the smooth_cutoff_g is equal to 1.'
          write(6,'(a)') ': read_d_param_mu = [real(dp), rank-1 array] : Note that this will be reshaped into a rank-3 array of dimension (order_mu_bf + 1, nctype, spin_nb).'
          write(6,'(a)') '                                               Note the order_mu_bf + 1 arises because of the definition of mu.'
          write(6,'(a)') '                                               The input for this option is the following: (for example order_mu_bf=1 nctype=2 spin_nb=2)'
          write(6,'(a)') '                                               read_d_param_mu (1,1,1)  (2,1,1)  (1,2,1)  (2,2,1) (1,1,2)  (2,1,2)  (1,2,2)  (2,2,2) end'
          write(6,'(a)') '                                               One final note: if mu_spin_dependent=.false. then set final index to 1.'
          write(6,'(a)') ': read_c_param_eta = [real(dp), rank-1 array] : Note that this will be reshaped into a rank-2 array of dimension (order_eta_bf + 1, eta_spin_dependence).'
          write(6,'(a)') '                                               Note the order_eta_bf + 1 arises because of the definition of eta.'
          write(6,'(a)') '                                               The input for this option is the following: (for example order_eta_bf=1, eta_spin_dependence=2 )'
          write(6,'(a)') '                                               read_d_param_mu (1,1)  (2,1)  (1,2)  (2,2) end'
          write(6,'(a)') ': read_a_param_phi = [real(dp), rank-1 array] : Note that this will be reshaped into a rank-3 array of dimension (nb_a_param_phi, ,nctype, phi_spin_dependence).'
          write(6,'(a)') '                                                nb_a_param_phi comes from the definition of phi and depends on order_phi_bf. This is identical to the case of the three body jastrow. '
          write(6,'(a)') '                                               The input for this option is the following: (for example order_phi_bf=3 => nb_a_param_phi=6, nctype=2 phi_spin_dependence=1 )'
          write(6,'(a)') '                                               read_a_param_phi (1,1,1)  (2,1,1)  (3,1,1)  (4,1,1)  (5,1,1) (6,1,1) (1,2,1)  (2,2,1) ... end'
          write(6,'(a)') ': read_b_param_theta = [real(dp), rank-1 array] : Note that this will be reshaped into a rank-3 array of dimension (nb_b_param_theta, ,nctype, theta_spin_dependence).'
          write(6,'(a)') '                                                nb_b_param_theta comes from the definition of theta and depends on order_theta_bf. This is identical to the case of the three body jastrow. '
          write(6,'(a)') '                                               The input for this option is the following: (for example order_theta_bf=3 => nb_b_param_theta=6, nctype=2 theta_spin_dependence=1 )'
          write(6,'(a)') '                                               read_b_param_theta (1,1,1)  (2,1,1)  (3,1,1)  (4,1,1)  (5,1,1) (6,1,1) (1,2,1)  (2,2,1) ... end'
          write(6,'(a)') ': end'


       case('order_eta_bf')
          call get_next_value(order_eta_bf)
          call require(lhere, 'order_eta_bf > 0', order_eta_bf .gt. 0)
          call object_modified('order_eta_bf')
          call object_associate('order_eta_bf',order_eta_bf)

       case('order_mu_bf')
          call get_next_value(order_mu_bf)
          call require(lhere, 'order_mu_bf > 0', order_mu_bf .gt. 0)
          call object_modified('order_mu_bf')
          call object_associate('order_mu_bf',order_mu_bf)

       case('order_phi_bf')
          call get_next_value(order_phi_bf)
          call require(lhere, 'order_phi_bf > 1', order_phi_bf .gt. 1)
          call object_modified('order_phi_bf')
          call object_associate('order_phi_bf',order_phi_bf)

       case('order_theta_bf')
          call get_next_value(order_theta_bf)
          call require(lhere, 'order_theta_bf > 1', order_theta_bf .gt. 1)
          call object_modified('order_theta_bf')
          call object_associate('order_theta_bf',order_theta_bf)

       case('scalek_bf')
          call get_next_value(scalek_bf)
          call require(lhere, 'scalek_bf > 0', scalek_bf .gt. 0.d0)
          call object_modified('scalek_bf')

       case('isc_bf')
          call get_next_value(isc_bf)
          call require(lhere, 'isc_bf > 0', isc_bf .gt. 0)
          call object_modified('isc_bf')

       case('mu_spin_dependent')
          call get_next_value(mu_spin_dependent)
          call object_modified('mu_spin_dependent')

       case('eta_spin_dependence')
          call get_next_value(eta_spin_dependence)
          call object_modified('eta_spin_dependence')

       case('phi_spin_dependence')
          call get_next_value(phi_spin_dependence)
          call object_modified('phi_spin_dependence')

       case('theta_spin_dependence')
          call get_next_value(theta_spin_dependence)
          call object_modified('theta_spin_dependence')

       case('en_bf')
          call get_next_value(en_bf)
          call object_modified('en_bf')

       case('ee_bf')
          call get_next_value(ee_bf)
          call object_modified('ee_bf')

       case('een_phi_bf')
          call get_next_value(een_phi_bf)
          call object_modified('een_phi_bf')

       case('een_theta_bf')
          call get_next_value(een_theta_bf)
          call object_modified('een_theta_bf')

       case('all_elec')
          call get_next_value(all_elec)
          call object_modified('all_elec')
          call object_associate('all_elec', all_elec)

       case ('read_d_param_mu')
          call get_next_value_list_object('read_d_param_mu', read_d_param_mu, read_d_param_mu_nb)
          call object_modified('read_d_param_mu')
          call object_modified('read_d_param_mu_nb')
          read_d_param=.true.
          call object_modified('read_d_param')

       case ('read_c_param_eta')
          call get_next_value_list_object('read_c_param_eta', read_c_param_eta, read_c_param_eta_nb)
          call object_modified('read_c_param_eta')
          call object_modified('read_c_param_eta_nb')
          read_c_param=.true.
          call object_modified('read_c_param')

       case ('read_a_param_phi')
          call get_next_value_list_object('read_a_param_phi', read_a_param_phi, read_a_param_phi_nb)
          call object_modified('read_a_param_phi')
          call object_modified('read_a_param_phi_nb')
          read_a_param=.true.
          call object_modified('read_a_param')

       case ('read_b_param_theta')
          call get_next_value_list_object('read_b_param_theta', read_b_param_theta, read_b_param_theta_nb)
          call object_modified('read_b_param_theta')
          call object_modified('read_b_param_theta_nb')
          read_b_param=.true.
          call object_modified('read_b_param')

       case ('threshold_l')
          call get_next_value_list_object('threshold_l', threshold_l, threshold_l_nb)
          call object_modified('threshold_l')
          call object_modified('threshold_l_nb')
          call object_associate('threshold_l_nb', threshold_l_nb)
          call object_provide('nctype')
          call require(lhere, 'Size of threshold_l = nctype', threshold_l_nb .eq. nctype)

       case ('end')
          exit

       case default
          call die (lhere, 'unknown keyword >'+trim(word)+'<.')

       end select

    enddo

    write(6,'(a)') 'End of backflow menu ---------------------------------------------------------------------------------'

  end subroutine backflow_menu

  !======================================================================================

  !==============================================================================================

  subroutine scaled_dist_en_wlk_bld

    !---------------------------------------------------------------------------
    ! Description : build scaled_dist_en_wlk(:,:,:). This is a rank three array.
    !               The first index is the electron.
    !               The second index is the nucleus
    !               The third index is the walker
    !               This returns the two body scaled distances between electrons and nuclei
    !               for each walker.
    ! Created     : F. Petruzielo, 18 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('scaled_dist_en_wlk')

       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed('nwalk')
       call object_needed ('dist_en_wlk')
       !isc_bf, scalek_bf are needed in a subroutine that is called by this build routine so make sure up to date.
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('scaled_dist_en_wlk', scaled_dist_en_wlk, nelec, ncent, nwalk)
    scaled_dist_en_wlk(:,:,:) = 0.d0

    call scaling_func_two_body(dist_en_wlk,scaled_dist_en_wlk)

  end subroutine scaled_dist_en_wlk_bld


  !==============================================================================================

  !==============================================================================================

  subroutine d_scaled_dist_en_wlk_d_r_bld

    !---------------------------------------------------------------------------
    ! Description : build d_scaled_dist_en_wlk_d_r(:,:,:). This is a rank three array.
    !               The first index is the electron.
    !               The second index is the nucleus
    !               The third index is the walker
    !               This returns the derivative of the scaled distances between electrons and nuclei
    !               for each walker.
    ! Created     : F. Petruzielo, 3 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('d_scaled_dist_en_wlk_d_r')

       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed('nwalk')
       call object_needed ('dist_en_wlk')
       !isc_bf, scalek_bf are needed in a subroutine that is called by this build routine so make sure up to date.
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('d_scaled_dist_en_wlk_d_r', d_scaled_dist_en_wlk_d_r, nelec, ncent, nwalk)
    d_scaled_dist_en_wlk_d_r(:,:,:) = 0.d0

    !calculate derivatives of scaled distances
    call d_scaling_func_two_body(dist_en_wlk, d_scaled_dist_en_wlk_d_r)

  end subroutine d_scaled_dist_en_wlk_d_r_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d2_scaled_dist_en_wlk_d_r_2_bld

    !---------------------------------------------------------------------------
    ! Description : build d2_scaled_dist_en_wlk_d_r_2(:,:,:). This is a rank three array.
    !               The first index is the electron.
    !               The second index is the center.
    !               The third index is the walker
    !               This returns the second derivative of the scaled distances between electrons and nuclei
    !               for each walker.
    ! Created     : F. Petruzielo, 20 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('d2_scaled_dist_en_wlk_d_r_2')

       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed('dist_en_wlk')
       call object_needed('ncent')
       !Make sure scalek_bf, isc_bf are up to date because needed in a subroutine that is called
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('d2_scaled_dist_en_wlk_d_r_2', d2_scaled_dist_en_wlk_d_r_2, nelec, ncent, nwalk)
    d2_scaled_dist_en_wlk_d_r_2(:,:,:) = 0.d0

    !calculate second derivative of scaled distances
    call d2_scaling_func_two_body(dist_en_wlk, d2_scaled_dist_en_wlk_d_r_2)

  end subroutine d2_scaled_dist_en_wlk_d_r_2_bld

  !==============================================================================================

  !==============================================================================================

  subroutine scaled_dist_ee_wlk_bld

    !---------------------------------------------------------------------------
    ! Description : build scaled_dist_ee_wlk(:,:,:). This is a rank three array.
    !               This is the two body scaled distance between electrons
    !               The first index is the electron.
    !               The second index is the electron
    !               The third index is the walker
    !
    ! Created     : F. Petruzielo, 22 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('scaled_dist_ee_wlk')

       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed ('dist_ee_wlk')
       !Make sure scalek_bf, isc_bf are up to date because needed in a subroutine that is called
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('scaled_dist_ee_wlk', scaled_dist_ee_wlk, nelec, nelec, nwalk)
    scaled_dist_ee_wlk(:,:,:) = 0.d0

    !calculate scaled distances
    call scaling_func_two_body(dist_ee_wlk,scaled_dist_ee_wlk)

  end subroutine scaled_dist_ee_wlk_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_scaled_dist_ee_wlk_d_r_bld

    !---------------------------------------------------------------------------
    ! Description : build d_scaled_dist_ee_wlk_d_r(:,:,:). This is a rank three array.
    !               The first index is the electron.
    !               The second index is the electron.
    !               The third index is the walker
    !               This returns the derivative of the scaled distances between electrons
    !               for each walker.
    ! Created     : F. Petruzielo, 5 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('d_scaled_dist_ee_wlk_d_r')

       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed ('dist_ee_wlk')
       !Make sure scalek_bf, isc_bf are up to date because needed in a subroutine that is called
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('d_scaled_dist_ee_wlk_d_r', d_scaled_dist_ee_wlk_d_r, nelec, nelec, nwalk)
    d_scaled_dist_ee_wlk_d_r(:,:,:) = 0.d0

    !calculate derivatives of scaled distances
    call d_scaling_func_two_body(dist_ee_wlk, d_scaled_dist_ee_wlk_d_r)

  end subroutine d_scaled_dist_ee_wlk_d_r_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d2_scaled_dist_ee_wlk_d_r_2_bld

    !---------------------------------------------------------------------------
    ! Description : build d2_scaled_dist_ee_wlk_d_r_2(:,:,:). This is a rank three array.
    !               The first index is the electron.
    !               The second index is the electron.
    !               The third index is the walker
    !               This returns the second derivative of the scaled distances between electrons
    !               for each walker.
    ! Created     : F. Petruzielo, 17 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('d2_scaled_dist_ee_wlk_d_r_2')

       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed ('dist_ee_wlk')
       !Make sure scalek_bf, isc_bf are up to date because needed in a subroutine that is called
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('d2_scaled_dist_ee_wlk_d_r_2', d2_scaled_dist_ee_wlk_d_r_2, nelec, nelec, nwalk)
    d2_scaled_dist_ee_wlk_d_r_2(:,:,:) = 0.d0

    !calculate second derivative of scaled distances
    call d2_scaling_func_two_body(dist_ee_wlk, d2_scaled_dist_ee_wlk_d_r_2)

  end subroutine d2_scaled_dist_ee_wlk_d_r_2_bld

  !==============================================================================================

  !==============================================================================================

  subroutine asymp_scaled_dist_two_body_bld

    !---------------------------------------------------------------------------
    ! Description : build asymp_scaled_dist_two_body. This is a real(dp).
    !               Calculate the value of the two body scaling functions evaluated at infinity
    !
    ! Created     : F. Petruzielo, 19 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !local
    character(len=max_string_len_rout), save :: lhere = 'asymp_scaled_dist_two_body_bld'

    ! header
    if (header_exe) then

       call object_create('asymp_scaled_dist_two_body')

       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    select case(isc_bf)
    case(4) ! r / (1 + r * k)
       asymp_scaled_dist_two_body = 1.d0 / scalek_bf
    case default
       call die(lhere,'Only isc = 4 is implemented.')
    end select

    !give scalar object a name
    call object_associate('asymp_scaled_dist_two_body', asymp_scaled_dist_two_body)

  end subroutine asymp_scaled_dist_two_body_bld

  !==============================================================================================

  !==============================================================================================

  subroutine scaling_func_two_body(distance_dummy, scaled_distance_dummy)
    !---------------------------------------------------------------------------
    ! Description : scaling_func_two_body. This subroutine calculates the scaled
    !               distances of rank-3 array of distances. The scaling form is specified
    !               by isc_bf and the scaling factor is specified by scalek_bf. Note that isc_bf and scalek_bf
    !               must be called with object_needed in any routine calling this routine such that they are valid
    !
    ! Created     : F. Petruzielo, 18 Jul 2008
    !---------------------------------------------------------------------------


    implicit none

    ! output
    real(dp), dimension(:,:,:), intent(out) :: scaled_distance_dummy

    ! input
    real(dp), dimension(:,:,:), intent(in) :: distance_dummy
    character(len=max_string_len_rout), save :: lhere = 'scaling_func_two_body'

    select case(isc_bf)
    case(4) ! r / (1 + r * k)
       scaled_distance_dummy = distance_dummy / (1.d0 + scalek_bf * distance_dummy)
    case default
       call die(lhere,'Only isc = 4 is implemented.')
    end select

  end subroutine scaling_func_two_body

  !==============================================================================================

  !==============================================================================================

  subroutine d_scaling_func_two_body(distance_dummy, d_scaled_distance_dummy)
    !---------------------------------------------------------------------------
    ! Description : d_scaling_func_two_body. This subroutine calculates the derivative of scaled
    !               distances of a rank-3 array of distances. The scaling form is specified
    !               by isc_bf and the scaling factor is specified by scalek_bf. Note that isc_bf and scalek_bf
    !               must be called with object_needed in any routine calling this routine such that they are valid
    !
    ! Created     : F. Petruzielo, 3 Feb 2009
    !---------------------------------------------------------------------------


    implicit none

    ! input
    real(dp), dimension(:,:,:), intent(in) :: distance_dummy

    ! output
    real(dp), dimension(:,:,:), intent(out) :: d_scaled_distance_dummy
    character(len=max_string_len_rout), save :: lhere = 'd_scaling_func_two_body'

    select case(isc_bf)
    case(4) ! d  [r / (1 + r * k)]
            !====================
            !     d r
       d_scaled_distance_dummy = 1.d0 / (1.d0 + scalek_bf * distance_dummy)**2
    case default
       call die(lhere,'Only isc = 4 is implemented.')
    end select

  end subroutine d_scaling_func_two_body

  !==============================================================================================

  !==============================================================================================

  subroutine d2_scaling_func_two_body(distance_dummy, d2_scaled_distance_dummy)
    !---------------------------------------------------------------------------
    ! Description : d2_scaling_func_two_body. This subroutine calculates the second derivative of scaled
    !               distances of a rank-3 array of distances. The scaling form is specified
    !               by isc_bf and the scaling factor is specified by scalek_bf. Note that isc_bf and scalek_bf
    !               must be called with object_needed in any routine calling this routine such that they are valid
    !
    ! Created     : F. Petruzielo, 3 Feb 2009
    !---------------------------------------------------------------------------

    implicit none

    ! input
    real(dp), dimension(:,:,:), intent(in) :: distance_dummy

    ! output
    real(dp), dimension(:,:,:), intent(out) :: d2_scaled_distance_dummy
    character(len=max_string_len_rout), save :: lhere = 'd2_scaling_func_two_body'

    select case(isc_bf)
    case(4) ! d^2  [r / (1 + r * k)]
            !====================
            !     d r^2
       d2_scaled_distance_dummy = -2.d0 * scalek_bf / (1.d0 + scalek_bf * distance_dummy)**3
    case default
       call die(lhere,'Only isc = 4 is implemented.')
    end select

  end subroutine d2_scaling_func_two_body

  !==============================================================================================

  !==============================================================================================

  subroutine scaled_dist_een_ee_wlk_bld

    !---------------------------------------------------------------------------
    ! Description : build scaled_dist_een_ee_wlk(:,:,:). Scaled distances between
    !               electrons calculated with 3-body scaling function
    !               This is a rank three array.
    !               The first index is the electron.
    !               The second index is the electron
    !               The third index is the walker
    !
    ! Created     : F. Petruzielo, 23 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('scaled_dist_een_ee_wlk')

       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed ('dist_ee_wlk')
       !needed in a subroutine that is called so make sure scalek_bf, isc_bf are up to date
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('scaled_dist_een_ee_wlk', scaled_dist_een_ee_wlk, nelec, nelec, nwalk)
    scaled_dist_een_ee_wlk(:,:,:) = 0.d0

    !calculate scaled distances
    call scaling_func_three_body(dist_ee_wlk,scaled_dist_een_ee_wlk)

  end subroutine scaled_dist_een_ee_wlk_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_scaled_dist_een_ee_wlk_d_r_bld

    !---------------------------------------------------------------------------
    ! Description : build d_scaled_dist_een_ee_wlk_d_r(:,:,:). This is a rank three array.
    !               The first index is an electron.
    !               The second index is another electron.
    !               The third index is the walker
    !               This returns the derivative of the 3-body scaled distances between electrons
    !               for each walker.
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('d_scaled_dist_een_ee_wlk_d_r')

       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed ('dist_ee_wlk')
       !needed in a subroutine that is called so make sure scalek_bf, isc_bf are up to date
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    ! allocation
    call object_alloc('d_scaled_dist_een_ee_wlk_d_r', d_scaled_dist_een_ee_wlk_d_r, nelec, nelec, nwalk)
    d_scaled_dist_een_ee_wlk_d_r(:,:,:) = 0.d0

    !calculate derivatives of scaled distances
    call d_scaling_func_three_body(dist_ee_wlk, d_scaled_dist_een_ee_wlk_d_r)

  end subroutine d_scaled_dist_een_ee_wlk_d_r_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d2_scaled_dist_een_ee_wlk_d_r_2_bld

    !---------------------------------------------------------------------------
    ! Description : build d2_scaled_dist_een_ee_wlk_d_r_2(:,:,:). This is a rank three array.
    !               The first index is an electron.
    !               The second index is another electron.
    !               The third index is the walker
    !               This returns the 2nd derivative of the 3-body scaled distances between electrons
    !               for each walker.
    ! Created     : F. Petruzielo, 26 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('d2_scaled_dist_een_ee_wlk_d_r_2')

       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed ('dist_ee_wlk')
       !needed in a subroutine that is called so make sure scalek_bf, isc_bf are up to date
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    ! allocation
    call object_alloc('d2_scaled_dist_een_ee_wlk_d_r_2', d2_scaled_dist_een_ee_wlk_d_r_2, nelec, nelec, nwalk)
    d2_scaled_dist_een_ee_wlk_d_r_2(:,:,:) = 0.d0

    !calculate 2nd derivatives of scaled distances
    call d2_scaling_func_three_body(dist_ee_wlk, d2_scaled_dist_een_ee_wlk_d_r_2)

  end subroutine d2_scaled_dist_een_ee_wlk_d_r_2_bld

  !==============================================================================================

  !==============================================================================================

  subroutine scaled_dist_een_en_wlk_bld

    !---------------------------------------------------------------------------
    ! Description : build scaled_dist_een_en_wlk(:,:,:). This is a rank three array.
    !               Calculate the scaled distances between electron and nucleus with the three body scaling functions.
    !               The first index is the electron.
    !               The second index is the electron
    !               The third index is the walker
    !
    ! Created     : F. Petruzielo, 23 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, nwalker, ncent
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('scaled_dist_een_en_wlk')

       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed('nwalk')
       call object_needed ('dist_en_wlk')
       !needed in a subroutine that is called so make sure scalek_bf, isc_bf are up to date
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('scaled_dist_een_en_wlk', scaled_dist_een_en_wlk, nelec, ncent, nwalk)
    scaled_dist_een_en_wlk(:,:,:) = 0.d0

    !calculate scaled distances
    call scaling_func_three_body(dist_en_wlk,scaled_dist_een_en_wlk)

  end subroutine scaled_dist_een_en_wlk_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_scaled_dist_een_en_wlk_d_r_bld

    !---------------------------------------------------------------------------
    ! Description : build d_scaled_dist_een_en_wlk_d_r(:,:,:). This is a rank three array.
    !               The first index is an electron.
    !               The second index is a center.
    !               The third index is the walker
    !               This returns the derivative of the 3-body scaled distances between electrons and nuclei
    !               for each walker.
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('d_scaled_dist_een_en_wlk_d_r')

       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed ('dist_en_wlk')
       !needed in a subroutine that is called so make sure scalek_bf, isc_bf are up to date
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('d_scaled_dist_een_en_wlk_d_r', d_scaled_dist_een_en_wlk_d_r, nelec, ncent, nwalk)
    d_scaled_dist_een_en_wlk_d_r(:,:,:) = 0.d0

    !calculate derivatives of scaled distances
    call d_scaling_func_three_body(dist_en_wlk, d_scaled_dist_een_en_wlk_d_r)

  end subroutine d_scaled_dist_een_en_wlk_d_r_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d2_scaled_dist_een_en_wlk_d_r_2_bld

    !---------------------------------------------------------------------------
    ! Description : build d2_scaled_dist_een_en_wlk_d_r_2(:,:,:). This is a rank three array.
    !               The first index is an electron.
    !               The second index is a center.
    !               The third index is the walker
    !               This returns the 2nd derivative of the 3-body scaled distances between electrons and nuclei
    !               for each walker.
    ! Created     : F. Petruzielo, 26 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('d2_scaled_dist_een_en_wlk_d_r_2')

       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed ('dist_en_wlk')
       !needed in a subroutine that is called so make sure scalek_bf, isc_bf are up to date
       call object_needed('scalek_bf')
       call object_needed('isc_bf')

       return
    endif

    call object_alloc('d2_scaled_dist_een_en_wlk_d_r_2', d2_scaled_dist_een_en_wlk_d_r_2, nelec, ncent, nwalk)
    d2_scaled_dist_een_en_wlk_d_r_2(:,:,:) = 0.d0

    !calculate 2nd derivatives of scaled distances
    call d2_scaling_func_three_body(dist_en_wlk, d2_scaled_dist_een_en_wlk_d_r_2)

  end subroutine d2_scaled_dist_een_en_wlk_d_r_2_bld

  !==============================================================================================

  !==============================================================================================

  subroutine scaling_func_three_body(distance_dummy, scaled_distance_dummy)
    !---------------------------------------------------------------------------
    ! Description : scaling_func_three_body. This subroutine calculates the scaled
    !               distances of rank-3 array of distances. The scaling form is specified
    !               by isc_bf and the scaling factor is specified by scalek_bf. Note that isc_bf and scalek_bf
    !               must be called with object_needed in any routine calling this routine such that they are valid
    !
    ! Created     : F. Petruzielo, 23 Jul 2008
    !---------------------------------------------------------------------------


    implicit none

    ! output
    real(dp), dimension(:,:,:), intent(out) :: scaled_distance_dummy

    ! input
    real(dp), dimension(:,:,:), intent(in) :: distance_dummy
    character(len=max_string_len_rout), save :: lhere = 'scaling_func_three_body'

    select case(isc_bf)
    case(4) ! 1 / (1 + r * k)
       scaled_distance_dummy = 1.d0 / (1.d0 + scalek_bf * distance_dummy)
    case default
       call die(lhere,'Only isc = 4 is implemented.')
    end select

  end subroutine scaling_func_three_body

  !==============================================================================================

  !==============================================================================================

  subroutine d_scaling_func_three_body(distance_dummy, d_scaled_distance_dummy)
    !---------------------------------------------------------------------------
    ! Description : d_scaling_func_three_body. This subroutine calculates the derivative of scaled
    !               distances (three-body scaling) of a rank-3 array of distances. The scaling form is specified
    !               by isc_bf and the scaling factor is specified by scalek_bf. Note that isc_bf and scalek_bf
    !               must be called with object_needed in any routine calling this routine such that they are valid
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------


    implicit none

    ! input
    real(dp), dimension(:,:,:), intent(in) :: distance_dummy

    ! output
    real(dp), dimension(:,:,:), intent(out) :: d_scaled_distance_dummy
    character(len=max_string_len_rout), save :: lhere = 'd_scaling_func_three_body'

    select case(isc_bf)
    case(4) ! d  [1 / (1 + r * k)]
            !====================
            !     d r
       d_scaled_distance_dummy = - scalek_bf / (1.d0 + scalek_bf * distance_dummy)**2
    case default
       call die(lhere,'Only isc = 4 is implemented.')
    end select

  end subroutine d_scaling_func_three_body

  !==============================================================================================

  !==============================================================================================

  subroutine d2_scaling_func_three_body(distance_dummy, d2_scaled_distance_dummy)
    !---------------------------------------------------------------------------
    ! Description : d2_scaling_func_three_body. This subroutine calculates the 2nd derivative of scaled
    !               distances (three-body scaling) of a rank-3 array of distances. The scaling form is specified
    !               by isc_bf and the scaling factor is specified by scalek_bf. Note that isc_bf and scalek_bf
    !               must be called with object_needed in any routine calling this routine such that they are valid
    !
    ! Created     : F. Petruzielo, 26 Mar 2009
    !---------------------------------------------------------------------------

    implicit none

    ! input
    real(dp), dimension(:,:,:), intent(in) :: distance_dummy

    ! output
    real(dp), dimension(:,:,:), intent(out) :: d2_scaled_distance_dummy
    character(len=max_string_len_rout), save :: lhere = 'd2_scaling_func_three_body'

    select case(isc_bf)
    case(4) ! d  [1 / (1 + r * k)]
            !====================
            !     d r
       d2_scaled_distance_dummy = 2 * scalek_bf**2 / (1.d0 + scalek_bf * distance_dummy)**3
    case default
       call die(lhere,'Only isc = 4 is implemented.')
    end select

  end subroutine d2_scaling_func_three_body

  !==============================================================================================

  !==============================================================================================

  subroutine smooth_cutoff_g_bld

    !---------------------------------------------------------------------------
    ! Description : build smooth_cutoff_g. .
    !               In all electron calculations, to impose the electron nuclear cusp, it is necessary to smoothly cutoff
    !               the backflow functions at the nucleus.
    !               The first index has dimension equal to the number of electrons.
    !               The second index has dimensions equal to the number of nuclei.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 20 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec
    include 'commons.h'

    !local
    integer  :: elec_i, cent_i, walk_i

    ! header
    if (header_exe) then

       call object_create('smooth_cutoff_g')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('threshold_l')
       call object_needed ('dist_en_wlk')

       return
    endif

    call object_alloc('smooth_cutoff_g', smooth_cutoff_g, nelec, ncent, nwalk)
    smooth_cutoff_g(:,:,:) = 0.d0

    do walk_i = 1, nwalk
       do cent_i = 1, ncent
          do elec_i = 1, nelec
             if( dist_en_wlk(elec_i, cent_i, walk_i) .ge. threshold_l(iwctype(cent_i))) then
                smooth_cutoff_g(elec_i, cent_i, walk_i) = 1.d0
             else
                smooth_cutoff_g(elec_i, cent_i, walk_i) = (dist_en_wlk(elec_i, cent_i, walk_i) / threshold_l(iwctype(cent_i)))**2 * (6.d0 - 8.d0 &
                     & * (dist_en_wlk(elec_i, cent_i, walk_i) / threshold_l(iwctype(cent_i))) + 3.d0 * (dist_en_wlk(elec_i, cent_i, walk_i) / threshold_l(iwctype(cent_i)))**2)
             end if
          enddo
       enddo
    enddo

  end subroutine smooth_cutoff_g_bld


  !==============================================================================================

  !==============================================================================================

  subroutine d_smooth_cutoff_g_d_r_bld

    !---------------------------------------------------------------------------
    ! Description : build d_smooth_cutoff_g_d_r. This is a rank three array.
    !               In all electron calculations, to impose the electron nuclear cusp, it is necessary to smoothly cutoff
    !               the backflow functions at the nucleus. This is the derivative of that function.
    !               The first index has dimension equal to the number of electrons, nelec.
    !               The second index has dimensions equal to the number of nuclei, ncent.
    !               The third index has dimension equal to the number of walkers, nwalk.
    !
    ! Created     : F. Petruzielo, 4 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec
    include 'commons.h'

    !local
    integer  :: elec_i, cent_i, walk_i

    ! header
    if (header_exe) then

       call object_create('d_smooth_cutoff_g_d_r')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('threshold_l')
       call object_needed ('dist_en_wlk')

       return
    endif

    ! allocation
    call object_alloc('d_smooth_cutoff_g_d_r', d_smooth_cutoff_g_d_r, nelec, ncent, nwalk)
    d_smooth_cutoff_g_d_r(:,:,:) = 0.d0

    do walk_i = 1, nwalk
       do cent_i = 1, ncent
          do elec_i = 1, nelec
             if( dist_en_wlk(elec_i, cent_i, walk_i) .ge. threshold_l(iwctype(cent_i))) then
                d_smooth_cutoff_g_d_r(elec_i, cent_i, walk_i) = 0.d0
             else
                d_smooth_cutoff_g_d_r(elec_i, cent_i, walk_i) = dist_en_wlk(elec_i, cent_i, walk_i) / ( threshold_l(iwctype(cent_i))**2 ) * ( 12.d0 - 24.d0 &
                     & * (dist_en_wlk(elec_i, cent_i, walk_i) / threshold_l(iwctype(cent_i))) + 12.d0 * (dist_en_wlk(elec_i, cent_i, walk_i) / threshold_l(iwctype(cent_i)))**2 )
             end if
          enddo
       enddo
    enddo

  end subroutine d_smooth_cutoff_g_d_r_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d2_smooth_cutoff_g_d_r_2_bld

    !---------------------------------------------------------------------------
    ! Description : build d2_smooth_cutoff_g_d_r_2. This is a rank three array.
    !               In all electron calculations, to impose the electron nuclear cusp, it is necessary to smoothly cutoff
    !               the backflow functions at the nucleus. This is the second derivative of that function.
    !               The first index has dimension equal to the number of electrons, nelec.
    !               The second index has dimensions equal to the number of nuclei, ncent.
    !               The third index has dimension equal to the number of walkers, nwalk.
    !
    ! Created     : F. Petruzielo, 17 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec
    include 'commons.h'

    !local
    integer  :: elec_i, cent_i, walk_i

    ! header
    if (header_exe) then

       call object_create('d2_smooth_cutoff_g_d_r_2')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('threshold_l')
       call object_needed ('dist_en_wlk')

       return
    endif

    ! allocation
    call object_alloc('d2_smooth_cutoff_g_d_r_2', d2_smooth_cutoff_g_d_r_2, nelec, ncent, nwalk)
    d2_smooth_cutoff_g_d_r_2(:,:,:) = 0.d0

    do walk_i = 1, nwalk
       do cent_i = 1, ncent
          do elec_i = 1, nelec
             if( dist_en_wlk(elec_i, cent_i, walk_i) .ge. threshold_l(iwctype(cent_i))) then
                d2_smooth_cutoff_g_d_r_2(elec_i, cent_i, walk_i) = 0.d0
             else
                d2_smooth_cutoff_g_d_r_2(elec_i, cent_i, walk_i) = ( 12.d0 - 48.d0 * dist_en_wlk(elec_i, cent_i, walk_i) / threshold_l(iwctype(cent_i)) + 36.d0 * (dist_en_wlk(elec_i, cent_i, walk_i) / threshold_l(iwctype(cent_i)))**2 ) /  threshold_l(iwctype(cent_i))**2
             end if
          enddo
       enddo
    enddo

  end subroutine d2_smooth_cutoff_g_d_r_2_bld


  !==============================================================================================

  !===================================================

  subroutine row_reduce(matrix_coeff)
    !---------------------------------------------------------------------------
    ! Description : This subroutine row reduces a matrix.
    !               The coefficients of the matrix are stored in array matrix_coeff
    !               Note the first index should be the row and the second index should be the column
    !
    ! Created     : F. Petruzielo, 10 August 2008
    !---------------------------------------------------------------------------

    implicit none

    !input
    real(dp),dimension(:,:),intent(inout) :: matrix_coeff

    !Local variables
    character(len=max_string_len_rout), save :: lhere = 'row_reduce'
    real(dp), allocatable :: temp_array(:)
    integer,dimension(1):: index_max_in_column
    integer :: index_max_in_column_temp
    integer :: row_i, column_i, pivot_row_i
    integer :: nb_row, nb_column, alloc_stat

    !get dimensions of matrix_coeff
    nb_column = size(matrix_coeff,2)
    nb_row = size(matrix_coeff,1)

    !allocate temporary array for switching rows
    allocate (temp_array(nb_column), stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       !allocation failed
       call require(lhere,'allocation successful',alloc_stat .eq. 0)
    endif
    !allocation successful

    !initialize pivot row
    pivot_row_i=1

    !loop over columns
    do column_i = 1, nb_column
       !Find row below pivot_row_i (inclusive) with largest value of |matrix_coeff(:,column_i)|
       !maxloc returns an array
       index_max_in_column = maxloc(abs(matrix_coeff(pivot_row_i:nb_row,column_i)))
       !correct position for cutting off rows above pivot_row_i
       index_max_in_column_temp = index_max_in_column(1) + pivot_row_i - 1

       !Interchange pivot_row_i with the row with max entry in column_i, if necessary
        if (index_max_in_column_temp .ne. pivot_row_i) then
           temp_array = matrix_coeff(pivot_row_i,:)
           matrix_coeff(pivot_row_i,:) = matrix_coeff(index_max_in_column_temp,:)
           matrix_coeff(index_max_in_column_temp,:) = temp_array
        end if

        if (abs(matrix_coeff(pivot_row_i,column_i)) .gt. 1.0d-10) then !at least one entry of the column is nonzero
           !Subtract multiples of pivot_row_i from subsequent rows to zero all subsequent coefficients of column_i (done more efficiently than simply subtracting)
           do row_i = pivot_row_i + 1, nb_row
              !note that in row_i the columns before column_i are already zero
              matrix_coeff(row_i,column_i + 1:nb_column) = matrix_coeff(row_i,column_i + 1:nb_column) - matrix_coeff(row_i, column_i) / matrix_coeff(pivot_row_i, column_i) * matrix_coeff(pivot_row_i,column_i + 1:nb_column)
              !zero out the rest of column_i
              matrix_coeff(row_i,column_i) = 0.d0
           end do

           !Subtract multiples of pivot_row_i from previous rows to zero all previous coefficients of column_i
           do row_i = 1, pivot_row_i - 1
              !note that in row_i the columns before column_i are already zero (or pivots)
              matrix_coeff(row_i,column_i + 1:nb_column) = matrix_coeff(row_i,column_i + 1:nb_column) - matrix_coeff(row_i, column_i) / matrix_coeff(pivot_row_i, column_i) * matrix_coeff(pivot_row_i,column_i + 1:nb_column)
              !zero out the rest of column_i
              matrix_coeff(row_i,column_i) = 0.d0
           enddo

           !normalize
           matrix_coeff(pivot_row_i,column_i+1:nb_column) = matrix_coeff(pivot_row_i,column_i+1:nb_column) / matrix_coeff(pivot_row_i,column_i)
           matrix_coeff(pivot_row_i,column_i)=1.d0
           !found a pivot in this row
           pivot_row_i=pivot_row_i +1

        else
           !every entry of column is zero
           cycle
        endif

        if (pivot_row_i >  nb_row) exit !found all the pivots
     end do

    deallocate (temp_array, stat = alloc_stat)
    if (alloc_stat .ne. 0) then
       !deallocation failed
       call require(lhere,'deallocation successful',alloc_stat .eq. 0)
    endif
    !deallocation successful

  end subroutine row_reduce

  !==============================================================================================

  !======================================================================================

  subroutine xi_bld
    !---------------------------------------------------------------------------
    ! Description : build object xi(:,:,:). This is a rank three array.
    !               It is the total backflow.
    !               The first index has dimension equal to the number of spatial dimensions, ndim.
    !               The second index has dimensions equal to the number of electrons, nelec.
    !               The third index has dimension equal to the number of walkers, nwalk.
    !
    ! Created     : F. Petruzielo, 10 Apr 2009
    !---------------------------------------------------------------------------
    implicit none

    !need ndim, nelec, nwalk
    include 'commons.h'

    !local
    character(len=max_string_len_rout), save :: lhere = 'xi_bld'

    ! header
    if (header_exe) then

       call object_create('xi')

       call object_needed('ndim')
       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed('en_bf')
       call object_needed('ee_bf')
       call object_needed('een_theta_bf')
       call object_needed('een_phi_bf')

       return
    endif

    call object_alloc('xi', xi, ndim, nelec, nwalk)
    xi(:,:,:) = 0.d0

    if (en_bf) then
       call object_provide_in_node(lhere, 'xi_en')
       xi(:, :, :) = xi(:, :, :) + xi_en(:, :, :)
    end if

    if (ee_bf) then
       call object_provide_in_node(lhere, 'xi_ee')
       xi(:, :, :) = xi(:, :, :) + xi_ee(:, :, :)
    end if

    if (een_theta_bf) then
       call object_provide_in_node(lhere, 'xi_een_theta')
       xi(:, :, :) = xi(:, :, :) + xi_een_theta(:, :, :)
    end if

    if (een_phi_bf) then
       call object_provide_in_node(lhere, 'xi_een_phi')
       xi(:, :, :) = xi(:, :, :) + xi_een_phi(:, :, :)
    end if

  end subroutine xi_bld

  !======================================================================================

  !==============================================================================================

  subroutine d_xi_j_beta_d_r_i_alpha_bld

    !---------------------------------------------------------------------------
    ! Description : build d_xi_j_beta_d_r_i_alpha. This is a rank five array.
    !               Derivative of xi with respect to electron position.
    !               It is necessary for the construction of the energy, drift velocity.
    !               The first and second indices are the number of spatial dimensions.
    !               The first index is for the components of xi (beta)
    !               The second index is for the components of the gradient (alpha)
    !               The third and fourth indices have dimensions equal to the number of electrons.
    !               The third index is the label of the backflow transformation of each electon
    !               The fourth index is for the electron that the gradient is taken with respect to
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 10 Apr 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim, ncent
    include 'commons.h'

    !local
    integer  :: elec_j

    !local
    character(len=max_string_len_rout), save :: lhere = 'd_xi_j_beta_d_r_i_alpha_bld'

    ! header
    if (header_exe) then

       call object_create('d_xi_j_beta_d_r_i_alpha')

       call object_needed('ndim')
       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed('en_bf')
       call object_needed('ee_bf')
       call object_needed('een_theta_bf')
       call object_needed('een_phi_bf')

       return
    endif

    call object_alloc('d_xi_j_beta_d_r_i_alpha', d_xi_j_beta_d_r_i_alpha, ndim, ndim, nelec, nelec, nwalk)
    d_xi_j_beta_d_r_i_alpha(:, :, :, :, :) = 0.d0

    if (en_bf) then
       call object_provide_in_node(lhere, 'd_xi_en_j_beta_d_r_i_alpha')
       do elec_j = 1, nelec
          d_xi_j_beta_d_r_i_alpha(:, :, elec_j, elec_j, :) = d_xi_j_beta_d_r_i_alpha(:, :, elec_j, elec_j, :) + d_xi_en_j_beta_d_r_i_alpha(:, :, elec_j, :)
       end do
    end if

    if (ee_bf) then
       call object_provide_in_node(lhere, 'd_xi_ee_j_beta_d_r_i_alpha')
       d_xi_j_beta_d_r_i_alpha(:, :, :, :, :) = d_xi_j_beta_d_r_i_alpha(:, :, :, :, :) + d_xi_ee_j_beta_d_r_i_alpha(:, :, :, :, :)
    end if

    if (een_theta_bf) then
       call object_provide_in_node(lhere, 'd_xi_een_theta_j_beta_d_r_i_alpha')
       d_xi_j_beta_d_r_i_alpha(:, :, :, :, :) = d_xi_j_beta_d_r_i_alpha(:, :, :, :, :) + d_xi_een_theta_j_beta_d_r_i_alpha(:, :, :, :, :)
    end if

    if (een_phi_bf) then
       call object_provide_in_node(lhere, 'd_xi_een_phi_j_beta_d_r_i_alpha')
       d_xi_j_beta_d_r_i_alpha(:, :, :, :, :) = d_xi_j_beta_d_r_i_alpha(:, :, :, :, :) + d_xi_een_phi_j_beta_d_r_i_alpha(:, :, :, :, :)
    end if

  end subroutine d_xi_j_beta_d_r_i_alpha_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_xi_j_beta_d_r_i_alpha_num_bld
    !---------------------------------------------------------------------------
    ! Description : build d_xi_j_beta_d_r_i_alpha__num. This is a rank five array.
    !               Numerical first derivative of the xi:
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               The first and second indices are the number of spatial dimensions.
    !               The first index is for the components of xi (beta)
    !               The second index is for the components of the gradient (alpha)
    !               The third and fourth indices have dimensions equal to the number of electrons.
    !               The third index is the label of the backflow transformation of each electon
    !               The fourth index is for the electron that the gradient is taken with respect to
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 10 Apr 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, dim_i, dim_j
    real(dp), allocatable :: xi_temp_plus(:,:,:), xi_temp_minus(:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'd_xi_j_beta_d_r_i_alpha_num_bld'

    ! header
    if (header_exe) then

       call object_create('d_xi_j_beta_d_r_i_alpha_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('xi')
       call object_needed('coord_elec_wlk')

       return
    endif

    ! allocation
    call object_alloc('d_xi_j_beta_d_r_i_alpha_num', d_xi_j_beta_d_r_i_alpha_num, ndim, ndim, nelec, nelec, nwalk)
    d_xi_j_beta_d_r_i_alpha_num(:,:,:,:,:) = 0.d0

    !allocate temporary objects
    allocate(xi_temp_plus(ndim, nelec, nwalk))
    allocate(xi_temp_minus(ndim, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

    do elec_j = 1, nelec
       do elec_i = 1, nelec
          do dim_j = 1, ndim  !beta
             do dim_i = 1, ndim  !alph
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) + epsilon
                call object_modified('coord_elec_wlk')
                !get updated xi
                call object_provide_in_node(lhere, 'xi')
                xi_temp_plus = xi
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) - epsilon
                call object_modified('coord_elec_wlk')
                !get updated xi
                call object_provide_in_node(lhere, 'xi')
                xi_temp_minus = xi
                !calculate derivative
                d_xi_j_beta_d_r_i_alpha_num(dim_j, dim_i, elec_j, elec_i, :) = (xi_temp_plus(dim_j, elec_j, :) - xi_temp_minus(dim_j, elec_j, :) ) / (2.d0 * epsilon)
             end do
          end do
       end do
    end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'xi')

  end subroutine d_xi_j_beta_d_r_i_alpha_num_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_j_beta_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_j_beta. This is a rank four array.
    !               Laplacian of the backflow function
    !               with respect to electron position.
    !               It is necessary for the construction of the energy.
    !               The first index is for the components of xi (beta)
    !               The second index is the label of the backflow transformed electon
    !               The third index is for the electron that the laplacian is taken with respect to
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 10 Apr 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j
    character(len=max_string_len_rout), save :: lhere = 'lap_i_xi_j_beta_bld'

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_j_beta')

       call object_needed('ndim')
       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed('en_bf')
       call object_needed('ee_bf')
       call object_needed('een_theta_bf')
       call object_needed('een_phi_bf')

       return
    endif

    call object_alloc('lap_i_xi_j_beta', lap_i_xi_j_beta, ndim, nelec, nelec, nwalk)
    lap_i_xi_j_beta(:,:,:,:) = 0.d0

    if (en_bf) then
       call object_provide_in_node(lhere, 'lap_i_xi_en_j_beta')
       do elec_j = 1, nelec
          lap_i_xi_j_beta(:, elec_j, elec_j, :) = lap_i_xi_j_beta(:, elec_j, elec_j, :) + lap_i_xi_en_j_beta(:, elec_j, :)
       end do
    end if

    if (ee_bf) then
       call object_provide_in_node(lhere, 'lap_i_xi_ee_j_beta')
       lap_i_xi_j_beta(:, :, :, :) = lap_i_xi_j_beta(:, :, :, :) + lap_i_xi_ee_j_beta(:, :, :, :)
    end if

    if (een_theta_bf) then
       call object_provide_in_node(lhere, 'lap_i_xi_een_theta_j_beta')
       lap_i_xi_j_beta(:, :, :, :) = lap_i_xi_j_beta(:, :, :, :) + lap_i_xi_een_theta_j_beta(:, :, :, :)
    end if

    if (een_phi_bf) then
       call object_provide_in_node(lhere, 'lap_i_xi_een_phi_j_beta')
       lap_i_xi_j_beta(:, :, :, :) = lap_i_xi_j_beta(:, :, :, :) + lap_i_xi_een_phi_j_beta(:, :, :, :)
    end if

  end subroutine lap_i_xi_j_beta_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_j_beta_num_bld
    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_j_beta_num. This is a rank four array.
    !               Numerical laplacian of xi:
    !
    !               (nabla_i)^2 xi_j^beta
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               The first index is for the components of xi (beta)
    !               The second index is the label of the backflow transformation of each electon
    !               The third index is for the electron that the laplacian is taken with respect to
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 10 Apr 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, dim_i, dim_j
    real(dp), allocatable :: d_xi_j_beta_d_r_i_alpha_plus(:,:,:,:,:), d_xi_j_beta_d_r_i_alpha_minus(:,:,:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'lap_i_xi_j_beta_num_bld'

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_j_beta_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('coord_elec_wlk')
       call object_needed('d_xi_j_beta_d_r_i_alpha')

       return
    endif

    ! allocation
    call object_alloc('lap_i_xi_j_beta_num', lap_i_xi_j_beta_num, ndim, nelec, nelec, nwalk)
    lap_i_xi_j_beta_num(:,:,:,:) = 0.d0

    !allocate temporary objects
    allocate(d_xi_j_beta_d_r_i_alpha_plus(ndim, ndim, nelec, nelec, nwalk))
    allocate(d_xi_j_beta_d_r_i_alpha_minus(ndim, ndim, nelec, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

     do elec_i = 1, nelec
        do elec_j = 1, nelec
           do dim_j = 1, ndim  !beta
              do dim_i = 1, ndim  !alpha
                 coord_elec_wlk = coord_elec_wlk_temp
                 !shift coordinates in a particular dimension, x,y,z, etc
                 coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) + epsilon
                 call object_modified('coord_elec_wlk')
                 !get updated d_xi_j_beta_d_r_i_alpha
                 call object_provide_in_node(lhere, 'd_xi_j_beta_d_r_i_alpha')
                 d_xi_j_beta_d_r_i_alpha_plus = d_xi_j_beta_d_r_i_alpha
                 coord_elec_wlk = coord_elec_wlk_temp
                 !shift coordinates in a particular dimension, x,y,z, etc
                 coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) - epsilon
                 call object_modified('coord_elec_wlk')
                 !get updated d_xi_j_beta_d_r_i_alpha
                 call object_provide_in_node(lhere, 'd_xi_j_beta_d_r_i_alpha')
                 d_xi_j_beta_d_r_i_alpha_minus = d_xi_j_beta_d_r_i_alpha
                 !calculate derivative
                 lap_i_xi_j_beta_num(dim_j, elec_j, elec_i, :) = lap_i_xi_j_beta_num(dim_j, elec_j, elec_i, :) + (d_xi_j_beta_d_r_i_alpha_plus(dim_j, dim_i, elec_j, elec_i, :) - d_xi_j_beta_d_r_i_alpha_minus(dim_j, dim_i, elec_j, elec_i, :) ) / (2.d0 * epsilon)
              end do
           end do
        end do
     end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'd_xi_j_beta_d_r_i_alpha')

  end subroutine lap_i_xi_j_beta_num_bld

  !==============================================================================================

  !======================================================================================

  subroutine xi_en_bld
    !---------------------------------------------------------------------------
    ! Description : build object xi_en(:,:,:). This is a rank three array.
    !               It is the electron-nuclear backflow.
    !               The first index has dimension equal to the number of spatial dimensions, ndim.
    !               The second index has dimensions equal to the number of electrons, nelec.
    !               The third index has dimension equal to the number of walkers, nwalk.
    !
    ! Created     : F. Petruzielo, 21 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need ndim, nelec, nwalk
    include 'commons.h'

    !local
    integer  :: dim_i

    ! header
    if (header_exe) then

       call object_create('xi_en')

       call object_needed('ndim')
       call object_needed('nelec')
       call object_needed('nwalk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('mu')

       return
    endif

    call object_alloc('xi_en', xi_en, ndim, nelec, nwalk)
    xi_en(:,:,:) = 0.d0

    do dim_i = 1, ndim
       xi_en(dim_i, :, :) = sum(mu(:, :, :) * vec_en_xyz_wlk(dim_i, :, :, :), 2)
    enddo

  end subroutine xi_en_bld

  !======================================================================================

  !======================================================================================

  subroutine mu_bld
    !---------------------------------------------------------------------------
    ! Description : build object mu(:,:,:). This is a rank three array.
    !               It is necessary for the construction of the electron-nuclear backflow terms
    !               The first index has dimension equal to the number of electrons, nelec.
    !               The second index has dimensions equal to the number of nuclei, ncent.
    !               The third index has dimension equal to the number of walkers, nwalk.
    !
    ! Created     : F. Petruzielo, 20 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: cent_i, order_p
    integer, parameter  :: up=1, down=2 !index for the spin
    character(len=max_string_len_rout), save :: lhere = 'mu_bld'

    ! header
    if (header_exe) then

       call object_create('mu')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('iwctype')
       call object_needed('order_mu_bf')
       call object_needed ('scaled_dist_en_wlk')
       call object_needed('d_param_mu')
       call object_needed('all_elec')

       return
    endif

    call object_alloc('mu', mu, nelec, ncent, nwalk)
    mu (:,:,:) = 0.d0

    do cent_i = 1, ncent
       !pade term of mu
       mu(1:nup, cent_i, :) = d_param_mu(1, iwctype(cent_i), up) * scaled_dist_en_wlk(1:nup, cent_i, :) / (1.d0 + d_param_mu(2, iwctype(cent_i), up) * scaled_dist_en_wlk(1:nup, cent_i, :))
       mu(nup+1:nelec, cent_i, :) = d_param_mu(1, iwctype(cent_i), down) * scaled_dist_en_wlk(nup+1:nelec, cent_i, :) / (1.d0 + d_param_mu(2, iwctype(cent_i), down) * scaled_dist_en_wlk(nup+1:nelec, cent_i, :))
       !power series terms of mu
       do order_p=2, order_mu_bf
          mu(1:nup, cent_i, :) = mu(1:nup, cent_i, :) + d_param_mu(order_p + 1, iwctype(cent_i), up) * scaled_dist_en_wlk(1:nup, cent_i, :)**order_p
          mu(nup+1:nelec, cent_i, :) = mu(nup+1:nelec, cent_i, :) + d_param_mu(order_p + 1, iwctype(cent_i), down) * scaled_dist_en_wlk(nup+1:nelec, cent_i, :)**order_p
       enddo
       if (all_elec) then
          !cutoff
          call object_provide_in_node(lhere, 'smooth_cutoff_g')
          mu(:,cent_i,:) = mu(:,cent_i,:) * product(smooth_cutoff_g(:,1:cent_i-1,:),2) * product(smooth_cutoff_g(:,cent_i+1:ncent,:),2)
       else
          !remove asymptotic value
          call object_provide_in_node(lhere, 'asymp_mu_pseudo')
          mu(1:nup,cent_i,:) = mu(1:nup,cent_i,:) - asymp_mu_pseudo(iwctype(cent_i), up)
          mu(nup+1:nelec,cent_i,:) = mu(nup+1:nelec,cent_i,:) - asymp_mu_pseudo(iwctype(cent_i), down)
       end if
    end do

  end subroutine mu_bld

  !==============================================================================================

  !==============================================================================================


  subroutine asymp_mu_pseudo_bld
    !---------------------------------------------------------------------------
    ! Description : build object asymp_mu_pseudo(:,:).
    !               This is a rank two array.
    !               The first index has dimension equal to number of center types, nctype.
    !               The second index has dimension equal to number of spins, spin_nb.
    !               Give the asymptotic value of the pade/power series expansion used in mu.
    !               This value is subtracted off from the pade/power series expansion such that mu goes to zero at infinity.
    !
    ! Created     : F. Petruzielo, 21 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_p

    ! header
    if (header_exe) then

       call object_create('asymp_mu_pseudo')

       call object_needed('nctype')
       call object_needed('order_mu_bf')
       call object_needed('d_param_mu')
       call object_needed('asymp_scaled_dist_two_body')

       return
    endif

    call object_alloc('asymp_mu_pseudo', asymp_mu_pseudo, nctype, spin_nb)
    asymp_mu_pseudo (:,:) = 0.d0

    !pade term of mu
    asymp_mu_pseudo(:, :) = d_param_mu(1, :, :) * asymp_scaled_dist_two_body / (1.d0 + d_param_mu(2, :, :) * asymp_scaled_dist_two_body)
    !power series terms of mu
    do order_p=2, order_mu_bf
       asymp_mu_pseudo(:, :) = asymp_mu_pseudo(:, :) + d_param_mu(order_p + 1, :, :) * asymp_scaled_dist_two_body**order_p
    enddo

  end subroutine asymp_mu_pseudo_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_param_mu_bld

    !---------------------------------------------------------------------------
    ! Description : build d_param_mu. This is a rank 3 array.
    !               List of d parameters for mu backflow function for each center type, and each spin
    !               The first index is which parameter.
    !               The second index is the type of center.
    !               The third index is the type of spin.
    !
    !               cusp condition constraints for all electron: mu(0)=0
    !               cusp condition constraints for pseudo: none
    !
    ! Created     : F. Petruzielo, 20 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_p, d_param_nb
    character(len=max_string_len_rout), save :: lhere = 'd_param_mu_bld'

    ! header
    if (header_exe) then

       call object_create('d_param_mu')

       call object_needed('all_elec')
       call object_needed('order_mu_bf')
       call object_needed('nctype')
       call object_needed('asymp_scaled_dist_two_body')
       call object_needed('mu_spin_dependent')
       call object_needed('read_d_param')
       call object_needed('read_d_param_mu_nb')

       return
    endif

    !need to have the option of using read_d_param_mu if it is read in
    if (read_d_param) then
       call object_provide_in_node(lhere, 'read_d_param_mu')
    end if

    call object_alloc('d_param_mu', d_param_mu, order_mu_bf + 1, nctype, spin_nb)
    d_param_mu = 0.d0

    !do we read in parameters?
    if (read_d_param) then
       if (mu_spin_dependent) then
          d_param_nb = (order_mu_bf + 1)  * nctype * spin_nb
          if (read_d_param_mu_nb == d_param_nb) then
             d_param_mu = reshape( read_d_param_mu, (/ order_mu_bf + 1, nctype, spin_nb /))
          else
             call require(lhere, 'number of parameters in read_d_param_mu equal to (order_mu_bf + 1) * nctype * spin_nb', .false.)
          end if
       else
          d_param_nb = (order_mu_bf + 1)  * nctype
          if (read_d_param_mu_nb == d_param_nb) then
             d_param_mu = reshape( read_d_param_mu, (/ order_mu_bf + 1, nctype, spin_nb /), pad=read_d_param_mu)
          else
             call require(lhere, 'number of parameters in read_d_param_mu equal to (order_mu_bf + 1) * nctype', .false.)
          end if
       end if
       !don't want to re-read
       read_d_param = .false.
    else
       call random_number(d_param_mu(:,:,:))
       if (.not. mu_spin_dependent) then
          d_param_mu(:,:,2) = d_param_mu(:,:,1)
       end if
    endif

    if (all_elec) then
       !impose electron nuclear cusp conditions assuming the two-body scaling functions evaluated at 0 are 0.
       d_param_mu(1,:,:) = 0.d0
       do order_p=2, order_mu_bf
          d_param_mu(1, :, :) = d_param_mu(1, :, :) - d_param_mu(order_p + 1, :, :) * asymp_scaled_dist_two_body**(order_p - 1)
       enddo
       d_param_mu(1, :, :) = d_param_mu(1, :, :) * ( 1.d0 +  d_param_mu(2, :, :) * asymp_scaled_dist_two_body)
    end if

  end subroutine d_param_mu_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_xi_en_j_beta_d_r_i_alpha_bld

    !---------------------------------------------------------------------------
    ! Description : build d_xi_en_j_beta_d_r_i_alpha. This is a rank four array.
    !               Derivative of the electron-nuclear backflow function with respect to electron position.
    !               It is necessary for the construction of the energy, drift velocity.
    !               The first and second indices are the number of spatial dimensions.
    !               The first index is for the components of xi_en (beta)
    !               The second index is for the components of the gradient (alpha)
    !               The third index has dimensions equal to the number of electrons.
    !               The third index is the label of the en backflow transformation of each electon
    !               and is also the label of the coordinate the gradient is taken wrt (because of delta_ij)
    !               The fourth index has dimension equal to the number of walkers.
    !               It is necessary for the construction of the energy, drift velocity.
    !
    ! Created     : F. Petruzielo, 3 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, nelec, ndim
    include 'commons.h'

    !local
    integer  :: dim_i, dim_j

    ! header
    if (header_exe) then

       call object_create('d_xi_en_j_beta_d_r_i_alpha')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('nelec')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('ndim')
       call object_needed('mu')
       call object_needed('d_mu_d_r_j_alpha')

       return
    endif

    call object_alloc('d_xi_en_j_beta_d_r_i_alpha', d_xi_en_j_beta_d_r_i_alpha, ndim, ndim, nelec, nwalk)
    d_xi_en_j_beta_d_r_i_alpha(:,:,:,:) = 0.d0


    !note the implicit delta function here since there is only one sum over electrons
    do dim_i = 1, ndim !alpha
       do dim_j = 1, ndim !beta
          d_xi_en_j_beta_d_r_i_alpha(dim_j, dim_i, :, :) = sum(d_mu_d_r_j_alpha(dim_i,:, :, :) * vec_en_xyz_wlk(dim_j, :, :, :), 2)
       end do
       d_xi_en_j_beta_d_r_i_alpha(dim_i, dim_i, :, :) = d_xi_en_j_beta_d_r_i_alpha(dim_i, dim_i, :, :) + sum(mu(:, :, :), 2)
    end do

  end subroutine d_xi_en_j_beta_d_r_i_alpha_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_xi_en_j_beta_d_r_i_alpha_num_bld
    !---------------------------------------------------------------------------
    ! Description : build d_xi_en_j_beta_d_r_i_alpha__num. This is a rank four array.
    !               Numerical first derivative of the xi_en:
    !
    !               d xi_en(r_jI)
    !               ----------
    !               d r_i_alpha
    !
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               Note alpha refers to cartesian coords. In 3-d, r_alpha takes on values x,y,z.
    !               The first index is the over spatial dimension.
    !               The second index is the over spatial dimension.
    !               The third index has dimensions equal to the number of electrons.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 4 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, dim_i, dim_j
    real(dp), allocatable :: xi_en_temp_plus(:,:,:), xi_en_temp_minus(:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'd_xi_en_j_beta_d_r_i_alpha_num_bld'

    ! header
    if (header_exe) then

       call object_create('d_xi_en_j_beta_d_r_i_alpha_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('xi_en')
       call object_needed('coord_elec_wlk')

       return
    endif

    call object_alloc('d_xi_en_j_beta_d_r_i_alpha_num', d_xi_en_j_beta_d_r_i_alpha_num, ndim, ndim, nelec, nwalk)
    d_xi_en_j_beta_d_r_i_alpha_num(:,:,:,:) = 0.d0

    allocate(xi_en_temp_plus(ndim, nelec, nwalk))
    allocate(xi_en_temp_minus(ndim, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

    do elec_i = 1, nelec
       do dim_i = 1, ndim  !alpha
          do dim_j = 1, ndim  !beta
             coord_elec_wlk = coord_elec_wlk_temp
             !shift coordinates in a particular dimension, x,y,z, etc
             coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) + epsilon
             call object_modified('coord_elec_wlk')
             !get updated xi_en
             call object_provide_in_node(lhere, 'xi_en')
             xi_en_temp_plus = xi_en
             coord_elec_wlk = coord_elec_wlk_temp
             !shift coordinates in a particular dimension, x,y,z, etc
             coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) - epsilon
             call object_modified('coord_elec_wlk')
             !get updated xi_en
             call object_provide_in_node(lhere, 'xi_en')
             xi_en_temp_minus = xi_en
             !calculate derivative
             d_xi_en_j_beta_d_r_i_alpha_num(dim_j, dim_i, elec_i, :) = (xi_en_temp_plus(dim_j, elec_i, :) - xi_en_temp_minus(dim_j, elec_i, :) ) / (2.d0*epsilon)
          end do
       end do
    end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'xi_en')

  end subroutine d_xi_en_j_beta_d_r_i_alpha_num_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_mu_d_r_j_alpha_bld

    !---------------------------------------------------------------------------
    ! Description : build d_mu_d_r_j_alpha. This is a rank four array.
    !               First derivative of the mu:
    !
    !               d mu(r_jI)
    !               ----------
    !               d r_j_alpha
    !
    !               Note alpha refers to cartesian coords. In 3-d, r_alpha takes on values x,y,z.
    !               The first index is the over spatial dimension.
    !               The second index has dimensions equal to the number of electrons.
    !               The third index has dimensions equal to the number of centers.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 4 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, nelec, ndim
    include 'commons.h'

    !local
    integer  :: cent_i, dim_i
    character(len=max_string_len_rout), save :: lhere = 'd_mu_d_r_j_alpha_bld'

    ! header
    if (header_exe) then

       call object_create('d_mu_d_r_j_alpha')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('nelec')
       call object_needed('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('d_mu_d_r_j_inuc')
       call object_needed('all_elec')

       return
    endif

    ! allocation
    call object_alloc('d_mu_d_r_j_alpha', d_mu_d_r_j_alpha, ndim, nelec, ncent, nwalk)
    d_mu_d_r_j_alpha(:,:,:,:) = 0.d0

    do dim_i = 1, ndim
       d_mu_d_r_j_alpha(dim_i, :, :, :) = d_mu_d_r_j_inuc(:, :, :) * vec_en_xyz_wlk(dim_i, :, :, :) / dist_en_wlk(:, :, :)
       if(all_elec) then
          call object_provide_in_node(lhere, 'mu')
          call object_provide_in_node(lhere, 'smooth_cutoff_g')
          call object_provide_in_node(lhere, 'd_smooth_cutoff_g_d_r')
          do cent_i = 1, ncent
             d_mu_d_r_j_alpha(dim_i, :, cent_i, :) = d_mu_d_r_j_alpha(dim_i, :, cent_i, :) * product(smooth_cutoff_g(:, 1:cent_i-1, :), 2) * product(smooth_cutoff_g(:, cent_i+1:ncent, :), 2) + mu(:,cent_i, :) * (sum(vec_en_xyz_wlk(dim_i, :, 1:cent_i-1, :)/ dist_en_wlk(:, 1:cent_i-1, :) / smooth_cutoff_g(:, 1:cent_i-1, :) * d_smooth_cutoff_g_d_r(:, 1:cent_i-1, :),2) + sum(vec_en_xyz_wlk(dim_i, :, cent_i+1:ncent, :)/ dist_en_wlk(:, cent_i+1:ncent, :) / smooth_cutoff_g(:, cent_i+1:ncent, :) * d_smooth_cutoff_g_d_r(:, cent_i+1:ncent, :),2))
          end do
       end if
    end do

  end subroutine d_mu_d_r_j_alpha_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_mu_d_r_j_inuc_bld

    !---------------------------------------------------------------------------
    ! Description : build d_mu_d_r_j_inuc. This is a rank three array.
    !
    !               d mu(r_jI)
    !               ----------    Without cutoff! (appears in all electron and pseudo case)
    !               d r_jI
    !
    !               The first index has dimensions equal to the number of electrons.
    !               The second index has dimensions equal to the number of centers.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 3 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, nelec, ndim
    include 'commons.h'

    !local
    integer  :: cent_i, order_p
    integer, parameter  :: up=1, down=2 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d_mu_d_r_j_inuc')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_mu_bf')
       call object_needed ('scaled_dist_en_wlk')
       call object_needed ('d_scaled_dist_en_wlk_d_r')
       call object_needed('d_param_mu')

       return
    endif

    call object_alloc('d_mu_d_r_j_inuc', d_mu_d_r_j_inuc, nelec, ncent, nwalk)
    d_mu_d_r_j_inuc(:,:,:) = 0.d0

    do cent_i = 1, ncent
       !pade terms
       d_mu_d_r_j_inuc(1:nup, cent_i, :) = d_param_mu(1, iwctype(cent_i), up) / (1.d0 + d_param_mu(2, iwctype(cent_i), up) * scaled_dist_en_wlk(1:nup, cent_i, :))**2
       d_mu_d_r_j_inuc(nup+1:nelec, cent_i, :) = d_param_mu(1, iwctype(cent_i), down) / (1.d0 + d_param_mu(2, iwctype(cent_i), down) * scaled_dist_en_wlk(nup+1:nelec, cent_i, :))**2
       !power series terms
       do order_p  = 2, order_mu_bf
          d_mu_d_r_j_inuc(1:nup, cent_i, :) = d_mu_d_r_j_inuc(1:nup, cent_i, :) + order_p * d_param_mu(order_p + 1, iwctype(cent_i), up) * scaled_dist_en_wlk(1:nup, cent_i, :)**(order_p - 1)
          d_mu_d_r_j_inuc(nup+1:nelec, cent_i, :) = d_mu_d_r_j_inuc(nup+1:nelec, cent_i, :) + order_p * d_param_mu(order_p + 1, iwctype(cent_i), down) * scaled_dist_en_wlk(nup+1:nelec, cent_i, :)**(order_p - 1)
       end do
    end do
    d_mu_d_r_j_inuc(:, :, :) = d_mu_d_r_j_inuc(:, :, :) * d_scaled_dist_en_wlk_d_r(:, :, :)

  end subroutine d_mu_d_r_j_inuc_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_en_j_beta_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_en_j_beta. This is a rank three array.
    !               Laplacian of the electron-electron backflow function
    !               with respect to electron position.
    !               It is necessary for the construction of the energy.
    !               The first index is for the components of xi_en (beta)
    !               The second index is the label of the ee backflow transformed electon and
    !               for the electron that the laplacian is taken with respect to (delta_ij)
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 20 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: dim_j

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_en_j_beta')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ndim')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('d_mu_d_r_j_alpha')
       call object_needed('lap_j_mu')

       return
    endif

    call object_alloc('lap_i_xi_en_j_beta', lap_i_xi_en_j_beta, ndim, nelec, nwalk)
    lap_i_xi_en_j_beta(:,:,:) = 0.d0

    do dim_j = 1, ndim
       lap_i_xi_en_j_beta(dim_j, :, :) =  sum(lap_j_mu(:, :, :) * vec_en_xyz_wlk(dim_j, :, :, :) + 2.d0 * d_mu_d_r_j_alpha(dim_j, :, :, :), 2)
    end do

  end subroutine lap_i_xi_en_j_beta_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_en_j_beta_num_bld
    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_en_j_beta_num. This is a rank three array.
    !               Numerical laplacian of xi_en:
    !
    !               (nabla_i)^2 xi_en_j^beta
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               The first index is for the components of xi_en (beta)
    !               The second index is the label of the ee backflow transformation of each electon
    !               and the electron that the laplacian is taken with respect to (delta_ij)
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 20 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j, dim_i, dim_j
    real(dp), allocatable :: d_xi_en_j_beta_d_r_i_alpha_plus(:,:,:,:), d_xi_en_j_beta_d_r_i_alpha_minus(:,:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'lap_i_xi_en_j_beta_num_bld'

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_en_j_beta_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('coord_elec_wlk')
       call object_needed('d_xi_en_j_beta_d_r_i_alpha')

       return
    endif

    ! allocation
    call object_alloc('lap_i_xi_en_j_beta_num', lap_i_xi_en_j_beta_num, ndim, nelec, nwalk)
    lap_i_xi_en_j_beta_num(:,:,:) = 0.d0

    !allocate temporary objects
    allocate(d_xi_en_j_beta_d_r_i_alpha_plus(ndim, ndim, nelec, nwalk))
    allocate(d_xi_en_j_beta_d_r_i_alpha_minus(ndim, ndim, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

    do elec_j = 1, nelec
       do dim_j = 1, ndim  !beta
          do dim_i = 1, ndim  !alpha
             coord_elec_wlk = coord_elec_wlk_temp
             !shift coordinates in a particular dimension, x,y,z, etc
             coord_elec_wlk(dim_i, elec_j, :) = coord_elec_wlk(dim_i, elec_j, :) + epsilon
             call object_modified('coord_elec_wlk')
             !get updated d_xi_en_j_beta_d_r_i_alpha
             call object_provide_in_node(lhere, 'd_xi_en_j_beta_d_r_i_alpha')
             d_xi_en_j_beta_d_r_i_alpha_plus = d_xi_en_j_beta_d_r_i_alpha
             coord_elec_wlk = coord_elec_wlk_temp
             !shift coordinates in a particular dimension, x,y,z, etc
             coord_elec_wlk(dim_i, elec_j, :) = coord_elec_wlk(dim_i, elec_j, :) - epsilon
             call object_modified('coord_elec_wlk')
             !get updated d_xi_en_j_beta_d_r_i_alpha
             call object_provide_in_node(lhere, 'd_xi_en_j_beta_d_r_i_alpha')
             d_xi_en_j_beta_d_r_i_alpha_minus = d_xi_en_j_beta_d_r_i_alpha
             !calculate derivative
             lap_i_xi_en_j_beta_num(dim_j, elec_j, :) = lap_i_xi_en_j_beta_num(dim_j, elec_j, :) + (d_xi_en_j_beta_d_r_i_alpha_plus(dim_j, dim_i, elec_j, :) - d_xi_en_j_beta_d_r_i_alpha_minus(dim_j, dim_i, elec_j, :) ) / (2.d0 * epsilon)
          end do
       end do
    end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'd_xi_en_j_beta_d_r_i_alpha')

  end subroutine lap_i_xi_en_j_beta_num_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_j_mu_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_j_mu. This is a rank three array.
    !               Laplacian of the mu:
    !
    !               (nabla_j)^2 mu(r_jI)       Note that this is derivative with respect to first index so derivatives of the cutoff
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of centers and labels the Ith center
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 20 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: cent_i, dim_i
    character(len=max_string_len_rout), save :: lhere = 'lap_j_mu_bld'

    ! header
    if (header_exe) then

       call object_create('lap_j_mu')

       call object_needed('nwalk')
       call object_needed('ndim')
       call object_needed('nelec')
       call object_needed('cent')
       call object_needed('all_elec')
       call object_needed('d_mu_d_r_j_inuc')
       call object_needed('d2_mu_d_r_j_inuc_2')
       call object_needed('dist_en_wlk')

       return
    endif

    call object_alloc('lap_j_mu', lap_j_mu, nelec, ncent, nwalk)
    lap_j_mu(:,:,:) = 0.d0

    lap_j_mu(:, :, :) = d2_mu_d_r_j_inuc_2(:, :, :) + (ndim - 1.d0) / dist_en_wlk(:, :, :) * d_mu_d_r_j_inuc(:,:,:)

    if (all_elec) then
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       call object_provide_in_node(lhere, 'd_smooth_cutoff_g_d_r')
       call object_provide_in_node(lhere, 'd2_smooth_cutoff_g_d_r_2')
       call object_provide_in_node(lhere, 'mu')
       call object_provide_in_node(lhere, 'vec_en_xyz_wlk')
       call object_provide_in_node(lhere, 'd_mu_d_r_j_alpha')

       do cent_i = 1, ncent
          lap_j_mu(:, cent_i, :) = lap_j_mu(:, cent_i, :) * product(smooth_cutoff_g(:, 1:cent_i-1, :), 2) * product(smooth_cutoff_g(:, cent_i + 1:ncent, :), 2)
          lap_j_mu(:, cent_i, :) = lap_j_mu(:, cent_i, :) + mu(:, cent_i, :) * ( sum((d2_smooth_cutoff_g_d_r_2(:, 1:cent_i-1, :) - d_smooth_cutoff_g_d_r(:, 1:cent_i-1, :) **2 / smooth_cutoff_g(:, 1:cent_i-1, :) + d_smooth_cutoff_g_d_r(:, 1:cent_i-1, :) * (ndim - 1.d0) / dist_en_wlk(:, 1:cent_i-1, :)) / smooth_cutoff_g(:, 1:cent_i-1, :), 2) + sum((d2_smooth_cutoff_g_d_r_2(:, cent_i+1:ncent, :) - d_smooth_cutoff_g_d_r(:, cent_i+1:ncent, :) **2 / smooth_cutoff_g(:, cent_i+1:ncent, :) + d_smooth_cutoff_g_d_r(:, cent_i+1:ncent, :) * (ndim - 1.d0) / dist_en_wlk(:, cent_i+1:ncent, :)) / smooth_cutoff_g(:, cent_i+1:ncent, :), 2) )
          do dim_i = 1, ndim
             lap_j_mu(:, cent_i, :) = lap_j_mu(:, cent_i, :) + (product(smooth_cutoff_g(:, 1:cent_i-1, :), 2) * product(smooth_cutoff_g(:, cent_i + 1:ncent, :), 2) * d_mu_d_r_j_inuc(:, cent_i, :) / dist_en_wlk(:, cent_i, :) * vec_en_xyz_wlk(dim_i, :, cent_i, :) + d_mu_d_r_j_alpha(dim_i, :, cent_i, :)) * ( sum(d_smooth_cutoff_g_d_r(:, 1:cent_i-1, :) / smooth_cutoff_g(:, 1:cent_i-1, :) / dist_en_wlk(:, 1:cent_i-1, :) * vec_en_xyz_wlk(dim_i, :, 1:cent_i-1, :), 2) + sum(d_smooth_cutoff_g_d_r(:, cent_i+1:ncent, :) / smooth_cutoff_g(:, cent_i+1:ncent, :) / dist_en_wlk(:, cent_i+1:ncent, :) * vec_en_xyz_wlk(dim_i, :, cent_i+1:ncent, :), 2) )
          end do
       end do
    end if

  end subroutine lap_j_mu_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d2_mu_d_r_j_inuc_2_bld

    !---------------------------------------------------------------------------
    ! Description : build d2_mu_d_r_j_inuc_2. This is a rank three array.
    !               Second derivative of mu without (without the cutoff function, used with the chain rule):
    !
    !               d^2 mu(r_j_inuc)
    !               ----------
    !               d r_j_inuc^2
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of centers and labels the inuc nucleus
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 20 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec
    include 'commons.h'

    !local
    integer  :: order_p, cent_i
    integer, parameter  :: up=1, down=2

    ! header
    if (header_exe) then

       call object_create('d2_mu_d_r_j_inuc_2')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed('nup')
       call object_needed('order_mu_bf')
       call object_needed ('scaled_dist_en_wlk')
       call object_needed ('d_scaled_dist_en_wlk_d_r')
       call object_needed ('d2_scaled_dist_en_wlk_d_r_2')
       call object_needed('d_param_mu')

       return
    endif

    call object_alloc('d2_mu_d_r_j_inuc_2', d2_mu_d_r_j_inuc_2, nelec, ncent, nwalk)
    d2_mu_d_r_j_inuc_2(:,:,:) = 0.d0

    do cent_i = 1, ncent
       !pade terms
       d2_mu_d_r_j_inuc_2(1:nup, cent_i, :) = d_param_mu(1, iwctype(cent_i), up) / (1.d0 + d_param_mu(2, iwctype(cent_i), up) * scaled_dist_en_wlk(1:nup, cent_i, :))**2
       d2_mu_d_r_j_inuc_2(nup+1:nelec, cent_i, :) = d_param_mu(1, iwctype(cent_i), down) / (1.d0 + d_param_mu(2, iwctype(cent_i), down) * scaled_dist_en_wlk(nup+1:nelec, cent_i, :))**2
       !power series terms
       do order_p  = 2, order_mu_bf
          d2_mu_d_r_j_inuc_2(1:nup, cent_i, :) = d2_mu_d_r_j_inuc_2(1:nup, cent_i, :) + order_p * d_param_mu(order_p + 1, iwctype(cent_i), up) * scaled_dist_en_wlk(1:nup, cent_i, :)**(order_p - 1)
          d2_mu_d_r_j_inuc_2(nup+1:nelec, cent_i, :) = d2_mu_d_r_j_inuc_2(nup+1:nelec, cent_i, :) + order_p * d_param_mu(order_p + 1, iwctype(cent_i), down) * scaled_dist_en_wlk(nup+1:nelec, cent_i, :)**(order_p - 1)
       end do
    end do
    d2_mu_d_r_j_inuc_2(:, :, :) = d2_mu_d_r_j_inuc_2(:, :, :) * d2_scaled_dist_en_wlk_d_r_2(:, :, :)

    do cent_i = 1, ncent
       !pade terms
       d2_mu_d_r_j_inuc_2(nup+1:nelec, cent_i, :) = d2_mu_d_r_j_inuc_2(nup+1:nelec, cent_i, :) - d_scaled_dist_en_wlk_d_r(nup+1:nelec, cent_i, :) **2 * (2.d0 * d_param_mu(1, iwctype(cent_i), down) * d_param_mu(2, iwctype(cent_i), down)/ (1.d0 + d_param_mu(2, iwctype(cent_i), down) * scaled_dist_en_wlk(nup+1:nelec, cent_i, :))**3)
       d2_mu_d_r_j_inuc_2( 1:nup, cent_i, :) = d2_mu_d_r_j_inuc_2( 1:nup, cent_i, :) - d_scaled_dist_en_wlk_d_r(1:nup, cent_i, :) **2 * (2.d0 * d_param_mu(1, iwctype(cent_i), up) * d_param_mu(2, iwctype(cent_i), up)/ (1.d0 + d_param_mu(2, iwctype(cent_i), up) * scaled_dist_en_wlk(1:nup, cent_i, :))**3)
       !power series terms
       do order_p  = 2, order_mu_bf
          d2_mu_d_r_j_inuc_2( nup+1:nelec, cent_i, :) = d2_mu_d_r_j_inuc_2( nup+1:nelec, cent_i, :) + order_p * (order_p - 1.d0) * d_param_mu(order_p + 1, iwctype(cent_i), down) * scaled_dist_en_wlk(nup+1:nelec, cent_i, :)**(order_p - 2) * d_scaled_dist_en_wlk_d_r(nup+1:nelec, cent_i, :) **2
          d2_mu_d_r_j_inuc_2( 1:nup, cent_i, :) = d2_mu_d_r_j_inuc_2( 1:nup, cent_i, :) + order_p * (order_p - 1.d0) * d_param_mu(order_p + 1, iwctype(cent_i), up) * scaled_dist_en_wlk(1:nup, cent_i, :)**(order_p - 2) * d_scaled_dist_en_wlk_d_r(1:nup, cent_i, :) **2
       end do
    end do

  end subroutine d2_mu_d_r_j_inuc_2_bld

  !==============================================================================================

  !======================================================================================

  subroutine xi_ee_bld
    !---------------------------------------------------------------------------
    ! Description : build object xi_ee(:,:). This is a rank three array.
    !               It is the electron-electron backflow.
    !               It is necessary for the construction of the total backflow transformation.
    !               The first index has dimension equal to the number of spatial dimensions, ndim.
    !               The second index has dimensions equal to the number of electrons, nelec.
    !               The third index has dimension equal to the number of walkers, nwalk.
    !
    ! Created     : F. Petruzielo, 22 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: dim_i

    ! header
    if (header_exe) then

       call object_create('xi_ee')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('ndim')
       call object_needed('eta')

       return
    endif

    call object_alloc('xi_ee', xi_ee, ndim, nelec, nwalk)
    xi_ee(:,:,:) = 0.d0

    do dim_i = 1, ndim
       xi_ee(dim_i, :, :) = sum(eta(:, :, :) * vec_ee_xyz_wlk(dim_i, :, :, :), 2)
    enddo

  end subroutine xi_ee_bld

  !======================================================================================

  !======================================================================================

  subroutine eta_bld
    !---------------------------------------------------------------------------
    ! Description : build object eta(:,:,:). This is a rank three array.
    !               It is necessary for the construction of the electron-electron backflow terms.
    !               The first index has dimension equal to the number of electrons, nelec.
    !               The second index has dimensions equal to the number of electrons, nelec.
    !               The third index has dimension equal to the number of walkers, nwalk.
    !
    ! Created     : F. Petruzielo, 22 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, order_p
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin
    character(len=max_string_len_rout), save :: lhere = 'eta_bld'

    ! header
    if (header_exe) then

       call object_create('eta')

       call object_needed('all_elec')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_eta_bf')
       call object_needed ('scaled_dist_ee_wlk')
       call object_needed('c_param_eta')
       call object_needed('asymp_eta')

       return
    endif

    call object_alloc('eta', eta, nelec, nelec, nwalk)
    eta (:,:,:) = 0.d0

    !pade term (none for parallel spin c_param_eta(1,alpha_alpha) = 0)
    eta(1:nup, nup+1:nelec, :) = c_param_eta(1, up_down) * scaled_dist_ee_wlk(1:nup, nup+1:nelec, :) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(1:nup, nup+1:nelec, :))
    eta(nup+1:nelec, 1:nup, :) = c_param_eta(1, up_down) * scaled_dist_ee_wlk(nup+1:nelec, 1:nup, :) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(nup+1:nelec, 1:nup, :))
    !power series terms
    do order_p=2, order_eta_bf
       eta(1:nup, 1:nup, :) = eta(1:nup, 1:nup, :) + c_param_eta(order_p + 1, up_up) * scaled_dist_ee_wlk(1:nup, 1:nup, :)**order_p
       eta(nup+1:nelec, nup+1:nelec, :) = eta(nup+1:nelec, nup+1:nelec, :) + c_param_eta(order_p + 1, down_down) * scaled_dist_ee_wlk(nup+1:nelec, nup+1:nelec, :)**order_p
       eta(nup+1:nelec, 1:nup, :) = eta(nup+1:nelec, 1:nup, :) + c_param_eta(order_p + 1, up_down) * scaled_dist_ee_wlk(nup+1:nelec, 1:nup, :)**order_p
       eta(1:nup, nup+1:nelec, :) = eta(1:nup, nup+1:nelec, :) + c_param_eta(order_p + 1, up_down) * scaled_dist_ee_wlk(1:nup, nup+1:nelec, :)**order_p
    enddo

    !subtract off asymptotic value
    eta(1:nup,1:nup,:) = eta(1:nup,1:nup,:) - asymp_eta(up_up)
    eta(nup + 1:nelec, nup + 1:nelec, :) = eta(nup + 1:nelec, nup + 1:nelec, :) - asymp_eta(down_down)
    eta(1:nup, nup + 1:nelec, :) = eta(1:nup, nup + 1:nelec, :) - asymp_eta(up_down)
    eta(nup + 1:nelec, 1:nup, :) = eta(nup + 1:nelec, 1:nup, :) - asymp_eta(up_down)

    if (all_elec) then
       !apply smooth cutoff function g
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       do elec_j = 1, nelec
          eta(:, elec_j, :) = eta(:, elec_j, :) * product(smooth_cutoff_g(:, :, :), 2)
       end do
    end if

    !no self backflow
    do elec_j = 1, nelec
       eta(elec_j, elec_j, :) = 0.d0
    enddo

  end subroutine eta_bld

  !=============================================================================================

  !==============================================================================================

  subroutine asymp_eta_bld
    !---------------------------------------------------------------------------
    ! Description : build object asymp_eta(:). This is a rank one array.
    !               The first index is the spin dependence.
    !               Give the asymptotic value of the pade/power series expansion used in eta.
    !               This value is subtracted off from the pade/power series expansion such that eta goes to zero at infinity.
    !
    ! Created     : F. Petruzielo, 22 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_p
    integer, parameter :: spin_dependencies_nb = 3
    integer, parameter  :: up_up=1, down_down=2, up_down=3

    ! header
    if (header_exe) then

       call object_create('asymp_eta')

       call object_needed('order_eta_bf')
       call object_needed('c_param_eta')
       call object_needed('asymp_scaled_dist_two_body')

       return
    endif

    call object_alloc('asymp_eta', asymp_eta, spin_dependencies_nb)
    asymp_eta(:) = 0.d0

    !power series terms
    do order_p=2, order_eta_bf
       asymp_eta(up_up) = asymp_eta(up_up) + c_param_eta(order_p + 1, up_up) * asymp_scaled_dist_two_body**order_p
       asymp_eta(down_down) = asymp_eta(down_down) + c_param_eta(order_p + 1, down_down) * asymp_scaled_dist_two_body**order_p
       asymp_eta(up_down) = asymp_eta(up_down) + c_param_eta(order_p + 1, up_down) * asymp_scaled_dist_two_body**order_p
    enddo

    !pade term, none for up_up or down_down because c_param_eta(1, alpha_alpha) = 0 for electron-electron cusp condition to be satisfied
    asymp_eta(up_down) = asymp_eta(up_down) + c_param_eta(1, up_down) * asymp_scaled_dist_two_body / (1.d0 + c_param_eta(2, up_down) * asymp_scaled_dist_two_body)

  end subroutine asymp_eta_bld

  !==============================================================================================

  !==============================================================================================

  subroutine c_param_eta_bld

    !---------------------------------------------------------------------------
    ! Description : build c_param_eta. This is a rank 2 array.
    !               List of c parameters for eta backflow function for each spin dependence
    !               The first index is which parameter.
    !               The second index is the type of spin relationship
    !
    !               cusp condition constraints:   d eta  |
    !                                             -----  |        = 0  for parallel spin electrons
    !                                             d r_ij |
    !                                                     r_ij=0
    !
    ! Created     : F. Petruzielo, 22 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: c_param_nb
    integer, parameter :: spin_dependencies_nb = 3
    character(len=max_string_len_rout), save :: lhere = 'c_param_eta_bld'

    ! header
    if (header_exe) then

       call object_create('c_param_eta')

       call object_needed('order_eta_bf')
       call object_needed('asymp_scaled_dist_two_body')
       call object_needed('eta_spin_dependence')
       call object_needed('read_c_param')
       call object_needed('read_c_param_eta_nb')

       return
    endif

    !need to have the option of using read_c_param_eta if it is read in
    if (read_c_param) then
       call object_provide_in_node(lhere, 'read_c_param_eta')
    end if

    call object_alloc('c_param_eta', c_param_eta, order_eta_bf + 1, spin_dependencies_nb)
    c_param_eta = 0.d0

    !do we read in parameters?
    if (read_c_param) then
       if (eta_spin_dependence == spin_dependencies_nb) then
          c_param_nb = (order_eta_bf + 1)  * spin_dependencies_nb
          if (read_c_param_eta_nb == c_param_nb) then
             c_param_eta = reshape( read_c_param_eta, (/ order_eta_bf + 1, spin_dependencies_nb /))
          else
             call require(lhere, 'number of parameters in read_c_param_eta equal to (order_eta_bf + 1) * 3 specified', .false.)
          end if
       elseif( eta_spin_dependence == spin_dependencies_nb -1 ) then
          c_param_nb = (order_eta_bf + 1)  * (spin_dependencies_nb-1)
          if (read_c_param_eta_nb == c_param_nb) then
             c_param_eta = reshape( read_c_param_eta, (/ order_eta_bf + 1, spin_dependencies_nb /), pad=read_c_param_eta )
             !want up_up in c_param_eta(:,1), up_up=down_down in c_param_eta(:,2), up_down in c_param_eta(:,3)
             c_param_eta(:,3)=c_param_eta(:,2)
             c_param_eta(:,2)=c_param_eta(:,1)
          else
             call require(lhere, 'number of parameters in read_c_param_eta equal to (order_eta_bf + 1) * 2 specified', .false.)
          end if
       else
          c_param_nb = (order_eta_bf + 1)
          if (read_c_param_eta_nb == c_param_nb) then
             c_param_eta = reshape( read_c_param_eta, (/ order_eta_bf + 1, spin_dependencies_nb /), pad=read_c_param_eta )
          else
             call require(lhere, 'number of parameters in read_c_param_eta equal to (order_eta_bf + 1) specified', .false.)
          end if
       end if
       !don't want to re-read
       read_c_param = .false.
    else
       call random_number(c_param_eta(:,:))
       if ( eta_spin_dependence .eq. 2) then
          c_param_eta(:,2) = c_param_eta(:,1)
       elseif( eta_spin_dependence .eq. 1) then
          c_param_eta(:,2) = c_param_eta(:,1)
          c_param_eta(:,3) = c_param_eta(:,1)
       end if
    endif

    !impose electron electron cusp conditions
    !Note that we assume the two-body scaling functions evaluated at 0 are 0 and have unit slope at 0.
    c_param_eta(1, 1) = 0.d0
    c_param_eta(1, 2) = 0.d0
    if (eta_spin_dependence .eq. 1) then
       c_param_eta(1, 3) = 0.d0
    endif

  end subroutine c_param_eta_bld

  !======================================================================================

  !==============================================================================================

  subroutine d_xi_ee_j_beta_d_r_i_alpha_bld

    !---------------------------------------------------------------------------
    ! Description : build d_xi_ee_j_beta_d_r_i_alpha. This is a rank five array.
    !               Derivative of the electron-electron backflow function
    !               with respect to electron position.
    !               It is necessary for the construction of the energy, drift velocity.
    !               The first and second indices are the number of spatial dimensions.
    !               The first index is for the components of xi_ee (beta)
    !               The second index is for the components of the gradient (alpha)
    !               The third and fourth indices have dimensions equal to the number of electrons.
    !               The third index is the label of the ee backflow transformed electon
    !               The fourth index is for the electron that the gradient is taken with respect to
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 4 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j, dim_i, dim_j

    ! header
    if (header_exe) then

       call object_create('d_xi_ee_j_beta_d_r_i_alpha')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('ndim')
       call object_needed('all_elec')
       call object_needed('d_eta_d_r_j_alpha_first')
       call object_needed('d_eta_d_r_i_alpha_second')
       call object_needed('eta')

       return
    endif

    call object_alloc('d_xi_ee_j_beta_d_r_i_alpha', d_xi_ee_j_beta_d_r_i_alpha, ndim, ndim, nelec, nelec, nwalk)
    d_xi_ee_j_beta_d_r_i_alpha(:,:,:,:,:) = 0.d0

    do elec_j = 1, nelec
       do dim_i = 1, ndim     !alpha
          !delta_{alpha,beta}
          d_xi_ee_j_beta_d_r_i_alpha(dim_i, dim_i, elec_j, elec_j, :) =  sum(eta(elec_j, :, :),1)
          do dim_j = 1, ndim   !beta
             d_xi_ee_j_beta_d_r_i_alpha(dim_j, dim_i, elec_j, elec_j, :) = d_xi_ee_j_beta_d_r_i_alpha(dim_j, dim_i, elec_j, elec_j, :) + sum(d_eta_d_r_j_alpha_first(dim_i, elec_j, :, :) * vec_ee_xyz_wlk(dim_j, elec_j, :, :), 1)
          end do
       end do
    end do

    do dim_i = 1, ndim     !alpha
       !delta_{alpha,beta}
       d_xi_ee_j_beta_d_r_i_alpha(dim_i, dim_i, :, :, :) = d_xi_ee_j_beta_d_r_i_alpha(dim_i, dim_i, :, :, :) - eta(:, :, :)
       do dim_j = 1, ndim   !beta
          d_xi_ee_j_beta_d_r_i_alpha(dim_j, dim_i, :, :, :) = d_xi_ee_j_beta_d_r_i_alpha(dim_j, dim_i, :, :, :) + d_eta_d_r_i_alpha_second(dim_i, :, :, :) * vec_ee_xyz_wlk(dim_j, :, :, :)
       end do
    end do

  end subroutine d_xi_ee_j_beta_d_r_i_alpha_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_xi_ee_j_beta_d_r_i_alpha_num_bld
    !---------------------------------------------------------------------------
    ! Description : build d_xi_ee_j_beta_d_r_i_alpha__num. This is a rank five array.
    !               Numerical first derivative of the xi_ee:
    !
    !               d xi_ee
    !               ----------
    !               d r_i_alpha
    !
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               The first and second indices are the number of spatial dimensions.
    !               The first index is for the components of xi_ee (beta)
    !               The second index is for the components of the gradient (alpha)
    !               The third and fourth indices have dimensions equal to the number of electrons.
    !               The third index is the label of the ee backflow transformation of each electon
    !               The fourth index is for the electron that the gradient is taken with respect to
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 4 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, dim_i, dim_j
    real(dp), allocatable :: xi_ee_temp_plus(:,:,:), xi_ee_temp_minus(:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'd_xi_ee_j_beta_d_r_i_alpha_num_bld'

    ! header
    if (header_exe) then

       call object_create('d_xi_ee_j_beta_d_r_i_alpha_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('xi_ee')
       call object_needed('coord_elec_wlk')

       return
    endif

    ! allocation
    call object_alloc('d_xi_ee_j_beta_d_r_i_alpha_num', d_xi_ee_j_beta_d_r_i_alpha_num, ndim, ndim, nelec, nelec, nwalk)
    d_xi_ee_j_beta_d_r_i_alpha_num(:,:,:,:,:) = 0.d0

    !allocate temporary objects
    allocate(xi_ee_temp_plus(ndim, nelec, nwalk))
    allocate(xi_ee_temp_minus(ndim, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

    do elec_i = 1, nelec
       do elec_j = 1, nelec
          do dim_j = 1, ndim  !beta
             do dim_i = 1, ndim  !alpha
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) + epsilon
                call object_modified('coord_elec_wlk')
                !get updated xi_ee
                call object_provide_in_node(lhere, 'xi_ee')
                xi_ee_temp_plus = xi_ee
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) - epsilon
                call object_modified('coord_elec_wlk')
                !get updated xi_ee
                call object_provide_in_node(lhere, 'xi_ee')
                xi_ee_temp_minus = xi_ee
                !calculate derivative
                d_xi_ee_j_beta_d_r_i_alpha_num(dim_j, dim_i, elec_j, elec_i, :) = (xi_ee_temp_plus(dim_j, elec_j, :) - xi_ee_temp_minus(dim_j, elec_j, :) ) / (2.d0* epsilon)
             end do
          end do
       end do
    end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'xi_ee')

  end subroutine d_xi_ee_j_beta_d_r_i_alpha_num_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_eta_d_r_i_alpha_second_bld

    !---------------------------------------------------------------------------
    ! Description : build d_eta_d_r_i_alpha_second. This is a rank four array.
    !               First derivative of the eta:
    !
    !               d eta(r_ji)       Note that this is derivative with respect to second index so no derivatives of the cutoff
    !               ----------
    !               d r_i_alpha
    !
    !               Note alpha refers to cartesian coords. In 3-d, r_alpha takes on values x,y,z.
    !               The first index is the over spatial dimension. This is the components of the gradient.
    !               The second index has dimensions equal to the number of electrons and labels the jth electron
    !               The third index has dimensions equal to the number of electrons and labels the ith electron
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 6 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j, dim_i
    character(len=max_string_len_rout), save :: lhere = 'd_eta_d_r_i_alpha_second_bld'

    ! header
    if (header_exe) then

       call object_create('d_eta_d_r_i_alpha_second')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed ('dist_ee_wlk')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('d_eta_d_r_ji')
       call object_needed('all_elec')

       return
    endif

    call object_alloc('d_eta_d_r_i_alpha_second', d_eta_d_r_i_alpha_second, ndim, nelec, nelec, nwalk)
    d_eta_d_r_i_alpha_second(:,:,:,:) = 0.d0

    do elec_j = 1, nelec
       do dim_i = 1, ndim
          d_eta_d_r_i_alpha_second(dim_i, elec_j, 1:elec_j-1, :) =  - d_eta_d_r_ji(elec_j, 1:elec_j-1, :) / dist_ee_wlk(elec_j, 1:elec_j-1, :) * vec_ee_xyz_wlk(dim_i, elec_j, 1:elec_j-1, :)
          d_eta_d_r_i_alpha_second(dim_i, elec_j, elec_j+1:nelec, :) =  - d_eta_d_r_ji(elec_j, elec_j+1:nelec, :) / dist_ee_wlk(elec_j, elec_j+1:nelec, :) * vec_ee_xyz_wlk(dim_i, elec_j, elec_j+1:nelec, :)
      end do
    end do

    if (all_elec) then
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       do elec_j = 1, nelec
          do dim_i = 1, ndim
             !cutoff
             d_eta_d_r_i_alpha_second(dim_i, :, elec_j, :) = d_eta_d_r_i_alpha_second(dim_i, :, elec_j, :) * product(smooth_cutoff_g(:, :, :), 2)
          end do
       end do
    end if

  end subroutine d_eta_d_r_i_alpha_second_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_eta_d_r_j_alpha_first_bld

    !---------------------------------------------------------------------------
    ! Description : build d_eta_d_r_j_alpha_first. This is a rank four array.
    !               First derivative of eta:
    !
    !               d eta_(r_jl)       Note that this is derivative with respect to first index so derivatives of the cutoff (in all electron case)
    !               ----------
    !               d r_j_alpha
    !
    !               Note alpha refers to cartesian coords. In 3-d, r_alpha takes on values x,y,z.
    !               The first index is the over spatial dimension. This is the components of the gradient.
    !               The second index has dimensions equal to the number of electrons and labels the jth electron
    !               The third index has dimensions equal to the number of electrons and labels the lth electron
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 9 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim, ncent
    include 'commons.h'

    !local
    integer  :: elec_l, dim_i
    character(len=max_string_len_rout), save :: lhere = 'd_eta_d_r_j_alpha_first_bld'

    ! header
    if (header_exe) then

       call object_create('d_eta_d_r_j_alpha_first')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('d_eta_d_r_i_alpha_second')
       call object_needed('all_elec')

       return
    endif

    call object_alloc('d_eta_d_r_j_alpha_first', d_eta_d_r_j_alpha_first, ndim, nelec, nelec, nwalk)
    d_eta_d_r_j_alpha_first(:,:,:,:) = 0.d0

    d_eta_d_r_j_alpha_first = - d_eta_d_r_i_alpha_second

    if (all_elec) then
       call object_provide_in_node(lhere, 'dist_en_wlk')
       call object_provide_in_node(lhere, 'vec_en_xyz_wlk')
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       call object_provide_in_node(lhere, 'd_smooth_cutoff_g_d_r')
       call object_provide_in_node(lhere, 'eta')

       do elec_l = 1, nelec
          do dim_i = 1, ndim
             d_eta_d_r_j_alpha_first(dim_i, :, elec_l, :) = d_eta_d_r_j_alpha_first(dim_i, :, elec_l, :) + eta(:, elec_l, :) * sum(d_smooth_cutoff_g_d_r(:, :, :) *  vec_en_xyz_wlk(dim_i, :, :, :) / smooth_cutoff_g(:, :, :) / dist_en_wlk(:, :, :), 2)
          end do
       end do

    end if

  end subroutine d_eta_d_r_j_alpha_first_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_eta_d_r_ji_bld

    !---------------------------------------------------------------------------
    ! Description : build d_eta_d_r_ji. This is a rank three array.
    !               First derivative of eta without (without the cutoff function, used with the chain rule):
    !
    !               d eta(r_ji)
    !               ----------
    !               d r_ji
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of electrons and labels the ith electron
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 6 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec
    include 'commons.h'

    !local
    integer  :: order_p
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d_eta_d_r_ji')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_eta_bf')
       call object_needed ('scaled_dist_ee_wlk')
       call object_needed ('d_scaled_dist_ee_wlk_d_r')
       call object_needed('c_param_eta')

       return
    endif

    call object_alloc('d_eta_d_r_ji', d_eta_d_r_ji, nelec, nelec, nwalk)
    d_eta_d_r_ji(:,:,:) = 0.d0

    !pade terms
    d_eta_d_r_ji( nup+1:nelec, 1:nup, :) = c_param_eta(1, up_down) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(nup+1:nelec, 1:nup, :))**2
    d_eta_d_r_ji( 1:nup, nup+1:nelec, :) = c_param_eta(1, up_down) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(1:nup, nup+1:nelec, :))**2
    !power series terms
    do order_p  = 2, order_eta_bf
       d_eta_d_r_ji( nup+1:nelec, 1:nup, :) = d_eta_d_r_ji( nup+1:nelec, 1:nup, :) + order_p * c_param_eta(order_p + 1, up_down) * scaled_dist_ee_wlk(nup+1:nelec, 1:nup, :)**(order_p - 1)
       d_eta_d_r_ji( 1:nup, nup+1:nelec, :) = d_eta_d_r_ji( 1:nup, nup+1:nelec, :) + order_p * c_param_eta(order_p + 1, up_down) * scaled_dist_ee_wlk(1:nup, nup+1:nelec, :)**(order_p - 1)
       d_eta_d_r_ji(1:nup, 1:nup, :) = d_eta_d_r_ji( 1:nup, 1:nup, :) + order_p * c_param_eta(order_p + 1, up_up) * scaled_dist_ee_wlk(1:nup, 1:nup, :)**(order_p - 1)
       d_eta_d_r_ji( nup+1:nelec, nup+1:nelec, :) = d_eta_d_r_ji( nup+1:nelec, nup+1:nelec, :) + order_p * c_param_eta(order_p + 1, down_down) * scaled_dist_ee_wlk(nup+1:nelec, nup+1:nelec, :)**(order_p - 1)
    end do

    !scale
    d_eta_d_r_ji(:, :, :) = d_eta_d_r_ji(:, :, :) * d_scaled_dist_ee_wlk_d_r(:, :, :)

  end subroutine d_eta_d_r_ji_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d2_eta_d_r_ji_2_bld

    !---------------------------------------------------------------------------
    ! Description : build d2_eta_d_r_ji_2. This is a rank three array.
    !               Second derivative of eta without (without the cutoff function, used with the chain rule):
    !
    !               d^2 eta(r_ji)
    !               ----------
    !               d r_ji^2
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of electrons and labels the ith electron
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 17 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec
    include 'commons.h'

    !local
    integer  :: order_p, elec_i
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d2_eta_d_r_ji_2')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_eta_bf')
       call object_needed ('scaled_dist_ee_wlk')
       call object_needed ('d_scaled_dist_ee_wlk_d_r')
       call object_needed ('d2_scaled_dist_ee_wlk_d_r_2')
       call object_needed('c_param_eta')

       return
    endif

    call object_alloc('d2_eta_d_r_ji_2', d2_eta_d_r_ji_2, nelec, nelec, nwalk)
    d2_eta_d_r_ji_2(:,:,:) = 0.d0

    !pade terms
    d2_eta_d_r_ji_2( nup+1:nelec, 1:nup, :) = c_param_eta(1, up_down) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(nup+1:nelec, 1:nup, :))**2
    d2_eta_d_r_ji_2( 1:nup, nup+1:nelec, :) = c_param_eta(1, up_down) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(1:nup, nup+1:nelec, :))**2
    !power series terms
    do order_p  = 2, order_eta_bf
       d2_eta_d_r_ji_2( nup+1:nelec, 1:nup, :) = d2_eta_d_r_ji_2( nup+1:nelec, 1:nup, :) + order_p * c_param_eta(order_p + 1, up_down) * scaled_dist_ee_wlk(nup+1:nelec, 1:nup, :)**(order_p - 1)
       d2_eta_d_r_ji_2( 1:nup, nup+1:nelec, :) = d2_eta_d_r_ji_2( 1:nup, nup+1:nelec, :) + order_p * c_param_eta(order_p + 1, up_down) * scaled_dist_ee_wlk(1:nup, nup+1:nelec, :)**(order_p - 1)
       d2_eta_d_r_ji_2(1:nup, 1:nup, :) = d2_eta_d_r_ji_2( 1:nup, 1:nup, :) + order_p * c_param_eta(order_p + 1, up_up) * scaled_dist_ee_wlk(1:nup, 1:nup, :)**(order_p - 1)
       d2_eta_d_r_ji_2( nup+1:nelec, nup+1:nelec, :) = d2_eta_d_r_ji_2( nup+1:nelec, nup+1:nelec, :) + order_p * c_param_eta(order_p + 1, down_down) * scaled_dist_ee_wlk(nup+1:nelec, nup+1:nelec, :)**(order_p - 1)
    end do

    !scale
    d2_eta_d_r_ji_2(:, :, :) = d2_eta_d_r_ji_2(:, :, :) * d2_scaled_dist_ee_wlk_d_r_2(:, :, :)

    !pade terms
    d2_eta_d_r_ji_2( nup+1:nelec, 1:nup, :) = d2_eta_d_r_ji_2( nup+1:nelec, 1:nup, :) - d_scaled_dist_ee_wlk_d_r(nup+1:nelec, 1:nup, :) **2 * (2.d0 * c_param_eta(1, up_down) * c_param_eta(2, up_down)/ (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(nup+1:nelec, 1:nup, :))**3)
    d2_eta_d_r_ji_2( 1:nup, nup+1:nelec, :) = d2_eta_d_r_ji_2( 1:nup, nup+1:nelec, :) - d_scaled_dist_ee_wlk_d_r(1:nup, nup+1:nelec, :) **2 * (2.d0 * c_param_eta(1, up_down) * c_param_eta(2, up_down)/ (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(1:nup, nup+1:nelec, :))**3)
    !power series terms
    do order_p  = 2, order_eta_bf
       d2_eta_d_r_ji_2( nup+1:nelec, 1:nup, :) = d2_eta_d_r_ji_2( nup+1:nelec, 1:nup, :) + order_p * (order_p - 1.d0) * c_param_eta(order_p + 1, up_down) * scaled_dist_ee_wlk(nup+1:nelec, 1:nup, :)**(order_p - 2) * d_scaled_dist_ee_wlk_d_r(nup+1:nelec, 1:nup, :) **2
       d2_eta_d_r_ji_2( 1:nup, 1:nup, :) = d2_eta_d_r_ji_2( 1:nup, 1:nup, :) + order_p * (order_p - 1.d0) * c_param_eta(order_p + 1, up_up) * scaled_dist_ee_wlk(1:nup, 1:nup, :)**(order_p - 2) * d_scaled_dist_ee_wlk_d_r(1:nup, 1:nup, :) **2
       d2_eta_d_r_ji_2( nup+1:nelec, nup+1:nelec, :) = d2_eta_d_r_ji_2( nup+1:nelec, nup+1:nelec, :) + order_p * (order_p - 1.d0) * c_param_eta(order_p + 1, down_down) * scaled_dist_ee_wlk(nup+1:nelec, nup+1:nelec, :)**(order_p - 2) * d_scaled_dist_ee_wlk_d_r(nup+1:nelec, nup+1:nelec, :) **2
       d2_eta_d_r_ji_2( 1:nup, nup+1:nelec, :) = d2_eta_d_r_ji_2( 1:nup, nup+1:nelec, :) + order_p * (order_p - 1.d0) * c_param_eta(order_p + 1, up_down) * scaled_dist_ee_wlk(1:nup, nup+1:nelec, :)**(order_p - 2) * d_scaled_dist_ee_wlk_d_r(1:nup, nup+1:nelec, :) **2
    end do

    !no self backflow
    do elec_i = 1, nelec
       d2_eta_d_r_ji_2(elec_i, elec_i, :) = 0.d0
    end do

  end subroutine d2_eta_d_r_ji_2_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_ee_j_beta_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_ee_j_beta. This is a rank four array.
    !               Laplacian of the electron-electron backflow function
    !               with respect to electron position.
    !               It is necessary for the construction of the energy.
    !               The first index is for the components of xi_ee (beta)
    !               The second index is the label of the ee backflow transformed electon
    !               The third index is for the electron that the laplacian is taken with respect to
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 17 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j, dim_j

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_ee_j_beta')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('ndim')
       call object_needed('d_eta_d_r_j_alpha_first')
       call object_needed('d_eta_d_r_i_alpha_second')
       call object_needed('lap_j_eta_first')
       call object_needed('lap_i_eta_second')

       return
    endif

    call object_alloc('lap_i_xi_ee_j_beta', lap_i_xi_ee_j_beta, ndim, nelec, nelec, nwalk)
    lap_i_xi_ee_j_beta(:,:,:,:) = 0.d0

    do elec_j = 1, nelec
       do dim_j = 1, ndim
          lap_i_xi_ee_j_beta(dim_j, elec_j, elec_j, :) =  sum(lap_j_eta_first(elec_j,:,:) * vec_ee_xyz_wlk(dim_j, elec_j, :, :) + 2.d0 * d_eta_d_r_j_alpha_first(dim_j, elec_j, :, :),1)

       end do
    end do

    do dim_j = 1, ndim
       lap_i_xi_ee_j_beta(dim_j, :, :, :) = lap_i_xi_ee_j_beta(dim_j, :, :, :) + lap_i_eta_second(:,:,:) * vec_ee_xyz_wlk(dim_j, :, :, :) - 2.d0 * d_eta_d_r_i_alpha_second(dim_j, :, :, :)
    end do

  end subroutine lap_i_xi_ee_j_beta_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_ee_j_beta_num_bld
    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_ee_j_beta_num. This is a rank four array.
    !               Numerical laplacian of xi_ee:
    !
    !               (nabla_i)^2 xi_ee_j^beta
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               The first index is for the components of xi_ee (beta)
    !               The second index is the label of the ee backflow transformation of each electon
    !               The third index is for the electron that the laplacian is taken with respect to
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 18 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, dim_i, dim_j
    real(dp), allocatable :: d_xi_ee_j_beta_d_r_i_alpha_plus(:,:,:,:,:), d_xi_ee_j_beta_d_r_i_alpha_minus(:,:,:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'lap_i_xi_ee_j_beta_num_bld'

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_ee_j_beta_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('coord_elec_wlk')
       call object_needed('d_xi_ee_j_beta_d_r_i_alpha')

       return
    endif

    ! allocation
    call object_alloc('lap_i_xi_ee_j_beta_num', lap_i_xi_ee_j_beta_num, ndim, nelec, nelec, nwalk)
    lap_i_xi_ee_j_beta_num(:,:,:,:) = 0.d0

    !allocate temporary objects
    allocate(d_xi_ee_j_beta_d_r_i_alpha_plus(ndim, ndim, nelec, nelec, nwalk))
    allocate(d_xi_ee_j_beta_d_r_i_alpha_minus(ndim, ndim, nelec, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

    do elec_i = 1, nelec
       do elec_j = 1, nelec
          do dim_j = 1, ndim  !beta
             do dim_i = 1, ndim  !alpha
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) + epsilon
                call object_modified('coord_elec_wlk')
                !get updated d_xi_ee_j_beta_d_r_i_alpha
                call object_provide_in_node(lhere, 'd_xi_ee_j_beta_d_r_i_alpha')
                d_xi_ee_j_beta_d_r_i_alpha_plus = d_xi_ee_j_beta_d_r_i_alpha
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) - epsilon
                call object_modified('coord_elec_wlk')
                !get updated d_xi_ee_j_beta_d_r_i_alpha
                call object_provide_in_node(lhere, 'd_xi_ee_j_beta_d_r_i_alpha')
                d_xi_ee_j_beta_d_r_i_alpha_minus = d_xi_ee_j_beta_d_r_i_alpha
                !calculate derivative
                lap_i_xi_ee_j_beta_num(dim_j, elec_j, elec_i, :) = lap_i_xi_ee_j_beta_num(dim_j, elec_j, elec_i, :) + (d_xi_ee_j_beta_d_r_i_alpha_plus(dim_j, dim_i, elec_j, elec_i, :) - d_xi_ee_j_beta_d_r_i_alpha_minus(dim_j, dim_i, elec_j, elec_i, :) ) / (2.d0 * epsilon)
             end do
          end do
       end do
    end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'd_xi_ee_j_beta_d_r_i_alpha')

  end subroutine lap_i_xi_ee_j_beta_num_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_eta_second_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_i_eta_second. This is a rank three array.
    !               Laplacian of the eta:
    !
    !               (nabla_i)^2 eta(r_ji)       Note that this is derivative with respect to second index so no derivatives of the cutoff
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of electrons and labels the ith electron
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 17 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j
    character(len=max_string_len_rout), save :: lhere = 'lap_i_eta_second_bld'

    ! header
    if (header_exe) then

       call object_create('lap_i_eta_second')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ndim')
       call object_needed('all_elec')
       call object_needed ('dist_ee_wlk')
       call object_needed('d_eta_d_r_ji')
       call object_needed('d2_eta_d_r_ji_2')

       return
    endif

    call object_alloc('lap_i_eta_second', lap_i_eta_second, nelec, nelec, nwalk)
    lap_i_eta_second(:,:,:) = 0.d0

    do elec_j = 1, nelec
       lap_i_eta_second(elec_j, 1:elec_j-1, :) = d2_eta_d_r_ji_2(elec_j, 1:elec_j-1, :) + (ndim - 1) * d_eta_d_r_ji(elec_j, 1:elec_j-1, :) / dist_ee_wlk(elec_j, 1:elec_j-1, :)
       lap_i_eta_second(elec_j, elec_j+1:nelec, :) = d2_eta_d_r_ji_2(elec_j, elec_j+1:nelec, :) + (ndim - 1) * d_eta_d_r_ji(elec_j, elec_j+1:nelec, :) / dist_ee_wlk(elec_j, elec_j+1:nelec, :)
    end do

    if (all_elec) then
       !cutoff
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       do elec_j = 1, nelec
          lap_i_eta_second(:, elec_j, :) = lap_i_eta_second(:, elec_j, :) * product(smooth_cutoff_g(:, :, :), 2)
       end do
    end if

  end subroutine lap_i_eta_second_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_j_eta_first_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_j_eta_first. This is a rank three array.
    !               Laplacian of the eta:
    !
    !               (nabla_j)^2 eta(r_jl)       Note that this is derivative with respect to first index so derivatives of the cutoff
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of electrons and labels the lth electron
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 18 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_l, dim_i
    character(len=max_string_len_rout), save :: lhere = 'lap_j_eta_first_bld'

    ! header
    if (header_exe) then

       call object_create('lap_j_eta_first')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('all_elec')
       call object_needed('lap_i_eta_second')
       call object_needed('ndim')

       return
    endif

    call object_alloc('lap_j_eta_first', lap_j_eta_first, nelec, nelec, nwalk)
    lap_j_eta_first(:,:,:) = 0.d0

    lap_j_eta_first = lap_i_eta_second

    if (all_elec) then
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       call object_provide_in_node(lhere, 'd_smooth_cutoff_g_d_r')
       call object_provide_in_node(lhere, 'd2_smooth_cutoff_g_d_r_2')
       call object_provide_in_node(lhere, 'eta')
       call object_provide_in_node(lhere, 'd_eta_d_r_j_alpha_first')
       call object_provide_in_node(lhere, 'd_eta_d_r_ji')
       call object_provide_in_node(lhere, 'dist_ee_wlk')
       call object_provide_in_node(lhere, 'dist_en_wlk')
       call object_provide_in_node(lhere, 'vec_en_xyz_wlk')
       call object_provide_in_node(lhere, 'vec_ee_xyz_wlk')

       do elec_l = 1, nelec
          lap_j_eta_first(:, elec_l, :) = lap_j_eta_first(:, elec_l, :) + eta(:, elec_l, :) * sum((d2_smooth_cutoff_g_d_r_2(:, :, :) - d_smooth_cutoff_g_d_r(:, :, :) ** 2 / smooth_cutoff_g(:, :, :) + d_smooth_cutoff_g_d_r(:, :, :) * (ndim - 1.d0)/ dist_en_wlk(:, :, :)) / smooth_cutoff_g(:, :, :), 2)
          do dim_i = 1, ndim
             lap_j_eta_first(:, elec_l, :) = lap_j_eta_first(:, elec_l, :) + d_eta_d_r_j_alpha_first(dim_i, :, elec_l, :) * sum(d_smooth_cutoff_g_d_r(:, :, :) / smooth_cutoff_g(:, :, :) / dist_en_wlk(:, :, :) * vec_en_xyz_wlk(dim_i, :, : ,:), 2)
             lap_j_eta_first(1:elec_l-1, elec_l, :) = lap_j_eta_first(1:elec_l-1, elec_l, :) + d_eta_d_r_ji(1:elec_l-1, elec_l, :) / dist_ee_wlk(1:elec_l-1, elec_l, :) * vec_ee_xyz_wlk(dim_i, 1:elec_l-1, elec_l, :) * product(smooth_cutoff_g(1:elec_l-1, :, :), 2) *  sum(d_smooth_cutoff_g_d_r(1:elec_l-1, :, :) / smooth_cutoff_g(1:elec_l-1, :, :) / dist_en_wlk(1:elec_l-1, :, :) * vec_en_xyz_wlk(dim_i, 1:elec_l-1, : ,:), 2)
             lap_j_eta_first(elec_l+1:nelec, elec_l, :) = lap_j_eta_first(elec_l+1:nelec, elec_l, :) + d_eta_d_r_ji(elec_l+1:nelec, elec_l, :) / dist_ee_wlk(elec_l+1:nelec, elec_l, :) * vec_ee_xyz_wlk(dim_i, elec_l+1:nelec, elec_l, :) * product(smooth_cutoff_g(elec_l+1:nelec, :, :), 2) *  sum(d_smooth_cutoff_g_d_r(elec_l+1:nelec, :, :) / smooth_cutoff_g(elec_l+1:nelec, :, :) / dist_en_wlk(elec_l+1:nelec, :, :) * vec_en_xyz_wlk(dim_i, elec_l+1:nelec, : ,:), 2)
          end do
       end do
    end if

  end subroutine lap_j_eta_first_bld

  !==============================================================================================

  !======================================================================================

  subroutine xi_een_phi_bld
    !---------------------------------------------------------------------------
    ! Description : build object xi_een_phi(:,:). This is a rank three array.
    !               It is part of electron-electron-nuclear backflow.
    !               The first index has dimension equal to the number of spatial dimensions.
    !               The second index has dimensions equal to the number of electrons.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 23 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, nelec, ndim
    include 'commons.h'

    !local
    integer  :: dim_i

    ! header
    if (header_exe) then

       call object_create('xi_een_phi')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('nelec')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('ndim')
       call object_needed('phi')

       return
    endif

    call object_alloc('xi_een_phi', xi_een_phi, ndim, nelec, nwalk)
    xi_een_phi(:,:,:) = 0.d0

    do dim_i = 1, ndim
       xi_een_phi(dim_i, :, :) = sum( sum(phi(:, :, :, :),3) * vec_ee_xyz_wlk(dim_i, :, :, :), 2)
    enddo

  end subroutine xi_een_phi_bld

  !======================================================================================

  !======================================================================================

  subroutine phi_bld
    !---------------------------------------------------------------------------
    ! Description : build object phi(:,:,:,:). This is a rank four array.
    !               It is necessary for the construction of the electron-electron-nuclear backflow terms.
    !               The first index has dimension equal to the number of electrons.
    !               The second index has dimension equal to the number of electrons.
    !               The third index has dimensions equal to the number of nuclei.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 23 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin
    character(len=max_string_len_rout), save :: lhere = 'phi_bld'

    ! header
    if (header_exe) then

       call object_create('phi')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_phi_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed('a_param_phi')
       call object_needed('all_elec')

       return
    endif

    call object_alloc('phi', phi, nelec, nelec, ncent, nwalk)
    phi (:,:,:,:) = 0.d0

    !power series terms of phi as seen in three-body jastrow
    param_i = 1
    do order_p = 2, order_phi_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_j = 1, nup
                      do elec_i = 1, nup
                         phi(elec_i, elec_j, cent_i, :) = phi(elec_i, elec_j, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_up) &
                              & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, :)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l &
                              & +  scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :) &
                              & *  scaled_dist_een_en_wlk(elec_j, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_i = nup + 1, nelec
                         phi(elec_i, elec_j, cent_i, :) = phi(elec_i, elec_j, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, :)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l +  scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l )&
                              & * ( scaled_dist_een_en_wlk(elec_i, cent_i, :) *  scaled_dist_een_en_wlk(elec_j, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                   do elec_j = nup + 1, nelec
                      do elec_i = 1, nup
                         phi(elec_i, elec_j, cent_i, :) = phi(elec_i, elec_j, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, :)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l +  scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l )&
                              & * ( scaled_dist_een_en_wlk(elec_i, cent_i, :) *  scaled_dist_een_en_wlk(elec_j, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_i = nup + 1, nelec
                         phi(elec_i, elec_j, cent_i, :) = phi(elec_i, elec_j, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), down_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, :)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l +  scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l )&
                              & * ( scaled_dist_een_en_wlk(elec_i, cent_i, :) *  scaled_dist_een_en_wlk(elec_j, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                end do
                param_i = param_i + 1
             end if
          enddo
       enddo
    enddo

    if (all_elec) then
       !apply smooth cutoff function
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       do cent_i=1, ncent
          do elec_j=1, nelec
             phi(:,elec_j,cent_i,:) = phi(:,elec_j,cent_i,:) * product(smooth_cutoff_g(:,1:cent_i-1,:),2) * product(smooth_cutoff_g(:,cent_i+1:ncent,:),2)
          end do
       enddo
    end if

    !no self backflow
    do elec_i=1, nelec
       phi(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine phi_bld

  !=============================================================================================

  !==============================================================================================

  subroutine a_param_phi_bld

    !---------------------------------------------------------------------------
    ! Description : build a_param_phi. This is a rank 3 array.
    !               List of a parameters for phi backflow function for each center type, and each spin dependence
    !               The first index is which parameter.
    !               The second index is the type of center.
    !               The third index is the type of spin dependence.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: dep_a_param_i, ctype_i, spin_dependency_i, a_param_nb
    character(len=max_string_len_rout), save :: lhere = 'a_param_phi_bld'
    integer, parameter :: spin_dependencies_nb = 3    !because three choices up_up, down_down, up_down

    ! header
    if (header_exe) then

       call object_create('a_param_phi')

       call object_needed('order_phi_bf')
       call object_needed('nctype')
       call object_needed('phi_spin_dependence')
       call object_needed('read_a_param')
       call object_needed('read_a_param_phi_nb')
       call object_needed('dep_a_param_phi_parallel')
       call object_needed('dep_a_param_phi_anti_parallel')
       call object_needed('a_param_phi_parallel_cond')
       call object_needed('a_param_phi_anti_parallel_cond')
       call object_needed('nb_a_param_phi')

       return
    endif

    !need to have the option of using read_a_param_phi if it is read in
    if (read_a_param) then
       call object_provide_in_node(lhere, 'read_a_param_phi')
    end if

    call object_alloc('a_param_phi', a_param_phi, nb_a_param_phi, nctype, spin_dependencies_nb)
    a_param_phi = 0.d0

    !do we read in parameters?
    if (read_a_param) then
       if (phi_spin_dependence == spin_dependencies_nb) then
          a_param_nb = nb_a_param_phi * nctype * spin_dependencies_nb
          if (read_a_param_phi_nb == a_param_nb) then
             a_param_phi = reshape( read_a_param_phi, (/ nb_a_param_phi, nctype, spin_dependencies_nb /))
          else
             call require(lhere, 'number of parameters in read_a_param_phi equal to nb_a_param_phi * nctype * 3', .false.)
          end if
       elseif( phi_spin_dependence == spin_dependencies_nb - 1) then
          a_param_nb = nb_a_param_phi * nctype * (spin_dependencies_nb - 1)
          if (read_a_param_phi_nb == a_param_nb) then
             a_param_phi = reshape( read_a_param_phi, (/ nb_a_param_phi, nctype, spin_dependencies_nb /), pad=read_a_param_phi)
             !want up_up in a_param_phi(:,:,1), up_up=down_down in a_param_phi(:,:,2), up_down in a_param_phi(:,:,3)
             a_param_phi(:,:,3)=a_param_phi(:,:,2)
             a_param_phi(:,:,2)=a_param_phi(:,:,1)
          else
             call require(lhere, 'number of parameters in read_a_param_phi equal to nb_a_param_phi * nctype * 2', .false.)
          end if
       else
          a_param_nb = nb_a_param_phi * nctype
          if (read_a_param_phi_nb == a_param_nb) then
             a_param_phi = reshape( read_a_param_phi, (/ nb_a_param_phi, nctype, spin_dependencies_nb /), pad=read_a_param_phi)
          else
             call require(lhere, 'number of parameters in read_a_param_phi equal to nb_a_param_phi * nctype', .false.)
          end if
       end if
       !don't want to re-read
       read_a_param = .false.
    else
       !start from random guess
       call random_number(a_param_phi(:,:,:))
    endif

    !impose cusp
    !To impose the cusp we take the pivot positions in our reduced row echelon form of the matrix of cusp constraints
    !and we set this dependent parameter equal to the negative of the dot product of the rest of the row with
    !the corresponding a parameters. This is nothing fancy, but simply a standard way of solving linear equations

    !impose cusp parallel
    do spin_dependency_i=1, 2
       do ctype_i=1, nctype
          do dep_a_param_i=1, size(dep_a_param_phi_parallel,1)
             a_param_phi(dep_a_param_phi_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                     & dot_product( -a_param_phi_parallel_cond(dep_a_param_i, dep_a_param_phi_parallel(dep_a_param_i)+1:),&
                     & a_param_phi(dep_a_param_phi_parallel(dep_a_param_i) + 1:, ctype_i, spin_dependency_i))
          end do
       end do
    end do
    !impose cusp anti_parallel
    spin_dependency_i = 3 !up_down or down_up
    do ctype_i=1, nctype
       do dep_a_param_i=1, size(dep_a_param_phi_anti_parallel,1)
          a_param_phi(dep_a_param_phi_anti_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
               & dot_product( -a_param_phi_anti_parallel_cond(dep_a_param_i, dep_a_param_phi_anti_parallel(dep_a_param_i)+1:),&
               & a_param_phi(dep_a_param_phi_anti_parallel(dep_a_param_i)+1:, ctype_i, spin_dependency_i))
       end do
    end do

    if (phi_spin_dependence == 2) then
       !require up-up = down-down
       a_param_phi(:,:,2) = a_param_phi(:,:,1)
    elseif (phi_spin_dependence == 1) then
       !require up-up = down-down=up-down
       a_param_phi(:,:,2) = a_param_phi(:,:,1)
       a_param_phi(:,:,3) = a_param_phi(:,:,1)
    end if

  end subroutine a_param_phi_bld

  !==============================================================================================

  !==============================================================================================

  subroutine nb_a_param_phi_bld

    !---------------------------------------------------------------------------
    ! Description : build nb_a_param_phi. This is an integer.
    !               Total number of a parameters in phi per nucleus type, per spin dependence.
    !               To get total number of a parameters multiply by number of nuclei types and
    !               spin dependence.
    !
    ! Created     : F. Petruzielo, 12 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_p, order_k, order_l, order_l_max

    ! header
    if (header_exe) then

       call object_create('nb_a_param_phi')

       !need order of expansion
       call object_needed('order_phi_bf')

       return
    endif

    call object_associate('nb_a_param_phi',nb_a_param_phi)

    !determine the number of coefficients for the three body bf
    nb_a_param_phi = 0 !intialize parameter counter
    do order_p = 2, order_phi_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                nb_a_param_phi = nb_a_param_phi + 1
             endif
          enddo
       enddo
    enddo

  end subroutine nb_a_param_phi_bld

  !==============================================================================================

  !==============================================================================================

  subroutine dep_a_param_phi_parallel_bld

    !---------------------------------------------------------------------------
    ! Description : build dep_a_param_phi_parallel. This is a rank one array.
    !               This produces a list of the dependent parameters from a_param_phi_parallel_cond.
    !               To do this we find the pivots in a row reduced matrix.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer :: nb_dep_a_param_phi_parallel
    integer :: row_i, pivot_column_i, column_i

    ! header
    if (header_exe) then

       call object_create('dep_a_param_phi_parallel')

       call object_needed('nb_a_param_phi')
       call object_needed('a_param_phi_parallel_cond')

       return
    endif



    !initialize the number of dependent parameters
    nb_dep_a_param_phi_parallel = 0
    do row_i = 1, size(a_param_phi_parallel_cond,1)
       if (sum(abs(a_param_phi_parallel_cond(row_i,:))) < 1.d-10) then
          !no pivot
          cycle
       end if
       !pivot found, so increment counter
       nb_dep_a_param_phi_parallel = nb_dep_a_param_phi_parallel + 1
    end do

    call object_alloc('dep_a_param_phi_parallel', dep_a_param_phi_parallel, nb_dep_a_param_phi_parallel)

    !intialize the pivot column
    pivot_column_i = 1
    do row_i = 1, nb_dep_a_param_phi_parallel
       do column_i = pivot_column_i, nb_a_param_phi
          !first non zero element  (should be a 1)
          if (a_param_phi_parallel_cond(row_i,column_i) > 1.d-10) then
             !pivot found in this row
             dep_a_param_phi_parallel(row_i) = column_i
             !next pivot column must be to the right of this one
             pivot_column_i=column_i + 1
             !so move on to next row
             exit
          end if
       end do
    end do

  end subroutine dep_a_param_phi_parallel_bld

  !==============================================================================================

  !==============================================================================================

  subroutine dep_a_param_phi_anti_parallel_bld

    !---------------------------------------------------------------------------
    ! Description : build dep_a_param_phi_anti_parallel. This is a rank one array.
    !               This produces a list of the dependent parameters from a_param_phi_anti_parallel_cond.
    !               To do this we find the pivots in a row reduced matrix.
    !
    ! Created     : F. Petruzielo, 13 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer :: nb_dep_a_param_phi_anti_parallel
    integer :: row_i, pivot_column_i, column_i

    ! header
    if (header_exe) then

       call object_create('dep_a_param_phi_anti_parallel')

       call object_needed('nb_a_param_phi')
       call object_needed('a_param_phi_anti_parallel_cond')

       return
    endif

    !initialize the number of dependent parameters
    nb_dep_a_param_phi_anti_parallel = 0
    do row_i = 1, size(a_param_phi_anti_parallel_cond,1)
       if (sum(abs(a_param_phi_anti_parallel_cond(row_i,:))) < 1.d-10) then
          !no pivot
          cycle
       end if
       !pivot found, so increment counter
       nb_dep_a_param_phi_anti_parallel = nb_dep_a_param_phi_anti_parallel + 1
    end do

    call object_alloc('dep_a_param_phi_anti_parallel', dep_a_param_phi_anti_parallel, nb_dep_a_param_phi_anti_parallel)

    !intialize the pivot column
    pivot_column_i = 1
    do row_i = 1, nb_dep_a_param_phi_anti_parallel
       do column_i = pivot_column_i, nb_a_param_phi
          !first non zero element  (should be a 1)
          if (a_param_phi_anti_parallel_cond(row_i,column_i) > 1.d-10) then
             !pivot found in this row
             dep_a_param_phi_anti_parallel(row_i) = column_i
             !next pivot column must be to the right of this one
             pivot_column_i=column_i + 1
             !so move on to next row
             exit
          end if
       end do
    end do

  end subroutine dep_a_param_phi_anti_parallel_bld

  !==============================================================================================

  !==============================================================================================

  subroutine a_param_phi_parallel_cond_bld

    !---------------------------------------------------------------------------
    ! Description : build a_param_phi_parallel_cond. This is a rank two array.
    !               The first index is the number of condition ( some are possibly redundant).
    !               The second index is the number a parameters.
    !               We use gaussian elimination to row reduce the conditions and check for redundancies.
    !
    ! Created     : F. Petruzielo, 13 Aug 2008
    !---------------------------------------------------------------------------

    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_p, order_k, order_l, order_l_max, param_i

    ! header
    if (header_exe) then

       call object_create('a_param_phi_parallel_cond')

       call object_needed('order_phi_bf')
       call object_needed('nb_a_param_phi')
       call object_needed('all_elec')

       return
    endif


    if (all_elec) then
       ! two derivative constraints evaluated at 0, one constraint on function, and one derivative constraint evaluated at another point so number of contraints is 2*(order_phi_bf -1 ) + 2 * order_phi_bf
       call object_alloc('a_param_phi_parallel_cond', a_param_phi_parallel_cond, 4 * order_phi_bf - 2, nb_a_param_phi)
    else
       ! two derivative constraints so number of contraints is 2*(order_phi_bf -1 )
       call object_alloc('a_param_phi_parallel_cond', a_param_phi_parallel_cond, 2*(order_phi_bf - 1), nb_a_param_phi)
    end if

    a_param_phi_parallel_cond(:,:) = 0.d0

    param_i = 0 !intialize parameter counter
    !setup matrix of constraints
    do order_p = 2, order_phi_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                param_i = param_i + 1
                !=================================
                ! d phi  |
                ! -----  |        = 0
                ! d r_iI |
                !         r_iI=0
                !=================================
                a_param_phi_parallel_cond((order_p+order_k-order_l)/2, param_i) = a_param_phi_parallel_cond( (order_p+order_k-order_l)/2, param_i) + (order_p - order_k + order_l) / 2
                !note that when p=ord_phi_bf and k+l = p then (p+k+l)/2 = ord_phi_bf but (p-k-l)/2 = 0 so this power does not contribute. Want to avoid out of bounds array
                if (order_p-order_k-order_l .ne. 0) then
                   a_param_phi_parallel_cond((order_p+order_k+order_l)/2, param_i) = a_param_phi_parallel_cond( (order_p+order_k+order_l)/2, param_i) + ( order_p - order_k - order_l) / 2
                end if
                !=================================
                ! d phi  |
                ! -----  |        = 0
                ! d r_ij |
                !         r_ij=0
                !=================================
                !note that the row needs to be shifted down by order_phi_bf - 1 because those rows hold the above constraints
                !note that when p=ord_phi_bf and k=0 then (p-k) = order_phi_bf but k=0 terms don't contribute. Want to avoid out of bounds array
                if (order_k .ne. 0) then
                   a_param_phi_parallel_cond(order_phi_bf - 1 + order_p-order_k, param_i) = a_param_phi_parallel_cond(order_phi_bf - 1 + order_p-order_k, param_i) + order_k
                end if
                if (all_elec) then
                   !================================================
                   ! phi  |
                   !      |        = 0
                   !      |
                   !       r_iI=0
                   !===============================================
                   !note that the row needs to be shifted down by 2*(order_phi_bf - 1) because those rows hold the above constraints
                   a_param_phi_parallel_cond(2 *(order_phi_bf - 1) + (order_p+order_k-order_l)/2, param_i) = a_param_phi_parallel_cond(2 * (order_phi_bf -1 ) + (order_p+order_k-order_l)/2, param_i) + 1.d0
                   a_param_phi_parallel_cond(2 *(order_phi_bf - 1) + (order_p+order_k+order_l)/2, param_i) = a_param_phi_parallel_cond(2 * (order_phi_bf -1 ) + (order_p+order_k+order_l)/2, param_i) + 1.d0
                   !================================================
                   !d phi  |
                   ! ----  |        = 0
                   !d r_ij |
                   !       r_iI=0
                   !===============================================
                   !note that the row needs to be shifted down by 2*(order_phi_bf - 1) + order_phi_bf because those rows hold the above constraints
                   !further we shift by an additional 1 because there can be 0th order terms.
                   a_param_phi_parallel_cond(3*order_phi_bf-2+(order_p+order_k-2-order_l)/2 + 1, param_i) = a_param_phi_parallel_cond(3*order_phi_bf-2+(order_p+order_k-2-order_l)/2+1, param_i) + order_k
                   a_param_phi_parallel_cond(3*order_phi_bf-2+(order_p+order_k-2+order_l)/2 + 1, param_i) = a_param_phi_parallel_cond(3*order_phi_bf-2+(order_p+order_k-2+order_l)/2+1, param_i) + order_k
                end if
             endif
          enddo
       enddo
    enddo

    !use gaussian elimintation on constraints to remove redundancies
    call row_reduce(a_param_phi_parallel_cond)

  end subroutine a_param_phi_parallel_cond_bld

  !==============================================================================================

  !==============================================================================================

  subroutine a_param_phi_anti_parallel_cond_bld

    !---------------------------------------------------------------------------
    ! Description : build a_param_phi_anti_parallel_spin_cond_bld. This is a rank two array.
    !               The first index is the number of condition ( some are possibly redundant).
    !               The second index is the number a parameters.
    !               We use gaussian elimination to row reduce the conditions and check for redundancies.
    !
    ! Created     : F. Petruzielo, 13 Aug 2008
    !---------------------------------------------------------------------------

    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_p, order_k, order_l, order_l_max, param_i

    ! header
    if (header_exe) then

       call object_create('a_param_phi_anti_parallel_cond')

       call object_needed('order_phi_bf')
       call object_needed('nb_a_param_phi')
       call object_needed('all_elec')

       return
    endif


    if (all_elec) then
       ! one derivative constraints evaluated at 0, one constraint on function, and one derivative constraint evaluated at another point so number of contraints is (order_phi_bf -1 ) + 2*order_phi_bf
       call object_alloc('a_param_phi_anti_parallel_cond', a_param_phi_anti_parallel_cond, 3 * order_phi_bf - 1, nb_a_param_phi)
    else
       ! one derivative constraints so number of contraints is (order_phi_bf -1 )
       call object_alloc('a_param_phi_anti_parallel_cond', a_param_phi_anti_parallel_cond, order_phi_bf - 1, nb_a_param_phi)
    end if

    a_param_phi_anti_parallel_cond(:,:) = 0.d0

    param_i = 0 !intialize parameter counter
    !setup matrix of constraints
    do order_p = 2, order_phi_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                param_i = param_i + 1
                !=================================
                ! d phi  |
                ! -----  |        = 0
                ! d r_iI |
                !         r_iI=0
                !=================================
                a_param_phi_anti_parallel_cond((order_p+order_k-order_l)/2, param_i) = a_param_phi_anti_parallel_cond( (order_p+order_k-order_l)/2, param_i) + (order_p - order_k + order_l) / 2
                !note that when p=ord_phi_bf and k+l = p then (p+k+l)/2 = ord_phi_bf but (p-k-l)/2 = 0 so this power does not contribute. Want to avoid out of bounds array
                if ((order_p + order_k + order_l)/2 .ne. order_phi_bf) then
                   a_param_phi_anti_parallel_cond((order_p+order_k+order_l)/2, param_i) = a_param_phi_anti_parallel_cond( (order_p+order_k+order_l)/2, param_i) + ( order_p - order_k - order_l) / 2
                end if
                if (all_elec) then
                   !================================================
                   ! phi  |
                   !      |        = 0
                   !      |
                   !       r_iI=0
                   !===============================================
                   !note that the row needs to be shifted down by (order_phi_bf - 1) because those rows hold the above constraints
                   a_param_phi_anti_parallel_cond((order_phi_bf - 1) + (order_p+order_k-order_l)/2, param_i) = a_param_phi_anti_parallel_cond((order_phi_bf -1 ) + (order_p+order_k-order_l)/2, param_i) + 1.d0
                   a_param_phi_anti_parallel_cond((order_phi_bf - 1) + (order_p+order_k+order_l)/2, param_i) = a_param_phi_anti_parallel_cond((order_phi_bf -1 ) + (order_p+order_k+order_l)/2, param_i) + 1.d0
                   !================================================
                   !d phi  |
                   ! ----  |        = 0
                   !d r_ij |
                   !       r_iI=0
                   !===============================================
                   !note that the row needs to be shifted down by (order_phi_bf - 1) + order_phi_bf because those rows hold the above constraints
                   a_param_phi_anti_parallel_cond(2*order_phi_bf-1+(order_p+order_k-2-order_l)/2+1, param_i) = a_param_phi_anti_parallel_cond(2*order_phi_bf-1+(order_p+order_k-2-order_l)/2+1, param_i) + order_k
                   a_param_phi_anti_parallel_cond(2*order_phi_bf-1+(order_p+order_k-2+order_l)/2+1, param_i) = a_param_phi_anti_parallel_cond(2*order_phi_bf-1+(order_p+order_k-2+order_l)/2+1, param_i) + order_k
                end if
             endif
          enddo
       enddo
    enddo

    !use gaussian elimintation on constraints to remove redundancies
    call row_reduce(a_param_phi_anti_parallel_cond)

  end subroutine a_param_phi_anti_parallel_cond_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_xi_een_phi_j_beta_d_r_i_alpha_bld

    !---------------------------------------------------------------------------
    ! Description : build d_xi_een_phi_j_beta_d_r_i_alpha. This is a rank five array.
    !               Derivative of the phi part of electron-electron-nuclear backflow function
    !               with respect to electron position.
    !               It is necessary for the construction of the energy, drift velocity.
    !               The first and second indices are the number of spatial dimensions.
    !               The first index is for the components of xi_een_phi (beta)
    !               The second index is for the components of the gradient (alpha)
    !               The third and fourth indices have dimensions equal to the number of electrons.
    !               The third index is the label of the een phi backflow transformation of each electon
    !               The fourth index is for the electron that the gradient is taken with respect to
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim, ncent
    include 'commons.h'

    !local
    integer  :: elec_j, dim_i, dim_j

    ! header
    if (header_exe) then

       call object_create('d_xi_een_phi_j_beta_d_r_i_alpha')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('ndim')
       call object_needed('d_phi_d_r_j_alpha_first')
       call object_needed('d_phi_d_r_i_alpha_second')
       call object_needed('phi')

       return
    endif

    call object_alloc('d_xi_een_phi_j_beta_d_r_i_alpha', d_xi_een_phi_j_beta_d_r_i_alpha, ndim, ndim, nelec, nelec, nwalk)
    d_xi_een_phi_j_beta_d_r_i_alpha(:,:,:,:,:) = 0.d0

    do elec_j = 1, nelec
       do dim_i = 1, ndim     !alpha
          !delta_{alpha,beta}
          d_xi_een_phi_j_beta_d_r_i_alpha(dim_i, dim_i, elec_j, elec_j, :) = sum(sum(phi(elec_j, :, :, :),1),1)
          do dim_j = 1, ndim   !beta
             d_xi_een_phi_j_beta_d_r_i_alpha(dim_j, dim_i, elec_j, elec_j, :) = d_xi_een_phi_j_beta_d_r_i_alpha(dim_j, dim_i, elec_j, elec_j, :) + sum(sum(d_phi_d_r_j_alpha_first(dim_i, elec_j, :, :, :),2) * vec_ee_xyz_wlk(dim_j, elec_j, :, :), 1)
          end do
       end do
    end do

    do dim_i = 1, ndim
       d_xi_een_phi_j_beta_d_r_i_alpha(dim_i, dim_i, :, :, :) = d_xi_een_phi_j_beta_d_r_i_alpha(dim_i, dim_i, :, :, :) - sum(phi(:, :, :, :),3)
       do dim_j = 1, ndim
          d_xi_een_phi_j_beta_d_r_i_alpha(dim_j, dim_i, :, :, :) = d_xi_een_phi_j_beta_d_r_i_alpha(dim_j, dim_i, :, :, :) + sum(d_phi_d_r_i_alpha_second(dim_i, :, :, :, :),3) * vec_ee_xyz_wlk(dim_j, :, :, :)
       end do
    end do

  end subroutine d_xi_een_phi_j_beta_d_r_i_alpha_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_xi_een_phi_j_beta_d_r_i_alpha_num_bld
    !---------------------------------------------------------------------------
    ! Description : build d_xi_een_phi_j_beta_d_r_i_alpha__num. This is a rank five array.
    !               Numerical first derivative of the xi_een_phi:
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               The first and second indices are the number of spatial dimensions.
    !               The first index is for the components of xi_een_phi (beta)
    !               The second index is for the components of the gradient (alpha)
    !               The third and fourth indices have dimensions equal to the number of electrons.
    !               The third index is the label of the een_phi backflow transformation of each electon
    !               The fourth index is for the electron that the gradient is taken with respect to
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 14 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, dim_i, dim_j
    real(dp), allocatable :: xi_een_phi_temp_plus(:,:,:), xi_een_phi_temp_minus(:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'd_xi_een_phi_j_beta_d_r_i_alpha_num_bld'

    ! header
    if (header_exe) then

       call object_create('d_xi_een_phi_j_beta_d_r_i_alpha_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('xi_een_phi')
       call object_needed('coord_elec_wlk')

       return
    endif

    ! allocation
    call object_alloc('d_xi_een_phi_j_beta_d_r_i_alpha_num', d_xi_een_phi_j_beta_d_r_i_alpha_num, ndim, ndim, nelec, nelec, nwalk)
    d_xi_een_phi_j_beta_d_r_i_alpha_num(:,:,:,:,:) = 0.d0

    !allocate temporary objects
    allocate(xi_een_phi_temp_plus(ndim, nelec, nwalk))
    allocate(xi_een_phi_temp_minus(ndim, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

    do elec_j = 1, nelec
       do elec_i = 1, nelec
          do dim_j = 1, ndim  !beta
             do dim_i = 1, ndim  !alph
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) + epsilon
                call object_modified('coord_elec_wlk')
                !get updated xi_een_phi
                call object_provide_in_node(lhere, 'xi_een_phi')
                xi_een_phi_temp_plus = xi_een_phi
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) - epsilon
                call object_modified('coord_elec_wlk')
                !get updated xi_een_phi
                call object_provide_in_node(lhere, 'xi_een_phi')
                xi_een_phi_temp_minus = xi_een_phi
                !calculate derivative
                d_xi_een_phi_j_beta_d_r_i_alpha_num(dim_j, dim_i, elec_j, elec_i, :) = (xi_een_phi_temp_plus(dim_j, elec_j, :) - xi_een_phi_temp_minus(dim_j, elec_j, :) ) / (2.d0 * epsilon)
             end do
          end do
       end do
    end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'xi_een_phi')

  end subroutine d_xi_een_phi_j_beta_d_r_i_alpha_num_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_phi_d_r_i_alpha_second_bld

    !---------------------------------------------------------------------------
    ! Description : build d_phi_d_r_i_alpha_second. This is a rank five array.
    !               First derivative of the phi:
    !
    !               d phi(r_jI,r_iI,r_ji)
    !               ----------
    !               d r_i_alpha
    !
    !               Note alpha refers to cartesian coords. In 3-d, r_alpha takes on values x,y,z.
    !               The first index is the over spatial dimension. This is the components of the gradient.
    !               The second index has dimensions equal to the number of electrons and labels the jth electron
    !               The third index has dimensions equal to the number of electrons and labels the ith electron
    !               The fourth index has dimension equal to the number of centers.
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim, ncent
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, dim_i
    character(len=max_string_len_rout), save :: lhere = 'd_phi_d_r_i_alpha_second_bld'


    ! header
    if (header_exe) then

       call object_create('d_phi_d_r_i_alpha_second')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed ('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed ('dist_ee_wlk')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('d_phi_d_r_j_inuc')
       call object_needed('d_phi_d_r_ji')
       call object_needed('all_elec')

       return
    endif

    call object_alloc('d_phi_d_r_i_alpha_second', d_phi_d_r_i_alpha_second, ndim, nelec, nelec, ncent, nwalk)
    d_phi_d_r_i_alpha_second(:,:,:,:,:) = 0.d0


    do elec_j = 1, nelec
       do dim_i = 1, ndim
          d_phi_d_r_i_alpha_second(dim_i, elec_j, :, :, :) = d_phi_d_r_j_inuc(:, elec_j, :, :) * vec_en_xyz_wlk(dim_i, :, :, :) / dist_en_wlk(:, :, :)
       end do
    end do

    do cent_i = 1, ncent
       do elec_j = 1, nelec
          do dim_i = 1, ndim
             d_phi_d_r_i_alpha_second(dim_i, elec_j, 1:elec_j-1, cent_i, :) = d_phi_d_r_i_alpha_second(dim_i, elec_j, 1:elec_j-1, cent_i, :) - d_phi_d_r_ji(elec_j, 1:elec_j-1, cent_i, :) * vec_ee_xyz_wlk(dim_i, elec_j, 1:elec_j-1, :) / dist_ee_wlk(elec_j, 1:elec_j-1, :)
             d_phi_d_r_i_alpha_second(dim_i, elec_j, elec_j+1:nelec, cent_i, :) = d_phi_d_r_i_alpha_second(dim_i, elec_j, elec_j+1:nelec, cent_i, :) - d_phi_d_r_ji(elec_j, elec_j+1:nelec, cent_i, :) * vec_ee_xyz_wlk(dim_i, elec_j, elec_j+1:nelec, :) / dist_ee_wlk(elec_j, elec_j+1:nelec, :)
          end do
       end do
    end do

    if (all_elec) then
       !apply smooth cutoff function
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       do cent_i=1, ncent
          do elec_i=1, nelec
             do dim_i = 1, ndim
                d_phi_d_r_i_alpha_second(dim_i, :, elec_i, cent_i, :) = d_phi_d_r_i_alpha_second(dim_i, :, elec_i, cent_i, :) * product(smooth_cutoff_g(:,1:cent_i-1,:),2) * product(smooth_cutoff_g(:,cent_i+1:ncent,:),2)
             end do
          enddo
       end do
    end if

  end subroutine d_phi_d_r_i_alpha_second_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_phi_d_r_j_alpha_first_bld

    !---------------------------------------------------------------------------
    ! Description : build d_phi_d_r_j_alpha_first. This is a rank five array.
    !               First derivative of the phi:
    !
    !               d phi(r_jI,r_lI,r_jl)
    !               ----------
    !               d r_j_alpha
    !
    !               Note alpha refers to cartesian coords. In 3-d, r_alpha takes on values x,y,z.
    !               The first index is the over spatial dimension. This is the components of the gradient.
    !               The second index has dimensions equal to the number of electrons and labels the jth electron
    !               The third index has dimensions equal to the number of electrons and labels the lth electron
    !               The fourth index has dimension equal to the number of centers.
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim, ncent
    include 'commons.h'

    !local
    integer  :: elec_j, elec_l, cent_i, dim_i
    character(len=max_string_len_rout), save :: lhere = 'd_phi_d_r_j_alpha_first_bld'

    ! header
    if (header_exe) then

       call object_create('d_phi_d_r_j_alpha_first')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed ('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed ('dist_ee_wlk')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('d_phi_d_r_j_inuc')
       call object_needed('d_phi_d_r_ji')
       call object_needed('all_elec')

       return
    endif

    call object_alloc('d_phi_d_r_j_alpha_first', d_phi_d_r_j_alpha_first, ndim, nelec, nelec, ncent, nwalk)
    d_phi_d_r_j_alpha_first(:,:,:,:,:) = 0.d0

    do elec_l = 1, nelec
       do dim_i = 1, ndim
          d_phi_d_r_j_alpha_first(dim_i, :, elec_l, :, :) = d_phi_d_r_j_inuc(:, elec_l, :, :) * vec_en_xyz_wlk(dim_i, :, :, :) / dist_en_wlk(:, :, :)
       end do
    end do

    do cent_i = 1, ncent
       do elec_j = 1, nelec
          do dim_i = 1, ndim
             d_phi_d_r_j_alpha_first(dim_i, elec_j, 1:elec_j-1, cent_i, :) = d_phi_d_r_j_alpha_first(dim_i, elec_j, 1:elec_j-1, cent_i, :) + d_phi_d_r_ji(elec_j, 1:elec_j-1, cent_i, :) * vec_ee_xyz_wlk(dim_i, elec_j, 1:elec_j-1, :) / dist_ee_wlk(elec_j, 1:elec_j-1, :)
             d_phi_d_r_j_alpha_first(dim_i, elec_j, elec_j+1:nelec, cent_i, :) = d_phi_d_r_j_alpha_first(dim_i, elec_j, elec_j+1:nelec, cent_i, :) + d_phi_d_r_ji(elec_j, elec_j+1:nelec, cent_i, :) * vec_ee_xyz_wlk(dim_i, elec_j, elec_j+1:nelec, :) / dist_ee_wlk(elec_j, elec_j+1:nelec, :)
          end do
       end do
    end do

    if (all_elec) then
       !apply smooth cutoff function
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       call object_provide_in_node(lhere, 'd_smooth_cutoff_g_d_r')
       call object_provide_in_node(lhere, 'phi')
       do cent_i=1, ncent
          do elec_l=1, nelec
             do dim_i = 1, ndim
                d_phi_d_r_j_alpha_first(dim_i, :, elec_l, cent_i, :) = d_phi_d_r_j_alpha_first(dim_i, :, elec_l, cent_i, :) * product(smooth_cutoff_g(:,1:cent_i-1,:),2) * product(smooth_cutoff_g(:,cent_i+1:ncent,:),2) + phi(:, elec_l, cent_i, :) * ( sum(vec_en_xyz_wlk(dim_i, :, 1:cent_i-1, :) / dist_en_wlk(:, 1:cent_i-1, :) / smooth_cutoff_g(:, 1:cent_i-1, :) * d_smooth_cutoff_g_d_r(:, 1:cent_i-1, :),2) + sum(vec_en_xyz_wlk(dim_i, :, cent_i+1:ncent, :)/ dist_en_wlk(:, cent_i+1:ncent, :) / smooth_cutoff_g(:, cent_i+1:ncent, :) * d_smooth_cutoff_g_d_r(:, cent_i+1:ncent, :),2))
             end do
          enddo
       end do
    end if

  end subroutine d_phi_d_r_j_alpha_first_bld

  !==============================================================================================

  !======================================================================================

  subroutine d_phi_d_r_ji_bld
    !---------------------------------------------------------------------------
    ! Description : build object d_phi_d_r_ji(:,:,:,:). This is a rank four array.
    !
    !               d phi(r_jI,r_iI,r_ji)    (without cutoff)
    !               ----------
    !               d r_ji
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d_phi_d_r_ji')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_phi_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_ee_wlk_d_r')
       call object_needed('a_param_phi')

       return
    endif

    call object_alloc('d_phi_d_r_ji', d_phi_d_r_ji, nelec, nelec, ncent, nwalk)
    d_phi_d_r_ji (:,:,:,:) = 0.d0


    param_i = 1 !intialize parameter counter
    do order_p = 2, order_phi_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d_phi_d_r_ji(elec_j, elec_i, cent_i, :) = d_phi_d_r_ji(elec_j, elec_i, cent_i, :) + order_k * a_param_phi(param_i, iwctype(cent_i), up_up)&
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) ) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l +  &
                              & scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_j = nup + 1, nelec
                         d_phi_d_r_ji(elec_j, elec_i, cent_i, :) = d_phi_d_r_ji(elec_j, elec_i, cent_i, :) + order_k * a_param_phi(param_i, iwctype(cent_i), up_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) ) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l +  &
                              & scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d_phi_d_r_ji(elec_j, elec_i, cent_i, :) = d_phi_d_r_ji(elec_j, elec_i, cent_i, :) + order_k * a_param_phi(param_i, iwctype(cent_i), up_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) ) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l +  &
                              & scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_j =  nup + 1, nelec
                         d_phi_d_r_ji(elec_j, elec_i, cent_i, :) = d_phi_d_r_ji(elec_j, elec_i, cent_i, :) + order_k * a_param_phi(param_i, iwctype(cent_i), down_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) ) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l +  &
                              & scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d_phi_d_r_ji(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d_phi_d_r_ji_bld

  !=============================================================================================

  !======================================================================================

  subroutine d_phi_d_r_j_inuc_bld
    !---------------------------------------------------------------------------
    ! Description : build object d_phi_d_r_j_inuc(:,:,:,:). This is a rank four array.
    !
    !               d phi(r_jI,r_iI,r_jI)    (without cutoff)
    !               ----------
    !               d r_jI
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d_phi_d_r_j_inuc')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_phi_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_en_wlk_d_r')
       call object_needed('a_param_phi')

       return
    endif

    call object_alloc('d_phi_d_r_j_inuc', d_phi_d_r_j_inuc, nelec, nelec, ncent, nwalk)
    d_phi_d_r_j_inuc (:,:,:,:) = 0.d0

    param_i = 1 !intialize parameter counter
    do order_p = 2, order_phi_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d_phi_d_r_j_inuc(elec_j, elec_i, cent_i, :) = d_phi_d_r_j_inuc(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_up) &
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k ) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) &
                              & * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) &
                              & + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                      do elec_j = nup+1, nelec
                         d_phi_d_r_j_inuc(elec_j, elec_i, cent_i, :) = d_phi_d_r_j_inuc(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_down) &
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k ) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) &
                              & * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) &
                              & + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d_phi_d_r_j_inuc(elec_j, elec_i, cent_i, :) = d_phi_d_r_j_inuc(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_down) &
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k ) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) &
                              & * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) &
                              & + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                      do elec_j = nup+1, nelec
                         d_phi_d_r_j_inuc(elec_j, elec_i, cent_i, :) = d_phi_d_r_j_inuc(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), down_down) &
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k ) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) &
                              & * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) &
                              & + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                   end do

                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d_phi_d_r_j_inuc(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d_phi_d_r_j_inuc_bld

  !=============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_een_phi_j_beta_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_een_phi_j_beta. This is a rank four array.
    !               Laplacian of the electron-electron-nuclear phi backflow function
    !               with respect to electron position.
    !               It is necessary for the construction of the energy.
    !               The first index is for the components of xi_een_phi (beta)
    !               The second index is the label of the ee backflow transformed electon
    !               The third index is for the electron that the laplacian is taken with respect to
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 31 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j, dim_j

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_een_phi_j_beta')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('ndim')
       call object_needed('d_phi_d_r_j_alpha_first')
       call object_needed('d_phi_d_r_i_alpha_second')
       call object_needed('lap_j_phi_first')
       call object_needed('lap_i_phi_second')

       return
    endif

    call object_alloc('lap_i_xi_een_phi_j_beta', lap_i_xi_een_phi_j_beta, ndim, nelec, nelec, nwalk)
    lap_i_xi_een_phi_j_beta(:,:,:,:) = 0.d0

    do elec_j = 1, nelec
       do dim_j = 1, ndim
          lap_i_xi_een_phi_j_beta(dim_j, elec_j, elec_j, :) =  sum(sum(lap_j_phi_first(elec_j,:,:,:), 2) * vec_ee_xyz_wlk(dim_j, elec_j, :, :) + 2.d0  * sum(d_phi_d_r_j_alpha_first(dim_j, elec_j, :, :, :), 2), 1)
       end do
    end do

    do dim_j = 1,ndim
       lap_i_xi_een_phi_j_beta(dim_j, :, :, :) = lap_i_xi_een_phi_j_beta(dim_j, :, :, :) + sum(lap_i_phi_second(:, :, :, :), 3) * vec_ee_xyz_wlk(dim_j, :, :, :) - 2.d0 * sum(d_phi_d_r_i_alpha_second(dim_j, :, :, :, :), 3)
    end do

  end subroutine lap_i_xi_een_phi_j_beta_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_een_phi_j_beta_num_bld
    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_een_phi_j_beta_num. This is a rank four array.
    !               Numerical laplacian of xi_een_phi:
    !
    !               (nabla_i)^2 xi_een_phi_j^beta
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               The first index is for the components of xi_een_phi (beta)
    !               The second index is the label of the een_phi backflow transformation of each electon
    !               The third index is for the electron that the laplacian is taken with respect to
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 31 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, dim_i, dim_j
    real(dp), allocatable :: d_xi_een_phi_j_beta_d_r_i_alpha_plus(:,:,:,:,:), d_xi_een_phi_j_beta_d_r_i_alpha_minus(:,:,:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'lap_i_xi_een_phi_j_beta_num_bld'

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_een_phi_j_beta_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('coord_elec_wlk')
       call object_needed('d_xi_een_phi_j_beta_d_r_i_alpha')

       return
    endif

    ! allocation
    call object_alloc('lap_i_xi_een_phi_j_beta_num', lap_i_xi_een_phi_j_beta_num, ndim, nelec, nelec, nwalk)
    lap_i_xi_een_phi_j_beta_num(:,:,:,:) = 0.d0

    !allocate temporary objects
    allocate(d_xi_een_phi_j_beta_d_r_i_alpha_plus(ndim, ndim, nelec, nelec, nwalk))
    allocate(d_xi_een_phi_j_beta_d_r_i_alpha_minus(ndim, ndim, nelec, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

     do elec_i = 1, nelec
        do elec_j = 1, nelec
           do dim_j = 1, ndim  !beta
              do dim_i = 1, ndim  !alpha
                 coord_elec_wlk = coord_elec_wlk_temp
                 !shift coordinates in a particular dimension, x,y,z, etc
                 coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) + epsilon
                 call object_modified('coord_elec_wlk')
                 !get updated d_xi_een_phi_j_beta_d_r_i_alpha
                 call object_provide_in_node(lhere, 'd_xi_een_phi_j_beta_d_r_i_alpha')
                 d_xi_een_phi_j_beta_d_r_i_alpha_plus = d_xi_een_phi_j_beta_d_r_i_alpha
                 coord_elec_wlk = coord_elec_wlk_temp
                 !shift coordinates in a particular dimension, x,y,z, etc
                 coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) - epsilon
                 call object_modified('coord_elec_wlk')
                 !get updated d_xi_een_phi_j_beta_d_r_i_alpha
                 call object_provide_in_node(lhere, 'd_xi_een_phi_j_beta_d_r_i_alpha')
                 d_xi_een_phi_j_beta_d_r_i_alpha_minus = d_xi_een_phi_j_beta_d_r_i_alpha
                 !calculate derivative
                 lap_i_xi_een_phi_j_beta_num(dim_j, elec_j, elec_i, :) = lap_i_xi_een_phi_j_beta_num(dim_j, elec_j, elec_i, :) + (d_xi_een_phi_j_beta_d_r_i_alpha_plus(dim_j, dim_i, elec_j, elec_i, :) - d_xi_een_phi_j_beta_d_r_i_alpha_minus(dim_j, dim_i, elec_j, elec_i, :) ) / (2.d0 * epsilon)
              end do
           end do
        end do
     end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'd_xi_een_phi_j_beta_d_r_i_alpha')

  end subroutine lap_i_xi_een_phi_j_beta_num_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_j_phi_first_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_j_phi_first. This is a rank four array.
    !               Laplacian of the phi:
    !
    !               (nabla_j)^2 phi(r_jI,r_lI,r_jl)       Note that this is derivative with respect to first index so derivatives of the cutoff
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of electrons and labels the lth electron
    !               The third index has dimensions equal to the number of centers and labels the Ith nucleus
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 31 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_l, dim_i, cent_i
    character(len=max_string_len_rout), save :: lhere = 'lap_j_phi_first_bld'

    ! header
    if (header_exe) then

       call object_create('lap_j_phi_first')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed('ndim')
       call object_needed('all_elec')
       call object_needed('d_phi_d_r_j_inuc')
       call object_needed('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('d_phi_d_r_ji')
       call object_needed('dist_ee_wlk')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('d2_phi_d_r_ji_2')
       call object_needed('d2_phi_d_r_j_inuc_2')
       call object_needed('d2_phi_d_r_j_inuc_d_r_ji')

       return
    endif

    call object_alloc('lap_j_phi_first', lap_j_phi_first, nelec, nelec, ncent, nwalk)
    lap_j_phi_first(:,:,:,:) = 0.d0

    do elec_l = 1, nelec
       lap_j_phi_first(:, elec_l, :, :) = d2_phi_d_r_j_inuc_2(:, elec_l, :, :) + (ndim - 1.d0) * d_phi_d_r_j_inuc(:, elec_l, :, :) / dist_en_wlk(:, :, :) + d2_phi_d_r_ji_2(:, elec_l, :, :)
    end do

    do cent_i = 1, ncent
       do elec_l = 1, nelec
          lap_j_phi_first(1:elec_l - 1, elec_l, cent_i, :) = lap_j_phi_first(1:elec_l - 1, elec_l, cent_i, :) + (ndim - 1.d0) * d_phi_d_r_ji(1:elec_l - 1, elec_l, cent_i, :) / dist_ee_wlk(1:elec_l - 1, elec_l, :) +  2.d0 * d2_phi_d_r_j_inuc_d_r_ji(1:elec_l - 1, elec_l, cent_i, :) / dist_en_wlk(1:elec_l - 1, cent_i, :) / dist_ee_wlk(1:elec_l - 1, elec_l, :) * sum(vec_en_xyz_wlk(:, 1:elec_l - 1, cent_i, :) * vec_ee_xyz_wlk(:, 1:elec_l - 1, elec_l, :), 1)
          lap_j_phi_first(elec_l + 1:nelec, elec_l, cent_i, :) = lap_j_phi_first(elec_l + 1:nelec, elec_l, cent_i, :) + (ndim - 1.d0) * d_phi_d_r_ji(elec_l + 1:nelec, elec_l, cent_i, :) / dist_ee_wlk(elec_l + 1:nelec, elec_l, :) + 2.d0 * d2_phi_d_r_j_inuc_d_r_ji(elec_l + 1:nelec, elec_l, cent_i, :) / dist_en_wlk(elec_l + 1:nelec, cent_i, :) / dist_ee_wlk(elec_l + 1:nelec, elec_l, :) * sum(vec_en_xyz_wlk(:, elec_l + 1:nelec, cent_i, :) * vec_ee_xyz_wlk(:, elec_l + 1:nelec, elec_l, :), 1)
       end do
    end do

    if (all_elec) then
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       call object_provide_in_node(lhere, 'd_smooth_cutoff_g_d_r')
       call object_provide_in_node(lhere, 'd2_smooth_cutoff_g_d_r_2')
       call object_provide_in_node(lhere, 'phi')
       call object_provide_in_node(lhere, 'd_phi_d_r_j_alpha_first')

       do cent_i = 1, ncent
          do elec_l = 1, nelec
             do dim_i = 1, ndim
                lap_j_phi_first(1:elec_l - 1, elec_l, cent_i, :) = lap_j_phi_first(1:elec_l - 1, elec_l, cent_i, :) + d_phi_d_r_ji(1:elec_l - 1, elec_l, cent_i, :) / dist_ee_wlk(1:elec_l - 1, elec_l, :) * vec_ee_xyz_wlk(dim_i, 1:elec_l - 1, elec_l, :) * ( sum(d_smooth_cutoff_g_d_r(1:elec_l - 1, 1:cent_i - 1, :) / smooth_cutoff_g(1:elec_l - 1, 1:cent_i - 1, :)  * vec_en_xyz_wlk(dim_i, 1:elec_l - 1, 1:cent_i - 1, :) / dist_en_wlk(1:elec_l - 1, 1:cent_i - 1, :), 2) + sum(d_smooth_cutoff_g_d_r(1:elec_l - 1, cent_i + 1:ncent, :) / smooth_cutoff_g(1:elec_l - 1, cent_i + 1:ncent, :)  * vec_en_xyz_wlk(dim_i, 1:elec_l - 1, cent_i + 1:ncent, :) / dist_en_wlk(1:elec_l - 1, cent_i + 1:ncent, :), 2) )
                lap_j_phi_first(elec_l + 1:nelec, elec_l, cent_i, :) = lap_j_phi_first(elec_l + 1:nelec, elec_l, cent_i, :) + d_phi_d_r_ji(elec_l + 1:nelec, elec_l, cent_i, :) / dist_ee_wlk(elec_l + 1:nelec, elec_l, :) * vec_ee_xyz_wlk(dim_i, elec_l + 1:nelec, elec_l, :) * ( sum(d_smooth_cutoff_g_d_r(elec_l + 1:nelec, 1:cent_i - 1, :) / smooth_cutoff_g(elec_l + 1:nelec, 1:cent_i - 1, :)  * vec_en_xyz_wlk(dim_i, elec_l + 1:nelec, 1:cent_i - 1, :) / dist_en_wlk(elec_l + 1:nelec, 1:cent_i - 1, :), 2) + sum(d_smooth_cutoff_g_d_r(elec_l + 1:nelec, cent_i + 1:ncent, :) / smooth_cutoff_g(elec_l + 1:nelec, cent_i + 1:ncent, :)  * vec_en_xyz_wlk(dim_i, elec_l + 1:nelec, cent_i + 1:ncent, :) / dist_en_wlk(elec_l + 1:nelec, cent_i + 1:ncent, :), 2) )
                lap_j_phi_first(:, elec_l, cent_i, :) = lap_j_phi_first(:, elec_l, cent_i, :) + d_phi_d_r_j_inuc(:, elec_l, cent_i, :) / dist_en_wlk(:, cent_i, :) * vec_en_xyz_wlk(dim_i, :, cent_i, :) * ( sum(d_smooth_cutoff_g_d_r(:, 1:cent_i - 1, :) / smooth_cutoff_g(:, 1:cent_i - 1, :)  * vec_en_xyz_wlk(dim_i, :, 1:cent_i - 1, :) / dist_en_wlk(:, 1:cent_i - 1, :), 2) + sum(d_smooth_cutoff_g_d_r(:, cent_i + 1:ncent, :) / smooth_cutoff_g(:, cent_i + 1:ncent, :)  * vec_en_xyz_wlk(dim_i, :, cent_i + 1:ncent, :) / dist_en_wlk(:, cent_i + 1:ncent, :), 2) )
             end do
          end do
       end do

       do cent_i = 1, ncent
          do elec_l = 1, nelec
             lap_j_phi_first(:, elec_l, cent_i, :) = lap_j_phi_first(:, elec_l, cent_i, :) * product(smooth_cutoff_g(:, 1:cent_i - 1, :), 2) * product(smooth_cutoff_g(:, cent_i + 1:ncent, :), 2)
          end do
       end do

       do cent_i = 1, ncent
          do elec_l = 1, nelec
             lap_j_phi_first(:, elec_l, cent_i, :) = lap_j_phi_first(:, elec_l, cent_i, :) + phi(:, elec_l, cent_i, :) * ( sum( (d2_smooth_cutoff_g_d_r_2(:, 1:cent_i - 1, :) - d_smooth_cutoff_g_d_r(:, 1:cent_i - 1, :) ** 2 / smooth_cutoff_g(:, 1:cent_i - 1, :) + d_smooth_cutoff_g_d_r(:, 1:cent_i - 1, :) * (ndim - 1.d0) / dist_en_wlk(:, 1:cent_i - 1, :) ) / smooth_cutoff_g(:, 1:cent_i - 1, :), 2) + sum( (d2_smooth_cutoff_g_d_r_2(:, cent_i + 1:ncent, :) - d_smooth_cutoff_g_d_r(:, cent_i + 1:ncent, :) ** 2 / smooth_cutoff_g(:, cent_i + 1:ncent, :) + d_smooth_cutoff_g_d_r(:, cent_i + 1:ncent, :) * (ndim - 1.d0) / dist_en_wlk(:, cent_i + 1:ncent, :) ) / smooth_cutoff_g(:, cent_i + 1:ncent, :), 2) )
             do dim_i = 1, ndim
                lap_j_phi_first(:, elec_l, cent_i, :) = lap_j_phi_first(:, elec_l, cent_i, :) + d_phi_d_r_j_alpha_first(dim_i, :, elec_l, cent_i, :) * ( sum(d_smooth_cutoff_g_d_r(:, 1:cent_i - 1, :) / smooth_cutoff_g(:, 1:cent_i - 1, :)  * vec_en_xyz_wlk(dim_i, :, 1:cent_i - 1, :) / dist_en_wlk(:, 1:cent_i - 1, :), 2) + sum(d_smooth_cutoff_g_d_r(:, cent_i + 1:ncent, :) / smooth_cutoff_g(:, cent_i + 1:ncent, :)  * vec_en_xyz_wlk(dim_i, :, cent_i + 1:ncent, :) / dist_en_wlk(:, cent_i + 1:ncent, :), 2) )
             end do
          end do
       end do

    end if

  end subroutine lap_j_phi_first_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_phi_second_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_i_phi_second. This is a rank four array.
    !               Laplacian of the phi:
    !
    !               (nabla_i)^2 phi(r_jI,r_iI,r_ji)       Note that this is derivative with respect to second index so no derivatives of the cutoff
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of electrons and labels the ith electron
    !               The third index has dimensions equal to the number of centers and labels the Ith nucleus
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 31 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j, cent_i
    character(len=max_string_len_rout), save :: lhere = 'lap_i_phi_second_bld'

    ! header
    if (header_exe) then

       call object_create('lap_i_phi_second')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed('ndim')
       call object_needed('all_elec')
       call object_needed('d_phi_d_r_j_inuc')
       call object_needed('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('d_phi_d_r_ji')
       call object_needed('dist_ee_wlk')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('d2_phi_d_r_ji_2')
       call object_needed('d2_phi_d_r_j_inuc_2')
       call object_needed('d2_phi_d_r_j_inuc_d_r_ji')

       return
    endif

    call object_alloc('lap_i_phi_second', lap_i_phi_second, nelec, nelec, ncent, nwalk)
    lap_i_phi_second(:,:,:,:) = 0.d0

    do elec_j = 1, nelec
       lap_i_phi_second(elec_j, :, :, :) = d2_phi_d_r_j_inuc_2(:, elec_j, :, :) + (ndim - 1.d0) * d_phi_d_r_j_inuc(:, elec_j, :, :) / dist_en_wlk(:, :, :) + d2_phi_d_r_ji_2(:, elec_j, :, :)
    end do

    do cent_i = 1, ncent
       do elec_j = 1, nelec
          lap_i_phi_second(elec_j, 1:elec_j - 1, cent_i, :) = lap_i_phi_second(elec_j, 1:elec_j - 1, cent_i, :) + (ndim - 1.d0) * d_phi_d_r_ji(1:elec_j - 1, elec_j, cent_i, :) / dist_ee_wlk(1:elec_j - 1, elec_j, :) +  2.d0 * d2_phi_d_r_j_inuc_d_r_ji(1:elec_j - 1, elec_j, cent_i, :) / dist_en_wlk(1:elec_j - 1, cent_i, :) / dist_ee_wlk(1:elec_j - 1, elec_j, :) * sum(vec_en_xyz_wlk(:, 1:elec_j - 1, cent_i, :) * vec_ee_xyz_wlk(:, 1:elec_j - 1, elec_j, :), 1)
          lap_i_phi_second(elec_j, elec_j + 1:nelec, cent_i, :) = lap_i_phi_second(elec_j, elec_j + 1:nelec, cent_i, :) + (ndim - 1.d0) * d_phi_d_r_ji(elec_j + 1:nelec, elec_j, cent_i, :) / dist_ee_wlk(elec_j + 1:nelec, elec_j, :) +  2.d0 * d2_phi_d_r_j_inuc_d_r_ji(elec_j + 1:nelec, elec_j, cent_i, :) / dist_en_wlk(elec_j + 1:nelec, cent_i, :) / dist_ee_wlk(elec_j + 1:nelec, elec_j, :) * sum(vec_en_xyz_wlk(:, elec_j + 1:nelec, cent_i, :) * vec_ee_xyz_wlk(:, elec_j + 1:nelec, elec_j, :), 1)
       end do
    end do

    if (all_elec) then
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       !cutoff
       do cent_i = 1, ncent
          do elec_j = 1, nelec
             lap_i_phi_second(:, elec_j, cent_i, :) = lap_i_phi_second(:, elec_j, cent_i, :) * product(smooth_cutoff_g(:, 1:cent_i - 1, :), 2) * product(smooth_cutoff_g(:, cent_i + 1:ncent, :), 2)
          end do
       end do
    end if

  end subroutine lap_i_phi_second_bld

  !==============================================================================================

  !======================================================================================

  subroutine d2_phi_d_r_ji_2_bld
    !---------------------------------------------------------------------------
    ! Description : build object d2_phi_d_r_ji_2(:,:,:,:). This is a rank four array.
    !
    !               d^2 phi(r_jI,r_iI,r_ji)    (without cutoff)
    !               ----------
    !               d r_(ji)^2
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 31 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d2_phi_d_r_ji_2')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_phi_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_ee_wlk_d_r')
       call object_needed ('d2_scaled_dist_een_ee_wlk_d_r_2')
       call object_needed('a_param_phi')

       return
    endif

    call object_alloc('d2_phi_d_r_ji_2', d2_phi_d_r_ji_2, nelec, nelec, ncent, nwalk)
    d2_phi_d_r_ji_2 (:,:,:,:) = 0.d0

    param_i = 1 !intialize parameter counter
    do order_p = 2, order_phi_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d2_phi_d_r_ji_2(elec_j, elec_i, cent_i, :) = d2_phi_d_r_ji_2(elec_j, elec_i, cent_i, :) + order_k * a_param_phi(param_i, iwctype(cent_i), up_up) * ( (order_k - 1) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-2) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) ** 2 + scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) * d2_scaled_dist_een_ee_wlk_d_r_2(elec_j, elec_i, :) ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l + scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_j = nup + 1, nelec
                         d2_phi_d_r_ji_2(elec_j, elec_i, cent_i, :) = d2_phi_d_r_ji_2(elec_j, elec_i, cent_i, :) + order_k * a_param_phi(param_i, iwctype(cent_i), up_down) * ( (order_k - 1) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-2) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) ** 2 + scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) * d2_scaled_dist_een_ee_wlk_d_r_2(elec_j, elec_i, :) ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l + scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d2_phi_d_r_ji_2(elec_j, elec_i, cent_i, :) = d2_phi_d_r_ji_2(elec_j, elec_i, cent_i, :) + order_k * a_param_phi(param_i, iwctype(cent_i), up_down) * ( (order_k - 1) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-2) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) ** 2 + scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) * d2_scaled_dist_een_ee_wlk_d_r_2(elec_j, elec_i, :) ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l + scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_j =  nup + 1, nelec
                         d2_phi_d_r_ji_2(elec_j, elec_i, cent_i, :) = d2_phi_d_r_ji_2(elec_j, elec_i, cent_i, :) + order_k * a_param_phi(param_i, iwctype(cent_i), down_down) * ( (order_k - 1) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-2) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) ** 2 + scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) * d2_scaled_dist_een_ee_wlk_d_r_2(elec_j, elec_i, :) ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l + scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d2_phi_d_r_ji_2(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d2_phi_d_r_ji_2_bld

  !=============================================================================================

  !======================================================================================

  subroutine d2_phi_d_r_j_inuc_d_r_ji_bld
    !---------------------------------------------------------------------------
    ! Description : build object d2_phi_d_r_j_inuc_d_r_ji(:,:,:,:). This is a rank four array.
    !
    !               d^2 phi(r_jI,r_iI,r_jI)    (without cutoff)
    !               ----------
    !               d r_jI d r_ji
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 31 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d2_phi_d_r_j_inuc_d_r_ji')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_phi_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_en_wlk_d_r')
       call object_needed ('d_scaled_dist_een_ee_wlk_d_r')
       call object_needed('a_param_phi')

       return
    endif

    call object_alloc('d2_phi_d_r_j_inuc_d_r_ji', d2_phi_d_r_j_inuc_d_r_ji, nelec, nelec, ncent, nwalk)
    d2_phi_d_r_j_inuc_d_r_ji (:,:,:,:) = 0.d0

    param_i = 1 !intialize parameter counter
    do order_p = 2, order_phi_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d2_phi_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) = d2_phi_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_up) * order_k * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k -1) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2)  * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                      do elec_j = nup+1, nelec
                         d2_phi_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) = d2_phi_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_down) * order_k * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k -1) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2)  * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d2_phi_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) = d2_phi_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_down) * order_k * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k -1) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2)  * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                      do elec_j = nup+1, nelec
                         d2_phi_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) = d2_phi_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), down_down) * order_k * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k -1) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2)  * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                   end do
                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d2_phi_d_r_j_inuc_d_r_ji(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d2_phi_d_r_j_inuc_d_r_ji_bld

  !=============================================================================================

  !======================================================================================

  subroutine d2_phi_d_r_j_inuc_2_bld
    !---------------------------------------------------------------------------
    ! Description : build object d2_phi_d_r_j_inuc_2(:,:,:,:). This is a rank four array.
    !
    !               d phi(r_jI,r_iI,r_jI)    (without cutoff)
    !               ----------
    !               d r_jI
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d2_phi_d_r_j_inuc_2')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_phi_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_en_wlk_d_r')
       call object_needed ('d2_scaled_dist_een_en_wlk_d_r_2')
       call object_needed('a_param_phi')

       return
    endif

    call object_alloc('d2_phi_d_r_j_inuc_2', d2_phi_d_r_j_inuc_2, nelec, nelec, ncent, nwalk)
    d2_phi_d_r_j_inuc_2 (:,:,:,:) = 0.d0

    param_i = 1 !intialize parameter counter
    do order_p = 2, order_phi_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d2_phi_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) = d2_phi_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_up) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k * scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) * ( d2_scaled_dist_een_en_wlk_d_r_2(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2)  + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2)) + d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) **2 * ( (order_p - order_k + order_l) / 2 * (order_p - order_k + order_l - 2) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-4)/2)  + (order_p - order_k - order_l) / 2  * (order_p - order_k - order_l - 2) / 2 * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-4)/2)) )
                      end do
                      do elec_j = nup+1, nelec
                         d2_phi_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) = d2_phi_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_down) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k * scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) * ( d2_scaled_dist_een_en_wlk_d_r_2(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2)  + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2)) + d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) **2 * ( (order_p - order_k + order_l) / 2 * (order_p - order_k + order_l - 2) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-4)/2)  + (order_p - order_k - order_l) / 2  * (order_p - order_k - order_l - 2) / 2 * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-4)/2)) )
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d2_phi_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) = d2_phi_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), up_down) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k * scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) * ( d2_scaled_dist_een_en_wlk_d_r_2(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2)  + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2)) + d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) **2 * ( (order_p - order_k + order_l) / 2 * (order_p - order_k + order_l - 2) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-4)/2)  + (order_p - order_k - order_l) / 2  * (order_p - order_k - order_l - 2) / 2 * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-4)/2)) )
                      end do
                      do elec_j = nup+1, nelec
                         d2_phi_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) = d2_phi_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) + a_param_phi(param_i, iwctype(cent_i), down_down) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k * scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) * ( d2_scaled_dist_een_en_wlk_d_r_2(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2)  + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2)) + d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) **2 * ( (order_p - order_k + order_l) / 2 * (order_p - order_k + order_l - 2) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-4)/2)  + (order_p - order_k - order_l) / 2  * (order_p - order_k - order_l - 2) / 2 * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-4)/2)) )
                      end do
                   end do
                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d2_phi_d_r_j_inuc_2(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d2_phi_d_r_j_inuc_2_bld

  !=============================================================================================

  !======================================================================================

  subroutine xi_een_theta_bld
    !---------------------------------------------------------------------------
    ! Description : build object xi_een_theta(:,:,:). This is a rank three array.
    !               It is part of electron-electron-nuclear backflow.
    !               It is necessary for the construction of the total backflow transformation.
    !               The first index has dimension equal to the number of spatial dimensions.
    !               The second index has dimensions equal to the number of electrons.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: dim_i

    ! header
    if (header_exe) then

       call object_create('xi_een_theta')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('ndim')
       call object_needed('theta')

       return
    endif

    call object_alloc('xi_een_theta', xi_een_theta, ndim, nelec, nwalk)
    xi_een_theta(:,:,:) = 0.d0

    do dim_i = 1, ndim
       xi_een_theta(dim_i, :, :) = sum(sum(theta(:, :, :, :), 2) * vec_en_xyz_wlk(dim_i, :, :, :),2)
    end do

  end subroutine xi_een_theta_bld

  !======================================================================================

  !======================================================================================

  subroutine theta_bld
    !---------------------------------------------------------------------------
    ! Description : build object theta(:,:,:,:). This is a rank four array.
    !               It is necessary for the construction of the electron-electron-nuclear backflow terms in all electron case.
    !               The first index has dimension equal to the number of electrons.
    !               The second index has dimension equal to the number of electrons.
    !               The third index has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin
    character(len=max_string_len_rout), save :: lhere = 'theta_bld'

    ! header
    if (header_exe) then

       call object_create('theta')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('nelec')
       call object_needed('iwctype')
       call object_needed('nup')
       call object_needed('order_theta_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed('b_param_theta')
       call object_needed('all_elec')

       return
    endif

    call object_alloc('theta', theta, nelec, nelec, ncent, nwalk)
    theta (:,:,:,:) = 0.d0

    param_i = 1
    do order_p = 2, order_theta_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_j = 1, nup
                      do elec_i = 1, nup
                         theta(elec_i, elec_j, cent_i, :) = theta(elec_i, elec_j, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_up) &
                              & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, :)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l &
                              & +  scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :) &
                              & *  scaled_dist_een_en_wlk(elec_j, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_i = nup + 1, nelec
                         theta(elec_i, elec_j, cent_i, :) = theta(elec_i, elec_j, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, :)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l +  scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l )&
                              & * ( scaled_dist_een_en_wlk(elec_i, cent_i, :) *  scaled_dist_een_en_wlk(elec_j, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                   do elec_j = nup + 1, nelec
                      do elec_i = 1, nup
                         theta(elec_i, elec_j, cent_i, :) = theta(elec_i, elec_j, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, :)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l +  scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l )&
                              & * ( scaled_dist_een_en_wlk(elec_i, cent_i, :) *  scaled_dist_een_en_wlk(elec_j, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_i = nup + 1, nelec
                         theta(elec_i, elec_j, cent_i, :) = theta(elec_i, elec_j, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), down_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, :)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l +  scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l )&
                              & * ( scaled_dist_een_en_wlk(elec_i, cent_i, :) *  scaled_dist_een_en_wlk(elec_j, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                end do
                param_i = param_i + 1
             end if
          enddo
       enddo
    enddo

    if (all_elec) then
       !apply smooth cutoff function
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       do cent_i=1, ncent
          do elec_j=1, nelec
             theta(:,elec_j,cent_i,:) = theta(:,elec_j,cent_i,:) * product(smooth_cutoff_g(:,1:cent_i-1,:),2) * product(smooth_cutoff_g(:,cent_i+1:ncent,:),2)
          end do
       enddo
    end if

    !no self backflow
    do elec_i=1, nelec
       theta(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine theta_bld

  !=============================================================================================

  !==============================================================================================

  subroutine b_param_theta_bld

    !---------------------------------------------------------------------------
    ! Description : build b_param_theta. This is a rank 3 array.
    !               List of a parameters for theta backflow function for each center type, and each spin dependence
    !               The first index is which parameter.
    !               The second index is the type of center.
    !               The third index is the type of spin dependence.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: b_param_nb, dep_b_param_i, ctype_i, spin_dependency_i
    character(len=max_string_len_rout), save :: lhere = 'b_param_theta_bld'
    integer, parameter :: spin_dependencies_nb = 3    !because three choices up_up, down_down, up_down

    ! header
    if (header_exe) then

       call object_create('b_param_theta')

       call object_needed('nctype')
       call object_needed('theta_spin_dependence')
       call object_needed('read_b_param')
       call object_needed('read_b_param_theta_nb')
       call object_needed('dep_b_param_theta')
       call object_needed('b_param_theta_cond')
       call object_needed('nb_b_param_theta')

       return
    endif

    !need to have the option of using read_b_param_theta if it is read in
    if (read_b_param) then
       call object_provide_in_node(lhere, 'read_b_param_theta')
    end if

    call object_alloc('b_param_theta', b_param_theta, nb_b_param_theta, nctype, spin_dependencies_nb)
    b_param_theta=0.d0

    !do we read in parameters?
    if (read_b_param) then
       if (theta_spin_dependence == spin_dependencies_nb) then
          b_param_nb = nb_b_param_theta * nctype * spin_dependencies_nb
          if (read_b_param_theta_nb == b_param_nb) then
             b_param_theta = reshape( read_b_param_theta, (/ nb_b_param_theta, nctype, spin_dependencies_nb /))
          else
             call require(lhere, 'number of parameters in read_b_param_theta equal to nb_b_param_theta * nctype * 3', .false.)
          end if
       elseif( theta_spin_dependence == spin_dependencies_nb - 1) then
          b_param_nb = nb_b_param_theta * nctype * (spin_dependencies_nb - 1)
          if (read_b_param_theta_nb == b_param_nb) then
             b_param_theta = reshape( read_b_param_theta, (/ nb_b_param_theta, nctype, spin_dependencies_nb /), pad=read_b_param_theta)
             !want up_up in b_param_theta(:,:,1), up_up=down_down in b_param_theta(:,:,2), up_down in b_param_theta(:,:,3)
             b_param_theta(:,:,3)=b_param_theta(:,:,2)
             b_param_theta(:,:,2)=b_param_theta(:,:,1)
          else
             call require(lhere, 'number of parameters in read_b_param_theta equal to nb_b_param_theta * nctype * 2', .false.)
          end if
       else
          b_param_nb = nb_b_param_theta * nctype
          if (read_b_param_theta_nb == b_param_nb) then
             b_param_theta = reshape( read_b_param_theta, (/ nb_b_param_theta, nctype, spin_dependencies_nb /), pad=read_b_param_theta)
          else
             call require(lhere, 'number of parameters in read_b_param_theta equal to nb_b_param_theta * nctype', .false.)
          end if
       end if
       !don't want to re-read
       read_b_param = .false.
    else
       !start from random guess
       call random_number(b_param_theta(:,:,:))
    endif

    !impose cusp
    !To impose the cusp we take the pivot positions in our reduced row echelon form of the matrix of cusp constraints
    !and we set this dependent parameter equal to the negative of the dot product of the rest of the row with
    !the corresponding a parameters. This is nothing fancy, but simply a standard way of solving linear equations

    do spin_dependency_i=1, theta_spin_dependence
       do ctype_i=1, nctype
          do dep_b_param_i=1, size(dep_b_param_theta,1)
             b_param_theta(dep_b_param_theta(dep_b_param_i), ctype_i, spin_dependency_i) = dot_product( -b_param_theta_cond(dep_b_param_i, dep_b_param_theta(dep_b_param_i)+1:), b_param_theta(dep_b_param_theta(dep_b_param_i) + 1:, ctype_i, spin_dependency_i))
          end do
       end do
    end do

    if (theta_spin_dependence .eq. 2) then
       !require up-up = down-down
       b_param_theta(:,:,2) = b_param_theta(:,:,1)
    elseif (theta_spin_dependence .eq. 1) then
       !require up-up = down-down=up-down
       b_param_theta(:,:,2) = b_param_theta(:,:,1)
       b_param_theta(:,:,3) = b_param_theta(:,:,1)
    end if

  end subroutine b_param_theta_bld

  !==============================================================================================

  !==============================================================================================

  subroutine nb_b_param_theta_bld

    !---------------------------------------------------------------------------
    ! Description : build nb_b_param_theta. This is an integer.
    !               Total number of a parameters in theta per nucleus type, per spin dependence.
    !               To get total number of a parameters multiply by number of nuclei types and
    !               spin dependence.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_p, order_k, order_l, order_l_max

    ! header
    if (header_exe) then

       call object_create('nb_b_param_theta')

       call object_needed('order_theta_bf')

       return
    endif

    call object_associate('nb_b_param_theta',nb_b_param_theta)

    !determine the number of coefficients for the three body bf
    nb_b_param_theta = 0 !initialize parameter counter
    do order_p = 2, order_theta_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                nb_b_param_theta = nb_b_param_theta + 1
             endif
          enddo
       enddo
    enddo

  end subroutine nb_b_param_theta_bld

  !==============================================================================================

  !==============================================================================================

  subroutine dep_b_param_theta_bld

    !---------------------------------------------------------------------------
    ! Description : build dep_b_param_theta. This is a rank one array.
    !               This produces a list of the dependent parameters from b_param_theta_cond.
    !               To do this we find the pivots in a row reduced matrix.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer :: nb_dep_b_param_theta
    integer :: row_i, pivot_column_i, column_i

    ! header
    if (header_exe) then

       call object_create('dep_b_param_theta')

       call object_needed('nb_b_param_theta')
       call object_needed('b_param_theta_cond')

       return
    endif

    !initialize the number of dependent parameters
    nb_dep_b_param_theta = 0
    do row_i = 1, size(b_param_theta_cond,1)
       if (sum(abs(b_param_theta_cond(row_i,:))) < 1.d-10) then
          !no pivot
          cycle
       end if
       !pivot found, so increment counter
       nb_dep_b_param_theta = nb_dep_b_param_theta + 1
    end do

    call object_alloc('dep_b_param_theta', dep_b_param_theta, nb_dep_b_param_theta)

    !intialize the pivot column
    pivot_column_i = 1
    do row_i = 1, nb_dep_b_param_theta
       do column_i = pivot_column_i, nb_b_param_theta
          !first non zero element  (should be a 1)
          if (b_param_theta_cond(row_i,column_i) > 1.d-10) then
             !pivot found in this row
             dep_b_param_theta(row_i) = column_i
             !next pivot column must be to the right of this one
             pivot_column_i=column_i + 1
             !so move on to next row
             exit
          end if
       end do
    end do

  end subroutine dep_b_param_theta_bld

  !==============================================================================================

  !==============================================================================================

  subroutine b_param_theta_cond_bld

    !---------------------------------------------------------------------------
    ! Description : build b_param_theta_cond_bld. This is a rank two array.
    !               The first index is the number of condition ( some are possibly redundant).
    !               The second index is the number a parameters.
    !               We use gaussian elimination to row reduce the conditions and check for redundancies.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------

    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_p, order_k, order_l, order_l_max, param_i

    ! header
    if (header_exe) then

       call object_create('b_param_theta_cond')

       call object_needed('order_theta_bf')
       call object_needed('nb_b_param_theta')
       call object_needed('all_elec')

       return
    endif

    if (all_elec) then
       ! one derivative constraints evaluated at 0, one constraint on function so number of contraints is 2 * order_theta_bf - 1
       call object_alloc('b_param_theta_cond', b_param_theta_cond, 2 * order_theta_bf - 1, nb_b_param_theta)
    else
       ! one derivative constraints so number of contraints is (order_theta_bf -1 )
       call object_alloc('b_param_theta_cond', b_param_theta_cond, order_theta_bf - 1, nb_b_param_theta)
    end if

    b_param_theta_cond(:,:) = 0.d0

    param_i = 0 !intialize parameter counter
    !setup matrix of constraints
    do order_p = 2, order_theta_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                param_i = param_i + 1
                !=================================
                ! d theta|
                ! -----  |        = 0
                ! d r_ij |
                !         r_ij=0
                !=================================
                !not necessary to have if condition for order_k = 0 (pseudo case) because there are rows below these
                if (order_k .ne. 0) then
                   b_param_theta_cond(order_p-order_k, param_i) = b_param_theta_cond(order_p-order_k, param_i) + order_k
                end if
                if (all_elec) then
                   !================================================
                   ! theta|
                   !      |        = 0
                   !      |
                   !       r_iI=0
                   !===============================================
                   !note that the row needs to be shifted down by (order_theta_bf - 1) because those rows hold the above constraints
                   b_param_theta_cond((order_theta_bf - 1) + (order_p+order_k-order_l)/2, param_i) = b_param_theta_cond((order_theta_bf -1 ) + (order_p+order_k-order_l)/2, param_i) + 1.d0
                   b_param_theta_cond((order_theta_bf - 1) + (order_p+order_k+order_l)/2, param_i) = b_param_theta_cond((order_theta_bf -1 ) + (order_p+order_k+order_l)/2, param_i) + 1.d0
                end if
             endif
          enddo
       enddo
    enddo

    !use gaussian elimintation on constraints to remove redundancies
    call row_reduce(b_param_theta_cond)

  end subroutine b_param_theta_cond_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_xi_een_theta_j_beta_d_r_i_alpha_bld

    !---------------------------------------------------------------------------
    ! Description : build d_xi_een_theta_j_beta_d_r_i_alpha. This is a rank five array.
    !               Derivative of the theta part of electron-electron-nuclear backflow function
    !               with respect to electron position.
    !               It is necessary for the construction of the energy, drift velocity.
    !               The first and second indices are the number of spatial dimensions.
    !               The first index is for the components of xi_een_theta (beta)
    !               The second index is for the components of the gradient (alpha)
    !               The third and fourth indices have dimensions equal to the number of electrons.
    !               The third index is the label of the een theta backflow transformation of each electon
    !               The fourth index is for the electron that the gradient is taken with respect to
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim, ncent
    include 'commons.h'

    !local
    integer  :: elec_j, dim_i, dim_j

    ! header
    if (header_exe) then

       call object_create('d_xi_een_theta_j_beta_d_r_i_alpha')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('ndim')
       call object_needed('d_theta_d_r_j_alpha_first')
       call object_needed('d_theta_d_r_i_alpha_second')
       call object_needed('theta')

       return
    endif

    call object_alloc('d_xi_een_theta_j_beta_d_r_i_alpha', d_xi_een_theta_j_beta_d_r_i_alpha, ndim, ndim, nelec, nelec, nwalk)
    d_xi_een_theta_j_beta_d_r_i_alpha(:,:,:,:,:) = 0.d0

    do elec_j = 1, nelec
       do dim_i = 1, ndim     !alpha
          !delta_{alpha,beta}
          d_xi_een_theta_j_beta_d_r_i_alpha(dim_i, dim_i, elec_j, elec_j, :) = sum(sum(theta(elec_j, :, :, :),1),1)
          do dim_j = 1, ndim   !beta
             d_xi_een_theta_j_beta_d_r_i_alpha(dim_j, dim_i, elec_j, elec_j, :) = d_xi_een_theta_j_beta_d_r_i_alpha(dim_j, dim_i, elec_j, elec_j, :) + sum(sum(d_theta_d_r_j_alpha_first(dim_i, elec_j, :, :, :),1) * vec_en_xyz_wlk(dim_j, elec_j, :, :), 1)
             d_xi_een_theta_j_beta_d_r_i_alpha(dim_j, dim_i, :, elec_j, :) = d_xi_een_theta_j_beta_d_r_i_alpha(dim_j, dim_i, :, elec_j, :) + sum(d_theta_d_r_i_alpha_second(dim_i, :, elec_j, :, :) * vec_en_xyz_wlk(dim_j, :, :, :),2)
          end do
       end do
    end do

  end subroutine d_xi_een_theta_j_beta_d_r_i_alpha_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_xi_een_theta_j_beta_d_r_i_alpha_num_bld
    !---------------------------------------------------------------------------
    ! Description : build d_xi_een_theta_j_beta_d_r_i_alpha__num. This is a rank five array.
    !               Numerical first derivative of the xi_een_theta:
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               The first and second indices are the number of spatial dimensions.
    !               The first index is for the components of xi_een_theta (beta)
    !               The second index is for the components of the gradient (alpha)
    !               The third and fourth indices have dimensions equal to the number of electrons.
    !               The third index is the label of the een_theta backflow transformation of each electon
    !               The fourth index is for the electron that the gradient is taken with respect to
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 14 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, dim_i, dim_j
    real(dp), allocatable :: xi_een_theta_temp_plus(:,:,:), xi_een_theta_temp_minus(:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'd_xi_een_theta_j_beta_d_r_i_alpha_num_bld'

    ! header
    if (header_exe) then

       call object_create('d_xi_een_theta_j_beta_d_r_i_alpha_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('xi_een_theta')
       call object_needed('coord_elec_wlk')

       return
    endif

    ! allocation
    call object_alloc('d_xi_een_theta_j_beta_d_r_i_alpha_num', d_xi_een_theta_j_beta_d_r_i_alpha_num, ndim, ndim, nelec, nelec, nwalk)
    d_xi_een_theta_j_beta_d_r_i_alpha_num(:,:,:,:,:) = 0.d0

    !allocate temporary objects
    allocate(xi_een_theta_temp_plus(ndim, nelec, nwalk))
    allocate(xi_een_theta_temp_minus(ndim, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

    do elec_j = 1, nelec
       do elec_i = 1, nelec
          do dim_j = 1, ndim  !beta
             do dim_i = 1, ndim  !alph
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) + epsilon
                call object_modified('coord_elec_wlk')
                !get updated xi_een_theta
                call object_provide_in_node(lhere, 'xi_een_theta')
                xi_een_theta_temp_plus = xi_een_theta
                coord_elec_wlk = coord_elec_wlk_temp
                !shift coordinates in a particular dimension, x,y,z, etc
                coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) - epsilon
                call object_modified('coord_elec_wlk')
                !get updated xi_een_theta
                call object_provide_in_node(lhere, 'xi_een_theta')
                xi_een_theta_temp_minus = xi_een_theta
                !calculate derivative
                d_xi_een_theta_j_beta_d_r_i_alpha_num(dim_j, dim_i, elec_j, elec_i, :) = (xi_een_theta_temp_plus(dim_j, elec_j, :) - xi_een_theta_temp_minus(dim_j, elec_j, :) ) / (2.d0 * epsilon)
             end do
          end do
       end do
    end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'xi_een_theta')

  end subroutine d_xi_een_theta_j_beta_d_r_i_alpha_num_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_theta_d_r_i_alpha_second_bld

    !---------------------------------------------------------------------------
    ! Description : build d_theta_d_r_i_alpha_second. This is a rank five array.
    !               First derivative of the theta:
    !
    !               d theta(r_jI,r_iI,r_ji)
    !               ----------
    !               d r_i_alpha
    !
    !               Note alpha refers to cartesian coords. In 3-d, r_alpha takes on values x,y,z.
    !               The first index is the over spatial dimension. This is the components of the gradient.
    !               The second index has dimensions equal to the number of electrons and labels the jth electron
    !               The third index has dimensions equal to the number of electrons and labels the ith electron
    !               The fourth index has dimension equal to the number of centers.
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim, ncent
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, dim_i
    character(len=max_string_len_rout), save :: lhere = 'd_theta_d_r_i_alpha_second_bld'


    ! header
    if (header_exe) then

       call object_create('d_theta_d_r_i_alpha_second')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed ('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed ('dist_ee_wlk')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('d_theta_d_r_j_inuc')
       call object_needed('d_theta_d_r_ji')
       call object_needed('all_elec')

       return
    endif

    call object_alloc('d_theta_d_r_i_alpha_second', d_theta_d_r_i_alpha_second, ndim, nelec, nelec, ncent, nwalk)
    d_theta_d_r_i_alpha_second(:,:,:,:,:) = 0.d0


    do elec_j = 1, nelec
       do dim_i = 1, ndim
          d_theta_d_r_i_alpha_second(dim_i, elec_j, :, :, :) = d_theta_d_r_j_inuc(:, elec_j, :, :) * vec_en_xyz_wlk(dim_i, :, :, :) / dist_en_wlk(:, :, :)
       end do
    end do

    do cent_i = 1, ncent
       do elec_j = 1, nelec
          do dim_i = 1, ndim
             d_theta_d_r_i_alpha_second(dim_i, elec_j, 1:elec_j-1, cent_i, :) = d_theta_d_r_i_alpha_second(dim_i, elec_j, 1:elec_j-1, cent_i, :) - d_theta_d_r_ji(elec_j, 1:elec_j-1, cent_i, :) * vec_ee_xyz_wlk(dim_i, elec_j, 1:elec_j-1, :) / dist_ee_wlk(elec_j, 1:elec_j-1, :)
             d_theta_d_r_i_alpha_second(dim_i, elec_j, elec_j+1:nelec, cent_i, :) = d_theta_d_r_i_alpha_second(dim_i, elec_j, elec_j+1:nelec, cent_i, :) - d_theta_d_r_ji(elec_j, elec_j+1:nelec, cent_i, :) * vec_ee_xyz_wlk(dim_i, elec_j, elec_j+1:nelec, :) / dist_ee_wlk(elec_j, elec_j+1:nelec, :)
          end do
       end do
    end do

    if (all_elec) then
       !apply smooth cutoff function
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       do cent_i=1, ncent
          do elec_i=1, nelec
             do dim_i = 1, ndim
                d_theta_d_r_i_alpha_second(dim_i, :, elec_i, cent_i, :) = d_theta_d_r_i_alpha_second(dim_i, :, elec_i, cent_i, :) * product(smooth_cutoff_g(:,1:cent_i-1,:),2) * product(smooth_cutoff_g(:,cent_i+1:ncent,:),2)
             end do
          enddo
       end do
    end if

  end subroutine d_theta_d_r_i_alpha_second_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_theta_d_r_j_alpha_first_bld

    !---------------------------------------------------------------------------
    ! Description : build d_theta_d_r_j_alpha_first. This is a rank five array.
    !               First derivative of the theta:
    !
    !               d theta(r_jI,r_lI,r_jl)
    !               ----------
    !               d r_j_alpha
    !
    !               Note alpha refers to cartesian coords. In 3-d, r_alpha takes on values x,y,z.
    !               The first index is the over spatial dimension. This is the components of the gradient.
    !               The second index has dimensions equal to the number of electrons and labels the jth electron
    !               The third index has dimensions equal to the number of electrons and labels the lth electron
    !               The fourth index has dimension equal to the number of centers.
    !               The fifth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim, ncent
    include 'commons.h'

    !local
    integer  :: elec_j, elec_l, cent_i, dim_i
    character(len=max_string_len_rout), save :: lhere = 'd_theta_d_r_j_alpha_first_bld'

    ! header
    if (header_exe) then

       call object_create('d_theta_d_r_j_alpha_first')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed ('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed ('dist_ee_wlk')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('d_theta_d_r_j_inuc')
       call object_needed('d_theta_d_r_ji')
       call object_needed('all_elec')

       return
    endif

    call object_alloc('d_theta_d_r_j_alpha_first', d_theta_d_r_j_alpha_first, ndim, nelec, nelec, ncent, nwalk)
    d_theta_d_r_j_alpha_first(:,:,:,:,:) = 0.d0

    do elec_l = 1, nelec
       do dim_i = 1, ndim
          d_theta_d_r_j_alpha_first(dim_i, :, elec_l, :, :) = d_theta_d_r_j_inuc(:, elec_l, :, :) * vec_en_xyz_wlk(dim_i, :, :, :) / dist_en_wlk(:, :, :)
       end do
    end do

    do cent_i = 1, ncent
       do elec_j = 1, nelec
          do dim_i = 1, ndim
             d_theta_d_r_j_alpha_first(dim_i, elec_j, 1:elec_j-1, cent_i, :) = d_theta_d_r_j_alpha_first(dim_i, elec_j, 1:elec_j-1, cent_i, :) + d_theta_d_r_ji(elec_j, 1:elec_j-1, cent_i, :) * vec_ee_xyz_wlk(dim_i, elec_j, 1:elec_j-1, :) / dist_ee_wlk(elec_j, 1:elec_j-1, :)
             d_theta_d_r_j_alpha_first(dim_i, elec_j, elec_j+1:nelec, cent_i, :) = d_theta_d_r_j_alpha_first(dim_i, elec_j, elec_j+1:nelec, cent_i, :) + d_theta_d_r_ji(elec_j, elec_j+1:nelec, cent_i, :) * vec_ee_xyz_wlk(dim_i, elec_j, elec_j+1:nelec, :) / dist_ee_wlk(elec_j, elec_j+1:nelec, :)
          end do
       end do
    end do

    if (all_elec) then
       !apply smooth cutoff function
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       call object_provide_in_node(lhere, 'd_smooth_cutoff_g_d_r')
       call object_provide_in_node(lhere, 'theta')
       do cent_i=1, ncent
          do elec_l=1, nelec
             do dim_i = 1, ndim
                d_theta_d_r_j_alpha_first(dim_i, :, elec_l, cent_i, :) = d_theta_d_r_j_alpha_first(dim_i, :, elec_l, cent_i, :) * product(smooth_cutoff_g(:,1:cent_i-1,:),2) * product(smooth_cutoff_g(:,cent_i+1:ncent,:),2) + theta(:, elec_l, cent_i, :) * ( sum(vec_en_xyz_wlk(dim_i, :, 1:cent_i-1, :) / dist_en_wlk(:, 1:cent_i-1, :) / smooth_cutoff_g(:, 1:cent_i-1, :) * d_smooth_cutoff_g_d_r(:, 1:cent_i-1, :),2) + sum(vec_en_xyz_wlk(dim_i, :, cent_i+1:ncent, :)/ dist_en_wlk(:, cent_i+1:ncent, :) / smooth_cutoff_g(:, cent_i+1:ncent, :) * d_smooth_cutoff_g_d_r(:, cent_i+1:ncent, :),2))
             end do
          enddo
       end do
    end if

  end subroutine d_theta_d_r_j_alpha_first_bld

  !==============================================================================================

  !======================================================================================

  subroutine d_theta_d_r_ji_bld
    !---------------------------------------------------------------------------
    ! Description : build object d_theta_d_r_ji(:,:,:,:). This is a rank four array.
    !
    !               d theta(r_jI,r_iI,r_ji)    (without cutoff)
    !               ----------
    !               d r_ji
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d_theta_d_r_ji')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_theta_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_ee_wlk_d_r')
       call object_needed('b_param_theta')

       return
    endif

    call object_alloc('d_theta_d_r_ji', d_theta_d_r_ji, nelec, nelec, ncent, nwalk)
    d_theta_d_r_ji (:,:,:,:) = 0.d0


    param_i = 1 !intialize parameter counter
    do order_p = 2, order_theta_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d_theta_d_r_ji(elec_j, elec_i, cent_i, :) = d_theta_d_r_ji(elec_j, elec_i, cent_i, :) + order_k * b_param_theta(param_i, iwctype(cent_i), up_up)&
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) ) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l +  &
                              & scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_j = nup + 1, nelec
                         d_theta_d_r_ji(elec_j, elec_i, cent_i, :) = d_theta_d_r_ji(elec_j, elec_i, cent_i, :) + order_k * b_param_theta(param_i, iwctype(cent_i), up_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) ) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l +  &
                              & scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d_theta_d_r_ji(elec_j, elec_i, cent_i, :) = d_theta_d_r_ji(elec_j, elec_i, cent_i, :) + order_k * b_param_theta(param_i, iwctype(cent_i), up_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) ) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l +  &
                              & scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_j =  nup + 1, nelec
                         d_theta_d_r_ji(elec_j, elec_i, cent_i, :) = d_theta_d_r_ji(elec_j, elec_i, cent_i, :) + order_k * b_param_theta(param_i, iwctype(cent_i), down_down)&
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) ) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l +  &
                              & scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d_theta_d_r_ji(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d_theta_d_r_ji_bld

  !=============================================================================================

  !======================================================================================

  subroutine d_theta_d_r_j_inuc_bld
    !---------------------------------------------------------------------------
    ! Description : build object d_theta_d_r_j_inuc(:,:,:,:). This is a rank four array.
    !
    !               d theta(r_jI,r_iI,r_jI)    (without cutoff)
    !               ----------
    !               d r_jI
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d_theta_d_r_j_inuc')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_theta_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_en_wlk_d_r')
       call object_needed('b_param_theta')

       return
    endif

    call object_alloc('d_theta_d_r_j_inuc', d_theta_d_r_j_inuc, nelec, nelec, ncent, nwalk)
    d_theta_d_r_j_inuc (:,:,:,:) = 0.d0

    param_i = 1 !intialize parameter counter
    do order_p = 2, order_theta_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d_theta_d_r_j_inuc(elec_j, elec_i, cent_i, :) = d_theta_d_r_j_inuc(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_up) &
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k ) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) &
                              & * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) &
                              & + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                      do elec_j = nup+1, nelec
                         d_theta_d_r_j_inuc(elec_j, elec_i, cent_i, :) = d_theta_d_r_j_inuc(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_down) &
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k ) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) &
                              & * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) &
                              & + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d_theta_d_r_j_inuc(elec_j, elec_i, cent_i, :) = d_theta_d_r_j_inuc(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_down) &
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k ) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) &
                              & * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) &
                              & + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                      do elec_j = nup+1, nelec
                         d_theta_d_r_j_inuc(elec_j, elec_i, cent_i, :) = d_theta_d_r_j_inuc(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), down_down) &
                              & * ( scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k ) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) &
                              & * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) &
                              & + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                   end do

                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d_theta_d_r_j_inuc(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d_theta_d_r_j_inuc_bld

  !=============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_een_theta_j_beta_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_een_theta_j_beta. This is a rank four array.
    !               Laplacian of the electron-electron-nuclear theta backflow function
    !               with respect to electron position.
    !               It is necessary for the construction of the energy.
    !               The first index is for the components of xi_een_theta (beta)
    !               The second index is the label of the ee backflow transformed electon
    !               The third index is for the electron that the laplacian is taken with respect to
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j, dim_j

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_een_theta_j_beta')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('ndim')
       call object_needed('d_theta_d_r_j_alpha_first')
       call object_needed('lap_j_theta_first')
       call object_needed('lap_i_theta_second')

       return
    endif

    call object_alloc('lap_i_xi_een_theta_j_beta', lap_i_xi_een_theta_j_beta, ndim, nelec, nelec, nwalk)
    lap_i_xi_een_theta_j_beta(:,:,:,:) = 0.d0

    do elec_j = 1, nelec
       do dim_j = 1, ndim
          lap_i_xi_een_theta_j_beta(dim_j, elec_j, elec_j, :) =  sum(sum(lap_j_theta_first(elec_j,:,:,:), 1) * vec_en_xyz_wlk(dim_j, elec_j, :, :) + 2.d0  * sum(d_theta_d_r_j_alpha_first(dim_j, elec_j, :, :, :), 1), 1)
          lap_i_xi_een_theta_j_beta(dim_j, :, elec_j, :) = lap_i_xi_een_theta_j_beta(dim_j, :, elec_j, :) + sum(lap_i_theta_second(:,elec_j,:,:) * vec_en_xyz_wlk(dim_j, :, :, :), 2)
       end do
    end do

  end subroutine lap_i_xi_een_theta_j_beta_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_xi_een_theta_j_beta_num_bld
    !---------------------------------------------------------------------------
    ! Description : build lap_i_xi_een_theta_j_beta_num. This is a rank four array.
    !               Numerical laplacian of xi_een_theta:
    !
    !               (nabla_i)^2 xi_een_theta_j^beta
    !
    !               NOTE THIS IS TO CHECK ANALYTIC DERIVATIVES
    !
    !               The first index is for the components of xi_een_theta (beta)
    !               The second index is the label of the een_theta backflow transformation of each electon
    !               The third index is for the electron that the laplacian is taken with respect to
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 26 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, dim_i, dim_j
    real(dp), allocatable :: d_xi_een_theta_j_beta_d_r_i_alpha_plus(:,:,:,:,:), d_xi_een_theta_j_beta_d_r_i_alpha_minus(:,:,:,:,:)
    real(dp), allocatable :: coord_elec_wlk_temp(:,:,:)
    real(dp) :: epsilon = 1.d-7
    character(len=max_string_len_rout), save :: lhere = 'lap_i_xi_een_theta_j_beta_num_bld'

    ! header
    if (header_exe) then

       call object_create('lap_i_xi_een_theta_j_beta_num')

       call object_needed('ndim')
       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('coord_elec_wlk')
       call object_needed('d_xi_een_theta_j_beta_d_r_i_alpha')

       return
    endif

    ! allocation
    call object_alloc('lap_i_xi_een_theta_j_beta_num', lap_i_xi_een_theta_j_beta_num, ndim, nelec, nelec, nwalk)
    lap_i_xi_een_theta_j_beta_num(:,:,:,:) = 0.d0

    !allocate temporary objects
    allocate(d_xi_een_theta_j_beta_d_r_i_alpha_plus(ndim, ndim, nelec, nelec, nwalk))
    allocate(d_xi_een_theta_j_beta_d_r_i_alpha_minus(ndim, ndim, nelec, nelec, nwalk))
    allocate(coord_elec_wlk_temp(ndim, nelec, nwalk))

    !store original values
    coord_elec_wlk_temp = coord_elec_wlk

     do elec_i = 1, nelec
        do elec_j = 1, nelec
           do dim_j = 1, ndim  !beta
              do dim_i = 1, ndim  !alpha
                 coord_elec_wlk = coord_elec_wlk_temp
                 !shift coordinates in a particular dimension, x,y,z, etc
                 coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) + epsilon
                 call object_modified('coord_elec_wlk')
                 !get updated d_xi_een_theta_j_beta_d_r_i_alpha
                 call object_provide_in_node(lhere, 'd_xi_een_theta_j_beta_d_r_i_alpha')
                 d_xi_een_theta_j_beta_d_r_i_alpha_plus = d_xi_een_theta_j_beta_d_r_i_alpha
                 coord_elec_wlk = coord_elec_wlk_temp
                 !shift coordinates in a particular dimension, x,y,z, etc
                 coord_elec_wlk(dim_i, elec_i, :) = coord_elec_wlk(dim_i, elec_i, :) - epsilon
                 call object_modified('coord_elec_wlk')
                 !get updated d_xi_een_theta_j_beta_d_r_i_alpha
                 call object_provide_in_node(lhere, 'd_xi_een_theta_j_beta_d_r_i_alpha')
                 d_xi_een_theta_j_beta_d_r_i_alpha_minus = d_xi_een_theta_j_beta_d_r_i_alpha
                 !calculate derivative
                 lap_i_xi_een_theta_j_beta_num(dim_j, elec_j, elec_i, :) = lap_i_xi_een_theta_j_beta_num(dim_j, elec_j, elec_i, :) + (d_xi_een_theta_j_beta_d_r_i_alpha_plus(dim_j, dim_i, elec_j, elec_i, :) - d_xi_een_theta_j_beta_d_r_i_alpha_minus(dim_j, dim_i, elec_j, elec_i, :) ) / (2.d0 * epsilon)
              end do
           end do
        end do
     end do

    !restore values
    coord_elec_wlk = coord_elec_wlk_temp
    call object_modified('coord_elec_wlk')
    call object_provide_in_node(lhere, 'd_xi_een_theta_j_beta_d_r_i_alpha')

  end subroutine lap_i_xi_een_theta_j_beta_num_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_j_theta_first_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_j_theta_first. This is a rank four array.
    !               Laplacian of the theta:
    !
    !               (nabla_j)^2 theta(r_jI,r_lI,r_jl)       Note that this is derivative with respect to first index so derivatives of the cutoff
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of electrons and labels the lth electron
    !               The third index has dimensions equal to the number of centers and labels the Ith nucleus
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_l, dim_i, cent_i
    character(len=max_string_len_rout), save :: lhere = 'lap_j_theta_first_bld'

    ! header
    if (header_exe) then

       call object_create('lap_j_theta_first')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed('ndim')
       call object_needed('all_elec')
       call object_needed('d_theta_d_r_j_inuc')
       call object_needed('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('d_theta_d_r_ji')
       call object_needed('dist_ee_wlk')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('d2_theta_d_r_ji_2')
       call object_needed('d2_theta_d_r_j_inuc_2')
       call object_needed('d2_theta_d_r_j_inuc_d_r_ji')

       return
    endif

    call object_alloc('lap_j_theta_first', lap_j_theta_first, nelec, nelec, ncent, nwalk)
    lap_j_theta_first(:,:,:,:) = 0.d0

    do elec_l = 1, nelec
       lap_j_theta_first(:, elec_l, :, :) = d2_theta_d_r_j_inuc_2(:, elec_l, :, :) + (ndim - 1.d0) * d_theta_d_r_j_inuc(:, elec_l, :, :) / dist_en_wlk(:, :, :) + d2_theta_d_r_ji_2(:, elec_l, :, :)
    end do

    do cent_i = 1, ncent
       do elec_l = 1, nelec
          lap_j_theta_first(1:elec_l - 1, elec_l, cent_i, :) = lap_j_theta_first(1:elec_l - 1, elec_l, cent_i, :) + (ndim - 1.d0) * d_theta_d_r_ji(1:elec_l - 1, elec_l, cent_i, :) / dist_ee_wlk(1:elec_l - 1, elec_l, :) +  2.d0 * d2_theta_d_r_j_inuc_d_r_ji(1:elec_l - 1, elec_l, cent_i, :) / dist_en_wlk(1:elec_l - 1, cent_i, :) / dist_ee_wlk(1:elec_l - 1, elec_l, :) * sum(vec_en_xyz_wlk(:, 1:elec_l - 1, cent_i, :) * vec_ee_xyz_wlk(:, 1:elec_l - 1, elec_l, :), 1)
          lap_j_theta_first(elec_l + 1:nelec, elec_l, cent_i, :) = lap_j_theta_first(elec_l + 1:nelec, elec_l, cent_i, :) + (ndim - 1.d0) * d_theta_d_r_ji(elec_l + 1:nelec, elec_l, cent_i, :) / dist_ee_wlk(elec_l + 1:nelec, elec_l, :) + 2.d0 * d2_theta_d_r_j_inuc_d_r_ji(elec_l + 1:nelec, elec_l, cent_i, :) / dist_en_wlk(elec_l + 1:nelec, cent_i, :) / dist_ee_wlk(elec_l + 1:nelec, elec_l, :) * sum(vec_en_xyz_wlk(:, elec_l + 1:nelec, cent_i, :) * vec_ee_xyz_wlk(:, elec_l + 1:nelec, elec_l, :), 1)
       end do
    end do

    if (all_elec) then
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       call object_provide_in_node(lhere, 'd_smooth_cutoff_g_d_r')
       call object_provide_in_node(lhere, 'd2_smooth_cutoff_g_d_r_2')
       call object_provide_in_node(lhere, 'theta')
       call object_provide_in_node(lhere, 'd_theta_d_r_j_alpha_first')

       do cent_i = 1, ncent
          do elec_l = 1, nelec
             do dim_i = 1, ndim
                lap_j_theta_first(1:elec_l - 1, elec_l, cent_i, :) = lap_j_theta_first(1:elec_l - 1, elec_l, cent_i, :) + d_theta_d_r_ji(1:elec_l - 1, elec_l, cent_i, :) / dist_ee_wlk(1:elec_l - 1, elec_l, :) * vec_ee_xyz_wlk(dim_i, 1:elec_l - 1, elec_l, :) * ( sum(d_smooth_cutoff_g_d_r(1:elec_l - 1, 1:cent_i - 1, :) / smooth_cutoff_g(1:elec_l - 1, 1:cent_i - 1, :)  * vec_en_xyz_wlk(dim_i, 1:elec_l - 1, 1:cent_i - 1, :) / dist_en_wlk(1:elec_l - 1, 1:cent_i - 1, :), 2) + sum(d_smooth_cutoff_g_d_r(1:elec_l - 1, cent_i + 1:ncent, :) / smooth_cutoff_g(1:elec_l - 1, cent_i + 1:ncent, :)  * vec_en_xyz_wlk(dim_i, 1:elec_l - 1, cent_i + 1:ncent, :) / dist_en_wlk(1:elec_l - 1, cent_i + 1:ncent, :), 2) )
                lap_j_theta_first(elec_l + 1:nelec, elec_l, cent_i, :) = lap_j_theta_first(elec_l + 1:nelec, elec_l, cent_i, :) + d_theta_d_r_ji(elec_l + 1:nelec, elec_l, cent_i, :) / dist_ee_wlk(elec_l + 1:nelec, elec_l, :) * vec_ee_xyz_wlk(dim_i, elec_l + 1:nelec, elec_l, :) * ( sum(d_smooth_cutoff_g_d_r(elec_l + 1:nelec, 1:cent_i - 1, :) / smooth_cutoff_g(elec_l + 1:nelec, 1:cent_i - 1, :)  * vec_en_xyz_wlk(dim_i, elec_l + 1:nelec, 1:cent_i - 1, :) / dist_en_wlk(elec_l + 1:nelec, 1:cent_i - 1, :), 2) + sum(d_smooth_cutoff_g_d_r(elec_l + 1:nelec, cent_i + 1:ncent, :) / smooth_cutoff_g(elec_l + 1:nelec, cent_i + 1:ncent, :)  * vec_en_xyz_wlk(dim_i, elec_l + 1:nelec, cent_i + 1:ncent, :) / dist_en_wlk(elec_l + 1:nelec, cent_i + 1:ncent, :), 2) )
                lap_j_theta_first(:, elec_l, cent_i, :) = lap_j_theta_first(:, elec_l, cent_i, :) + d_theta_d_r_j_inuc(:, elec_l, cent_i, :) / dist_en_wlk(:, cent_i, :) * vec_en_xyz_wlk(dim_i, :, cent_i, :) * ( sum(d_smooth_cutoff_g_d_r(:, 1:cent_i - 1, :) / smooth_cutoff_g(:, 1:cent_i - 1, :)  * vec_en_xyz_wlk(dim_i, :, 1:cent_i - 1, :) / dist_en_wlk(:, 1:cent_i - 1, :), 2) + sum(d_smooth_cutoff_g_d_r(:, cent_i + 1:ncent, :) / smooth_cutoff_g(:, cent_i + 1:ncent, :)  * vec_en_xyz_wlk(dim_i, :, cent_i + 1:ncent, :) / dist_en_wlk(:, cent_i + 1:ncent, :), 2) )
             end do
          end do
       end do

       do cent_i = 1, ncent
          do elec_l = 1, nelec
             lap_j_theta_first(:, elec_l, cent_i, :) = lap_j_theta_first(:, elec_l, cent_i, :) * product(smooth_cutoff_g(:, 1:cent_i - 1, :), 2) * product(smooth_cutoff_g(:, cent_i + 1:ncent, :), 2)
          end do
       end do

       do cent_i = 1, ncent
          do elec_l = 1, nelec
             lap_j_theta_first(:, elec_l, cent_i, :) = lap_j_theta_first(:, elec_l, cent_i, :) + theta(:, elec_l, cent_i, :) * ( sum( (d2_smooth_cutoff_g_d_r_2(:, 1:cent_i - 1, :) - d_smooth_cutoff_g_d_r(:, 1:cent_i - 1, :) ** 2 / smooth_cutoff_g(:, 1:cent_i - 1, :) + d_smooth_cutoff_g_d_r(:, 1:cent_i - 1, :) * (ndim - 1.d0) / dist_en_wlk(:, 1:cent_i - 1, :) ) / smooth_cutoff_g(:, 1:cent_i - 1, :), 2) + sum( (d2_smooth_cutoff_g_d_r_2(:, cent_i + 1:ncent, :) - d_smooth_cutoff_g_d_r(:, cent_i + 1:ncent, :) ** 2 / smooth_cutoff_g(:, cent_i + 1:ncent, :) + d_smooth_cutoff_g_d_r(:, cent_i + 1:ncent, :) * (ndim - 1.d0) / dist_en_wlk(:, cent_i + 1:ncent, :) ) / smooth_cutoff_g(:, cent_i + 1:ncent, :), 2) )
             do dim_i = 1, ndim
                lap_j_theta_first(:, elec_l, cent_i, :) = lap_j_theta_first(:, elec_l, cent_i, :) + d_theta_d_r_j_alpha_first(dim_i, :, elec_l, cent_i, :) * ( sum(d_smooth_cutoff_g_d_r(:, 1:cent_i - 1, :) / smooth_cutoff_g(:, 1:cent_i - 1, :)  * vec_en_xyz_wlk(dim_i, :, 1:cent_i - 1, :) / dist_en_wlk(:, 1:cent_i - 1, :), 2) + sum(d_smooth_cutoff_g_d_r(:, cent_i + 1:ncent, :) / smooth_cutoff_g(:, cent_i + 1:ncent, :)  * vec_en_xyz_wlk(dim_i, :, cent_i + 1:ncent, :) / dist_en_wlk(:, cent_i + 1:ncent, :), 2) )
             end do
          end do
       end do

    end if

  end subroutine lap_j_theta_first_bld

  !==============================================================================================

  !==============================================================================================

  subroutine lap_i_theta_second_bld

    !---------------------------------------------------------------------------
    ! Description : build lap_i_theta_second. This is a rank four array.
    !               Laplacian of the theta:
    !
    !               (nabla_i)^2 theta(r_jI,r_iI,r_ji)       Note that this is derivative with respect to second index so no derivatives of the cutoff
    !
    !               The first index has dimensions equal to the number of electrons and labels the jth electron
    !               The second index has dimensions equal to the number of electrons and labels the ith electron
    !               The third index has dimensions equal to the number of centers and labels the Ith nucleus
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_j, cent_i
    character(len=max_string_len_rout), save :: lhere = 'lap_i_theta_second_bld'

    ! header
    if (header_exe) then

       call object_create('lap_i_theta_second')

       call object_needed('nwalk')
       call object_needed('nelec')
       call object_needed('ncent')
       call object_needed('ndim')
       call object_needed('all_elec')
       call object_needed('d_theta_d_r_j_inuc')
       call object_needed('dist_en_wlk')
       call object_needed('vec_en_xyz_wlk')
       call object_needed('d_theta_d_r_ji')
       call object_needed('dist_ee_wlk')
       call object_needed('vec_ee_xyz_wlk')
       call object_needed('d2_theta_d_r_ji_2')
       call object_needed('d2_theta_d_r_j_inuc_2')
       call object_needed('d2_theta_d_r_j_inuc_d_r_ji')

       return
    endif

    call object_alloc('lap_i_theta_second', lap_i_theta_second, nelec, nelec, ncent, nwalk)
    lap_i_theta_second(:,:,:,:) = 0.d0

    do elec_j = 1, nelec
       lap_i_theta_second(elec_j, :, :, :) = d2_theta_d_r_j_inuc_2(:, elec_j, :, :) + (ndim - 1.d0) * d_theta_d_r_j_inuc(:, elec_j, :, :) / dist_en_wlk(:, :, :) + d2_theta_d_r_ji_2(:, elec_j, :, :)
    end do

    do cent_i = 1, ncent
       do elec_j = 1, nelec
          lap_i_theta_second(elec_j, 1:elec_j - 1, cent_i, :) = lap_i_theta_second(elec_j, 1:elec_j - 1, cent_i, :) + (ndim - 1.d0) * d_theta_d_r_ji(1:elec_j - 1, elec_j, cent_i, :) / dist_ee_wlk(1:elec_j - 1, elec_j, :) +  2.d0 * d2_theta_d_r_j_inuc_d_r_ji(1:elec_j - 1, elec_j, cent_i, :) / dist_en_wlk(1:elec_j - 1, cent_i, :) / dist_ee_wlk(1:elec_j - 1, elec_j, :) * sum(vec_en_xyz_wlk(:, 1:elec_j - 1, cent_i, :) * vec_ee_xyz_wlk(:, 1:elec_j - 1, elec_j, :), 1)
          lap_i_theta_second(elec_j, elec_j + 1:nelec, cent_i, :) = lap_i_theta_second(elec_j, elec_j + 1:nelec, cent_i, :) + (ndim - 1.d0) * d_theta_d_r_ji(elec_j + 1:nelec, elec_j, cent_i, :) / dist_ee_wlk(elec_j + 1:nelec, elec_j, :) +  2.d0 * d2_theta_d_r_j_inuc_d_r_ji(elec_j + 1:nelec, elec_j, cent_i, :) / dist_en_wlk(elec_j + 1:nelec, cent_i, :) / dist_ee_wlk(elec_j + 1:nelec, elec_j, :) * sum(vec_en_xyz_wlk(:, elec_j + 1:nelec, cent_i, :) * vec_ee_xyz_wlk(:, elec_j + 1:nelec, elec_j, :), 1)
       end do
    end do

    if (all_elec) then
       call object_provide_in_node(lhere, 'smooth_cutoff_g')
       !cutoff
       do cent_i = 1, ncent
          do elec_j = 1, nelec
             lap_i_theta_second(:, elec_j, cent_i, :) = lap_i_theta_second(:, elec_j, cent_i, :) * product(smooth_cutoff_g(:, 1:cent_i - 1, :), 2) * product(smooth_cutoff_g(:, cent_i + 1:ncent, :), 2)
          end do
       end do
    end if

  end subroutine lap_i_theta_second_bld

  !==============================================================================================

  !======================================================================================

  subroutine d2_theta_d_r_ji_2_bld
    !---------------------------------------------------------------------------
    ! Description : build object d2_theta_d_r_ji_2(:,:,:,:). This is a rank four array.
    !
    !               d^2 theta(r_jI,r_iI,r_ji)    (without cutoff)
    !               ----------
    !               d r_(ji)^2
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 26 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d2_theta_d_r_ji_2')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_theta_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_ee_wlk_d_r')
       call object_needed ('d2_scaled_dist_een_ee_wlk_d_r_2')
       call object_needed('b_param_theta')

       return
    endif

    call object_alloc('d2_theta_d_r_ji_2', d2_theta_d_r_ji_2, nelec, nelec, ncent, nwalk)
    d2_theta_d_r_ji_2 (:,:,:,:) = 0.d0

    param_i = 1 !intialize parameter counter
    do order_p = 2, order_theta_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d2_theta_d_r_ji_2(elec_j, elec_i, cent_i, :) = d2_theta_d_r_ji_2(elec_j, elec_i, cent_i, :) + order_k * b_param_theta(param_i, iwctype(cent_i), up_up) * ( (order_k - 1) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-2) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) ** 2 + scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) * d2_scaled_dist_een_ee_wlk_d_r_2(elec_j, elec_i, :) ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l + scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_j = nup + 1, nelec
                         d2_theta_d_r_ji_2(elec_j, elec_i, cent_i, :) = d2_theta_d_r_ji_2(elec_j, elec_i, cent_i, :) + order_k * b_param_theta(param_i, iwctype(cent_i), up_down) * ( (order_k - 1) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-2) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) ** 2 + scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) * d2_scaled_dist_een_ee_wlk_d_r_2(elec_j, elec_i, :) ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l + scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d2_theta_d_r_ji_2(elec_j, elec_i, cent_i, :) = d2_theta_d_r_ji_2(elec_j, elec_i, cent_i, :) + order_k * b_param_theta(param_i, iwctype(cent_i), up_down) * ( (order_k - 1) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-2) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) ** 2 + scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) * d2_scaled_dist_een_ee_wlk_d_r_2(elec_j, elec_i, :) ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l + scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                      do elec_j =  nup + 1, nelec
                         d2_theta_d_r_ji_2(elec_j, elec_i, cent_i, :) = d2_theta_d_r_ji_2(elec_j, elec_i, cent_i, :) + order_k * b_param_theta(param_i, iwctype(cent_i), down_down) * ( (order_k - 1) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-2) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) ** 2 + scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k-1) * d2_scaled_dist_een_ee_wlk_d_r_2(elec_j, elec_i, :) ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :)**order_l + scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l ) * ( scaled_dist_een_en_wlk(elec_j, cent_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :) ) ** ((order_p - order_k - order_l) / 2)
                      end do
                   end do
                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d2_theta_d_r_ji_2(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d2_theta_d_r_ji_2_bld

  !=============================================================================================

  !======================================================================================

  subroutine d2_theta_d_r_j_inuc_d_r_ji_bld
    !---------------------------------------------------------------------------
    ! Description : build object d2_theta_d_r_j_inuc_d_r_ji(:,:,:,:). This is a rank four array.
    !
    !               d^2 theta(r_jI,r_iI,r_jI)    (without cutoff)
    !               ----------
    !               d r_jI d r_ji
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 26 Mar 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d2_theta_d_r_j_inuc_d_r_ji')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_theta_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_en_wlk_d_r')
       call object_needed ('d_scaled_dist_een_ee_wlk_d_r')
       call object_needed('b_param_theta')

       return
    endif

    call object_alloc('d2_theta_d_r_j_inuc_d_r_ji', d2_theta_d_r_j_inuc_d_r_ji, nelec, nelec, ncent, nwalk)
    d2_theta_d_r_j_inuc_d_r_ji (:,:,:,:) = 0.d0

    param_i = 1 !intialize parameter counter
    do order_p = 2, order_theta_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d2_theta_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) = d2_theta_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_up) * order_k * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k -1) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2)  * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                      do elec_j = nup+1, nelec
                         d2_theta_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) = d2_theta_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_down) * order_k * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k -1) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2)  * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d2_theta_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) = d2_theta_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_down) * order_k * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k -1) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2)  * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                      do elec_j = nup+1, nelec
                         d2_theta_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) = d2_theta_d_r_j_inuc_d_r_ji(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), down_down) * order_k * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**(order_k -1) * d_scaled_dist_een_ee_wlk_d_r(elec_j, elec_i, :) *  scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2)  * d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2) + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2))
                      end do
                   end do
                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d2_theta_d_r_j_inuc_d_r_ji(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d2_theta_d_r_j_inuc_d_r_ji_bld

  !=============================================================================================

  !======================================================================================

  subroutine d2_theta_d_r_j_inuc_2_bld
    !---------------------------------------------------------------------------
    ! Description : build object d2_theta_d_r_j_inuc_2(:,:,:,:). This is a rank four array.
    !
    !               d theta(r_jI,r_iI,r_jI)    (without cutoff)
    !               ----------
    !               d r_jI
    !
    !               The first index labels j and has dimension equal to the number of electrons.
    !               The second index labels i and has dimension equal to the number of electrons.
    !               The third index labels I and has dimensions equal to the number of nuclei.
    !               The fourth index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 25 Feb 2009
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_j, elec_i, cent_i, order_p, order_k, order_l, param_i
    integer  :: order_l_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('d2_theta_d_r_j_inuc_2')

       call object_needed('nwalk')
       call object_needed('ncent')
       call object_needed('iwctype')
       call object_needed('nelec')
       call object_needed('nup')
       call object_needed('order_theta_bf')
       call object_needed ('scaled_dist_een_en_wlk')
       call object_needed ('scaled_dist_een_ee_wlk')
       call object_needed ('d_scaled_dist_een_en_wlk_d_r')
       call object_needed ('d2_scaled_dist_een_en_wlk_d_r_2')
       call object_needed('b_param_theta')

       return
    endif

    call object_alloc('d2_theta_d_r_j_inuc_2', d2_theta_d_r_j_inuc_2, nelec, nelec, ncent, nwalk)
    d2_theta_d_r_j_inuc_2 (:,:,:,:) = 0.d0

    param_i = 1 !intialize parameter counter
    do order_p = 2, order_theta_bf
       do order_k = 0, order_p - 1
          if (order_k .ne. 0) then
             order_l_max = order_p - order_k
          else
             order_l_max = order_p - 2
          end if
          do order_l = 0, order_l_max
             if (modulo(order_p - order_k - order_l, 2) .ne. 0) then
                cycle
             else
                do cent_i = 1, ncent
                   do elec_i = 1, nup
                      do elec_j = 1, nup
                         d2_theta_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) = d2_theta_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_up) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k * scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) * ( d2_scaled_dist_een_en_wlk_d_r_2(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2)  + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2)) + d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) **2 * ( (order_p - order_k + order_l) / 2 * (order_p - order_k + order_l - 2) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-4)/2)  + (order_p - order_k - order_l) / 2  * (order_p - order_k - order_l - 2) / 2 * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-4)/2)) )
                      end do
                      do elec_j = nup+1, nelec
                         d2_theta_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) = d2_theta_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_down) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k * scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) * ( d2_scaled_dist_een_en_wlk_d_r_2(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2)  + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2)) + d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) **2 * ( (order_p - order_k + order_l) / 2 * (order_p - order_k + order_l - 2) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-4)/2)  + (order_p - order_k - order_l) / 2  * (order_p - order_k - order_l - 2) / 2 * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-4)/2)) )
                      end do
                   end do
                   do elec_i = nup + 1, nelec
                      do elec_j = 1, nup
                         d2_theta_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) = d2_theta_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), up_down) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k * scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) * ( d2_scaled_dist_een_en_wlk_d_r_2(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2)  + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2)) + d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) **2 * ( (order_p - order_k + order_l) / 2 * (order_p - order_k + order_l - 2) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-4)/2)  + (order_p - order_k - order_l) / 2  * (order_p - order_k - order_l - 2) / 2 * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-4)/2)) )
                      end do
                      do elec_j = nup+1, nelec
                         d2_theta_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) = d2_theta_d_r_j_inuc_2(elec_j, elec_i, cent_i, :) + b_param_theta(param_i, iwctype(cent_i), down_down) * scaled_dist_een_ee_wlk(elec_j, elec_i, :)**order_k * scaled_dist_een_en_wlk(elec_i, cent_i, :)  ** ((order_p - order_k - order_l) / 2) * ( d2_scaled_dist_een_en_wlk_d_r_2(elec_j, cent_i, :) * ( (order_p - order_k + order_l) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-2)/2)  + (order_p - order_k - order_l) / 2  * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-2)/2)) + d_scaled_dist_een_en_wlk_d_r(elec_j, cent_i, :) **2 * ( (order_p - order_k + order_l) / 2 * (order_p - order_k + order_l - 2) / 2  * scaled_dist_een_en_wlk(elec_j, cent_i, :)**((order_p - order_k + order_l-4)/2)  + (order_p - order_k - order_l) / 2  * (order_p - order_k - order_l - 2) / 2 * scaled_dist_een_en_wlk(elec_i, cent_i, :)**order_l * scaled_dist_een_en_wlk(elec_j, cent_i, :) ** ((order_p - order_k - order_l-4)/2)) )
                      end do
                   end do
                end do
                param_i = param_i + 1
             endif
          enddo
       enddo
    enddo

    !no self backflow
    do elec_i=1, nelec
       d2_theta_d_r_j_inuc_2(elec_i, elec_i, :, :) = 0.d0
    end do

  end subroutine d2_theta_d_r_j_inuc_2_bld

  !=============================================================================================

end module backflow_mod
