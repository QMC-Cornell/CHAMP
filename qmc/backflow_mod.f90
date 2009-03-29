module backflow_mod
  use parser_tools_mod, only: word, get_next_word, get_next_value, get_next_value_list_object
  use constants_mod    !give access to various common blocks
  use basic_tools_mod, only: require, die
  use strings_tools_mod !require operator overloading of + to concatenate strings
  use objects_mod, only: object_modified, object_create, object_needed, object_alloc, object_provide, object_associate, object_provide_in_node
  use nodes_mod, only: header_exe
  use variables_mod, only: spin_nb
  use electrons_mod, only: dist_en_wlk, vec_en_xyz_wlk, vec_ee_xyz_wlk, dist_ee_wlk

  !Declaration of global variables and default values

  integer                                        :: threshold_l_nb
  real(dp), allocatable                          :: threshold_l(:)
  real(dp), allocatable                          :: smooth_cutoff_g(:,:,:)

  integer                                        :: isc_bf
  real(dp)                                       :: scalek_bf

  logical                                        :: all_elec

  logical                                        :: en_bf
  real(dp), allocatable                          :: xi_en(:,:,:)
  integer                                        :: order_mu_bf
  integer                                        :: read_d_param_mu_nb
  logical                                        :: read_d_param
  real(dp), allocatable                          :: read_d_param_mu(:)
  real(dp), allocatable                          :: mu_pseudo(:,:,:)
  real(dp), allocatable                          :: d_param_mu_pseudo(:,:,:)
  real(dp), allocatable                          :: mu_all_elec(:,:,:)
  real(dp), allocatable                          :: d_param_mu_all_elec(:,:,:)
  logical                                        :: mu_spin_dependent
  real(dp), allocatable                          :: asymp_mu_pseudo(:,:)

  logical                                        :: ee_bf
  real(dp), allocatable                          :: xi_ee(:,:,:)
  integer                                        :: order_eta_bf
  integer                                        :: read_c_param_eta_nb
  logical                                        :: read_c_param
  real(dp), allocatable                          :: read_c_param_eta(:)
  real(dp), allocatable                          :: eta_pseudo(:,:,:)
  real(dp), allocatable                          :: eta_all_elec(:,:,:)
  real(dp), allocatable                          :: c_param_eta(:,:)
  integer                                        :: eta_spin_dependence
  real(dp), allocatable                          :: asymp_eta(:)

  logical                                        :: een_phi_bf
  real(dp), allocatable                          :: xi_een_phi(:,:,:)
  integer                                        :: order_phi_bf
  integer                                        :: read_a_param_phi_nb
  logical                                        :: read_a_param
  real(dp), allocatable                          :: read_a_param_phi(:)
  real(dp), allocatable                          :: phi_pseudo(:,:,:,:)
  real(dp), allocatable                          :: phi_all_elec(:,:,:,:)
  real(dp), allocatable                          :: a_param_phi_all_elec(:,:,:)
  real(dp), allocatable                          :: a_param_phi_pseudo(:,:,:)
  real(dp), allocatable                          :: a_param_phi_pseudo_anti_parallel_cond(:,:)
  real(dp), allocatable                          :: a_param_phi_pseudo_parallel_cond(:,:)
  real(dp), allocatable                          :: a_param_phi_all_elec_parallel_cond(:,:)
  real(dp), allocatable                          :: a_param_phi_all_elec_anti_parallel_cond(:,:)
  integer                                        :: phi_spin_dependence
  integer                                        :: nb_a_param_phi
  integer, allocatable                           :: dep_a_param_phi_pseudo_parallel(:)
  integer, allocatable                           :: dep_a_param_phi_pseudo_anti_parallel(:)
  integer, allocatable                           :: dep_a_param_phi_all_elec_parallel(:)
  integer, allocatable                           :: dep_a_param_phi_all_elec_anti_parallel(:)

  logical                                        :: een_theta_bf
  real(dp), allocatable                          :: xi_een_theta(:,:,:)
  integer                                        :: order_theta_bf
  integer                                        :: read_b_param_theta_nb
  logical                                        :: read_b_param
  real(dp), allocatable                          :: read_b_param_theta(:)
  real(dp), allocatable                          :: theta_pseudo(:,:,:,:)
  real(dp), allocatable                          :: theta_all_elec(:,:,:,:)
  real(dp), allocatable                          :: b_param_theta_all_elec(:,:,:)
  real(dp), allocatable                          :: b_param_theta_pseudo(:,:,:)
  real(dp), allocatable                          :: b_param_theta_pseudo_cond(:,:)
  real(dp), allocatable                          :: b_param_theta_all_elec_cond(:,:)
  integer                                        :: theta_spin_dependence
  integer                                        :: nb_b_param_theta
  integer, allocatable                           :: dep_b_param_theta_pseudo(:)
  integer, allocatable                           :: dep_b_param_theta_all_elec(:)

  real(dp), allocatable                          :: scaled_dist_een_ee_wlk(:,:,:)
  real(dp), allocatable                          :: scaled_dist_een_en_wlk(:,:,:)
  real(dp), allocatable                          :: scaled_dist_ee_wlk(:,:,:)
  real(dp), allocatable                          :: scaled_dist_en_wlk(:,:,:)
  real(dp)                                       :: asymp_scaled_dist_two_body

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
    !object associate will add a pointer to this object so that it knows its name.
    !Note that this is called automatically when you call object_alloc for an array.
    !When we have a real number we must call explicitly

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

    mu_spin_dependent = .true.                    !  up != down
    call object_modified('mu_spin_dependent')
    call object_associate('mu_spin_dependent', mu_spin_dependent)

    eta_spin_dependence = 3                       !  up_up != down_down != up_down
    call object_modified('eta_spin_dependence')
    call object_associate('eta_spin_dependence', eta_spin_dependence)

    phi_spin_dependence = 3                       !  up_up != down_down != up_down
    call object_modified('phi_spin_dependence')
    call object_associate('phi_spin_dependence', phi_spin_dependence)

    theta_spin_dependence = 3                     !  up_up != down_down != up_down
    call object_modified('theta_spin_dependence')
    call object_associate('theta_spin_dependence', theta_spin_dependence)

    !by default we assume that the backflow parameters are not read in. However the parameters can be read in. See menu.
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

    !by default we assume that the b parameters for theta are not read in. However the parameters can be read in. See menu.
    read_b_param=.false.
    call object_modified('read_b_param')
    call object_associate('read_b_param', read_b_param)

    read_b_param_theta_nb = 0
    call object_modified('read_b_param_theta_nb')
    call object_associate('read_b_param_theta_nb', read_b_param_theta_nb)

    !set default values for scalek_bf, isc_bf
    !by default use the jastrow values
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

       !Find next word in input file
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
          !read in order_eta_bf
          call get_next_value(order_eta_bf)
          !check that order_eta_bf is non-negative
          call require(lhere, 'order_eta_bf > 0', order_eta_bf .gt. 0)
          !Want order_eta_bf as part of tree.
          !Note if object_needed('order_eta_bf') is invoked without the object_modified() having been executed the program will stop.
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
          !read in list of parameters. Note that read_d_param_mu_nb returns the length of the list
          call get_next_value_list_object('read_d_param_mu', read_d_param_mu, read_d_param_mu_nb)
          call object_modified('read_d_param_mu')
          call object_modified('read_d_param_mu_nb')
          !d parameters were read in
          read_d_param=.true.
          call object_modified('read_d_param')


       case ('read_c_param_eta')
          !read in list of parameters. Note that read_c_param_eta_nb returns the length of the list
          call get_next_value_list_object('read_c_param_eta', read_c_param_eta, read_c_param_eta_nb)
          call object_modified('read_c_param_eta')
          call object_modified('read_c_param_eta_nb')
          !c parameters were read in
          read_c_param=.true.
          call object_modified('read_c_param')

       case ('read_a_param_phi')
          !read in list of parameters. Note that read_a_param_phi_nb returns the length of the list
          call get_next_value_list_object('read_a_param_phi', read_a_param_phi, read_a_param_phi_nb)
          call object_modified('read_a_param_phi')
          call object_modified('read_a_param_phi_nb')
          !a parameters were read in
          read_a_param=.true.
          call object_modified('read_a_param')

       case ('read_b_param_theta')
          !read in list of parameters. Note that read_b_param_theta_nb returns the length of the list
          call get_next_value_list_object('read_b_param_theta', read_b_param_theta, read_b_param_theta_nb)
          call object_modified('read_b_param_theta')
          call object_modified('read_b_param_theta_nb')
          !b parameters were read in
          read_b_param=.true.
          call object_modified('read_b_param')

       case ('threshold_l')
          call get_next_value_list_object('threshold_l', threshold_l, threshold_l_nb)
          call object_modified('threshold_l')
          call object_modified('threshold_l_nb')
          call object_associate('threshold_l_nb', threshold_l_nb)
          call object_provide('nctype')
          !specify l value for each type of nucleus
          call require(lhere, 'Size of threshold_l = nctype', threshold_l_nb .eq. nctype)

       case ('end')
          !done reading
          exit

       case default
          !do not recognize word
          !note that + operator is overloaded in strings_tools_mod.f90
          call die (lhere, 'unknown keyword >'+trim(word)+'<.')

       end select

    enddo ! end loop over menu lines

!testing     !testing will need to add ee een
!     !make sure we are doing backflow
!     call require(lhere, 'If you have a backflow menu you should be doing backflow',en_bf) ee_bf somewhere here too

    write(6,'(a)') 'End of backflow menu ---------------------------------------------------------------------------------'

!testing     !testing
!     ! launch backflow
!     ! need a flag here


  end subroutine backflow_menu

!======================================================================================

!======================================================================================

  subroutine xi_en_bld
    !---------------------------------------------------------------------------
    ! Description : build object xi_en(:,:). This is a rank three array.
    !               It is the electron-nuclear backflow.
    !               It is necessary for the construction of the total backflow transformation.
    !               The first index has dimension equal to the number of spatial dimensions.
    !               The second index has dimensions equal to the number of electrons.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 21 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, cent_i, walk_i,dim_i
    character(len=max_string_len_rout), save :: lhere = 'xi_en_bld'

    ! header
    if (header_exe) then

       call object_create('xi_en')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       ! number of electrons
       call object_needed('nelec')
       !cartesian vector between electrons and nuclei (for each walker)
       call object_needed('vec_en_xyz_wlk')
       !dimension of space
       call object_needed('ndim')
       !is this an all-electron calculation?
       call object_needed('all_elec')

       return
    endif

    !Note that we need to have the option of either all electron or pseudo.
    !However, we want only one object xi_en. Thus we don't call object_needed('mu_all_elec')
    !and  object_needed('mu_pseudo') for instance.
    ! This would result in evaluating both everytime we want xi_en.
    !To accomplish this correctly we need to use the utility object_provide_in_node
    !We cant simply use object_provide because we are inside the tree.

    if (all_elec) then
       !all electron case
       call object_provide_in_node(lhere, 'mu_all_elec')

       !allocation

       call object_alloc('xi_en', xi_en, ndim, nelec, nwalk)

       !build xi_en
       xi_en(:,:,:) = 0.d0

       do walk_i = 1, nwalk
          do cent_i = 1, ncent
             do elec_i = 1, nelec  !note that the first nup electrons have spin-up
                do dim_i = 1, ndim
                   xi_en(dim_i, elec_i, walk_i) = xi_en(dim_i, elec_i, walk_i) + mu_all_elec(elec_i, cent_i, walk_i) * vec_en_xyz_wlk(dim_i, elec_i, cent_i, walk_i)
                enddo
             enddo
          enddo
       enddo
    else
       !pseudo case
       call object_provide_in_node(lhere, 'mu_pseudo')
       !allocation

       call object_alloc('xi_en', xi_en, ndim, nelec, nwalk)

       !build xi_en
       xi_en(:,:,:) = 0.d0

       do walk_i = 1, nwalk
          do cent_i = 1, ncent
             do elec_i = 1, nelec  !note that the first nup electrons have spin-up
                do dim_i = 1, ndim
                   xi_en(dim_i, elec_i, walk_i) = xi_en(dim_i, elec_i, walk_i) + mu_pseudo(elec_i, cent_i, walk_i) * vec_en_xyz_wlk(dim_i, elec_i, cent_i, walk_i)
                enddo
             enddo
          enddo
       enddo
    endif


  end subroutine xi_en_bld

  !======================================================================================

  !======================================================================================

  subroutine mu_all_elec_bld
    !---------------------------------------------------------------------------
    ! Description : build object mu_all_elec(:,:,:). This is a rank three array.
    !               It is necessary for the construction of the electron-nuclear backflow terms in all electron case.
    !               The first index has dimension equal to the number of electrons.
    !               The second index has dimensions equal to the number of nuclei.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 20 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_i, cent_i, walk_i, order_i, prod_i
    integer, parameter  :: up=1, down=2 !index for the spin

    ! header
    if (header_exe) then

       call object_create('mu_all_elec')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       !specify atom-type for each atom
       call object_needed('iwctype')
       ! number of electrons
       call object_needed('nelec')
       !number of up-spin electrons
       call object_needed('nup')
       !order of expansion
       call object_needed('order_mu_bf')
       !scaled distances of each electron of each walker to each nucleus
       call object_needed ('scaled_dist_en_wlk')
       !list of d parameters (satisfying all-elec conditions) for mu backflow function for each center type, and each spin
       call object_needed('d_param_mu_all_elec')
       !In all electron calculations, to impose the electron nuclear cusp, it is necessary to smoothly cutoff the backflow functions at the nucleus.
       call object_needed('smooth_cutoff_g')

       return
    endif

    ! allocation

    !this can be made more efficient
    call object_alloc('mu_all_elec', mu_all_elec, nelec, ncent, nwalk)

    ! build mu_all_elec

    !initialize
    mu_all_elec (:,:,:) = 0.d0

    !spin-up electrons
    !fortran arrays are in column-major storage order. So, left most index is fastest, use in the innermost loop.
    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !pade term of mu_all_elec
             mu_all_elec(elec_i, cent_i, walk_i) = d_param_mu_all_elec(1, iwctype(cent_i), up) * scaled_dist_en_wlk(elec_i, cent_i, walk_i) / (1.d0 + d_param_mu_all_elec(2, iwctype(cent_i), up) * scaled_dist_en_wlk(elec_i, cent_i, walk_i))
             !power series terms of mu_all_elec
             do order_i=2, order_mu_bf
                mu_all_elec(elec_i, cent_i, walk_i) = mu_all_elec(elec_i, cent_i, walk_i) + d_param_mu_all_elec(order_i + 1, iwctype(cent_i), up) * scaled_dist_en_wlk(elec_i, cent_i, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo


    !spin-down electrons
    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons: nup+1 to nelec
          do elec_i = nup+1, nelec
             !pade term of mu_all_elec
             mu_all_elec(elec_i, cent_i, walk_i) = d_param_mu_all_elec(1, iwctype(cent_i), down) * scaled_dist_en_wlk(elec_i, cent_i, walk_i) / (1.d0 + d_param_mu_all_elec(2, iwctype(cent_i), down) * scaled_dist_en_wlk(elec_i, cent_i, walk_i))
             !power series terms of mu_all_elec
             do order_i=2, order_mu_bf
                mu_all_elec(elec_i, cent_i, walk_i) = mu_all_elec(elec_i, cent_i, walk_i) + d_param_mu_all_elec(order_i + 1, iwctype(cent_i), down) * scaled_dist_en_wlk(elec_i, cent_i, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo

    !note we don't have to remove asymptotic value from mu because it is forced to be zero by satisfying the electron-nuclear cusp conditions.

    !apply smooth cutoff function g
    do cent_i=1, ncent
       do prod_i=1, cent_i-1
          mu_all_elec(:,cent_i,:) = mu_all_elec(:,cent_i,:) * smooth_cutoff_g(:,prod_i,:)
       end do
       do prod_i=cent_i+1, ncent
          mu_all_elec(:,cent_i,:) = mu_all_elec(:,cent_i,:) * smooth_cutoff_g(:,prod_i,:)
       end do
    end do

  end subroutine mu_all_elec_bld

  !=============================================================================================

  !=============================================================================================

  subroutine mu_pseudo_bld
    !---------------------------------------------------------------------------
    ! Description : build object mu_pseudo(:,:,:). This is a rank three array.
    !               It is necessary for the construction of the electron-nuclear backflow terms in pseudo case.
    !               The first index has dimension equal to the number of electrons.
    !               The second index has dimensions equal to the number of nuclei.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 21 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_i, cent_i, walk_i, order_i
    integer, parameter  :: up=1, down=2 !index for the spin

    ! header
    if (header_exe) then

       call object_create('mu_pseudo')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       !specify atom-type for each atom
       call object_needed('iwctype')
       ! number of electrons
       call object_needed('nelec')
       !number of up-spin electrons
       call object_needed('nup')
       !order of expansion
       call object_needed('order_mu_bf')
       !scaled distances of each electron of each walker to each nucleus
       call object_needed ('scaled_dist_en_wlk')
       !list of d parameters for mu backflow function for each center type, and each spin
       call object_needed('d_param_mu_pseudo')
       !asymptotic value of mu for each nucleus type, spin type
       call object_needed('asymp_mu_pseudo')

       return
    endif

    ! allocation

    !this can be made more efficient
    call object_alloc('mu_pseudo', mu_pseudo, nelec, ncent, nwalk)

    ! build mu_pseudo

    !initialize
    mu_pseudo (:,:,:) = 0.d0

    !spin-up electrons
    !fortran arrays are in column-major storage order. So, left most index is fastest, use in the innermost loop.
    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !pade term of mu_pseudo
             mu_pseudo(elec_i, cent_i, walk_i) = d_param_mu_pseudo(1, iwctype(cent_i), up) * scaled_dist_en_wlk(elec_i, cent_i, walk_i) / (1.d0 + d_param_mu_pseudo(2, iwctype(cent_i), up) * scaled_dist_en_wlk(elec_i, cent_i, walk_i))
             !power series terms of mu_pseudo
             do order_i=2, order_mu_bf
                mu_pseudo(elec_i, cent_i, walk_i) = mu_pseudo(elec_i, cent_i, walk_i) + d_param_mu_pseudo(order_i + 1, iwctype(cent_i), up) * scaled_dist_en_wlk(elec_i, cent_i, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo


    !remove asymptotic value for spin-up electrons
    !loop over number of centers
    do cent_i = 1, ncent
       mu_pseudo(1:nup,cent_i,:) = mu_pseudo(1:nup,cent_i,:) - asymp_mu_pseudo(iwctype(cent_i), up)
    enddo

    !spin-down electrons
    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons: nup+1 to nelec
          do elec_i = nup+1, nelec
             !pade term of mu_pseudo
             mu_pseudo(elec_i, cent_i, walk_i) = d_param_mu_pseudo(1, iwctype(cent_i), down) * scaled_dist_en_wlk(elec_i, cent_i, walk_i) / (1.d0 + d_param_mu_pseudo(2, iwctype(cent_i), down) * scaled_dist_en_wlk(elec_i, cent_i, walk_i))
             !power series terms of mu_pseudo
             do order_i=2, order_mu_bf
                mu_pseudo(elec_i, cent_i, walk_i) = mu_pseudo(elec_i, cent_i, walk_i) + d_param_mu_pseudo(order_i + 1, iwctype(cent_i), down) * scaled_dist_en_wlk(elec_i, cent_i, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo


    !remove asymptotic value for spin-down electrons
    !loop over number of centers
    do cent_i = 1, ncent
       mu_pseudo(nup+1:nelec,cent_i,:) = mu_pseudo(nup+1:nelec,cent_i,:) - asymp_mu_pseudo(iwctype(cent_i), down)
    enddo

  end subroutine mu_pseudo_bld

  !==============================================================================================

  !==============================================================================================


  subroutine asymp_mu_pseudo_bld
    !---------------------------------------------------------------------------
    ! Description : build object asymp_mu_pseudo(:,:). This is a rank two array.
    !               The first index is the nucleus type.
    !               The second index is the spin type.
    !               Give the asymptotic value of the pade/power series expansion used in mu.
    !               This value is subtracted off from the pade/power series expansion such that mu goes to zero at infinity.
    !
    ! Created     : F. Petruzielo, 21 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: ctype_i, spin_i, order_i

    ! header
    if (header_exe) then

       call object_create('asymp_mu_pseudo')

       !how many types of atoms
       call object_needed('nctype')
       !order of expansion
       call object_needed('order_mu_bf')
       !list of d parameters for mu_pseudo backflow function for each center type, and each spin
       call object_needed('d_param_mu_pseudo')
       !asymptotic scaled distance for electron nucleus scaling
       call object_needed('asymp_scaled_dist_two_body')

       return
    endif

    ! allocation
    !spin_nb is number of spins. Two because spin 1/2
    call object_alloc('asymp_mu_pseudo', asymp_mu_pseudo, nctype, spin_nb)

    ! build asymp_mu
    !initialize
    asymp_mu_pseudo (:,:) = 0.d0

    !calculate asymptotic values for each type a nucleus and each type of spin
    !loop over spin
    do spin_i = 1, spin_nb
       !loop over center types
       do ctype_i = 1, nctype
          !pade term of mu
          asymp_mu_pseudo(ctype_i, spin_i) = d_param_mu_pseudo(1, ctype_i, spin_i) * asymp_scaled_dist_two_body / (1.d0 + d_param_mu_pseudo(2, ctype_i, spin_i) * asymp_scaled_dist_two_body)
          !power series terms of mu
          do order_i=2, order_mu_bf
             asymp_mu_pseudo(ctype_i, spin_i) = asymp_mu_pseudo(ctype_i, spin_i) + d_param_mu_pseudo(order_i + 1, ctype_i, spin_i) * asymp_scaled_dist_two_body**order_i
          enddo
       enddo
    enddo

  end subroutine asymp_mu_pseudo_bld

  !==============================================================================================

  !==============================================================================================

  subroutine scaled_dist_en_wlk_bld

    !---------------------------------------------------------------------------
    ! Description : build scaled_dist_en_wlk(:,:,:). This is a rank three array.
    !               The first index is the electron.
    !               The second index is the nucleus
    !               The third index is the walker
    !
    ! Created     : F. Petruzielo, 18 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nelec, ncent, nwalker
    include 'commons.h'

    ! header
    if (header_exe) then

       call object_create('scaled_dist_en_wlk')

       ! number of electrons
       call object_needed('nelec')
       !number of atoms/centers
       call object_needed('ncent')
       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !distances of each electron of each walker to each nucleus
       call object_needed ('dist_en_wlk')

       !Make sure scalek_bf, isc_bf are up to date even though the values are not
       !needed directly in this build routine. They are, however, needed in a subroutine that is called
       !by this build routine.
       !scale factor for the scaling functions need to bring an up to date version into tree!!!!!
       call object_needed('scalek_bf')
       !scaling form for the scaling functions
       call object_needed('isc_bf')

       return
    endif

    ! allocation

    !this can be made more efficient
    call object_alloc('scaled_dist_en_wlk', scaled_dist_en_wlk, nelec, ncent, nwalk)

    ! build scaled_dist_en_wlk
    !initialize
    scaled_dist_en_wlk(:,:,:) = 0.d0

    !calculate scaled distances
    call scaling_func_two_body(dist_en_wlk,scaled_dist_en_wlk)

  end subroutine scaled_dist_en_wlk_bld

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

       !scale factor for the scaling functions
       call object_needed('scalek_bf')
       !scaling form for the scaling functions
       call object_needed('isc_bf')

       return
    endif



    select case(isc_bf)
    case(4) ! r / (1 + r * k)
       asymp_scaled_dist_two_body = 1.d0 / scalek_bf
    case default
       call die(lhere,'Only isc = 4 is implemented.')
    end select
    !testing
    call object_associate('asymp_scaled_dist_two_body', asymp_scaled_dist_two_body)



  end subroutine asymp_scaled_dist_two_body_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_param_mu_all_elec_bld

    !---------------------------------------------------------------------------
    ! Description : build d_param_mu_all_elec. This is a rank 3 array.
    !               List of d parameters for mu backflow function for each center type, and each spin
    !               The first index is which parameter.
    !               The second index is the type of center.
    !               The third index is the type of spin.
    !
    !               cusp condition constraints: mu(0)=0
    !
    ! Created     : F. Petruzielo, 20 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'


    !local
    integer  :: order_i, d_param_spin_nb, d_param_no_spin_nb
    character(len=max_string_len_rout), save :: lhere = 'd_param_mu_all_elec_bld'

    ! header
    if (header_exe) then

       call object_create('d_param_mu_all_elec')

       !order of expansion
       call object_needed('order_mu_bf')
       !how many types of atoms
       call object_needed('nctype')
       !asymptotic scaled distance for electron nucleus scaling
       call object_needed('asymp_scaled_dist_two_body')
       !are the parameters spin dependent?
       call object_needed('mu_spin_dependent')
       !have the parameters been read in?
       call object_needed('read_d_param')
       !how many parameters were read in?
       call object_needed('read_d_param_mu_nb')

       return
    endif



    !Note that we need to have the option of using read_d_param_mu if it is read in
    !However, we don't know apriori if it will be read in. 
    !Thus we don't call object_needed('read_d_param_mu')
    !To accomplish this correctly we need to use the utility object_provide_in_node
    !We cant simply use object_provide because we are inside the tree.

    if (read_d_param) then
       !d parameters were read in
       call object_provide_in_node(lhere, 'read_d_param_mu')
    end if



    ! allocation
    !spin_nb is number of spins. Two because spin 1/2
    call object_alloc('d_param_mu_all_elec', d_param_mu_all_elec, order_mu_bf + 1, nctype, spin_nb)

    !do we read in parameters?
    if (read_d_param) then
       !read in parameters
       d_param_spin_nb = (order_mu_bf + 1)  * nctype * spin_nb
       d_param_no_spin_nb = (order_mu_bf + 1)  * nctype
       if (read_d_param_mu_nb == d_param_spin_nb) then
          !enough parameters were read in for a spin dependent mu
          if (.not. mu_spin_dependent) then
             !too many parameters were read in
             call require(lhere, 'number of parameters in read_d_param_mu is equal to (order_mu_bf + 1) * nctype', mu_spin_dependent)
          else
             !correct number of parameters so reshape into correct dimensions for d_param_mu_all_elec
             d_param_mu_all_elec = reshape( read_d_param_mu, (/ order_mu_bf + 1, nctype, spin_nb /))
             !don't want to re-read
             read_d_param = .false.
          end if
       elseif( read_d_param_mu_nb == d_param_no_spin_nb) then
          !number of parameters read in is only sufficient for no spin dependence
          if ( mu_spin_dependent) then
             !not enough parameters read in
             call require(lhere, 'number of parameters in read_d_param_mu is equal to (order_mu_bf + 1) * nctype * spin_nb',.not. mu_spin_dependent)
          else
             !correct number of parameters so reshape into correct dimensions for d_param_mu_all_elec
             write(*,*) "debug fp"
             d_param_mu_all_elec = reshape( read_d_param_mu, (/ order_mu_bf + 1, nctype, spin_nb /), pad=read_d_param_mu)
             !note that only the first index has any meaning from the previous statement so duplicate the meaningful values

!debug             !d_param_mu_all_elec(:,:,2)=d_param_mu_all_elec(:,:,1)

             !don't want to re-read
             read_d_param = .false.
          end if
       else
          !wrong number of parameters
          call die(lhere,'read_d_param_mu has an incorrect number of parameters.')
       end if
    else
       !don't read in parameters
       !start from random guess
       call random_number(d_param_mu_all_elec(:,:,:))

       !are the parameters spin dependent?
       if (.not. mu_spin_dependent) then
          !no spin dependence
          d_param_mu_all_elec(:,:,2) = d_param_mu_all_elec(:,:,1)
       end if
    endif

    !impose electron nuclear cusp conditions
    !Note that we assume the two-body scaling functions evaluated at 0 are 0.
    d_param_mu_all_elec(1,:,:) = 0.d0
    do order_i=2, order_mu_bf
       d_param_mu_all_elec(1, :, :) = d_param_mu_all_elec(1, :, :) - d_param_mu_all_elec(order_i + 1, :, :) * asymp_scaled_dist_two_body**(order_i - 1) * ( 1.d0 +  d_param_mu_all_elec(2, :, :) * asymp_scaled_dist_two_body)
    enddo

  end subroutine d_param_mu_all_elec_bld

  !==============================================================================================

  !==============================================================================================

  subroutine d_param_mu_pseudo_bld

    !---------------------------------------------------------------------------
    ! Description : build d_param_mu_pseudo. This is a rank 3 array.
    !               List of d parameters for mu backflow function for each center type, and each spin
    !               The first index is which parameter.
    !               The second index is the type of center.
    !               The third index is the type of spin.
    !
    !               cusp condition constraints: none
    !
    ! Created     : F. Petruzielo, 21 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: d_param_spin_nb, d_param_no_spin_nb
    character(len=max_string_len_rout), save :: lhere = 'd_param_mu_pseudo_bld'

    ! header
    if (header_exe) then

       call object_create('d_param_mu_pseudo')

       !order of expansion
       call object_needed('order_mu_bf')
       !how many types of atoms
       call object_needed('nctype')
       !are the parameters spin dependent?
       call object_needed('mu_spin_dependent')
       !have the parameters been read in?
       call object_needed('read_d_param')
       !how many parameters were read in?
       call object_needed('read_d_param_mu_nb')

       return
    endif



    !Note that we need to have the option of using read_d_param_mu if it is read in
    !However, we don't know apriori if it will be read in. 
    !Thus we don't call object_needed('read_d_param_mu')
    !To accomplish this correctly we need to use the utility object_provide_in_node
    !We cant simply use object_provide because we are inside the tree.

    if (read_d_param) then
       !d parameters were read in
       call object_provide_in_node(lhere, 'read_d_param_mu')
    end if


    ! allocation
    !spin_nb is number of spins. Two because spin 1/2
    call object_alloc('d_param_mu_pseudo', d_param_mu_pseudo, order_mu_bf + 1, nctype, spin_nb)

    !do we read in parameters?
    if (read_d_param) then
       !read in parameters
       d_param_spin_nb = (order_mu_bf + 1)  * nctype * spin_nb
       d_param_no_spin_nb = (order_mu_bf + 1)  * nctype
       if (read_d_param_mu_nb == d_param_spin_nb) then
          !enough parameters were read in for a spin dependent mu
          if (.not. mu_spin_dependent) then
             !too many parameters were read in
             call require(lhere, 'number of parameters in read_d_param_mu is equal to (order_mu_bf + 1) * nctype', mu_spin_dependent)
          else
             !correct number of parameters so reshape into correct dimensions for d_param_mu_pseudo
             d_param_mu_pseudo = reshape( read_d_param_mu, (/ order_mu_bf + 1, nctype, spin_nb /))
             !don't want to re-read
             read_d_param = .false.
          end if
       elseif( read_d_param_mu_nb == d_param_no_spin_nb) then
          !number of parameters read in is only sufficient for no spin dependence in mu
          if ( mu_spin_dependent) then
             !not enough parameters read in
             call require(lhere, 'number of parameters in read_d_param_mu is equal to (order_mu_bf + 1) * nctype * spin_nb',.not. mu_spin_dependent)
          else
             !correct number of parameters so reshape into correct dimensions for d_param_mu_pseudo
             d_param_mu_pseudo = reshape( read_d_param_mu, (/ order_mu_bf + 1, nctype, spin_nb /), pad=read_d_param_mu)
             !note that only the first index has any meaning from the previous assignment so duplicate the meaningful values
!debug             d_param_mu_pseudo(:,:,2)=d_param_mu_pseudo(:,:,1)
             !don't want to re-read
             read_d_param = .false.
          end if
       else
          !wrong number of parameters read in
          call die(lhere,'read_d_param_mu has an incorrect number of parameters.')
       end if
    else
       !don't read in parameters
       !start from random guess
       call random_number(d_param_mu_pseudo(:,:,:))

       !are the parameters spin dependent?
       if (.not. mu_spin_dependent) then
          !no spin dependence
          d_param_mu_pseudo(:,:,2) = d_param_mu_pseudo(:,:,1)
       end if
    endif

    !no electron nuclear cusp conditions to impose because pseudo case

  end subroutine d_param_mu_pseudo_bld

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

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       !specify atom-type for each atom
       call object_needed('iwctype')
       ! number of electrons
       call object_needed('nelec')
       !threshhold for each nucleus type
       call object_needed('threshold_l')
       !distances of each electron of each walker to each nucleus
       call object_needed ('dist_en_wlk')

       return
    endif

    ! allocation
    call object_alloc('smooth_cutoff_g', smooth_cutoff_g, nelec, ncent, nwalk)

    smooth_cutoff_g(:,:,:) = 0.d0
    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
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

  !======================================================================================

  subroutine xi_ee_bld
    !---------------------------------------------------------------------------
    ! Description : build object xi_ee(:,:). This is a rank three array.
    !               It is the electron-electron backflow.
    !               It is necessary for the construction of the total backflow transformation.
    !               The first index has dimension equal to the number of spatial dimensions.
    !               The second index has dimensions equal to the number of electrons.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 22 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, walk_i, dim_i
    character(len=max_string_len_rout), save :: lhere = 'xi_ee_bld'

    ! header
    if (header_exe) then

       call object_create('xi_ee')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       ! number of electrons
       call object_needed('nelec')
       !cartesian vector between different electrons (for each walker)
       call object_needed('vec_ee_xyz_wlk')
       !dimension of space
       call object_needed('ndim')
       !is this an all-electron calculation?
       call object_needed('all_elec')

       return
    endif

    if (all_elec) then
       !all electron case
       call object_provide_in_node(lhere, 'eta_all_elec')

       !allocation

       call object_alloc('xi_ee', xi_ee, ndim, nelec, nwalk)

       !build xi_ee
       xi_ee(:,:,:) = 0.d0

       do walk_i = 1, nwalk
          do elec_i= 1, nelec
             do elec_j = 1, elec_i - 1
                do dim_i = 1, ndim
                   xi_ee(dim_i, elec_i, walk_i) = xi_ee(dim_i, elec_i, walk_i) + eta_all_elec(elec_i, elec_j, walk_i) * vec_ee_xyz_wlk(dim_i, elec_i, elec_j, walk_i)
                enddo
             enddo
             do elec_j = elec_i + 1, nelec
                do dim_i = 1, ndim
                   xi_ee(dim_i, elec_i, walk_i) = xi_ee(dim_i, elec_i, walk_i) + eta_all_elec(elec_i, elec_j, walk_i) * vec_ee_xyz_wlk(dim_i, elec_i, elec_j, walk_i)
                enddo
             enddo
          enddo
       enddo
    else
       !pseudo case
       call object_provide_in_node(lhere, 'eta_pseudo')
       !allocation

       call object_alloc('xi_ee', xi_ee, ndim, nelec, nwalk)

       !build xi_ee
       xi_ee(:,:,:) = 0.d0

       do walk_i = 1, nwalk
          do elec_i = 1, nelec
             do elec_j = 1, elec_i - 1
                do dim_i = 1, ndim
                   xi_ee(dim_i, elec_i, walk_i) = xi_ee(dim_i, elec_i, walk_i) + eta_pseudo(elec_i, elec_j, walk_i) * vec_ee_xyz_wlk(dim_i, elec_i, elec_j, walk_i)
                enddo
             enddo
             do elec_j = elec_i + 1, nelec
                do dim_i = 1, ndim
                   xi_ee(dim_i, elec_i, walk_i) = xi_ee(dim_i, elec_i, walk_i) + eta_pseudo(elec_i, elec_j, walk_i) * vec_ee_xyz_wlk(dim_i, elec_i, elec_j, walk_i)
                enddo
             enddo
          enddo
       enddo
    endif


  end subroutine xi_ee_bld



  !======================================================================================

  !======================================================================================

  subroutine eta_all_elec_bld
    !---------------------------------------------------------------------------
    ! Description : build object eta_all_elec(:,:,:). This is a rank three array.
    !               It is necessary for the construction of the electron-electron backflow terms in all electron case.
    !               The first index has dimension equal to the number of electrons.
    !               The second index has dimensions equal to the number of electrons.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 22 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, walk_i, order_i
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('eta_all_elec')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       ! number of electrons
       call object_needed('nelec')
       !number of up-spin electrons
       call object_needed('nup')
       !order of expansion
       call object_needed('order_eta_bf')
       !scaled distances of each electron of each walker to every other electron of that walker
       call object_needed ('scaled_dist_ee_wlk')
       !list of c parameters for eta backflow function for each spin relationship: up_up, down_down, up_down
       call object_needed('c_param_eta')
       !In all electron calculations, to impose the electron nuclear cusp, it is necessary to smoothly cutoff the backflow functions at the nucleus.
       call object_needed('smooth_cutoff_g')
       !Asymptotic part of the pade/power series expansion for each spin relationship: up_up, down_down, up_down
       call object_needed('asymp_eta')

       return
    endif

    ! allocation

    call object_alloc('eta_all_elec', eta_all_elec, nelec, nelec, nwalk)

    ! build eta_all_elec

    !initialize
    eta_all_elec (:,:,:) = 0.d0

    !calculate the up_up terms

    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over up-spin of electrons
       do elec_i = 1, nup
          !loop over up-spin electrons
          do elec_j = 1, elec_i - 1
             !no pade term because c_param_eta(1,up_up) = 0
             !power series terms of mu_all_elec
             do order_i=2, order_eta_bf
                eta_all_elec(elec_i, elec_j, walk_i) = eta_all_elec(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, up_up) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
          do elec_j = elec_i + 1, nup
             !no pade term because c_param_eta(1,up_up) = 0
             !power series terms of eta_all_elec
             do order_i=2, order_eta_bf
                eta_all_elec(elec_i, elec_j, walk_i) = eta_all_elec(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, up_up) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo

    !subtract off asymptotic value
    eta_all_elec(1:nup,1:nup,:) = eta_all_elec(1:nup,1:nup,:) - asymp_eta(up_up)

    !calculate the down_down terms

    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over down-spin electrons
       do elec_i = nup + 1, nelec
          !loop over down-spin electrons
          do elec_j = nup + 1, elec_i - 1
             !no pade term because c_param_eta(1,down_down) = 0
             !power series terms of eta_all_elec
             do order_i=2, order_eta_bf
                eta_all_elec(elec_i, elec_j, walk_i) = eta_all_elec(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, down_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
          do elec_j = elec_i + 1, nelec
             !no pade term because c_param_eta(1,down_down) = 0
             !power series terms of eta_all_elec
             do order_i=2, order_eta_bf
                eta_all_elec(elec_i, elec_j, walk_i) = eta_all_elec(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, down_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo

    !subtract off asymptotic value
    eta_all_elec(nup + 1:nelec, nup + 1:nelec, :) = eta_all_elec(nup + 1:nelec, nup + 1:nelec, :) - asymp_eta(down_down)

    !calculate the up_down terms

    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over up-spin electrons
       do elec_i = 1, nup
          !loop over down-spin electrons
          do elec_j = nup + 1, nelec
             !pade term
             eta_all_elec(elec_i, elec_j, walk_i) = c_param_eta(1, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i))
             !power series terms of eta_all_elec
             do order_i=2, order_eta_bf
                eta_all_elec(elec_i, elec_j, walk_i) = eta_all_elec(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo

    !subtract off asymptotic value
    eta_all_elec(1:nup, nup + 1:nelec, :) = eta_all_elec(1:nup, nup + 1:nelec, :) - asymp_eta(up_down)

    !calculate the down_up terms. Note that c_param_eta, and hence asymp_eta, are the same for up_down and down_up

    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over down-spin electrons
       do elec_i = nup + 1, nelec
          !loop over up-spin electrons
          do elec_j = 1, nup
             !pade term
             eta_all_elec(elec_i, elec_j, walk_i) = c_param_eta(1, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i))
             !power series terms of eta_all_elec
             do order_i=2, order_eta_bf
                eta_all_elec(elec_i, elec_j, walk_i) = eta_all_elec(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo

    !subtract off asymptotic value
    eta_all_elec(nup + 1:nelec, 1:nup, :) = eta_all_elec(nup + 1:nelec, 1:nup, :) - asymp_eta(up_down)

    !apply smooth cutoff function g
    do elec_i = 1, nelec
       do elec_j = 1, nelec
          do walk_i = 1, nwalk
             eta_all_elec(elec_i, elec_j, walk_i) = eta_all_elec(elec_i, elec_j, walk_i) * product(smooth_cutoff_g(elec_i, : , walk_i))
          enddo
       enddo
    enddo


  end subroutine eta_all_elec_bld

  !=============================================================================================


  !======================================================================================

  subroutine eta_pseudo_bld
    !---------------------------------------------------------------------------
    ! Description : build object eta_pseudo(:,:,:). This is a rank three array.
    !               It is necessary for the construction of the electron-electron backflow terms in the pseudo case.
    !               The first index has dimension equal to the number of electrons.
    !               The second index has dimensions equal to the number of electrons.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 22 Jul 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, walk_i, order_i
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('eta_pseudo')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       ! number of electrons
       call object_needed('nelec')
       !number of up-spin electrons
       call object_needed('nup')
       !order of expansion
       call object_needed('order_eta_bf')
       !scaled distances of each electron of each walker to every other electron of that walker
       call object_needed ('scaled_dist_ee_wlk')
       !list of c parameters for eta backflow function for each spin relationship: up_up, down_down, up_down
       call object_needed('c_param_eta')
       !Asymptotic part of the pade/power series expansion for each spin relationship: up_up, down_down, up_down
       call object_needed('asymp_eta')

       return
    endif

    ! allocation

    call object_alloc('eta_pseudo', eta_pseudo, nelec, nelec, nwalk)

    ! build eta_pseudo

    !initialize
    eta_pseudo (:,:,:) = 0.d0

    !calculate the up_up terms

    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over up-spin of electrons
       do elec_i = 1, nup
          !loop over up-spin electrons
          do elec_j = 1, elec_i - 1
             !no pade term because c_param_eta(1,up_up) = 0
             !power series terms of eta_pseudo
             do order_i=2, order_eta_bf
                eta_pseudo(elec_i, elec_j, walk_i) = eta_pseudo(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, up_up) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
          do elec_j = elec_i + 1, nup
             !no pade term because c_param_eta(1,up_up) = 0
             !power series terms of eta_pseudo
             do order_i=2, order_eta_bf
                eta_pseudo(elec_i, elec_j, walk_i) = eta_pseudo(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, up_up) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo

    !subtract off asymptotic value
    eta_pseudo(1:nup,1:nup,:) = eta_pseudo(1:nup,1:nup,:) - asymp_eta(up_up)

    !calculate the down_down terms

    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over down-spin electrons
       do elec_i = nup + 1, nelec
          !loop over down-spin electrons
          do elec_j = nup + 1, elec_i - 1
             !no pade term because c_param_eta(1,down_down) = 0
             !power series terms of eta_pseudo
             do order_i=2, order_eta_bf
                eta_pseudo(elec_i, elec_j, walk_i) = eta_pseudo(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, down_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
          do elec_j = elec_i + 1, nelec
             !no pade term because c_param_eta(1,down_down) = 0
             !power series terms of eta_pseudo
             do order_i=2, order_eta_bf
                eta_pseudo(elec_i, elec_j, walk_i) = eta_pseudo(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, down_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo

    !subtract off asymptotic value
    eta_pseudo(nup + 1:nelec, nup + 1:nelec, :) = eta_pseudo(nup + 1:nelec, nup + 1:nelec, :) - asymp_eta(down_down)

    !calculate the up_down terms

    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over up-spin electrons
       do elec_i = 1, nup
          !loop over down-spin electrons
          do elec_j = nup + 1, nelec
             !pade term
             eta_pseudo(elec_i, elec_j, walk_i) = c_param_eta(1, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i))
             !power series terms of eta_pseudo
             do order_i=2, order_eta_bf
                eta_pseudo(elec_i, elec_j, walk_i) = eta_pseudo(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo

    !subtract off asymptotic value
    eta_pseudo(1:nup, nup + 1:nelec, :) = eta_pseudo(1:nup, nup + 1:nelec, :) - asymp_eta(up_down)

    !calculate the down_down terms. Note that c_param_eta, and hence asymp_eta, are the same for up_down and down_up

    !loop over number of walkers
    do walk_i = 1, nwalk
       !loop over down-spin electrons
       do elec_i = nup + 1, nelec
          !loop over up-spin electrons
          do elec_j = 1, nup
             !pade term
             eta_pseudo(elec_i, elec_j, walk_i) = c_param_eta(1, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i) / (1.d0 + c_param_eta(2, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i))
             !power series terms of eta_pseudo
             do order_i=2, order_eta_bf
                eta_pseudo(elec_i, elec_j, walk_i) = eta_pseudo(elec_i, elec_j, walk_i) + c_param_eta(order_i + 1, up_down) * scaled_dist_ee_wlk(elec_i, elec_j, walk_i)**order_i
             enddo
          enddo
       enddo
    enddo

    !subtract off asymptotic value
    eta_pseudo(nup + 1:nelec, 1:nup, :) = eta_pseudo(nup + 1:nelec, 1:nup, :) - asymp_eta(up_down)

    !no cutoff because pseudo

  end subroutine eta_pseudo_bld

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
    integer  :: order_i
    integer, parameter :: spin_dependencies_nb = 3    !because three choices up_up, down_down, up_down
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('asymp_eta')

       !order of expansion
       call object_needed('order_eta_bf')
       !list of c parameters for eta backflow function for each spin dependence
       call object_needed('c_param_eta')
       !asymptotic scaled distance for electron electron scaling
       call object_needed('asymp_scaled_dist_two_body')

       return
    endif

    ! allocation
    call object_alloc('asymp_eta', asymp_eta, spin_dependencies_nb)

    ! build asymp_eta
    !initialize
    asymp_eta(:) = 0.d0

    !up_up
    !no pade term because c_param_eta(1, up_up) = 0 for electron-electron cusp condition to be satisfied
    !power series terms
    do order_i=2, order_eta_bf
       asymp_eta(up_up) = asymp_eta(up_up) + c_param_eta(order_i + 1, up_up) * asymp_scaled_dist_two_body**order_i
    enddo

    !down_down
    !no pade term because c_param_eta(1, down_down) = 0 for electron-electron cusp condition to be satisfied
    !power series terms
    do order_i=2, order_eta_bf
       asymp_eta(down_down) = asymp_eta(down_down) + c_param_eta(order_i + 1, down_down) * asymp_scaled_dist_two_body**order_i
    enddo

    !up_down (or down_up)
    !pade term
    asymp_eta(up_down) = c_param_eta(1, up_down) * asymp_scaled_dist_two_body / (1.d0 + c_param_eta(2, up_down) * asymp_scaled_dist_two_body)
    !power series terms
    do order_i=2, order_eta_bf
       asymp_eta(up_down) = asymp_eta(up_down) + c_param_eta(order_i + 1, up_down) * asymp_scaled_dist_two_body**order_i
    enddo

  end subroutine asymp_eta_bld

  !==============================================================================================


  !==============================================================================================

  subroutine scaled_dist_ee_wlk_bld

    !---------------------------------------------------------------------------
    ! Description : build scaled_dist_ee_wlk(:,:,:). This is a rank three array.
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

       ! number of electrons
       call object_needed('nelec')
       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !distances of each electron of each walker to every other electron of that walker
       call object_needed ('dist_ee_wlk')

       !Make sure scalek_bf, isc_bf are up to date even though the values are not
       !needed directly in this build routine. They are, however, needed in a subroutine that is called
       !by this build routine.
       !scale factor for the scaling functions
       call object_needed('scalek_bf')
       !scaling form for the scaling functions
       call object_needed('isc_bf')

       return
    endif

    ! allocation

    !this can be made more efficient
    call object_alloc('scaled_dist_ee_wlk', scaled_dist_ee_wlk, nelec, nelec, nwalk)

    ! build scaled_dist_ee_wlk
    !initialize
    scaled_dist_ee_wlk(:,:,:) = 0.d0

    !calculate scaled distances
    call scaling_func_two_body(dist_ee_wlk,scaled_dist_ee_wlk)

  end subroutine scaled_dist_ee_wlk_bld

  !==============================================================================================

  !==============================================================================================

  subroutine c_param_eta_bld

    !---------------------------------------------------------------------------
    ! Description : build c_param_eta. This is a rank 2 array.
    !               List of c parameters for eta backflow function for each spin dependence
    !               The first index is which parameter.
    !               The second index is the type of spin dependence: 3(default), 2, or 1  (as explained in help section of backflow menu)
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
    integer  :: c_param_spin_3_nb, c_param_spin_2_nb, c_param_spin_1_nb
    integer, parameter :: spin_dependencies_nb = 3    !because three choices up_up, down_down, up_down
    character(len=max_string_len_rout), save :: lhere = 'c_param_eta_bld'

    ! header
    if (header_exe) then

       call object_create('c_param_eta')

       !order of expansion
       call object_needed('order_eta_bf')
       !asymptotic scaled distance for electron electron scaling
       call object_needed('asymp_scaled_dist_two_body')
       !spin dependence of the parameters
       call object_needed('eta_spin_dependence')
       !have the parameters been read in?
       call object_needed('read_c_param')
       !how many parameters were read in?
       call object_needed('read_c_param_eta_nb')

       return
    endif



    !Note that we need to have the option of using read_c_param_eta if it is read in
    !However, we don't know apriori if it will be read in. 
    !Thus we don't call object_needed('read_c_param_eta')
    !To accomplish this correctly we need to use the utility object_provide_in_node
    !We cant simply use object_provide because we are inside the tree.

    if (read_c_param) then
       !d parameters were read in
       call object_provide_in_node(lhere, 'read_c_param_eta')
    end if


    ! allocation
    call object_alloc('c_param_eta', c_param_eta, order_eta_bf + 1, spin_dependencies_nb)

    !do we read in parameters?
    if (read_c_param) then
       !read in parameters
       c_param_spin_3_nb = (order_eta_bf + 1)  * 3         !  up_up != down_down != up_down
       c_param_spin_2_nb = (order_eta_bf + 1)  * 2         !  up_up = down_down != up_down
       c_param_spin_1_nb = order_eta_bf + 1                !  up_up = down_down = up_down
       if (read_c_param_eta_nb == c_param_spin_3_nb) then
          !enough parameters were read in for:  up_up != down_down != up_down
          if (eta_spin_dependence .ne. 3) then
             !too many parameters were read in
             call require(lhere, 'number of parameters correct for eta_spin_dependence specified', eta_spin_dependence .eq. 3)
          else
             !correct number of parameters so reshape into correct dimensions for c_param_eta
             c_param_eta = reshape( read_c_param_eta, (/ order_eta_bf + 1, spin_dependencies_nb /)) !  up_up != down_down != up_down
             !don't want to re-read
             read_c_param = .false.
          end if
       elseif( read_c_param_eta_nb == c_param_spin_2_nb) then
          !enough parameters read in for:  up_up = down_down != up_down
          if ( eta_spin_dependence .ne. 2) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for eta_spin_dependence specified', eta_spin_dependence .eq. 2)
          else
             !correct number of parameters so reshape into correct dimensions for c_param_eta
             c_param_eta = reshape( read_c_param_eta, (/ order_eta_bf + 1, spin_dependencies_nb /), pad=read_c_param_eta )   !up_up = down_down != up_down
             ! up_up is in c_param_eta(:,1), up_down is in c_param_eta(:,2) and c_param_eta(:,3) has garbage in it
             !want up_up in c_param_eta(:,1), up_up=down_down in c_param_eta(:,2), up_down in c_param_eta(:,3)
             c_param_eta(:,3)=c_param_eta(:,2)
             c_param_eta(:,2)=c_param_eta(:,1)
             !don't want to re-read
             read_c_param = .false.
          end if
       elseif( read_c_param_eta_nb == c_param_spin_1_nb) then
          !enough parameters read in for:  up_up = down_down = up_down
          if ( eta_spin_dependence .ne. 1) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for eta_spin_dependence specified', eta_spin_dependence .eq. 1)
          else
             !correct number of parameters so reshape into correct dimensions for c_param_eta
             c_param_eta = reshape( read_c_param_eta, (/ order_eta_bf + 1, spin_dependencies_nb /), pad=read_c_param_eta)   !up_up = down_down = up_down
             ! up_up is in c_param_eta(:,1), c_param_eta(:,2) has garbage in it and c_param_eta(:,3) has garbage in it
             !want up_up in c_param_eta(:,1), up_up=down_down in c_param_eta(:,2), up_up=up_down in c_param_eta(:,3)
             c_param_eta(:,2)=c_param_eta(:,1)
             c_param_eta(:,3)=c_param_eta(:,1)
             !don't want to re-read
             read_c_param = .false.
          end if
       else
          !wrong number of parameters
          call die(lhere,'read_c_param_eta has an incorrect number of parameters.')
       end if
    else
       !don't read in parameters
       !start from random guess
       call random_number(c_param_eta(:,:))

       !What is the spin dependence?
       if ( eta_spin_dependence .eq. 2) then
          !        up_up = down_down != up_down
          c_param_eta(:,2) = c_param_eta(:,1)
       elseif( eta_spin_dependence .eq. 1) then
          !        up_up = down_down = up_down
          c_param_eta(:,2) = c_param_eta(:,1)
          c_param_eta(:,3) = c_param_eta(:,1)
       end if
    endif

    !impose electron electron cusp conditions
    !Note that we assume the two-body scaling functions evaluated at 0 are 0 and have unit slope at 0.
    c_param_eta(1, 1) = 0.d0 ! parallel spin
    c_param_eta(1, 2) = 0.d0 ! parallel spin
    if (eta_spin_dependence .eq. 1) then
       c_param_eta(1, 3) = 0.d0
    endif

  end subroutine c_param_eta_bld

  !==============================================================================================


  !======================================================================================

  subroutine xi_een_phi_bld
    !---------------------------------------------------------------------------
    ! Description : build object xi_een_phi(:,:). This is a rank three array.
    !               It is part of electron-electron-nuclear backflow.
    !               It is necessary for the construction of the total backflow transformation.
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
    integer  :: elec_i, elec_j, cent_i, walk_i, dim_i
    character(len=max_string_len_rout), save :: lhere = 'xi_een_phi_bld'

    ! header
    if (header_exe) then

       call object_create('xi_een_phi')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       ! number of electrons
       call object_needed('nelec')
       !cartesian vector between electrons and other electrons (for each walker)
       call object_needed('vec_ee_xyz_wlk')
       !dimension of space
       call object_needed('ndim')
       !is this an all-electron calculation?
       call object_needed('all_elec')

       return
    endif

    if (all_elec) then
       !all electron case
       call object_provide_in_node(lhere, 'phi_all_elec')

       !allocation

       call object_alloc('xi_een_phi', xi_een_phi, ndim, nelec, nwalk)

       !build xi_een_phi
       xi_een_phi(:,:,:) = 0.d0

       !loop over walkers
       do walk_i = 1, nwalk
          !loop over nuclei
          do cent_i = 1, ncent
             !loop over electrons
             do elec_i = 1, nelec
                !loop over other electrons before elec_i
                do elec_j = 1, elec_i -1
                   !loop over spatial dimensions
                   do dim_i = 1, ndim
                      xi_een_phi(dim_i, elec_i, walk_i) = xi_een_phi(dim_i, elec_i, walk_i) +  phi_all_elec(elec_i, elec_j, cent_i, walk_i) * vec_ee_xyz_wlk(dim_i, elec_i, elec_j, walk_i)
                   enddo
                enddo
                !loop over other electrons after elec_i
                do elec_j = elec_i + 1, nelec
                   !loop over spatial dimensions
                   do dim_i = 1, ndim
                      xi_een_phi(dim_i, elec_i, walk_i) = xi_een_phi(dim_i, elec_i, walk_i) +  phi_all_elec(elec_i, elec_j, cent_i, walk_i) * vec_ee_xyz_wlk(dim_i, elec_i, elec_j, walk_i)
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       !pseudo case
       call object_provide_in_node(lhere, 'phi_pseudo')

       !allocation

       call object_alloc('xi_een_phi', xi_een_phi, ndim, nelec, nwalk)

       !build xi_en
       xi_een_phi(:,:,:) = 0.d0

       !loop over walkers
       do walk_i = 1, nwalk
          !loop over nuclei
          do cent_i = 1, ncent
             !loop over electrons
             do elec_i = 1, nelec
                !loop over other electrons before elec_i
                do elec_j = 1, elec_i -1
                   !loop over spatial dimensions
                   do dim_i = 1, ndim
                      xi_een_phi(dim_i, elec_i, walk_i) = xi_een_phi(dim_i, elec_i, walk_i) +  phi_pseudo(elec_i, elec_j, cent_i, walk_i) * vec_ee_xyz_wlk(dim_i, elec_i, elec_j, walk_i)
                   enddo
                enddo
                !loop over other electrons after elec_i
                do elec_j = elec_i + 1, nelec
                   !loop over spatial dimensions
                   do dim_i = 1, ndim
                      xi_een_phi(dim_i, elec_i, walk_i) = xi_een_phi(dim_i, elec_i, walk_i) +  phi_pseudo(elec_i, elec_j, cent_i, walk_i) * vec_ee_xyz_wlk(dim_i, elec_i, elec_j, walk_i)
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif


  end subroutine xi_een_phi_bld

  !======================================================================================

  !======================================================================================

  subroutine phi_all_elec_bld
    !---------------------------------------------------------------------------
    ! Description : build object phi_all_elec(:,:,:,:). This is a rank four array.
    !               It is necessary for the construction of the electron-electron-nuclear backflow terms in all electron case.
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
    integer  :: elec_i, elec_j, cent_i, walk_i, order_i, order_j, order_k, prod_i, param_i
    integer  :: order_k_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('phi_all_elec')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       !specify atom-type for each atom
       call object_needed('iwctype')
       ! number of electrons
       call object_needed('nelec')
       !number of up-spin electrons
       call object_needed('nup')
       !order of expansion
       call object_needed('order_phi_bf')
       !scaled distances of each electron of each walker to each nucleus (three body scaling functions)
       call object_needed ('scaled_dist_een_en_wlk')
       !scaled distances of each electron of each walker to every other electron of that walker (three body scaling functions)
       call object_needed ('scaled_dist_een_ee_wlk')
       !list of a parameters (satisfying all-elec conditions) for phi backflow function for each center type, and each spin dependence
       call object_needed('a_param_phi_all_elec')
       !In all electron calculations, to impose the electron nuclear cusp, it is necessary to smoothly cutoff the backflow functions at the nucleus.
       call object_needed('smooth_cutoff_g')


       return
    endif

    ! allocation

    call object_alloc('phi_all_elec', phi_all_elec, nelec, nelec, ncent, nwalk)

    ! build phi_all_elec

    !initialize
    phi_all_elec (:,:,:,:) = 0.d0

    !up_up terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !loop over other up-spin electrons below elec_i
             do elec_j=1, elec_i - 1
                !power series terms of phi_all_elec as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_all_elec(elec_i, elec_j, cent_i, walk_i) = phi_all_elec(elec_i, elec_j, cent_i, walk_i) + a_param_phi_all_elec(param_i, iwctype(cent_i), up_up) &
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k &
                                 & +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) &
                                 & *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !loop over other up-spin electrons above elec_i
             do elec_j= elec_i + 1, nup
                !power series terms of phi_all_elec
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_all_elec(elec_i, elec_j, cent_i, walk_i) = phi_all_elec(elec_i, elec_j, cent_i, walk_i) + a_param_phi_all_elec(param_i, iwctype(cent_i), up_up) &
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k&
                                 & +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)&
                                 & *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo



    !down_down terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons
          do elec_i = nup + 1, nelec
             !loop over other down-spin electrons below elec_i
             do elec_j=nup + 1, elec_i - 1
                !power series terms of phi_all_elec as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_all_elec(elec_i, elec_j, cent_i, walk_i) = phi_all_elec(elec_i, elec_j, cent_i, walk_i) + a_param_phi_all_elec(param_i, iwctype(cent_i), down_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !loop over other down-spin electrons above elec_i
             do elec_j= elec_i + 1, nelec
                !power series terms of phi_all_elec
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_all_elec(elec_i, elec_j, cent_i, walk_i) = phi_all_elec(elec_i, elec_j, cent_i, walk_i) + a_param_phi_all_elec(param_i, iwctype(cent_i), down_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k&
                                 & +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)&
                                 & *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


    !up_down terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !loop over down-spin electrons
             do elec_j= nup + 1, nelec
                !power series terms of phi_all_elec as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_all_elec(elec_i, elec_j, cent_i, walk_i) = phi_all_elec(elec_i, elec_j, cent_i, walk_i) + a_param_phi_all_elec(param_i, iwctype(cent_i), up_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


    !down_up terms (note that a_param_phi_all_elec are same for up_down and down_up)
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons
          do elec_i = nup + 1, nelec
             !loop over up-spin electrons
             do elec_j = 1, nup
                !power series terms of phi_all_elec as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_all_elec(elec_i, elec_j, cent_i, walk_i) = phi_all_elec(elec_i, elec_j, cent_i, walk_i) + a_param_phi_all_elec(param_i, iwctype(cent_i), up_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    !apply smooth cutoff function g
    do cent_i=1, ncent
       do elec_j=1, nelec
          do prod_i=1, cent_i-1
             phi_all_elec(:,elec_j,cent_i,:) = phi_all_elec(:,elec_j,cent_i,:) * smooth_cutoff_g(:,prod_i,:)
          end do
          do prod_i=cent_i+1, ncent
             phi_all_elec(:,elec_j, cent_i,:) = phi_all_elec(:,elec_j, cent_i,:) * smooth_cutoff_g(:,prod_i,:)
          end do
       end do
    enddo

  end subroutine phi_all_elec_bld

  !=============================================================================================

  !======================================================================================

  subroutine phi_pseudo_bld
    !---------------------------------------------------------------------------
    ! Description : build object phi_pseudo(:,:,:,:). This is a rank four array.
    !               It is necessary for the construction of the electron-electron-nuclear backflow terms in the pseudo case.
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
    integer  :: elec_i, elec_j, cent_i, walk_i, order_i, order_j, order_k, param_i
    integer  :: order_k_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('phi_pseudo')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       !specify atom-type for each atom
       call object_needed('iwctype')
       ! number of electrons
       call object_needed('nelec')
       !number of up-spin electrons
       call object_needed('nup')
       !order of expansion
       call object_needed('order_phi_bf')
       !scaled distances of each electron of each walker to each nucleus (three body scaling functions)
       call object_needed ('scaled_dist_een_en_wlk')
       !scaled distances of each electron of each walker to every other electron of that walker (three body scaling functions)
       call object_needed ('scaled_dist_een_ee_wlk')
       !list of a parameters (satisfying all-elec conditions) for phi backflow function for each center type, and each spin dependence
       call object_needed('a_param_phi_pseudo')

       return
    endif

    ! allocation

    call object_alloc('phi_pseudo', phi_pseudo, nelec, nelec, ncent, nwalk)

    ! build phi_pseudo

    !initialize
    phi_pseudo (:,:,:,:) = 0.d0

    !up_up terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !loop over other up-spin electrons below elec_i
             do elec_j=1, elec_i - 1
                !power series terms of phi_pseudo as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_pseudo(elec_i, elec_j, cent_i, walk_i) = phi_pseudo(elec_i, elec_j, cent_i, walk_i) + a_param_phi_pseudo(param_i, iwctype(cent_i), up_up)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !loop over other up-spin electrons above elec_i
             do elec_j= elec_i + 1, nup
                !power series terms of phi_pseudo
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_pseudo(elec_i, elec_j, cent_i, walk_i) = phi_pseudo(elec_i, elec_j, cent_i, walk_i) + a_param_phi_pseudo(param_i, iwctype(cent_i), up_up)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo



    !down_down terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons
          do elec_i = nup + 1, nelec
             !loop over other down-spin electrons below elec_i
             do elec_j=nup + 1, elec_i - 1
                !power series terms of phi_pseudo as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_pseudo(elec_i, elec_j, cent_i, walk_i) = phi_pseudo(elec_i, elec_j, cent_i, walk_i) + a_param_phi_pseudo(param_i, iwctype(cent_i), down_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !loop over other down-spin electrons above elec_i
             do elec_j= elec_i + 1, nelec
                !power series terms of phi_pseudo
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_pseudo(elec_i, elec_j, cent_i, walk_i) = phi_pseudo(elec_i, elec_j, cent_i, walk_i) + a_param_phi_pseudo(param_i, iwctype(cent_i), down_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


    !up_down terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !loop over down-spin electrons
             do elec_j= nup + 1, nelec
                !power series terms of phi_pseudo as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_pseudo(elec_i, elec_j, cent_i, walk_i) = phi_pseudo(elec_i, elec_j, cent_i, walk_i) + a_param_phi_pseudo(param_i, iwctype(cent_i), up_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


    !down_up terms (note that a_param_phi_pseudo are same for up_down and down_up)
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons
          do elec_i = nup + 1, nelec
             !loop over up-spin electrons
             do elec_j = 1, nup
                !power series terms of phi_pseudo as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_phi_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            phi_pseudo(elec_i, elec_j, cent_i, walk_i) = phi_pseudo(elec_i, elec_j, cent_i, walk_i) + a_param_phi_pseudo(param_i, iwctype(cent_i), up_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    !no cutoff because pseudo

  end subroutine phi_pseudo_bld

  !=============================================================================================



  !==============================================================================================

  subroutine scaled_dist_een_ee_wlk_bld

    !---------------------------------------------------------------------------
    ! Description : build scaled_dist_een_ee_wlk(:,:,:). This is a rank three array.
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

       ! number of electrons
       call object_needed('nelec')
       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !distances of each electron of each walker to every other electron of that walker
       call object_needed ('dist_ee_wlk')

       !Make sure scalek_bf, isc_bf are up to date even though the values are not
       !needed directly in this build routine. They are, however, needed in a subroutine that is called
       !by this build routine.
       !scale factor for the scaling functions
       call object_needed('scalek_bf')
       !scaling form for the scaling functions
       call object_needed('isc_bf')

       return
    endif

    ! allocation

    call object_alloc('scaled_dist_een_ee_wlk', scaled_dist_een_ee_wlk, nelec, nelec, nwalk)

    ! build scaled_dist_een_ee_wlk
    !initialize
    scaled_dist_een_ee_wlk(:,:,:) = 0.d0

    !calculate scaled distances
    call scaling_func_three_body(dist_ee_wlk,scaled_dist_een_ee_wlk)

  end subroutine scaled_dist_een_ee_wlk_bld

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

       ! number of electrons
       call object_needed('nelec')
       !number of nucei
       call object_needed('ncent')
       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !distances of each electron of each walker to each nucleus
       call object_needed ('dist_en_wlk')

       !Make sure scalek_bf, isc_bf are up to date even though the values are not
       !needed directly in this build routine. They are, however, needed in a subroutine that is called
       !by this build routine.
       !scale factor for the scaling functions
       call object_needed('scalek_bf')
       !scaling form for the scaling functions
       call object_needed('isc_bf')

       return
    endif

    ! allocation

    call object_alloc('scaled_dist_een_en_wlk', scaled_dist_een_en_wlk, nelec, ncent, nwalk)

    ! build scaled_dist_een_en_wlk
    !initialize
    scaled_dist_een_en_wlk(:,:,:) = 0.d0

    !calculate scaled distances
    call scaling_func_three_body(dist_en_wlk,scaled_dist_een_en_wlk)

  end subroutine scaled_dist_een_en_wlk_bld

  !==============================================================================================



  !==============================================================================================

  subroutine a_param_phi_pseudo_bld

    !---------------------------------------------------------------------------
    ! Description : build a_param_phi_pseudo. This is a rank 3 array.
    !               List of a parameters for phi backflow function for each center type, and each spin dependence
    !               The first index is which parameter.
    !               The second index is the type of center.
    !               The third index is the type of spin dependence.
    !
    ! Created     : F. Petruzielo, 13 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: a_param_spin_3_nb, a_param_spin_2_nb, a_param_spin_1_nb, dep_a_param_i, ctype_i, spin_dependency_i
    character(len=max_string_len_rout), save :: lhere = 'a_param_phi_pseudo_bld'
    integer, parameter :: spin_dependencies_nb = 3    !because three choices up_up, down_down, up_down

    ! header
    if (header_exe) then

       call object_create('a_param_phi_pseudo')

       !order of expansion
       call object_needed('order_phi_bf')
       !how many types of atoms
       call object_needed('nctype')
       !what is the spin dependence?
       call object_needed('phi_spin_dependence')
       !have the parameters been read in?
       call object_needed('read_a_param')
       !how many parameters were read in?
       call object_needed('read_a_param_phi_nb')
       !what are the dependent parameters for parallel spin?
       call object_needed('dep_a_param_phi_pseudo_parallel')
       !what are the dependent parameters for anti parallel spin?
       call object_needed('dep_a_param_phi_pseudo_anti_parallel')
       !what are the conditions for parallel spin?
       call object_needed('a_param_phi_pseudo_parallel_cond')
       !what are the conditions for anti parallel spin?
       call object_needed('a_param_phi_pseudo_anti_parallel_cond')
       !How many a parameters
       call object_needed('nb_a_param_phi')

       return
    endif


    !Note that we need to have the option of using read_a_param_phi if it is read in
    !However, we don't know apriori if it will be read in.
    !Thus we don't call object_needed('read_a_param_phi')
    !To accomplish this correctly we need to use the utility object_provide_in_node
    !We cant simply use object_provide because we are inside the tree.

    if (read_a_param) then
       !a parameters were read in
       call object_provide_in_node(lhere, 'read_a_param_phi')
    end if

    ! allocation
    call object_alloc('a_param_phi_pseudo', a_param_phi_pseudo, nb_a_param_phi, nctype, spin_dependencies_nb)

    !do we read in parameters?
    if (read_a_param) then
       !read in parameters
       a_param_spin_3_nb = (nb_a_param_phi) * nctype * 3         !  up_up != down_down != up_down
       a_param_spin_2_nb = (nb_a_param_phi) * nctype * 2         !  up_up = down_down != up_down
       a_param_spin_1_nb = (nb_a_param_phi) * nctype             !  up_up = down_down = up_down
       if (read_a_param_phi_nb == a_param_spin_3_nb) then
          !enough parameters were read in for:  up_up != down_down != up_down
          if (phi_spin_dependence .ne. 3) then
             !too many parameters were read in
             call require(lhere, 'number of parameters correct for phi_spin_dependence specified', phi_spin_dependence .eq. 3)
          else
             !correct number of parameters so reshape into correct dimensions for a_param_phi
             a_param_phi_pseudo = reshape( read_a_param_phi, (/ nb_a_param_phi, nctype, spin_dependencies_nb /)) !  up_up != down_down != up_down
             !don't want to re-read
             read_a_param = .false.
          end if
       elseif( read_a_param_phi_nb == a_param_spin_2_nb) then
          !enough parameters read in for:  up_up = down_down != up_down
          if ( phi_spin_dependence .ne. 2) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for phi_spin_dependence specified', phi_spin_dependence .eq. 2)
          else
             !correct number of parameters so reshape into correct dimensions for a_param_phi_pseudo
             a_param_phi_pseudo = reshape( read_a_param_phi, (/ nb_a_param_phi, nctype, spin_dependencies_nb /), pad=read_a_param_phi)   !up_up = down_down != up_down
             ! up_up is in a_param_phi_pseudo(:,:,1), up_down is in a_param_phi_pseudo(:,:,2) and a_param_phi_pseudo(:,:,3) has garbage in it
             !want up_up in a_param_phi_pseudo(:,:,1), up_up=down_down in a_param_phi_pseudo(:,:,2), up_down in a_param_phi_pseudo(:,:,3)
             a_param_phi_pseudo(:,:,3)=a_param_phi_pseudo(:,:,2)
             a_param_phi_pseudo(:,:,2)=a_param_phi_pseudo(:,:,1)
             !don't want to re-read
             read_a_param = .false.
          end if
       elseif( read_a_param_phi_nb == a_param_spin_1_nb) then
          !enough parameters read in for:  up_up = down_down = up_down
          if ( phi_spin_dependence .ne. 1) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for phi_spin_dependence specified', phi_spin_dependence .eq. 1)
          else
             !correct number of parameters so reshape into correct dimensions for a_param_phi_pseudo
             a_param_phi_pseudo = reshape( read_a_param_phi, (/ nb_a_param_phi, nctype, spin_dependencies_nb /), pad=read_a_param_phi)   !up_up = down_down = up_down
             ! up_up is in a_param_phi_pseudo(:,:,1), a_param_phi_pseudo(:,:,2) has garbage in it and a_param_phi_pseudo(:,:,3) has garbage in it
             !want up_up in a_param_phi_pseudo(:,:,1), up_up=down_down in a_param_phi_pseudo(:,:,2), up_up=up_down in a_param_phi_pseudo(:,:,3)
             a_param_phi_pseudo(:,:,2)=a_param_phi_pseudo(:,:,1)
             a_param_phi_pseudo(:,:,3)=a_param_phi_pseudo(:,:,1)
             !don't want to re-read
             read_a_param = .false.
          end if
       else
          !wrong number of parameters
          call die(lhere,'read_a_param_phi has an incorrect number of parameters.')
       end if
    else
       !don't read in parameters
       !start from random guess
       call random_number(a_param_phi_pseudo(:,:,:))
    endif



    !impose cusp
    !To impose the cusp we take the pivot positions (dependent parameters) in our reduced row echelon form of the matrix of cusp constraints
    !and we set this dependent parameter equal to the negative of the dot product of the rest of the row with
    !the corresponding a parameters. This is nothing fancy, but simply a standard way of solving linear equations

    !What is the spin dependence?
    if ( phi_spin_dependence .eq. 3) then
       !impose cusp parallel
       !loop over depend parameters
       do dep_a_param_i=1, size(dep_a_param_phi_pseudo_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             !loop of parallel spins: up-up and down-down
             do spin_dependency_i=1, 2
                a_param_phi_pseudo(dep_a_param_phi_pseudo_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                     & dot_product( -a_param_phi_pseudo_parallel_cond(dep_a_param_i, dep_a_param_phi_pseudo_parallel(dep_a_param_i)+1:),&
                     & a_param_phi_pseudo(dep_a_param_phi_pseudo_parallel(dep_a_param_i) + 1:, ctype_i, spin_dependency_i))
             end do
          end do
       end do
       !impose cusp anti_parallel
       spin_dependency_i = 3 !up_down or down_up
       !loop over depend parameters
       do dep_a_param_i=1, size(dep_a_param_phi_pseudo_anti_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             a_param_phi_pseudo(dep_a_param_phi_pseudo_anti_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -a_param_phi_pseudo_anti_parallel_cond(dep_a_param_i, dep_a_param_phi_pseudo_anti_parallel(dep_a_param_i)+1:),&
                  & a_param_phi_pseudo(dep_a_param_phi_pseudo_anti_parallel(dep_a_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
    elseif ( phi_spin_dependence .eq. 2) then
       !        up_up = down_down != up_down
       !impose cusp parallel
       spin_dependency_i = 1 !up-up
       !loop over depend parameters
       do dep_a_param_i=1, size(dep_a_param_phi_pseudo_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             a_param_phi_pseudo(dep_a_param_phi_pseudo_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -a_param_phi_pseudo_parallel_cond(dep_a_param_i, dep_a_param_phi_pseudo_parallel(dep_a_param_i)+1:),&
                  & a_param_phi_pseudo(dep_a_param_phi_pseudo_parallel(dep_a_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !impose cusp anti_parallel
       spin_dependency_i = 3 !up_down or down_up
       !loop over depend parameters
       do dep_a_param_i=1, size(dep_a_param_phi_pseudo_anti_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             a_param_phi_pseudo(dep_a_param_phi_pseudo_anti_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -a_param_phi_pseudo_anti_parallel_cond(dep_a_param_i, dep_a_param_phi_pseudo_anti_parallel(dep_a_param_i)+1:),&
                  & a_param_phi_pseudo(dep_a_param_phi_pseudo_anti_parallel(dep_a_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !require up-up = down-down
       a_param_phi_pseudo(:,:,2) = a_param_phi_pseudo(:,:,1)
    elseif( phi_spin_dependence .eq. 1) then
       !        up_up = down_down = up_down
       !impose cusp parallel
       spin_dependency_i = 1 !up-up
       !loop over depend parameters
       do dep_a_param_i=1, size(dep_a_param_phi_pseudo_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             a_param_phi_pseudo(dep_a_param_phi_pseudo_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -a_param_phi_pseudo_parallel_cond(dep_a_param_i, dep_a_param_phi_pseudo_parallel(dep_a_param_i)+1:),&
                  & a_param_phi_pseudo(dep_a_param_phi_pseudo_parallel(dep_a_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !require up-up = down-down=up-down
       a_param_phi_pseudo(:,:,2) = a_param_phi_pseudo(:,:,1)
       a_param_phi_pseudo(:,:,3) = a_param_phi_pseudo(:,:,1)
    end if




  end subroutine a_param_phi_pseudo_bld

  !==============================================================================================





  !==============================================================================================

  subroutine a_param_phi_all_elec_bld

    !---------------------------------------------------------------------------
    ! Description : build a_param_phi_all_elec. This is a rank 3 array.
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
    integer  :: a_param_spin_3_nb, a_param_spin_2_nb, a_param_spin_1_nb, dep_a_param_i, ctype_i, spin_dependency_i
    character(len=max_string_len_rout), save :: lhere = 'a_param_phi_all_elec_bld'
    integer, parameter :: spin_dependencies_nb = 3    !because three choices up_up, down_down, up_down

    ! header
    if (header_exe) then

       call object_create('a_param_phi_all_elec')

       !order of expansion
       call object_needed('order_phi_bf')
       !how many types of atoms
       call object_needed('nctype')
       !what is the spin dependence?
       call object_needed('phi_spin_dependence')
       !have the parameters been read in?
       call object_needed('read_a_param')
       !how many parameters were read in?
       call object_needed('read_a_param_phi_nb')
       !what are the dependent parameters for parallel spin?
       call object_needed('dep_a_param_phi_all_elec_parallel')
       !what are the dependent parameters for anti parallel spin?
       call object_needed('dep_a_param_phi_all_elec_anti_parallel')
       !what are the conditions for parallel spin?
       call object_needed('a_param_phi_all_elec_parallel_cond')
       !what are the conditions for anti parallel spin?
       call object_needed('a_param_phi_all_elec_anti_parallel_cond')
       !How many a parameters
       call object_needed('nb_a_param_phi')

       return
    endif


    !Note that we need to have the option of using read_a_param_phi if it is read in
    !However, we don't know apriori if it will be read in.
    !Thus we don't call object_needed('read_a_param_phi')
    !To accomplish this correctly we need to use the utility object_provide_in_node
    !We cant simply use object_provide because we are inside the tree.

    if (read_a_param) then
       !a parameters were read in
       call object_provide_in_node(lhere, 'read_a_param_phi')
    end if

    ! allocation
    call object_alloc('a_param_phi_all_elec', a_param_phi_all_elec, nb_a_param_phi, nctype, spin_dependencies_nb)

    !do we read in parameters?
    if (read_a_param) then
       !read in parameters
       a_param_spin_3_nb = (nb_a_param_phi) * nctype * 3         !  up_up != down_down != up_down
       a_param_spin_2_nb = (nb_a_param_phi) * nctype * 2         !  up_up = down_down != up_down
       a_param_spin_1_nb = (nb_a_param_phi) * nctype             !  up_up = down_down = up_down
       if (read_a_param_phi_nb == a_param_spin_3_nb) then
          !enough parameters were read in for:  up_up != down_down != up_down
          if (phi_spin_dependence .ne. 3) then
             !too many parameters were read in
             call require(lhere, 'number of parameters correct for phi_spin_dependence specified', phi_spin_dependence .eq. 3)
          else
             !correct number of parameters so reshape into correct dimensions for a_param_phi
             a_param_phi_all_elec = reshape( read_a_param_phi, (/ nb_a_param_phi, nctype, spin_dependencies_nb /)) !  up_up != down_down != up_down
             !don't want to re-read
             read_a_param = .false.
          end if
       elseif( read_a_param_phi_nb == a_param_spin_2_nb) then
          !enough parameters read in for:  up_up = down_down != up_down
          if ( phi_spin_dependence .ne. 2) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for phi_spin_dependence specified', phi_spin_dependence .eq. 2)
          else
             !correct number of parameters so reshape into correct dimensions for a_param_phi_all_elec
             a_param_phi_all_elec = reshape( read_a_param_phi, (/ nb_a_param_phi, nctype, spin_dependencies_nb /), pad=read_a_param_phi)   !up_up = down_down != up_down
             ! up_up is in a_param_phi_all_elec(:,:,1), up_down is in a_param_phi_all_elec(:,:,2) and a_param_phi_all_elec(:,:,3) has garbage in it
             !want up_up in a_param_phi_all_elec(:,:,1), up_up=down_down in a_param_phi_all_elec(:,:,2), up_down in a_param_phi_all_elec(:,:,3)
             a_param_phi_all_elec(:,:,3)=a_param_phi_all_elec(:,:,2)
             a_param_phi_all_elec(:,:,2)=a_param_phi_all_elec(:,:,1)
             !don't want to re-read
             read_a_param = .false.
          end if
       elseif( read_a_param_phi_nb == a_param_spin_1_nb) then
          !enough parameters read in for:  up_up = down_down = up_down
          if ( phi_spin_dependence .ne. 1) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for phi_spin_dependence specified', phi_spin_dependence .eq. 1)
          else
             !correct number of parameters so reshape into correct dimensions for a_param_phi_all_elec
             a_param_phi_all_elec = reshape( read_a_param_phi, (/ nb_a_param_phi, nctype, spin_dependencies_nb /), pad=read_a_param_phi)   !up_up = down_down = up_down
             ! up_up is in a_param_phi_all_elec(:,:,1), a_param_phi_all_elec(:,:,2) has garbage in it and a_param_phi_all_elec(:,:,3) has garbage in it
             !want up_up in a_param_phi_all_elec(:,:,1), up_up=down_down in a_param_phi_all_elec(:,:,2), up_up=up_down in a_param_phi_all_elec(:,:,3)
             a_param_phi_all_elec(:,:,2)=a_param_phi_all_elec(:,:,1)
             a_param_phi_all_elec(:,:,3)=a_param_phi_all_elec(:,:,1)
             !don't want to re-read
             read_a_param = .false.
          end if
       else
          !wrong number of parameters
          call die(lhere,'read_a_param_phi has an incorrect number of parameters.')
       end if
    else
       !don't read in parameters
       !start from random guess
       call random_number(a_param_phi_all_elec(:,:,:))
    endif



    !impose cusp
    !To impose the cusp we take the pivot positions in our reduced row echelon form of the matrix of cusp constraints
    !and we set this dependent parameter equal to the negative of the dot product of the rest of the row with
    !the corresponding a parameters. This is nothing fancy, but simply a standard way of solving linear equations

    !What is the spin dependence?
    if ( phi_spin_dependence .eq. 3) then
       !impose cusp parallel
       do dep_a_param_i=1, size(dep_a_param_phi_all_elec_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             !loop of parallel spins: up-up and down-down
             do spin_dependency_i=1, 2
                a_param_phi_all_elec(dep_a_param_phi_all_elec_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                     & dot_product( -a_param_phi_all_elec_parallel_cond(dep_a_param_i, dep_a_param_phi_all_elec_parallel(dep_a_param_i)+1:),&
                     & a_param_phi_all_elec(dep_a_param_phi_all_elec_parallel(dep_a_param_i) + 1:, ctype_i, spin_dependency_i))
             end do
          end do
       end do
       !impose cusp anti_parallel
       spin_dependency_i = 3 !up_down or down_up
       !loop over dependent parameters
       do dep_a_param_i=1, size(dep_a_param_phi_all_elec_anti_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             a_param_phi_all_elec(dep_a_param_phi_all_elec_anti_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -a_param_phi_all_elec_anti_parallel_cond(dep_a_param_i, dep_a_param_phi_all_elec_anti_parallel(dep_a_param_i)+1:),&
                  & a_param_phi_all_elec(dep_a_param_phi_all_elec_anti_parallel(dep_a_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
    elseif ( phi_spin_dependence .eq. 2) then
       !        up_up = down_down != up_down
       !impose cusp parallel
       spin_dependency_i = 1 !up-up
       !loop over dependent parameters
       do dep_a_param_i=1, size(dep_a_param_phi_all_elec_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             a_param_phi_all_elec(dep_a_param_phi_all_elec_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -a_param_phi_all_elec_parallel_cond(dep_a_param_i, dep_a_param_phi_all_elec_parallel(dep_a_param_i)+1:),&
                  & a_param_phi_all_elec(dep_a_param_phi_all_elec_parallel(dep_a_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !impose cusp anti_parallel
       spin_dependency_i = 3 !up_down or down_up
       !loop over dependent parameters
       do dep_a_param_i=1, size(dep_a_param_phi_all_elec_anti_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             a_param_phi_all_elec(dep_a_param_phi_all_elec_anti_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -a_param_phi_all_elec_anti_parallel_cond(dep_a_param_i, dep_a_param_phi_all_elec_anti_parallel(dep_a_param_i)+1:),&
                  & a_param_phi_all_elec(dep_a_param_phi_all_elec_anti_parallel(dep_a_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !require up-up = down-down
       a_param_phi_all_elec(:,:,2) = a_param_phi_all_elec(:,:,1)
    elseif( phi_spin_dependence .eq. 1) then
       !        up_up = down_down = up_down
       !impose cusp parallel
       spin_dependency_i = 1 !up-up
       !loop over dependent parameters
       do dep_a_param_i=1, size(dep_a_param_phi_all_elec_parallel,1)
          !loop over center types
          do ctype_i=1, nctype
             a_param_phi_all_elec(dep_a_param_phi_all_elec_parallel(dep_a_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -a_param_phi_all_elec_parallel_cond(dep_a_param_i, dep_a_param_phi_all_elec_parallel(dep_a_param_i)+1:),&
                  & a_param_phi_all_elec(dep_a_param_phi_all_elec_parallel(dep_a_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !require up-up = down-down=up-down
       a_param_phi_all_elec(:,:,2) = a_param_phi_all_elec(:,:,1)
       a_param_phi_all_elec(:,:,3) = a_param_phi_all_elec(:,:,1)
    end if




  end subroutine a_param_phi_all_elec_bld

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
    integer  :: order_i, order_j, order_k, order_k_max

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
    do order_i = 2, order_phi_bf
       do order_j = 0, order_i - 1
          if (order_j .ne. 0) then
             order_k_max = order_i - order_j
          else
             order_k_max = order_i - 2
          end if
          do order_k = 0, order_k_max
             if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
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

  subroutine dep_a_param_phi_pseudo_parallel_bld

    !---------------------------------------------------------------------------
    ! Description : build dep_a_param_phi_pseudo_parallel. This is a rank one array.
    !               This produces a list of the dependent parameters from a_param_phi_pseudo_parallel_cond.
    !               To do this we find the pivots in a row reduced matrix.
    !
    ! Created     : F. Petruzielo, 12 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer :: nb_dep_a_param_phi_pseudo_parallel
    integer :: row_i, pivot_column_i, column_i

    ! header
    if (header_exe) then

       call object_create('dep_a_param_phi_pseudo_parallel')

       !need number of a parameters
       call object_needed('nb_a_param_phi')
       !need parallel spin conditions
       call object_needed('a_param_phi_pseudo_parallel_cond')

       return
    endif


    !initialize the number of dependent parameters
    nb_dep_a_param_phi_pseudo_parallel = 0
    do row_i = 1, size(a_param_phi_pseudo_parallel_cond,1)
       if (sum(abs(a_param_phi_pseudo_parallel_cond(row_i,:))) < 1.d-10) then
          !no pivot
          cycle
       end if
       !pivot found, so increment counter
       nb_dep_a_param_phi_pseudo_parallel = nb_dep_a_param_phi_pseudo_parallel + 1
    end do

    ! allocation
    call object_alloc('dep_a_param_phi_pseudo_parallel', dep_a_param_phi_pseudo_parallel, nb_dep_a_param_phi_pseudo_parallel)


    !intialize the pivot column
    pivot_column_i = 1
    do row_i = 1, nb_dep_a_param_phi_pseudo_parallel
       do column_i = pivot_column_i, nb_a_param_phi
          !first non zero element  (should be a 1)
          if (a_param_phi_pseudo_parallel_cond(row_i,column_i) > 1.d-10) then
             !pivot found in this row
             dep_a_param_phi_pseudo_parallel(row_i) = column_i
             !next pivot column must be to the right of this one
             pivot_column_i=column_i + 1
             !so move on to next row
             exit
          end if
       end do
    end do

  end subroutine dep_a_param_phi_pseudo_parallel_bld

  !==============================================================================================




  !==============================================================================================

  subroutine dep_a_param_phi_all_elec_parallel_bld

    !---------------------------------------------------------------------------
    ! Description : build dep_a_param_phi_all_elec_parallel. This is a rank one array.
    !               This produces a list of the dependent parameters from a_param_phi_all_elec_parallel_cond.
    !               To do this we find the pivots in a row reduced matrix.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer :: nb_dep_a_param_phi_all_elec_parallel
    integer :: row_i, pivot_column_i, column_i

    ! header
    if (header_exe) then

       call object_create('dep_a_param_phi_all_elec_parallel')

       !need number of a parameters
       call object_needed('nb_a_param_phi')
       !need parallel spin conditions
       call object_needed('a_param_phi_all_elec_parallel_cond')

       return
    endif


    !initialize the number of dependent parameters
    nb_dep_a_param_phi_all_elec_parallel = 0
    do row_i = 1, size(a_param_phi_all_elec_parallel_cond,1)
       if (sum(abs(a_param_phi_all_elec_parallel_cond(row_i,:))) < 1.d-10) then
          !no pivot
          cycle
       end if
       !pivot found, so increment counter
       nb_dep_a_param_phi_all_elec_parallel = nb_dep_a_param_phi_all_elec_parallel + 1
    end do

    ! allocation
    call object_alloc('dep_a_param_phi_all_elec_parallel', dep_a_param_phi_all_elec_parallel, nb_dep_a_param_phi_all_elec_parallel)


    !intialize the pivot column
    pivot_column_i = 1
    do row_i = 1, nb_dep_a_param_phi_all_elec_parallel
       do column_i = pivot_column_i, nb_a_param_phi
          !first non zero element  (should be a 1)
          if (a_param_phi_all_elec_parallel_cond(row_i,column_i) > 1.d-10) then
             !pivot found in this row
             dep_a_param_phi_all_elec_parallel(row_i) = column_i
             !next pivot column must be to the right of this one
             pivot_column_i=column_i + 1
             !so move on to next row
             exit
          end if
       end do
    end do

  end subroutine dep_a_param_phi_all_elec_parallel_bld

  !==============================================================================================




  !==============================================================================================

  subroutine dep_a_param_phi_pseudo_anti_parallel_bld

    !---------------------------------------------------------------------------
    ! Description : build dep_a_param_phi_pseudo_anti_parallel. This is a rank one array.
    !               This produces a list of the dependent parameters from a_param_phi_pseudo_anti_parallel_cond.
    !               To do this we find the pivots in a row reduced matrix.
    !
    ! Created     : F. Petruzielo, 13 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer :: nb_dep_a_param_phi_pseudo_anti_parallel
    integer :: row_i, pivot_column_i, column_i

    ! header
    if (header_exe) then

       call object_create('dep_a_param_phi_pseudo_anti_parallel')

       !need number of a parameters
       call object_needed('nb_a_param_phi')
       !need anti parallel spin conditions
       call object_needed('a_param_phi_pseudo_anti_parallel_cond')

       return
    endif


    !initialize the number of dependent parameters
    nb_dep_a_param_phi_pseudo_anti_parallel = 0
    do row_i = 1, size(a_param_phi_pseudo_anti_parallel_cond,1)
       if (sum(abs(a_param_phi_pseudo_anti_parallel_cond(row_i,:))) < 1.d-10) then
          !no pivot
          cycle
       end if
       !pivot found, so increment counter
       nb_dep_a_param_phi_pseudo_anti_parallel = nb_dep_a_param_phi_pseudo_anti_parallel + 1
    end do

    ! allocation
    call object_alloc('dep_a_param_phi_pseudo_anti_parallel', dep_a_param_phi_pseudo_anti_parallel, nb_dep_a_param_phi_pseudo_anti_parallel)


    !intialize the pivot column
    pivot_column_i = 1
    do row_i = 1, nb_dep_a_param_phi_pseudo_anti_parallel
       do column_i = pivot_column_i, nb_a_param_phi
          !first non zero element  (should be a 1)
          if (a_param_phi_pseudo_anti_parallel_cond(row_i,column_i) > 1.d-10) then
             !pivot found in this row
             dep_a_param_phi_pseudo_anti_parallel(row_i) = column_i
             !next pivot column must be to the right of this one
             pivot_column_i=column_i + 1
             !so move on to next row
             exit
          end if
       end do
    end do

  end subroutine dep_a_param_phi_pseudo_anti_parallel_bld

  !==============================================================================================





  !==============================================================================================

  subroutine dep_a_param_phi_all_elec_anti_parallel_bld

    !---------------------------------------------------------------------------
    ! Description : build dep_a_param_phi_all_elec_anti_parallel. This is a rank one array.
    !               This produces a list of the dependent parameters from a_param_phi_all_elec_anti_parallel_cond.
    !               To do this we find the pivots in a row reduced matrix.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer :: nb_dep_a_param_phi_all_elec_anti_parallel
    integer :: row_i, pivot_column_i, column_i

    ! header
    if (header_exe) then

       call object_create('dep_a_param_phi_all_elec_anti_parallel')

       !need number of a parameters
       call object_needed('nb_a_param_phi')
       !need anti parallel spin conditions
       call object_needed('a_param_phi_all_elec_anti_parallel_cond')

       return
    endif


    !initialize the number of dependent parameters
    nb_dep_a_param_phi_all_elec_anti_parallel = 0
    do row_i = 1, size(a_param_phi_all_elec_anti_parallel_cond,1)
       if (sum(abs(a_param_phi_all_elec_anti_parallel_cond(row_i,:))) < 1.d-10) then
          !no pivot
          cycle
       end if
       !pivot found, so increment counter
       nb_dep_a_param_phi_all_elec_anti_parallel = nb_dep_a_param_phi_all_elec_anti_parallel + 1
    end do

    ! allocation
    call object_alloc('dep_a_param_phi_all_elec_anti_parallel', dep_a_param_phi_all_elec_anti_parallel, nb_dep_a_param_phi_all_elec_anti_parallel)


    !intialize the pivot column
    pivot_column_i = 1
    do row_i = 1, nb_dep_a_param_phi_all_elec_anti_parallel
       do column_i = pivot_column_i, nb_a_param_phi
          !first non zero element  (should be a 1)
          if (a_param_phi_all_elec_anti_parallel_cond(row_i,column_i) > 1.d-10) then
             !pivot found in this row
             dep_a_param_phi_all_elec_anti_parallel(row_i) = column_i
             !next pivot column must be to the right of this one
             pivot_column_i=column_i + 1
             !so move on to next row
             exit
          end if
       end do
    end do

  end subroutine dep_a_param_phi_all_elec_anti_parallel_bld

  !==============================================================================================





  !==============================================================================================

  subroutine a_param_phi_pseudo_anti_parallel_cond_bld

    !---------------------------------------------------------------------------
    ! Description : build a_param_phi_pseudo_anti_parallel_spin_cond_bld. This is a rank two array.
    !               The first index is the number of condition ( some are possibly redundant).
    !               The second index is the number a parameters.
    !               We use gaussian elimination to row reduce the conditions and check for redundancies.
    !
    !               cusp condition constraints:   d phi  |
    !                                             -----  |        = 0
    !                                             d r_iI |
    !                                                     r_iI=0
    !
    ! Created     : F. Petruzielo, 12 Aug 2008
    !---------------------------------------------------------------------------

    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_i, order_j, order_k, order_k_max, param_i

    ! header
    if (header_exe) then

       call object_create('a_param_phi_pseudo_anti_parallel_cond')

       !need order of expansion
       call object_needed('order_phi_bf')
       !need number of a parameters
       call object_needed('nb_a_param_phi')

       return
    endif

    ! allocation
    ! one derivative constraints so number of contraints is (order_phi_bf -1 )
    call object_alloc('a_param_phi_pseudo_anti_parallel_cond', a_param_phi_pseudo_anti_parallel_cond, order_phi_bf - 1, nb_a_param_phi)

    !initialize
    a_param_phi_pseudo_anti_parallel_cond(:,:) = 0.d0

    param_i = 0 !intialize parameter counter
    !setup matrix of constraints
    do order_i = 2, order_phi_bf
       do order_j = 0, order_i - 1
          if (order_j .ne. 0) then
             order_k_max = order_i - order_j
          else
             order_k_max = order_i - 2
          end if
          do order_k = 0, order_k_max
             if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                cycle
             else
                param_i = param_i + 1
                !=================================
                ! d phi  |
                ! -----  |        = 0
                ! d r_iI |
                !         r_iI=0
                !=================================
                a_param_phi_pseudo_anti_parallel_cond((order_i+order_j-order_k)/2, param_i) = a_param_phi_pseudo_anti_parallel_cond( (order_i+order_j-order_k)/2, param_i) + (order_k + order_i - order_j) / 2
                a_param_phi_pseudo_anti_parallel_cond((order_i+order_j+order_k)/2, param_i) = a_param_phi_pseudo_anti_parallel_cond( (order_i+order_j+order_k)/2, param_i) + ( order_i - order_j - order_k) / 2
             endif
          enddo
       enddo
    enddo

    !use gaussian elimintation on constraints to remove redundancies
    call row_reduce(a_param_phi_pseudo_anti_parallel_cond)

  end subroutine a_param_phi_pseudo_anti_parallel_cond_bld

  !==============================================================================================

  !==============================================================================================

  subroutine a_param_phi_pseudo_parallel_cond_bld

    !---------------------------------------------------------------------------
    ! Description : build a_param_phi_pseudo_parallel_spin_cond_bld. This is a rank two array.
    !               The first index is the number of condition ( some are possibly redundant).
    !               The second index is the number a parameters.
    !               We use gaussian elimination to row reduce the conditions and check for redundancies.
    !
    !               cusp condition constraints:   d phi  |
    !                                             -----  |        = 0
    !                                             d r_iI |
    !                                                     r_iI=
    !                                             d phi  |
    !                                             -----  |        = 0
    !                                             d r_ij |
    !                                                     r_ij=0
    !
    ! Created     : F. Petruzielo, 12 Aug 2008
    !---------------------------------------------------------------------------

    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_i, order_j, order_k, order_k_max, param_i

    ! header
    if (header_exe) then

       call object_create('a_param_phi_pseudo_parallel_cond')

       !need order of expansion
       call object_needed('order_phi_bf')
       !need number of a parameters
       call object_needed('nb_a_param_phi')

       return
    endif

    ! allocation
    ! two derivative constraints so number of contraints is 2*(order_phi_bf -1 )
    call object_alloc('a_param_phi_pseudo_parallel_cond', a_param_phi_pseudo_parallel_cond, 2*(order_phi_bf - 1), nb_a_param_phi)


    !initialize
    a_param_phi_pseudo_parallel_cond(:,:) = 0.d0

    param_i = 0 !intialize parameter counter
    !setup matrix of constraints
    do order_i = 2, order_phi_bf
       do order_j = 0, order_i - 1
          if (order_j .ne. 0) then
             order_k_max = order_i - order_j
          else
             order_k_max = order_i - 2
          end if
          do order_k = 0, order_k_max
             if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                cycle
             else
                param_i = param_i + 1
                !=================================
                ! d phi  |
                ! -----  |        = 0
                ! d r_iI |
                !         r_iI=0
                !=================================
                a_param_phi_pseudo_parallel_cond((order_i+order_j-order_k)/2, param_i) = a_param_phi_pseudo_parallel_cond( (order_i+order_j-order_k)/2, param_i) + (order_k + order_i - order_j) / 2
                a_param_phi_pseudo_parallel_cond((order_i+order_j+order_k)/2, param_i) = a_param_phi_pseudo_parallel_cond( (order_i+order_j+order_k)/2, param_i) + ( order_i - order_j - order_k) / 2
                !=================================
                ! d phi  |
                ! -----  |        = 0
                ! d r_ij |
                !         r_ij=0
                !=================================
                !note that the row needs to be shifted down by order_phi_bf - 1 because those rows hold the above constraints
                a_param_phi_pseudo_parallel_cond(order_phi_bf - 1 + order_i-order_j, param_i) = a_param_phi_pseudo_parallel_cond(order_phi_bf - 1 + order_i-order_j, param_i) + order_j
             endif
          enddo
       enddo
    enddo

    !use gaussian elimintation on constraints to remove redundancies
    call row_reduce(a_param_phi_pseudo_parallel_cond)

  end subroutine a_param_phi_pseudo_parallel_cond_bld

  !==============================================================================================


  !==============================================================================================

  subroutine a_param_phi_all_elec_parallel_cond_bld

    !---------------------------------------------------------------------------
    ! Description : build a_param_phi_all_elec_parallel_spin_cond_bld. This is a rank two array.
    !               The first index is the number of condition ( some are possibly redundant).
    !               The second index is the number a parameters.
    !               We use gaussian elimination to row reduce the conditions and check for redundancies.
    !
    !               cusp condition constraints:   d phi  |
    !                                             -----  |        = 0
    !                                             d r_iI |
    !                                                     r_iI=0
    !                                             d phi  |
    !                                             -----  |        = 0
    !                                             d r_ij |
    !                                                     r_ij=0
    !
    !                                              phi  |
    !                                                   |        = 0
    !                                                   |
    !                                                     r_iI=0
    !
    !
    !                                             d phi  |
    !                                             -----  |        = 0      some of these are redundant conditions
    !                                             d r_ij |
    !                                                     r_iI=0

    !
    !
    !
    !
    ! Created     : F. Petruzielo, 13 Aug 2008
    !---------------------------------------------------------------------------

    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_i, order_j, order_k, order_k_max, param_i

    ! header
    if (header_exe) then

       call object_create('a_param_phi_all_elec_parallel_cond')

       !need order of expansion
       call object_needed('order_phi_bf')
       !need number of a parameters
       call object_needed('nb_a_param_phi')

       return
    endif

    ! allocation
    ! two derivative constraints evaluated at 0, one constraint on function
    ! and one derivative constraint evaluated at another point so number of contraints is 2*(order_phi_bf -1 ) + order_phi_bf + order_phi_bf
    call object_alloc('a_param_phi_all_elec_parallel_cond', a_param_phi_all_elec_parallel_cond, 4 * order_phi_bf - 2, nb_a_param_phi)

    !initialize
    a_param_phi_all_elec_parallel_cond(:,:) = 0.d0

    param_i = 0 !intialize parameter counter
    !setup matrix of constraints
    do order_i = 2, order_phi_bf
       do order_j = 0, order_i - 1
          if (order_j .ne. 0) then
             order_k_max = order_i - order_j
          else
             order_k_max = order_i - 2
          end if
          do order_k = 0, order_k_max
             if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                cycle
             else
                param_i = param_i + 1
                !=================================
                ! d phi  |
                ! -----  |        = 0
                ! d r_iI |
                !         r_iI=0
                !=================================
                a_param_phi_all_elec_parallel_cond((order_i+order_j-order_k)/2, param_i) = a_param_phi_all_elec_parallel_cond( (order_i+order_j-order_k)/2, param_i) + (order_k + order_i - order_j) / 2
                a_param_phi_all_elec_parallel_cond((order_i+order_j+order_k)/2, param_i) = a_param_phi_all_elec_parallel_cond( (order_i+order_j+order_k)/2, param_i) + ( order_i - order_j - order_k) / 2
                !=================================
                ! d phi  |
                ! -----  |        = 0
                ! d r_ij |
                !         r_ij=0
                !=================================
                !note that the row needs to be shifted down by order_phi_bf - 1 because those rows hold the above constraints
                a_param_phi_all_elec_parallel_cond(order_phi_bf - 1 + order_i-order_j, param_i) = a_param_phi_all_elec_parallel_cond(order_phi_bf - 1 + order_i-order_j, param_i) + order_j
                !================================================
                ! phi  |
                !      |        = 0
                !      |
                !       r_iI=0
                !===============================================
                !note that the row needs to be shifted down by 2*(order_phi_bf - 1) because those rows hold the above constraints
                a_param_phi_all_elec_parallel_cond(2 *(order_phi_bf - 1) + (order_i+order_j-order_k)/2, param_i) = a_param_phi_all_elec_parallel_cond(2 * (order_phi_bf -1 ) + (order_i+order_j-order_k)/2, param_i) + 1.d0
                a_param_phi_all_elec_parallel_cond(2 *(order_phi_bf - 1) + (order_i+order_j+order_k)/2, param_i) = a_param_phi_all_elec_parallel_cond(2 * (order_phi_bf -1 ) + (order_i+order_j+order_k)/2, param_i) + 1.d0
                !================================================
                !d phi  |
                ! ----  |        = 0
                !d r_ij |
                !       r_iI=0
                !===============================================
                !note that the row needs to be shifted down by 2*(order_phi_bf - 1) + order_phi_bf because those rows hold the above constraints
                a_param_phi_all_elec_parallel_cond(3*order_phi_bf-2+(order_i+order_j-2-order_k)/2, param_i) = a_param_phi_all_elec_parallel_cond(3*order_phi_bf-2+(order_i+order_j-2-order_k)/2, param_i) + order_j
                a_param_phi_all_elec_parallel_cond(3*order_phi_bf-2+(order_i+order_j-2+order_k)/2, param_i) = a_param_phi_all_elec_parallel_cond(3*order_phi_bf-2+(order_i+order_j-2+order_k)/2, param_i) + order_j
             endif
          enddo
       enddo
    enddo

    !use gaussian elimintation on constraints to remove redundancies
    call row_reduce(a_param_phi_all_elec_parallel_cond)

  end subroutine a_param_phi_all_elec_parallel_cond_bld

  !==============================================================================================



  !==============================================================================================

  subroutine a_param_phi_all_elec_anti_parallel_cond_bld

    !---------------------------------------------------------------------------
    ! Description : build a_param_phi_all_elec_anti_parallel_spin_cond_bld. This is a rank two array.
    !               The first index is the number of condition ( some are possibly redundant).
    !               The second index is the number a parameters.
    !               We use gaussian elimination to row reduce the conditions and check for redundancies.
    !
    !               cusp condition constraints:   d phi  |
    !                                             -----  |        = 0
    !                                             d r_iI |
    !                                                     r_iI=0
    !
    !                                              phi  |
    !                                                   |        = 0
    !                                                   |
    !                                                     r_iI=0
    !
    !
    !                                             d phi  |
    !                                             -----  |        = 0      some of these are redundant conditions
    !                                             d r_ij |
    !                                                     r_iI=0

    !
    !
    !
    !
    ! Created     : F. Petruzielo, 13 Aug 2008
    !---------------------------------------------------------------------------

    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_i, order_j, order_k, order_k_max, param_i

    ! header
    if (header_exe) then

       call object_create('a_param_phi_all_elec_anti_parallel_cond')

       !need order of expansion
       call object_needed('order_phi_bf')
       !need number of a parameters
       call object_needed('nb_a_param_phi')

       return
    endif

    ! allocation
    ! one derivative constraints evaluated at 0, one constraint on function
    ! and one derivative constraint evaluated at another point so number of contraints is (order_phi_bf -1 ) + order_phi_bf + order_phi_bf
    call object_alloc('a_param_phi_all_elec_anti_parallel_cond', a_param_phi_all_elec_anti_parallel_cond, 3 * order_phi_bf - 1, nb_a_param_phi)

    !initialize
    a_param_phi_all_elec_anti_parallel_cond(:,:) = 0.d0

    param_i = 0 !intialize parameter counter
    !setup matrix of constraints
    do order_i = 2, order_phi_bf
       do order_j = 0, order_i - 1
          if (order_j .ne. 0) then
             order_k_max = order_i - order_j
          else
             order_k_max = order_i - 2
          end if
          do order_k = 0, order_k_max
             if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                cycle
             else
                param_i = param_i + 1
                !=================================
                ! d phi  |
                ! -----  |        = 0
                ! d r_iI |
                !         r_iI=0
                !=================================
                a_param_phi_all_elec_anti_parallel_cond((order_i+order_j-order_k)/2, param_i) = a_param_phi_all_elec_anti_parallel_cond( (order_i+order_j-order_k)/2, param_i) + (order_k + order_i - order_j) / 2
                a_param_phi_all_elec_anti_parallel_cond((order_i+order_j+order_k)/2, param_i) = a_param_phi_all_elec_anti_parallel_cond( (order_i+order_j+order_k)/2, param_i) + ( order_i - order_j - order_k) / 2
                !================================================
                ! phi  |
                !      |        = 0
                !      |
                !       r_iI=0
                !===============================================
                !note that the row needs to be shifted down by (order_phi_bf - 1) because those rows hold the above constraints
                a_param_phi_all_elec_anti_parallel_cond((order_phi_bf - 1) + (order_i+order_j-order_k)/2, param_i) = a_param_phi_all_elec_anti_parallel_cond((order_phi_bf -1 ) + (order_i+order_j-order_k)/2, param_i) + 1.d0
                a_param_phi_all_elec_anti_parallel_cond((order_phi_bf - 1) + (order_i+order_j+order_k)/2, param_i) = a_param_phi_all_elec_anti_parallel_cond((order_phi_bf -1 ) + (order_i+order_j+order_k)/2, param_i) + 1.d0
                !================================================
                !d phi  |
                ! ----  |        = 0
                !d r_ij |
                !       r_iI=0
                !===============================================
                !note that the row needs to be shifted down by (order_phi_bf - 1) + order_phi_bf because those rows hold the above constraints
                a_param_phi_all_elec_anti_parallel_cond(2*order_phi_bf-1+(order_i+order_j-2-order_k)/2, param_i) = a_param_phi_all_elec_anti_parallel_cond(2*order_phi_bf-1+(order_i+order_j-2-order_k)/2, param_i) + order_j
                a_param_phi_all_elec_anti_parallel_cond(2*order_phi_bf-1+(order_i+order_j-2+order_k)/2, param_i) = a_param_phi_all_elec_anti_parallel_cond(2*order_phi_bf-1+(order_i+order_j-2+order_k)/2, param_i) + order_j
             endif
          enddo
       enddo
    enddo

    !use gaussian elimintation on constraints to remove redundancies
    call row_reduce(a_param_phi_all_elec_anti_parallel_cond)

  end subroutine a_param_phi_all_elec_anti_parallel_cond_bld

  !==============================================================================================


  !======================================================================================

  subroutine xi_een_theta_bld
    !---------------------------------------------------------------------------
    ! Description : build object xi_een_theta(:,:). This is a rank three array.
    !               It is part of electron-electron-nuclear backflow.
    !               It is necessary for the construction of the total backflow transformation.
    !               The first index has dimension equal to the number of spatial dimensions.
    !               The second index has dimensions equal to the number of electrons.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, nelec, ndim
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, cent_i, walk_i, dim_i
    character(len=max_string_len_rout), save :: lhere = 'xi_een_theta_bld'

    ! header
    if (header_exe) then

       call object_create('xi_een_theta')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       ! number of electrons
       call object_needed('nelec')
       !cartesian vector between electrons and nuclei (for each walker)
       call object_needed('vec_en_xyz_wlk')
       !dimension of space
       call object_needed('ndim')
       !is this an all-electron calculation?
       call object_needed('all_elec')

       return
    endif

    if (all_elec) then
       !all electron case
       call object_provide_in_node(lhere, 'theta_all_elec')

       !allocation

       call object_alloc('xi_een_theta', xi_een_theta, ndim, nelec, nwalk)

       !build xi_een_theta
       xi_een_theta(:,:,:) = 0.d0

       !loop over walkers
       do walk_i = 1, nwalk
          !loop over nuclei
          do cent_i = 1, ncent
             !loop over electrons
             do elec_i = 1, nelec
                !loop over other electrons before elec_i
                do elec_j = 1, elec_i -1
                   !loop over spatial dimensions
                   do dim_i = 1, ndim
                      xi_een_theta(dim_i, elec_i, walk_i) = xi_een_theta(dim_i, elec_i, walk_i) +  theta_all_elec(elec_i, elec_j, cent_i, walk_i) * vec_en_xyz_wlk(dim_i, elec_i, cent_i, walk_i)
                   enddo
                enddo
                !loop over other electrons after elec_i
                do elec_j = elec_i + 1, nelec
                   !loop over spatial dimensions
                   do dim_i = 1, ndim
                      xi_een_theta(dim_i, elec_i, walk_i) = xi_een_theta(dim_i, elec_i, walk_i) +  theta_all_elec(elec_i, elec_j, cent_i, walk_i) * vec_en_xyz_wlk(dim_i, elec_i, cent_i, walk_i)
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       !pseudo case
       call object_provide_in_node(lhere, 'theta_pseudo')

       !allocation

       call object_alloc('xi_een_theta', xi_een_theta, ndim, nelec, nwalk)

       !build xi_en
       xi_een_theta(:,:,:) = 0.d0

       !loop over walkers
       do walk_i = 1, nwalk
          !loop over nuclei
          do cent_i = 1, ncent
             !loop over electrons
             do elec_i = 1, nelec
                !loop over other electrons before elec_i
                do elec_j = 1, elec_i -1
                   !loop over spatial dimensions
                   do dim_i = 1, ndim
                      xi_een_theta(dim_i, elec_i, walk_i) = xi_een_theta(dim_i, elec_i, walk_i) +  theta_pseudo(elec_i, elec_j, cent_i, walk_i) * vec_en_xyz_wlk(dim_i, elec_i, cent_i, walk_i)
                   enddo
                enddo
                !loop over other electrons after elec_i
                do elec_j = elec_i + 1, nelec
                   !loop over spatial dimensions
                   do dim_i = 1, ndim
                      xi_een_theta(dim_i, elec_i, walk_i) = xi_een_theta(dim_i, elec_i, walk_i) +  theta_pseudo(elec_i, elec_j, cent_i, walk_i) * vec_en_xyz_wlk(dim_i, elec_i, cent_i, walk_i)
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif


  end subroutine xi_een_theta_bld

  !======================================================================================



  !======================================================================================

  subroutine theta_all_elec_bld
    !---------------------------------------------------------------------------
    ! Description : build object theta_all_elec(:,:,:,:). This is a rank four array.
    !               It is necessary for the construction of the electron-electron-nuclear backflow terms in all electron case.
    !               The first index has dimension equal to the number of electrons.
    !               The second index has dimension equal to the number of electrons.
    !               The third index has dimensions equal to the number of nuclei.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, cent_i, walk_i, order_i, order_j, order_k, prod_i, param_i
    integer  :: order_k_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('theta_all_elec')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       !specify atom-type for each atom
       call object_needed('iwctype')
       ! number of electrons
       call object_needed('nelec')
       !number of up-spin electrons
       call object_needed('nup')
       !order of expansion
       call object_needed('order_theta_bf')
       !scaled distances of each electron of each walker to each nucleus (three body scaling functions)
       call object_needed ('scaled_dist_een_en_wlk')
       !scaled distances of each electron of each walker to every other electron of that walker (three body scaling functions)
       call object_needed ('scaled_dist_een_ee_wlk')
       !list of a parameters (satisfying all-elec conditions) for theta backflow function for each center type, and each spin dependence
       call object_needed('b_param_theta_all_elec')
       !In all electron calculations, to impose the electron nuclear cusp, it is necessary to smoothly cutoff the backflow functions at the nucleus.
       call object_needed('smooth_cutoff_g')


       return
    endif

    ! allocation

    call object_alloc('theta_all_elec', theta_all_elec, nelec, nelec, ncent, nwalk)

    ! build theta_all_elec

    !initialize
    theta_all_elec (:,:,:,:) = 0.d0

    !up_up terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !loop over other up-spin electrons below elec_i
             do elec_j=1, elec_i - 1
                !power series terms of theta_all_elec as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_all_elec(elec_i, elec_j, cent_i, walk_i) = theta_all_elec(elec_i, elec_j, cent_i, walk_i) + b_param_theta_all_elec(param_i, iwctype(cent_i), up_up) &
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k &
                                 & +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) &
                                 & *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !loop over other up-spin electrons above elec_i
             do elec_j= elec_i + 1, nup
                !power series terms of theta_all_elec
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_all_elec(elec_i, elec_j, cent_i, walk_i) = theta_all_elec(elec_i, elec_j, cent_i, walk_i) + b_param_theta_all_elec(param_i, iwctype(cent_i), up_up) &
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k&
                                 & +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)&
                                 & *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo



    !down_down terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons
          do elec_i = nup + 1, nelec
             !loop over other down-spin electrons below elec_i
             do elec_j=nup + 1, elec_i - 1
                !power series terms of theta_all_elec as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_all_elec(elec_i, elec_j, cent_i, walk_i) = theta_all_elec(elec_i, elec_j, cent_i, walk_i) + b_param_theta_all_elec(param_i, iwctype(cent_i), down_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !loop over other down-spin electrons above elec_i
             do elec_j= elec_i + 1, nelec
                !power series terms of theta_all_elec
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_all_elec(elec_i, elec_j, cent_i, walk_i) = theta_all_elec(elec_i, elec_j, cent_i, walk_i) + b_param_theta_all_elec(param_i, iwctype(cent_i), down_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k&
                                 & +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)&
                                 & *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


    !up_down terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !loop over down-spin electrons
             do elec_j= nup + 1, nelec
                !power series terms of theta_all_elec as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_all_elec(elec_i, elec_j, cent_i, walk_i) = theta_all_elec(elec_i, elec_j, cent_i, walk_i) + b_param_theta_all_elec(param_i, iwctype(cent_i), up_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


    !down_up terms (note that b_param_theta_all_elec are same for up_down and down_up)
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons
          do elec_i = nup + 1, nelec
             !loop over up-spin electrons
             do elec_j = 1, nup
                !power series terms of theta_all_elec as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_all_elec(elec_i, elec_j, cent_i, walk_i) = theta_all_elec(elec_i, elec_j, cent_i, walk_i) + b_param_theta_all_elec(param_i, iwctype(cent_i), up_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    !apply smooth cutoff function g
    do cent_i=1, ncent
       do elec_j=1, nelec
          do prod_i=1, cent_i-1
             theta_all_elec(:,elec_j,cent_i,:) = theta_all_elec(:,elec_j,cent_i,:) * smooth_cutoff_g(:,prod_i,:)
          end do
          do prod_i=cent_i+1, ncent
             theta_all_elec(:,elec_j, cent_i,:) = theta_all_elec(:,elec_j, cent_i,:) * smooth_cutoff_g(:,prod_i,:)
          end do
       end do
    enddo

  end subroutine theta_all_elec_bld

  !=============================================================================================

  !======================================================================================

  subroutine theta_pseudo_bld
    !---------------------------------------------------------------------------
    ! Description : build object theta_pseudo(:,:,:,:). This is a rank four array.
    !               It is necessary for the construction of the electron-electron-nuclear backflow terms in the pseudo case.
    !               The first index has dimension equal to the number of electrons.
    !               The second index has dimension equal to the number of electrons.
    !               The third index has dimensions equal to the number of nuclei.
    !               The third index has dimension equal to the number of walkers.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nwalk, ncent, iwctype, nelec, nup
    include 'commons.h'

    !local
    integer  :: elec_i, elec_j, cent_i, walk_i, order_i, order_j, order_k, param_i
    integer  :: order_k_max
    integer, parameter  :: up_up=1, down_down=2, up_down=3 !index for the spin

    ! header
    if (header_exe) then

       call object_create('theta_pseudo')

       !number of walkers. nwalk=1 for vmc, but can be larger for dmc
       call object_needed('nwalk')
       !number of atoms/centers
       call object_needed('ncent')
       !specify atom-type for each atom
       call object_needed('iwctype')
       ! number of electrons
       call object_needed('nelec')
       !number of up-spin electrons
       call object_needed('nup')
       !order of expansion
       call object_needed('order_theta_bf')
       !scaled distances of each electron of each walker to each nucleus (three body scaling functions)
       call object_needed ('scaled_dist_een_en_wlk')
       !scaled distances of each electron of each walker to every other electron of that walker (three body scaling functions)
       call object_needed ('scaled_dist_een_ee_wlk')
       !list of a parameters (satisfying all-elec conditions) for theta backflow function for each center type, and each spin dependence
       call object_needed('b_param_theta_pseudo')

       return
    endif

    ! allocation

    call object_alloc('theta_pseudo', theta_pseudo, nelec, nelec, ncent, nwalk)

    ! build theta_pseudo

    !initialize
    theta_pseudo (:,:,:,:) = 0.d0

    !up_up terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !loop over other up-spin electrons below elec_i
             do elec_j=1, elec_i - 1
                !power series terms of theta_pseudo as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_pseudo(elec_i, elec_j, cent_i, walk_i) = theta_pseudo(elec_i, elec_j, cent_i, walk_i) + b_param_theta_pseudo(param_i, iwctype(cent_i), up_up)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !loop over other up-spin electrons above elec_i
             do elec_j= elec_i + 1, nup
                !power series terms of theta_pseudo
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_pseudo(elec_i, elec_j, cent_i, walk_i) = theta_pseudo(elec_i, elec_j, cent_i, walk_i) + b_param_theta_pseudo(param_i, iwctype(cent_i), up_up)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo



    !down_down terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons
          do elec_i = nup + 1, nelec
             !loop over other down-spin electrons below elec_i
             do elec_j=nup + 1, elec_i - 1
                !power series terms of theta_pseudo as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_pseudo(elec_i, elec_j, cent_i, walk_i) = theta_pseudo(elec_i, elec_j, cent_i, walk_i) + b_param_theta_pseudo(param_i, iwctype(cent_i), down_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !loop over other down-spin electrons above elec_i
             do elec_j= elec_i + 1, nelec
                !power series terms of theta_pseudo
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_pseudo(elec_i, elec_j, cent_i, walk_i) = theta_pseudo(elec_i, elec_j, cent_i, walk_i) + b_param_theta_pseudo(param_i, iwctype(cent_i), down_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


    !up_down terms
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over up-spin electrons
          do elec_i = 1, nup
             !loop over down-spin electrons
             do elec_j= nup + 1, nelec
                !power series terms of theta_pseudo as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_pseudo(elec_i, elec_j, cent_i, walk_i) = theta_pseudo(elec_i, elec_j, cent_i, walk_i) + b_param_theta_pseudo(param_i, iwctype(cent_i), up_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


    !down_up terms (note that b_param_theta_pseudo are same for up_down and down_up)
    do walk_i = 1, nwalk
       !loop over number of centers
       do cent_i = 1, ncent
          !loop over down-spin electrons
          do elec_i = nup + 1, nelec
             !loop over up-spin electrons
             do elec_j = 1, nup
                !power series terms of theta_pseudo as seen in three-body jastrow
                param_i = 1 !intialize parameter counter
                do order_i = 2, order_theta_bf
                   do order_j = 0, order_i - 1
                      if (order_j .ne. 0) then
                         order_k_max = order_i - order_j
                      else
                         order_k_max = order_i - 2
                      end if
                      do order_k = 0, order_k_max
                         if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                            cycle
                         else
                            theta_pseudo(elec_i, elec_j, cent_i, walk_i) = theta_pseudo(elec_i, elec_j, cent_i, walk_i) + b_param_theta_pseudo(param_i, iwctype(cent_i), up_down)&
                                 & * ( scaled_dist_een_ee_wlk(elec_i, elec_j, walk_i)**order_j ) * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i)**order_k +  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i)**order_k )&
                                 & * ( scaled_dist_een_en_wlk(elec_i, cent_i, walk_i) *  scaled_dist_een_en_wlk(elec_j, cent_i, walk_i) ) ** ((order_i - order_j - order_k) / 2)
                            param_i = param_i + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    !no cutoff because pseudo

  end subroutine theta_pseudo_bld

  !=============================================================================================



  !==============================================================================================

  subroutine b_param_theta_pseudo_bld

    !---------------------------------------------------------------------------
    ! Description : build b_param_theta_pseudo. This is a rank 3 array.
    !               List of b parameters for theta backflow function for each center type, and each spin dependence
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
    integer  :: b_param_spin_3_nb, b_param_spin_2_nb, b_param_spin_1_nb, dep_b_param_i, ctype_i, spin_dependency_i
    character(len=max_string_len_rout), save :: lhere = 'b_param_theta_pseudo_bld'
    integer, parameter :: spin_dependencies_nb = 3    !because three choices up_up, down_down, up_down

    ! header
    if (header_exe) then

       call object_create('b_param_theta_pseudo')

       !order of expansion
       call object_needed('order_theta_bf')
       !how many types of atoms
       call object_needed('nctype')
       !what is the spin dependence?
       call object_needed('theta_spin_dependence')
       !have the parameters been read in?
       call object_needed('read_b_param')
       !how many parameters were read in?
       call object_needed('read_b_param_theta_nb')
       !what are the dependent parameters?
       call object_needed('dep_b_param_theta_pseudo')
       !what are the conditions?
       call object_needed('b_param_theta_pseudo_cond')
       !How many a parameters
       call object_needed('nb_b_param_theta')

       return
    endif


    !Note that we need to have the option of using read_b_param_theta if it is read in
    !However, we don't know apriori if it will be read in.
    !Thus we don't call object_needed('read_b_param_theta')
    !To accomplish this correctly we need to use the utility object_provide_in_node
    !We cant simply use object_provide because we are inside the tree.

    if (read_b_param) then
       !a parameters were read in
       call object_provide_in_node(lhere, 'read_b_param_theta')
    end if

    ! allocation
    call object_alloc('b_param_theta_pseudo', b_param_theta_pseudo, nb_b_param_theta, nctype, spin_dependencies_nb)

    !do we read in parameters?
    if (read_b_param) then
       !read in parameters
       b_param_spin_3_nb = (nb_b_param_theta) * nctype * 3         !  up_up != down_down != up_down
       b_param_spin_2_nb = (nb_b_param_theta) * nctype * 2         !  up_up = down_down != up_down
       b_param_spin_1_nb = (nb_b_param_theta) * nctype             !  up_up = down_down = up_down
       if (read_b_param_theta_nb == b_param_spin_3_nb) then
          !enough parameters were read in for:  up_up != down_down != up_down
          if (theta_spin_dependence .ne. 3) then
             !too many parameters were read in
             call require(lhere, 'number of parameters correct for theta_spin_dependence specified', theta_spin_dependence .eq. 3)
          else
             !correct number of parameters so reshape into correct dimensions for b_param_theta
             b_param_theta_pseudo = reshape( read_b_param_theta, (/ nb_b_param_theta, nctype, spin_dependencies_nb /)) !  up_up != down_down != up_down
             !don't want to re-read
             read_b_param = .false.
          end if
       elseif( read_b_param_theta_nb == b_param_spin_2_nb) then
          !enough parameters read in for:  up_up = down_down != up_down
          if ( theta_spin_dependence .ne. 2) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for theta_spin_dependence specified', theta_spin_dependence .eq. 2)
          else
             !correct number of parameters so reshape into correct dimensions for b_param_theta_pseudo
             b_param_theta_pseudo = reshape( read_b_param_theta, (/ nb_b_param_theta, nctype, spin_dependencies_nb /), pad=read_b_param_theta)   !up_up = down_down != up_down
             ! up_up is in b_param_theta_pseudo(:,:,1), up_down is in b_param_theta_pseudo(:,:,2) and b_param_theta_pseudo(:,:,3) has garbage in it
             !want up_up in b_param_theta_pseudo(:,:,1), up_up=down_down in b_param_theta_pseudo(:,:,2), up_down in b_param_theta_pseudo(:,:,3)
             b_param_theta_pseudo(:,:,3)=b_param_theta_pseudo(:,:,2)
             b_param_theta_pseudo(:,:,2)=b_param_theta_pseudo(:,:,1)
             !don't want to re-read
             read_b_param = .false.
          end if
       elseif( read_b_param_theta_nb == b_param_spin_1_nb) then
          !enough parameters read in for:  up_up = down_down = up_down
          if ( theta_spin_dependence .ne. 1) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for theta_spin_dependence specified', theta_spin_dependence .eq. 1)
          else
             !correct number of parameters so reshape into correct dimensions for b_param_theta_pseudo
             b_param_theta_pseudo = reshape( read_b_param_theta, (/ nb_b_param_theta, nctype, spin_dependencies_nb /), pad=read_b_param_theta)   !up_up = down_down = up_down
             ! up_up is in b_param_theta_pseudo(:,:,1), b_param_theta_pseudo(:,:,2) has garbage in it and b_param_theta_pseudo(:,:,3) has garbage in it
             !want up_up in b_param_theta_pseudo(:,:,1), up_up=down_down in b_param_theta_pseudo(:,:,2), up_up=up_down in b_param_theta_pseudo(:,:,3)
             b_param_theta_pseudo(:,:,2)=b_param_theta_pseudo(:,:,1)
             b_param_theta_pseudo(:,:,3)=b_param_theta_pseudo(:,:,1)
             !don't want to re-read
             read_b_param = .false.
          end if
       else
          !wrong number of parameters
          call die(lhere,'read_b_param_theta has an incorrect number of parameters.')
       end if
    else
       !don't read in parameters
       !start from random guess
       call random_number(b_param_theta_pseudo(:,:,:))
    endif



    !impose cusp
    !To impose the cusp we take the pivot positions (dependent parameters) in our reduced row echelon form of the matrix of cusp constraints
    !and we set this dependent parameter equal to the negative of the dot product of the rest of the row with
    !the corresponding a parameters. This is nothing fancy, but simply a standard way of solving linear equations

    !What is the spin dependence?
    if ( theta_spin_dependence .eq. 3) then
       !loop over depend parameters
       do dep_b_param_i=1, size(dep_b_param_theta_pseudo,1)
          !loop over center types
          do ctype_i=1, nctype
             !loop of spin_dependeny
             do spin_dependency_i=1, theta_spin_dependence
                b_param_theta_pseudo(dep_b_param_theta_pseudo(dep_b_param_i), ctype_i, spin_dependency_i) =&
                     & dot_product( -b_param_theta_pseudo_cond(dep_b_param_i, dep_b_param_theta_pseudo(dep_b_param_i)+1:),&
                     & b_param_theta_pseudo(dep_b_param_theta_pseudo(dep_b_param_i) + 1:, ctype_i, spin_dependency_i))
             end do
          end do
       end do
    elseif ( theta_spin_dependence .eq. 2) then
       !        up_up = down_down != up_down
       spin_dependency_i = 1 !up-up
       !loop over depend parameters
       do dep_b_param_i=1, size(dep_b_param_theta_pseudo,1)
          !loop over center types
          do ctype_i=1, nctype
             b_param_theta_pseudo(dep_b_param_theta_pseudo(dep_b_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -b_param_theta_pseudo_cond(dep_b_param_i, dep_b_param_theta_pseudo(dep_b_param_i)+1:),&
                  & b_param_theta_pseudo(dep_b_param_theta_pseudo(dep_b_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       spin_dependency_i = 3 !up_down or down_up
       !loop over depend parameters
       do dep_b_param_i=1, size(dep_b_param_theta_pseudo,1)
          !loop over center types
          do ctype_i=1, nctype
             b_param_theta_pseudo(dep_b_param_theta_pseudo(dep_b_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -b_param_theta_pseudo_cond(dep_b_param_i, dep_b_param_theta_pseudo(dep_b_param_i)+1:),&
                  & b_param_theta_pseudo(dep_b_param_theta_pseudo(dep_b_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !require up-up = down-down
       b_param_theta_pseudo(:,:,2) = b_param_theta_pseudo(:,:,1)
    elseif( theta_spin_dependence .eq. 1) then
       !        up_up = down_down = up_down
       spin_dependency_i = 1 !up-up
       !loop over depend parameters
       do dep_b_param_i=1, size(dep_b_param_theta_pseudo,1)
          !loop over center types
          do ctype_i=1, nctype
             b_param_theta_pseudo(dep_b_param_theta_pseudo(dep_b_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -b_param_theta_pseudo_cond(dep_b_param_i, dep_b_param_theta_pseudo(dep_b_param_i)+1:),&
                  & b_param_theta_pseudo(dep_b_param_theta_pseudo(dep_b_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !require up-up = down-down=up-down
       b_param_theta_pseudo(:,:,2) = b_param_theta_pseudo(:,:,1)
       b_param_theta_pseudo(:,:,3) = b_param_theta_pseudo(:,:,1)
    end if




  end subroutine b_param_theta_pseudo_bld

  !==============================================================================================




  !==============================================================================================

  subroutine b_param_theta_all_elec_bld

    !---------------------------------------------------------------------------
    ! Description : build b_param_theta_all_elec. This is a rank 3 array.
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
    integer  :: b_param_spin_3_nb, b_param_spin_2_nb, b_param_spin_1_nb, dep_b_param_i, ctype_i, spin_dependency_i
    character(len=max_string_len_rout), save :: lhere = 'b_param_theta_all_elec_bld'
    integer, parameter :: spin_dependencies_nb = 3    !because three choices up_up, down_down, up_down

    ! header
    if (header_exe) then

       call object_create('b_param_theta_all_elec')

       !order of expansion
       call object_needed('order_theta_bf')
       !how many types of atoms
       call object_needed('nctype')
       !what is the spin dependence?
       call object_needed('theta_spin_dependence')
       !have the parameters been read in?
       call object_needed('read_b_param')
       !how many parameters were read in?
       call object_needed('read_b_param_theta_nb')
       !what are the dependent parameters?
       call object_needed('dep_b_param_theta_all_elec')
       !what are the conditions?
       call object_needed('b_param_theta_all_elec_cond')
       !How many a parameters
       call object_needed('nb_b_param_theta')

       return
    endif


    !Note that we need to have the option of using read_b_param_theta if it is read in
    !However, we don't know apriori if it will be read in.
    !Thus we don't call object_needed('read_b_param_theta')
    !To accomplish this correctly we need to use the utility object_provide_in_node
    !We cant simply use object_provide because we are inside the tree.

    if (read_b_param) then
       !a parameters were read in
       call object_provide_in_node(lhere, 'read_b_param_theta')
    end if

    ! allocation
    call object_alloc('b_param_theta_all_elec', b_param_theta_all_elec, nb_b_param_theta, nctype, spin_dependencies_nb)

    !do we read in parameters?
    if (read_b_param) then
       !read in parameters
       b_param_spin_3_nb = (nb_b_param_theta) * nctype * 3         !  up_up != down_down != up_down
       b_param_spin_2_nb = (nb_b_param_theta) * nctype * 2         !  up_up = down_down != up_down
       b_param_spin_1_nb = (nb_b_param_theta) * nctype             !  up_up = down_down = up_down
       if (read_b_param_theta_nb == b_param_spin_3_nb) then
          !enough parameters were read in for:  up_up != down_down != up_down
          if (theta_spin_dependence .ne. 3) then
             !too many parameters were read in
             call require(lhere, 'number of parameters correct for theta_spin_dependence specified', theta_spin_dependence .eq. 3)
          else
             !correct number of parameters so reshape into correct dimensions for b_param_theta
             b_param_theta_all_elec = reshape( read_b_param_theta, (/ nb_b_param_theta, nctype, spin_dependencies_nb /)) !  up_up != down_down != up_down
             !don't want to re-read
             read_b_param = .false.
          end if
       elseif( read_b_param_theta_nb == b_param_spin_2_nb) then
          !enough parameters read in for:  up_up = down_down != up_down
          if ( theta_spin_dependence .ne. 2) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for theta_spin_dependence specified', theta_spin_dependence .eq. 2)
          else
             !correct number of parameters so reshape into correct dimensions for b_param_theta_all_elec
             b_param_theta_all_elec = reshape( read_b_param_theta, (/ nb_b_param_theta, nctype, spin_dependencies_nb /), pad=read_b_param_theta)   !up_up = down_down != up_down
             ! up_up is in b_param_theta_all_elec(:,:,1), up_down is in b_param_theta_all_elec(:,:,2) and b_param_theta_all_elec(:,:,3) has garbage in it
             !want up_up in b_param_theta_all_elec(:,:,1), up_up=down_down in b_param_theta_all_elec(:,:,2), up_down in b_param_theta_all_elec(:,:,3)
             b_param_theta_all_elec(:,:,3)=b_param_theta_all_elec(:,:,2)
             b_param_theta_all_elec(:,:,2)=b_param_theta_all_elec(:,:,1)
             !don't want to re-read
             read_b_param = .false.
          end if
       elseif( read_b_param_theta_nb == b_param_spin_1_nb) then
          !enough parameters read in for:  up_up = down_down = up_down
          if ( theta_spin_dependence .ne. 1) then
             !incorrect number of parameters
             call require(lhere, 'number of parameters correct for theta_spin_dependence specified', theta_spin_dependence .eq. 1)
          else
             !correct number of parameters so reshape into correct dimensions for b_param_theta_all_elec
             b_param_theta_all_elec = reshape( read_b_param_theta, (/ nb_b_param_theta, nctype, spin_dependencies_nb /), pad=read_b_param_theta)   !up_up = down_down = up_down
             ! up_up is in b_param_theta_all_elec(:,:,1), b_param_theta_all_elec(:,:,2) has garbage in it and b_param_theta_all_elec(:,:,3) has garbage in it
             !want up_up in b_param_theta_all_elec(:,:,1), up_up=down_down in b_param_theta_all_elec(:,:,2), up_up=up_down in b_param_theta_all_elec(:,:,3)
             b_param_theta_all_elec(:,:,2)=b_param_theta_all_elec(:,:,1)
             b_param_theta_all_elec(:,:,3)=b_param_theta_all_elec(:,:,1)
             !don't want to re-read
             read_b_param = .false.
          end if
       else
          !wrong number of parameters
          call die(lhere,'read_b_param_theta has an incorrect number of parameters.')
       end if
    else
       !don't read in parameters
       !start from random guess
       call random_number(b_param_theta_all_elec(:,:,:))
    endif



    !impose cusp
    !To impose the cusp we take the pivot positions in our reduced row echelon form of the matrix of cusp constraints
    !and we set this dependent parameter equal to the negative of the dot product of the rest of the row with
    !the corresponding a parameters. This is nothing fancy, but simply a standard way of solving linear equations

    !What is the spin dependence?
    if ( theta_spin_dependence .eq. 3) then
       do dep_b_param_i=1, size(dep_b_param_theta_all_elec,1)
          !loop over center types
          do ctype_i=1, nctype
             !loop over spins
             do spin_dependency_i=1, theta_spin_dependence
                b_param_theta_all_elec(dep_b_param_theta_all_elec(dep_b_param_i), ctype_i, spin_dependency_i) =&
                     & dot_product( -b_param_theta_all_elec_cond(dep_b_param_i, dep_b_param_theta_all_elec(dep_b_param_i)+1:),&
                     & b_param_theta_all_elec(dep_b_param_theta_all_elec(dep_b_param_i) + 1:, ctype_i, spin_dependency_i))
             end do
          end do
       end do
    elseif ( theta_spin_dependence .eq. 2) then
       !        up_up = down_down != up_down
       spin_dependency_i = 1 !up-up
       !loop over dependent parameters
       do dep_b_param_i=1, size(dep_b_param_theta_all_elec,1)
          !loop over center types
          do ctype_i=1, nctype
             b_param_theta_all_elec(dep_b_param_theta_all_elec(dep_b_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -b_param_theta_all_elec_cond(dep_b_param_i, dep_b_param_theta_all_elec(dep_b_param_i)+1:),&
                  & b_param_theta_all_elec(dep_b_param_theta_all_elec(dep_b_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       spin_dependency_i = 3 !up_down or down_up
       !loop over dependent parameters
       do dep_b_param_i=1, size(dep_b_param_theta_all_elec,1)
          !loop over center types
          do ctype_i=1, nctype
             b_param_theta_all_elec(dep_b_param_theta_all_elec(dep_b_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -b_param_theta_all_elec_cond(dep_b_param_i, dep_b_param_theta_all_elec(dep_b_param_i)+1:),&
                  & b_param_theta_all_elec(dep_b_param_theta_all_elec(dep_b_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !require up-up = down-down
       b_param_theta_all_elec(:,:,2) = b_param_theta_all_elec(:,:,1)
    elseif( theta_spin_dependence .eq. 1) then
       !        up_up = down_down = up_down
       spin_dependency_i = 1 !up-up
       !loop over dependent parameters
       do dep_b_param_i=1, size(dep_b_param_theta_all_elec,1)
          !loop over center types
          do ctype_i=1, nctype
             b_param_theta_all_elec(dep_b_param_theta_all_elec(dep_b_param_i), ctype_i, spin_dependency_i) =&
                  & dot_product( -b_param_theta_all_elec_cond(dep_b_param_i, dep_b_param_theta_all_elec(dep_b_param_i)+1:),&
                  & b_param_theta_all_elec(dep_b_param_theta_all_elec(dep_b_param_i)+1:, ctype_i, spin_dependency_i))
          end do
       end do
       !require up-up = down-down=up-down
       b_param_theta_all_elec(:,:,2) = b_param_theta_all_elec(:,:,1)
       b_param_theta_all_elec(:,:,3) = b_param_theta_all_elec(:,:,1)
    end if




  end subroutine b_param_theta_all_elec_bld

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
    integer  :: order_i, order_j, order_k, order_k_max

    ! header
    if (header_exe) then

       call object_create('nb_b_param_theta')

       !need order of expansion
       call object_needed('order_theta_bf')

       return
    endif

    call object_associate('nb_b_param_theta',nb_b_param_theta)

    !determine the number of coefficients for the three body bf
    nb_b_param_theta = 0 !intialize parameter counter
    do order_i = 2, order_theta_bf
       do order_j = 0, order_i - 1
          if (order_j .ne. 0) then
             order_k_max = order_i - order_j
          else
             order_k_max = order_i - 2
          end if
          do order_k = 0, order_k_max
             if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
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

  subroutine dep_b_param_theta_pseudo_bld

    !---------------------------------------------------------------------------
    ! Description : build dep_b_param_theta_pseudo. This is a rank one array.
    !               This produces a list of the dependent parameters from b_param_theta_pseudo_cond.
    !               To do this we find the pivots in a row reduced matrix.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer :: nb_dep_b_param_theta_pseudo
    integer :: row_i, pivot_column_i, column_i


    ! header
    if (header_exe) then

       call object_create('dep_b_param_theta_pseudo')

       !need number of a parameters
       call object_needed('nb_b_param_theta')
       !need parallel spin conditions
       call object_needed('b_param_theta_pseudo_cond')

       return
    endif


    !initialize the number of dependent parameters
    nb_dep_b_param_theta_pseudo = 0
    do row_i = 1, size(b_param_theta_pseudo_cond,1)
       if (sum(abs(b_param_theta_pseudo_cond(row_i,:))) < 1.d-10) then
          !no pivot
          cycle
       end if
       !pivot found, so increment counter
       nb_dep_b_param_theta_pseudo = nb_dep_b_param_theta_pseudo + 1
    end do

    ! allocation
    call object_alloc('dep_b_param_theta_pseudo', dep_b_param_theta_pseudo, nb_dep_b_param_theta_pseudo)


    !intialize the pivot column
    pivot_column_i = 1
    do row_i = 1, nb_dep_b_param_theta_pseudo
       do column_i = pivot_column_i, nb_b_param_theta
          !first non zero element  (should be a 1)
          if (b_param_theta_pseudo_cond(row_i,column_i) > 1.d-10) then
             !pivot found in this row
             dep_b_param_theta_pseudo(row_i) = column_i
             !next pivot column must be to the right of this one
             pivot_column_i=column_i + 1
             !so move on to next row
             exit
          end if
       end do
    end do

  end subroutine dep_b_param_theta_pseudo_bld

  !==============================================================================================




  !==============================================================================================

  subroutine dep_b_param_theta_all_elec_bld

    !---------------------------------------------------------------------------
    ! Description : build dep_b_param_theta_all_elec. This is a rank one array.
    !               This produces a list of the dependent parameters from b_param_theta_all_elec_cond.
    !               To do this we find the pivots in a row reduced matrix.
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------
    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer :: nb_dep_b_param_theta_all_elec
    integer :: row_i, pivot_column_i, column_i

    ! header
    if (header_exe) then

       call object_create('dep_b_param_theta_all_elec')

       !need number of a parameters
       call object_needed('nb_b_param_theta')
       !need anti parallel spin conditions
       call object_needed('b_param_theta_all_elec_cond')

       return
    endif


    !initialize the number of dependent parameters
    nb_dep_b_param_theta_all_elec = 0
    do row_i = 1, size(b_param_theta_all_elec_cond,1)
       if (sum(abs(b_param_theta_all_elec_cond(row_i,:))) < 1.d-10) then
          !no pivot
          cycle
       end if
       !pivot found, so increment counter
       nb_dep_b_param_theta_all_elec = nb_dep_b_param_theta_all_elec + 1
    end do

    ! allocation
    call object_alloc('dep_b_param_theta_all_elec', dep_b_param_theta_all_elec, nb_dep_b_param_theta_all_elec)


    !intialize the pivot column
    pivot_column_i = 1
    do row_i = 1, nb_dep_b_param_theta_all_elec
       do column_i = pivot_column_i, nb_b_param_theta
          !first non zero element  (should be a 1)
          if (b_param_theta_all_elec_cond(row_i,column_i) > 1.d-10) then
             !pivot found in this row
             dep_b_param_theta_all_elec(row_i) = column_i
             !next pivot column must be to the right of this one
             pivot_column_i=column_i + 1
             !so move on to next row
             exit
          end if
       end do
    end do

  end subroutine dep_b_param_theta_all_elec_bld


  !===================================================



  !==============================================================================================

  subroutine b_param_theta_pseudo_cond_bld

    !---------------------------------------------------------------------------
    ! Description : build b_param_theta_pseudo_cond_bld. This is a rank two array.
    !               The first index is the number of conditions ( some are possibly redundant).
    !               The second index is the number a parameters.
    !               We use gaussian elimination to row reduce the conditions and check for redundancies.
    !
    !               cusp condition constraints:   d theta|
    !                                             -----  |        = 0
    !                                             d r_ij |
    !                                                     r_ij=0
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------

    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_i, order_j, order_k, order_k_max, param_i

    ! header
    if (header_exe) then

       call object_create('b_param_theta_pseudo_cond')

       !need order of expansion
       call object_needed('order_theta_bf')
       !need number of a parameters
       call object_needed('nb_b_param_theta')

       return
    endif

    ! allocation
    ! one derivative constraints so number of contraints is (order_theta_bf -1 )
    call object_alloc('b_param_theta_pseudo_cond', b_param_theta_pseudo_cond, order_theta_bf - 1, nb_b_param_theta)

    !initialize
    b_param_theta_pseudo_cond(:,:) = 0.d0

    param_i = 0 !intialize parameter counter
    !setup matrix of constraints
    do order_i = 2, order_theta_bf
       do order_j = 0, order_i - 1
          if (order_j .ne. 0) then
             order_k_max = order_i - order_j
          else
             order_k_max = order_i - 2
          end if
          do order_k = 0, order_k_max
             if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                cycle
             else
                param_i = param_i + 1
                !=================================
                ! d theta|
                ! -----  |        = 0
                ! d r_ij |
                !         r_ij=0
                !=================================
                b_param_theta_pseudo_cond(order_i-order_j, param_i) = b_param_theta_pseudo_cond(order_i-order_j, param_i) + order_j
             endif
          enddo
       enddo
    enddo

    !use gaussian elimintation on constraints to remove redundancies
    call row_reduce(b_param_theta_pseudo_cond)

  end subroutine b_param_theta_pseudo_cond_bld

  !==============================================================================================



  !==============================================================================================

  subroutine b_param_theta_all_elec_cond_bld

    !---------------------------------------------------------------------------
    ! Description : build b_param_theta_all_elec_cond_bld. This is a rank two array.
    !               The first index is the number of condition ( some are possibly redundant).
    !               The second index is the number a parameters.
    !               We use gaussian elimination to row reduce the conditions and check for redundancies.
    !
    !               cusp condition constraints:   d theta|
    !                                             -----  |        = 0
    !                                             d r_ij |
    !                                                     r_ij=0
    !
    !                                              theta|
    !                                                   |        = 0
    !                                                   |
    !                                                     r_iI=0
    !
    !
    !
    !
    !
    ! Created     : F. Petruzielo, 14 Aug 2008
    !---------------------------------------------------------------------------

    implicit none

    !need nctype
    include 'commons.h'

    !local
    integer  :: order_i, order_j, order_k, order_k_max, param_i

    ! header
    if (header_exe) then

       call object_create('b_param_theta_all_elec_cond')

       !need order of expansion
       call object_needed('order_theta_bf')
       !need number of a parameters
       call object_needed('nb_b_param_theta')

       return
    endif

    ! allocation
    ! one derivative constraints evaluated at 0, one constraint on function so number of contraints is 2 * order_theta_bf - 1
    call object_alloc('b_param_theta_all_elec_cond', b_param_theta_all_elec_cond, 2 * order_theta_bf - 1, nb_b_param_theta)

    !initialize
    b_param_theta_all_elec_cond(:,:) = 0.d0

    param_i = 0 !intialize parameter counter
    !setup matrix of constraints
    do order_i = 2, order_theta_bf
       do order_j = 0, order_i - 1
          if (order_j .ne. 0) then
             order_k_max = order_i - order_j
          else
             order_k_max = order_i - 2
          end if
          do order_k = 0, order_k_max
             if (modulo(order_i - order_j - order_k, 2) .ne. 0) then
                cycle
             else
                param_i = param_i + 1
                !=================================
                ! d theta|
                ! -----  |        = 0
                ! d r_ij |
                !         r_ij=0
                !=================================
                b_param_theta_all_elec_cond(order_i-order_j, param_i) = b_param_theta_all_elec_cond(order_i-order_j, param_i) + order_j
                !================================================
                ! theta|
                !      |        = 0
                !      |
                !       r_iI=0
                !===============================================
                !note that the row needs to be shifted down by (order_theta_bf - 1) because those rows hold the above constraints
                b_param_theta_all_elec_cond((order_theta_bf - 1) + (order_i+order_j-order_k)/2, param_i) = b_param_theta_all_elec_cond((order_theta_bf -1 ) + (order_i+order_j-order_k)/2, param_i) + 1.d0
                b_param_theta_all_elec_cond((order_theta_bf - 1) + (order_i+order_j+order_k)/2, param_i) = b_param_theta_all_elec_cond((order_theta_bf -1 ) + (order_i+order_j+order_k)/2, param_i) + 1.d0
             endif
          enddo
       enddo
    enddo

    !use gaussian elimintation on constraints to remove redundancies
    call row_reduce(b_param_theta_all_elec_cond)

  end subroutine b_param_theta_all_elec_cond_bld

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
    real(dp),dimension(:,:),intent(INOUT) :: matrix_coeff

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

End module backflow_mod
