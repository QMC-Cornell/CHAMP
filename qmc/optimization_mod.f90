module optimization_mod

  use all_tools_mod
  use opt_nwt_mod
  use opt_lin_mod
  use opt_ptb_mod
  use opt_ovlp_fn_mod, only: delta_ovlp_fn, delta_ovlp_fn_linear, opt_ovlp_fn_menu, l_opt_ovlp_fn_nonsemiortho, l_opt_ovlp_fn_linear, ovlp_trial_fn, ovlp_trial_fn_sav, ovlp_trial_fn_over_ovlp_trial, ovlp_trial_fn_over_ovlp_trial_sav, l_opt_ovlp_branching, gradient_ovlp, gradient_ovlp_norm, delta_ovlp_fn_nonsemiortho
  use opt_common_mod
  use nuclei_mod
  use orbitals_mod
  use periodic_jastrow_mod
  use deriv_mod
  use montecarlo_mod
  use control_mod
  use vmc_mod
  use dmc_mod
  use cusp_mod
  use optim_mod
  use dets_mod

! Declaration of global variables and default values
  character(len=max_string_len)   :: opt_method = 'linear'

  integer                 :: parameter_type_nb = 0
# if defined (PATHSCALE)
   character(len=max_string_len) :: parameter_type (max_string_array_len) ! for pathscale compiler
# else
   character(len=max_string_len), allocatable :: parameter_type (:)
# endif

  integer                 :: iter_opt_min_nb = 0
  logical                 :: l_check_convergence = .true.
  integer                 :: check_convergence_nb = 1

  logical                 :: l_increase_accuracy = .true.
  logical                 :: l_decrease_error      = .true.
  logical                 :: l_decrease_error_adaptative  = .false.
  real(dp)                :: decrease_error_factor = 2.d0
  real(dp)                :: decrease_error_limit  = 0.d0
  logical                 :: l_increase_blocks     = .false.
  real(dp)                :: increase_blocks_factor = 2.d0
  logical                 :: l_decrease_p_var = .false.
  logical                 :: l_ortho_orb_vir_to_orb_occ = .false.
  logical                 :: l_approx_orb_rot = .false.
  logical                 :: l_reweight = .false.
!  logical                 :: l_reset_walker_weights_sum_block = .false.
  integer                 :: reweight_power = 1
  real(dp)                :: reweight_scale = 10.d0

  real(dp)                :: energy_threshold          = 1.d-3
  real(dp)                :: gradient_norm_threshold   = 0.1d0

  real(dp)                :: delta_param_norm
  real(dp)                :: delta_csf_norm
  real(dp)                :: delta_jas_norm
  real(dp)                :: delta_param_norm_max_3d = 10.d0
  real(dp)                :: delta_param_norm_max_2d = 500.d0
  real(dp), allocatable   :: delta_param (:)
  real(dp), allocatable   :: delta_csf (:)
  real(dp), allocatable   :: delta_csf_rot (:)
  real(dp), allocatable   :: delta_jas (:)
  real(dp), allocatable   :: delta_pjas (:)
  real(dp), allocatable   :: delta_coef_ex (:)
  real(dp), allocatable   :: delta_mat_rot_1st_order (:,:)
  real(dp), allocatable   :: delta_mat_rot_real (:,:)
  real(dp), allocatable   :: delta_exp (:)
  real(dp), allocatable   :: delta_geo (:)
  real(dp), allocatable   :: delta_c_rp (:,:)
  real(dp), allocatable   :: delta_c_rm (:,:)
  real(dp), allocatable   :: delta_c_ip (:,:)
  real(dp), allocatable   :: delta_c_im (:,:)
  real(dp)                :: add_diag_max  = 1.d10
  logical                 :: l_reset_add_diag = .true.
  real(dp)                :: add_diag_reset_value  = 1.d-8
  logical                 :: do_add_diag_mult_exp = .false.
  logical                 :: l_last_run = .true.

  contains

!===========================================================================
  subroutine optimization_menu
!---------------------------------------------------------------------------
! Description : menu for opt
!
! Created     : J. Toulouse, 12 Oct 2005
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'optimization_menu'
  integer nterms4, param_type_i, param_i, na1, na2, ia, isp, it
  logical l_launch_opt

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of optimization menu ---------------------------------------------------------------------------'

! initialization
  l_launch_opt = .true.
  l_opt = .true.
  l_stab = .true.
  target_state = 0
  l_opt_jas = .false.
  l_opt_pjas = .false.
  l_opt_pjasen = .false.
  l_opt_pjasee = .false.
  l_opt_csf = .false.
  l_opt_orb = .false.
  l_opt_exp = .false.
  l_opt_geo = .false.
  l_approx_orb_rot = .false.
  l_active_orb_ortho_constraint = .true.

! temporary error message
  if (nopt_iter <= 0 .and. .not. use_parser) then
   call die (lhere, 'nopt_iter should be > 0 for optimization.')
  endif

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)') 'HELP for optimization menu:'
   write(6,'(a)') 'optimization'
   write(6,'(a)') ' optimize = [logical] : launch wave function optimization? (default=true)'
   write(6,'(a)') ' parameters jastrow csfs orbitals exponents geometry end : list of parameter types to optimized'
   write(6,'(a)') ' method                    : energy optimization method'
   write(6,'(a)') '        = linear (default) : linear energy optimization method'
   write(6,'(a)') '        = newton           : Newton energy optimization method'
   write(6,'(a)') '        = perturbative     : Perturbative energy optimization method'
   write(6,'(a)') '        = overlap_fn       : Maximize overlap with fixed-nodes wavefunction method'
   write(6,'(a)') '        = gradient         : compute only the energy gradient'
   write(6,'(a)') ' newton ... end : menu for Newton optimization method'
   write(6,'(a)') ' linear ... end : menu for linear optimization method'
   write(6,'(a)') ' perturbative ... end : menu for perturbative optimization method'
   write(6,'(a)') ' overlap_fn ... end : menu for overlap optimization method'
   write(6,'(a)') ' p_var = [real] : fraction of variance to minimize (default=0)'
   write(6,'(a)') ' stabilize  = [logical] stabilize the minimization? (default=true)'
   write(6,'(a)') ' stabilization : choice of stabilization of the minimization'
   write(6,'(a)') '              = identity: add multiple of identity matrix to Hessian or Hamiltonian (default)'
   write(6,'(a)') '              = overlap: add multiple of overlap matrix to Hamiltonian (only for linear method)'
   write(6,'(a)') ' add_diag     = [real] : stabilization constant (default=1.d-8)'
   write(6,'(a)') ' add_diag_max = [real] : maximum allowed value of add_diag (default=1.d10)'
   write(6,'(a)') ' reset_add_diag = [bool] : reset add_diag to add_diag_reset_value at each step before adjustment (default=true)'
   write(6,'(a)') ' add_diag_reset_value = [real] : value to which add_diag will be reset to at each step (default=1.d-12)'
   write(6,'(a)') ' iter_opt_min_nb = [integer] : minimum number of optimization iterations (default=0)'
   write(6,'(a)') ' iter_opt_max_nb = [integer] : maximum number of optimization iterations (default=100)'
   write(6,'(a)') ' last_run = [logical] : perform a last run with the last predicted parameters? (default=true)'
   write(6,'(a)') ' increase_accuracy = [logical] : default=true, increase statistical accuracy at each step?'
   write(6,'(a)') ' decrease_error = [logical] : default=true, decrease statistical error at each step?'
   write(6,'(a)') ' decrease_error_adaptative = [logical] : default=true, decrease statistical error according to energy difference'
   write(6,'(a)') ' decrease_error_factor = [real] : default=2, decrease statistical error at each step by this factor'
   write(6,'(a)') ' decrease_error_limit = [real] : default=energy_threshold/2, decrease statistical error at each step until this limit'
   write(6,'(a)') ' increase_blocks = [logical] : default=false, increase number of blocks at each step?'
   write(6,'(a)') ' increase_blocks_factor = [real] : default=2, increase number of blocks at each step by this factor'
   write(6,'(a)') ' increase_blocks_limit = [integer] : default=10000, increase number of blocks at each step until this limit'
   write(6,'(a)') ' check_convergence = [logical] : default=true, check and stop if energy difference is less than tolerance?'
   write(6,'(a)') ' check_convergence_nb = [integer] : default=1, number of consecutive steps energy difference must be less than tolerance before stopping'
   write(6,'(a)') ' energy_threshold = [real] : threshold on energy for stopping criterium (default=1.d-3)'
   write(6,'(a)') ' casscf = [logical] : default=false, is it a CASSCF wave function? If true, help the program by removing redundant active-active excitations for orbital optimization'
   write(6,'(a)') ' check_redundant_orbital_derivative = [logical] : default=true, check for additional redundancies in orbital derivatives'
   write(6,'(a)') ' approx_orb_rot = [logical] : approximate first-order orbital rotation (default=false)'
   write(6,'(a)') ' active_orb_ortho_constraint = [logical] : impose active orbitals orthogonality constraint during optimization (default=true)'
   write(6,'(a)') ' deriv2nd = [logical] : default=true, compute second-order wave function derivatives if necessary?'
   write(6,'(a)') ' delta_param_norm_max_3d = [real] : maximum parameter variation norm allowed for 3-dim systems (default=10)'
   write(6,'(a)') ' delta_param_norm_max_2d = [real] : maximum parameter variation norm allowed for 2,1-dim systems (default=10)'
   write(6,'(a)') ' hessian_variance = [linear|levenberg_marquardt|levenberg_marquardt_cov] : choice of variance hessian (default=linear)'
   write(6,'(a)') ' decrease_p_var= [bool] : decrease progressively proportion of variance (default=false)'
   write(6,'(a)') ' print_orbital_excitations = [bool] print orbital excitation information? (default=false)'
   write(6,'(a)') ' orthonormalize_orbitals = [bool] orthonormalize orbitals at each optimization step? (default=false)'
   write(6,'(a)') ' ortho_orb_vir_to_orb_occ = [bool] : orthogonalize virtual orbitals to occupied orbitals (default=false)'
   write(6,'(a)') ' do_add_diag_mult_exp = [bool] : allow the basis exponents to change only by 20% by increasing add_diag_mult_exp if change > 5% (default=false)'
   write(6,'(a)') ' exp_opt_restrict = [bool] : restriction on exponent parameters to optimize according to basis function types? (default=true)'
   write(6,'(a)') ' deriv_bound = [bool] : applying a bound on the log derivatives of the wave function wrt parameters? (default=false)'
   write(6,'(a)') ' deriv_bound_value = [real] : value of the bound mentioned above (default=10.d0)'
   write(6,'(a)') ' reweight = [bool] : reweight expectation values in optimization with weights depending on distance to node (default=false)'
   write(6,'(a)') ' reweight_power = [integer] : value of power in reweighting expression (default=1)'
   write(6,'(a)') ' reweight_scale = [real] : value of scaling factor in reweighting expression (default=10.d0)'
   write(6,'(a)') ' csfs ... end: list of csfs to optimize (default=all csfs)'
   write(6,'(a)') ' slater_mat_ex_trans_inv_sm = [bool] Sherman-Morison for inverses of single-excited determinants? It shoud be set to true for large electron numbers (default=true)'
   write(6,'(a)') 'end'

  case ('optimize')
   call get_next_value (l_launch_opt)

  case ('parameters')
# if defined (PATHSCALE)
   call get_next_value_list_string ('parameter_type', parameter_type, parameter_type_nb) ! for pathscale compiler
# else
   call get_next_value_list ('parameter_type', parameter_type, parameter_type_nb)
# endif

  case ('method')
   call get_next_value (opt_method)

  case ('newton')
   call opt_nwt_menu

  case ('linear')
   call opt_lin_menu

  case ('perturbative')
   call opt_ptb_menu

  case ('overlap_fn')
   call opt_ovlp_fn_menu

  case ('p_var')
   call get_next_value (p_var)
   call object_modified ('p_var')  
   call require (lhere, 'p_var >= 0 and p_var =< 1', p_var >= 0.d0 .and. p_var <= 1.d0)

  case ('iter_opt_min_nb')
   call get_next_value (iter_opt_min_nb)

  case ('iter_opt_max_nb')
   call get_next_value (iter_opt_max_nb)
   nopt_iter = iter_opt_max_nb

  case ('last_run')
   call get_next_value (l_last_run)

  case ('stabilize')
   call get_next_value (l_stab)

  case ('add_diag')
   call get_next_value (add_diag(1))
   diag_stab = add_diag(1)
   call object_modified ('diag_stab')

  case ('reset_add_diag')
   call get_next_value (l_reset_add_diag)

  case ('add_diag_reset_value')
   call get_next_value (add_diag_reset_value)

  case ('stabilization')
   call get_next_value (stabilization)

  case ('add_diag_max')
   call get_next_value (add_diag_max)

  case ('increase_accuracy')
   call get_next_value (l_increase_accuracy)

  case ('decrease_error')
   call get_next_value (l_decrease_error)
   if(l_decrease_error) then
    l_increase_blocks = .false.
   endif

  case ('decrease_error_adaptative')
   call get_next_value (l_decrease_error_adaptative)

  case ('decrease_error_factor')
   call get_next_value (decrease_error_factor)
   call require (lhere, 'decrease_error_factor > 0', decrease_error_factor > 0) !fp

  case ('decrease_error_limit')
   call get_next_value (decrease_error_limit)
   call require (lhere, 'decrease_error_limit > 0', decrease_error_limit > 0) !fp

  case ('increase_blocks')
   call get_next_value (l_increase_blocks)
   if (l_increase_blocks) then
     l_decrease_error = .false.
   endif

  case ('increase_blocks_factor')
   call get_next_value (increase_blocks_factor)
   call require (lhere, 'increase_blocks_factor > 0', increase_blocks_factor > 0) !fp

  case ('increase_blocks_limit')
   call get_next_value (increase_blocks_limit)
   call require (lhere, 'increase_blocks_limit > 0', increase_blocks_limit > 0) !fp
   nblk_max = nint(increase_blocks_limit)

  case ('check_convergence')
   call get_next_value (l_check_convergence)

  case ('check_convergence_nb')
   call get_next_value (check_convergence_nb)
   call require (lhere, 'check_convergence_nb > 0', check_convergence_nb > 0) !fp

  case ('energy_threshold')
   call get_next_value (energy_threshold)
   call require (lhere, 'energy_threshold > 0', energy_threshold > 0) !fp
   tol_energy = energy_threshold
   call object_modified ('energy_threshold')

  case ('casscf')
   call get_next_value (l_casscf)

  case ('check_redundant_orbital_derivative')
   call get_next_value (l_check_redundant_orbital_derivative)

  case ('approx_orb_rot')
   call get_next_value (l_approx_orb_rot)

  case ('active_orb_ortho_constraint')
   call get_next_value (l_active_orb_ortho_constraint)

  case ('deriv2nd')
   call get_next_value (l_deriv2nd)

  case ('delta_param_norm_max_3d')
   call get_next_value (delta_param_norm_max_3d)

  case ('delta_param_norm_max_2d')
   call get_next_value (delta_param_norm_max_2d)

  case ('hessian_variance')
   call get_next_value (hessian_variance_type)

  case ('decrease_p_var')
   call get_next_value (l_decrease_p_var)

  case ('print_orbital_excitations')
   call get_next_value (l_print_orbital_excitations)

  case ('orthonormalize_orbitals')
   call get_next_value (l_ortho_orb_opt)

  case ('ortho_orb_vir_to_orb_occ')
   call get_next_value (l_ortho_orb_vir_to_orb_occ)

  case ('do_add_diag_mult_exp')
   call get_next_value (do_add_diag_mult_exp)

  case ('exp_opt_restrict')
   call get_next_value (l_exp_opt_restrict)

  case ('deriv_bound')
   call get_next_value (l_deriv_bound)

  case ('deriv_bound_value')
   call get_next_value (deriv_bound_value)
   call require (lhere, 'deriv_bound_value > 1', deriv_bound_value > 1) !fp

  case ('reweight')
   call get_next_value (l_reweight)

  case ('reweight_power')
   call get_next_value (reweight_power)

  case ('reweight_scale')
   call get_next_value (reweight_scale)

  case ('csfs')
   call get_next_value_list_object ('iwcsf', iwcsf, nparmcsf)
   call object_modified ('nparmcsf')

  case ('slater_mat_ex_trans_inv_sm')
   call get_next_value (l_slater_mat_ex_trans_inv_sm)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown keyword >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines

  write(6,'(a,i4)') ' maximum number of iterations = ', iter_opt_max_nb
  write(6,'(a,1pd9.1)') ' stabilization constant: add_diag=',add_diag(1)
  write(6,'(a,f12.4)') ' fraction of variance: p_var=',p_var
  write(6,'(a,es12.4)') ' energy threshold for convergence =',energy_threshold

! if(iter_opt_max_nb /= 0) igradhess=1 ! possibly needed to do optimization with new input but using Cyrus' optimization routines
  igradhess=0 ! turn off Cyrus' optimization routines
  if(iter_opt_max_nb /= 0 .and. nforce > 1) stop 'nforce > 1 not allowed in optimization. At present can optim 1 wf only'
  if(iter_opt_max_nb /= 0 .and. (nwf.lt.3)) stop 'for optimization nwf should be >=3'

! set some default values
  if (decrease_error_limit == 0.d0 ) then
    decrease_error_limit = energy_threshold/2.d0
  endif

  if (.not. l_increase_accuracy) then
    l_decrease_error = .false.
    l_increase_blocks = .false.
  endif

  add_diag_mult_exp = 1.d0
  call object_modified ('add_diag_mult_exp')

  if (.not. l_active_orb_ortho_constraint) then
    l_approx_orb_rot = .true.
    write(6,'(a)') " Warning: switch on approximate first-order orbital rotation because it is the only one available for option active_orb_ortho_constraint = false"
  endif

! choice of optimization method
  select case(trim(opt_method))
   case ('linear')
    l_opt_lin = .true.
   case ('newton')
    l_opt_nwt = .true.
   case ('perturbative')
    l_opt_ptb = .true.
   case ('overlap_fn')
    l_opt_ovlp_fn = .true.
    call require (lhere, 'l_mode_dmc = true', l_mode_dmc)
   case ('gradient')
    l_opt_grad = .true.
    l_last_run = .false.
   case default
    call die (lhere, 'unknown optimization method >'+trim(opt_method)+'<.')
  end select

! parameters to optimize
  write(6,'(a,10a10)') ' Requested parameter types: ',parameter_type(:)
  do param_type_i = 1, parameter_type_nb
   select case(trim(parameter_type(param_type_i)))
   case ('jastrow')
    l_opt_jas = .true.
   case ('pjasen')
    l_opt_pjasen = .true.
    l_opt_pjas = .true.
   case ('pjasee')
    l_opt_pjasee = .true.
    l_opt_pjas = .true.
   case ('csfs')
    l_opt_csf = .true.
   case ('orbitals')
    l_opt_orb = .true.
   case ('exponents')
    l_opt_exp = .true.
   case ('geometry')
    l_opt_geo = .true.
   case default
    call die (lhere, 'unknown parameter type >'+trim(parameter_type(param_type_i))+'<.')
   end select
  enddo

  if (use_parser) then

  call object_provide ('ncsf')
  call object_provide ('nparmcsf')
  call object_provide ('iwcsf')
  if(nparmcsf == ncsf) then
     write(6,'(a,i5,a,i5)') ' Warning: since normalization of wavefn. is arb. nparmcsf=',nparmcsf,' should be <= ncsf-1=',ncsf-1
  endif
  if (l_opt_ovlp_fn_linear .and. l_opt_csf) then
    if (iwcsf(1) /= 1) then
     write(6,*) " Warning: first csf not optimized in ovlp_fn_linear method?"
!     call die (lhere, 'with ovlp_fn_linear method the first csf must be optimized.')
    endif
  endif


! default jastrow parameters to optimize
! For the e-n parameters we are assuming that a(1) and a(2) are not optimized, which is often not true for all-electron calculations.
  if(ijas.le.3) then
   na1=nspin1
   na2=nspin2
  else
   na1=1
   na2=nctype
  endif
  nparmot=0
  call alloc ('nparma', nparma, na2-na1+1)
  do ia=na1,na2
   if(norda==0) then
    nparma(ia)=0
   elseif(norda>=1) then
    nparma(ia)=norda-1
   else
    call die (lhere, 'norda must be >= 0 norda='+norda)
   endif
  enddo
  call alloc ('nparmb', nparmb, nspin2b-nspin1+1)
  do isp=nspin1,nspin2b
   if(nordb==0) then
    nparmb(isp)=0
   elseif(nordb>=1) then
    nparmb(isp)=nordb
   else
    call die (lhere, 'nordb must be >= 0 nordb='+nordb)
   endif
  enddo
  call alloc ('nparmc', nparmc, nctype)
  do it=1,nctype
   if(nordc==0) then
    nparmc(it)=0
   elseif(nordc>=1) then
    nparmc(it)=15
    nparmc(it)=nterms4(nordc)-2*(nordc-1)
   else
    call die (lhere, 'nordc must be >= 0 nordc='+nordc)
   endif
  enddo
  if(ijas.ge.4.and.ijas.le.6) then
    do it=1,nctype
          if(nloc.eq.0) then
! All-electron with analytic slater basis
            if((norda.eq.0.and.nparma(it).gt.0).or.(norda.gt.0 .and. nparma(it).gt.norda+1)) then
              write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
              stop 'nparma too large for norda in all-electron calculation'
            endif
           else
! Pseudopotential with numerical basis (cannot vary a(1) or a(2)
            if(norda.eq.1) stop 'makes no sense to have norda=1 for nloc!=0, i.e. for psp. atoms or Je spheres.'
            if((norda.eq.0.and.nparma(it).gt.0).or.(norda.gt.0 .and. nparma(it).gt.norda-1)) then
              write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
              stop 'nparma too large for norda in pseudopot calculation'
            endif
          endif
          if(isc.le.10 .and.((nordc.le.2.and.nparmc(it).gt.0)           &
          .or.(nordc.eq.3.and.nparmc(it).gt.2).or.(nordc.eq.4.and.nparmc(it).gt.7)  &
          .or.(nordc.eq.5.and.nparmc(it).gt.15).or.(nordc.eq.6.and.nparmc(it).gt.27)&
          .or.(nordc.eq.7.and.nparmc(it).gt.43))) then
            write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
            stop 'nparmc too large for nordc in J_een with cusp conds'
          endif
          if(isc.gt.10 .and.((nordc.le.1.and.nparmc(it).gt.0).or.(nordc.eq.2.and.nparmc(it).gt.2) &
          .or.(nordc.eq.3.and.nparmc(it).gt.6).or.(nordc.eq.4.and.nparmc(it).gt.13)               &
          .or.(nordc.eq.5.and.nparmc(it).gt.23).or.(nordc.eq.6.and.nparmc(it).gt.37)              &
          .or.(nordc.eq.7.and.nparmc(it).gt.55))) then
            write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
            stop 'nparmc too large for nordc without cusp conds'
          endif
     enddo
! For the b coefs. we assume that b(1) is fixed by the cusp-cond.
        do isp=1,nspin1,nspin2b
            if((nordb.eq.0.and.nparmb(isp).gt.0).or.(nordb.gt.0 .and. nparmb(isp).gt.nordb)) then
              write(6,'(''isp,nordb,nparmb(isp)'',3i5)') isp,nordb,nparmb(isp)
              stop 'nparmb too large for nordb'
            endif
        enddo
      endif
! compute nparmj
      nparmj=0
      call alloc ('npoint', npoint, nctype)
      call alloc ('npointa', npointa, na2)
      npointa(1)=0
      do ia=na1,na2
        if(ia.gt.1) npointa(ia)=npointa(ia-1)+nparma(ia-1)
        nparmj=nparmj+nparma(ia)
      enddo
      do isp=nspin1,nspin2b
        nparmj=nparmj+nparmb(isp)
      enddo
      npoint(1)=nparmj
      do it=1,nctype
        if(it.gt.1) npoint(it)=npoint(it-1)+nparmc(it-1)
        nparmj=nparmj+nparmc(it)
      enddo
      nparmjs=nparmj+nparms
      call object_modified ('nparmjs')

      if(ijas.ge.4.and.ijas.le.6) then
!       call alloc ('iwjasa', iwjasa, nparmj, nctype)
        call alloc ('iwjasa', iwjasa, max(maxval(nparma),1), nctype)
        do it=1,nctype 
!         iwjasa (1:nparma(it),it) = (/ 3, 4, 5, 6/)
          do param_i = 1, nparma(it)
            iwjasa(param_i,it) = param_i + 2
          enddo
        enddo
!       call alloc ('iwjasb', iwjasb, nparmj, nspin2b-nspin1+1)
        call alloc ('iwjasb', iwjasb, max(maxval(nparmb),1), nspin2b-nspin1+1)
        do isp=nspin1,nspin2b
!         iwjasb(1:nparmb(isp),isp) = (/2, 3, 4, 5, 6/)
          do param_i = 1, nparmb(isp)
            iwjasb(param_i,isp) = param_i + 1
          enddo
        enddo
!       call alloc ('iwjasc', iwjasc, nparmj, nctype)
        call alloc ('iwjasc', iwjasc, max(maxval(nparmc),1), nctype)
        if(nordc==3) then
          do  it=1,nctype
            iwjasc(1:nparmc(it),it) = (/ 3,   5/)
          enddo
        elseif(nordc==4) then
          do  it=1,nctype
            iwjasc(1:nparmc(it),it) = (/ 3,   5,   7, 8, 9,    11,    13/)
          enddo
        elseif(nordc==5) then
          do  it=1,nctype
            iwjasc(1:nparmc(it),it) = (/ 3,   5,   7, 8, 9,    11,    13, 14, 15, 16, 17, 18,    20, 21,    23/)
          enddo
        elseif(nordc==6) then
          do  it=1,nctype
            iwjasc(1:nparmc(it),it) = (/ 3,   5,   7, 8, 9,    11,    13, 14, 15, 16, 17, 18,    20, 21,    23, 24, 25, 26, 27, 28, 29, 30, 31,    33, 34,    36, 37/)
          enddo
        elseif(nordc==7) then
          do  it=1,nctype
            iwjasc(1:nparmc(it),it) = (/ 3,   5,   7, 8, 9,    11,    13, 14, 15, 16, 17, 18,    20, 21,    23, 24, 25, 26, 27, 28, 29, 30, 31,    33, 34,    36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,    50, 51, 52, 54, 55/)
          enddo
        elseif(nordc>=8) then
          call die (lhere, 'iwjasc is not implemented for nordc >= 8 nordc='+nordc)
        endif
      endif

      if(icusp2.ge.1 .and. ijas.eq.3 .and. isc.le.7) call cuspinit3(1)
      if(icusp2.ge.1 .and. ijas.eq.4 .and. isc.le.10) call cuspinit4(0)
      call object_modified ('nparma')
      call object_modified ('nparmb')
      call object_modified ('nparmc')
      call object_modified ('iwjasa')
      call object_modified ('iwjasb')
      call object_modified ('iwjasc')

  endif ! if use_parser

! set numbers parameters to zero if not optimized
!JT  nparmj = nparmj_input
!JT  nparmcsf = nparmcsf_input
  call object_modified ('nparmj')
  call object_modified ('nparmcsf')
  if (.not. l_opt_jas) then
    nparmj=0
    call object_modified ('nparmj')
  endif
  if (.not. l_opt_pjas) then
     param_pjas_nb  = 0
     call object_modified ('param_pjas_nb')
  endif
  if (.not. l_opt_pjasen) then
     param_pjasen_nb  = 0
     call object_modified ('param_pjasen_nb')
  endif
  if (.not. l_opt_pjasee) then
     param_pjasee_nb  = 0
     call object_modified ('param_pjasee_nb')
  endif
  if (.not. l_opt_csf) then
    nparmcsf=0
    call object_modified ('nparmcsf')
  endif
  if (.not. l_opt_orb) then
    param_orb_nb  = 0
    call object_modified ('param_orb_nb')
  endif
  if (.not. l_opt_exp) then
    param_exp_nb  = 0
    call object_modified ('param_exp_nb')
  endif
  if (.not. l_opt_geo) then
    param_geo_nb  = 0
    call object_modified ('param_geo_nb')
  endif
  
  if (ibasis.le.3) nparmd=nparmcsf
  call object_modified ('nparmd')

! check consistency of options
  if (l_opt_ptb) then
   if (l_opt_jas .or. l_opt_csf) then
    if (l_opt_orb_eig) then
      call die (lhere, 'perturbative method with option use_orbital_eigenvalues=true is allowed when optimizing orbitals only.')
    endif
    if (l_diagonal_overlap) then
      call die (lhere, 'perturbative method with option overlap_diagonal=true is allowed when optimizing orbitals only.')
    endif
   endif
  endif

  if (l_opt_nwt .and. l_opt_pjas) then
     call  die (lhere, 'Optimization of periodic Jastrow parameters is done with linear method only')
  endif

  if (l_opt_ovlp_fn_linear .and. (l_opt_jas .or. l_opt_exp)) then
     call  die (lhere, 'Overlap-fn optimization method with linear scaling is only for CSF and orbital parameters')
  endif

! Warnings for pertubative method
  if (l_opt_ptb) then
    if(l_opt_jas) then
      write(6,'(a)') 'Warning: the perturbative method is usually very bad for optimizing the Jastrow parameters.'
      l_warning = .true.
    endif
    if (l_opt_csf) then
      write(6,'(a)') 'Warning: the perturbative method is usally not very good for optimizing of the CSF parameters.'
      l_warning = .true.
    endif
  endif

! compute or not energies associated to orbitals derivatives
  l_opt_orb_energy = l_opt_orb .and. .not. l_opt_orb_eig

! compute or not 2nd derivatives of Jastrow
  if (l_opt_jas .and. l_opt_nwt) then
    l_opt_jas_2nd_deriv = .true.
  endif

! Orbital optimization
  if (l_opt_orb) then
!  print information for orbital optimization
   write(6,*)
   write(6,'(3a)') ' Orbital optimization information:'
   call object_provide ('param_orb_nb')
   call object_provide ('det_ex_unq_up_nb')

   call object_provide ('orb_opt_last_lab')
   norb = orb_opt_last_lab
   write(6,'(a,i8)') ' Number of computed orbitals will be ', norb
  endif

! Exponent optimization
  if (l_opt_exp) then
!  print information for exponent optimization
   write(6,*)
   write(6,'(3a)') ' Exponent optimization information:'
   call object_provide ('param_exp_nb')

!  compute orbital coefficients on orthonormalized basis functions
   if (trim(basis_functions_varied) == 'orthonormalized') then
    call coef_orb_on_ortho_basis_from_coef (1)
   endif
  endif ! l_opt_exp

! Print number of parameters to optimized
  call object_provide ('nparm')
  call object_provide ('nparmcsf')
  call object_provide ('param_pjas_nb')
  call object_provide ('param_orb_nb')
  call object_provide ('param_exp_nb')
  call object_provide ('param_geo_nb')
  call object_provide ('param_nb')
  write(6,*)
  write(6,'(a,i5)') ' Number of Jastrow parameters:   ', nparmj
  write(6,'(a,i5)') ' Number of periodic Jastrow parameters: ', param_pjas_nb
  write(6,'(a,i5)') ' Number of CSF parameters:       ', nparmcsf
  write(6,'(a,i5)') ' Number of orbital parameters:   ', param_orb_nb
  write(6,'(a,i5)') ' Number of exponent parameters:  ', param_exp_nb
  write(6,'(a,i5)') ' Number of geometry parameters:  ', param_geo_nb
  write(6,'(a,i5)') ' Total number of parameters:     ', param_nb
  write(6,*)

! set nparm if new input
  if (use_parser) then
   nparm = nparmj+nparmcsf
   call object_modified ('nparm')
  endif

  write(6,'(a)') 'End of optimization menu ---------------------------------------------------------------------------------'

! launch optimization
  if (l_launch_opt) then
   call optimization
   run_done = .true.
  endif

  end subroutine optimization_menu

!===========================================================================
  subroutine optimization
!---------------------------------------------------------------------------
! Description : routine for wave function optimization
!
! Created     : J. Toulouse, 18 Feb 2006
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'optimization'
  integer iter, parm_i, parm_j, param_i, iter_best
  real(dp) eloc_av_previous
  real(dp) d_eloc_av
  logical  l_convergence_reached
  integer convergence_reached_nb
  integer :: move_rejected = 0
  real(dp) energy_plus_err, energy_plus_err_best, energy_best

! begin
  write(6,*)
  write(6,'(a)') '************************* WAVE FUNCTION OPTIMIZATION *************************'

! Initializations
  if (l_mode_vmc) then 
    call vmc_init
   elseif (l_mode_dmc) then
    call dmc_init
   else
    call die (lhere, 'If doing new optimization then mode must be vmc or dmc.')
  endif
  energy_plus_err_best=1.d99
  eloc_av_previous =  0.d0
  d_eloc_av = 0.d0
  l_convergence_reached = .false.
  convergence_reached_nb  = 0
  add_diag_mult_exp = 1.d0
  call object_modified ('add_diag_mult_exp')

! orthonormalization the orbitals
  if (l_ortho_orb_opt) then
    call ortho_orb
  endif

! Optimization method
  write(6,*)
  write(6,'(3a)') 'Optimization will be done with the ',trim(opt_method),' method.'
  write(6,*)

  if (l_opt_ovlp_fn .and. .not. l_opt_ovlp_branching) then
   l_branching = .false.
   write(6,'(a)') 'Warning: turn off branching in DMC for overlap fixed-node optimization method.'
  endif

  if (l_opt_ovlp_fn_linear) then
   l_stab = .false.
   write(6,'(a)') 'Warning: turn off stabilization for overlap-fn method with linear scaling.'
  endif

! choice of stabilization
  if (l_stab) then
   select case (trim(stabilization))
    case ('identity')
     write (6,'(a)') 'optimization will be stabilized by adding multiple of identity matrix'
    case ('overlap')
     write (6,'(a)') 'optimization will be stabilized by adding multiple of overlap matrix'
     if (.not. l_opt_lin) then
      call die (lhere, 'stabilization = overlap is only implemented for the linear optimization method')
     endif
    case ('symmetrize')
     write (6,'(a)') 'optimization will be stabilized by symmetrizing the hamiltonian matrix'
    case default
     call die (lhere, 'unknown stabilization choice >'+trim(stabilization)+'<.')
   end select
  else
   write(6,'(a)') 'stabilization turned off in optimization.'
  endif


! Nice printing
  write(6,'(a,i5,a,i5,a,i7,a,i5,a,i5,3a)') 'OPT: optimization of',nparmj+param_pjas_nb,' Jastrow,', nparmcsf,' CSF,',param_orb_nb,' orbital,', param_exp_nb, ' exponent and ',param_geo_nb," geometry parameters with ",trim(opt_method)," method:"


! Start optimization loop ------------------------------------------------------------------------
  do iter = 1, iter_opt_max_nb
   iter_global_optimization = iter

   write(6,'()')
   write(6,'(a,i3)') 'Beginning optimization iteration # ',iter

!  Define averages and statistical errors to compute:

!  For energy gradient
   call object_average_request ('dpsi_av')
   call object_average_request ('dpsi_eloc_av')
   call object_error_request ('gradient_norm_err')

!!!!! temporary for testing gradient
   if (l_opt_grad) then
     call object_average_request ('deloc_av') 

     call object_error_request ('deloc_av_err') 
     call object_error_request ('dpsi_eloc_covar_err') 
     call object_error_request ('dpsi_eloc_covar_deloc_av_err')
     call object_error_request ('dpsi_eloc_covar_err_2') 
     call object_error_request ('dpsi_eloc_covar_deloc_av_err_2')

!    more rigorous estimators of the errors
     call object_variance_request ('dpsi_eloc_av_var')
     call object_variance_request ('dpsi_av_var')
     call object_covariance_request ('dpsi_eloc_av_dpsi_av_covar')
     call object_covariance_request ('dpsi_eloc_av_eloc_av_covar')
     call object_covariance_request ('dpsi_av_eloc_av_covar')
     call object_variance_request ('deloc_av_var')
     call object_covariance_request ('deloc_av_dpsi_eloc_av_covar')
     call object_covariance_request ('deloc_av_dpsi_av_covar')
     call object_covariance_request ('deloc_av_eloc_av_covar')
   endif

!  For variance gradient and hessian
   if (p_var /= 0.d0) then
    call object_average_request ('deloc_eloc_av')
    call object_average_request ('dpsi_eloc_sq_av')
    call object_average_request ('deloc_deloc_av')

!   for linear Hessian
    if (trim(hessian_variance_type) == 'linear') then
     call object_average_request ('dpsi_deloc_eloc_av')
     call object_average_request ('dpsi_dpsi_eloc_sq_av')
    endif
   endif

!  Newton method
   if (l_opt_nwt) then
    call object_average_request ('dpsi_dpsi_av')
    call object_average_request ('deloc_av')
    call object_average_request ('dpsi_deloc_av')
    call object_average_request ('dpsi_dpsi_eloc_av')
    call object_average_request ('d2psi_av')
    call object_average_request ('d2psi_eloc_av')
   endif

!  linear method
   if (l_opt_lin) then
    call object_average_request ('dpsi_dpsi_av')
    call object_average_request ('deloc_av')
    call object_average_request ('dpsi_deloc_av')
    call object_average_request ('dpsi_dpsi_eloc_av')
   endif

! perturbative method
   if (l_opt_ptb) then

!   compute full overlap or just diagonal?
    if (l_diagonal_overlap) then
     call object_average_request ('dpsi_sq_av')
    else
     call object_average_request ('dpsi_dpsi_av')
    endif

!   compute energy denominators in QMC?
    if (iter == 1 .and. .not. l_opt_orb_eig) then
      call object_average_request ('deloc_av')
      call object_average_request ('dpsi_eloc_av')
      call object_average_request ('dpsi_deloc_av')
      call object_average_request ('dpsi_sq_eloc_av')
    else
      l_opt_orb_energy = .false.
    endif

   endif

! for overlap_fn method
    if (l_opt_ovlp_fn) then
      call object_average_request ('dpsi_uwav')
      call object_average_request ('dpsi_dpsi_uwav')
    endif

! for overlap_fn method with linear scaling
    if (l_opt_ovlp_fn_linear) then
      call object_average_request ('dpsi_over_jas2_av')
      call object_average_request ('dpsi_over_jas2_uwav')
      call object_average_request ('dpsi2_over_jas2_uwav')
      call object_average_request ('dpsi_dpsi_over_jas2_uwav') !temp
      call object_average_request ('first_csf_over_jas2_av')
      call object_average_request ('first_csf_over_jas2_uwav')
      call object_average_request ('first_csf2_over_jas2_uwav')
      call object_average_request ('one_over_jas2_uwav')
      call object_average_request ('one_over_jas2_av') !temp
    endif

!  request additional averages for bounds on dpsi and deloc
   if (l_deriv_bound) then
      call object_average_request ('dpsi_sq_av')
      call object_average_request ('deloc_sq_av')
   endif

!  set norb
   if (l_opt_orb) then
    call object_provide ('orb_opt_last_lab')
    norb = orb_opt_last_lab

!   orbital overlap
    if (l_ortho_orb_vir_to_orb_occ) then
     call ortho_orb_vir_to_orb_occ
     call object_provide ('orb_ovlp')
     call object_write ('orb_ovlp')
    endif

   endif

!  request vb_weights
   if (l_vb_weights) then
      call object_block_average_request ('product_csf_over_psid_bav')
      call object_average_request ('product_csf_over_psid_av')
      call object_block_average_request ('csf_over_psid_bav')
      call object_average_request ('vb_weights_chirgwin_coulson_av')
      call object_error_request ('vb_weights_chirgwin_coulson_av_err')
      call object_average_request ('vb_weights_lowdin_av')
      call object_error_request ('vb_weights_lowdin_av_err')
      call routine_write_final_request ('vb_weights_wrt')
   endif

!  QMC run
   nforce=1
   nwftype=1
   if(l_mode_vmc) then
     call vmc
    elseif(l_mode_dmc) then
     call dmc
    else
     call die(lhere, 'If doing new optimization then mode must be vmc or dmc.')
   endif

!  set back norb
   if (l_opt_orb) then
    call object_provide ('orb_occ_last_in_wf_lab')
    norb = orb_occ_last_in_wf_lab
   endif

!  Calculate and print gradient
   call object_provide ('gradient')
   call object_provide ('gradient_norm')
   call object_provide ('gradient_norm_err')
   call object_provide ('param_type')
   write(6,*)
   write(6,'(a)') 'Gradient with respect to the parameters:'
   do param_i = 1, param_nb
     write(6,'(a,i5,a,es15.8,3a)') 'gradient component # ',param_i,' : ', gradient (param_i), ' (',trim(param_type (param_i)),')'
   enddo
   write(6,'(a,es15.8,a,es15.8)') 'gradient norm :              ',gradient_norm, ' +- ',gradient_norm_err
   write(6,*)
   if (l_opt_ovlp_fn) then
     call object_provide ('gradient_ovlp')
     call object_provide ('gradient_ovlp_norm')
      write(6,*)
      write(6,'(a)') 'Gradient of overlap with respect to the parameters:'
      do param_i = 1, param_nb
        write(6,'(a,i5,a,es15.8,3a)') 'gradient ovlp component # ',param_i,' : ', gradient_ovlp (param_i), ' (',trim(param_type (param_i)),')'
      enddo
      write(6,'(a,es15.8,a,es15.8)') 'gradient ovlp norm :              ',gradient_ovlp_norm  !, ' +- ',gradient_ovlp_norm_err
      write(6,*)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  temporary printing for testing DMC gradient
   if (l_opt_grad) then
    write(6,'(a)') 'Gradient components (second error more rigorously calculated from covariances):'
    call object_provide ('dpsi_eloc_covar')
    call object_provide ('dpsi_eloc_covar_err')
    call object_provide ('dpsi_eloc_covar_err_2')
    call object_provide ('deloc_av')
    call object_provide ('deloc_av_err')
    call object_provide ('dpsi_eloc_covar_deloc_av')
    call object_provide ('dpsi_eloc_covar_deloc_av_err')
    call object_provide ('dpsi_eloc_covar_deloc_av_err_2')
    do  parm_i=1,param_nb
     write(6,'(i5,a,f10.6,a,i6,a,i6,a, a,f10.6,a,i6,a, a,f10.6,a,i6,a,i6,a)') &
       parm_i, ' 2*dpsi_eloc_covar=',2*dpsi_eloc_covar(parm_i), '(',nint(2*dpsi_eloc_covar_err(parm_i)*10**6),'|',nint(2*dpsi_eloc_covar_err_2(parm_i)*10**6),') ', &
      '2*deloc_av=',2*deloc_av(parm_i), '(',nint(2*deloc_av_err(parm_i)*10**6),') ',  &
      '2*dpsi_eloc_covar_deloc_av=',2*dpsi_eloc_covar_deloc_av(parm_i), '(',nint(2*dpsi_eloc_covar_deloc_av_err(parm_i)*10**6),'|',nint(2*dpsi_eloc_covar_deloc_av_err_2(parm_i)*10**6),')'
    enddo
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  check vanishing components or linear dependencies in gradient
   do parm_i = 1, param_nb
      if (abs(gradient (parm_i)) < 1.d-10) then
       l_warning = .true.
       write(6,'(a)') 'Warning: zero or very small gradient component:'
       write(6,'(a,i3,a,es15.8,a,i3,a,es15.8)') 'Warning: gradient (',parm_i,')=',gradient (parm_i)
       cycle
      endif
     do parm_j = parm_i+1, param_nb
      if (abs(gradient (parm_i) - gradient (parm_j)) < 1.d-10) then
       l_warning = .true.
       write(6,'(a)') 'Warning: possible linear dependency:'
       write(6,'(a,i3,a,es15.8,a,i3,a,es15.8)') 'Warning: gradient (',parm_i,')=',gradient (parm_i),' is identical or very close to gradient (',parm_j,')=',gradient (parm_j)
      endif
     enddo
   enddo

   if (l_opt_grad) then
    exit
   endif


!  calculate and print deloc_av_norm
   if (l_opt_lin .or. l_opt_nwt) then

    call object_provide ('deloc_av_abs_max')
    write(6,*)
    write(6,'(a,es15.8,a)') 'Maximum absolute value of local energy derivatives :', deloc_av_abs_max, ' (must be zero in VMC within statistical noise except for geometry optimization)'
   endif

!  calculate and print hessian
   if (l_opt_nwt) then
    call object_provide ('hess_nwt_eigval')
    write(6,*)
    write(6,'(a)') 'Hessian eigenvalues:'
    do param_i = 1, param_nb
     write(6,'(a,i5,a,es15.8)') 'eigenvalue # ',param_i,' : ', hess_nwt_eigval (param_i)
    enddo
   endif

!  variation of energy
   if (iter > 1) then
    d_eloc_av = energy(1) - eloc_av_previous
   endif

!  save current energy
   eloc_av_previous = energy(1)

!  If this is the best yet, save it.  Since we are primarily interested in the energy we always use
!  that as part of the criterion.  By adding in energy_err we favor those iterations where the energy
!  has a smaller error, either because of a reduction in sigma and Tcorr or because nblk is increasing.
!  If p_var!=0 then we add that to the criterion too.
   energy_plus_err=energy(1)+3*energy_err(1)+p_var*(energy_sigma(1)+3*error_sigma)
   if(energy_plus_err.lt.energy_plus_err_best) then
    iter_best = iter
    energy_plus_err_best=energy_plus_err
    energy_best = energy(1)
    call wf_best_save
   endif

!  check convergence
   if (l_check_convergence .and. iter > 1 .and. iter >= iter_opt_min_nb) then

!    convergence reached?
     if (dabs(d_eloc_av) <= energy_threshold .and. energy_err(1) <= energy_threshold/2.d0) then
       convergence_reached_nb = convergence_reached_nb + 1
     else
       convergence_reached_nb = 0
     endif

     if (convergence_reached_nb == check_convergence_nb) then
       l_convergence_reached = .true.
       if (.not. l_last_run) exit
     endif

   endif ! convergence

!  decrease error threshold provided previous move was not rejected
   if(move_rejected == 0) then
     if (l_increase_accuracy) then
       if (l_decrease_error) then
         if(l_decrease_error_adaptative) then
           if (d_eloc_av /= 0.d0) then
             error_threshold = min(energy_err(1),max(dabs(d_eloc_av)/decrease_error_factor,energy_err(1)/decrease_error_factor,decrease_error_limit))
           else
             error_threshold = energy_err(1)
           endif
         else
           error_threshold = max(energy_err(1)/decrease_error_factor,decrease_error_limit)
         endif
       endif
       if (l_increase_blocks) then
         nblk = nint(min(nblk*increase_blocks_factor,increase_blocks_limit))
       endif
     endif
   endif

!  If wave function got significantly worse, go back to previous wave function and increase diag_stab.
!  Unlike the criterion for selecting the best wavefn., during the optim. we want to reject the move only if it got much worse.
!  I have tried 2 criteria.  In both of them large increases in energy are penalized even if p_var=1.
!  The reason is that otherwise with variance minimization one can optimize to a very diffuse state
!  that may be an approximation to an excited state.  A possible solution to this is to optimize about
!  a fixed guessed energy, close to the ground-state energy, rather than about the sample average,
!  which is what is done in fit.
!  In the 1st criterion we penalize sigma even if p_var=0.
!  The 1st criterion has separate conditions for the energy and sigma, the 2nd criterion combines them.
!  In the first criterion:
!  For variance minimization we allow sigma to be at most 1.5 times worse
!  For energy   minimization we allow sigma to be at most 2.5 times worse
!  and
!  For variance minimization we allow the energy to be at most 6 std dev. worse
!  For energy   minimization we allow the energy to be at most 3 std dev. worse
!  The 2nd criterion is a bit closer to the objective function being optimized.
!  Since we we want to reject the move only if it got much worse, for the proposed move we
!  add 3*err and for the saved move 6*err and 10*err.  We add 10*error_sigma_sav because sigma
!  and particularly the error in sigma have large uncertainty, if the cusp is not imposed.
!  Since the 2nd criterion does not seem decisively better than the 1st, I am still using the
!  1st criterion, but leave the 2nd in the code, commented out.
   if(l_stab .and. iter > 1) then
!  1st criterion
   if(energy_sigma(1) > (2.5d0-p_var)*energy_sigma_sav .or.  &
     energy(1)-energy_sav > 3*(1+p_var)*(sqrt(energy_err(1)**2 + energy_err_sav**2))) then
!  2nd criterion
!  if((2-p_var)*(energy(1) +3*energy_err(1)) +p_var*(energy_sigma(1) +3*error_sigma)**2 > &
!     (2-p_var)*(energy_sav+6*energy_err_sav)+p_var*(energy_sigma_sav+10*error_sigma_sav)**2) then

     move_rejected=move_rejected+1

     call wf_restore
     if (l_opt_pjas) call restore_pjas
     call object_restore ('gradient')
     if(l_opt_nwt) then
      call object_restore ('hess_nwt')
     endif
     if(l_opt_lin) then
      call object_restore ('ham_lin')
      call object_restore ('ovlp_lin')
      call object_restore ('renorm_vector')
     endif
     if(l_opt_ptb) then
     if (l_diagonal_overlap) then
      call object_restore ('dpsi_sq_covar')
     else
      call object_restore ('dpsi_dpsi_covar_inv')
     endif
!        call object_restore ('delta_e_ptb') !!
     endif
     if(l_opt_ovlp_fn) then
      call object_restore ('dpsi_av')
      call object_restore ('dpsi_uwav')
      call object_restore ('dpsi_dpsi_uwcovar')
     endif
     diag_stab = min(100.d0 * diag_stab, add_diag_max)
     if(diag_stab == add_diag_max) call die (lhere, 'diag_stab too large')
     call object_modified ('diag_stab')
     write(6,'(a,1pd9.1)') 'Wave function got worse, increase add_diag up to ',diag_stab
     call wf_update_and_check_and_stab
!    write new wave function
     write(6,'(a)') ''
     write(6,'(''For next iteration #'',i4,'' new wave function:'')') iter+1
     call write_wf_new

!    just in case mc config is in crazy place, reset mc_configs by calling sites, only in VMC for now. Need to do it for DMC as well?
     if (l_mode_vmc) then
      isite = 1; call mc_configs_read
     endif

     if (l_decrease_error) then
      write(6,'(a,i3,t10,f12.7,a,f11.7,f10.5,f9.5,a,f9.5,f12.5,a,f9.5,f6.3,f9.5,1pd9.1)') 'OPT:',iter, energy(1),' +-', &
      energy_err(1), d_eloc_av, energy_sigma(1), ' +-', error_sigma, gradient_norm, ' +-', gradient_norm_err, p_var, error_threshold, diag_stab
     else
      write(6,'(a,i3,t10,f12.7,a,f11.7,f10.5,f9.5,a,f9.5,f12.5,a,f9.5,f6.3,i9,1pd9.1)') 'OPT:',iter, energy(1),' +-', &
      energy_err(1), d_eloc_av, energy_sigma(1), ' +-', error_sigma, gradient_norm, ' +-', gradient_norm_err, p_var, nblk, diag_stab
     endif
     if (l_opt_ovlp_fn) then
       write(6,'(a,i3,a,f13.6,a,f13.6,a,f12.5)') 'OPTd:',iter,' ovlp1=',ovlp_trial_fn, ' ovlp2=',ovlp_trial_fn_over_ovlp_trial, ' gradient_ovlp_norm= ', gradient_ovlp_norm
     endif
     cycle
    endif
   endif

   move_rejected=0

!  wave function got better or at least not significantly worse so save various quantities
   energy_sav=energy(1)
   energy_sigma_sav=energy_sigma(1)
   energy_err_sav=energy_err(1)
   error_sigma_sav = error_sigma
   ene_var_sav=(1-p_var)*energy(1)+p_var*energy_sigma(1)**2
   if (l_opt_ovlp_fn) then
    ovlp_trial_fn_sav = ovlp_trial_fn
    ovlp_trial_fn_over_ovlp_trial_sav = ovlp_trial_fn_over_ovlp_trial
   endif

!  save wavefunction, gradient, Hamiltonian and overlap
   call wf_save
   if (do_pjas) call save_pjas
   call object_save ('gradient')
   if(l_opt_nwt) then
    call object_save ('hess_nwt')
   endif
   if(l_opt_lin) then
      call object_save ('ham_lin')
      call object_save ('ovlp_lin')
      call object_save ('renorm_vector')
   endif
   if(l_opt_ptb) then
      if (l_diagonal_overlap) then
       call object_save ('dpsi_sq_covar')
      else
       call object_save ('dpsi_dpsi_covar_inv')
      endif
!      call object_save ('delta_e_ptb')  !
   endif
   if(l_opt_ovlp_fn) then
      call object_save ('dpsi_av')
      call object_save ('dpsi_uwav')
      call object_save ('dpsi_dpsi_uwcovar')
   endif

!  adjust diag_stab
   if (l_stab) then
     call adjust_diag_stab
   endif

!  update parameters
   call wf_update_and_check_and_stab

!  pretty printing
   write(6,*)
!   call object_provide ('sigma')
   call object_provide ('gradient_norm')
   call object_provide ('gradient_norm_err')
   if (l_opt_ovlp_fn) then
    call object_provide ('gradient_ovlp_norm')
   endif

   if (iter == 1) then
    if (l_decrease_error) then
     write(6,'(a)') 'OPT: iter    energy         error      diff          sigma                grad norm        p_var  nxt err  nxt stab'
    else
     write(6,'(a)') 'OPT: iter    energy         error      diff          sigma                grad norm        p_var  nxt nblk nxt stab'
    endif
   endif
   if (l_decrease_error) then
    if (.not. l_convergence_reached) then
     write(6,'(a,i3,t10,f12.7,a,f11.7,f10.5,f9.5,a,f9.5,f12.5,a,f9.5,f6.3,f9.5,1pd9.1)') 'OPT:',iter, energy_sav,' +-', &
     energy_err_sav, d_eloc_av, energy_sigma_sav, ' +-', error_sigma_sav, gradient_norm, ' +-', gradient_norm_err, p_var, error_threshold, diag_stab
    else
     write(6,'(a,i3,t10,f12.7,a,f11.7,f10.5,f9.5,a,f9.5,f12.5,a,f9.5,f6.3,f9.5,1pd9.1,a)') 'OPT:',iter, energy_sav,' +-', &
     energy_err_sav, d_eloc_av, energy_sigma_sav, ' +-', error_sigma_sav, gradient_norm, ' +-', gradient_norm_err, p_var, error_threshold, diag_stab,' converged'
    endif
   else
    if(.not. l_convergence_reached) then
     write(6,'(a,i3,t10,f12.7,a,f11.7,f10.5,f9.5,a,f9.5,f12.5,a,f9.5,f6.3,i9,1pd9.1)') 'OPT:',iter, energy_sav,' +-', &
     energy_err_sav, d_eloc_av, energy_sigma_sav, ' +-', error_sigma_sav, gradient_norm, ' +-', gradient_norm_err, p_var, nblk, diag_stab
    else
     write(6,'(a,i3,t10,f12.7,a,f11.7,f10.5,f9.5,a,f9.5,f12.5,a,f9.5,f6.3,i9,1pd9.1,a)') 'OPT:',iter, energy_sav,' +-', &
     energy_err_sav, d_eloc_av, energy_sigma_sav, ' +-', error_sigma_sav, gradient_norm, ' +-', gradient_norm_err, p_var, nblk, diag_stab,' converged'
    endif
   endif
   if (l_opt_ovlp_fn) then
      write(6,'(a,i3,a,f13.6,a,f13.6,a,f12.5)') 'OPTd:',iter,' ovlp1=',ovlp_trial_fn_sav, ' ovlp2=',ovlp_trial_fn_over_ovlp_trial_sav, ' gradient_ovlp_norm= ', gradient_ovlp_norm
   endif

!  decrease p_var
   if (l_decrease_p_var) then
     p_var = p_var/2.d0
     if (p_var < 0.01d0) p_var = 0.d0
     call object_modified ('p_var')
   endif

!  write new wave function
   write(6,'(a)') ''
   write(6,'(''For next iteration #'',i4,'' new wave function:'')') iter+1
   call write_wf_new

   if (l_convergence_reached) exit

  enddo ! end optimization loop

  if (l_opt_grad) then
   write(6,*)
   write(6,'(a)') 'End of gradient calculation'
   return
  endif

! final printing
  write(6,*) ''
  write(6,'(a,i3,a)') 'Optimization ended after ',iter,' iterations.'

  if (l_convergence_reached) then
   write(6,'(a)') 'Convergence reached.'
   write(6,'(a,f12.7,a,i2,a)') 'Threshold on energy ', energy_threshold,' reached for ', check_convergence_nb,' consecutive steps.'
   if (.not. l_last_run) then
    write(6,*)
    write(6,'(a,i3,t10,f12.7,a,f11.7,f10.5,f9.5,a,f9.5,f12.5,a,f9.5,f6.3,a)') 'OPT:',iter,energy(1),' +-',energy_err(1), d_eloc_av, energy_sigma(1), ' +-', error_sigma, gradient_norm, ' +-', gradient_norm_err, p_var, '      converged'
    if (l_opt_ovlp_fn) then
      write(6,'(a,i3,a,f13.6,a,f13.6,a,f12.5)') 'OPTd:',iter,' ovlp1=',ovlp_trial_fn, ' ovlp2=',ovlp_trial_fn_over_ovlp_trial, ' gradient_ovlp_norm= ', gradient_ovlp_norm
    endif
   endif
   iter = iter + 1
  else
   l_warning = .true.
   write(6,'(a)') 'Warning: Convergence not reached.'
   write(6,'(2a,i3,a)') trim(lhere),': Maximun number of iterations ',  iter_opt_max_nb,' reached.'
  endif

! write final wave function
  if (.not. l_last_run) then
   write(6,*)
   write(6,'(a)') 'Wave function at final iteration:'
   call write_wf_new
  endif

! do a last vmc with the last predicted parameters without calculating the derivatives
  if (l_last_run) then
  write(6,*)
  write(6,'(a)') 'Performing last qmc run'

! request vb_weights
  if (l_vb_weights) then
      l_lastrun_global_optimization = l_last_run
      call object_block_average_request ('product_csf_over_psid_bav')
      call object_block_average_request ('csf_over_psid_bav')
      call object_average_request ('vb_weights_chirgwin_coulson_av')
      call object_error_request ('vb_weights_chirgwin_coulson_av_err')
      call object_average_request ('vb_weights_lowdin_av')
      call object_error_request ('vb_weights_lowdin_av_err')
      call routine_write_final_request ('vb_weights_wrt')
  endif

  nforce=1
  nwftype=1
  if(l_mode_vmc) then
    call vmc
   elseif(l_mode_dmc) then
    call dmc
   else
    call die(lhere, 'If doing new optimization then mode must be vmc or dmc.')
  endif

  d_eloc_av = energy(1) - eloc_av_previous
  write(6,*)
  write(6,'(a,i3,t10,f12.7,a,f11.7,f10.5,f9.5,a,f9.5,a)') 'OPT:',iter,energy(1),' +-',energy_err(1), d_eloc_av, energy_sigma(1), ' +-', error_sigma,'                                    last run'
  if (l_opt_ovlp_fn) then
   write(6,'(a,i3,a,f13.6,a,f13.6)') 'OPTd:',iter,' ovlp1=',ovlp_trial_fn, ' ovlp2=',ovlp_trial_fn_over_ovlp_trial
  endif

! If this is the best yet, save it.  Since we are primarily interested in the energy we use
! that as part of the criterion.  By adding in energy_err we favor those iterations where the energy
! has a smaller error, either because of a reduction in sigma and Tcorr or because nblk is increasing.
! If p_var!=0 then we add that to the criterion too.
  energy_plus_err=energy(1)+3*energy_err(1)+p_var*(energy_sigma(1)+3*error_sigma)
  if(energy_plus_err.lt.energy_plus_err_best) then
    iter_best = iter
    energy_plus_err_best=energy_plus_err
    energy_best = energy(1)
    call wf_best_save
  endif
  endif !l_last_run

!

! Print best wave function
  write(6,*)
! write(6,'(a,i3)') 'OPT: the best wave function was found at iteration # ',iter_best
  write(6,'(a,i3,f13.7,a,f11.7)') 'OPT: the best wave function was found at iteration # ',iter_best,energy_best,' +-',(energy_plus_err_best-energy_best)/3
  write(6,'(a)') 'Best wave function:'
  call write_wf_best

! restore best wave function and set new etrial for subsequent calculations
  call wf_best_restore
  etrial = energy_best

  end subroutine optimization

!===========================================================================
  subroutine wf_update_and_check (is_bad_move, is_bad_move_exp, exp_move_big)
!---------------------------------------------------------------------------
! Description : update and check parameters of wave function
!
! Output      : is_bad_move = 1 if move is bad
!
! Created     : J. Toulouse, 18 Jan 2006
!---------------------------------------------------------------------------
  include 'modules.h'
  use dim_mod
  implicit none

! output
  integer, intent(out) :: is_bad_move, is_bad_move_exp
  logical, intent(out) :: exp_move_big

! local
  integer bas_i
  integer i, ict, isp, iparmcsf, iparm
  integer dexp_i, dexp_to_all_bas_i
  real(dp), parameter :: AMAX_NONLIN = 100.d0
  integer exponent_negative_nb
  real(dp) parm2min
  integer iparmpjase
  integer force_i, cent_i, dim_i

! begin

! Update parameters:

! CSF parameters
  if (l_opt_csf) then
   call object_provide ('ncsf')
   call object_provide ('iwcsf')
   call object_provide ('delta_csf')
   do iparmcsf = 1, nparmcsf
     csf_coef(iwcsf(iparmcsf),iwf)=csf_coef(iwcsf(iparmcsf),1) + delta_csf (iparmcsf)
   enddo
   call object_modified ('csf_coef')
  endif ! l_opt_csf

! Jastrow parameters
  if (l_opt_jas) then

    call object_provide ('nctype')
    call object_provide ('nparma')
    call object_provide ('nparmb')
    call object_provide ('nparmc')
    call object_provide ('iwjasa')
    call object_provide ('iwjasb')
    call object_provide ('iwjasc')
    call object_provide ('a4')
    call object_provide ('b')
    call object_provide ('c')
    call object_provide ('delta_jas')

    iparm = 0

!   e-N term
    do ict=1,nctype
       do i=1,nparma(ict)
          iparm=iparm+1
          a4(iwjasa(i,ict),ict,iwf)=a4(iwjasa(i,ict),ict,1) + delta_jas (iparm)
! warning: remove this stop.
!          if(iwjasa(i,ict) == 2 .and. a4(iwjasa(i,ict),ict,iwf) > AMAX_NONLIN) then
!           call die (lhere, 'probably do not want a(2) > AMAX_NONLIN')
!          endif
       enddo
     enddo

!    e-e term
     do isp=nspin1,nspin2b
        do i=1,nparmb(isp)
          iparm=iparm+1
          b(iwjasb(i,1),isp,iwf)=b(iwjasb(i,1),isp,1) + delta_jas (iparm)
! warning: remove this stop. Is it good?
!          if(iwjasb(i,1) == 2 .and. b(iwjasb(i,1),isp,1) > AMAX_NONLIN) then
!            call die (lhere, 'probably do not want b(2) > AMAX_NONLIN')
!          endif
        enddo
     enddo

!     e-e-N term
      do ict=1,nctype
        do i=1,nparmc(ict)
          iparm=iparm+1
          c(iwjasc(i,ict),ict,iwf)=c(iwjasc(i,ict),ict,1) + delta_jas (iparm)
        enddo
      enddo

!     find not optimized e-e-N coefficients by nullifying the contribution of the e-e-N
!     to e-e and e-N cusps
      if(ijas == 4 .and. isc <= 9) call cuspexact4(1,iwf)

  endif ! l_opt_jas

! Periodic Jastrow  parameters
  if (l_opt_pjas) then
     call object_provide ('delta_pjas')
     do iparmpjase = 1, param_pjas_nb
        pjas_parms (iparmpjase,iwf)=pjas_parms (iparmpjase,1) + delta_pjas (iparmpjase)
     enddo
     call object_modified ('pjas_parms')
  endif ! l_opt_pjas

! Exponent parameters
  if (l_opt_exp) then

   call object_provide ('param_exp_nb')
   call object_provide ('dexp_to_all_bas_nb')
   call object_provide ('dexp_to_all_bas')
   call object_provide ('zex')
   call object_provide ('delta_exp')

   do dexp_i = 1, param_exp_nb
     do dexp_to_all_bas_i = 1, dexp_to_all_bas_nb (dexp_i)
        bas_i = dexp_to_all_bas (dexp_i)%row (dexp_to_all_bas_i)
        if (l_optimize_log_exp) then
          zex (bas_i, iwf) = exp(log(zex (bas_i, 1)) + delta_exp (dexp_i))
!          zex (bas_i, iwf) = zex (bas_i, 1) * (1.d0 + delta_exp (dexp_i))
        else
          zex (bas_i, iwf) = zex (bas_i, 1) + delta_exp (dexp_i)
        endif
     enddo ! dexp_to_all_bas_i
   enddo ! dexp_i

   call object_modified ('zex')

!  recompute exponents zex2
   call copy_zex_zex2

!  update orbital coefficients due to exponent variations
   select case (trim(basis_functions_varied))
   case ('unnormalized')
    call coef_orb_on_norm_basis_from_coef (1)
   case ('normalized')
    call coef_from_coef_orb_on_norm_basis (1)
   case ('orthonormalized')
    call coef_from_coef_orb_on_ortho_basis (1)
    call coef_orb_on_norm_basis_from_coef (iwf)
   case default
    call die (here, 'unknown case >'+trim(basis_functions_varied)+'< for basis_functions_varied.')
   end select

  endif ! l_opt_exp

! Orbital parameters
  if (l_opt_orb) then

!  non-periodic case
   if (iperiodic == 0) then

     if (l_opt_exp) then
     select case (trim(basis_functions_varied))
     case ('unnormalized')
      call update_coef_by_rot
      call coef_orb_on_norm_basis_from_coef (iwf)
     case ('normalized')
      call update_coef_orb_on_norm_basis_by_rot
      call coef_from_coef_orb_on_norm_basis (iwf)
     case ('orthonormalized')
      call update_coef_orb_on_ortho_basis_by_rot
      call coef_from_coef_orb_on_ortho_basis (iwf)
      call coef_orb_on_norm_basis_from_coef (iwf)
     case default
      call die (here, 'unknown case >'+trim(basis_functions_varied)+'< for basis_functions_varied.')
     end select
     else
      call update_coef_by_rot
      call coef_orb_on_norm_basis_from_coef (iwf)
     endif

!    checking or imposing e-n cusp conditions
!    Warning: I should check if something is needed to be done for exponent optimization (coef_from_coef_orb_on_norm_basis?)
     if (l_check_cusp_en .or. l_impose_cusp_en_opt) then
      call cusp_en_orb
     endif

!  periodic case
   else

!      write(6,*) "ngvec_orb = ", ngvec_orb

!    provide needed objects
     call object_provide ('ngvec_orb')
     call object_provide ('orb_tot_nb')
     call object_provide ('c_rp')
     call object_provide ('c_rm')
     call object_provide ('c_ip')
     call object_provide ('c_im')
     call object_provide ('delta_c_rp')
     call object_provide ('delta_c_rm')
     call object_provide ('delta_c_ip')
     call object_provide ('delta_c_im')

!    update orbital coefficients
     c_rp (1:ngvec_orb, 1:orb_tot_nb) = c_rp (1:ngvec_orb, 1:orb_tot_nb) + delta_c_rp (:,:)
     c_rm (1:ngvec_orb, 1:orb_tot_nb) = c_rm (1:ngvec_orb, 1:orb_tot_nb) + delta_c_rm (:,:)
     c_ip (1:ngvec_orb, 1:orb_tot_nb) = c_ip (1:ngvec_orb, 1:orb_tot_nb) + delta_c_ip (:,:)
     c_im (1:ngvec_orb, 1:orb_tot_nb) = c_im (1:ngvec_orb, 1:orb_tot_nb) + delta_c_im (:,:)

!    coefficients have been modified
     call object_modified ('c_rp')
     call object_modified ('c_rm')
     call object_modified ('c_ip')
     call object_modified ('c_im')

   endif ! if iperiodic == 0

  endif ! l_opt_orb

! Geometry parameters
  if (l_opt_geo) then

   call object_provide ('ncent')
   call object_provide ('nwf')
   call object_provide ('znuc')
   call object_provide ('iwctype')
   call object_provide ('ndim')
   call object_provide ('param_geo_nb')
   call object_provide ('cent')
   call object_provide ('delta_geo')

   call object_alloc ('cent2', cent2, ndim, ncent, 3)
   cent2 (1:ndim,1:ncent,iwf) = cent (1:ndim,1:ncent)
   do force_i = 1, param_geo_nb
    cent_i = forces_cent (force_i) 
    dim_i = forces_direct (force_i) 
    cent2 (dim_i, cent_i, iwf) = cent (dim_i, cent_i) + delta_geo (force_i)
   enddo
   call object_modified ('cent2')

!  recalculating nuclear potential energy
   call alloc ('pecentn', pecentn, nwf)
   call pot_nn(cent2(:,:,iwf),znuc,iwctype,ncent,pecentn(iwf))
   call object_modified ('pecentn')

   if (iwf == 1) then
    cent (1:ndim,1:ncent) = cent2 (1:ndim,1:ncent,1)
    pecent=pecentn(1) 
    call object_modified ('cent')
    call object_modified ('pecent')
   endif ! iwf == 1
   
  endif ! l_opt_geo

! check if really correct to call after each update
  call set_scale_dist(ipr,iwf)

! print new parameters
  call write_wf

! Check move:
  is_bad_move = 0
  is_bad_move_exp = 0
  exp_move_big = .false.

! test norm of csf coefficient variations
  if (l_opt_csf) then
    call object_provide ('delta_csf_norm')
    if (delta_csf_norm > 1.d0) then
      is_bad_move = 1
      write (6,'(a,es15.8,a)') 'This is a bad move because the norm of the csf coefficient variations is too large: delta_csf_norm=',delta_csf_norm,' > 1.'
    endif
  endif

! test norm of jastrow parameter variations
  if (l_opt_jas) then
    call object_provide ('delta_jas_norm')
    if (ndim == 3) then
      if(delta_jas_norm > max(10.d0,5.d0/scalek(iwf))) then
        is_bad_move = 1
        write (6,'(a,es15.8,a)') 'This is a bad move because the norm of the jastrow parameter variations is too large: delta_jas_norm=',delta_jas_norm,' > 10 or 5/scalek'
      endif
    else
      if(delta_jas_norm > max(500.d0,50.d0/scalek(iwf))) then
        is_bad_move = 1
        write (6,'(a,es15.8,a)') 'This is a bad move because the norm of the jastrow parameter variations is too large: delta_jas_norm=',delta_jas_norm,' > 500 or 50/scalek'
      endif
    endif
  endif

! test nonlinear Jastrow parameters
  if (l_opt_jas) then
   if (isc.ne.8 .and. isc.ne.10) then
    parm2min=-scalek(1)
   else
    parm2min=-1.d0
   endif
  do ict = 1, nctype
    if (a4(2,ict,iwf) < parm2min) then
     is_bad_move = 2
     write (6,'(a,es15.8,a,es15.8)') 'This is a bad move because a2=',a4(2,ict,iwf),' < parm2min=',parm2min
    endif
    if (a4(2,ict,iwf) > AMAX_NONLIN) then
     is_bad_move = 1
     write (6,'(a,es15.8,a,es15.8)') 'This is a bad move because a2=',a4(2,ict,iwf),' > AMAX_NONLIN=',AMAX_NONLIN
    endif
  enddo

  do isp = nspin1, nspin2b
   if (b(2,isp,iwf) < parm2min) then
     is_bad_move = 2
     write (6,'(a,es15.8,a,es15.8)') 'This is a bad move because b2=',b(2,isp,iwf),' < parm2min=',parm2min
   endif
   if (b(2,isp,iwf) > AMAX_NONLIN) then
     is_bad_move = 1
     write (6,'(a,es15.8,a,es15.8)') 'This is a bad move because b2=',b(2,isp,iwf),' > AMAX_NONLIN=',AMAX_NONLIN
   endif
  enddo
  endif ! l_opt_jas

! test exponents
  if (l_opt_exp) then
    exponent_negative_nb = 0
    do dexp_i = 1, param_exp_nb
      do dexp_to_all_bas_i = 1, dexp_to_all_bas_nb (dexp_i)
        bas_i = dexp_to_all_bas (dexp_i)%row (dexp_to_all_bas_i)
        if (zex (bas_i, iwf) < 0.d0) then
          is_bad_move = 2
          exponent_negative_nb = exponent_negative_nb + 1
        endif
        if (do_add_diag_mult_exp) then
           if (abs (delta_exp (dexp_i)) > 0.05 * zex (bas_i, 1)) then
              exp_move_big = .true.
           endif
           if (abs (delta_exp (dexp_i)) > 0.2 * zex (bas_i, 1)) then
              is_bad_move_exp = 1
              write(6,'(a,f10.3,a,I4,a,f10.1,a)') "This is a bad move because change in exponent ", zex (bas_i, 1), &
                   &  " of basis ", bas_i , " is ", 100 * abs (delta_exp (dexp_i))/ zex (bas_i, 1), " percent"
           endif
        endif ! do_add_diag_mult_exp
      enddo ! dexp_to_all_bas_i
    enddo ! dexp_i
    if (exponent_negative_nb > 0) then
       write (6,'(a,i3,a)') 'This is a bad move because ',exponent_negative_nb,' exponents are negative.'
    endif
  endif ! l_opt_exp

! check total parameter variations norm
  call object_provide ('delta_param_norm')
  if (ndim == 3) then
    if (delta_param_norm > delta_param_norm_max_3d) then
      is_bad_move = max(is_bad_move,1)
      write (6,'(a,es15.8,a,es15.8)') 'This is a bad move because the norm of the parameter variations is too large: delta_param_norm=',delta_param_norm,' > ',delta_param_norm_max_3d
    endif
   else
    if (delta_param_norm > delta_param_norm_max_2d) then
      is_bad_move = max(is_bad_move,1)
      write (6,'(a,es15.8,a,es15.8)') 'This is a bad move because the norm of the parameter variations is too large: delta_param_norm=',delta_param_norm,' > ',delta_param_norm_max_2d
    endif
  endif

! test norm of linear wave function variation
  if (l_opt_lin) then
   if (psi_lin_var_norm > psi_lin_var_norm_max) then
    is_bad_move = max(is_bad_move,1)
    write (6,'(a,es15.8,a,es15.8)') 'This is a bad move because the norm of the linear wave function variation is too large: psi_lin_var_norm=',psi_lin_var_norm,' > ',psi_lin_var_norm_max
   endif
  endif

!  call write_wf

  end subroutine wf_update_and_check

!===========================================================================
  subroutine wf_update_and_check_and_stab
!---------------------------------------------------------------------------
! Description : update, check parameters of wave function
! Description : and stabilize if necessary
!
! Created     : J. Toulouse, 18 Apr 2007
!---------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'wf_update_and_check_and_stab'
  integer is_bad_move, loop
  integer is_bad_move_exp
  logical exp_move_big

! begin

! loop over updates of wave function
  loop = 0
  do
   loop = loop + 1
   if (loop >= 50) then
    call die (lhere, 'move is still rejected after increasing add_diag 50 times!')
   endif

!  update and check wave function
   call wf_update_and_check (is_bad_move, is_bad_move_exp, exp_move_big)

!  if move for the exponents was not big, decrease add_diag_mult_exp
   if (l_stab .and. l_opt_exp .and. do_add_diag_mult_exp .and. .not. exp_move_big) then
       if (add_diag_mult_exp  > 1.0_dp) then
         add_diag_mult_exp = max (1.0_dp, add_diag_mult_exp/10.d0)
         write(6,'(a,f10.1)') "add_diag_mult_exp is decreased to", add_diag_mult_exp
         call object_modified ('add_diag_mult_exp')
       endif
   endif

!  if the move is bad, increase add_diag and retry
   if ((l_stab .and. is_bad_move == 1) .or. (is_bad_move == 2)) then !warning: is_bad_move == 2  forces stabilization
     call wf_restore
     if (l_opt_pjas) call restore_pjas
     diag_stab = min(diag_stab * 10.d0, add_diag_max)
     if(diag_stab == add_diag_max) call die (lhere, 'diag_stab too large')
     call object_modified ('diag_stab')
     write(6,'(a,1pd9.1)') 'increasing add_diag up to', diag_stab
     cycle
   endif

!  if move is bad only for the exponents, increase add_diag_mult_exp and retry
   if (l_stab .and. l_opt_exp .and. do_add_diag_mult_exp .and. is_bad_move_exp == 1) then
       call wf_restore
       if (l_opt_pjas) call restore_pjas
       add_diag_mult_exp = add_diag_mult_exp * 10
       write(6,'(a,f10.1)') "add_diag_mult_exp is increased to ", add_diag_mult_exp
       call object_modified ('add_diag_mult_exp')
     cycle
    endif

!   if the move was good, exit loop
    exit

  enddo ! end loop

! orthonormalization the orbitals
  if (l_ortho_orb_opt) then
    call ortho_orb
  endif

  end subroutine wf_update_and_check_and_stab

!===========================================================================
  subroutine adjust_diag_stab
!---------------------------------------------------------------------------
! Description : find optimal value of diag_stab
!
! Created     : J. Toulouse, 21 Jan 2006
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'adjust_diag_stab'
  integer is_bad_move_1, is_bad_move_2, is_bad_move_3
  integer iadd_diag
  real(dp) error_threshold_save
  real(dp) :: ene_var (3)
  integer loop
  integer nblk_sav, nblk_small
  logical :: calculation_reliable (3)
  integer :: calculation_reliable_nb
  integer :: is_bad_move_exp1, is_bad_move_exp2, is_bad_move_exp3
  logical :: exp_move_big1, exp_move_big2, exp_move_big3

! begin
  write(6,*)
  write(6,'(a)') 'Searching for optimal stabilizing add_diag...'

! reset add_diag to (small) fixed value, so that the found optimal add_diag will tend to be closer to this (small) fixed value
  if (l_reset_add_diag) then
    diag_stab = add_diag_reset_value
    call object_modified ('diag_stab')
    write(6,'(a,1pd9.1)') 'Resetting add_diag to add_diag_reset_value = ',diag_stab
  endif

! smaller value of number of blocks
  nblk_small=min(nblk_max,max(10,block_nb/10))

! Duplicate 3 times the wave function for use in correlated sampling
  nforce=3
  nwftype=3
  call wf_copy2
  call wf_copy
  if (do_pjas) call copy_pjas  ! WAS
! initialize asymptotic values related to scalek(2) and scalek(3)
  call set_scale_dist(ipr,2)
  call set_scale_dist(ipr,3)

! make sure that diag_stab is not tiny compared to the smallest eigenvalue of Hessian
  if (l_opt_nwt) then
    call object_provide ('hess_nwt_eigval_min')
    diag_stab = min(max(diag_stab,1.d-1*hess_nwt_eigval_min),add_diag_max)
    if(diag_stab == add_diag_max) call die (lhere, 'diag_stab too large')
    call object_modified ('diag_stab')
  endif

  add_diag (1) = diag_stab

! Compute 3 optimized wavefunctions with different values of diag_stab
  loop = 0
  do

   loop = loop + 1

   if (loop >= 15) then
    call die (lhere, 'add_diag has been increased 15 times but moves are still too bad')
   endif

   iwf = 2
   add_diag (2) = 0.1d0 * add_diag (1)
   diag_stab = add_diag (2)
   call object_modified ('diag_stab')
   write(6,*)
   write(6,'(a,i1,a,1pd9.1)') 'Trying add_diag (', iwf , ')=', diag_stab
   call wf_update_and_check (is_bad_move_2, is_bad_move_exp2, exp_move_big2)

   iwf = 3
   add_diag (3) = 10.d0 * add_diag (1)
   diag_stab = add_diag (3)
   call object_modified ('diag_stab')
   write(6,*)
   write(6,'(a,i1,a,1pd9.1)') 'Trying add_diag (', iwf , ')=', diag_stab
   call wf_update_and_check (is_bad_move_3, is_bad_move_exp3, exp_move_big3)

   iwf = 1
   diag_stab =  add_diag (1)
   call object_modified ('diag_stab')
   write(6,*)
   write(6,'(a,i1,a,1pd9.1)') 'Trying add_diag (', iwf , ')=', diag_stab
   call wf_update_and_check (is_bad_move_1, is_bad_move_exp1, exp_move_big1)

!  update nuclear coordinates
   if (l_opt_geo) then
    call object_provide ('ndim')
    call object_provide ('ncent')
    call object_provide ('nwf')
    call object_provide ('cent2')
    call alloc ('cent_ref', cent_ref, 3, ncent)
    call alloc ('delc', delc, 3, ncent, nwf)
    cent_ref (1:ndim,1:ncent) = cent2 (1:ndim,1:ncent,1)
    delc (1:ndim,1:ncent,2) = cent2 (1:ndim,1:ncent,2) - cent_ref (1:ndim,1:ncent)
    delc (1:ndim,1:ncent,3) = cent2 (1:ndim,1:ncent,3) - cent_ref (1:ndim,1:ncent)
    call object_modified ('cent_ref')
    call object_modified ('delc')
   endif

!  if move too large, restore wave function and increase diag_stab_ref, otherwise exit loop
   if (is_bad_move_1 /= 0 .or. is_bad_move_2 /= 0 .or. is_bad_move_3 /= 0) then
    call wf_restore
    if (l_opt_pjas) call restore_pjas
    add_diag (1) = add_diag (1) * 10**(is_bad_move_1 + is_bad_move_2 + is_bad_move_3)
    write(6,'(a,1pd9.1)') 'Increasing add_diag to ',add_diag (1)
    cycle
   endif

   if (l_opt_exp .and. do_add_diag_mult_exp) then

      if (.not. exp_move_big1 .and. .not. exp_move_big2  .and.  .not. exp_move_big3) then
         if (add_diag_mult_exp  > 1.0_dp) then
            add_diag_mult_exp = max (1.0_dp, add_diag_mult_exp /10.d0)
            write(6,'(a,f10.1)') "add_diag_mult_exp is decreased to ", add_diag_mult_exp
            call object_modified ('add_diag_mult_exp')
         endif
      endif

      if (is_bad_move_exp1 == 1 .or. is_bad_move_exp2 == 1 .or. is_bad_move_exp3 == 1) then
         call wf_restore
         if (l_opt_pjas) call restore_pjas
         add_diag_mult_exp = add_diag_mult_exp * 10
         write(6,'(a,f10.1)') "add_diag_mult_exp is increased to ", add_diag_mult_exp
         call object_modified ('add_diag_mult_exp')
         cycle
      endif

   endif ! l_opt_exp .and. do_add_diag_mult_exp

!  3 small correlated sampling without computing gradient
!   igradhess=0
   error_threshold_save = error_threshold
!   error_threshold = 10.d0 * error_threshold
   error_threshold = 1.d30
   nblk_sav=nblk
   nblk=nblk_small
   if(l_mode_vmc) then
     call vmc
    elseif(l_mode_dmc) then
     call dmc
    else
    call die(lhere, 'If doing new optimization then mode must be vmc or dmc.')
   endif
   nblk=nblk_sav
   error_threshold = error_threshold_save

!  check whether correlated calculations are reliable by looking at their autocorrelation times
   write (6,*)
   calculation_reliable_nb = 3
   calculation_reliable (1:3) = .true.
   call object_provide ('eloc_tc')
   do iadd_diag = 1, 3
    if (eloc_tc (iadd_diag) >= 50.d0) then
     calculation_reliable_nb = calculation_reliable_nb - 1
     calculation_reliable (iadd_diag) = .false.
     l_warning = .true.
     write (6,'(a,i1,a)') 'Warning: correlated calculation # ',iadd_diag,' has Tc > 50 and thus will be ignored.'
    endif
   enddo

!  if the central calculation is not reliable, increase add_diag and cycle
   if (.not. calculation_reliable (1)) then
    add_diag (1) =add_diag (1) * 10.d0
    write(6,'(a,1pd9.1)') 'Correlated calculation # 1 is not reliable, increase add_diag to ',add_diag (1)
!   just in case mc config is in crazy place, reset mc_configs by calling sites
    if (l_mode_vmc) then
      isite=1; call mc_configs_read
    endif
    cycle
   endif

! If wave function got significantly worse, go back to previous wave function and increase diag_stab.
! Unlike the criterion for selecting the best wavefn., during the optim. we want to reject the move only if it got much worse.
! I have tried 2 criteria.  In both of them large increases in energy are penalized even if p_var=1.
! The reason is that otherwise with variance minimization one can optimize to a very diffuse state
! that may be an approximation to an excited state.  A possible solution to this is to optimize about
! a fixed guessed energy, close to the ground-state energy, rather than about the sample average,
! which is what is done in fit.
! In the 1st criterion we penalize sigma even if p_var=0.
! The 1st criterion has separate conditions for the energy and sigma, the 2nd criterion combines them.
! In the first criterion:
! For variance minimization we allow sigma to be at most 1.5 times worse
! For energy   minimization we allow sigma to be at most 2.5 times worse
! and
! For variance minimization we allow the energy to be at most 6 std dev. worse
! For energy   minimization we allow the energy to be at most 3 std dev. worse
! The 2nd criterion is a bit closer to the objective function being optimized.
! Since we we want to reject the move only if it got much worse, for the proposed move we
! add 3*err and for the saved move 6*err and 10*err.  We add 10*error_sigma_sav because sigma
! and particularly the error in sigma have large uncertainty, if the cusp is not imposed.
! Since the correlated sampling runs are 10 times shorter than the others, include an extra factor of 3.3
! Since the 2nd criterion does not seem decisively better than the 1st, I am still using the
! 1st criterion, but leave the 2nd in the code, commented out.
! I used to reject if the central calculation (1) or the one with larger add_diag (3) were bad.
! Then we started to reject only if (1) is bad.
! Now we reject only if all 3 are bad.

! 1st criterion
!  if(energy_sigma(1) > (2.5d0-p_var)*energy_sigma_sav .or.     &
!     energy_sigma(3) > (2.5d0-p_var)*energy_sigma_sav .or.     &
!     energy(1)-energy_sav > 3*(1+p_var)*(sqrt(energy_err(1)**2+energy_err_sav**2)) .or. &
!     energy(3)-energy_sav > 3*(1+p_var)*(sqrt(energy_err(3)**2+energy_err_sav**2))) then

! 1st criterion
!  check only central calculation
!  if(energy_sigma(1) > (2.5d0-p_var)*energy_sigma_sav .or.     &
!    energy(1)-energy_sav > 3*(1+p_var)*(sqrt(energy_err(1)**2+energy_err_sav**2))) then
!  check all 3 and reject only if all are bad
   if((energy_sigma(1) > (2.5d0-p_var)*energy_sigma_sav .or.     &
     energy(1)-energy_sav > 3*(1+p_var)*(sqrt(energy_err(1)**2+energy_err_sav**2))) .and. &
     (energy_sigma(2) > (2.5d0-p_var)*energy_sigma_sav .or.     &
     energy(2)-energy_sav > 3*(1+p_var)*(sqrt(energy_err(2)**2+energy_err_sav**2))) .and. &
     (energy_sigma(3) > (2.5d0-p_var)*energy_sigma_sav .or.     &
     energy(3)-energy_sav > 3*(1+p_var)*(sqrt(energy_err(3)**2+energy_err_sav**2)))) then
! 2nd criterion
!  if((2-p_var)*(energy(1) +3*energy_err(1))     +p_var*(energy_sigma(1) +3*error_sigma)**2 > &
!     (2-p_var)*(energy_sav+6*3.3*energy_err_sav)+p_var*(energy_sigma_sav+10*3.3*error_sigma_sav)**2) then

    add_diag (1) = add_diag (1) * 100.d0
    write(6,'(a,1pd9.1)') 'Energy or sigma of correlation calculation # 1 went up too much, increase add_diag to ',add_diag (1)

    call wf_restore
    if (l_opt_pjas) call restore_pjas
    call object_restore ('gradient')
    if(l_opt_nwt) then
     call object_restore ('hess_nwt')
    endif
    if(l_opt_lin) then
     call object_restore ('ham_lin')
     call object_restore ('ovlp_lin')
     call object_restore ('renorm_vector')
    endif
    if(l_opt_ptb) then
     if (l_diagonal_overlap) then
      call object_restore ('dpsi_sq_covar')
     else
      call object_restore ('dpsi_dpsi_covar_inv')
     endif
!       call object_restore ('delta_e_ptb')  !
    endif
    if(l_opt_ovlp_fn) then
      call object_restore ('dpsi_av')
      call object_restore ('dpsi_uwav')
      call object_restore ('dpsi_dpsi_uwcovar')
    endif

!   just in case mc config is in crazy place, reset mc_configs by calling sites
    if (l_mode_vmc) then
      isite=1; call mc_configs_read
    endif

    cycle
   endif

   exit
  enddo

! if correlated calculations # 2 or 3 not reliable, but correlated calculation # 1 is reliable, then accept the current value of add_diag
! this may be dangerous, maybe need to be modified

! if the 3 correlated calculations are reliable, find optimal add_diag by parabolic interpolation
  if (calculation_reliable_nb == 3) then
!  This is the objective function being optimized
   do iadd_diag=1,3
     ene_var(iadd_diag)=(1.d0-p_var)*energy(iadd_diag)+p_var*energy_sigma(iadd_diag)**2
   enddo

!  Find optimal add_diag
!  call quad_min(energy_sav,energy_err_sav,energy,energy_err,force,force_err,ene_var,add_diag,3,0.d0,0.d0,p_var)
   call quad_min(ene_var,3)
  endif

! final value of diag_stab
  iwf = 1
  call wf_restore
  if (l_opt_pjas) call restore_pjas
  call object_restore ('gradient')
  if(l_opt_nwt) then
    call object_restore ('hess_nwt')
  endif
  if(l_opt_lin) then
     call object_restore ('ham_lin')
     call object_restore ('ovlp_lin')
     call object_restore ('renorm_vector')
  endif
  if(l_opt_ptb) then
     if (l_diagonal_overlap) then
      call object_restore ('dpsi_sq_covar')
     else
      call object_restore ('dpsi_dpsi_covar_inv')
     endif
!     call object_restore ('delta_e_ptb') !!
  endif  
  if(l_opt_ovlp_fn) then
      call object_restore ('dpsi_av')
      call object_restore ('dpsi_uwav')
      call object_restore ('dpsi_dpsi_uwcovar')
  endif
  diag_stab = min(add_diag (1), add_diag_max)
  if(diag_stab == add_diag_max) call die (lhere, 'diag_stab too large')
  call object_modified ('diag_stab')

  end subroutine adjust_diag_stab

!===========================================================================
  subroutine write_wf
!---------------------------------------------------------------------------
! Description : write optimized parameters of current wave function
!
! Created     : J. Toulouse, 23 Jan 2006
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'write_wf'
  integer orb_i, cent_i, dim_i, bas_i, bas_cent_i
  integer i, ict, isp
  character(len=80) fmt

! begin

! print CSFs coefficients
  if (l_opt_csf) then
    call object_provide ('ncsf')
    call object_provide ('csf_coef')
    write(6,'(a)') 'CSFs coefficients:'
    write(6,'(a)') 'csfs'
    write(6,'(a)') ' csf_coef'
# if defined (PATHSCALE)
    write(6,'(1000f15.8)') csf_coef(1:ncsf,iwf) ! for pathscale compiler
# else
    if (use_parser) then
     write(6,'(1000f15.8)') csf_coef(1:ncsf,iwf)
    else
     write(fmt,"( '(', i8, 'f15.8,a)' )") ncsf
     write(6,fmt) csf_coef(1:ncsf,iwf),' (csf_coef(icsf),icsf=1,ncsf)'
    endif
# endif
    write(6,'(a)') ' end'
    write(6,'(a)') 'end'
  endif ! l_opt_csf

! print Jastrow parameters
  if (l_opt_jas) then

    call object_provide ('nctype')
    call object_provide ('a4')
    call object_provide ('b')
    call object_provide ('c')

      write(6,'(a)') 'Jastrow parameters:'
      write(6,'(a)') 'jastrow'
      write(6,'(a)') ' parameters'
      if (nparma_read > 0) then
        write(fmt,'(''('',i2,''g20.12,a)'')') nparma_read
      else
        write(fmt,'(''(a)'')')
      endif
      do ict=1,nctype
        write(6,fmt) (a4(i,ict,iwf),i=1,nparma_read),' (a(iparmj),iparmj=1,nparma)'
      enddo

      if(nparmb_read > 0) then
        write(fmt,'(''('',i2,''g20.12,a)'')') nparmb_read
      else
        write(fmt,'(''(a)'')')
      endif
      do isp=nspin1,nspin2b
        write(6,fmt) (b(i,isp,iwf),i=1,nparmb_read),' (b(iparmj),iparmj=1,nparmb)'
      enddo

      if(nparmc_read > 0) then
        write(fmt,'(''('',i2,''g20.12,a)'')') nparmc_read
        do ict=1,nctype
          write(6,fmt) (c(i,ict,iwf),i=1,nparmc_read),' (c(iparmj),iparmj=1,nparmc)'
        enddo
      endif

      write(6,'(a)') ' end'
      write(6,'(a)') 'end'

  endif ! l_opt_jas

! print pjas parameters
  if (l_opt_pjas) then
     if ( l_opt_pjasen) then
        write(6,'(''periodic Jastrow parameters (pjas_en)'',1000f10.6)') pjas_parms (1:param_pjasen_nb, iwf)
        if (.not. inversion) then
           write(6,'(a,1000f10.6)') "cosine", (pjas_parms (i, iwf), i=1,param_pjasen_nb-1,2)
           write(6,'(a,1000f10.6)') "sine", (pjas_parms (i+1, iwf), i=1,param_pjasen_nb-1,2)
        endif
     endif
     if (l_opt_pjasee) then
        write(6,'(''periodic Jastrow parameters (pjas_ee)'',1000f10.6)') pjas_parms (param_pjasen_nb+1:param_pjas_nb, iwf)
     endif
  endif

! print orbitals coefficients
  if (l_opt_orb .or. (l_opt_exp .and. trim(basis_functions_varied) /= 'normalized')) then

    if (iperiodic == 0) then
      call object_provide ('nbasis')
      call object_provide ('orb_tot_nb')
      call object_provide ('coef_orb_on_norm_basis')

      write(6,'(a)') 'Orbital coefficients:'
      write(6,'(a)') 'orbitals'
      write(6,'(a)') ' coefficients'
      do orb_i = 1, orb_tot_nb
        if(orb_i==1) then
# if defined (PATHSCALE)
          write(6,'(1000es16.8)') coef_orb_on_norm_basis (1:nbasis, orb_i, iwf) ! for pathscale compiler
# else
          write(fmt,"( '(', i8, 'es16.8,a)' )") nbasis
          write(6,fmt) coef_orb_on_norm_basis(1:nbasis, orb_i, iwf), ' (coef(i,j),j=1,nbasis)'
# endif
        else
          write(fmt,"( '(', i8, 'es16.8)' )") nbasis
          write(6,fmt) coef_orb_on_norm_basis(1:nbasis, orb_i, iwf)
        endif ! orb_i==1
      enddo ! orb_i

    else
      call write_orbitals_pw_real
      write(6,'(3a)') 'Orbital coefficients written in file >',trim(file_orbitals_pw_out),'<'
    endif ! iperiodic

    write(6,'(a)') ' end'
    write(6,'(a)') 'end'
  endif ! l_opt_orb

! print basis exponents
  if (l_opt_exp) then
    call object_provide ('nbasis')
    call object_provide ('zex')
    if (use_parser) then
      write(6,'(a)') 'basis'
      write(6,'(a)') ' basis_functions'
      bas_i=0
      do ict=1,nctype
        write(6,'(i3)') ict
        do bas_cent_i = 1, basis_fns_by_center_type_nb (ict)
          bas_i = bas_i + 1
          select case (trim(radial_basis_type))
          case ('analytical')
            write(6,'(a6,3x,f12.8)') basis_fns_type_by_center_type (ict)%row(bas_cent_i), zex (bas_i, iwf)
          case ('numerical')
           call die (lhere, 'numerical basis functions cannot be optimized.')
          case ('numerical-analytical')
            write(6,'(a6,3x,f12.8,3x,i3)') basis_fns_type_by_center_type (ict)%row(bas_cent_i), zex (bas_i, iwf), basis_fns_rad_by_center_type(ict)%row(bas_cent_i)
          case default
           call die (lhere, 'unknown keyword >'+trim(radial_basis_type)+'<')
          end select
        enddo
      enddo
      write(6,'(a)') ' end'
      write(6,'(a)') 'end'
    else
      write(6,'(a)') 'Basis exponents:'
# if defined (PATHSCALE)
      write(6,'(100f10.6)') zex (1:nbasis, iwf) ! for pathscale compiler
# else
      write(fmt,"( '(', i8, 'f10.6,a)' )") nbasis
      write(6,fmt) zex(1:nbasis, iwf), ' (zex(i),i=1,nbasis)'
# endif
    endif ! use_parser
  endif ! l_opt_exp

! print geometry
  if (l_opt_geo) then
   call object_provide ('ncent')
   call object_provide ('ndim')
   call object_provide ('cent2')
   write(6,'(a)') 'Geometry:'
   do cent_i = 1, ncent
    write(6,'(a,i3,a,3f9.5)') 'center # ',cent_i, ' :',(cent2 (dim_i, cent_i,iwf), dim_i = 1, ndim)
   enddo ! cent_i
  endif ! l_opt_geo

  write(6,*)

  end subroutine write_wf

!===========================================================================
  subroutine write_wf_new
!---------------------------------------------------------------------------
! Description : write optimized parameters of new wave function (at each optimization step)
! Description : at each optimization step
!
! Created     : J. Toulouse, 04 Dec 2007
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'write_wf_new'
  integer orb_i, cent_i, dim_i, bas_i, bas_cent_i
  integer i, ict, isp
  character(len=80) fmt

! begin

! print CSFs coefficients
  if (l_opt_csf) then
    call object_provide ('ncsf')
    call object_provide ('csf_coef')
    write(6,'(a)') 'CSFs coefficients:'
    write(6,'(a)') 'csfs'
    write(6,'(a)') ' csf_coef'
# if defined (PATHSCALE)
    write(6,'(1000f15.8)') csf_coef(1:ncsf,iwf) ! for pathscale compiler
# else
    if (use_parser) then
     write(6,'(1000f15.8)') csf_coef(1:ncsf,iwf)
    else
     write(fmt,"( '(', i8, 'f15.8,a)' )") ncsf
     write(6,fmt) csf_coef(1:ncsf,iwf), ' (csf_coef_new(icsf),icsf=1,ncsf)'
    endif
# endif
    write(6,'(a)') ' end'
    write(6,'(a)') 'end'
  endif ! l_opt_csf

! print Jastrow parameters
  if (l_opt_jas) then

    call object_provide ('nctype')
    call object_provide ('a4')
    call object_provide ('b')
    call object_provide ('c')

      write(6,'(a)') 'Jastrow parameters:'
      write(6,'(a)') 'jastrow'
      write(6,'(a)') ' parameters'
      if (nparma_read > 0) then
        write(fmt,'(''(1p,'',i2,''g20.12,a)'')') nparma_read
       else
        write(fmt,'(''(a)'')')
      endif
      do ict=1,nctype
        write(6,fmt) (a4(i,ict,iwf),i=1,nparma_read),' (a_new(iparmj),iparmj=1,nparma)'
      enddo

      if(nparmb_read > 0) then
        write(fmt,'(''(1p,'',i2,''g20.12,a)'')') nparmb_read
       else
        write(fmt,'(''(a)'')')
      endif
      do isp=nspin1,nspin2b
        write(6,fmt) (b(i,isp,iwf),i=1,nparmb_read),' (b_new(iparmj),iparmj=1,nparmb)'
      enddo

      if(nparmc_read > 0) then
        write(fmt,'(''(1p,'',i2,''g20.12,a)'')') nparmc_read
       else
        write(fmt,'(''(a)'')')
      endif
      do ict=1,nctype
        write(6,fmt) (c(i,ict,iwf),i=1,nparmc_read),' (c_new(iparmj),iparmj=1,nparmc)'
      enddo

      write(6,'(a)') ' end'
      write(6,'(a)') 'end'
  endif ! l_opt_jas

! print orbitals coefficients
  if (l_opt_orb .or. (l_opt_exp .and. trim(basis_functions_varied) /= 'normalized')) then

    if (iperiodic == 0) then
      call object_provide ('nbasis')
      call object_provide ('orb_tot_nb')
      call object_provide ('coef_orb_on_norm_basis')

      write(6,'(a)') 'Orbital coefficients:'
      write(6,'(a)') 'orbitals'
      write(6,'(a)') ' coefficients'
      do orb_i = 1, orb_tot_nb
        if(orb_i==1) then
# if defined (PATHSCALE)
          write(6,'(1000es16.8)') coef_orb_on_norm_basis (1:nbasis, orb_i, iwf) ! for pathscale compiler
# else
          write(fmt,"( '(', i8, 'es16.8,a)' )") nbasis
          write(6,fmt) coef_orb_on_norm_basis(1:nbasis, orb_i, iwf), ' (coef_new(i,j),j=1,nbasis)'
# endif
        else
          write(fmt,"( '(', i8, 'es16.8)' )") nbasis
          write(6,fmt) coef_orb_on_norm_basis(1:nbasis, orb_i, iwf)
        endif ! iorb_i==1
      enddo ! orb_i

    else
      call write_orbitals_pw_real
      write(6,'(3a)') 'Orbital coefficients written in file >',trim(file_orbitals_pw_out),'<'
    endif ! iperiodic

    write(6,'(a)') ' end'
    write(6,'(a)') 'end'

  endif ! l_opt_orb

! print basis exponents
  if (l_opt_exp) then
    call object_provide ('nbasis')
    call object_provide ('zex')
    if (use_parser) then
      write(6,'(a)') 'basis'
      write(6,'(a)') ' basis_functions'
      bas_i=0
      do ict=1,nctype
        write(6,'(i3)') ict
        do bas_cent_i = 1, basis_fns_by_center_type_nb (ict)
          bas_i = bas_i + 1
          select case (trim(radial_basis_type))
          case ('analytical')
            write(6,'(a6,3x,f12.8)') basis_fns_type_by_center_type (ict)%row(bas_cent_i), zex (bas_i, iwf)
          case ('numerical')
           call die (lhere, 'numerical basis functions cannot be optimized.')
          case ('numerical-analytical')
            write(6,'(a6,3x,f12.8,3x,i3)') basis_fns_type_by_center_type (ict)%row(bas_cent_i), zex (bas_i, iwf), basis_fns_rad_by_center_type(ict)%row(bas_cent_i)
          case default
           call die (lhere, 'unknown keyword >'+trim(radial_basis_type)+'<')
          end select
        enddo
      enddo
      write(6,'(a)') ' end'
      write(6,'(a)') 'end'
    else
      write(6,'(a)') 'Basis exponents:'
# if defined (PATHSCALE)
      write(6,'(100f10.6)') zex (1:nbasis, iwf) ! for pathscale compiler
# else
      write(fmt,"( '(', i8, 'f12.8,a)' )") nbasis
      write(6,fmt) zex(1:nbasis, iwf), ' (zex_new(i),i=1,nbasis)'
# endif
    endif ! use_parser
  endif ! l_opt_exp

! print geometry
  if (l_opt_geo) then
    call object_provide ('ncent')
    call object_provide ('ndim')
    call object_provide ('cent2')
    write(6,'(a)') 'Geometry:'
    do cent_i = 1, ncent
      write(6,'(a,i3,a,3f9.5)') 'center # ',cent_i, ' :',(cent2 (dim_i, cent_i,iwf), dim_i = 1, ndim)
    enddo ! cent_i
  endif ! l_opt_geo

  write(6,*)
  call systemflush(6)

  end subroutine write_wf_new

!===========================================================================
  subroutine write_wf_best
!---------------------------------------------------------------------------
! Description : write parameters of best wave function
!
! Created     : J. Toulouse, 20 Nov 2007
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'write_wf_best'
  integer orb_i, cent_i, dim_i, bas_i, bas_cent_i
  integer i, ict, isp
  character(len=80) fmt

! begin

! print CSFs coefficients
  if (l_opt_csf) then
    call object_provide ('ncsf')
    call object_provide ('csf_coef_best')
    write(6,'(a)') 'CSFs coefficients:'
    write(6,'(a)') 'csfs'
    write(6,'(a)') ' csf_coef'
# if defined (PATHSCALE)
    write(6,'(1000f15.8)') csf_coef_best(1:ncsf) ! for pathscale compiler
# else
    write(fmt,"( '(', i8, 'f15.8,a)' )") ncsf
    write(6,fmt) csf_coef_best(1:ncsf), ' (csf_coef_best(icsf),icsf=1,ncsf)'
# endif
    write(6,'(a)') ' end'
    write(6,'(a)') 'end'
  endif ! l_opt_csf

! print Jastrow parameters
  if (l_opt_jas) then

    call object_provide ('nctype')
    call object_provide ('a4_best')
    call object_provide ('b_best')
    call object_provide ('c_best')

    write(6,'(a)') 'Jastrow parameters:'
    write(6,'(a)') 'jastrow'
    write(6,'(a)') ' parameters'
    if (nparma_read > 0) then
      write(fmt,'(''(1p,'',i2,''g20.12,a)'')') nparma_read
    else
      write(fmt,'(''(a)'')')
    endif
    do ict=1,nctype
      write(6,fmt) (a4_best(i,ict),i=1,nparma_read),' (a_best(iparmj),iparmj=1,nparma)'
    enddo

    if(nparmb_read > 0) then
      write(fmt,'(''(1p,'',i2,''g20.12,a)'')') nparmb_read
    else
      write(fmt,'(''(a)'')')
    endif
    do isp=nspin1,nspin2b
      write(6,fmt) (b_best(i,isp),i=1,nparmb_read),' (b_best(iparmj),iparmj=1,nparmb)'
    enddo

    if(nparmc_read > 0) then
      write(fmt,'(''(1p,'',i2,''g20.12,a)'')') nparmc_read
    else
      write(fmt,'(''(a)'')')
    endif
    do ict=1,nctype
      write(6,fmt) (c_best(i,ict),i=1,nparmc_read),' (c_best(iparmj),iparmj=1,nparmc)'
    enddo

    write(6,'(a)') ' end'
    write(6,'(a)') 'end'
  endif ! l_opt_jas

! print orbitals coefficients
  if (l_opt_orb .or. (l_opt_exp .and. trim(basis_functions_varied) /= 'normalized')) then

    if (iperiodic == 0) then
    call object_provide ('nbasis')
    call object_provide ('orb_tot_nb')
    call object_provide ('coef_orb_on_norm_basis_best')

    write(6,'(a)') 'Orbital coefficients:'
    write(6,'(a)') 'orbitals'
    write(6,'(a)') ' coefficients'
    do orb_i = 1, orb_tot_nb
      if(orb_i==1) then
# if defined (PATHSCALE)
        write(6,'(1000es16.8)') coef_orb_on_norm_basis_best (1:nbasis, orb_i) ! for pathscale compiler
# else
        write(fmt,"( '(', i8, 'es16.8,a)' )") nbasis
        write(6,fmt) coef_orb_on_norm_basis_best(1:nbasis, orb_i), ' (coef_best(i,j),j=1,nbasis)'
# endif
      else
        write(fmt,"( '(', i8, 'es16.8)' )") nbasis
        write(6,fmt) coef_orb_on_norm_basis_best(1:nbasis, orb_i)
      endif ! orb_i==1
    enddo ! orb_i

    else
      call write_orbitals_pw_real
      write(6,'(3a)') 'Best orbital coefficients written in file >',trim(file_orbitals_pw_out),'<'
    endif ! iperiodic

    write(6,'(a)') ' end'
    write(6,'(a)') 'end'

  endif ! l_opt_orb

! print basis exponents
  if (l_opt_exp) then
    call object_provide ('nbasis')
    call object_provide ('zex_best')
    if (use_parser) then
      write(6,'(a)') 'basis'
      write(6,'(a)') ' basis_functions'
      bas_i=0
      do ict=1,nctype
        write(6,'(i3)') ict
        do bas_cent_i = 1, basis_fns_by_center_type_nb (ict)
          bas_i = bas_i + 1
          select case (trim(radial_basis_type))
          case ('analytical')
            write(6,'(a6,3x,f12.8)') basis_fns_type_by_center_type (ict)%row(bas_cent_i), zex_best (bas_i)
          case ('numerical')
           call die (lhere, 'numerical basis functions cannot be optimized.')
          case ('numerical-analytical')
            write(6,'(a6,3x,f12.8,3x,i3)') basis_fns_type_by_center_type (ict)%row(bas_cent_i), zex_best (bas_i), basis_fns_rad_by_center_type(ict)%row(bas_cent_i)
          case default
           call die (lhere, 'unknown keyword >'+trim(radial_basis_type)+'<')
          end select
        enddo
      enddo
      write(6,'(a)') ' end'
      write(6,'(a)') 'end'
    else
      write(6,'(a)') 'Basis exponents:'
# if defined (PATHSCALE)
      write(6,'(1000f10.6)') zex_best (1:nbasis) ! for pathscale compiler
# else
      write(fmt,"( '(', i8, 'f12.8,a)' )") nbasis
      write(6,fmt) zex_best(1:nbasis), ' (zex_best(i),i=1,nbasis)'
# endif
    endif ! use_parser
  endif ! l_opt_exp

! print geometry
  if (l_opt_geo) then
    call object_provide ('ncent')
    call object_provide ('ndim')
    call object_provide ('cent_best')
    write(6,'(a)') 'Geometry:'
    do cent_i = 1, ncent
      write(6,'(a,i3,a,3f9.5)') 'center # ',cent_i, ' :',(cent_best (dim_i, cent_i), dim_i = 1, ndim)
    enddo ! cent_i
  endif ! l_opt_geo

  write(6,*)

  end subroutine write_wf_best

! ==============================================================================
  subroutine delta_param_bld
! ------------------------------------------------------------------------------
! Description   : parameter changes after one optimization step
!
! Created       : J. Toulouse, 17 Oct 2005
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len), save :: lhere = 'delta_param_bld'
  integer shift


! header
  if (header_exe) then

   call object_create ('delta_param')
   call object_create ('delta_param_norm')
   call object_create ('delta_csf')
   call object_create ('delta_csf_norm')
   call object_create ('delta_jas')
   call object_create ('delta_jas_norm')
   call object_create ('delta_pjas')
   call object_create ('delta_coef_ex')
   call object_create ('delta_exp')
   call object_create ('delta_geo')

   call object_needed ('param_nb')

   return

  endif

! begin

! allocations
  call object_alloc ('delta_param', delta_param, param_nb)

! parameter variations from one optimization method
  if (l_opt_nwt) then
      call object_provide_in_node (lhere, 'delta_nwt')
      delta_param (:) = delta_nwt (:)
  elseif (l_opt_lin) then
      call object_provide_in_node (lhere, 'delta_lin')
      delta_param (:) = delta_lin (:)
  elseif (l_opt_ptb) then
      call object_provide_in_node (lhere, 'delta_ptb')
      delta_param (:) = delta_ptb (:)
  elseif (l_opt_ovlp_fn_linear) then
      call object_provide_in_node (lhere, 'delta_ovlp_fn_linear')
      delta_param (:) = delta_ovlp_fn_linear (:)
  elseif (l_opt_ovlp_fn_nonsemiortho) then
      call object_provide_in_node (lhere, 'delta_ovlp_fn_nonsemiortho')
      delta_param (:) = delta_ovlp_fn_nonsemiortho (:)
  elseif (l_opt_ovlp_fn) then
      call object_provide_in_node (lhere, 'delta_ovlp_fn')
      delta_param (:) = delta_ovlp_fn (:)
  else
      call die (lhere, 'No parameter variations available.')
  endif

! compute the reduced norm of delta_param
  delta_param_norm = dsqrt(dot_product(delta_param,delta_param)/param_nb)

  write(6,*)

! split parameter variations
  shift = 0
  if (l_opt_csf) then
   call object_alloc ('delta_csf', delta_csf, nparmcsf)
   delta_csf (:) = delta_param (shift+1: nparmcsf)
   shift = shift + nparmcsf
   write(6,'(a,500f12.7)') 'CSF parameters variations=', delta_csf (:)
   delta_csf_norm = dsqrt(dot_product(delta_csf,delta_csf)/nparmcsf)
  endif

  if (l_opt_jas) then
   call object_alloc ('delta_jas', delta_jas, nparmj)
   delta_jas (:)= delta_param (shift+1: shift+nparmj)
   shift = shift + nparmj
   write(6,'(a,500f12.7)') 'Jastrow parameters variations=', delta_jas (:)
   delta_jas_norm = dsqrt(dot_product(delta_jas,delta_jas)/nparmj)
  endif

  if (l_opt_pjas) then
     call object_alloc ('delta_pjas', delta_pjas, param_pjas_nb)
     delta_pjas (:)= delta_param (shift+1: shift+param_pjas_nb)
     shift = shift + param_pjas_nb
     write(6,'(a,600f12.7)') 'Periodic Jastrow parameters variations=', delta_pjas (:)
  endif

  if (l_opt_exp) then
   call object_alloc ('delta_exp', delta_exp, param_exp_nb)
   delta_exp (:)= delta_param (shift+1: shift+param_exp_nb)
   shift = shift + param_exp_nb
   write(6,'(a,600f12.7)') 'Exponent parameters variations=', delta_exp (:)
  endif

  if (l_opt_orb) then
   call object_alloc ('delta_coef_ex', delta_coef_ex, param_orb_nb)
   delta_coef_ex (:)= delta_param (shift+1: shift+param_orb_nb)
   shift = shift + param_orb_nb
   write(6,'(a,600f12.7)') 'Orbital rotations parameters variations=', delta_coef_ex (:)
  endif

  if (l_opt_geo) then
   call object_alloc ('delta_geo', delta_geo, param_geo_nb)
   delta_geo (:)= delta_param (shift+1: shift+param_geo_nb)
   shift = shift + param_geo_nb
   write(6,'(a,600f12.7)') 'Geometry parameters variations=', delta_geo (:)
  endif

 end subroutine delta_param_bld

! ==============================================================================
  subroutine delta_csf_rot_bld
! ------------------------------------------------------------------------------
! Description   : CSF parameter changes by rotation (without Jastrow)
!
! Created       : J. Toulouse, 29 Mar 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer i
  real(dp) d, cosd, sindoverd

! header
  if (header_exe) then

   call object_create ('delta_csf_rot')

   call object_needed ('nparmcsf')
   call object_needed ('ncsf')
   call object_needed ('iwcsf')
   call object_needed ('delta_csf')
   call object_needed ('csf_coef')

   return

  endif

! begin
! allocations
  call object_alloc ('delta_csf_rot', delta_csf_rot, ncsf)

  d = 0.d0
  do i = 1, nparmcsf
      d = d  + delta_csf (i) * delta_csf (i)
  enddo

  d = dsqrt(d)
  cosd = dcos(d)
  sindoverd = dsin(d)/d

  delta_csf_rot(1) = (cosd - 1.d0) * csf_coef(1,1)

  do i = 1, nparmcsf
     delta_csf_rot(i+1)= (cosd - 1.d0) * csf_coef(i+1,1) + sindoverd * delta_csf (i)
  enddo

 end subroutine delta_csf_rot_bld

! ==============================================================================
  subroutine delta_csf_rot_bld_old2
! ------------------------------------------------------------------------------
! Description   : CSF parameter changes by rotation (without Jastrow)
!
! Created       : J. Toulouse, 29 Mar 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer iparmcsf, jparmcsf
  real(dp) d, sum_temp, sindoverd

! header
  if (header_exe) then

   call object_create ('delta_csf_rot')

   call object_needed ('nparmcsf')
   call object_needed ('ncsf')
   call object_needed ('iwcsf')
   call object_needed ('delta_csf')
   call object_needed ('csf_coef')
   call object_needed ('wfdet_ovlp') !!
   call object_needed ('csfs_wfdet_ovlp') !!
   call object_needed ('csfs_ovlp') !!!

   return

  endif

! begin
! allocations
  call object_alloc ('delta_csf_rot', delta_csf_rot, nparmcsf)

  d = 0.d0
  do iparmcsf = 1, nparmcsf
   do jparmcsf = 1, nparmcsf
      d = d  + delta_csf (iparmcsf) * delta_csf (jparmcsf) * (csfs_ovlp(iwcsf(iparmcsf),iwcsf(jparmcsf)) - csfs_wfdet_ovlp(iwcsf(iparmcsf)) * csfs_wfdet_ovlp(iwcsf(jparmcsf)) )/ wfdet_ovlp
   enddo
  enddo
  d = dsqrt(d)

  sum_temp = 0.d0
  do iparmcsf = 1, nparmcsf
     sum_temp = sum_temp + delta_csf (iparmcsf) * csfs_wfdet_ovlp(iwcsf(iparmcsf))/wfdet_ovlp
  enddo

  sindoverd = dsin(d)/d

  do iparmcsf=1,nparmcsf
     delta_csf_rot(iparmcsf)= (dcos(d) - sindoverd * sum_temp - 1.d0) * csf_coef(iwcsf(iparmcsf),1)  &
                                 + sindoverd * delta_csf (iparmcsf)
  enddo


 end subroutine delta_csf_rot_bld_old2

! ==============================================================================
  subroutine delta_csf_rot_bld_old
! ------------------------------------------------------------------------------
! Description   : CSF parameter changes by rotation (incorporating the Jastrow)
!
! Created       : J. Toulouse, 27 Mar 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer iparmcsf, jparmcsf
  real(dp) d, sum_temp, sindoverd

! header
  if (header_exe) then

   call object_create ('delta_csf_rot')

   call object_needed ('nparmcsf')
   call object_needed ('ncsf')
   call object_needed ('iwcsf')
   call object_needed ('delta_csf')
   call object_needed ('csf_coef')
   call object_needed ('dpsi_csf_av')
   call object_needed ('dpsi_csf_dpsi_csf_covar')

   return

  endif

! begin
! allocations
  call object_alloc ('delta_csf_rot', delta_csf_rot, nparmcsf)

  d = 0.d0
  do iparmcsf = 1, nparmcsf
   do jparmcsf = 1, nparmcsf
      d = d  + delta_csf (iparmcsf) * delta_csf (jparmcsf) * dpsi_csf_dpsi_csf_covar (iparmcsf, jparmcsf)
   enddo
  enddo
  d = dsqrt(d)

  sum_temp = 0.d0
  do iparmcsf = 1, nparmcsf
     sum_temp = sum_temp + delta_csf (iparmcsf) * dpsi_csf_av (iparmcsf)
  enddo

  sindoverd = dsin(d)/d

  do iparmcsf=1,nparmcsf
     delta_csf_rot(iparmcsf)= (dcos(d) - sindoverd * sum_temp - 1.d0) * csf_coef(iwcsf(iparmcsf),1)  &
                                 + sindoverd * delta_csf (iparmcsf)
  enddo

 end subroutine delta_csf_rot_bld_old

! ==============================================================================
  subroutine delta_mat_rot_1st_order_bld
! ------------------------------------------------------------------------------
! Description   : approximate first-order orbital rotation matrix = kappa matrix
!
! Created       : J. Toulouse, 02 May 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer ex_i, dorb_i, orb_1st, orb_2nd

! header
  if (header_exe) then

   call object_create ('delta_mat_rot_1st_order')

   call object_needed ('orb_tot_nb')
   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_1st_lab')
   call object_needed ('ex_orb_2nd_lab')
   call object_needed ('delta_coef_ex')

   return

  endif

! set up kappa matrix
! this matrix is real anti-symmetric if all orthogonality constraints kij=-kji are imposed
! for active-active excitations the orthogonality constraints may or may not be imposed
  call object_alloc ('delta_mat_rot_1st_order', delta_mat_rot_1st_order, orb_tot_nb, orb_tot_nb)
  delta_mat_rot_1st_order (:,:) =  0.d0
  do dorb_i = 1, param_orb_nb
    ex_i = ex_orb_ind (dorb_i)
    orb_1st = ex_orb_1st_lab (ex_i)
    orb_2nd = ex_orb_2nd_lab (ex_i)
    delta_mat_rot_1st_order (orb_1st, orb_2nd) = delta_coef_ex (dorb_i)
!   anti-symmetrize elements that are zero (but do not necessarily anti-symmetrize active-active block for which 
!   we can have independent kij and kji if orthogonality constraint is not imposed
    if (delta_mat_rot_1st_order (orb_2nd, orb_1st) == 0.d0) then
      delta_mat_rot_1st_order (orb_2nd, orb_1st) = - delta_mat_rot_1st_order (orb_1st, orb_2nd)
    endif
  enddo ! dorb_i

  end subroutine delta_mat_rot_1st_order_bld

! ==============================================================================
  subroutine delta_mat_rot_real_bld
! ------------------------------------------------------------------------------
! Description   : orbital rotation matrix
!
! Created       : J. Toulouse, 07 Feb 2006
! Modified      : J. Toulouse, 01 May 2010: begin case of no orthogonality constraint for active orbitals
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_mat_rot_real_bld'
  integer ex_i, i, j, k, dorb_i
  integer orb_1st, orb_2nd
  integer lwork
  integer info
  integer eig
  real(dp), allocatable :: kappa (:,:)
  real(dp), allocatable :: mat_rot_real (:,:)
  real(dp), allocatable :: mat_a (:,:)
  real(dp), allocatable :: mat_vr (:,:)
  real(dp), allocatable :: mat_vl (:,:)
  real(dp), allocatable :: mat_wr (:)
  real(dp), allocatable :: mat_wi (:)
  real(dp), allocatable :: work (:)
  real(dp), allocatable :: identity (:,:)
  double complex, allocatable :: mat_w (:)
!  double complex, allocatable :: identity_complex (:,:)
  double complex, allocatable :: mat_u (:,:), mat_u_inv (:,:)
  double complex, allocatable :: kappa_check (:,:)
  double complex, allocatable :: mat_exp_w (:)
  double complex, allocatable :: delta_mat_rot (:,:)
  double complex, allocatable :: mat_rot (:,:)


! header
  if (header_exe) then

   call object_create ('delta_mat_rot_real')

   call object_needed ('orb_tot_nb')
   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_1st_lab')
   call object_needed ('ex_orb_2nd_lab')
   call object_needed ('delta_coef_ex')

   return

  endif

! set up kappa matrix
! this matrix is real anti-symmetric if all orthogonality constraints kij=-kji are imposed
! for active-active excitations the orthogonality constraints may or may not be imposed
  call alloc ('kappa', kappa, orb_tot_nb, orb_tot_nb)
  kappa (:,:) =  0.d0
  do dorb_i = 1, param_orb_nb
    ex_i = ex_orb_ind (dorb_i)
    orb_1st = ex_orb_1st_lab (ex_i)
    orb_2nd = ex_orb_2nd_lab (ex_i)
    kappa (orb_1st, orb_2nd) = delta_coef_ex (dorb_i)
!   anti-symmetrize elements that are zero (but do not necessarily anti-symmetrize active-active block for which 
!   we can have independent kij and kji if orthogonality constraint is not imposed
    if (kappa (orb_2nd, orb_1st) == 0.d0) then
      kappa (orb_2nd, orb_1st) = - kappa (orb_1st, orb_2nd)
    endif
  enddo ! dorb_i

! diagonalize kappa
  call alloc ('mat_a', mat_a, orb_tot_nb, orb_tot_nb)
  call alloc ('mat_vr', mat_vr, orb_tot_nb, orb_tot_nb)
  call alloc ('mat_vl', mat_vl, orb_tot_nb, orb_tot_nb)
  call alloc ('mat_wr', mat_wr, orb_tot_nb)
  call alloc ('mat_wi', mat_wi, orb_tot_nb)
  lwork = 1
  call alloc ('work', work, lwork)
  call dgeev('V','V', orb_tot_nb, mat_a, orb_tot_nb, mat_wr, mat_wi, mat_vl, orb_tot_nb, mat_vr, orb_tot_nb, work, -1, info)
  if (info /= 0) then
   call die (here, 'problem in dgeev')
  endif
  lwork = nint(work(1))
  call alloc ('work', work, lwork)
  mat_a (:,:) = kappa (:,:)
  call dgeev('V','V', orb_tot_nb, mat_a, orb_tot_nb, mat_wr, mat_wi, mat_vl, orb_tot_nb, mat_vr, orb_tot_nb, work, lwork, info)
  call release ('mat_a', mat_a)
  call release ('work', work)
  if (info /= 0) then
   call die (here, 'problem in dgeev')
  endif

! complex eigenvalues
  call alloc ('mat_w', mat_w, orb_tot_nb)
  mat_w (:) = dcmplx(mat_wr(:),mat_wi(:))
!  write(6,'(2a,f,a,f,a)') trim(here),': (complex) eigenvalues of kappa:'
!  do i = 1, orb_tot_nb
!   write(6,'(2a,i3,a,e,a,e,a)') trim(here),': ',i,' eig=',mat_wr(i), ' +', mat_wi(i),' i'
!  enddo

! check that all eigenvalues are purely imaginary if orthogonality constraint is imposed
  if (l_active_orb_ortho_constraint) then
   do i = 1, orb_tot_nb
    if (dabs(mat_wr(i)) > 1.d-8) then
      write(6,'(2a,i3,a,es15.8,a,es15.8,a)') trim(here),': eigenvalue # ',i,' : ' ,mat_wr(i), ' +', mat_wi(i),' i  is not purely imaginary'
      call die (here)
    endif
   enddo
  endif

! set up matrix of eigenvectors. It is unitary if orthogonality constraint is imposed.
  call alloc ('mat_u', mat_u, orb_tot_nb, orb_tot_nb)
  eig = 0
  do
   eig = eig + 1
   if (eig > orb_tot_nb) exit
   mat_u (:,eig) = dcmplx(mat_vr(:,eig), 0.d0)
   if (eig+1 > orb_tot_nb) exit
   if (mat_w (eig) == dconjg(mat_w(eig+1))) then
!     write(6,*) trim(here),': eigenvalue # ',eig,'and # ',eig+1,' are congugate'
     mat_u (:,eig)   = dcmplx(mat_vr(:,eig), mat_vr(:,eig+1))
     mat_u (:,eig+1) = dcmplx(mat_vr(:,eig), -mat_vr(:,eig+1))
     eig = eig + 1
   endif
   cycle
  enddo

! compute inverse of eigenvector matrix if not the matrix is not unitary
! unfinished...
  if (.not. l_active_orb_ortho_constraint) then
   call die (lhere, 'routine not yet finished for active_orb_ortho_constraint = false')
   call alloc ('mat_u_inv', mat_u_inv, orb_tot_nb, orb_tot_nb)
   call inverse_by_svd (mat_u, mat_u_inv, orb_tot_nb, 1.d-8)
  endif

!  write(6,'(2a)') trim(here), 'unitary eigenvectors:'
!  do j = 1, orb_tot_nb
!    write(6,*) trim(here),': eigenvector # ',j, (mat_u (i,j),i=1,orb_tot_nb)
!  enddo
!

! check that U*U^dag = 1: this is not true numerically in the subspace of zero eigenvalue
!  call alloc ('identity_complex', identity_complex, orb_tot_nb, orb_tot_nb)
!  identity_complex (:,:) = matmul(mat_u,transpose(dconjg(mat_u)))
!
!  do i = 1,orb_tot_nb
!   do j = 1, orb_tot_nb
!     write(6,'(2a,i3,a,i3,a,f,f)') trim(here),': identity(',i,',',j,')=',identity_complex(i,j)
!   enddo
!  enddo

! check the recovering of kappa
  call alloc ('kappa_check', kappa_check, orb_tot_nb, orb_tot_nb)
  kappa_check (:,:) = 0.d0
  do i = 1,orb_tot_nb
   do j = 1,orb_tot_nb
    do k = 1, orb_tot_nb
      kappa_check (i,j) = kappa_check (i,j) +  mat_u (i, k) * mat_w (k) * dconjg(mat_u (j,k))
    enddo
!   check absolute accuracy
    if (dabs(real(kappa_check(i,j))-kappa(i,j)) > 10d-10 .or. aimag(kappa_check(i,j)) > 10d-10) then
!   check relative accuracy
     if (dabs(1.d0 - real(kappa_check(i,j))/kappa(i,j)) > 10d-10) then
         write(6,'(2a,i3,a,i3,a,es15.8,es15.8)') trim(here),': kappa(',i,',',j,')=',kappa(i,j)
         write(6,'(2a,i3,a,i3,a,es15.8,es15.8)') trim(here),': kappa_check(',i,',',j,')=',kappa_check(i,j)
         write (6,'(2a)') trim(here),': problem in the diagonalization of orbital rotation generator matrix kappa.'
         write (6,'(2a)') trim(here),': most likely, this happened because there are redundant orbital parameters.'
         call die (here)
     endif
    endif
   enddo
  enddo

! set up orbital (real orthogonal) rotation matrix:  X = exp (kappa) = U exp (eigenvalues) U^dag
  call alloc ('mat_exp_w', mat_exp_w, orb_tot_nb)
  call alloc ('delta_mat_rot', delta_mat_rot, orb_tot_nb, orb_tot_nb)
  mat_exp_w (:) = exp (mat_w (:))
  delta_mat_rot (:,:) = 0.d0
  do i = 1,orb_tot_nb
   do j = 1,orb_tot_nb
    do k = 1, orb_tot_nb
      delta_mat_rot (i,j) = delta_mat_rot (i,j) +  mat_u (i, k) * (mat_exp_w (k) - 1.d0) * dconjg(mat_u (j,k))
    enddo
   enddo
  enddo

  call alloc ('mat_rot', mat_rot, orb_tot_nb, orb_tot_nb)
  mat_rot (:,:) = delta_mat_rot(:,:)
  do i = 1, orb_tot_nb
   mat_rot (i,i) = 1.d0 + delta_mat_rot (i,i)
  enddo

! check that rotation matrix is real
  do i = 1,orb_tot_nb
   do j = 1,orb_tot_nb
      if (dabs(aimag(mat_rot (i,j))) > 10.d-10) then
         write(6,'(2a,i3,a,i3,a,es15.8,es15.8)') trim(here),': mat_rot(',i,',',j,')=',mat_rot(i,j)
         call die (here, 'orbital rotation matrix is not real')
      endif
   enddo
  enddo

  call alloc ('mat_rot_real', mat_rot_real, orb_tot_nb, orb_tot_nb)
  mat_rot_real = real(mat_rot)

! check orthogonality of rotation matrix
  call alloc ('identity', identity, orb_tot_nb, orb_tot_nb)
  identity (:,:) = matmul(mat_rot_real,transpose(mat_rot_real))
  do i = 1,orb_tot_nb
   do j = 1,orb_tot_nb
     if (i == j) then
       if(dabs(identity(i,j) -1.d0) > 10.d-10) then
           write(6,'(2a,i3,a,i3,a,es15.8,a)') trim(here),': identity(',i,',',j,')=',identity(i,j),' /= 1'
           call die (here, 'orbital rotation matrix is not orthogonal')
       endif
     else
       if(dabs(identity(i,j)) > 10.d-10) then
           write(6,'(2a,i3,a,i3,a,es15.8,a)') trim(here),': identity(',i,',',j,')=',identity(i,j), ' /= 0'
           call die (here, 'orbital rotation matrix is not orthogonal')
       endif
     endif
   enddo
  enddo

  call release ('kappa', kappa)
  call release ('mat_vr', mat_vr)
  call release ('mat_vl', mat_vl)
  call release ('mat_wr', mat_wr)
  call release ('mat_wi', mat_wi)
  call release ('mat_w', mat_w)
  call release ('mat_u', mat_u)
  call release ('mat_exp_w', mat_exp_w)
  call release ('mat_rot', mat_rot)
  call release ('mat_rot_real', mat_rot_real)

  call object_alloc ('delta_mat_rot_real', delta_mat_rot_real, orb_tot_nb, orb_tot_nb)
  delta_mat_rot_real = real(delta_mat_rot)
  call release ('delta_mat_rot', delta_mat_rot)

  end subroutine delta_mat_rot_real_bld

! ==============================================================================
  subroutine update_coef_by_rot
! ------------------------------------------------------------------------------
! Description   : update orbital coefficient by rotation
!
! Created       : J. Toulouse, 31 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer orb_i, orb_j
  real(dp), allocatable :: coef_new (:,:)

! begin

! provide needed objects
  call object_provide ('nbasis')
  call object_provide ('orb_tot_nb')
  call object_provide ('coef')

  call alloc ('coef_new', coef_new, nbasis, orb_tot_nb)
  coef_new (1:nbasis, 1:orb_tot_nb) = coef (1:nbasis, 1:orb_tot_nb, 1)

! update orbital coefficients by first-order rotation
  if (l_approx_orb_rot) then
    call object_provide ('delta_mat_rot_1st_order')
    do orb_i = 1, orb_tot_nb
      do orb_j = 1, orb_tot_nb
        coef_new (1:nbasis, orb_i) = coef_new (1:nbasis, orb_i) + delta_mat_rot_1st_order (orb_i, orb_j) * coef (1:nbasis, orb_j, 1)
      enddo ! orb_j
     enddo ! orb_i

! update orbital coefficients by exact rotation
  else
    call object_provide ('delta_mat_rot_real')
    do orb_i = 1, orb_tot_nb
      do orb_j = 1, orb_tot_nb
        coef_new (1:nbasis, orb_i) = coef_new (1:nbasis, orb_i) + delta_mat_rot_real (orb_i, orb_j) * coef (1:nbasis, orb_j, 1)
      enddo ! orb_j
     enddo ! orb_i
  endif

  coef (1:nbasis, 1:orb_tot_nb, iwf) = coef_new (1:nbasis, 1:orb_tot_nb)
  call object_modified ('coef')

  end subroutine update_coef_by_rot

! ==============================================================================
  subroutine update_coef_orb_on_norm_basis_by_rot
! ------------------------------------------------------------------------------
! Description   : update orbital coefficient on normalized basis by rotation
!
! Created       : J. Toulouse, 31 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer orb_i, orb_j
  real(dp), allocatable :: coef_new (:,:)

! begin

! provide needed objects
  call object_provide ('nbasis')
  call object_provide ('orb_tot_nb')
  call object_provide ('coef_orb_on_norm_basis')

  call alloc ('coef_new', coef_new, nbasis, orb_tot_nb)
  coef_new (1:nbasis, 1:orb_tot_nb) = coef_orb_on_norm_basis (1:nbasis, 1:orb_tot_nb, 1)

! update orbital coefficients by first-order rotation
  if (l_approx_orb_rot) then
    call object_provide ('delta_mat_rot_1st_order')
    do orb_i = 1, orb_tot_nb
      do orb_j = 1, orb_tot_nb
        coef_new (1:nbasis, orb_i) = coef_new (1:nbasis, orb_i) + delta_mat_rot_1st_order (orb_i, orb_j) * coef_orb_on_norm_basis (1:nbasis, orb_j, 1)
      enddo ! orb_j
     enddo ! orb_i

! update orbital coefficients by exact rotation
  else
    call object_provide ('delta_mat_rot_real')
    do orb_i = 1, orb_tot_nb
      do orb_j = 1, orb_tot_nb
        coef_new (1:nbasis, orb_i) = coef_new (1:nbasis, orb_i) + delta_mat_rot_real (orb_i, orb_j) * coef_orb_on_norm_basis (1:nbasis, orb_j, 1)
      enddo ! orb_j
     enddo ! orb_i
   endif

  coef_orb_on_norm_basis (1:nbasis, 1:orb_tot_nb, iwf) = coef_new (1:nbasis, 1:orb_tot_nb)
  call object_modified ('coef_orb_on_norm_basis')

  end subroutine update_coef_orb_on_norm_basis_by_rot

! ==============================================================================
  subroutine update_coef_orb_on_ortho_basis_by_rot
! ------------------------------------------------------------------------------
! Description   : update orbital coefficient on orthonormalized basis by rotation
!
! Created       : J. Toulouse, 31 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer orb_i, orb_j
  real(dp), allocatable :: coef_new (:,:)

! begin

! provide needed objects
  call object_provide ('nbasis')
  call object_provide ('orb_tot_nb')
  call object_provide ('coef_orb_on_ortho_basis')

  call alloc ('coef_new', coef_new, nbasis, orb_tot_nb)
  coef_new (1:nbasis, 1:orb_tot_nb) = coef_orb_on_ortho_basis (1:nbasis, 1:orb_tot_nb, 1)

! update orbital coefficients by first-order rotation
  if (l_approx_orb_rot) then
    call object_provide ('delta_mat_rot_1st_order')
    do orb_i = 1, orb_tot_nb
      do orb_j = 1, orb_tot_nb
        coef_new (1:nbasis, orb_i) = coef_new (1:nbasis, orb_i) + delta_mat_rot_1st_order (orb_i, orb_j) * coef_orb_on_ortho_basis (1:nbasis, orb_j, 1)
      enddo ! orb_j
     enddo ! orb_i

! update orbital coefficients by exact rotation
  else
    call object_provide ('delta_mat_rot_real')
    do orb_i = 1, orb_tot_nb
      do orb_j = 1, orb_tot_nb
        coef_new (1:nbasis, orb_i) = coef_new (1:nbasis, orb_i) + delta_mat_rot_real (orb_i, orb_j) * coef_orb_on_ortho_basis (1:nbasis, orb_j, 1)
      enddo ! orb_j
     enddo ! orb_i
   endif

  coef_orb_on_ortho_basis (1:nbasis, 1:orb_tot_nb, iwf) = coef_new (1:nbasis, 1:orb_tot_nb)
  call object_modified ('coef_orb_on_ortho_basis')

  end subroutine update_coef_orb_on_ortho_basis_by_rot

! ==============================================================================
  subroutine delta_coef_pw_bld
! ------------------------------------------------------------------------------
! Description   : Variation of the orbital coefficients by approximate rotation
! Description   : for plane wave basis
!
! Created       : J. Toulouse, 26 Mar 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer ex_i, orb_i
  integer orb_1st, orb_2nd


! header
  if (header_exe) then

   call object_create ('delta_c_rp')
   call object_create ('delta_c_rm')
   call object_create ('delta_c_ip')
   call object_create ('delta_c_im')

   call object_needed ('orb_tot_nb')
   call object_needed ('ngvec_orb')
   call object_needed ('param_orb_nb')
   call object_needed ('ex_orb_ind')
   call object_needed ('ex_orb_1st_lab')
   call object_needed ('ex_orb_2nd_lab')
   call object_needed ('delta_coef_ex')
   call object_needed ('c_rp')
   call object_needed ('c_rm')
   call object_needed ('c_ip')
   call object_needed ('c_im')

   return

  endif

! begin

! allocations
  call object_alloc ('delta_c_rp', delta_c_rp, ngvec_orb, orb_tot_nb)
  call object_alloc ('delta_c_rm', delta_c_rm, ngvec_orb, orb_tot_nb)
  call object_alloc ('delta_c_ip', delta_c_ip, ngvec_orb, orb_tot_nb)
  call object_alloc ('delta_c_im', delta_c_im, ngvec_orb, orb_tot_nb)

  delta_c_rp (:,:) = 0.d0
  delta_c_rm (:,:) = 0.d0
  delta_c_ip (:,:) = 0.d0
  delta_c_im (:,:) = 0.d0

  do orb_i = 1, param_orb_nb

    ex_i = ex_orb_ind (orb_i)
    orb_1st = ex_orb_1st_lab (ex_i)
    orb_2nd = ex_orb_2nd_lab (ex_i)

!   excitation between orbitals both coming from the real part or imaginary part of complex orbitals
    if (ireal_imag(orb_1st) == ireal_imag(orb_2nd)) then
      delta_c_rp (:, orb_1st) = delta_c_rp (:, orb_1st) + delta_coef_ex (ex_i) * c_rp (1:ngvec_orb, orb_2nd)
      delta_c_rm (:, orb_1st) = delta_c_rm (:, orb_1st) + delta_coef_ex (ex_i) * c_rm (1:ngvec_orb, orb_2nd)
      delta_c_ip (:, orb_1st) = delta_c_ip (:, orb_1st) + delta_coef_ex (ex_i) * c_ip (1:ngvec_orb, orb_2nd)
      delta_c_im (:, orb_1st) = delta_c_im (:, orb_1st) + delta_coef_ex (ex_i) * c_im (1:ngvec_orb, orb_2nd)

!   excitation between orbitals coming from the real part and imaginary part of complex orbitals
    else
      delta_c_rp (:, orb_1st) = delta_c_rp (:, orb_1st) + delta_coef_ex (ex_i) * c_ip (1:ngvec_orb, orb_2nd)
      delta_c_rm (:, orb_1st) = delta_c_rm (:, orb_1st) + delta_coef_ex (ex_i) * c_im (1:ngvec_orb, orb_2nd)
      delta_c_ip (:, orb_1st) = delta_c_ip (:, orb_1st) - delta_coef_ex (ex_i) * c_rp (1:ngvec_orb, orb_2nd)
      delta_c_im (:, orb_1st) = delta_c_im (:, orb_1st) - delta_coef_ex (ex_i) * c_rm (1:ngvec_orb, orb_2nd)
    endif

  enddo ! orb_i

  end subroutine delta_coef_pw_bld

! ==============================================================================
  subroutine iwcsf_bld
! ------------------------------------------------------------------------------
! Description : default csfs to optimize if not read
!
! Created     : J. Toulouse, 30 Jul 2012
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i

! header
  if (header_exe) then

   call object_create ('nparmcsf')
   call object_create ('iwcsf')

   call object_needed ('ncsf')

   return

  endif

! begin
  if (l_opt_ovlp_fn_linear) then
   nparmcsf=ncsf
   call object_alloc ('iwcsf', iwcsf, nparmcsf)
   do param_i = 1, nparmcsf
    iwcsf (param_i) = param_i
   enddo
  else
   nparmcsf=ncsf-1
   call object_alloc ('iwcsf', iwcsf, nparmcsf)
   do param_i = 1, nparmcsf
    iwcsf (param_i) = param_i + 1
   enddo
  endif

  end subroutine iwcsf_bld

end module optimization_mod
