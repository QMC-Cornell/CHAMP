module opt_lin_mod

  use all_tools_mod
  use opt_common_mod
  use deriv_mod

! Declaration of global variables and default values
  character(len=max_string_len)   :: update_nonlinear = 'semiorthogonal'
  real(dp)                        :: xi = 1.d0
  logical                         :: l_sym_ham = .false.
  logical                         :: l_opt_orb_orb_eig = .false.
  logical                         :: l_opt_orb_orb_diag = .false.
  logical                         :: l_mixed_wf_geo = .false.
  logical                         :: l_sym_ham_geo = .true.
  logical                         :: l_sym_ham_first_row_column_geo = .false.
  logical                         :: l_renormalize = .false.
  logical                         :: l_print_eigenvector_norm = .false.
  logical                         :: l_opt_lin_left_eigvec = .false.
  logical                         :: l_ham_1st_col_eq_1st_row = .false.


  real(dp), allocatable           :: ovlp_lin(:,:)
  real(dp), allocatable           :: ovlp_lin_renorm(:,:)
  real(dp), allocatable           :: ovlp_lin_eigvec(:,:)
  real(dp), allocatable           :: ovlp_lin_eigval(:)
  real(dp), allocatable           :: renorm_vector(:)
  real(dp), allocatable           :: ham_lin_energy(:,:)
  real(dp), allocatable           :: ham_lin_variance(:,:)
  real(dp), allocatable           :: ham_lin(:,:)
  real(dp), allocatable           :: ham_lin_renorm(:,:)
  real(dp), allocatable           :: ham_lin_renorm_stab(:,:)
  real(dp), allocatable           :: ham_ovlp_lin(:,:)
  real(dp), allocatable           :: ham_ovlp_lin_eigvec(:,:)
  real(dp), allocatable           :: ham_ovlp_lin_eigval(:)
  real(dp), allocatable           :: ovlp_lin_inv(:,:)

  real(dp), allocatable           :: delta_lin(:)
  real(dp)                        :: psi_lin_norm_sq

  real(dp), allocatable           :: ham_eigval_av(:)
  real(dp), allocatable           :: ham_eigval_av_err(:)
  real(dp), allocatable           :: ovlp_lin_av(:,:)
  real(dp), allocatable           :: ham_lin_energy_av(:,:)
  real(dp), allocatable           :: ham_lin_energy_av_err(:,:)

  logical                         :: l_select_eigvec_lowest = .false.
  logical                         :: l_select_eigvec_largest_1st_coef = .false.
  logical                         :: l_select_eigvec_smallest_norm = .false.

  integer                         :: target_state = 0
  integer                         :: target_state_above_groundstate = 0
  integer                         :: target_state_above_groundstate_or_target_smallest_norm = 0
  real(dp)                        :: add_diag_mult_exp = 1.d0
  real(dp)                        :: eigval_lower_bound
  real(dp)                        :: eigval_upper_bound
  logical                         :: l_eigval_lower_bound_fixed = .false.
  logical                         :: l_eigval_upper_bound_fixed = .false.
  logical                         :: l_print_eigval_errors
  real(dp)                        :: add_diag_ovlp  = 1.d-8

  contains

!===========================================================================
  subroutine opt_lin_menu
!---------------------------------------------------------------------------
! Description : menu for linear optimization
!
! Created     : J. Toulouse, 24 Apr 2006
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'opt_lin_menu'

! begin
  target_state = 0
  target_state_above_groundstate = 0
  target_state_above_groundstate_or_target_smallest_norm = 0
  l_eigval_lower_bound_fixed = .false.
  l_eigval_upper_bound_fixed = .false.
  l_opt_lin_left_eigvec = .false.
  l_ham_1st_col_eq_1st_row = .false.
  l_print_eigval_errors = .false.

! loop over menu lines
  do
  call get_next_word(word)

  select case (trim(word))
  case ('help')
   write(6,'(a)') ' HELP for linear optimization menu'
   write(6,'(a)') '  linear'
   write(6,'(a)') '   update_nonlinear = [original|semiorthogonal] : default=semiorthogonal, choice of update of nonlinear paramaters'
   write(6,'(a)') '   xi = [real] : update of nonlinear paramaters by orthogonalization to xi Psi_0 +(1-xi) Psi_lin'
   write(6,'(a)') '                - xi=1: orthogonalization to Psi_0 (default)'
   write(6,'(a)') '                - xi=0: orthogonalization to Psi_lin, ensures min |Psi_lin-Psi_0|'
   write(6,'(a)') '                - xi=0.5: orthogonalization to Psi_0 + Psi_lin, ensures |Psi_0|=|Psi_lin|'
   write(6,'(a)') '    use_orbital_eigenvalues = [logical] : approximate orbital part of Hamiltonian using orbital eigenvalues? (default=false)'
   write(6,'(a)') '    symmetrize_hamiltonian = [logical] : symmetrize Hamiltonian(default=false)'
   write(6,'(a)') '    approx_orb_orb = [logical] : approximate orbital-orbital part of Hamiltonian only (default=false)'
   write(6,'(a)') '    approx_orb_orb_diag = [logical] : diagonal only approximation for orbital-orbital block (default=false)'
   write(6,'(a)') '    mixed_wf_geo  = [logical] : use mixed wave function-geometry terms? (default=false)'
   write(6,'(a)') '    sym_ham_geo = [logical] : symmetrize geometry part of Hamiltonian(default=true)'
   write(6,'(a)') '    sym_ham_first_row_column_geo = [logical] : symmetrize first row/column for geometry part of Hamiltonian(default=false)'
   write(6,'(a)') '    renormalize = [logical] : renormalize generalized eigenvalue equation with square root of overlap matrix diagonal(default=false)'
   write(6,'(a)') '    select_eigvec_lowest = [bool] : select lowest reasonable eigenvector for ground state optimization (default=false)'
   write(6,'(a)') '    select_eigvec_largest_1st_coef = [bool] : select eigenvector with largest first coefficient for ground state optimization (default=false)'
   write(6,'(a)') '    select_eigvec_smallest_norm = [bool] : select eigenvector with smallest norm(Psi_lin-Psi_)) for nonlinear params for ground state optimization (default=false)'
   write(6,'(a)') '    eigval_lower_bound = [real] : lower bound for reasonable eigenvalue window (default: internally calculated)'
   write(6,'(a)') '    eigval_upper_bound = [real] : upper bound for reasonable eigenvalue window (default: internally calculated)'
   write(6,'(a)') '    target_state = [integer] : index of target state to optimize (default is the most reasonable ground-state)'
   write(6,'(a)') '    target_state_above_groundstate = [integer] : index of target state above the ground state to optimize (default=0 for ground-state)'
   write(6,'(a)') '    target_state_above_groundstate_or_target_smallest_norm = [integer] : target that state above ground state or target smallest norm (default=0)'
   write(6,'(a)') '    print_eigenvector_norm = [bool] : print norm of all eigenvectors? (default=false)'
   write(6,'(a)') '    left_eigvec = [bool]: use left eigenvector instead of right one? (default=false)'
   write(6,'(a)') '    ham_1st_col_eq_1st_row = [bool]: set 1st column equal to 1st row in Hamiltonian(default=false)'
   write(6,'(a)') '    print_eigval_errors = [bool]: calculate and print statistical errors on energy eigenvalues (default=false)'
   write(6,'(a)') '    add_diag_ovlp = [real] : value added to diagonal of overlap matrix in linear method (default=1.d-8)'
   write(6,'(a)') ' end'

  case ('update_nonlinear')
   call get_next_value(update_nonlinear)

  case ('xi')
   call get_next_value(xi)

  case ('use_orbital_eigenvalues')
   call get_next_value(l_opt_orb_eig)

  case ('symmetrize_hamiltonian')
   call get_next_value(l_sym_ham)

  case ('approx_orb_orb')
   call get_next_value(l_opt_orb_orb_eig)

  case ('approx_orb_orb_diag')
   call get_next_value(l_opt_orb_orb_diag)

  case ('mixed_wf_geo')
   call get_next_value(l_mixed_wf_geo)

  case ('sym_ham_geo')
   call get_next_value(l_sym_ham_geo)

  case ('sym_ham_first_row_column_geo')
   call get_next_value(l_sym_ham_first_row_column_geo)

  case ('renormalize')
   call get_next_value(l_renormalize)

  case ('select_eigvec_lowest')
   call get_next_value(l_select_eigvec_lowest)
   if(l_select_eigvec_lowest) then
     l_select_eigvec_largest_1st_coef = .false.
     l_select_eigvec_smallest_norm = .false.
   endif

  case ('select_eigvec_largest_1st_coef')
   call get_next_value(l_select_eigvec_largest_1st_coef)
   if(l_select_eigvec_largest_1st_coef) then
     l_select_eigvec_lowest = .false.
     l_select_eigvec_smallest_norm = .false.
   endif

  case ('select_eigvec_smallest_norm')
   call get_next_value(l_select_eigvec_smallest_norm)
   if(l_select_eigvec_smallest_norm) then
     l_select_eigvec_lowest = .false.
     l_select_eigvec_largest_1st_coef = .false.
   endif

  case ('eigval_lower_bound')
   call get_next_value (eigval_lower_bound)
   l_eigval_lower_bound_fixed = .true.

  case ('eigval_upper_bound')
   call get_next_value (eigval_upper_bound)
   l_eigval_upper_bound_fixed = .true.

  case ('target_state')
   call get_next_value(target_state)
   call require (lhere, 'target_state >= 0', target_state >= 0)

  case ('target_state_above_groundstate')
   call get_next_value(target_state_above_groundstate)
   call require (lhere, 'target_state_above_groundstate >= 0', target_state_above_groundstate >= 0)

  case ('target_state_above_groundstate_or_target_smallest_norm')
   call get_next_value(target_state_above_groundstate_or_target_smallest_norm)
   call require (lhere, 'target_state_above_groundstate_or_target_smallest_norm >= 0', target_state_above_groundstate_or_target_smallest_norm >= 0)

  case ('print_eigenvector_norm')
   call get_next_value (l_print_eigenvector_norm)

  case ('left_eigvec')
   call get_next_value (l_opt_lin_left_eigvec)

  case ('ham_1st_col_eq_1st_row')
   call get_next_value (l_ham_1st_col_eq_1st_row)

  case ('print_eigval_errors')
   call get_next_value (l_print_eigval_errors)

  case ('add_diag_ovlp')
   call get_next_value (add_diag_ovlp)

  case ('end')
   exit

  case default
   call die(lhere, 'unknown word >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  if(trim(update_nonlinear) == 'original') then
   write(6,'(a)') ' update of nonlinear parameters in linear optimization method will be done using the original derivatives'
   write(6,'(a)') ' Warning: this update choice is very bad for the Jastrow parameters!'
   l_warning = .true.
  else
   write(6,'(a)') ' update of nonlinear parameters in linear optimization method will be done using semiorthogonal derivatives:'
   write(6,'(a,es15.8)') ' the derivatives will be orthogonalized to [xi Psi_0 +(1-xi) Psi_lin], with xi=',xi
  endif

  if (target_state /= 0 .and. target_state_above_groundstate /= 0) then
   target_state = 0
   write(6,'(a)') ' Warning: target_state_above_groundstate /= 0, thus target_state is ignored (reset to 0)'
   l_warning = .true.
  endif
! when targeting an excited state above the ground state, select the ground state by lowest energy criterium
! (since smaller norm criterium may not be relevant)
  if (target_state_above_groundstate /= 0) then
   l_select_eigvec_lowest = .true.
   if (.not. l_eigval_lower_bound_fixed .or. .not. l_eigval_upper_bound_fixed) then
    call die (lhere, 'eigval_lower_bound and eigval_upper_bound for the ground state must be specified when using target_state_above_groundstate') 
   endif
  endif
  if (target_state_above_groundstate_or_target_smallest_norm /= 0) then
   l_select_eigvec_lowest = .true.
  endif
  if(target_state_above_groundstate/=0 .and. target_state_above_groundstate_or_target_smallest_norm/=0) then
   call die (lhere, 'only one of target_state_above_groundstate or target_state_above_groundstate_or_target_smallest_norm should be used')
  endif

  end subroutine opt_lin_menu

! ==============================================================================
  subroutine ovlp_lin_bld
! ------------------------------------------------------------------------------
! Description   : overlap matrix for linear method
! Description   : < Psi(i) | Psi(j) > / < Psi(0) | Psi(0) >
!
! Created       : J. Toulouse, 10 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j


! header
  if(header_exe) then

   call object_create('ovlp_lin')
   call object_create('ovlp_lin_eigvec')
   call object_create('ovlp_lin_eigval')

   call object_needed('param_nb')
   call object_needed('param_aug_nb')
   call object_needed('dpsi_av')
   call object_needed('dpsi_dpsi_covar')
   call object_needed ('is_param_type_orb')
   call object_needed ('is_param_type_geo')

   return

  endif

! begin

! allocations
  call object_alloc('ovlp_lin', ovlp_lin, param_aug_nb, param_aug_nb)

! first element
  ovlp_lin(1,1) = 1.d0

! first row and first column
  do i = 1, param_nb
   ovlp_lin(1,i+1) = 0.d0
   ovlp_lin(i+1,1) = 0.d0
  enddo

! derivative-derivative part
  do i = 1, param_nb
   do j = i, param_nb

!   diagonal-only approximation for orbital-orbital part
    if (l_opt_orb_orb_diag .and. is_param_type_orb(i) .and. is_param_type_orb(j) .and. i /= j) then
     ovlp_lin (i+1,j+1) = 0.d0

!   no mixed wave function-geometry terms
    elseif (l_opt_geo .and. (.not. l_mixed_wf_geo) .and.  &
    ((is_param_type_geo (i) .and. .not. is_param_type_geo (j)).or.(is_param_type_geo (j) .and. .not. is_param_type_geo (i)) ) ) then
     ovlp_lin (i+1,j+1) = 0.d0

!   normal overlap
    else
     ovlp_lin(i+1,j+1) = dpsi_dpsi_covar(i,j)

    endif
   
!   force symmetrization of overlap matrix (important for numerics?)
    if (i /= j) then
     ovlp_lin(j+1,i+1) = ovlp_lin(i+1,j+1)
    endif

   enddo
  enddo

! Warning: tmp: add to diagonal of overlap
  write(6,'(''Adding to ovlp_lin'',es12.4)') add_diag_ovlp
  do i = 1, param_aug_nb
!   ovlp_lin(i,i) = ovlp_lin(i,i)+diag_stab
    ovlp_lin(i,i) = ovlp_lin(i,i)+add_diag_ovlp
  enddo

!  do i = 1, param_nb
!    write(6,'(2a,100f12.4)') trim(here),': ovlp_lin=',(ovlp_lin(i,j),j=1,param_nb)
!  enddo

! check eigenvalues
  call object_alloc('ovlp_lin_eigvec', ovlp_lin_eigvec, param_aug_nb, param_aug_nb)
  call object_alloc('ovlp_lin_eigval', ovlp_lin_eigval, param_aug_nb)
  call eigensystem(ovlp_lin, ovlp_lin_eigvec, ovlp_lin_eigval, param_aug_nb)

  write(6,*)
  write(6,'(a)') 'Eigenvalues of overlap matrix of current wave function and its first-order derivatives:'
  do i = 1, param_aug_nb
    write(6,'(a,i4,a,es15.8)') 'overlap eigenvalue # ',i,': ',ovlp_lin_eigval(i)
  enddo

  call object_modified ('ovlp_lin')
! call object_write ('ovlp_lin')

  end subroutine ovlp_lin_bld

! ==============================================================================
  subroutine ovlp_lin_renorm_bld
! ------------------------------------------------------------------------------
! Description   : renormalized overlap matrix for linear method
! Description   : seems to be not good for exponent optimization
! Description   : when there are some very small overlap elements
! Description   : needs something better -> e.g., Lowdin's canonical orthogonalizatoin
! Helps for dots, counter-productive for exponents
!
! Created       : J. Toulouse, 13 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j

! header
  if(header_exe) then

   call object_create('ovlp_lin_renorm')
   call object_create('renorm_vector')

   call object_needed('param_aug_nb')
   call object_needed('ovlp_lin')

   return

  endif

! begin

! allocations
  call object_alloc('renorm_vector', renorm_vector, param_aug_nb)
  call object_alloc('ovlp_lin_renorm', ovlp_lin_renorm, param_aug_nb, param_aug_nb)

! computing renormalization matrix
  if(l_renormalize) then
   do i = 1, param_aug_nb
    renorm_vector(i) = dsqrt(ovlp_lin(i,i))
   enddo
  else
    renorm_vector(:) = 1.d0
  endif

! renormalized overlap matrix
  do i = 1, param_aug_nb
     do j = 1, param_aug_nb
       ovlp_lin_renorm(i,j) = ovlp_lin(i,j) /(renorm_vector(i) * renorm_vector(j))
     enddo
  enddo

  end subroutine ovlp_lin_renorm_bld

! ==============================================================================
  subroutine ham_lin_energy_bld
! ------------------------------------------------------------------------------
! Description   : Hamiltonian matrix over basis Psi(i) for linear method
!
! Created       : J. Toulouse, 10 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'ham_lin_energy_bld'
  integer i, j, pair

! header
  if(header_exe) then

   call object_create('ham_lin_energy')

   call object_needed('param_nb')
   call object_needed('param_aug_nb')
   call object_needed('param_pairs')
   call object_needed('eloc_av')
   call object_needed('deloc_av')
   call object_needed('dpsi_eloc_av')
   call object_needed('dpsi_eloc_covar')
   call object_needed('dpsi_deloc_covar')
   call object_needed('dpsi_dpsi_eloc_av')
   call object_needed('dpsi_av')
   call object_needed('is_param_type_orb')
   call object_needed('is_param_type_geo')

   return

  endif

! begin

! For approximate Hamitonian for orbitals
  if(l_opt_orb_eig .or. l_opt_orb_orb_eig) then
   call object_provide(lhere, 'ovlp_lin')
   call object_provide(lhere, 'delta_eps')
  endif

! allocations
  call object_alloc('ham_lin_energy', ham_lin_energy, param_aug_nb, param_aug_nb)

! first element
  ham_lin_energy(1,1) = eloc_av

! first row and first column
  do i = 1, param_nb

!   approximate Hamiltonian for orbitals
    if(l_opt_orb_eig .and. is_param_type_orb(i)) then
     ham_lin_energy(1+i,1) = dpsi_eloc_covar(i) 
!     ham_lin_energy(1,i+1) = ovlp_lin(1,i+1) *(eloc_av + delta_eps(i-nparmcsf-nparmj))
     ham_lin_energy(1,1+i) = ham_lin_energy(1+i,1)  ! symmetric for orbitals

!   symmetric Hamiltonian for geometry parameters?
    elseif(l_opt_geo .and. (l_sym_ham_geo .or. l_sym_ham_first_row_column_geo) .and. is_param_type_geo(i)) then
     ham_lin_energy(1+i,1) = dpsi_eloc_covar(i) + deloc_av(i)/2.d0
     ham_lin_energy(1,1+i) = ham_lin_energy(1+i,1)

!   1st column equal 1st row (mainly for DMC optimization)
    elseif(l_ham_1st_col_eq_1st_row) then
     ham_lin_energy(1+i,1) = dpsi_eloc_covar(i) + deloc_av(i)
     ham_lin_energy(1,1+i) = dpsi_eloc_covar(i) + deloc_av(i)

!   normal Hamiltoniam
    else
     ham_lin_energy(1+i,1) = dpsi_eloc_covar(i)
     ham_lin_energy(1,1+i) = dpsi_eloc_covar(i) + deloc_av(i)
    endif

  enddo ! i

! derivative-derivative part
  do j = 1, param_nb
   do i = 1, param_nb
     pair = param_pairs(i,j)

!   approximate Hamiltonian for Jastrow-orbital, CSF-orbital mixed terms(swap i and j)
    if(l_opt_orb_eig .and. .not. is_param_type_orb(i) .and. is_param_type_orb(j)) then
     ham_lin_energy(i+1,j+1) =  ham_lin_energy(j+1,i+1)

!   diagonal-only approximation for orbital-orbital block
    elseif(l_opt_orb_orb_diag .and. is_param_type_orb(i) .and. is_param_type_orb(j) .and. i /= j ) then
     ham_lin_energy(i+1,j+1) = 0.d0

!   approximate Hamiltonian for orbital-orbital terms
    elseif((l_opt_orb_eig .or. l_opt_orb_orb_eig) .and. is_param_type_orb(i) .and. is_param_type_orb(j)) then
!     ham_lin_energy(i+1,j+1) = ovlp_lin(i+1,j+1) *(eloc_av + delta_eps(j-nparmcsf-nparmj))
     ham_lin_energy(i+1,j+1) = ovlp_lin(i+1,j+1) *((eloc_av + delta_eps(j-nparmcsf-nparmj)) +(eloc_av + delta_eps(i-nparmcsf-nparmj)))/2.d0

!   no mixed wave function-geometry terms?
    elseif (l_opt_geo .and. (.not. l_mixed_wf_geo) .and.  &
    ((is_param_type_geo (i) .and. .not. is_param_type_geo (j)).or.(is_param_type_geo (j) .and. .not. is_param_type_geo (i)) ) ) then
     ham_lin_energy (i+1,j+1) =  0.d0

!   symmetric Hamiltonian for geometry-geometry block?
    elseif (l_opt_geo .and. l_sym_ham_geo .and. (is_param_type_geo (i) .and. is_param_type_geo (j))) then
     ham_lin_energy (i+1,j+1) =  dpsi_dpsi_eloc_av (pair)                                        &
                               - dpsi_av (j) * dpsi_eloc_av (i) - dpsi_av (i) * dpsi_eloc_av (j) &
                               + dpsi_av (i) * dpsi_av (j) * eloc_av                             &
                               + (dpsi_deloc_covar (i, j) + dpsi_deloc_covar (j, i))/2.d0
!JT                               + (dpsi_deloc_covar (i, j) + dpsi_deloc_covar (j, i))/1.d0

!   normal Hamiltonian (Eq. 54d of 2007 JCP)
    else
     ham_lin_energy(i+1,j+1) =  dpsi_dpsi_eloc_av(pair)                                     &
                              - dpsi_av(j) * dpsi_eloc_av(i) - dpsi_av(i) * dpsi_eloc_av(j) &
                              + dpsi_av(i) * dpsi_av(j) * eloc_av                           &
                              + dpsi_deloc_covar(i, j)
    endif

!   if(i /= j) then
!
!!   approximate Hamiltonian for Jastrow-orbital, CSF-orbital and orbital-orbital terms
!    if(l_opt_orb_eig .and. i > nparmcsf+nparmj) then
!     ham_lin_energy(j+1,i+1) = ovlp_lin(j+1,i+1) *(eloc_av + delta_eps(i-nparmcsf-nparmj))
!
!!   normal Hamiltoniam
!    else
!     ham_lin_energy(j+1,i+1) = dpsi_deloc_covar(j, i) + dpsi_dpsi_eloc_av(pair)     &
!                       - dpsi_av(i) * dpsi_eloc_av(j)                     &
!                       - dpsi_av(j) *( dpsi_eloc_av(i)) &
!                       + dpsi_av(i) * dpsi_av(j) * eloc_av
!    endif
!
!   endif !i /= j

   enddo
  enddo

! symmetrize Hamiltonian
  if(l_sym_ham) then
   ham_lin_energy =(ham_lin_energy + transpose(ham_lin_energy))/2.d0
  endif


!  write(6,*)
!  write(6,'(a)') 'Hamiltonian matrix:'
  call object_modified ('ham_lin_energy')
! call object_write ('ham_lin_energy')
!  do i = 1, param_aug_nb
!    write(6,'(100e16.8)')(ham_lin_energy(i,j),j=1,param_aug_nb)
!  enddo

!  call object_provide('ovlp_lin')
!  write(6,'(2a,100f12.4)') trim(here),': ham_lin_energy diagonal=',(ham_lin_energy(i,i)/ovlp_lin(i,i),i=1,param_aug_nb)

  end subroutine ham_lin_energy_bld

!! ==============================================================================
!  subroutine ham_lin_energy_av_bld
!! ------------------------------------------------------------------------------
!! Description   : same as ham_lin_energy but running average over the block
!! Description   : for calculating statistical error
!!
!! Created       : J. Toulouse, 22 Jan 2016
!! ------------------------------------------------------------------------------
!  use all_modules_mod
!  implicit none
!
!! local
!  character(len=max_string_len_rout), save :: lhere = 'ham_lin_energy_av_bld'
!  integer i, j, pair
!
!! header
!  if(header_exe) then
!
!   call object_create('ham_lin_energy_av')
!   call object_error_define ('ham_lin_energy_av', 'ham_lin_energy_av_err')
!
!   call object_needed('param_nb')
!   call object_needed('param_aug_nb')
!   call object_needed('param_pairs')
!   call object_needed('eloc_av')
!   call object_needed('deloc_av')
!   call object_needed('dpsi_eloc_av')
!   call object_needed('dpsi_eloc_covar')
!   call object_needed('dpsi_deloc_covar')
!   call object_needed('dpsi_dpsi_eloc_av')
!   call object_needed('dpsi_av')
!
!   return
!
!  endif
!
!! begin
!
!! allocations
!  call object_alloc('ham_lin_energy_av', ham_lin_energy_av, param_aug_nb, param_aug_nb)
!  call object_alloc('ham_lin_energy_av_err', ham_lin_energy_av_err, param_aug_nb, param_aug_nb)
!
!! first element
!  ham_lin_energy_av(1,1) = eloc_av
!
!! first row and first column
!  do i = 1, param_nb
!     ham_lin_energy_av(1+i,1) = dpsi_eloc_covar(i)
!     ham_lin_energy_av(1,1+i) = dpsi_eloc_covar(i) + deloc_av(i)
!  enddo ! i
!
!! derivative-derivative part
!  do j = 1, param_nb
!   do i = 1, param_nb
!     pair = param_pairs(i,j)
!     ham_lin_energy_av(i+1,j+1) =  dpsi_dpsi_eloc_av(pair)                                     &
!                              - dpsi_av(j) * dpsi_eloc_av(i) - dpsi_av(i) * dpsi_eloc_av(j) &
!                              + dpsi_av(i) * dpsi_av(j) * eloc_av                           &
!                              + dpsi_deloc_covar(i, j)
!   enddo
!  enddo
!
!  end subroutine ham_lin_energy_av_bld

! ==============================================================================
  subroutine ham_lin_variance_bld
! ------------------------------------------------------------------------------
! Description   : Hamiltonian matrix over basis Psi(i) for linear method
! Description   : for variance minimization
!
! Created       : J. Toulouse, 25 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j

! header
  if(header_exe) then

   call object_create('ham_lin_variance')

   call object_needed('param_nb')
   call object_needed('param_aug_nb')
   call object_needed('eloc_var')
   call object_needed('gradient_variance')
   call object_needed('hessian_variance')
   call object_needed('dpsi_dpsi_covar')

   return

  endif

! begin

! allocations
  call object_alloc('ham_lin_variance', ham_lin_variance, param_aug_nb, param_aug_nb)

! first element
  ham_lin_variance(1,1) = eloc_var

! first row and first column
  do i = 1, param_nb
    ham_lin_variance(i+1,1) = gradient_variance(i)/2.d0
    ham_lin_variance(1,i+1) = ham_lin_variance(i+1,1)
  enddo

! derivative-derivative part
  do j = 1, param_nb
   do i = 1, param_nb
     ham_lin_variance(i+1,j+1) = hessian_variance(i,j)/2.d0 + eloc_var * dpsi_dpsi_covar(i, j)
   enddo
  enddo

  end subroutine ham_lin_variance_bld

! ==============================================================================
  subroutine ham_lin_bld
! ------------------------------------------------------------------------------
! Description   : Total hamiltonian matrix over basis Psi(i) for linear method
!
! Created       : J. Toulouse, 25 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'ham_lin_bld'

! header
  if(header_exe) then

   call object_create('ham_lin')

   call object_needed('param_nb')
   call object_needed('param_aug_nb')
   call object_needed('ham_lin_energy')
   call object_needed('p_var')

   return

  endif

! begin

! allocations
  call object_alloc('ham_lin', ham_lin, param_aug_nb, param_aug_nb)

  if(p_var /= 0) then
    call object_provide(lhere, 'ham_lin_variance')
    ham_lin(:,:) =(1.d0 - p_var) * ham_lin_energy(:,:) + p_var * ham_lin_variance(:,:)
  else
    ham_lin(:,:) = ham_lin_energy(:,:)
  endif

  end subroutine ham_lin_bld

! ==============================================================================
  subroutine ham_lin_renorm_bld
! ------------------------------------------------------------------------------
! Description   : Renormalized hamiltonian matrix for linear method
!
! Created       : J. Toulouse, 13 Jul 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i, j

! header
  if(header_exe) then

   call object_create('ham_lin_renorm')

   call object_needed('param_aug_nb')
   call object_needed('ham_lin')
   call object_needed('renorm_vector')

   return

  endif

! begin

! allocations
  call object_alloc('ham_lin_renorm', ham_lin_renorm, param_aug_nb, param_aug_nb)

! renormalizing Hamiltonian matrix
  do i = 1, param_aug_nb
    do j = 1, param_aug_nb
      ham_lin_renorm(i,j) = ham_lin(i,j) /(renorm_vector(i) * renorm_vector(j))
    enddo
  enddo

  end subroutine ham_lin_renorm_bld

! ==============================================================================
  subroutine ham_lin_renorm_stab_bld
! ------------------------------------------------------------------------------
! Description   : Hamiltonian matrix over basis Psi(i) for linear method
! Description   : stabilized
!
! Created       : J. Toulouse, 24 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer i

! header
  if(header_exe) then

   call object_create('ham_lin_renorm_stab')

   call object_needed('param_nb')
   call object_needed('param_aug_nb')
   call object_needed('ham_lin_renorm')
   call object_needed('diag_stab')
   call object_needed('add_diag_mult_exp')

   return

  endif

! begin
! allocations
  call object_alloc('ham_lin_renorm_stab', ham_lin_renorm_stab, param_aug_nb, param_aug_nb)

  ham_lin_renorm_stab(:,:) = ham_lin_renorm(:,:)

  select case (trim(stabilization))

! stabilization by adding identity matrix
  case ('identity')
   do i = 1, param_nb
     if(i > nparmcsf+nparmj .and. i <= nparmcsf+nparmj+param_exp_nb) then
       ham_lin_renorm_stab(i+1,i+1) = ham_lin_renorm(i+1,i+1) + diag_stab * add_diag_mult_exp ! multiplicative factor for exponent parameters
     else
       ham_lin_renorm_stab(i+1,i+1) = ham_lin_renorm(i+1,i+1) + diag_stab
     endif
   enddo

! stabilization by adding overlap matrix
  case ('overlap')
   call object_provide('ovlp_lin')
   do i = 1, param_nb
       ham_lin_renorm_stab(i+1,i+1) = ham_lin_renorm(i+1,i+1) + diag_stab * ovlp_lin(i+1,i+1)
   enddo

! stabilization by symmetrizing the hamiltonian
  case ('symmetrize')
       ham_lin_renorm_stab = ham_lin_renorm + diag_stab/(1.d0+diag_stab) * (transpose(ham_lin_renorm) - ham_lin_renorm)/2.d0

  case default
   call die (here, 'unknown stabilization choice >'+trim(stabilization)+'<.')
  end select

  end subroutine ham_lin_renorm_stab_bld

! ==============================================================================
  subroutine delta_lin_bld
! ------------------------------------------------------------------------------
! Description   : variation of parameters for linear method
!
! Created       : J. Toulouse and Cyrus Umrigar, 10 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_lin_bld'
  integer i, j
  integer iparmcsf, iparm, jparm
  integer lwork, info
  real(dp), allocatable :: mat_a(:,:), mat_b(:,:)
  real(dp), allocatable :: eigvec(:,:) !, eigvec_left(:,:)
  real(dp), allocatable :: eigval_r(:), eigval_i(:), eigval_denom(:)
  real(dp), allocatable :: work(:)
  real(dp) eigvec_max_1st_compon, eigvec_first_coef, lowest_eigval, smallest_norm, tmp
  integer eig_ind, eigvec_max_1st_compon_ind, eigvec_lowest_eigval_ind, eigvec_smallest_norm_ind

  integer, allocatable :: eigval_srt_ind_to_eigval_ind(:), eigval_ind_to_eigval_srt_ind(:)
  integer temp, eig_excited_ind_test


! header
  if(header_exe) then

   call object_create('delta_lin')

   call object_needed('param_aug_nb')
   call object_needed('param_nb')
   call object_needed('ovlp_lin')
   call object_needed('ham_lin_renorm_stab')
   call object_needed('dpsi_av')
   call object_needed('ovlp_lin_renorm')
   call object_needed('renorm_vector')
   call object_needed('p_var')

   return

  endif

! begin

! allocation
  call object_alloc('delta_lin', delta_lin, param_nb)

! temporary arrays
  call alloc('mat_a', mat_a, param_aug_nb, param_aug_nb)
  call alloc('mat_b', mat_b, param_aug_nb, param_aug_nb)
  call alloc('eigvec', eigvec, param_aug_nb, param_aug_nb)
! call alloc('eigvec_left', eigvec_left, param_aug_nb, param_aug_nb)
  call alloc('eigval_r', eigval_r, param_aug_nb)
  call alloc('eigval_i', eigval_i, param_aug_nb)
  call alloc('eigval_denom', eigval_denom, param_aug_nb)

  mat_a(:,:) = ham_lin_renorm_stab(:,:)
  mat_b(:,:) = ovlp_lin_renorm(:,:)
! write(6,*) 'mat_a=',mat_a
! write(6,*) 'mat_b=',mat_b

! calculate optimal value of lwork
  lwork = 1
  call alloc('work', work, lwork)
  call dggev('N','V',param_aug_nb, mat_a, param_aug_nb, mat_b, param_aug_nb, eigval_r, eigval_i,  &
             eigval_denom, eigvec, param_aug_nb, eigvec, param_aug_nb, work, -1, info)
! call dggev('V','V',param_aug_nb, mat_a, param_aug_nb, mat_b, param_aug_nb, eigval_r, eigval_i,  &
!            eigval_denom, eigvec_left, param_aug_nb, eigvec, param_aug_nb, work, -1, info)
  if(info /= 0) then
    call die(lhere, 'problem in dggev(while calculating optimal value of lwork): info='+info+' /= 0')
  endif
  lwork =  nint(work(1))
  call alloc('work', work, lwork)

! solve generalized eigenvalue problem A*x = lambda*B*x
  write(6,*)
  write(6,'(a,1pd9.1)') 'Solving generalized eigenvalue equation of linear method with a_diag =', diag_stab
! note: large jobs can die in this routine because of lack of memory
  
! dggev scales eigenvectors so that the largest component of each eigenvector has abs(real)+abs(imag)=1.
  if(.not.l_opt_lin_left_eigvec) then
!   compute right eigenvectors (standard)
    call dggev('N','V',param_aug_nb, mat_a, param_aug_nb, mat_b, param_aug_nb, eigval_r, eigval_i,  &
              eigval_denom, eigvec, param_aug_nb, eigvec, param_aug_nb, work, lwork, info)
  else
!   compute left eigenvectors (this is much more noisy since the strong zero-variance principle is lost) 
    write(6,'(a)') 'Warning: use left eigenvectors instead of right eigenvectors'
    call dggev('V','N',param_aug_nb, mat_a, param_aug_nb, mat_b, param_aug_nb, eigval_r, eigval_i,  &
              eigval_denom, eigvec, param_aug_nb, eigvec, param_aug_nb, work, lwork, info)
  endif
! call dggev('V','V',param_aug_nb, mat_a, param_aug_nb, mat_b, param_aug_nb, eigval_r, eigval_i,  &
!            eigval_denom, eigvec_left, param_aug_nb, eigvec, param_aug_nb, work, lwork, info)
! write(6,*) 'after dggev'
  call release ('work', work)
  call release ('mat_a', mat_a)
  call release ('mat_b', mat_b)
  if(info /= 0) then
    call die(lhere, 'problem in dggev: info='+info+' /= 0')
  endif

! write(6,'(2a,100d10.2)') trim(lhere),': eigval_denom=', eigval_denom(:)

! calculate eigenvalue
  do i = 1, param_aug_nb
    eigval_r(i) = eigval_r(i) / eigval_denom(i)
    eigval_i(i) = eigval_i(i) / eigval_denom(i)
  enddo

  call release ('eigval_denom', eigval_denom)

! print eigenvalues
!  write(6,'(a)') 'Unsorted (complex) eigenvalues:'
!  do i = 1, param_aug_nb
!    write(6,'(a,i5,a,2(f10.6,a))') 'eigenvalue #',i,': ',eigval_r(i), ' +', eigval_i(i),' i'
!  enddo

! print eigenvectors
!  write(6,'(a)') 'Unsorted eigenvectors:'
!  do j = 1, param_aug_nb
!    write(6,'(a,i3,a,100f12.6)') 'right eigenvector # ', j,' :',(eigvec(i, j), i = 1, param_aug_nb)
!    write(6,'(a,i3,a,100f12.6)') 'left  eigenvector # ', j,' :',(eigvec_left(i, j), i = 1, param_aug_nb)
!  enddo

! Sorting out eigenvalues
! eigval_srt_ind_to_eigval_ind is the map from sorted eigenvalues to original eigenvalues
  call alloc('eigval_srt_ind_to_eigval_ind', eigval_srt_ind_to_eigval_ind, param_aug_nb)
  do i = 1, param_aug_nb
    eigval_srt_ind_to_eigval_ind(i) = i
  enddo
  write(6,'(''eigval_srt_ind_to_eigval_ind='',10000i6)') eigval_srt_ind_to_eigval_ind(1:param_aug_nb)
  do i = 1, param_aug_nb
    do j = i+1, param_aug_nb
      if(eigval_r(eigval_srt_ind_to_eigval_ind(j)) < eigval_r(eigval_srt_ind_to_eigval_ind(i))) then
        temp = eigval_srt_ind_to_eigval_ind(i)
        eigval_srt_ind_to_eigval_ind(i) = eigval_srt_ind_to_eigval_ind(j)
        eigval_srt_ind_to_eigval_ind(j) = temp
      endif
    enddo
  enddo
  write(6,'(''eigval_srt_ind_to_eigval_ind='',10000i6)') eigval_srt_ind_to_eigval_ind(1:param_aug_nb)
! eigval_ind_to_eigval_srt_ind is the map from original eigenvalues to sorted eigenvalues
  call alloc('eigval_ind_to_eigval_srt_ind', eigval_ind_to_eigval_srt_ind, param_aug_nb)
  do i = 1, param_aug_nb
    eigval_ind_to_eigval_srt_ind(eigval_srt_ind_to_eigval_ind(i)) = i
  enddo
  write(6,'(''eigval_ind_to_eigval_srt_ind='',10000i6)') eigval_ind_to_eigval_srt_ind(1:param_aug_nb)

! print eigenvalues
  write(6,'(a)') 'Sorted (complex) eigenvalues:'
  do i = 1, param_aug_nb
    write(6,'(a,i5,a,2(f12.6,a))') 'eigenvalue #',i,': ',eigval_r(eigval_srt_ind_to_eigval_ind(i)), ' +', eigval_i(eigval_srt_ind_to_eigval_ind(i)),' i'
  enddo

! print eigenvectors
! write(6,'(a)') 'Sorted right eigenvectors:'
! do j = 1, param_aug_nb
!   write(6,'(a,i3,a,100f12.6)') 'eigenvector # ', j,' :',(eigvec(i, eigval_srt_ind_to_eigval_ind(j)), i = 1, param_aug_nb)
! enddo

! Find eigenvector with smallest norm(Psi_lin-Psi_0) for nonlinear parameters
  if (l_print_eigenvector_norm) then
    write(6,'(/,a)') 'Norm of linear wave function variation for nonlinear parameters for unsorted eigenvectors:'
  endif
  smallest_norm=9.d99
  eigvec_smallest_norm_ind=1
  do i = 1, param_aug_nb
    psi_lin_var_norm = 0
    do iparm = nparmcsf+1, param_nb
      do jparm = nparmcsf+1, param_nb
        psi_lin_var_norm = psi_lin_var_norm + eigvec(1+iparm,i)*eigvec(1+jparm,i)*ovlp_lin(1+iparm,1+jparm)/(renorm_vector(1+iparm)*renorm_vector(1+jparm))
      enddo
    enddo
    psi_lin_var_norm = psi_lin_var_norm/eigvec(1,i)**2  ! Normalize so first component is 1
    if (l_print_eigenvector_norm) then
!     write(6,'(a,i5,a,f20.8)') 'eigenvector #',i, ': norm =', psi_lin_var_norm
      write(6,'(a,2i5,a,f18.8,2(f12.6,a))') 'eigenvector unsorted and sorted #',i, eigval_ind_to_eigval_srt_ind(i),': norm, eigenvalue =', psi_lin_var_norm, eigval_r(i), ' +', eigval_i(i),' i'
    endif
    if(psi_lin_var_norm < smallest_norm) then
      if(psi_lin_var_norm >= 0.d0) then
        smallest_norm = psi_lin_var_norm
        eigvec_smallest_norm_ind = i
      else
       if (l_print_eigenvector_norm) then
        write(6,'(a,es12.4,a)') 'Warning: psi_lin_var_norm=',psi_lin_var_norm ,' < 0 because of roundoff' ! This can happen only because of numerical roundoff errors
       endif
      endif
    endif
    if(psi_lin_var_norm == smallest_norm) then ! When several eigenvectors have the same psi_lin_var_norm pick the one with lowest energy
      if(eigval_r(i) < eigval_r(eigvec_smallest_norm_ind)) then
        smallest_norm = psi_lin_var_norm
        eigvec_smallest_norm_ind = i
      endif
    endif
  enddo
 
  if (l_print_eigenvector_norm) then
    write(6,'(a,i4,a,g16.8)') 'Unsorted eigenvector # ',eigvec_smallest_norm_ind,' has smallest norm of linear wave function variation for nonlinear parameters =', smallest_norm
  endif

! JT: increase upper bound by 0.01 to avoid abusive rejection
  if (.not. l_eigval_upper_bound_fixed) then
    eigval_upper_bound = (1-p_var)*energy_sav+p_var*energy_sigma_sav**2 + 0.01d0
  endif
  if (.not. l_eigval_lower_bound_fixed) then
    eigval_lower_bound = (1-p_var)*(energy_sav-energy_sigma_sav) + 0.25*p_var*energy_sigma_sav**2
  endif
  write(6,'(/,a,2f12.5)') 'Reasonable eigenvalue window is:', eigval_lower_bound, eigval_upper_bound

  write(6,'(a,t64,i4,a,2(f10.4,a),f10.6)') '(sorted) eigenvector with smallest norm of wavefn variation is #',eigval_ind_to_eigval_srt_ind(eigvec_smallest_norm_ind), ': ',eigval_r(eigvec_smallest_norm_ind), ' +', eigval_i(eigvec_smallest_norm_ind),' i, norm change=', smallest_norm
! write(6,'(a,f10.6)') 'Norm of linear wave function variation for nonlin. params for this eigenvector =', smallest_norm

! Find eigenvector with largest first coefficient
! When there are more than one with largest first coeff., pick out of these the one with the lowest eigenvalue
  eigvec_max_1st_compon = dabs(eigvec(1,1))
  eigvec_max_1st_compon_ind = 1
  do i = 1, param_aug_nb
    if(dabs(eigvec(1,i)) > eigvec_max_1st_compon) then
      eigvec_max_1st_compon = dabs(eigvec(1,i))
      eigvec_max_1st_compon_ind = i
    endif
    if(dabs(eigvec(1,i)) == eigvec_max_1st_compon) then ! When several eigenvectors have the same eigvec_max_1st_compon pick the one with lowest energy
      if(eigval_r(i) < eigval_r(eigvec_max_1st_compon_ind)) then
       eigvec_max_1st_compon = dabs(eigvec(1,i))
       eigvec_max_1st_compon_ind = i
      endif
    endif
  enddo
  psi_lin_var_norm = 0
  do iparm = nparmcsf+1, param_nb
    do jparm = nparmcsf+1, param_nb
      psi_lin_var_norm = psi_lin_var_norm + eigvec(1+iparm,eigvec_max_1st_compon_ind)*eigvec(1+jparm,eigvec_max_1st_compon_ind)*ovlp_lin(1+iparm,1+jparm)/(renorm_vector(1+iparm)*renorm_vector(1+jparm))
    enddo
  enddo
  psi_lin_var_norm = psi_lin_var_norm/eigvec(1,eigvec_max_1st_compon_ind)**2 ! Normalize so first component is 1
  write(6,'(a,t64,i4,a,2(f10.4,a),f10.6)') '(sorted) eigenvector with largest first coefficient is #',eigval_ind_to_eigval_srt_ind(eigvec_max_1st_compon_ind), ': ',eigval_r(eigvec_max_1st_compon_ind), ' +', eigval_i(eigvec_max_1st_compon_ind),' i, norm change=', psi_lin_var_norm

! Find eigenvector with lowest eigenvalue that is in a reasonable window
! If eigenvec with lowest eigenvalue is degenerate, choose the one with largest first coef
  lowest_eigval = 9.d99
  eigvec_lowest_eigval_ind = 0
  do i = 1, param_aug_nb
    if(eigval_r(i) < lowest_eigval .and. eigval_r(i) < eigval_upper_bound .and. eigval_r(i) > eigval_lower_bound) then
! Old criteria
!   if((p_var < 1.d0 .and. eigval_r(i) < lowest_eigval .and. dabs(eigval_r(i)-(1-p_var)*etrial) < 10.d0) .or. &
!      (p_var == 1.d0 .and. eigval_r(i) < lowest_eigval .and. eigval_r(i) > 0.d0)) then
      lowest_eigval = eigval_r(i)
      eigvec_lowest_eigval_ind = i
    endif
    if(eigval_r(i) == lowest_eigval) then
      if(dabs(eigvec(1,i)) > dabs(eigvec(1, eigvec_lowest_eigval_ind))) then
        lowest_eigval = eigval_r(i)
        eigvec_lowest_eigval_ind = i
      endif
    endif
  enddo

! Find the norm of the change for the eigenvector with the lowest reasonable eigenvalue
  if(eigvec_lowest_eigval_ind /= 0) then
    psi_lin_var_norm = 0
    do iparm = nparmcsf+1, param_nb
      do jparm = nparmcsf+1, param_nb
        psi_lin_var_norm = psi_lin_var_norm + eigvec(1+iparm,eigvec_lowest_eigval_ind)*eigvec(1+jparm,eigvec_lowest_eigval_ind)*ovlp_lin(1+iparm,1+jparm)/(renorm_vector(1+iparm)*renorm_vector(1+jparm))
      enddo
    enddo
    psi_lin_var_norm = psi_lin_var_norm/eigvec(1,eigvec_lowest_eigval_ind)**2 ! Normalize so first component is 1
    if(psi_lin_var_norm < 0.d0) psi_lin_var_norm= 1.d99
    write(6,'(a,t64,i4,a,2(f10.4,a),f10.6)') '(sorted) eigenvector with lowest reasonable eigenvalue is #',eigval_ind_to_eigval_srt_ind(eigvec_lowest_eigval_ind), ': ',eigval_r(eigvec_lowest_eigval_ind), ' +', eigval_i(eigvec_lowest_eigval_ind),' i, norm change=', psi_lin_var_norm
!   write(6,'(a,f10.6)') 'Norm of linear wave function variation for nonlin. params for this eigenvector =', psi_lin_var_norm
  else
    write(6,'(a)') 'Warning: all the eigenvalues are outside the reasonable eigenvalue windows so select eigenvector with smallest norm!'
  endif
  call systemflush(6)

! if target_state = 0, select eigenvector from one of the 3 criteria:
! 1) lowest reasonable eigenvalue,
!    if no reasonable lowest eigenvalue found, then take the one with smallest norm(Psi_lin-Psi_0)
! 2) largest 1st component relative and from these the one with the lowest eigenvalue
! 3) smallest norm(Psi_lin-Psi_0) for nonlinear parameters
  if(target_state == 0) then

    if(l_select_eigvec_lowest .and. eigvec_lowest_eigval_ind /= 0) then
      eig_ind = eigvec_lowest_eigval_ind
    elseif(l_select_eigvec_largest_1st_coef) then
      eig_ind = eigvec_max_1st_compon_ind
    elseif(l_select_eigvec_smallest_norm .or. (l_select_eigvec_lowest .and. eigvec_lowest_eigval_ind == 0)) then
      eig_ind = eigvec_smallest_norm_ind

!   default selection criterion:
    else
      if(eigvec_lowest_eigval_ind /= 0 .and. psi_lin_var_norm < 100*smallest_norm) then
        eig_ind = eigvec_lowest_eigval_ind
      else
        eig_ind = eigvec_smallest_norm_ind
      endif

    endif

! if target_state >= 1, simply select the corresponding wanted target state
  else

    if(target_state > param_aug_nb) then
      call die(lhere, 'target_state='+target_state+' > param_aug_nb='+param_aug_nb)
    endif
    eig_ind = eigval_srt_ind_to_eigval_ind(target_state)

  endif

!! If we are targeting the ground state ideally all 3 criteria for selecting the eigenvector should coincide
! Ideally all 3 criteria for selecting the eigenvector should coincide
! if((target_state == 0 .and. target_state_above_groundstate == 0 .and. target_state_above_groundstate_or_target_smallest_norm == 0) .and. (eigvec_lowest_eigval_ind /= eigvec_max_1st_compon_ind .or. eigvec_max_1st_compon_ind /= eigvec_smallest_norm_ind)) then
  if(eigvec_lowest_eigval_ind == 0) then
    write(6,'(a,2(i4,a,2(f8.3,a)))') 'Warning: lowest_eigval, largest_1st_coef, smallest_norm eigenvecs are:, no reasonable eigenvalue', &
    eigval_ind_to_eigval_srt_ind(eigvec_max_1st_compon_ind), ': ',eigval_r(eigvec_max_1st_compon_ind), ' +', eigval_i(eigvec_max_1st_compon_ind),' i', &
    eigval_ind_to_eigval_srt_ind(eigvec_smallest_norm_ind), ': ',eigval_r(eigvec_smallest_norm_ind), ' +', eigval_i(eigvec_smallest_norm_ind),' i'
    l_warning = .true.
  else
    if((target_state == 0) .and. (eigval_srt_ind_to_eigval_ind (eigval_ind_to_eigval_srt_ind (eigvec_lowest_eigval_ind) + target_state_above_groundstate + target_state_above_groundstate_or_target_smallest_norm) /= eigvec_max_1st_compon_ind .or. eigvec_max_1st_compon_ind /= eigvec_smallest_norm_ind)) then
      write(6,'(a,3(i4,a,2(f8.3,a)))') 'Warning: lowest_eigval, largest_1st_coef, smallest_norm eigenvecs are:', &
      eigval_ind_to_eigval_srt_ind(eigvec_lowest_eigval_ind), ': ',eigval_r(eigvec_lowest_eigval_ind), ' +', eigval_i(eigvec_lowest_eigval_ind),' i', &
      eigval_ind_to_eigval_srt_ind(eigvec_max_1st_compon_ind), ': ',eigval_r(eigvec_max_1st_compon_ind), ' +', eigval_i(eigvec_max_1st_compon_ind),' i', &
      eigval_ind_to_eigval_srt_ind(eigvec_smallest_norm_ind), ': ',eigval_r(eigvec_smallest_norm_ind), ' +', eigval_i(eigvec_smallest_norm_ind),' i'
      l_warning = .true.
    endif
  endif
  call systemflush(6)

! possibility of selecting an eigenvector above the selected ground state for excited states
! Using target_state_above_groundstate is less robust than using target_state_above_groundstate_or_target_smallest_norm
  if (target_state_above_groundstate /= 0) then
    if (eigvec_lowest_eigval_ind ==0) then
      call die (lhere, 'target_state_above_groundstate /= 0 and no ground state found in the energy window')
    endif
    write(6,'(a,t64,i4,a,2(f10.4,a))') 'selected (sorted) ground state eigenvector is  #',eigval_ind_to_eigval_srt_ind(eig_ind), ': ',eigval_r(eig_ind), ' +', eigval_i(eig_ind),' i'
    eig_ind = eigval_srt_ind_to_eigval_ind (eigval_ind_to_eigval_srt_ind (eig_ind) + target_state_above_groundstate)
    write(6,'(a,i4,a)') 'The excited state # ',target_state_above_groundstate,' above this ground state will be selected'
  endif

! At this point eig_ind could either point to the guessed ground state (if eigvec_lowest_eigval_ind /= 0) or the guessed excited state (the state with the change that has the smallest norm).
  if (target_state_above_groundstate_or_target_smallest_norm /= 0) then
    if(eigvec_lowest_eigval_ind /= 0) then
      write(6,'(''psi_lin_var_norm, smallest_norm'',9es12.4)') psi_lin_var_norm, smallest_norm
      eig_excited_ind_test = eigval_srt_ind_to_eigval_ind (eigval_ind_to_eigval_srt_ind (eig_ind) + target_state_above_groundstate_or_target_smallest_norm)
      psi_lin_var_norm = 0.d0
      do iparm = nparmcsf+1, param_nb
        do jparm = nparmcsf+1, param_nb
! Warning: Shouldn't the next line have /(renorm_vector(1+iparm)*renorm_vector(1+jparm))
          psi_lin_var_norm = psi_lin_var_norm + eigvec(1+iparm,eig_excited_ind_test)*eigvec(1+jparm,eig_excited_ind_test)*ovlp_lin(1+iparm,1+jparm)
        enddo
      enddo
      psi_lin_var_norm = psi_lin_var_norm/eigvec(1,eig_excited_ind_test)**2 ! Normalize so first component is 1
      write(6,'(a,f10.6)') 'Norm of linear wave function variation for nonlin. params for this excited trial eigenvector =', psi_lin_var_norm
    endif
    if(eigvec_lowest_eigval_ind /=0 .and. psi_lin_var_norm < 100*smallest_norm) then
      write(6,'(a,t64,i4,a,2(f10.4,a))') 'selected ground state eigenvector is  #',eigval_ind_to_eigval_srt_ind(eig_ind), ': ',eigval_r(eig_ind), ' +', eigval_i(eig_ind),' i'
!     eig_ind = eigval_srt_ind_to_eigval_ind (eigval_ind_to_eigval_srt_ind (eig_ind) + target_state_above_groundstate_or_target_smallest_norm)
      eig_ind = eig_excited_ind_test
      write(6,'(a,i4,a)') 'The excited state # ',target_state_above_groundstate_or_target_smallest_norm,' above this ground state will be selected'
    else
      eig_ind = eigvec_smallest_norm_ind
      write(6,'(a)') 'The state with the change that has the smallest norm will be selected'
! This seems like a good check, but it can hurt at times because the ground state can get shifted by add_diag to a higher energy than the excited state.
      if(eigval_ind_to_eigval_srt_ind (eig_ind) < 1+target_state_above_groundstate_or_target_smallest_norm) then
        eig_ind = eigval_srt_ind_to_eigval_ind (1+target_state_above_groundstate_or_target_smallest_norm)
        write(6,'(a,i6)') 'Warning: state with the change that has the smallest norm is not sufficiently excited, so choose instead (sorted) eigenvector #',eigval_ind_to_eigval_srt_ind(eig_ind)
      endif
    endif
  endif

! undo renormalization
! warning: only for selected eigenvector
  eigvec(:, eig_ind) = eigvec(:, eig_ind) / renorm_vector(:)

! Warning: tmp
  write(6,'(/,''eigvec(1, eig_ind)'',9es12.4)') eigvec(1,eig_ind)

! normalize eigenvector so that first component is 1
  tmp=eigvec(1,eig_ind)
! eigvec(:,eig_ind) = eigvec(:,eig_ind)/eigvec(1,eig_ind)
  eigvec(1:param_aug_nb, eig_ind) = eigvec(1:param_aug_nb, eig_ind)/tmp

! norm of linear wave function variation for nonlinear parameter
  psi_lin_var_norm = 0.d0
  do iparm = nparmcsf+1, param_nb
   do jparm = nparmcsf+1, param_nb
     psi_lin_var_norm = psi_lin_var_norm + eigvec(1+iparm,eig_ind)*eigvec(1+jparm,eig_ind)*ovlp_lin(1+iparm,1+jparm)
   enddo
  enddo
  psi_lin_var_norm = psi_lin_var_norm/eigvec(1,eig_ind)**2 ! Normalize so first component is 1
  write(6,'(a,t64,i4,a,2(f10.4,a),f10.6)') 'selected (sorted) eigenvector is #',eigval_ind_to_eigval_srt_ind(eig_ind), ': ',eigval_r(eig_ind), ' +', eigval_i(eig_ind),' i, norm change=', psi_lin_var_norm
! write(6,'(a,f10.6)') 'Norm of linear wavefunction variation for nonlin. params for this eigenvector =', psi_lin_var_norm
  if(psi_lin_var_norm > 100*smallest_norm) write(6,'(''Warning: psi_lin_var_norm > 100*smallest_norm'',2es12.4)')  psi_lin_var_norm,smallest_norm

  call release ('eigval_r', eigval_r)
  call release ('eigval_i', eigval_i)

! Warning: tmp  This is truly bizarre that when the original value of eigvec(1,eig_ind) is not 1, eigvec(1,eig_ind) and eigvec_first_coef are different!
! write(6,'(/,''eigvec(1, eig_ind)'',9es12.4)') eigvec(1,eig_ind)

! calculate the actual parameter variations
  eigvec_first_coef = eigvec(1,eig_ind)

! write(6,'(/,''Warning eigvec_first_coef='',9es12.4)') eigvec_first_coef
  if(eigvec_first_coef.ne.eigvec(1,eig_ind)) then
    write(6,'(/,''This happens because of compiler bug in gfortran 4.6.3-1ubuntu5.  It is fixed in 4.8.4-2ubuntu1~14.04.1'')')
    write(6,'(''eigvec(1, eig_ind), eigvec_first_coef'',9es12.4)') eigvec(1,eig_ind), eigvec_first_coef
    call die(lhere,'compiler bug in gfortran version 4.6.3-1ubuntu5, fixed in 4.8.4-2ubuntu1~14.04.1')
  endif

  select case (trim(update_nonlinear))

! original: come back to original derivatives basis for all parameters
   case ('original')
    do i = 1, param_nb
      eigvec_first_coef = eigvec_first_coef - eigvec(1+i,eig_ind) * dpsi_av(i)
    enddo

! semiorthogonal: use semiorthognal derivatives for nonlinear parameters
   case ('semiorthogonal')

!   come back to original derivatives for the CSFs only
    do iparmcsf = 1, nparmcsf
      eigvec_first_coef = eigvec_first_coef - eigvec(1+iparmcsf,eig_ind) * dpsi_av(iparmcsf)
    enddo

!   nonlinear parameter
    eigvec_first_coef = eigvec_first_coef +(1.d0-xi)*psi_lin_var_norm/((1.d0-xi) + xi*(1.d0+psi_lin_var_norm))

  case default
    call die(lhere, 'unknown update choice >'+trim(update_nonlinear)+'<')
  end select

! Warning: tmp
  write(6,'(/,''add_diag, delta_lin1='',es9.1,1000es12.4)') diag_stab, eigvec(1:param_nb+1, eig_ind)

! final parameter variations
  do iparm = 1, param_nb
    delta_lin(iparm) = eigvec(1+iparm, eig_ind) / eigvec_first_coef
  enddo

! Warning: tmp
  write(6,'(/,''add_diag, delta_lin2='',es9.1,12x,1000es12.4)') diag_stab, delta_lin(1:param_nb)

  call release ('eigvec', eigvec)

  end subroutine delta_lin_bld

! ==============================================================================
  subroutine psi_lin_norm_sq_bld
! ------------------------------------------------------------------------------
! Description   : square of norm of linear wave function(in the original basis)
!
! Created       : J. Toulouse, 18 May 2008
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer param_i, param_j, pair

! header
  if(header_exe) then

   call object_create('psi_lin_norm_sq')

   call object_needed('param_nb')
   call object_needed('delta_lin')
   call object_needed('dpsi_av')
   call object_needed('dpsi_dpsi_av')
   call object_needed('param_pairs')

   return

  endif

! begin

! allocation
  call object_associate('psi_lin_norm_sq', psi_lin_norm_sq)

  psi_lin_norm_sq = 1.d0

  do param_i = 1, param_nb
    psi_lin_norm_sq = psi_lin_norm_sq + 2.d0 * delta_lin(param_i) * dpsi_av(param_i)
    do param_j = 1, param_nb
      pair = param_pairs(param_i, param_j)
      psi_lin_norm_sq = psi_lin_norm_sq + delta_lin(param_i) * delta_lin(param_j) * dpsi_dpsi_av(pair)
    enddo ! param_j
  enddo ! param_i

  end subroutine psi_lin_norm_sq_bld

! ==============================================================================
  subroutine ham_eigval_av_bld
! ------------------------------------------------------------------------------
! Description   : running average of eigenvalues of Hamiltonian
! Description   : for calculating statistical errors
!
! Created       : J. Toulouse, 18 Feb 2016
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'ham_eigval_av_bld'
  integer i, j, pair
  integer lwork, info
  real(dp), allocatable :: mat_a(:,:), mat_b(:,:)
  real(dp), allocatable :: eigvec(:,:)
  real(dp), allocatable :: eigval_r(:), eigval_i(:), eigval_denom(:)
  real(dp), allocatable :: work(:)
  integer, allocatable :: eigval_srt_ind_to_eigval_ind(:), eigval_ind_to_eigval_srt_ind(:)
  integer temp


! header
  if(header_exe) then

   call object_create ('ham_eigval_av')
   call object_create ('ham_lin_energy_av')
   call object_create ('ovlp_lin_av')
   call object_error_define ('ham_lin_energy_av', 'ham_lin_energy_av_err')
   call object_error_define ('ham_eigval_av', 'ham_eigval_av_err')

   call object_needed('param_nb')
   call object_needed('param_aug_nb')
   call object_needed('param_pairs')
   call object_needed('eloc_av')
   call object_needed('deloc_av')
   call object_needed('dpsi_eloc_av')
   call object_needed('dpsi_eloc_covar')
   call object_needed('dpsi_deloc_covar')
   call object_needed('dpsi_dpsi_eloc_av')
   call object_needed('dpsi_av')
   call object_needed('dpsi_dpsi_covar')

   return

  endif

! begin

! allocation
  call object_alloc('ham_eigval_av', ham_eigval_av, param_aug_nb)
  call object_alloc('ham_eigval_av_err', ham_eigval_av_err, param_aug_nb)
  call object_alloc('ham_lin_energy_av', ham_lin_energy_av, param_aug_nb, param_aug_nb)
  call object_alloc('ham_lin_energy_av_err', ham_lin_energy_av_err, param_aug_nb, param_aug_nb)
  call object_alloc('ovlp_lin_av', ovlp_lin_av, param_aug_nb, param_aug_nb)

! Calculate Hamiltonian matrix:
! first element
  ham_lin_energy_av(1,1) = eloc_av

! first row and first column
  do i = 1, param_nb
     ham_lin_energy_av(1+i,1) = dpsi_eloc_covar(i)
     ham_lin_energy_av(1,1+i) = dpsi_eloc_covar(i) + deloc_av(i)
  enddo ! i

! derivative-derivative part
  do j = 1, param_nb
   do i = 1, param_nb
     pair = param_pairs(i,j)
     ham_lin_energy_av(i+1,j+1) =  dpsi_dpsi_eloc_av(pair)                                     &
                              - dpsi_av(j) * dpsi_eloc_av(i) - dpsi_av(i) * dpsi_eloc_av(j) &
                              + dpsi_av(i) * dpsi_av(j) * eloc_av                           &
                              + dpsi_deloc_covar(i, j)
   enddo
  enddo

! Calculate overlap matrix:
! first element
  ovlp_lin_av(1,1) = 1.d0

! first row and first column
  do i = 1, param_nb
   ovlp_lin_av(1,i+1) = 0.d0
   ovlp_lin_av(i+1,1) = 0.d0
  enddo

! derivative-derivative part
  do i = 1, param_nb
   do j = i, param_nb
     ovlp_lin_av(i+1,j+1) = dpsi_dpsi_covar(i,j)
   enddo
!   force symmetrization of overlap matrix (important for numerics?)
    if (i /= j) then
     ovlp_lin_av(j+1,i+1) = ovlp_lin_av(i+1,j+1)
    endif
  enddo

! Solve generalized eigenvalue equation

! temprorary arrays
  call alloc('mat_a', mat_a, param_aug_nb, param_aug_nb)
  call alloc('mat_b', mat_b, param_aug_nb, param_aug_nb)
  call alloc('eigvec', eigvec, param_aug_nb, param_aug_nb)
  call alloc('eigval_r', eigval_r, param_aug_nb)
  call alloc('eigval_i', eigval_i, param_aug_nb)
  call alloc('eigval_denom', eigval_denom, param_aug_nb)

  mat_a(:,:) = ham_lin_energy_av(:,:)
  mat_b(:,:) = ovlp_lin_av(:,:)

! calculate optimal value of lwork
  lwork = 1
  call alloc('work', work, lwork)
  call dggev('N','V',param_aug_nb, mat_a, param_aug_nb, mat_b, param_aug_nb, eigval_r, eigval_i,  &
             eigval_denom, eigvec, param_aug_nb, eigvec, param_aug_nb, work, -1, info)
  if(info /= 0) then
   call die(lhere, 'problem in dggev(while calculating optimal value of lwork): info='+info+' /= 0')
  endif
  lwork =  nint(work(1))
  call alloc('work', work, lwork)

  call dggev('N','V',param_aug_nb, mat_a, param_aug_nb, mat_b, param_aug_nb, eigval_r, eigval_i,  &
              eigval_denom, eigvec, param_aug_nb, eigvec, param_aug_nb, work, lwork, info)
  call release ('work', work)
  call release ('mat_a', mat_a)
  call release ('mat_b', mat_b)
  if(info /= 0) then
   call die(lhere, 'problem in dggev: info='+info+' /= 0')
  endif

! calculate eigenvalue
  do i = 1, param_aug_nb
    eigval_r(i) = eigval_r(i) / eigval_denom(i)
    eigval_i(i) = eigval_i(i) / eigval_denom(i)
  enddo

  call release ('eigval_denom', eigval_denom)

! print eigenvalues
!  write(6,'(a)') 'Unsorted (complex) eigenvalues:'
!  do i = 1, param_aug_nb
!    write(6,'(a,i5,a,2(f10.6,a))') 'eigenvalue #',i,': ',eigval_r(i), ' +', eigval_i(i),' i'
!  enddo

! print eigenvectors
!  write(6,'(a)') 'Unsorted eigenvectors:'
!  do j = 1, param_aug_nb
!    write(6,'(a,i3,a,100f12.6)') 'right eigenvector # ', j,' :',(eigvec(i, j), i = 1, param_aug_nb)
!  enddo

! Sorting out eigenvalues
! eigval_srt_ind_to_eigval_ind is the map from sorted eigenvalues to original eigenvalues
  call alloc('eigval_srt_ind_to_eigval_ind', eigval_srt_ind_to_eigval_ind, param_aug_nb)
  do i = 1, param_aug_nb
    eigval_srt_ind_to_eigval_ind(i) = i
  enddo
  do i = 1, param_aug_nb
    do j = i+1, param_aug_nb
      if(eigval_r(eigval_srt_ind_to_eigval_ind(j)) < eigval_r(eigval_srt_ind_to_eigval_ind(i))) then
        temp = eigval_srt_ind_to_eigval_ind(i)
        eigval_srt_ind_to_eigval_ind(i) = eigval_srt_ind_to_eigval_ind(j)
        eigval_srt_ind_to_eigval_ind(j) = temp
      endif
    enddo
  enddo
! eigval_ind_to_eigval_srt_ind is the map from original eigenvalues to sorted eigenvalues
  call alloc('eigval_ind_to_eigval_srt_ind', eigval_ind_to_eigval_srt_ind, param_aug_nb)
  do i = 1, param_aug_nb
   eigval_ind_to_eigval_srt_ind(eigval_srt_ind_to_eigval_ind(i)) = i
  enddo

! print eigenvalues
  write(6,'(a)') 'Sorted (complex) eigenvalues:'
  do i = 1, param_aug_nb
    write(6,'(a,i5,a,2(f12.6,a))') 'eigenvalue #',i,': ',eigval_r(eigval_srt_ind_to_eigval_ind(i)), ' +', eigval_i(eigval_srt_ind_to_eigval_ind(i)),' i'
  enddo

! save sorted eigenvalues in ham_eigval_av
  do i = 1, param_aug_nb
   ham_eigval_av (i) = eigval_r(eigval_srt_ind_to_eigval_ind(i))
  enddo

  call release ('eigval_r', eigval_r)
  call release ('eigval_i', eigval_i)
  call release ('eigvec', eigvec)

  end subroutine ham_eigval_av_bld

end module opt_lin_mod
