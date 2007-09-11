module opt_lin_mod

  use all_tools_mod
  use deriv_mod

! Declaration of global variables and default values
  character(len=max_string_len)   :: update_nonlinear = 'semiorthogonal'
  real(dp)                        :: xi = 1.d0
  logical                         :: l_sym_ham = .false.
  logical                         :: l_opt_orb_orb_eig = .false.
  logical                         :: l_opt_orb_orb_diag = .false.
  logical                         :: l_renormalize = .false.

  integer                         :: param_aug_nb

  real(dp), allocatable           :: ovlp_lin (:,:)
  real(dp), allocatable           :: ovlp_lin_renorm (:,:)
  real(dp), allocatable           :: ovlp_lin_eigvec (:,:)
  real(dp), allocatable           :: ovlp_lin_eigval (:)
  real(dp), allocatable           :: renorm_vector (:)
  real(dp), allocatable           :: ham_lin_energy (:,:)
  real(dp), allocatable           :: ham_lin_variance (:,:)
  real(dp), allocatable           :: ham_lin (:,:)
  real(dp), allocatable           :: ham_lin_renorm (:,:)
  real(dp), allocatable           :: ham_lin_renorm_stab (:,:)
  real(dp), allocatable           :: ham_ovlp_lin (:,:)
  real(dp), allocatable           :: ham_ovlp_lin_eigvec (:,:)
  real(dp), allocatable           :: ham_ovlp_lin_eigval (:)
  real(dp), allocatable           :: ovlp_lin_inv (:,:)

  real(dp), allocatable           :: delta_lin (:)
  real(dp)                        :: psi_lin_var_norm = 0.d0
  real(dp)                        :: psi_lin_var_norm_max = 10.d0

  logical                         :: l_select_eigvec_lowest = .true.
  logical                         :: l_select_eigvec_largest_1st_coef = .false.
  integer                         :: target_state = 0

  contains

!===========================================================================
  subroutine opt_lin_menu
!---------------------------------------------------------------------------
! Description : menu for linear optimization
!
! Created     : J. Toulouse, 24 Apr 2006
!---------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'opt_lin_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case(trim(word))
  case ('help')
   write(6,'(a)') 'HELP for linear optimization menu'
   write(6,'(a)') ' linear'
   write(6,'(a)') '  update_nonlinear = [original|semiorthogonal] : default=semiorthogonal, choice of update of nonlinear paramaters'
   write(6,'(a)') '  xi = [real] : update of nonlinear paramaters by orthogonalization to xi Psi_0 + (1-xi) Psi_lin'
   write(6,'(a)') '               - xi=1: orthogonalization to Psi_0 (default)'
   write(6,'(a)') '               - xi=0: orthogonalization to Psi_lin, ensures min |Psi_lin-Psi_0|'
   write(6,'(a)') '               - xi=0.5: orthogonalization to Psi_0 + Psi_lin, ensures |Psi_0|=|Psi_lin|'
   write(6,'(a)') '   use_orbital_eigenvalues = [logical] : approximate orbital part of Hamiltonian using orbital eigenvalues? (default=false)'
   write(6,'(a)') '   symmetrize_hamiltonian = [logical] : symmetrize Hamiltonian (default=false)'
   write(6,'(a)') '   approx_orb_orb = [logical] : approximate orbital-orbital part of Hamiltonian only (default=false)'
   write(6,'(a)') '   approx_orb_orb_diag = [logical] : diagonal only approximation for orbital-orbital block (default=false)'
   write(6,'(a)') '   renormalize = [logical] : renormalize generalized eigenvalue equation with square root of overlap matrix diagonal (default=false)'
   write(6,'(a)') '   select_eigvec_lowest = [bool] : select lowest reasonable eigenvector for ground state optimization (default=true)'
   write(6,'(a)') '   select_eigvec_largest_1st_coef = [bool] : select eigenvector with largest first coefficient for ground state optimization (default=false)'
   write(6,'(a)') '   target_state = [integer] : index of target state to optimize (default is ground-state)'
   write(6,'(a)') ' end'

  case ('update_nonlinear')
   call get_next_value (update_nonlinear)

  case ('xi')
   call get_next_value (xi)

  case ('use_orbital_eigenvalues')
   call get_next_value (l_opt_orb_eig)

  case ('symmetrize_hamiltonian')
   call get_next_value (l_sym_ham)

  case ('approx_orb_orb')
   call get_next_value (l_opt_orb_orb_eig)

  case ('approx_orb_orb_diag')
   call get_next_value (l_opt_orb_orb_diag)

  case ('renormalize')
   call get_next_value (l_renormalize)

  case ('select_eigvec_lowest')
   call get_next_value (l_select_eigvec_lowest)
   l_select_eigvec_largest_1st_coef = .false.

  case ('select_eigvec_largest_1st_coef')
   call get_next_value (l_select_eigvec_largest_1st_coef)
   l_select_eigvec_lowest = .false.

  case ('target_state')
   call get_next_value (target_state)
   call require ('target_state >= 0', target_state >= 0)

  case ('end')
   exit

  case default
   call die (lhere, 'unknown word >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  if (trim(update_nonlinear) == 'original') then
   write(6,'(2a)') trim(lhere),': update of nonlinear parameters in linear optimization method will be done using the original derivatives'
   write(6,'(2a)') trim(lhere),': WARNING: this is very bad for the Jastrow parameters!'
  else
   write(6,'(2a)') trim(lhere),': update of nonlinear parameters in linear optimization method will be done using semiorthogonal derivatives'
   write(6,'(2a,f)') trim(lhere),': the derivatives will be orthogonalized to [xi Psi_0 + (1-xi) Psi_lin], with xi=',xi
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
  implicit none
  include 'commons.h'

! local
  integer i, j


! header
  if (header_exe) then

   call object_create ('param_aug_nb')
   call object_create ('ovlp_lin')
   call object_create ('ovlp_lin_eigvec')
   call object_create ('ovlp_lin_eigval')

   call object_needed ('param_nb')
   call object_needed ('dpsi_av')
   call object_needed ('dpsi_dpsi_c_av')

   return

  endif

! begin

! allocations
  param_aug_nb = param_nb + 1
  call object_alloc ('ovlp_lin', ovlp_lin, param_aug_nb, param_aug_nb)

! first element
  ovlp_lin (1,1) = 1.d0

! first row and first column
  do i = 1, param_nb
   ovlp_lin (1,i+1) = 0.d0
   ovlp_lin (i+1,1) = 0.d0
  enddo

! derivative-derivative part
  do i = 1, param_nb
   do j = 1, param_nb

!   diagonal-only approximation for orbital-orbital part
    if (l_opt_orb_orb_diag .and. i > nparmcsf+nparmj .and. j > nparmcsf+nparmj .and. i /= j) then

     ovlp_lin (i+1,j+1) = 0.d0

!   normal overlap
    else

     ovlp_lin (i+1,j+1) = dpsi_dpsi_c_av (i,j)

    endif

!    if (i /= j) then
!     ovlp_lin (j+1,i+1) = ovlp_lin (i+1,j+1)
!    endif

   enddo
  enddo

!  do i = 1, param_nb
!    write(6,'(2a,100f12.4)') trim(here),': ovlp_lin=',(ovlp_lin(i,j),j=1,param_nb)
!  enddo

! check eigenvalues
  call object_alloc ('ovlp_lin_eigvec', ovlp_lin_eigvec, param_aug_nb, param_aug_nb)
  call object_alloc ('ovlp_lin_eigval', ovlp_lin_eigval, param_aug_nb)
  call eigensystem (ovlp_lin, ovlp_lin_eigvec, ovlp_lin_eigval, param_aug_nb)

  write(6,*)
  write(6,'(a)') 'Eigenvalues of overlap matrix of current wave function and its first-order derivatives:'
  do i = 1, param_aug_nb
    write(6,'(a,i3,a,e)') 'eigenvalue # ',i,': ',ovlp_lin_eigval(i)
  enddo

  end subroutine ovlp_lin_bld

! ==============================================================================
  subroutine ovlp_lin_renorm_bld
! ------------------------------------------------------------------------------
! Description   : renormalized overlap matrix for linear method
! Description   : seems to be not good for exponent optimization
! Description   : when there are some very small overlap elements
! Description   : needs something better -> e.g., Lowdin's canonical orthogonalizatoin
!
! Created       : J. Toulouse, 13 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer i, j

! header
  if (header_exe) then

   call object_create ('ovlp_lin_renorm')
   call object_create ('renorm_vector')

   call object_needed ('param_aug_nb')
   call object_needed ('ovlp_lin')

   return

  endif

! begin

! allocations
  call object_alloc ('renorm_vector', renorm_vector, param_aug_nb)
  call object_alloc ('ovlp_lin_renorm', ovlp_lin_renorm, param_aug_nb, param_aug_nb)

! computing renormalization matrix
  if (l_renormalize) then
   do i = 1, param_aug_nb
    renorm_vector (i) = dsqrt(ovlp_lin (i,i))
   enddo
  else
    renorm_vector (:) = 1.d0
  endif

! renormalized overlap matrix
  do i = 1, param_aug_nb
     do j = 1, param_aug_nb
       ovlp_lin_renorm (i,j) = ovlp_lin (i,j) / (renorm_vector (i) * renorm_vector (j))
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
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'ham_lin_energy_bld'
  integer i, j, pair

! header
  if (header_exe) then

   call object_create ('ham_lin_energy')

   call object_needed ('param_nb')
   call object_needed ('param_aug_nb')
   call object_needed ('param_pairs')
   call object_needed ('eloc_av')
   call object_needed ('deloc_av')
   call object_needed ('dpsi_eloc_av')
   call object_needed ('dpsi_eloc_c_av')
   call object_needed ('dpsi_deloc_c_av')
   call object_needed ('dpsi_dpsi_eloc_av')
   call object_needed ('dpsi_av')

   return

  endif

! begin

! For approximate Hamitonian for orbitals
  if (l_opt_orb_eig .or. l_opt_orb_orb_eig) then
   call object_provide (lhere, 'ovlp_lin')
   call object_provide (lhere, 'delta_eps')
  endif

! allocations
  call object_alloc ('ham_lin_energy', ham_lin_energy, param_aug_nb, param_aug_nb)

! first element
  ham_lin_energy (1,1) = eloc_av

! first row and first column
  do i = 1, param_nb
    ham_lin_energy (i+1,1) = dpsi_eloc_c_av (i)

!   approximate Hamiltonian for orbitals
    if (l_opt_orb_eig .and. i > nparmcsf+nparmj) then
!     ham_lin_energy (1,i+1) = ovlp_lin(1,i+1) * (eloc_av + delta_eps (i-nparmcsf-nparmj))
     ham_lin_energy (1,i+1) = ham_lin_energy (i+1,1)  ! symmetric for orbitals

!   normal Hamiltoniam
    else
     ham_lin_energy (1,i+1) = dpsi_eloc_c_av (i) + deloc_av (i)
    endif

  enddo


! derivative-derivative part
  do j = 1, param_nb
   do i = 1, param_nb
     pair = param_pairs (i,j)

!   approximate Hamiltonian for Jastrow-orbital, CSF-orbital mixed terms (swap i and j)
    if (l_opt_orb_eig .and. i <= nparmcsf+nparmj .and. j > nparmcsf+nparmj) then
     ham_lin_energy (i+1,j+1) =  ham_lin_energy (j+1,i+1)

!   diagonal-only approximation for orbital-orbital block
    elseif (l_opt_orb_orb_diag .and. i > nparmcsf+nparmj .and. j > nparmcsf+nparmj .and. i /= j ) then
     ham_lin_energy (i+1,j+1) = 0.d0

!   approximate Hamiltonian for orbital-orbital terms
    elseif ((l_opt_orb_eig .or. l_opt_orb_orb_eig) .and. i > nparmcsf+nparmj .and. j > nparmcsf+nparmj) then
!     ham_lin_energy (i+1,j+1) = ovlp_lin (i+1,j+1) * (eloc_av + delta_eps (j-nparmcsf-nparmj))
     ham_lin_energy (i+1,j+1) = ovlp_lin (i+1,j+1) * ((eloc_av + delta_eps (j-nparmcsf-nparmj)) + (eloc_av + delta_eps (i-nparmcsf-nparmj)))/2.d0

!   normal Hamiltoniam
    else

     ham_lin_energy (i+1,j+1) = dpsi_deloc_c_av (i, j) + dpsi_dpsi_eloc_av (pair)     &
                       - dpsi_av (j) * dpsi_eloc_av (i) - dpsi_av (i) * dpsi_eloc_av (j) &
                       + dpsi_av (i) * dpsi_av (j) * eloc_av
    endif

!   if (i /= j) then
!
!!   approximate Hamiltonian for Jastrow-orbital, CSF-orbital and orbital-orbital terms
!    if (l_opt_orb_eig .and. i > nparmcsf+nparmj) then
!     ham_lin_energy (j+1,i+1) = ovlp_lin (j+1,i+1) * (eloc_av + delta_eps (i-nparmcsf-nparmj))
!
!!   normal Hamiltoniam
!    else
!     ham_lin_energy (j+1,i+1) = dpsi_deloc_c_av (j, i) + dpsi_dpsi_eloc_av (pair)     &
!                       - dpsi_av (i) * dpsi_eloc_av (j)                     &
!                       - dpsi_av (j) * ( dpsi_eloc_av (i)) &
!                       + dpsi_av (i) * dpsi_av (j) * eloc_av
!    endif
!
!   endif !i /= j

   enddo
  enddo

! symmetrize Hamiltonian
  if (l_sym_ham) then
   ham_lin_energy = (ham_lin_energy + transpose(ham_lin_energy))/2.d0
  endif


!  write(6,*)
!  write(6,'(a)') 'Hamiltonian matrix:'
!  do i = 1, param_aug_nb
!    write(6,'(100e16.8)') (ham_lin_energy(i,j),j=1,param_aug_nb)
!  enddo

!  call object_provide ('ovlp_lin')
!  write(6,'(2a,100f12.4)') trim(here),': ham_lin_energy diagonal=',(ham_lin_energy(i,i)/ovlp_lin(i,i),i=1,param_aug_nb)

  end subroutine ham_lin_energy_bld

! ==============================================================================
  subroutine ham_lin_variance_bld
! ------------------------------------------------------------------------------
! Description   : Hamiltonian matrix over basis Psi(i) for linear method
! Description   : for variance minimization
!
! Created       : J. Toulouse, 25 Apr 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  integer i, j

! header
  if (header_exe) then

   call object_create ('ham_lin_variance')

   call object_needed ('param_nb')
   call object_needed ('param_aug_nb')
   call object_needed ('eloc_var')
   call object_needed ('gradient_variance')
   call object_needed ('hessian_variance')
   call object_needed ('dpsi_dpsi_c_av')

   return

  endif

! begin

! allocations
  call object_alloc ('ham_lin_variance', ham_lin_variance, param_aug_nb, param_aug_nb)

! first element
  ham_lin_variance (1,1) = eloc_var

! first row and first column
  do i = 1, param_nb
    ham_lin_variance (i+1,1) = gradient_variance (i)/2.d0
    ham_lin_variance (1,i+1) = ham_lin_variance (i+1,1)
  enddo

! derivative-derivative part
  do j = 1, param_nb
   do i = 1, param_nb
     ham_lin_variance (i+1,j+1) = hessian_variance (i,j)/2.d0 + eloc_var * dpsi_dpsi_c_av (i, j)
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
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'ham_lin_bld'
  integer i, j

! header
  if (header_exe) then

   call object_create ('ham_lin')

   call object_needed ('param_nb')
   call object_needed ('param_aug_nb')
   call object_needed ('ham_lin_energy')
   call object_needed ('p_var')

   return

  endif

! begin

! allocations
  call object_alloc ('ham_lin', ham_lin, param_aug_nb, param_aug_nb)

  if (p_var /= 0) then
    call object_provide (lhere, 'ham_lin_variance')
    ham_lin (:,:) = (1.d0 - p_var) * ham_lin_energy (:,:) + p_var * ham_lin_variance (:,:)
  else
    ham_lin (:,:) = ham_lin_energy (:,:)
  endif

  end subroutine ham_lin_bld

! ==============================================================================
  subroutine ham_lin_renorm_bld
! ------------------------------------------------------------------------------
! Description   : Renormalized hamiltonian matrix for linear method
!
! Created       : J. Toulouse, 13 Jul 2007
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'ham_lin_renorm_bld'
  integer i, j

! header
  if (header_exe) then

   call object_create ('ham_lin_renorm')

   call object_needed ('param_aug_nb')
   call object_needed ('ham_lin')
   call object_needed ('renorm_vector')

   return

  endif

! begin

! allocations
  call object_alloc ('ham_lin_renorm', ham_lin_renorm, param_aug_nb, param_aug_nb)

! renormalizing Hamiltonian matrix
  do i = 1, param_aug_nb
    do j = 1, param_aug_nb
      ham_lin_renorm (i,j) = ham_lin (i,j) / (renorm_vector (i) * renorm_vector (j))
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
  implicit none

! local
  character(len=max_string_len_rout), save :: here = 'ham_lin_renorm_stab_bld'
  integer i

! header
  if (header_exe) then

   call object_create ('ham_lin_renorm_stab')

   call object_needed ('param_nb')
   call object_needed ('param_aug_nb')
   call object_needed ('ham_lin_renorm')
   call object_needed ('diag_stab')

   return

  endif

! begin
! allocations
  call object_alloc ('ham_lin_renorm_stab', ham_lin_renorm_stab, param_aug_nb, param_aug_nb)

  ham_lin_renorm_stab (:,:) = ham_lin_renorm (:,:)

! stabilization by adding overlap matrix
  if (trim(stabilization) == 'overlap') then
  call object_provide ('ovlp_lin')
  do i = 1, param_nb
      ham_lin_renorm_stab (i+1,i+1) = ham_lin_renorm (i+1,i+1) + diag_stab * ovlp_lin (i+1,i+1)
  enddo

! stabilization by adding identity matrix
  else
  do i = 1, param_nb
    ham_lin_renorm_stab (i+1,i+1) = ham_lin_renorm (i+1,i+1) + diag_stab
  enddo
  endif

  end subroutine ham_lin_renorm_stab_bld

! ==============================================================================
  subroutine delta_lin_bld
! ------------------------------------------------------------------------------
! Description   : variation of parameters for linear method
!
! Created       : J. Toulouse, 10 Jan 2006
! ------------------------------------------------------------------------------
  implicit none
  include 'commons.h'

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_lin_bld'
  integer i, j
  integer ex_i, iparmcsf, iparmj, iparm, jparm
  integer lwork, info
  real(dp), allocatable :: mat_a (:,:), mat_b (:,:)
  real(dp), allocatable :: eigvec (:,:)
  real(dp), allocatable :: eigval_r (:), eigval_i (:), eigval_denom (:)
  real(dp), allocatable :: work (:)
  real(dp) eigvec_1st_component_max, eigval_r_min, eigvec_first_coef, eigval_lowest
  integer eig_ind, eig_1st_component_max_ind, eigval_lowest_ind
  integer, allocatable :: eigval_srt_ind_to_eigval_ind (:), eigval_ind_to_eigval_srt_ind (:)
  integer temp
  logical target_state_found


! header
  if (header_exe) then

   call object_create ('delta_lin')

   call object_needed ('param_aug_nb')
   call object_needed ('param_nb')
   call object_needed ('ovlp_lin')
   call object_needed ('ham_lin_renorm_stab')
   call object_needed ('dpsi_av')
   call object_needed ('ovlp_lin_renorm')
   call object_needed ('renorm_vector')

   return

  endif

! begin

! allocation
  call object_alloc ('delta_lin', delta_lin, param_nb)

! temprorary arrays
  call alloc ('mat_a', mat_a, param_aug_nb, param_aug_nb)
  call alloc ('mat_b', mat_b, param_aug_nb, param_aug_nb)
  call alloc ('eigvec', eigvec, param_aug_nb, param_aug_nb)
  call alloc ('eigval_r', eigval_r, param_aug_nb)
  call alloc ('eigval_i', eigval_i, param_aug_nb)
  call alloc ('eigval_denom', eigval_denom, param_aug_nb)

  mat_a (:,:) = ham_lin_renorm_stab (:,:)
  mat_b (:,:) = ovlp_lin_renorm (:,:)

! calculate optimal value of lwork
  lwork = 1
  call alloc ('work', work, lwork)
  call dggev('N','V',param_aug_nb, mat_a, param_aug_nb, mat_b, param_aug_nb, eigval_r, eigval_i,  &
             eigval_denom, eigvec, param_aug_nb, eigvec, param_aug_nb, work, -1, info)
  if (info /= 0) then
   call die (lhere, 'problem in dggev (while calculating optimal value of lwork): info='+info+' /= 0')
  endif
  lwork =  work(1)
  call alloc ('work', work, lwork)

! solve generalized eigenvalue problem A*x = lambda*B*x
  write(6,*)
  write(6,'(a,1pd9.1)') 'Solving generalized eigenvalue equation of linear method with a_diag =', diag_stab
  call dggev('N','V',param_aug_nb, mat_a, param_aug_nb, mat_b, param_aug_nb, eigval_r, eigval_i,  &
             eigval_denom, eigvec, param_aug_nb, eigvec, param_aug_nb, work, lwork, info)
  if (info /= 0) then
   call die (lhere, 'problem in dggev: info='+info+' /= 0')
  endif

!   write(6,'(2a,100d10.2)') trim(lhere),': eigval_denom=', eigval_denom (:)

! calculate eigenvalue
  do i = 1, param_aug_nb
    eigval_r (i) = eigval_r (i) / eigval_denom (i)
    eigval_i (i) = eigval_i (i) / eigval_denom (i)
  enddo

! Sorting out eigenvalues
! eigval_srt_ind_to_eigval_ind is the map from sorted eigenvalues to original eigenvalues
  call alloc ('eigval_srt_ind_to_eigval_ind', eigval_srt_ind_to_eigval_ind, param_aug_nb)
  do i = 1, param_aug_nb
    eigval_srt_ind_to_eigval_ind (i) = i
  enddo
  do i = 1, param_aug_nb
    do j = i+1, param_aug_nb
      if (eigval_r (eigval_srt_ind_to_eigval_ind (j)) < eigval_r (eigval_srt_ind_to_eigval_ind (i))) then
        temp = eigval_srt_ind_to_eigval_ind (i)
        eigval_srt_ind_to_eigval_ind (i) = eigval_srt_ind_to_eigval_ind (j)
        eigval_srt_ind_to_eigval_ind (j) = temp
      endif
    enddo
  enddo
! eigval_ind_to_eigval_srt_ind is the map from original eigenvalues to sorted eigenvalues
  call alloc ('eigval_ind_to_eigval_srt_ind', eigval_ind_to_eigval_srt_ind, param_aug_nb)
  do i = 1, param_aug_nb
   eigval_ind_to_eigval_srt_ind (eigval_srt_ind_to_eigval_ind (i)) = i
  enddo

! print eigenvalues
  write(6,'(a)') 'Sorted (complex) eigenvalues:'
  do i = 1, param_aug_nb
    write(6,'(a,i5,a,f12.6,a,f12.6,a)') 'eigenvalue #',i,': ',eigval_r (eigval_srt_ind_to_eigval_ind (i)), ' +', eigval_i (eigval_srt_ind_to_eigval_ind (i)),' i'
  enddo

! print eigenvectors
!  write(6,'(a)') 'Eigenvectors:'
!  do j = 1, param_aug_nb
!    write(6,'(a,i3,a,100f12.6)') 'eigenvector # ', j,' :', (eigvec (i, eigval_srt_ind_to_eigval_ind (j)), i = 1, param_aug_nb)
!  enddo

! Find eigenvector with largest first coefficient
  eigvec_1st_component_max = dabs(eigvec (1,1))
  eig_1st_component_max_ind = 1
  do i = 1, param_aug_nb
    if (dabs(eigvec (1,i)) > eigvec_1st_component_max) then
      eigvec_1st_component_max = dabs(eigvec (1,i))
      eig_1st_component_max_ind = i
    endif
    if (dabs(eigvec (1,i)) == eigvec_1st_component_max) then
      if (eigval_r (i) < eigval_r (eig_1st_component_max_ind)) then
       eigvec_1st_component_max = dabs(eigvec (1,i))
       eig_1st_component_max_ind = i
      endif
    endif
  enddo
  write(6,'(a,i5,a,f12.6,a,f12.6,a)') 'The (sorted) eigenvector with largest first coefficient is #',eigval_ind_to_eigval_srt_ind (eig_1st_component_max_ind), ': ',eigval_r (eig_1st_component_max_ind), ' +', eigval_i (eig_1st_component_max_ind),' i'

! Find eigenvector with lowest eigenvalue that is not crazy
  eigval_lowest = 9.d99
  eigval_lowest_ind = 0
  do i = 1, param_aug_nb
    if (eigval_r (i) < eigval_lowest .and. dabs(eigval_r (i)-etrial) < 10.d0) then
      eigval_lowest = eigval_r (i)
      eigval_lowest_ind = i
    endif
    if (eigval_r (i) == eigval_lowest) then
      if (dabs(eigvec (1,i)) > dabs(eigvec (1, eigval_lowest_ind))) then
        eigval_lowest = eigval_r (i)
        eigval_lowest_ind = i
      endif
    endif
  enddo

  if (eigval_lowest_ind /= 0) then
    write(6,'(a,i5,a,f12.6,a,f12.6,a)') 'The (sorted) eigenvector with lowest reasonable eigenvalue is #',eigval_ind_to_eigval_srt_ind (eigval_lowest_ind), ': ',eigval_r (eigval_lowest_ind), ' +', eigval_i (eigval_lowest_ind),' i'
  else
!   if no reasonable lowest eigenvalue found, then just take the lowest one
    eigval_lowest = 9.d99
    eigval_lowest_ind = 0
    do i = 1, param_aug_nb
      if (eigval_r (i) < eigval_lowest) then
        eigval_lowest = eigval_r (i)
        eigval_lowest_ind = i
      endif
      if (eigval_r (i) == eigval_lowest) then
        if (dabs(eigvec (1,i)) > dabs(eigvec (1, eigval_lowest_ind))) then
          eigval_lowest = eigval_r (i)
          eigval_lowest_ind = i
        endif
      endif
    enddo
    write(6,'(a,i5,a,f12.6,a,f12.6,a)') 'The (sorted) eigenvector with lowest eigenvalue is #',eigval_ind_to_eigval_srt_ind (eigval_lowest_ind), ': ',eigval_r (eigval_lowest_ind), ' +', eigval_i (eigval_lowest_ind),' i'
    write(6,'(a)') 'Warning: all the eigenvalues are outside the resonable energy windows!'
  endif

! if target_state = 0, select eigenvector with lowest reasonable eigenvalue or largest 1st components
  if (target_state == 0) then

   if (l_select_eigvec_lowest) then
     eig_ind = eigval_lowest_ind
   elseif (l_select_eigvec_largest_1st_coef) then
     eig_ind = eig_1st_component_max_ind
   else
     call die (lhere, 'Both select_eigvec_lowest and select_eigvec_largest_1st_coef are false.')
   endif

! if target_state >= 1, simply select the corresponding the wanted target state
  else

    if (target_state > param_aug_nb) then
     call die (lhere, 'target_state='+target_state+' > param_aug_nb='+param_aug_nb)
    endif
    eig_ind = eigval_srt_ind_to_eigval_ind (target_state)

  endif

  write(6,'(a,i5,a,f12.6,a,f12.6,a)') 'The selected (sorted) eigenvector is #',eigval_ind_to_eigval_srt_ind (eig_ind), ': ',eigval_r (eig_ind), ' +', eigval_i (eig_ind),' i'
!  write(6,'(2a,100f12.7)') trim(lhere),': parameters variations =', eigvec_lin (:)

! undo renormalization
! warning: only for selected eigenvector
  eigvec (:, eig_ind) = eigvec (:, eig_ind) / renorm_vector (:)

! normalize eigenvector so that first component is 1
  eigvec(:,eig_ind) = eigvec(:,eig_ind)/eigvec(1,eig_ind)

! norm of linear wave function variation for nonlinear parameter
  psi_lin_var_norm = 0.d0
  do iparm = nparmcsf+1, param_nb
   do jparm = nparmcsf+1, param_nb
     psi_lin_var_norm = psi_lin_var_norm + eigvec(1+iparm,eig_ind)*eigvec(1+jparm,eig_ind)*ovlp_lin(1+iparm,1+jparm)
   enddo
  enddo
  write (6,'(a,f)') 'Norm of linear wave function variation for nonlinear parameters =', psi_lin_var_norm


! calculate the actual parameter variations
  eigvec_first_coef = eigvec(1,eig_ind)

  select case(trim(update_nonlinear))

! original: come back to original derivatives basis for all parameters
   case ('original')
    do i = 1, param_nb
       eigvec_first_coef = eigvec_first_coef - eigvec(1+i,eig_ind) * dpsi_av (i)
    enddo

! semiorthogonal: use semiorthognal derivatives for nonlinear parameters
   case ('semiorthogonal')

!   come back to original derivatives for the CSFs only
    do iparmcsf = 1, nparmcsf
       eigvec_first_coef = eigvec_first_coef - eigvec(1+iparmcsf,eig_ind) * dpsi_av (iparmcsf)
    enddo

!   nonlinear paramater
    eigvec_first_coef = eigvec_first_coef + (1.d0-xi)*psi_lin_var_norm/((1.d0-xi) + xi*(1.d0+psi_lin_var_norm))

  case default
    call die (lhere, 'unknown update choice >'+trim(update_nonlinear)+'<')
  end select

! final paramter variations
  do iparm = 1, param_nb
    delta_lin (iparm) = eigvec (1+iparm, eig_ind) / eigvec_first_coef
  enddo

  end subroutine delta_lin_bld

end module opt_lin_mod
