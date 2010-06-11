module opt_ovlp_fn_mod

  use all_tools_mod
  use jastrow_mod
  use csfs_mod
!  use opt_common_mod
  use deriv_mod, only: dpsi_dpsi, param_pairs

! Declaration of global variables and default values
  logical                         :: l_weighted_overlap = .true.
  real(dp), allocatable           :: csf_over_psit_j (:)
  real(dp), allocatable           :: csf_over_psit_j_av (:)
  real(dp), allocatable           :: delta_ovlp_fn (:)
  real(dp), allocatable           :: delta_ovlp_fn_exact (:)
  real(dp), allocatable           :: ovlp_ovlp_fn (:,:)
  real(dp), allocatable           :: ovlp_ovlp_fn_av (:,:)
  real(dp)                        :: wt_lambda = 1.d0

  contains

!===========================================================================
  subroutine opt_ovlp_fn_menu
!---------------------------------------------------------------------------
! Description : menu for overlap with FN wavefn. optimization
!
! Created     : Cyrus Umrigar and Frank Petruzielo 7 Jun 2010
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'opt_ovlp_fn_menu'

! begin

! loop over menu lines
  do
  call get_next_word(word)

  select case (trim(word))
  case ('help')
   write(6,'(a)') ' HELP for overlap fixed-node optimization menu'
   write(6,'(a)') '  overlap_fn'
   write(6,'(a)') '   weighted_overlap = [logical] : have zero variance estimator by calculating the difference between the FN and the trial wavefn. (default=true)'
   write(6,'(a)') '   wt_lambda = [real] : Power to which DMC wts are raised unit physical time later. (default=1.d0)'
   write(6,'(a)') ' end'

  case ('weighted_overlap')
   call get_next_value(l_weighted_overlap)

  case ('wt_lambda')
   call get_next_value(wt_lambda)

  case ('end')
   exit

  case default
   call die(lhere, 'unknown word >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  call require (lhere, 'wt_lambda >= 0 and wt_lambda =< 1', wt_lambda >= 0.d0 .and. wt_lambda <= 1.d0)

  end subroutine opt_ovlp_fn_menu

! ==============================================================================
  subroutine csf_over_psit_j_bld
! ------------------------------------------------------------------------------
! Description   : csf divided by Psi_T times Jastrow^2
!
! Created       : Cyrus Umrigar and Frank Petruzielo, 7 Jun 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer csf_i

! header
  if(header_exe) then

   call object_create('csf_over_psit_j')
   call object_average_define ('csf_over_psit_j', 'csf_over_psit_j_av')

   call object_needed('ncsf')
   call object_needed('csf_over_psid')
   call object_needed('psi_jas')

   return

  endif

! begin

! allocations
  call object_alloc('csf_over_psit_j', csf_over_psit_j, ncsf)
  call object_alloc('csf_over_psit_j_av', csf_over_psit_j_av, ncsf)

  do csf_i = 1, ncsf
    csf_over_psit_j(csf_i) = csf_over_psid (csf_i)/psi_jas**2
  enddo ! csf_i

  end subroutine csf_over_psit_j_bld

! ==============================================================================
  subroutine delta_ovlp_fn_bld
! ------------------------------------------------------------------------------
! Description   : variation of parameters for overlap FN method
!
! Created       : Cyrus Umrigar and Frank Petruzielo, 7 Jun 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_bld'
  integer iparm

! header
  if(header_exe) then

   call object_create('delta_ovlp_fn')

   call object_needed('csf_over_psit_j_av')
   call object_needed ('param_nb')

   return

  endif

! begin

! allocation
  call object_alloc('delta_ovlp_fn', delta_ovlp_fn, param_nb)

  call object_provide ('csf_coef')
  call object_provide ('delta_ovlp_fn_exact')

  csf_over_psit_j_av=csf_over_psit_j_av/csf_over_psit_j_av(1)
  do iparm = 1, param_nb
     write(6,*) "temp: csf_over_psit_j_av", csf_over_psit_j_av(iparm)
     delta_ovlp_fn(iparm) = csf_over_psit_j_av(iparm) - csf_coef(iparm, 1)
     write(6,*) "temp: :",csf_coef(iparm, 1), delta_ovlp_fn(iparm)
  enddo

  end subroutine delta_ovlp_fn_bld

! ==============================================================================
  subroutine delta_ovlp_fn_exact_bld
! ------------------------------------------------------------------------------
! Description   : variation of parameters for overlap FN method
!
! Created       : Cyrus Umrigar and Frank Petruzielo, 7 Jun 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_exact_bld'
  integer iparm, jparm, ierr
  real(dp), allocatable   :: ovlp_ovlp_fn_tmp(:,:)

! header
  if(header_exe) then

   call object_create('delta_ovlp_fn_exact')

   call object_needed('csf_over_psit_j_av')
   call object_needed('param_nb')
   call object_needed('ovlp_ovlp_fn_av')

   return

  endif

! begin

! allocation
  call object_alloc('delta_ovlp_fn_exact', delta_ovlp_fn_exact, param_nb)
  allocate(ovlp_ovlp_fn_tmp(param_nb,param_nb),stat=ierr)
  if(ierr.ne.0) call die(lhere,'failed to allocate ovlp_ovlp_fn_tmp')

  ovlp_ovlp_fn_tmp=ovlp_ovlp_fn_av
  delta_ovlp_fn_exact=csf_over_psit_j_av

  call object_provide ('csf_coef')

  do iparm=1,param_nb
     write(6,'(''ovlp='', 9d12.4)') (ovlp_ovlp_fn_tmp(iparm,jparm),jparm=1,param_nb)
   enddo

! Do cholesky decomposition
  call chlsky(ovlp_ovlp_fn_tmp,param_nb,param_nb,ierr)
  if(ierr.ne.0) call die(lhere,'cholesky decomposition failed')

! Symmetrize decomposed matrix (needs to be done before calling uxb
! or need to modify uxb)
  do iparm=1,param_nb
    do jparm=iparm+1,param_nb
      ovlp_ovlp_fn_tmp(iparm,jparm)=ovlp_ovlp_fn_tmp(jparm,iparm)
    enddo
  enddo

! Solve linear equations
  call lxb(ovlp_ovlp_fn_tmp,param_nb,param_nb,delta_ovlp_fn_exact)
  call uxb(ovlp_ovlp_fn_tmp,param_nb,param_nb,delta_ovlp_fn_exact)

  delta_ovlp_fn_exact=delta_ovlp_fn_exact/delta_ovlp_fn_exact(1)
  do iparm = 1, param_nb
     write(6,*) "temp: delta_ovlp_fn_exact", delta_ovlp_fn_exact(iparm)
     delta_ovlp_fn_exact(iparm) = delta_ovlp_fn_exact(iparm) - csf_coef(iparm, 1)
     write(6,*) "temp: :",csf_coef(iparm, 1), delta_ovlp_fn_exact(iparm)
  enddo

  end subroutine delta_ovlp_fn_exact_bld

! ==============================================================================
  subroutine ovlp_ovlp_fn_bld
! ------------------------------------------------------------------------------
! Description   : overlap matrix for ovlp_fn method
!                 Note that unlike the linear method, here the first element is not for the trial wavefn.
! Description   : < Psi(i)  Psi(j) / J J > / <Psi(0) Psi(0) >   \approx   < Psi(i)  Psi(j) / Psi(0) Psi(0) J J >_Psi(0)^2
!
! Created       : Cyrus Umrigar and Frank Petruzielo by modifying ovlp_lin_bld of J. Toulouse. 8 Jun 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer i, j, ierr
  real(dp), allocatable           :: ovlp_ovlp_fn_eigvec(:,:)
  real(dp), allocatable           :: ovlp_ovlp_fn_eigval(:)

! header
  if(header_exe) then

   call object_create('ovlp_ovlp_fn')
   call object_create('ovlp_ovlp_fn_av')
   call object_average_define ('ovlp_ovlp_fn','ovlp_ovlp_fn_av')
   call object_create('ovlp_ovlp_fn_eigvec')
   call object_create('ovlp_ovlp_fn_eigval')

   call object_needed('param_pairs')
   call object_needed('param_nb')
   call object_needed('dpsi_dpsi')
   call object_needed('psi_jas')

   return

  endif

! begin

! allocations
  call object_alloc('ovlp_ovlp_fn', ovlp_ovlp_fn, param_nb, param_nb)
  call object_alloc('ovlp_ovlp_fn_av', ovlp_ovlp_fn_av, param_nb, param_nb)

! derivative-derivative part
  do i = 1, param_nb
   do j = i, param_nb

     ovlp_ovlp_fn(i,j) = dpsi_dpsi(param_pairs(i,j)) / psi_jas**2

!   force symmetrization of overlap matrix (important for numerics?)
    if (i /= j) then
     ovlp_ovlp_fn(j,i) = ovlp_ovlp_fn(i,j)
    endif

   enddo
  enddo

!   do i=1,param_nb
!     write(6,'(''ovlp='', 9d12.4)') (ovlp_ovlp_fn(i,j),j=1,param_nb)
!   enddo

! ! check eigenvalues
!   allocate(ovlp_ovlp_fn_eigvec(param_nb, param_nb), stat=ierr)
!   allocate(ovlp_ovlp_fn_eigval(param_nb), stat=ierr)
!   call eigensystem(ovlp_ovlp_fn, ovlp_ovlp_fn_eigvec, ovlp_ovlp_fn_eigval, param_nb)

!   write(6,*)
!   write(6,'(a)') 'Eigenvalues of overlap matrix of current wave function and its first-order derivatives:'
!   do i = 1, param_nb
!     write(6,'(a,i4,a,es15.8)') 'overlap eigenvalue # ',i,': ',ovlp_ovlp_fn_eigval(i)
!   enddo

  end subroutine ovlp_ovlp_fn_bld

end module opt_ovlp_fn_mod

