module opt_ovlp_fn_mod

  use all_tools_mod
  use jastrow_mod
  use csfs_mod
  use opt_common_mod
  use deriv_mod

! Declaration of global variables and default values
  character(len=max_string_len)   :: update_nonlinear = 'semiorthogonal'
  real(dp)                        :: xi = 1.d0
  logical                         :: l_weighted_overlap = .true.
!  real(dp), allocatable           :: csf_over_psit_j (:)
!  real(dp), allocatable           :: csf_over_psit_j_av (:)
  real(dp), allocatable           :: delta_ovlp_fn (:)
!  real(dp), allocatable           :: delta_ovlp_fn_exact (:)
  real(dp), allocatable           :: ovlp_ovlp_fn (:,:)
  real(dp), allocatable           :: ovlp_ovlp_fn_inv (:,:)
!JT  real(dp), allocatable           :: ovlp_ovlp_fn_av (:,:)
  real(dp)                        :: wt_lambda = 1.d0
  real(dp)                        :: ovlp_trial_fn
  real(dp)                        :: ovlp_trial_fn_over_ovlp_trial

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
   write(6,'(a)') '   update_nonlinear = [original|semiorthogonal] : default=semiorthogonal, choice of update of nonlinear paramaters'
   write(6,'(a)') '   xi = [real] : update of nonlinear paramaters by orthogonalization to xi Psi_0 +(1-xi) Psi_lin'
   write(6,'(a)') '                - xi=1: orthogonalization to Psi_0 (default)'
   write(6,'(a)') '                - xi=0: orthogonalization to Psi_lin, ensures min |Psi_lin-Psi_0|'
   write(6,'(a)') '                - xi=0.5: orthogonalization to Psi_0 + Psi_lin, ensures |Psi_0|=|Psi_lin|'
   write(6,'(a)') ' end'

  case ('weighted_overlap')
   call get_next_value (l_weighted_overlap)

  case ('wt_lambda')
   call get_next_value( wt_lambda)

  case ('update_nonlinear')
   call get_next_value (update_nonlinear)

  case ('xi')
   call get_next_value (xi)

  case ('end')
   exit

  case default
   call die(lhere, 'unknown word >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  call require (lhere, 'wt_lambda >= 0 and wt_lambda =< 1', wt_lambda >= 0.d0 .and. wt_lambda <= 1.d0)

  end subroutine opt_ovlp_fn_menu


!!===========================================================================
!  subroutine delta_ovlp_fn_bld_old
!! ------------------------------------------------------------------------------
!! Description   : variation of parameters for overlap FN method
!!
!! Created       : Cyrus Umrigar and Frank Petruzielo, 7 Jun 2010
!! ------------------------------------------------------------------------------
!  include 'modules.h'
!  implicit none
!
!! local
!  character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_bld'
!  integer iparm, jparm, ierr, iparmcsf, jparmcsf
!  real(dp), allocatable   :: ovlp_ovlp_fn_tmp(:,:), dpsi_tmp(:)
!  real(dp) :: psi_lin_var_norm, dpsi_tmp_first_coef, d, normalization, denominator
!
!! header
!  if(header_exe) then
!
!   call object_create('delta_ovlp_fn')
!
!   call object_needed('dpsi_av')
!   call object_needed('param_nb')
!   call object_needed('ovlp_ovlp_fn_av')
!
!   return
!
!  endif
!
!! begin
!
!! allocation
!  call object_alloc('delta_ovlp_fn', delta_ovlp_fn, param_nb)
!  call alloc('ovlp_ovlp_fn_tmp', ovlp_ovlp_fn_tmp, param_nb+1, param_nb+1)
!  call alloc('dpsi_tmp', dpsi_tmp, param_nb+1)
!!  allocate(ovlp_ovlp_fn_tmp(param_nb,param_nb),stat=ierr)
!!  if(ierr.ne.0) call die(lhere,'failed to allocate ovlp_ovlp_fn_tmp')
!
!  do iparm=1,param_nb+1
!    write(6,'(''ovlp='',9d12.4)') (ovlp_ovlp_fn_av(iparm,jparm),jparm=1,param_nb+1)
!  enddo
!
!  ovlp_ovlp_fn_tmp = ovlp_ovlp_fn_av
!  dpsi_tmp(1) = wgcum(1)/(nstep*nblk*nwalk*nproc)
!  dpsi_tmp(2:) = dpsi_av
!  write(6,'(''dpsi_tmp before='',100es12.4)') dpsi_tmp
!
!!   do iparm=1,param_nb
!!      write(6,'(''ovlp='', 9d12.4)') (ovlp_ovlp_fn_tmp(iparm,jparm),jparm=1,param_nb)
!!    enddo
!
!! Do cholesky decomposition
!  call chlsky(ovlp_ovlp_fn_tmp,param_nb+1,param_nb+1,ierr)
!  if(ierr.ne.0) call die(lhere,'cholesky decomposition failed')
!
!! Symmetrize decomposed matrix (needs to be done before calling uxb
!! or need to modify uxb)
!  do iparm=1,param_nb+1
!    do jparm=iparm+1,param_nb+1
!      ovlp_ovlp_fn_tmp(iparm,jparm)=ovlp_ovlp_fn_tmp(jparm,iparm)
!    enddo
!  enddo
!
!! Solve linear equations
!  call lxb(ovlp_ovlp_fn_tmp,param_nb+1,param_nb+1,dpsi_tmp)
!  call uxb(ovlp_ovlp_fn_tmp,param_nb+1,param_nb+1,dpsi_tmp)
!
!! normalize so that first component is 1
!  write(6,'(''dpsi_tmp after='',100es12.4)') dpsi_tmp
!  write(6,'(''dpsi_tmp(1)='',9es12.4)') dpsi_tmp(1)
!  dpsi_tmp = dpsi_tmp/dpsi_tmp(1)
!  delta_ovlp_fn = dpsi_tmp(2:)
!
!! Eqn.(8) 2007PRL
!  d = 1.d0
!  do iparm = nparmcsf+1, param_nb
!     d = d + 2 * ovlp_ovlp_fn_av(1, iparm+1) * delta_ovlp_fn(iparm)
!     do jparm = nparmcsf+1, param_nb
!        d = d +  delta_ovlp_fn(iparm) * ovlp_ovlp_fn_av(iparm+1, jparm+1) * delta_ovlp_fn(jparm)
!     enddo
!  enddo
!  d = sqrt(d)
!
!! Eqn.(7) 2007PRL
!  denominator = 1.d0
!
!  do iparm = nparmcsf+1, param_nb
!     normalization = - ( xi * d * ovlp_ovlp_fn_av(1,iparm+1) + (1 - xi)*(ovlp_ovlp_fn_av(1,iparm+1) &
!          + sum(ovlp_ovlp_fn_av(iparm+1,nparmcsf+2:)*delta_ovlp_fn(nparmcsf+1:)))) &
!          / ( xi * d + (1 - xi)*(1 + sum(ovlp_ovlp_fn_av(1,nparmcsf+2:)*delta_ovlp_fn(nparmcsf+1:))))
!     denominator = denominator - normalization * delta_ovlp_fn(iparm)
!  enddo
!
!  delta_ovlp_fn = delta_ovlp_fn / denominator
!
!! ! norm of linear wave function variation for nonlinear parameter
!!   psi_lin_var_norm = 0.d0
!!   do iparm = nparmcsf+1, param_nb
!!    do jparm = nparmcsf+1, param_nb
!!      psi_lin_var_norm = psi_lin_var_norm + dpsi_tmp(1+iparm)*dpsi_tmp(1+jparm)*ovlp_ovlp_fn(1+iparm,1+jparm)
!!    enddo
!!   enddo
!!   write(6,'(a,f10.6)') 'Norm of linear wave function variation for nonlin. params =', psi_lin_var_norm
!! !  if(psi_lin_var_norm > 10*smallest_norm) write(6,'(''Warning: psi_lin_var_norm > 10*smallest_norm'',2d12.4)')  psi_lin_var_norm,smallest_norm
!
!! ! calculate the actual parameter variations
!!   dpsi_tmp_first_coef = dpsi_tmp(1)
!
!!   select case (trim(update_nonlinear))
!
!! ! original: come back to original derivatives basis for all parameters
!!    case ('original')
!!     do iparm = 1, param_nb
!!        dpsi_tmp_first_coef = dpsi_tmp_first_coef - dpsi_tmp(1+iparm) * dpsi_av(iparm)
!!     enddo
!
!! ! semiorthogonal: use semiorthognal derivatives for nonlinear parameters
!!    case ('semiorthogonal')
!
!! !   come back to original derivatives for the CSFs only
!
!!     write (6,*) "dpsi_tmp_first_coef: ", dpsi_tmp_first_coef
!!     do iparmcsf = 1, nparmcsf
!!        dpsi_tmp_first_coef = dpsi_tmp_first_coef - dpsi_tmp(1+iparmcsf) * dpsi_av(iparmcsf)
!!     enddo
!!     write (6,*) "dpsi_tmp_first_coef: ", dpsi_tmp_first_coef
!! !   nonlinear paramater
!!     dpsi_tmp_first_coef = dpsi_tmp_first_coef +(1.d0-xi)*psi_lin_var_norm/((1.d0-xi) + xi*(1.d0+psi_lin_var_norm))
!!     write (6,*) "dpsi_tmp_first_coef: ", dpsi_tmp_first_coef
!!     write (6,*) " "
!!    case default
!!     call die(lhere, 'unknown update choice >'+trim(update_nonlinear)+'<')
!
!!   end select
!
!! final parameter variations
!!  delta_ovlp_fn = dpsi_tmp(2:) / dpsi_tmp_first_coef
!
!  !temporary
!
!
!!  call release ('dpsi_tmp', dpsi_tmp)
!
!  end subroutine delta_ovlp_fn_bld_old
!
!! ==============================================================================
!  subroutine ovlp_ovlp_fn_bld_old
!! ------------------------------------------------------------------------------
!! Description   : overlap matrix for ovlp_fn method
!!                 The local object is the same as in the linear method but the average here is performed over Psi(0)^2
!! Description   : < Psi(i)  Psi(j) > / <Psi(0) Psi(0) >_|Psi(0)|^2
!!
!! Created       : Cyrus Umrigar and Frank Petruzielo by modifying ovlp_lin_bld of J. Toulouse. 8 Jun 2010
!! ------------------------------------------------------------------------------
!  include 'modules.h'
!  implicit none
!
!! local
!  integer i, j
!!  real(dp), allocatable           :: ovlp_ovlp_fn_eigvec(:,:)
!!  real(dp), allocatable           :: ovlp_ovlp_fn_eigval(:)
!
!! header
!  if(header_exe) then
!
!   call object_create('ovlp_ovlp_fn')
!   call object_create('ovlp_ovlp_fn_av')
!!   call object_average_define ('ovlp_ovlp_fn','ovlp_ovlp_fn_av')
!!   call object_create('ovlp_ovlp_fn_eigvec')
!!   call object_create('ovlp_ovlp_fn_eigval')
!
!   call object_needed('param_pairs')
!   call object_needed('param_nb')
!   call object_needed('dpsi_dpsi')
!   call object_needed('dpsi')
!!   call object_needed('psi_jas')
!
!   return
!
!  endif
!
!! begin
!
!! allocations
!  call object_alloc('ovlp_ovlp_fn', ovlp_ovlp_fn, param_nb+1, param_nb+1)
!  call object_alloc('ovlp_ovlp_fn_av', ovlp_ovlp_fn_av, param_nb+1, param_nb+1)
!
!  ovlp_ovlp_fn(1,1) = 1.d0
!  ovlp_ovlp_fn(1,2:) = dpsi
!  ovlp_ovlp_fn(2:,1) = dpsi
!!   do i = 1, param_nb
!!   ovlp_ovlp_fn(1,i+1) = 0.d0
!!   ovlp_ovlp_fn(i+1,1) = 0.d0
!!  enddo
!
!! derivative-derivative part
!  do i = 1, param_nb
!   do j = i, param_nb
!
!     ovlp_ovlp_fn(i+1,j+1) = dpsi_dpsi(param_pairs(i,j))
!
!!   symmetrize overlap matrix
!    if (i /= j) then
!     ovlp_ovlp_fn(j+1,i+1) = ovlp_ovlp_fn(i+1,j+1)
!    endif
!
!   enddo
!  enddo
!
!!   do i=1,param_nb+1
!!      write(6,'(''ovlp_local='',9d12.4)') (ovlp_ovlp_fn(i,j),j=1,param_nb+1)
!!   enddo
!
!  end subroutine ovlp_ovlp_fn_bld_old
!
! ==============================================================================
  subroutine ovlp_trial_fn_bld
! ------------------------------------------------------------------------------
! Description   : overlap of trial wavefunction with fn wavefunction
!
! Description   : < Psi(0) Psi_FN > / sqrt(<Psi(0) Psi(0) > <Psi_FN Psi_FN>
!
! Created       : Cyrus Umrigar and Frank Petruzielo  1 Jul 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer i, j

! header
  if(header_exe) then

   call object_create('ovlp_trial_fn')

   call object_needed('wgcum1')
   call object_needed('wgcm21')

   return

  endif

! begin
  ovlp_trial_fn = wgcum1(1) / sqrt(wgcm21(1) *(nstep*nblk*nwalk*nproc))

  end subroutine ovlp_trial_fn_bld

! ==============================================================================
  subroutine ovlp_trial_fn_over_ovlp_trial_bld
! ------------------------------------------------------------------------------
! Description   : overlap of trial wavefunction with fn wavefunction
! Description   : normalized to < Psi(0)|Psi(0) >
!
! Description   : < Psi(0)|Psi_FN > / < Psi(0)|Psi(0) >
!
! Created       : J. Toulouse, 8 Jul 2010
! ------------------------------------------------------------------------------
  implicit none

! header
  if(header_exe) then

   call object_create ('ovlp_trial_fn_over_ovlp_trial')

   call object_needed ('walker_weights_sum')
   call object_needed ('total_iterations_nb')

   return

  endif

! begin
  ovlp_trial_fn_over_ovlp_trial = walker_weights_sum / total_iterations_nb

  end subroutine ovlp_trial_fn_over_ovlp_trial_bld

!===========================================================================
  subroutine delta_ovlp_fn_bld
! ------------------------------------------------------------------------------
! Description   : variation of parameters for overlap FN method
!
! Created       : J. Toulouse, 7 Jun 2010
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_bld'
  integer param_i, param_j
  real(dp) transform_factor

! header
  if (header_exe) then

   call object_create ('delta_ovlp_fn')

   call object_needed ('param_nb')
   call object_needed ('dpsi_av')
   call object_needed ('dpsi_uwav')
   call object_needed ('dpsi_dpsi_uwcovar')
   call object_needed ('dpsi_dpsi_uwcovar_inv')

   return

  endif

! allocation
  call object_alloc ('delta_ovlp_fn', delta_ovlp_fn, param_nb)

  write(6,'(a)') "Projection of wave function derivatives unto fixed-node and trial wave functions:"
  do param_i = 1, param_nb
   write(6,'(a,i5,a,f12.6,a,f12.6,a,f12.6)') "parameter # ",param_i, ": <Psi_i|Psi_FN>= ", dpsi_av (param_i), ", <Psi_i|Psi_0>= ", dpsi_uwav (param_i), ", difference= ", dpsi_av (param_i) - dpsi_uwav (param_i)
  enddo


! find parameter variations in semiorthogonalized basis
  do param_i = 1, param_nb
   delta_ovlp_fn (param_i) = 0.d0
   do param_j = 1, param_nb
    delta_ovlp_fn (param_i) = delta_ovlp_fn (param_i) + dpsi_dpsi_uwcovar_inv (param_i, param_j) * (dpsi_av (param_j) - dpsi_uwav (param_j))
   enddo
  enddo

! norm of linear wave function variation for nonlinear parameter
  call object_provide ('nparmcsf')
  psi_lin_var_norm = 0.d0
  do param_i = nparmcsf+1, param_nb
   do param_j = nparmcsf+1, param_nb
     psi_lin_var_norm = psi_lin_var_norm + delta_ovlp_fn (param_i) * delta_ovlp_fn (param_j) * dpsi_dpsi_uwcovar (param_i, param_j)
   enddo
  enddo
  write(6,'(a,f10.6)') 'Norm of linear wave function variation for nonlinear parameters =', psi_lin_var_norm

! transformation of parameter variations to desired basis
  transform_factor = 1.d0
  select case (trim(update_nonlinear))

! original: come back to original derivatives basis for all parameters
   case ('original')
    do param_i = 1, param_nb
       transform_factor = transform_factor - dpsi_av (param_i) * delta_ovlp_fn (param_i)
    enddo

! semiorthogonal: use semiorthognal derivatives for nonlinear parameters
   case ('semiorthogonal')

!   come back to original derivatives for the CSF parameters only
    do param_i = 1, nparmcsf
       transform_factor = transform_factor - dpsi_av (param_i) * delta_ovlp_fn (param_i)
    enddo

!   nonlinear parameters
    transform_factor = transform_factor + (1.d0-xi)*psi_lin_var_norm/((1.d0-xi) + xi*(1.d0+psi_lin_var_norm))

  case default
    call die (lhere, 'unknown update choice >'+trim(update_nonlinear)+'<')
  end select

! final parameter variations
  delta_ovlp_fn (:) = delta_ovlp_fn (:) / transform_factor

  end subroutine delta_ovlp_fn_bld


! ! ==============================================================================
!   subroutine csf_over_psit_j_bld
! ! ------------------------------------------------------------------------------
! ! Description   : csf divided by Psi_T times Jastrow^2
! !
! ! Created       : Cyrus Umrigar and Frank Petruzielo, 7 Jun 2010
! ! ------------------------------------------------------------------------------
!   include 'modules.h'
!   implicit none

! ! local
!   integer csf_i

! ! header
!   if(header_exe) then

!    call object_create('csf_over_psit_j')
!    call object_average_define ('csf_over_psit_j', 'csf_over_psit_j_av')

!    call object_needed('ncsf')
!    call object_needed('csf_over_psid')
!    call object_needed('psi_jas')

!    return

!   endif

! ! begin

! ! allocations
!   call object_alloc('csf_over_psit_j', csf_over_psit_j, ncsf)
!   call object_alloc('csf_over_psit_j_av', csf_over_psit_j_av, ncsf)

!   do csf_i = 1, ncsf
!     csf_over_psit_j(csf_i) = csf_over_psid (csf_i)/psi_jas**2
!   enddo ! csf_i

!   end subroutine csf_over_psit_j_bld

! ! ==============================================================================
!   subroutine delta_ovlp_fn_bld
! ! ------------------------------------------------------------------------------
! ! Description   : variation of parameters for overlap FN method
! !
! ! Created       : Cyrus Umrigar and Frank Petruzielo, 7 Jun 2010
! ! ------------------------------------------------------------------------------
!   include 'modules.h'
!   implicit none

! ! local
!   character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_bld'
!   integer iparm

! ! header
!   if(header_exe) then

!    call object_create('delta_ovlp_fn')

!    call object_needed('csf_over_psit_j_av')
!    call object_needed ('param_nb')

!    return

!   endif

! ! begin

! ! allocation
!   call object_alloc('delta_ovlp_fn', delta_ovlp_fn, param_nb)

!   call object_provide ('csf_coef')
!   call object_provide ('delta_ovlp_fn_exact')

!   csf_over_psit_j_av=csf_over_psit_j_av/csf_over_psit_j_av(1)
!   do iparm = 1, param_nb
!      write(6,*) "temp: csf_over_psit_j_av", csf_over_psit_j_av(iparm)
!      delta_ovlp_fn(iparm) = csf_over_psit_j_av(iparm) - csf_coef(iparm, 1)
!      write(6,*) "temp: :",csf_coef(iparm, 1), delta_ovlp_fn(iparm)
!   enddo

!   end subroutine delta_ovlp_fn_bld

! ! ==============================================================================
!   subroutine delta_ovlp_fn_exact_bld
! ! ------------------------------------------------------------------------------
! ! Description   : variation of parameters for overlap FN method
! !
! ! Created       : Cyrus Umrigar and Frank Petruzielo, 7 Jun 2010
! ! ------------------------------------------------------------------------------
!   include 'modules.h'
!   implicit none

! ! local
!   character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_exact_bld'
!   integer iparm, jparm, ierr
!   real(dp), allocatable   :: ovlp_ovlp_fn_tmp(:,:)

! ! header
!   if(header_exe) then

!    call object_create('delta_ovlp_fn_exact')

!    call object_needed('csf_over_psit_j_av')
!    call object_needed('param_nb')
!    call object_needed('ovlp_ovlp_fn_av')

!    return

!   endif

! ! begin

! ! allocation
!   call object_alloc('delta_ovlp_fn_exact', delta_ovlp_fn_exact, param_nb)
!   allocate(ovlp_ovlp_fn_tmp(param_nb,param_nb),stat=ierr)
!   if(ierr.ne.0) call die(lhere,'failed to allocate ovlp_ovlp_fn_tmp')

!   ovlp_ovlp_fn_tmp=ovlp_ovlp_fn_av
!   delta_ovlp_fn_exact=csf_over_psit_j_av

!   call object_provide ('csf_coef')

!   do iparm=1,param_nb
!      write(6,'(''ovlp='', 9d12.4)') (ovlp_ovlp_fn_tmp(iparm,jparm),jparm=1,param_nb)
!    enddo

! ! Do cholesky decomposition
!   call chlsky(ovlp_ovlp_fn_tmp,param_nb,param_nb,ierr)
!   if(ierr.ne.0) call die(lhere,'cholesky decomposition failed')

! ! Symmetrize decomposed matrix (needs to be done before calling uxb
! ! or need to modify uxb)
!   do iparm=1,param_nb
!     do jparm=iparm+1,param_nb
!       ovlp_ovlp_fn_tmp(iparm,jparm)=ovlp_ovlp_fn_tmp(jparm,iparm)
!     enddo
!   enddo

! ! Solve linear equations
!   call lxb(ovlp_ovlp_fn_tmp,param_nb,param_nb,delta_ovlp_fn_exact)
!   call uxb(ovlp_ovlp_fn_tmp,param_nb,param_nb,delta_ovlp_fn_exact)

!   delta_ovlp_fn_exact=delta_ovlp_fn_exact/delta_ovlp_fn_exact(1)
!   do iparm = 1, param_nb
!      write(6,*) "temp: delta_ovlp_fn_exact", delta_ovlp_fn_exact(iparm)
!      delta_ovlp_fn_exact(iparm) = delta_ovlp_fn_exact(iparm) - csf_coef(iparm, 1)
!      write(6,*) "temp: :",csf_coef(iparm, 1), delta_ovlp_fn_exact(iparm)
!   enddo

!   end subroutine delta_ovlp_fn_exact_bld


! ! ==============================================================================
!   subroutine ovlp_ovlp_fn_bld
! ! ------------------------------------------------------------------------------
! ! Description   : overlap matrix for ovlp_fn method
! !  
! ! Description   : < Psi(i)  Psi(j) / J J > / <Psi(0) Psi(0) >   \approx   < Psi(i)  Psi(j) / Psi(0) Psi(0) J J >_Psi(0)^2
! !
! ! Created       : Cyrus Umrigar and Frank Petruzielo by modifying ovlp_lin_bld of J. Toulouse. 8 Jun 2010
! ! ------------------------------------------------------------------------------
!   include 'modules.h'
!   implicit none

! ! local
!   integer i, j, ierr
!   real(dp), allocatable           :: ovlp_ovlp_fn_eigvec(:,:)
!   real(dp), allocatable           :: ovlp_ovlp_fn_eigval(:)

! ! header
!   if(header_exe) then

!    call object_create('ovlp_ovlp_fn')
!    call object_create('ovlp_ovlp_fn_av')
!    call object_average_define ('ovlp_ovlp_fn','ovlp_ovlp_fn_av')
!    call object_create('ovlp_ovlp_fn_eigvec')
!    call object_create('ovlp_ovlp_fn_eigval')

!    call object_needed('param_pairs')
!    call object_needed('param_nb')
!    call object_needed('dpsi_dpsi')
!    call object_needed('psi_jas')

!    return

!   endif

! ! begin

! ! allocations
!   call object_alloc('ovlp_ovlp_fn', ovlp_ovlp_fn, param_nb, param_nb)
!   call object_alloc('ovlp_ovlp_fn_av', ovlp_ovlp_fn_av, param_nb, param_nb)

! ! derivative-derivative part
!   do i = 1, param_nb
!    do j = i, param_nb

!      ovlp_ovlp_fn(i,j) = dpsi_dpsi(param_pairs(i,j)) / psi_jas**2

! !   force symmetrization of overlap matrix (important for numerics?)
!     if (i /= j) then
!      ovlp_ovlp_fn(j,i) = ovlp_ovlp_fn(i,j)
!     endif

!    enddo
!   enddo

! !   do i=1,param_nb
! !     write(6,'(''ovlp='', 9d12.4)') (ovlp_ovlp_fn(i,j),j=1,param_nb)
! !   enddo

! ! ! check eigenvalues
! !   allocate(ovlp_ovlp_fn_eigvec(param_nb, param_nb), stat=ierr)
! !   allocate(ovlp_ovlp_fn_eigval(param_nb), stat=ierr)
! !   call eigensystem(ovlp_ovlp_fn, ovlp_ovlp_fn_eigvec, ovlp_ovlp_fn_eigval, param_nb)

! !   write(6,*)
! !   write(6,'(a)') 'Eigenvalues of overlap matrix of current wave function and its first-order derivatives:'
! !   do i = 1, param_nb
! !     write(6,'(a,i4,a,es15.8)') 'overlap eigenvalue # ',i,': ',ovlp_ovlp_fn_eigval(i)
! !   enddo

!   end subroutine ovlp_ovlp_fn_bld

end module opt_ovlp_fn_mod
