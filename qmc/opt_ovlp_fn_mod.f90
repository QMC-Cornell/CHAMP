module opt_ovlp_fn_mod

  use all_tools_mod
  use jastrow_mod
  use csfs_mod
  use opt_common_mod
  use deriv_mod
  use optim_mod

! Declaration of global variables and default values
  character(len=max_string_len)   :: update_nonlinear = 'semiorthogonal'
  real(dp)                        :: xi = 1.d0
  real(dp)                        :: weight_trial = 1.0
  logical                         :: l_opt_ovlp_fn_nonsemiortho = .false.
  logical                         :: l_opt_ovlp_fn_linear = .false.
  logical                         :: l_ovlp_identity = .false.
  logical                         :: l_ovlp_diagonal = .true.
  logical                         :: l_opt_ovlp_branching = .false.
  real(dp), allocatable           :: delta_ovlp_fn (:)
  real(dp), allocatable           :: delta_ovlp_fn_nonsemiortho (:)
  real(dp), allocatable           :: delta_ovlp_fn_linear (:)
  real(dp)                        :: ovlp_trial_fn
  real(dp)                        :: ovlp_trial_fn_sav
  real(dp)                        :: ovlp_trial_fn_over_ovlp_trial
  real(dp)                        :: ovlp_trial_fn_over_ovlp_trial_sav
  real(dp), allocatable           :: dpsi_over_jas2 (:)
  real(dp), allocatable           :: dpsi_over_jas2_av (:)
  real(dp), allocatable           :: dpsi_over_jas2_uwav (:)
  real(dp), allocatable           :: dpsi2_over_jas2 (:)
  real(dp), allocatable           :: dpsi2_over_jas2_uwav (:)
  real(dp), allocatable           :: dpsi_dpsi_over_jas2 (:,:)  !temp
  real(dp), allocatable           :: dpsi_dpsi_over_jas2_uwav (:,:) !temp
  real(dp)                        :: first_csf_over_jas2
  real(dp)                        :: first_csf_over_jas2_av
  real(dp)                        :: first_csf_over_jas2_uwav
  real(dp)                        :: first_csf2_over_jas2
  real(dp)                        :: first_csf2_over_jas2_uwav
  real(dp)                        :: one_over_jas2
  real(dp)                        :: one_over_jas2_uwav
  real(dp)                        :: one_over_jas2_av
  real(dp), allocatable           :: gradient_ovlp (:)
  real(dp)                        :: gradient_ovlp_norm 
  real(dp), allocatable           :: ovlp_fn (:,:)

  contains

!===========================================================================
  subroutine opt_ovlp_fn_menu
!---------------------------------------------------------------------------
! Description : menu for overlap with FN wavefn. optimization
!
! Created     : Cyrus Umrigar and Frank Petruzielo 7 Jun 2010
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'opt_ovlp_fn_menu'

! begin

  write(6,'(''Overlap-FN optimization'')')

! loop over menu lines
  do
  call get_next_word(word)

  select case (trim(word))
  case ('help')
   write(6,'(a)') ' HELP for overlap fixed-node optimization menu'
   write(6,'(a)') '  overlap_fn'
   write(6,'(a)') '   weight_trial = [real] : Psi_FN = weight_trial *  Psi_T + sum_{j=1}^{NCSF} delta_j Psi_j'
   write(6,'(a)') '                - weight_trial=1: zero-variance algorithm (default)'
   write(6,'(a)') '                - weight_trial=0: Reboredo method'
   write(6,'(a)') '   opt_ovlp_fn_nonsemiortho = [logical] : Solve linear system in non-semiorthogonal basis. (default=false)'
   write(6,'(a)') '   opt_ovlp_fn_linear = [logical] : Perform overlap-fn method with weight of 1/J^2. Assume overlap matrix is diagonal. (default=false)'
   write(6,'(a)') '   ovlp_identity = [logical] : Assume overlap matrix is the identiry matrix in opt_ovlp_fn_linear. (default=false)'
   write(6,'(a)') '   ovlp_diagonal = [logical] : Use diagonal overlap matrix only in opt_ovlp_fn_linear. (default=true)'
   write(6,'(a)') '   wt_lambda = [real] : Power to which DMC wts are raised unit physical time later. (default=1.d0)'
   write(6,'(a)') '   l_opt_ovlp_branching = [logical] : do branching (default=false)'
   write(6,'(a)') '   update_nonlinear = [original|semiorthogonal] : default=semiorthogonal, choice of update of nonlinear paramaters'
   write(6,'(a)') '   xi = [real] : update of nonlinear paramaters by orthogonalization to xi Psi_0 +(1-xi) Psi_lin'
   write(6,'(a)') '                - xi=1: orthogonalization to Psi_0 (default)'
   write(6,'(a)') '                - xi=0: orthogonalization to Psi_lin, ensures min |Psi_lin-Psi_0|'
   write(6,'(a)') '                - xi=0.5: orthogonalization to Psi_0 + Psi_lin, ensures |Psi_0|=|Psi_lin|'
   write(6,'(a)') ' end'

  case ('weight_trial')
   call get_next_value (weight_trial)

  case ('opt_ovlp_fn_nonsemiortho')
   call get_next_value (l_opt_ovlp_fn_nonsemiortho)

  case ('opt_ovlp_fn_linear')
   call get_next_value (l_opt_ovlp_fn_linear)

  case ('ovlp_identity')
   call get_next_value (l_ovlp_identity)

  case ('ovlp_diagonal')
   call get_next_value (l_ovlp_diagonal)

  case ('wt_lambda')
   call get_next_value( wt_lambda)
   write(6,'(''wt_lambda='',es10.3)') wt_lambda

  case ('l_opt_ovlp_branching')
   call get_next_value(l_opt_ovlp_branching)

  case ('update_nonlinear')
   call get_next_value (update_nonlinear)

  case ('xi')
   call get_next_value (xi)
   write(6,'(''xi='',f6.3)') xi

  case ('end')
   exit

  case default
   call die(lhere, 'unknown word >'+trim(word)+'<')
  end select

  enddo ! end loop over menu lines

  call require (lhere, 'wt_lambda >= 0 and wt_lambda =< 1', wt_lambda >= 0.d0 .and. wt_lambda <= 1.d0)

  if (l_opt_ovlp_fn_linear) then
   write(6,'(a)') ' overlap-fn-linear optimization method (with or without full overlap matrix)'
   if (l_ovlp_identity) then
    write(6,'(a)') ' assume identity matrix for the overlap matrix'
   endif
   if (l_ovlp_diagonal) then
    write(6,'(a)') ' assume diagonal matrix for the overlap matrix'
   endif
   if (l_ovlp_identity .and. l_ovlp_diagonal) then
    call die (lhere, 'ovlp_identity and ovlp_diagonal are both true.')
   endif
  endif

  end subroutine opt_ovlp_fn_menu

! ==============================================================================
  subroutine ovlp_trial_fn_bld
! ------------------------------------------------------------------------------
! Description   : overlap of trial wavefunction with fn wavefunction
! Description   : < Psi(0) Psi_FN > / sqrt(<Psi(0) Psi(0) > <Psi_FN Psi_FN>)
! Description   : <Psi_FN Psi_FN> is approximated as the sum of the squares of the weights
!
! Created       : Cyrus Umrigar and Frank Petruzielo  1 Jul 2010
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if(header_exe) then

   call object_create('ovlp_trial_fn')

   call object_needed ('walker_weights_sum')
   call object_needed ('walker_weights_sq_sum')
   call object_needed ('total_iterations_nb')

   return

  endif

! begin
  ovlp_trial_fn = walker_weights_sum / sqrt(total_iterations_nb * walker_weights_sq_sum)

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

 ! ==============================================================================
   subroutine dpsi_over_jas2_bld
 ! ------------------------------------------------------------------------------
 ! Description   : derivatives (Psi_i/Psi_0) (1/Jastrow^2)
 ! Description   : for overlap_fn method with linear scaling
 !
 ! Created       : J. Toulouse, 21 Jul 2010
 ! ------------------------------------------------------------------------------
   implicit none

 ! local
   integer param_i

 ! header
   if(header_exe) then

    call object_create ('dpsi_over_jas2')
    call object_average_define ('dpsi_over_jas2', 'dpsi_over_jas2_av')
    call object_average_unweighted_define ('dpsi_over_jas2', 'dpsi_over_jas2_uwav')

    call object_needed ('param_nb')
    call object_needed ('dpsi')
    call object_needed ('psi_jas')

    return

   endif

 ! begin

 ! allocations
   call object_alloc ('dpsi_over_jas2', dpsi_over_jas2, param_nb)
   call object_alloc ('dpsi_over_jas2_av', dpsi_over_jas2_av, param_nb)
   call object_alloc ('dpsi_over_jas2_uwav', dpsi_over_jas2_uwav, param_nb)

   do param_i = 1, param_nb
     dpsi_over_jas2 (param_i) = dpsi (param_i) / (psi_jas**2)
   enddo

   end subroutine dpsi_over_jas2_bld

 ! ==============================================================================
   subroutine dpsi2_over_jas2_bld
 ! ------------------------------------------------------------------------------
 ! Description   :  (Psi_i/Psi_0) (Psi_i/Psi_0) (1/Jastrow^2)
 ! Description   : for overlap_fn method with linear scaling
 !
 ! Created       : F. Petruzielo, 10 Aug 2010
 ! ------------------------------------------------------------------------------
   implicit none

 ! local
   integer param_i
   integer param_j !temp

 ! header
   if(header_exe) then

    call object_create ('dpsi2_over_jas2')
    call object_create ('dpsi_dpsi_over_jas2') !temp
    call object_average_unweighted_define ('dpsi2_over_jas2', 'dpsi2_over_jas2_uwav')
    call object_average_unweighted_define ('dpsi_dpsi_over_jas2', 'dpsi_dpsi_over_jas2_uwav') !temp

    call object_needed ('param_nb')
    call object_needed ('dpsi')
    call object_needed ('psi_jas')

    return

   endif

 ! begin

 ! allocations
   call object_alloc ('dpsi2_over_jas2', dpsi2_over_jas2, param_nb)
   call object_alloc ('dpsi2_over_jas2_uwav', dpsi2_over_jas2_uwav, param_nb)
   call object_alloc ('dpsi_dpsi_over_jas2', dpsi_dpsi_over_jas2, param_nb, param_nb) !temp
   call object_alloc ('dpsi_dpsi_over_jas2_uwav', dpsi_dpsi_over_jas2_uwav, param_nb, param_nb) !temp

   do param_i = 1, param_nb
     dpsi2_over_jas2 (param_i) = dpsi (param_i)**2 / (psi_jas**2)
     do param_j =1, param_nb !temp
        dpsi_dpsi_over_jas2 (param_i,param_j) = dpsi (param_i) * dpsi (param_j) / (psi_jas**2) !temp
     end do !temp
   enddo


   end subroutine dpsi2_over_jas2_bld


 ! ==============================================================================
   subroutine first_csf_over_jas2_bld
 ! ------------------------------------------------------------------------------
 ! Description   : for overlap_fn method with linear scaling
 !
 ! Created       : J. Toulouse, 21 Jul 2010
 ! ------------------------------------------------------------------------------
   implicit none

 ! header
   if(header_exe) then

    call object_create ('first_csf_over_jas2')
    call object_average_define ('first_csf_over_jas2', 'first_csf_over_jas2_av')
    call object_average_unweighted_define ('first_csf_over_jas2', 'first_csf_over_jas2_uwav')

    call object_needed ('first_csf_over_psid')
    call object_needed ('psi_jas')

    return

   endif

 ! begin
   call object_associate ('first_csf_over_jas2', first_csf_over_jas2)
   call object_associate ('first_csf_over_jas2_av', first_csf_over_jas2_av)
   call object_associate ('first_csf_over_jas2_uwav', first_csf_over_jas2_uwav)

   first_csf_over_jas2 = first_csf_over_psid / (psi_jas**2)

   end subroutine first_csf_over_jas2_bld

 ! ==============================================================================
   subroutine first_csf2_over_jas2_bld
 ! ------------------------------------------------------------------------------
 ! Description   : for overlap_fn method with linear scaling
 !
 ! Created       : J. Toulouse, 21 Jul 2010
 ! ------------------------------------------------------------------------------
   implicit none

 ! header
   if(header_exe) then

    call object_create ('first_csf2_over_jas2')
    call object_average_unweighted_define ('first_csf2_over_jas2', 'first_csf2_over_jas2_uwav')

    call object_needed ('first_csf_over_psid')
    call object_needed ('psi_jas')

    return

   endif

 ! begin
   call object_associate ('first_csf2_over_jas2', first_csf2_over_jas2)
   call object_associate ('first_csf2_over_jas2_uwav', first_csf2_over_jas2_uwav)

   first_csf2_over_jas2 = first_csf_over_psid **2 / (psi_jas**2)

   end subroutine first_csf2_over_jas2_bld

 ! ==============================================================================
   subroutine one_over_jas2_bld
 ! ------------------------------------------------------------------------------
 ! Description   : for ovlp_fn optimization method with linear scaling
 !
 ! Created       : F. Petruzielo, 10 Aug 2010
 ! ------------------------------------------------------------------------------
   implicit none

 ! header
   if(header_exe) then

    call object_create ('one_over_jas2')
    call object_average_unweighted_define ('one_over_jas2', 'one_over_jas2_uwav')
    call object_average_define ('one_over_jas2', 'one_over_jas2_av') !temp

    call object_needed ('psi_jas')

    return

   endif

 ! begin
   call object_associate ('one_over_jas2', one_over_jas2)
   call object_associate ('one_over_jas2_uwav', one_over_jas2_uwav)
   call object_associate ('one_over_jas2_av', one_over_jas2_av) !temp

   one_over_jas2 = 1.d0 / (psi_jas**2)

   end subroutine one_over_jas2_bld

!===========================================================================
  subroutine gradient_ovlp_bld
! ------------------------------------------------------------------------------
! Description   : gradient of overlap for overlap FN method
!
! Created       : J. Toulouse, 26 Jul 2012
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i

! header
  if (header_exe) then

   call object_create ('gradient_ovlp')

   call object_needed ('param_nb')
   call object_needed ('dpsi_av')
   call object_needed ('dpsi_uwav')

   return

  endif

! allocation
  call object_alloc ('gradient_ovlp', gradient_ovlp, param_nb)

  do param_i = 1, param_nb
   gradient_ovlp (param_i) = dpsi_av (param_i) - dpsi_uwav (param_i)
  enddo

  end subroutine gradient_ovlp_bld

!===========================================================================
  subroutine gradient_ovlp_norm_bld
! ------------------------------------------------------------------------------
! Description   : norm of gradient of overlap for overlap FN method
!
! Created       : J. Toulouse, 26 Jul 2012
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i

! header
  if (header_exe) then

   call object_create ('gradient_ovlp_norm')

   call object_needed ('param_nb')
   call object_needed ('gradient_ovlp')

   return

  endif

  gradient_ovlp_norm = 0.d0
  do param_i = 1, param_nb
   gradient_ovlp_norm =  gradient_ovlp_norm + gradient_ovlp (param_i)**2
  enddo
  gradient_ovlp_norm = dsqrt(gradient_ovlp_norm)

  end subroutine gradient_ovlp_norm_bld

 ! ==============================================================================
   subroutine delta_ovlp_fn_linear_bld
 ! ------------------------------------------------------------------------------
 ! Description   : variation of parameters for overlap FN method
 ! Description   : with diagonal overlap matrix and weighting of 1/J^2
 !
 ! Created       : F. Petruzielo, 10 Aug 2010
 ! ------------------------------------------------------------------------------
   use all_modules_mod
   implicit none

 ! local
   character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_linear_bld'
   integer param_i

 ! header
   if(header_exe) then

    call object_create ('delta_ovlp_fn_linear')

    call object_needed ('param_nb')
    call object_needed ('dpsi_over_jas2_av')
    call object_needed ('dpsi_over_jas2_uwav')
    call object_needed ('dpsi2_over_jas2_uwav')
    call object_needed ('first_csf_over_jas2_av')
    call object_needed ('first_csf_over_jas2_uwav')
    call object_needed ('first_csf2_over_jas2_uwav')
    call object_needed ('one_over_jas2_uwav')
    call object_needed ('one_over_jas2_av') !temp
    call object_needed ('ovlp_trial_fn')
    call object_needed ('ovlp_trial_fn_over_ovlp_trial')
    return

   endif

 ! begin

 ! allocation
   call object_alloc ('delta_ovlp_fn_linear', delta_ovlp_fn_linear, param_nb)

!  rescaling so that first CSF coefficient is unchanged
   call object_provide ('csf_coef')
   call object_provide ('nparmcsf')
   call object_provide ('iwcsf')
   write(6,*) "S=1", one_over_jas2_uwav / sum(csf_coef(:,1)**2)
   write(6,*) "S=S_{ii}=", dpsi2_over_jas2_uwav(:)
   do param_i = 1, param_nb !temp
      write(6,*) "S=S_{ij}=", dpsi_dpsi_over_jas2_uwav (param_i,:)  !temp
   end do !temp

   if (l_ovlp_identity) then
      delta_ovlp_fn_linear (:) = (dpsi_over_jas2_av (:) * ovlp_trial_fn  - weight_trial * dpsi_over_jas2_uwav (:)) * sum(csf_coef(:,1)**2) / one_over_jas2_uwav
      do param_i = 1, nparmcsf
         delta_ovlp_fn_linear (param_i) = (weight_trial * csf_coef (iwcsf(param_i), 1) + delta_ovlp_fn_linear (param_i)) * csf_coef(1,1) / (weight_trial * csf_coef(1,1) + (first_csf_over_jas2_av * ovlp_trial_fn - weight_trial * first_csf_over_jas2_uwav) * sum(csf_coef(:,1)**2) / one_over_jas2_uwav) - csf_coef (iwcsf(param_i), 1)
      enddo

   elseif (l_ovlp_diagonal) then
!JT      delta_ovlp_fn_linear (:) = (dpsi_over_jas2_av (:) * ovlp_trial_fn - weight_trial * dpsi_over_jas2_uwav (:)) / dpsi2_over_jas2_uwav(:)
      delta_ovlp_fn_linear (:) = (dpsi_over_jas2_av (:) * ovlp_trial_fn_over_ovlp_trial - weight_trial * dpsi_over_jas2_uwav (:)) / dpsi2_over_jas2_uwav(:)
     
!     rescaling of CSF parameters 
      if (weight_trial /= 1.d0) then
       call die (lhere, 'case weight_trial /= 1 not implemented')
      endif
      
      call object_provide ('is_param_type_csf')
      call object_provide ('is_param_type_jas')
      call object_provide ('is_param_type_pjas')
      call object_provide ('is_param_type_exp')
      call object_provide ('is_param_type_orb')
      call object_provide ('is_param_type_geo')
      do param_i = 1, param_nb
       if (is_param_type_csf(param_i)) then
        delta_ovlp_fn_linear (param_i) = (csf_coef(1,1)/(csf_coef(1,1)+delta_ovlp_fn_linear(1)))*(csf_coef (iwcsf(param_i), 1) + delta_ovlp_fn_linear(param_i)) - csf_coef (iwcsf(param_i), 1)
       elseif (is_param_type_jas(param_i)) then
        call die (lhere, 'not implemented for Jastrow parameters')
       elseif (is_param_type_pjas(param_i)) then
        call die (lhere, 'not implemented for Jastrow parameters')
       elseif (is_param_type_exp(param_i)) then
        call die (lhere, 'not implemented for exponent parameters')
       elseif (is_param_type_orb(param_i)) then
        delta_ovlp_fn_linear (param_i) = (csf_coef(1,1)/(csf_coef(1,1)+delta_ovlp_fn_linear(1)))*(delta_ovlp_fn_linear(param_i))
       elseif (is_param_type_geo(param_i)) then
        call die (lhere, 'not implemented for geometry parameters')
       else
        call die (lhere, 'unknown parameter type')
       endif
      enddo
      

!JT      do param_i = 1, nparmcsf
!JT         delta_ovlp_fn_linear (param_i) = (weight_trial * csf_coef (iwcsf(param_i), 1) + delta_ovlp_fn_linear (param_i)) * csf_coef(1,1) / (weight_trial * csf_coef(1,1) + (first_csf_over_jas2_av * ovlp_trial_fn - weight_trial * first_csf_over_jas2_uwav) / first_csf2_over_jas2_uwav) - csf_coef (iwcsf(param_i), 1)
!JT      end do
      ! do param_i = nparmcsf+1, param_nb !temp
      !    delta_ovlp_fn_linear (param_i) = (delta_ovlp_fn_linear (param_i)) * csf_coef(1,1) / (weight_trial * csf_coef(1,1) + (first_csf_over_jas2_av * ovlp_trial_fn - weight_trial * first_csf_over_jas2_uwav) / first_csf2_over_jas2_uwav) !temp
      ! enddo !temp

    else
!   with full overlap matrix:
      call die (lhere, 'full overlap matrix not implemented')
!      dpsi_dpsi_over_jas2_uwav
     
   end if

   ! if (l_ovlp_identity) then   !temp
   !    delta_ovlp_fn_linear (:) = (dpsi_over_jas2_av (:) * one_over_jas2_uwav / one_over_jas2_av  - weight_trial * dpsi_over_jas2_uwav (:)) * sum(csf_coef(:,1)**2) / one_over_jas2_uwav   !temp
   !    do param_i = 1, nparmcsf   !temp
   !       delta_ovlp_fn_linear (param_i) = (weight_trial * csf_coef (iwcsf(param_i), 1) + delta_ovlp_fn_linear (param_i)) * csf_coef(1,1) / (weight_trial * csf_coef(1,1) + (first_csf_over_jas2_av * one_over_jas2_uwav / one_over_jas2_av - weight_trial * first_csf_over_jas2_uwav) * sum(csf_coef(:,1)**2) / one_over_jas2_uwav) - csf_coef (iwcsf(param_i), 1)   !temp
   !    enddo   !temp
   ! else   !temp
   !    delta_ovlp_fn_linear (:) = (dpsi_over_jas2_av (:) * one_over_jas2_uwav / one_over_jas2_av - weight_trial * dpsi_over_jas2_uwav (:)) / dpsi2_over_jas2_uwav(:)   !temp
   !    do param_i = 1, nparmcsf   !temp
   !       delta_ovlp_fn_linear (param_i) = (weight_trial * csf_coef (iwcsf(param_i), 1) + delta_ovlp_fn_linear (param_i)) * csf_coef(1,1) / (weight_trial * csf_coef(1,1) + (first_csf_over_jas2_av * one_over_jas2_uwav / one_over_jas2_av - weight_trial * first_csf_over_jas2_uwav) / first_csf2_over_jas2_uwav) - csf_coef (iwcsf(param_i), 1)   !temp
   !    end do   !temp
   !    do param_i = nparmcsf+1, param_nb !temp   !temp
   !       delta_ovlp_fn_linear (param_i) = (delta_ovlp_fn_linear (param_i)) * csf_coef(1,1) / (weight_trial * csf_coef(1,1) + (first_csf_over_jas2_av * one_over_jas2_uwav / one_over_jas2_av - weight_trial * first_csf_over_jas2_uwav) / first_csf2_over_jas2_uwav) !temp   !temp
   !    enddo !temp   !temp
   ! end if   !temp


   end subroutine delta_ovlp_fn_linear_bld

!===========================================================================
  subroutine delta_ovlp_fn_bld
! ------------------------------------------------------------------------------
! Description   : variation of parameters for overlap FN method
! Description   : with overlap matrix in semiorthogonalized basis
!
! Created       : J. Toulouse, 7 Jun 2010
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_bld'
  integer param_i, param_j, info
  real(dp) transform_factor
  real(dp), allocatable :: dpsi_dpsi_uwcovar_stab (:,:)
  real(dp), allocatable :: dpsi_dpsi_uwcovar_stab_inv (:,:)
  real(dp), allocatable :: vect (:)
  real(dp) threshold

! header
  if (header_exe) then

   call object_create ('delta_ovlp_fn')

   call object_needed ('param_nb')
   call object_needed ('dpsi_av')
   call object_needed ('dpsi_uwav')
   call object_needed ('dpsi_dpsi_uwcovar')
   call object_needed ('gradient_ovlp')
   call object_needed ('diag_stab')

   return

  endif

! allocation
  call object_alloc ('delta_ovlp_fn', delta_ovlp_fn, param_nb)

  write(6,'(a)') "Projection of wave function derivatives unto fixed-node and trial wave functions:"
  do param_i = 1, param_nb
   write(6,'(a,i5,a,f12.6,a,f12.6,a,f12.6)') "parameter # ",param_i, ": <Psi_i|Psi_FN>= ", dpsi_av (param_i), ", <Psi_i|Psi_0>= ", dpsi_uwav (param_i), ", difference= ", dpsi_av (param_i) - dpsi_uwav (param_i)
  enddo

! stabilize overlap matrix
  call alloc ('dpsi_dpsi_uwcovar_stab', dpsi_dpsi_uwcovar_stab, param_nb, param_nb)
  dpsi_dpsi_uwcovar_stab (:,:) = dpsi_dpsi_uwcovar (:,:)
  do param_i = 1, param_nb
    dpsi_dpsi_uwcovar_stab (param_i, param_i) = dpsi_dpsi_uwcovar_stab (param_i, param_i) + diag_stab
  enddo

! find parameter variations by inverting overlap matrix
!  call alloc ('dpsi_dpsi_uwcovar_stab_inv', dpsi_dpsi_uwcovar_stab_inv, param_nb, param_nb)
!  threshold = 1.d-10
!  threshold = 0.d0 !!!!!
!  call inverse_by_svd (dpsi_dpsi_uwcovar_stab, dpsi_dpsi_uwcovar_stab_inv, param_nb, threshold)
!  do param_i = 1, param_nb
!   delta_ovlp_fn (param_i) = 0.d0
!   do param_j = 1, param_nb
!    delta_ovlp_fn (param_i) = delta_ovlp_fn (param_i) + dpsi_dpsi_uwcovar_stab_inv (param_i, param_j) * gradient_ovlp (param_j)
!   enddo
!  enddo

! find parameter variations by solving the linear equation problem
  call alloc ('vect', vect, param_nb)
  vect(:) = gradient_ovlp (:)
  call dposv('U',param_nb,1,dpsi_dpsi_uwcovar_stab,param_nb,vect,param_nb,info)
!  do param_i = 1, param_nb
!   write(6,'(a,i4,a,f12.6,a,f12.6)') 'param_i=',param_i,' delta_ovlp_fn=',delta_ovlp_fn (param_i), ' vect=',vect (param_i)
!  enddo
  delta_ovlp_fn (:) = vect (:)
  

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

!===========================================================================
  subroutine ovlp_fn_bld
! ------------------------------------------------------------------------------
! Description   : overlap matrix for overlap-fn optimization method in
! Description   : non-semiorthogonalized basis
!
! Created       : J. Toulouse, 30 Jul 2012
! ------------------------------------------------------------------------------
  implicit none

! local
  integer param_i, param_j
  real(dp), allocatable :: ovlp_fn_eigvec (:,:)
  real(dp), allocatable :: ovlp_fn_eigval (:)


! header
  if (header_exe) then

   call object_create ('ovlp_fn')

   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('dpsi_uwav')
   call object_needed ('dpsi_dpsi_uwav')
   call object_needed ('diag_stab')
   return

  endif

! allocation
  call object_alloc ('ovlp_fn', ovlp_fn, param_nb+1,param_nb+1)
  ovlp_fn (1,1) = 1.d0
  do param_i = 1, param_nb
   ovlp_fn (1,param_i+1) = dpsi_uwav (param_i)
   ovlp_fn (param_i+1,1) = dpsi_uwav (param_i)
   do param_j = 1, param_nb
     ovlp_fn (param_i+1,param_j+1) = dpsi_dpsi_uwav (param_pairs(param_i,param_j))
   enddo
  enddo

! check eigenvalues
  call alloc('ovlp_fn_eigvec', ovlp_fn_eigvec, param_nb+1, param_nb+1)
  call alloc('ovlp_fn_eigval', ovlp_fn_eigval, param_nb+1)
  call eigensystem(ovlp_fn, ovlp_fn_eigvec, ovlp_fn_eigval, param_nb+1)

  write(6,*)
  write(6,'(a)') 'Eigenvalues of overlap matrix of current wave function and its first-order derivatives:'
  do param_i = 1, param_nb+1
    write(6,'(a,i4,a,es15.8)') 'overlap eigenvalue # ',param_i,': ',ovlp_fn_eigval(param_i)
  enddo

  end subroutine ovlp_fn_bld

!===========================================================================
  subroutine delta_ovlp_fn_nonsemiortho_bld
! ------------------------------------------------------------------------------
! Description   : variation of parameters for overlap FN method
! Description   : with overlap matrix in non-semiorthogonalized basis
!
! Created       : J. Toulouse, 30 Jul 2012
! ------------------------------------------------------------------------------
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'delta_ovlp_fn_nonsemiortho_bld'
  integer param_i, param_j
  real(dp) transform_factor
  real(dp), allocatable :: ovlp_fn_stab (:,:)
  real(dp), allocatable :: ovlp_fn_stab_inv (:,:)
  real(dp), allocatable :: vector (:)
  real(dp), allocatable :: delta_ovlp_fn_plus_one (:)
  real(dp) threshold

! header
  if (header_exe) then

   call object_create ('delta_ovlp_fn_nonsemiortho')

   call object_needed ('param_nb')
   call object_needed ('ovlp_fn')
   call object_needed ('dpsi_av')
   call object_needed ('ovlp_trial_fn_over_ovlp_trial')
   call object_needed ('diag_stab')

   return

  endif

! allocation
  call object_alloc ('delta_ovlp_fn_nonsemiortho', delta_ovlp_fn_nonsemiortho, param_nb)

! stabilize overlap matrix
  call alloc ('ovlp_fn_stab', ovlp_fn_stab, param_nb+1, param_nb+1)
  ovlp_fn_stab (:,:) = ovlp_fn (:,:)
  do param_i = 2, param_nb
    ovlp_fn_stab (param_i, param_i) = ovlp_fn (param_i, param_i) + diag_stab
  enddo

! inverse overlap matrix
  call alloc ('ovlp_fn_stab_inv', ovlp_fn_stab_inv, param_nb+1, param_nb+1)
  threshold = 1.d-10
  call inverse_by_svd (ovlp_fn_stab, ovlp_fn_stab_inv, param_nb+1, threshold)

! right-hand-side vector
  call alloc ('vector', vector, param_nb+1)
  vector (1) = ovlp_trial_fn_over_ovlp_trial
  do param_i = 1,  param_nb
    vector (param_i+1) = dpsi_av (param_i) * ovlp_trial_fn_over_ovlp_trial
  enddo

! find parameter variations
  call alloc ('delta_ovlp_fn_plus_one', delta_ovlp_fn_plus_one, param_nb+1)
  do param_i = 1, param_nb+1
   delta_ovlp_fn_plus_one (param_i) = 0.d0
   do param_j = 1, param_nb+1
    delta_ovlp_fn_plus_one (param_i) = delta_ovlp_fn_plus_one (param_i) + ovlp_fn_stab_inv (param_i, param_j) * vector (param_j)
   enddo
  enddo
  do param_i = 1, param_nb
   delta_ovlp_fn_nonsemiortho (param_i) = delta_ovlp_fn_plus_one (param_i+1) / delta_ovlp_fn_plus_one (1)
  enddo

! norm of linear wave function variation for nonlinear parameter
  call object_provide ('nparmcsf')
  psi_lin_var_norm = 0.d0
  do param_i = nparmcsf+1, param_nb
   do param_j = nparmcsf+1, param_nb
     psi_lin_var_norm = psi_lin_var_norm + delta_ovlp_fn_nonsemiortho (param_i) * delta_ovlp_fn_nonsemiortho (param_j) * ovlp_fn (param_i+1, param_j+1)
   enddo
  enddo
  write(6,'(a,f10.6)') 'Norm of linear wave function variation for nonlinear parameters =', psi_lin_var_norm

! transformation of parameter variations to desired basis
  transform_factor = 1.d0
  select case (trim(update_nonlinear))

! original: 
   case ('original')

! semiorthogonal: use semiorthognal derivatives for nonlinear parameters
   case ('semiorthogonal')
    call die (lhere, 'semiorthogonal not implemented')

  case default
    call die (lhere, 'unknown update choice >'+trim(update_nonlinear)+'<')
  end select

! final parameter variations
!  delta_ovlp_fn (:) = delta_ovlp_fn (:) / transform_factor

  end subroutine delta_ovlp_fn_nonsemiortho_bld

end module opt_ovlp_fn_mod
