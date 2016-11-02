module linearresponse_mod

  use all_tools_mod
  use opt_common_mod
  use dets_mod
  use vmc_mod
  use deriv_mod
  use opt_lin_mod

! "tda_only"          : will not calculate the full linresp equations, but only the A matrix
! "hessian"           : will calculate the Hessian of the SCF (A B \\ B A)
! "real_hessian"      : will calculate the real Hessian of the SCF (A+B)
! "compare_to_tda"    : will calculate the linresp equations (A 0 \\ 0 A)
! "compare_to_optlin" : will give an output comparable to that obtained with
!                       an "optimization" command together with "compare_to_linresp"

! Declaration of global variables and default values

! TDA
  logical                         :: l_tda_only         =.false.
  real(dp), allocatable           :: tda_av_eigenval_r(:)
  real(dp), allocatable           :: tda_av_eigenval_i(:)
  real(dp), allocatable           :: tda_av_eigenval_r_err(:)
  real(dp), allocatable           :: tda_av_eigenval_i_err(:)
  real(dp), allocatable           :: tda_av_eigenval_via_super_r(:)
  real(dp), allocatable           :: tda_av_eigenval_via_super_i(:)
  real(dp), allocatable           :: tda_av_eigenval_via_super_r_err(:)
  real(dp), allocatable           :: tda_av_eigenval_via_super_i_err(:)
! prints and debug
  logical                         :: l_print            =.false.
  logical                         :: l_print_eigenvec   =.false.
  logical                         :: l_print_every_block=.false.
  logical                         :: l_compare_to_tda   =.false.
  logical                         :: l_symm_amat        =.false.
! hessian
  logical                         :: l_hessian          =.false.
  real(dp), allocatable           :: hessian_av_eigenval_r(:)
  real(dp), allocatable           :: hessian_av_eigenval_i(:)
  real(dp), allocatable           :: hessian_av_eigenval_r_err(:)
  real(dp), allocatable           :: hessian_av_eigenval_i_err(:)
! real hessian
  logical                         :: l_real_hessian     =.false.
  real(dp), allocatable           :: real_hessian_av_eigenval_r(:)
  real(dp), allocatable           :: real_hessian_av_eigenval_i(:)
  real(dp), allocatable           :: real_hessian_av_eigenval_r_err(:)
  real(dp), allocatable           :: real_hessian_av_eigenval_i_err(:)
! "normal" linresp
  real(dp), allocatable           :: amat_av(:,:)
  real(dp), allocatable           :: super_amat_av(:,:)
  real(dp), allocatable           :: bmat_av(:,:)
  real(dp), allocatable           :: ovlp_psii_psij_av(:,:)
  real(dp), allocatable           :: super_ovlp_psii_psij_av(:,:)
  real(dp), allocatable           :: linresp_mat(:,:)
  real(dp), allocatable           :: ovlp_mat(:,:)
  real(dp), allocatable           :: linresp_av_eigenval_r(:)
  real(dp), allocatable           :: linresp_av_eigenval_i(:)
  real(dp), allocatable           :: linresp_av_eigenval_r_err(:)
  real(dp), allocatable           :: linresp_av_eigenval_i_err(:)

  contains

!===========================================================================
  subroutine linearresponse_menu
!---------------------------------------------------------------------------
! Description : menu for linear-response calculations
!
! Created     : J. Toulouse, 18 May 2016
! Modified    : B. Mussard, June 2016
!---------------------------------------------------------------------------
  use all_modules_mod
  use derivatives_fast_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'linearresponse_menu'
  integer                 :: i
  integer                 :: parameter_type_nb = 0
# if defined (PATHSCALE)
   character(len=max_string_len) :: parameter_type (max_string_array_len) ! for pathscale compiler
# else
   character(len=max_string_len), allocatable :: parameter_type (:)
# endif

! begin
  write(6,*)
  write(6,'(a)') 'Beginning of linearresponse menu -------------------------------------------------------------------------'

! initialization

  itest=0
! loop over menu lines
  do
    call get_next_word (word)

    select case(trim(word))
    case ('help')
      write(6,'(a)') 'HELP for linearresponse menu:'
      write(6,'(a)') ' parameters jastrow csfs orbitals exponents end : list of parameter types to take into account'
      write(6,'(a)') 'linearresponse'
      write(6,'(a)') 'end'

    case ('parameters')
# if defined (PATHSCALE)
      call get_next_value_list_string ('parameter_type', parameter_type, parameter_type_nb) ! for pathscale compiler
# else
      call get_next_value_list ('parameter_type', parameter_type, parameter_type_nb)
# endif

!   triplet
    case ('triplet')
      call get_next_value (l_triplet)

!   TDA
    case ('tda_only')
      call get_next_value (l_tda_only)

!   Hessians
    case ('hessian')
      call get_next_value (l_hessian)
    case ('real_hessian')
      call get_next_value (l_real_hessian)

!   prints and debug
    case ('print_eigenvec')
      call get_next_value (l_print_eigenvec)
    case ('print_every_block')
      call get_next_value (l_print_every_block)
    case ('print')
      call get_next_value (l_print)
    case ('symm_amat')
      call get_next_value (l_symm_amat)
    case ('compare_to_tda')
      call get_next_value (l_compare_to_tda)
    case ('itest')
      call get_next_value (itest)
    case ('compare_to_optlin')
      call get_next_value (l_compare_linresp_and_optlin)
      if (l_compare_linresp_and_optlin) then
        l_tda_only=.true.
        l_print=.true.
      endif

    case ('end')
      exit

    case default
      call die (lhere, 'unknown keyword >'+trim(word)+'<.')
    end select

  enddo ! end loop over menu lines

! initialize
  l_opt=.true.
  l_opt_jas=.false.
  l_opt_csf=.false.
  l_opt_orb=.false.
  l_opt_exp=.false.
  igradhess=0

! subroutine gathering "old" commands needed to properly get "nparmj"
  call get_nparmj

! parameters type
  write(6,'(a,10a10)') ' Requested parameter types: ',parameter_type(:)
  do i = 1, parameter_type_nb
   select case(trim(parameter_type(i)))
   case ('jastrow')
    l_opt_jas = .true.
   case ('csfs')
    l_opt_csf = .true.
   case ('orbitals')
    l_opt_orb = .true.
   case ('exponents')
    l_opt_exp = .true.
   case default
    call die (lhere, 'unknown parameter type >'+trim(parameter_type(i))+'<.')
   end select
  enddo

! ORB parameters
  if (l_opt_orb) then
    write(6,*)
    write(6,'(3a)') ' Orbital parameter information:'
    call object_provide ('param_orb_nb')
    call object_provide ('det_ex_unq_up_nb')
    call object_provide ('orb_opt_last_lab')
    norb = orb_opt_last_lab
    write(6,'(a,i8)') ' Number of computed orbitals will be ', norb
    if ((.not.l_compare_linresp_and_optlin).and.(.not.l_compare_to_tda).and.(.not.l_tda_only)) then
      call object_provide ('double_ex_nb')
    endif
  else
    param_orb_nb  =  0
    call object_modified ('param_orb_nb')
  endif

! EXP parameters
  if (l_opt_exp) then
    write(6,*)
    write(6,'(3a)') ' Exponent parameter information:'
    call object_provide ('param_exp_nb')
  else
    param_exp_nb  =  0
    call object_modified ('param_exp_nb')
  endif

! JAS parameters
  if (l_opt_jas) then
    l_opt_jas_2nd_deriv=.true.
  else
    nparmj=0
    call object_modified ('nparmj')
  endif

! CSF parameters
  if (l_opt_csf) then
    call object_provide ('ncsf')
    call object_provide ('nparmcsf')
  else
    nparmcsf=0
    call object_modified ('nparmcsf')
  endif
  if (ibasis.le.3) nparmd=nparmcsf
  call object_modified ('nparmd')

! final messages 
  nparm = nparmj+nparmcsf
  call object_modified ('nparm')
  call object_provide ('nparmcsf')
  call object_provide ('param_orb_nb')
  call object_provide ('param_exp_nb')
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

  write(6,'(a)') 'End of linearresponse menu -------------------------------------------------------------------------------'

  call linearresponse

  end subroutine linearresponse_menu

!===========================================================================
  subroutine linearresponse
!---------------------------------------------------------------------------
! Description : routine for linear-response calculations
!
! Created     : J. Toulouse, 18 May 2016
! Modified    : B. Mussard, June 2016
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'linearresponse'

! begin
  write(6,*)
  write(6,'(a)') '************************************** LINEAR-RESPONSE CALCULATION ***************************************'

! initializations
  if (.not. l_mode_vmc) then 
    call die (lhere, 'mode must be vmc for linear-response calculations')
  endif 
  call vmc_init

! average needed by the linresp calculations
  call object_average_request('dpsi_av')
  call object_average_request('dpsi_eloc_av')
  call object_average_request('dpsi_dpsi_av')
  call object_average_request('dpsi_dpsi_eloc_av')
  call object_average_request('deloc_av')
  call object_average_request('dpsi_deloc_av')
  if ((.not.l_compare_linresp_and_optlin).and.(.not.l_compare_to_tda).and.(.not.l_tda_only)) then
    call object_average_request('d2psi_av')
    call object_average_request('d2psi_eloc_av')
  endif

! error needed by the linresp calculations
  if (l_hessian) then
    call object_error_request('hessian_av_eigenval_r_err')
    call object_error_request('hessian_av_eigenval_i_err')
  endif
  if (l_real_hessian) then
    call object_error_request('real_hessian_av_eigenval_r_err')
    call object_error_request('real_hessian_av_eigenval_i_err')
  endif
  call object_error_request('tda_av_eigenval_via_super_r_err')
  call object_error_request('tda_av_eigenval_via_super_i_err')
  call object_error_request('tda_av_eigenval_r_err')
  call object_error_request('tda_av_eigenval_i_err')
  if (.not.l_tda_only) then
    call object_error_request('linresp_av_eigenval_r_err')
    call object_error_request('linresp_av_eigenval_i_err')
  endif

! VMC run
  nforce=1
  nwftype=1
  call vmc
  run_done=.true.

! test-phase of FastDet
!   call object_provide('gamma_up')

! final calculation
  if (l_hessian) then
    call hessian_av_eigenval_bld
  endif
  if (l_real_hessian) then
    call real_hessian_av_eigenval_bld
  endif
  call tda_av_eigenval_bld
  call tda_av_eigenval_via_super_bld
  if (.not.l_tda_only) then
    call linresp_av_eigenval_bld
  endif

  end subroutine linearresponse

!===========================================================================
  subroutine linresp_av_eigenval_bld
!---------------------------------------------------------------------------
! Description : Solve the generalized eigenvalue equation
! Description :   ABBA . evec =  eval S . evec
! Description : The eigenvalues are sorted and printed.
!
! Created     : B. Mussard, June 2016 (from "ham_lin_av_bld")
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none
  
! local
  character(len=max_string_len_rout), save :: lhere = 'linresp_av_eigenval_bld'
  integer                         :: i,info,lwork
  real(dp), allocatable           :: linresp_local(:,:)
  real(dp), allocatable           :: ovlp_local(:,:)
  real(dp), allocatable           :: eigvec(:,:)
  real(dp), allocatable           :: eigval_r(:)
  real(dp), allocatable           :: eigval_i(:)
  real(dp), allocatable           :: eigval_denom(:)
  real(dp), allocatable           :: work(:)
  integer,  allocatable           :: backward_sort(:)
  integer,  allocatable           :: forward_sort(:)

! begin
  if (header_exe) then
    call object_create('linresp_av_eigenval_r')
    call object_create('linresp_av_eigenval_i')
    call object_error_define ('linresp_av_eigenval_r','linresp_av_eigenval_r_err')
    call object_error_define ('linresp_av_eigenval_i','linresp_av_eigenval_i_err')

    call object_needed('param_nb')

    call object_needed('linresp_mat')
    call object_needed('ovlp_mat')

    return
  endif

  call object_alloc ('linresp_av_eigenval_r',linresp_av_eigenval_r,2*param_nb)
  call object_alloc ('linresp_av_eigenval_i',linresp_av_eigenval_i,2*param_nb)
  call object_alloc ('linresp_av_eigenval_r_err',linresp_av_eigenval_r_err,2*param_nb)
  call object_alloc ('linresp_av_eigenval_i_err',linresp_av_eigenval_i_err,2*param_nb)

! prepare arrays and output message
  call alloc('eigvec',       eigvec,       2*param_nb, 2*param_nb)
  call alloc('eigval_r',     eigval_r,     2*param_nb)
  call alloc('eigval_i',     eigval_i,     2*param_nb)
  call alloc('eigval_denom', eigval_denom, 2*param_nb)
  if (run_done.or.l_print_every_block) then
  write(6,*)
  write(6,'(a,1pd9.1)') 'Solving the LINRESP generalized eigenvalue equation'
  endif

! calculate optimal value of lwork
  lwork = 1
  call alloc('work', work, lwork)
  call alloc('linresp_local',linresp_local, 2*param_nb, 2*param_nb)
  call alloc('ovlp_local',   ovlp_local   , 2*param_nb, 2*param_nb)
  linresp_local=linresp_mat
  ovlp_local=ovlp_mat
  call dggev('N','V', &
             2*param_nb, linresp_local,  &
             2*param_nb, ovlp_local,     &
             2*param_nb, eigval_r, eigval_i, eigval_denom, &
             eigvec, 2*param_nb, eigvec, 2*param_nb, &
             work, -1, info)
  if(info /= 0) call die(lhere, 'problem in dggev_alloc: info='+info+' /= 0')
  lwork =  nint(work(1))
  call alloc('work', work, lwork)

! generalized eigenvalue problem
  call dggev('N','V', &
             2*param_nb, linresp_local,  &
             2*param_nb, ovlp_local,     &
             2*param_nb, eigval_r, eigval_i, eigval_denom, &
             eigvec, 2*param_nb, eigvec, 2*param_nb, &
             work, lwork, info)
  if(info /= 0) call die(lhere, 'problem in dggev: info='+info+' /= 0 (compare to 2*param_nb='+2*param_nb+')')

! scale eigenvalue
  do i = 1, 2*param_nb
    eigval_r(i) = eigval_r(i) / eigval_denom(i)
    eigval_i(i) = eigval_i(i) / eigval_denom(i)
  enddo

! sort eigenvalues
  call sorting_real_im(eigval_r,eigval_i,forward_sort,backward_sort,2*param_nb) 

! the sorted eigenvalues are "linresp_av_eigenval"
  do i = 1, 2*param_nb
   linresp_av_eigenval_r (i) = eigval_r(backward_sort(i))
   linresp_av_eigenval_i (i) = eigval_i(backward_sort(i))
  enddo

! output analysis
  call analysis_eigval_eigvec(eigvec,backward_sort, &
    linresp_av_eigenval_r,linresp_av_eigenval_r_err,&
    linresp_av_eigenval_i,linresp_av_eigenval_i_err,2*param_nb)

! release statements
  call release ('linresp_local', linresp_local)
  call release ('ovlp_local', ovlp_local)
  call release ('work', work)
  call release ('eigvec', eigvec)
  call release ('eigval_r', eigval_r)
  call release ('eigval_i', eigval_i)
  call release ('eigval_denom', eigval_denom)
  call release ('backward_sort', backward_sort)
  call release ('forward_sort', forward_sort)

  end subroutine linresp_av_eigenval_bld

!===========================================================================
  subroutine tda_av_eigenval_bld
!---------------------------------------------------------------------------
! Description : Solve the TDA eigenvalue equation
! Description :   A . evec =  eval evec
! Description : The eigenvalues are sorted and printed.
!
! Created     : B. Mussard, June 2016 (from "ham_lin_av_bld")
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none
  
! local
  character(len=max_string_len_rout), save :: lhere = 'tda_av_eigenval_bld'
  integer                         :: i,info,lwork
  real(dp), allocatable           :: amat_local(:,:)
  real(dp), allocatable           :: ovlp_local(:,:)
  real(dp), allocatable           :: eigvec(:,:)
  real(dp), allocatable           :: eigval_r(:)
  real(dp), allocatable           :: eigval_i(:)
  real(dp), allocatable           :: eigval_denom(:)
  real(dp), allocatable           :: work(:)
  integer,  allocatable           :: backward_sort(:)
  integer,  allocatable           :: forward_sort(:)

! begin
  if (header_exe) then
    call object_create('tda_av_eigenval_r')
    call object_create('tda_av_eigenval_i')
    call object_error_define ('tda_av_eigenval_r','tda_av_eigenval_r_err')
    call object_error_define ('tda_av_eigenval_i','tda_av_eigenval_i_err')

    call object_needed('param_nb')

    call object_needed('amat_av')
    call object_needed('ovlp_psii_psij_av')

    return
  endif

  call object_alloc ('tda_av_eigenval_r',tda_av_eigenval_r,param_nb)
  call object_alloc ('tda_av_eigenval_i',tda_av_eigenval_i,param_nb)
  call object_alloc ('tda_av_eigenval_r_err',tda_av_eigenval_r_err,param_nb)
  call object_alloc ('tda_av_eigenval_i_err',tda_av_eigenval_i_err,param_nb)

! prepare arrays and output message
  call alloc('eigvec',       eigvec,       param_nb, param_nb)
  call alloc('eigval_r',     eigval_r,     param_nb)
  call alloc('eigval_i',     eigval_i,     param_nb)
  call alloc('eigval_denom', eigval_denom, param_nb)
  if (run_done.or.l_print_every_block) then
  write(6,*)
  write(6,'(a,1pd9.1)') 'Solving the TDA eigenvalue equation'
  endif

! calculate optimal value of lwork
  lwork = 1
  call alloc('work', work, lwork)
  call alloc('amat_local', amat_local, param_nb, param_nb)
  call alloc('ovlp_local', ovlp_local, param_nb, param_nb)
  amat_local=amat_av
  ovlp_local=ovlp_psii_psij_av
  call dggev('N','V', &
             param_nb, amat_local,         &
             param_nb, ovlp_local,         &
             param_nb, eigval_r, eigval_i, eigval_denom, &
             eigvec, param_nb, eigvec, param_nb, &
             work, -1, info)
  if(info /=  0) call die(lhere, 'problem in dggev_alloc: info='+info+' /=  0')
  lwork=nint(work(1))
  call alloc('work', work, lwork)

! generalized eigenvalue problem
  call dggev('N','V', &
             param_nb, amat_local,         &
             param_nb, ovlp_local,         &
             param_nb, eigval_r, eigval_i, eigval_denom, &
             eigvec, param_nb, eigvec, param_nb, &
             work, lwork, info)
  if(info /=  0) call die(lhere, 'problem in dggev: info='+info+' /=  0 (compare to  param_nb='+param_nb+')')

! scale eigenvalue
  do i = 1, param_nb
    eigval_r(i) = eigval_r(i) / eigval_denom(i)
    eigval_i(i) = eigval_i(i) / eigval_denom(i)
  enddo

! sort eigenvalues
  call sorting_real_im(eigval_r,eigval_i,forward_sort,backward_sort,param_nb) 

! the sorted eigenvalues are "linresp_av_eigenval"
  do i = 1, param_nb
   tda_av_eigenval_r (i) = eigval_r(backward_sort(i))
   tda_av_eigenval_i (i) = eigval_i(backward_sort(i))
  enddo

! output analysis
  call analysis_eigval_eigvec(eigvec,backward_sort, &
    tda_av_eigenval_r,tda_av_eigenval_r_err, &
    tda_av_eigenval_i,tda_av_eigenval_i_err,param_nb)

! release statements
  call release ('amat_local', amat_local)
  call release ('ovlp_local', ovlp_local)
  call release ('work', work)
  call release ('eigvec', eigvec)
  call release ('eigval_r', eigval_r)
  call release ('eigval_i', eigval_i)
  call release ('eigval_denom', eigval_denom)
  call release ('backward_sort', backward_sort)
  call release ('forward_sort', forward_sort)

  end subroutine tda_av_eigenval_bld

!===========================================================================
  subroutine tda_av_eigenval_via_super_bld
!---------------------------------------------------------------------------
! Description : Compute the TDA energies via the super-CI equation
! Description : The eigenvalues are sorted and printed.
!
! Created     : B. Mussard, June 2016 (from "ham_lin_av_bld")
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none
  
! local
  character(len=max_string_len_rout), save :: lhere = 'tda_av_eigenval_via_super_bld'
  integer                         :: i,info,lwork
  real(dp), allocatable           :: amat_local(:,:)
  real(dp), allocatable           :: ovlp_local(:,:)
  real(dp), allocatable           :: eigvec(:,:)
  real(dp), allocatable           :: eigval_r(:)
  real(dp), allocatable           :: eigval_i(:)
  real(dp), allocatable           :: eigval_denom(:)
  real(dp), allocatable           :: work(:)
  integer,  allocatable           :: backward_sort(:)
  integer,  allocatable           :: forward_sort(:)

! begin
  if (header_exe) then
    call object_create('tda_av_eigenval_via_super_r')
    call object_create('tda_av_eigenval_via_super_i')
    call object_error_define ('tda_av_eigenval_via_super_r','tda_av_eigenval_via_super_r_err')
    call object_error_define ('tda_av_eigenval_via_super_i','tda_av_eigenval_via_super_i_err')

    call object_needed('param_nb')

    call object_needed('super_amat_av')
    call object_needed('super_ovlp_psii_psij_av')

    return
  endif

  call object_alloc ('tda_av_eigenval_via_super_r',tda_av_eigenval_via_super_r,param_nb+1)
  call object_alloc ('tda_av_eigenval_via_super_i',tda_av_eigenval_via_super_i,param_nb+1)
  call object_alloc ('tda_av_eigenval_via_super_r_err',tda_av_eigenval_via_super_r_err,param_nb+1)
  call object_alloc ('tda_av_eigenval_via_super_i_err',tda_av_eigenval_via_super_i_err,param_nb+1)

! prepare arrays and output message
  call alloc('eigvec',       eigvec,       param_nb+1, param_nb+1)
  call alloc('eigval_r',     eigval_r,     param_nb+1)
  call alloc('eigval_i',     eigval_i,     param_nb+1)
  call alloc('eigval_denom', eigval_denom, param_nb+1)
  if (run_done.or.l_print_every_block) then
  write(6,*)
  write(6,'(a,1pd9.1)') 'Solving the TDA eigenvalue equation (via SUPER)'
  endif

! calculate optimal value of lwork
  lwork = 1
  call alloc('work', work, lwork)
  call alloc('amat_local', amat_local, param_nb+1, param_nb+1)
  call alloc('ovlp_local', ovlp_local, param_nb+1, param_nb+1)
  amat_local=super_amat_av
  ovlp_local=super_ovlp_psii_psij_av
  call dggev('N','V', &
             param_nb+1, amat_local,         &
             param_nb+1, ovlp_local,         &
             param_nb+1, eigval_r, eigval_i, eigval_denom, &
             eigvec, param_nb+1, eigvec, param_nb+1, &
             work, -1, info)
  if(info /=  0) call die(lhere, 'problem in dggev_alloc: info='+info+' /=  0')
  lwork=nint(work(1))
  call alloc('work', work, lwork)

! generalized eigenvalue problem
  call dggev('N','V', &
             param_nb+1, amat_local,         &
             param_nb+1, ovlp_local,         &
             param_nb+1, eigval_r, eigval_i, eigval_denom, &
             eigvec, param_nb+1, eigvec, param_nb+1, &
             work, lwork, info)
  if(info /=  0) call die(lhere, 'problem in dggev: info='+info+' /=  0 (compare to  param_nb='+param_nb+'+1)')

! scale eigenvalue
  do i = 1, param_nb+1
    eigval_r(i) = eigval_r(i) / eigval_denom(i)
    eigval_i(i) = eigval_i(i) / eigval_denom(i)
  enddo

! sort eigenvalues
  call sorting_real_im(eigval_r,eigval_i,forward_sort,backward_sort,param_nb+1) 

! the sorted eigenvalues are "linresp_av_eigenval"
  do i = 1, param_nb+1
   tda_av_eigenval_via_super_r (i) = eigval_r(backward_sort(i))
   tda_av_eigenval_via_super_i (i) = eigval_i(backward_sort(i))
  enddo

! output analysis
  call analysis_eigval_eigvec(eigvec,backward_sort, &
    tda_av_eigenval_via_super_r,tda_av_eigenval_via_super_r_err, &
    tda_av_eigenval_via_super_i,tda_av_eigenval_via_super_i_err,param_nb+1)

! release statements
  call release ('amat_local', amat_local)
  call release ('ovlp_local', ovlp_local)
  call release ('work', work)
  call release ('eigvec', eigvec)
  call release ('eigval_r', eigval_r)
  call release ('eigval_i', eigval_i)
  call release ('eigval_denom', eigval_denom)
  call release ('backward_sort', backward_sort)
  call release ('forward_sort', forward_sort)

  end subroutine tda_av_eigenval_via_super_bld

!===========================================================================
  subroutine hessian_av_eigenval_bld
!---------------------------------------------------------------------------
! Description : Solve the HESSIAN eigenvalue equation.
! Description : The eigenvalues are sorted and printed.
!
! Created     : B. Mussard, June 2016 (from "ham_lin_av_bld")
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none
  
! local
  character(len=max_string_len_rout), save :: lhere = 'hessian_av_eigenval_bld'
  integer                         :: i,info,lwork
  real(dp), allocatable           :: linresp_local(:,:)
  real(dp), allocatable           :: eigvec(:,:)
  real(dp), allocatable           :: eigval_r(:)
  real(dp), allocatable           :: eigval_i(:)
  integer , allocatable           :: backward_sort(:)
  integer , allocatable           :: forward_sort(:)
  real(dp), allocatable           :: work(:)

! begin
  if (header_exe) then
    call object_create('hessian_av_eigenval_r')
    call object_create('hessian_av_eigenval_i')
    call object_error_define ('hessian_av_eigenval_r','hessian_av_eigenval_r_err')
    call object_error_define ('hessian_av_eigenval_i','hessian_av_eigenval_i_err')

    call object_needed('param_nb')

    call object_needed('linresp_mat')

    return
  endif

  call object_alloc ('hessian_av_eigenval_r',hessian_av_eigenval_r,2*param_nb)
  call object_alloc ('hessian_av_eigenval_i',hessian_av_eigenval_i,2*param_nb)
  call object_alloc ('hessian_av_eigenval_r_err',hessian_av_eigenval_r_err,2*param_nb)
  call object_alloc ('hessian_av_eigenval_i_err',hessian_av_eigenval_i_err,2*param_nb)

! prepare arrays and output message
  call alloc('eigvec',       eigvec,       2*param_nb, 2*param_nb)
  call alloc('eigval_r',     eigval_r,     2*param_nb)
  call alloc('eigval_i',     eigval_i,     2*param_nb)
  if (run_done.or.l_print_every_block) then
  write(6,*)
  write(6,'(a,1pd9.1)') 'Solving the HESSIAN eigenvalue equation'
  endif

! calculate optimal value of lwork
  lwork = 1
  call alloc('work', work, lwork)
  call alloc('linresp_local', linresp_local, 2*param_nb, 2*param_nb)
  linresp_local=linresp_mat
  call dgeev('N','V', &
             2*param_nb, linresp_local,      &
             2*param_nb, eigval_r, eigval_i, &
             eigvec, 2*param_nb, eigvec, 2*param_nb, &
             work, -1, info)
  if(info /=  0) call die(lhere, 'problem in dgeev_alloc: info='+info+' /=  0')
  lwork=nint(work(1))
  call alloc('work', work, lwork)

! generalized eigenvalue problem
  call dgeev('N','V', &
             2*param_nb, linresp_local,      &
             2*param_nb, eigval_r, eigval_i, &
             eigvec, 2*param_nb, eigvec, 2*param_nb, &
             work, lwork, info)
  if(info /=  0) call die(lhere, 'problem in dgeev: info='+info+' /=  0 (compare to  2*param_nb='+2*param_nb+')')

! sort eigenvalues
  call sorting_real_im(eigval_r,eigval_i,forward_sort,backward_sort,2*param_nb) 

! the sorted eigenvalues are "linresp_av_eigenval"
  do i = 1, 2*param_nb
   hessian_av_eigenval_r (i) = eigval_r(backward_sort(i))
   hessian_av_eigenval_i (i) = eigval_i(backward_sort(i))
  enddo

! output analysis
  call analysis_eigval_eigvec(eigvec,backward_sort, &
    hessian_av_eigenval_r,hessian_av_eigenval_r_err, &
    hessian_av_eigenval_i,hessian_av_eigenval_i_err,2*param_nb)

! release statements
  call release ('linresp_local', linresp_local)
  call release ('work', work)
  call release ('eigvec', eigvec)
  call release ('eigval_r', eigval_r)
  call release ('eigval_i', eigval_i)
  call release ('backward_sort', backward_sort)
  call release ('forward_sort', forward_sort)

  end subroutine hessian_av_eigenval_bld

!===========================================================================
  subroutine real_hessian_av_eigenval_bld
!---------------------------------------------------------------------------
! Description : Solve the REAL HESSIAN eigenvalue equation.
! Description : The eigenvalues are sorted and printed.
!
! Created     : B. Mussard, June 2016 (from "ham_lin_av_bld")
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none
  
! local
  character(len=max_string_len_rout), save :: lhere = 'real_hessian_av_eigenval_bld'
  integer                         :: i,info,lwork
  real(dp), allocatable           :: abmat_local(:,:)
  real(dp), allocatable           :: eigvec(:,:)
  real(dp), allocatable           :: eigval_r(:)
  real(dp), allocatable           :: eigval_i(:)
  integer , allocatable           :: backward_sort(:)
  integer , allocatable           :: forward_sort(:)
  real(dp), allocatable           :: work(:)

! begin
  if (header_exe) then
    call object_create('real_hessian_av_eigenval_r')
    call object_create('real_hessian_av_eigenval_i')
    call object_error_define ('real_hessian_av_eigenval_r','real_hessian_av_eigenval_r_err')
    call object_error_define ('real_hessian_av_eigenval_i','real_hessian_av_eigenval_i_err')

    call object_needed('param_nb')

    call object_needed('amat_av')
    call object_needed('bmat_av')

    return
  endif

  call object_alloc ('real_hessian_av_eigenval_r',real_hessian_av_eigenval_r,param_nb)
  call object_alloc ('real_hessian_av_eigenval_i',real_hessian_av_eigenval_i,param_nb)
  call object_alloc ('real_hessian_av_eigenval_r_err',real_hessian_av_eigenval_r_err,param_nb)
  call object_alloc ('real_hessian_av_eigenval_i_err',real_hessian_av_eigenval_i_err,param_nb)

! prepare arrays and output message
  call alloc('eigvec',       eigvec,       param_nb, param_nb)
  call alloc('eigval_r',     eigval_r,     param_nb)
  call alloc('eigval_i',     eigval_i,     param_nb)
  if (run_done.or.l_print_every_block) then
  write(6,*)
  write(6,'(a,1pd9.1)') 'Solving the REAL HESSIAN eigenvalue equation'
  endif

! calculate optimal value of lwork
  lwork = 1
  call alloc('work', work, lwork)
  call alloc('abmat_local', abmat_local, param_nb, param_nb)
  abmat_local=amat_av+bmat_av
  call dgeev('N','V', &
             param_nb, abmat_local,    &
             param_nb, eigval_r, eigval_i, &
             eigvec, param_nb, eigvec, param_nb, &
             work, -1, info)
  if(info /=  0) call die(lhere, 'problem in dgeev_alloc: info='+info+' /=  0')
  lwork=nint(work(1))
  call alloc('work', work, lwork)

! generalized eigenvalue problem
  call dgeev('N','V', &
             param_nb, abmat_local,    &
             param_nb, eigval_r, eigval_i, &
             eigvec, param_nb, eigvec, param_nb, &
             work, lwork, info)
  if(info /=  0) call die(lhere, 'problem in dgeev: info='+info+' /=  0 (compare to  param_nb='+param_nb+')')

! sort eigenvalues
  call sorting_real_im(eigval_r,eigval_i,forward_sort,backward_sort,param_nb) 

! the sorted eigenvalues are "linresp_av_eigenval"
  do i = 1, param_nb
   real_hessian_av_eigenval_r (i) = eigval_r(backward_sort(i))
   real_hessian_av_eigenval_i (i) = eigval_i(backward_sort(i))
  enddo

! output analysis
  call analysis_eigval_eigvec(eigvec,backward_sort, &
    real_hessian_av_eigenval_r,real_hessian_av_eigenval_r_err,&
    real_hessian_av_eigenval_i,real_hessian_av_eigenval_i_err,param_nb)

! release statements
  call release ('abmat_local', abmat_local)
  call release ('work', work)
  call release ('eigvec', eigvec)
  call release ('eigval_r', eigval_r)
  call release ('eigval_i', eigval_i)
  call release ('backward_sort', backward_sort)
  call release ('forward_sort', forward_sort)

  end subroutine real_hessian_av_eigenval_bld

!===========================================================================
  subroutine ovlp_mat_bld
!---------------------------------------------------------------------------
! Description : Construct the OVERLAP 2*n-matrix from the overlap matrix
! Description :
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! begin
  if (header_exe) then
    call object_create('ovlp_mat')

    call object_needed('param_nb')
    call object_needed('ovlp_psii_psij_av')

    return
  endif

  call object_alloc('ovlp_mat',ovlp_mat,2*param_nb,2*param_nb)

  ovlp_mat=0.d0
  ovlp_mat(1:param_nb,1:param_nb)=ovlp_psii_psij_av
  if (l_compare_linresp_and_optlin) then
    ovlp_mat(param_nb+1:2*param_nb,param_nb+1:2*param_nb)=+ovlp_psii_psij_av
  else
    ovlp_mat(param_nb+1:2*param_nb,param_nb+1:2*param_nb)=-ovlp_psii_psij_av
  endif

  end subroutine ovlp_mat_bld

!===========================================================================
  subroutine super_ovlp_psii_psij_av_bld
!---------------------------------------------------------------------------
! Description : Construct the OVERLAP (n+1)-matrix from the overlap matrix
! Description :
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer :: i,j

! begin
  if (header_exe) then
    call object_create('super_ovlp_psii_psij_av')

    call object_needed('param_nb')
    call object_needed('ovlp_psii_psij_av')

    return
  endif

  call object_alloc('super_ovlp_psii_psij_av',super_ovlp_psii_psij_av,param_nb+1,param_nb+1)

  super_ovlp_psii_psij_av(1,1)=1.d0
  super_ovlp_psii_psij_av(2:,2:)=ovlp_psii_psij_av

  if (run_done.or.l_print_every_block) then
  if (l_print) then
  do i=1,param_nb+1
  do j=1,param_nb+1
    write(6,*) '/print_too_much/super_ovlp',i,j,super_ovlp_psii_psij_av(i,j)
  enddo
  enddo
  endif
  endif

  end subroutine super_ovlp_psii_psij_av_bld

!===========================================================================
  subroutine linresp_mat_bld
!---------------------------------------------------------------------------
! Description : Construct the ABBA 2*n-matrix from A and B matrices
! Description :
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! begin
  if (header_exe) then
    call object_create('linresp_mat')

    call object_needed('param_nb')
    call object_needed('amat_av')

    return
  endif

  call object_alloc('linresp_mat',linresp_mat,2*param_nb,2*param_nb)
  linresp_mat=0.d0
  linresp_mat(1:param_nb,1:param_nb)=amat_av
  linresp_mat(param_nb+1:2*param_nb,param_nb+1:2*param_nb)=amat_av
  if ((.not.l_compare_linresp_and_optlin).and.(.not.l_compare_to_tda)) then
    call object_provide('bmat_av')
    linresp_mat(param_nb+1:2*param_nb,1:param_nb)=bmat_av
    linresp_mat(1:param_nb,param_nb+1:2*param_nb)=bmat_av
  endif

  end subroutine linresp_mat_bld

!===========================================================================
  subroutine ovlp_psii_psij_av_bld
!---------------------------------------------------------------------------
! Description : <Psi_i|Psi_j> overlap 
! Description : (this is inspired by "ham_eigval_av_bld",
! Description :  which does many other things useless in this context)
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer :: i,j
  real(dp), allocatable           :: ovlp_eigval(:)
  real(dp), allocatable           :: ovlp_eigvec(:,:)

! begin
  if (header_exe) then
    call object_create('ovlp_psii_psij_av')

    call object_needed('param_nb')
    call object_needed('dpsi_dpsi_covar')

    return
  endif

  call object_alloc('ovlp_psii_psij_av',ovlp_psii_psij_av,param_nb,param_nb)

! Eq(53d) of JCP 126 084102 (2007)
  do i = 1, param_nb
    do j = i, param_nb
      ovlp_psii_psij_av(i,j) = dpsi_dpsi_covar(i,j)
      if (i.ne.j) then
       ovlp_psii_psij_av(j,i) = ovlp_psii_psij_av(i,j)
      endif
    enddo
  enddo

  if (run_done.or.l_print_every_block) then
  if (l_print) then
  do i= 1,param_nb
  write(6,*) '/print_too_much/ovlp',i,ovlp_psii_psij_av(i,:)
  enddo
  call alloc('ovlp_eigvec', ovlp_eigvec, param_nb, param_nb)
  call alloc('ovlp_eigval', ovlp_eigval, param_nb)
  call eigensystem(ovlp_psii_psij_av, ovlp_eigvec, ovlp_eigval, param_nb)
  write(6,*)
  write(6,'(a)') '/print_too_much/eigenvalues:'
  do i = 1, param_nb
    write(6,'(a,i4,a,es15.8)') '/print_too_much/eigenvalue # ',i,': ',ovlp_eigval(i)
  enddo
  call release('ovlp_eigvec',ovlp_eigvec)
  call release('ovlp_eigval',ovlp_eigval)
  endif
  endif

  end subroutine ovlp_psii_psij_av_bld

!===========================================================================
  subroutine amat_av_bld
!---------------------------------------------------------------------------
! Description : <Psi_i|H|Psi_j>
! Description : (this is inspired by "ham_eigval_av_bld",
! Description :  which does many other things useless in this context)
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer :: i,j,ij

! begin
  if (header_exe) then
    call object_create('amat_av')

    call object_needed('param_nb')
    call object_needed('param_pairs')

    call object_needed('dpsi_dpsi_av')
    call object_needed('dpsi_dpsi_eloc_av')
    call object_needed('dpsi_av')
    call object_needed('dpsi_eloc_av')
    call object_needed('eloc_av')
    call object_needed('dpsi_deloc_covar')

    return
  endif

  call object_alloc('amat_av',amat_av,param_nb,param_nb)

! Eq(54d) of JCP 126 084102 (2007)
  do j=1,param_nb
    do i=1,param_nb
      ij = param_pairs(i,j)
      if (l_compare_linresp_and_optlin) then
        amat_av(i,j)=dpsi_dpsi_eloc_av(ij)         &
                   - dpsi_av(j)*dpsi_eloc_av(i)    &
                   - dpsi_av(i)*dpsi_eloc_av(j)    &
                   + dpsi_av(i)*dpsi_av(j)*eloc_av &
                   + dpsi_deloc_covar(i,j)
      else
        amat_av(i,j)=dpsi_dpsi_eloc_av(ij)         &
               -     dpsi_av(j)*dpsi_eloc_av(i)    &
               -     dpsi_av(i)*dpsi_eloc_av(j)    &
               +2.d0*dpsi_av(i)*dpsi_av(j)*eloc_av &
               +     dpsi_deloc_covar(i,j)         &
               -     eloc_av*dpsi_dpsi_av(param_pairs(i,j))
      endif
    enddo
  enddo

  if (l_symm_amat) then
    amat_av=0.5d0*(amat_av+transpose(amat_av))
  endif

  if (run_done.or.l_print_every_block) then
  if (l_print) then
  do i=1,param_nb
  do j=1,param_nb
    write(6,*) '/print_too_much/amat',i,j,amat_av(i,j)
  enddo
  enddo
  endif
  endif

  end subroutine amat_av_bld

!===========================================================================
  subroutine super_amat_av_bld
!---------------------------------------------------------------------------
! Description : <Psi_i|H|Psi_j>
! Description : (this is inspired by "ham_eigval_av_bld",
! Description :  which does many other things useless in this context)
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer :: i,j,ij

! begin
  if (header_exe) then
    call object_create('super_amat_av')

    call object_needed('param_nb')
    call object_needed('param_pairs')

    call object_needed('eloc_av')
    call object_needed('dpsi_eloc_covar')
    call object_needed('deloc_av')
    call object_needed('amat_av')

    return
  endif

  call object_alloc('super_amat_av',super_amat_av,param_nb+1,param_nb+1)

  super_amat_av(1,1)=eloc_av

  do i=1,param_nb
    super_amat_av(i+1,1)=dpsi_eloc_covar(i)
    super_amat_av(1,i+1)=dpsi_eloc_covar(i) + deloc_av(i)
  enddo

  super_amat_av(2:,2:)=amat_av
  do j=1,param_nb
    do i=1,param_nb
      ij = param_pairs(i,j)
      super_amat_av(i+1,j+1)=dpsi_dpsi_eloc_av(ij)         &
                           - dpsi_av(j)*dpsi_eloc_av(i)    &
                           - dpsi_av(i)*dpsi_eloc_av(j)    &
                           + dpsi_av(i)*dpsi_av(j)*eloc_av &
                           + dpsi_deloc_covar(i,j)
    enddo
  enddo
  if (l_symm_amat) then
    super_amat_av(2:,2:)=0.5d0*(super_amat_av(2:,2:)+transpose(super_amat_av(2:,2:)))
  endif

  if (run_done.or.l_print_every_block) then
  if (l_print) then
  do i=1,param_nb+1
  do j=1,param_nb+1
    write(6,*) '/print_too_much/super_amat',i,j,super_amat_av(i,j)
  enddo
  enddo
  endif
  endif

  end subroutine super_amat_av_bld

!===========================================================================
  subroutine bmat_av_bld
!---------------------------------------------------------------------------
! Description : <Psi_ij|H|Psi_0>
! Description : 
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  integer :: i,j,ij

! begin
  if (header_exe) then
    call object_create('bmat_av')

    call object_needed('param_nb')
    call object_needed('param_pairs')

    call object_needed('dpsi_av')
    call object_needed('dpsi_eloc_av')
    call object_needed('eloc_av')
    call object_needed('d2psi_av')
    call object_needed('d2psi_eloc_av')

    return
  endif

  call object_alloc('bmat_av',bmat_av,param_nb,param_nb)

  do j=1,param_nb
    do i=1,param_nb
      ij = param_pairs(i,j)
      bmat_av(i,j)=(2*dpsi_av(i)*dpsi_av(j)-d2psi_av(ij))*eloc_av &
                  - dpsi_av(j)*dpsi_eloc_av(i)                    &
                  - dpsi_av(i)*dpsi_eloc_av(j)                    &
                  + d2psi_eloc_av(ij)
    enddo
  enddo

  if (run_done.or.l_print_every_block) then
  if (l_print) then
  do i=1,param_nb
  do j=1,param_nb
    write(6,*) '/print_too_much/bmat',i,j,bmat_av(i,j)
  enddo
  enddo
  endif
  endif

  end subroutine bmat_av_bld

!===========================================================================
  subroutine sorting_real_im(eigval_r,eigval_i,forward_sort,backward_sort,n)
!---------------------------------------------------------------------------
! Description : Sorting out eigenvalues
! Description : 
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input/output
  real(dp), allocatable           :: eigval_r(:)
  real(dp), allocatable           :: eigval_i(:)
  integer,  allocatable           :: backward_sort(:)
  integer,  allocatable           :: forward_sort(:)
  integer                         :: n

! local
  integer :: i,j,temp

! begin
! "backward_sort"
! is the map from sorted eigenvalues to original eigenvalues
  call alloc('backward_sort', backward_sort, n)
  do i = 1, n
    backward_sort(i) = i
  enddo
  do i = 1, n
    do j = i+1, n
      if    ((eigval_r(backward_sort(j))-eigval_r(backward_sort(i))<-0.0001) &
       .or. ((eigval_r(backward_sort(j))-eigval_r(backward_sort(i))<+0.0001) &
       .and. (eigval_i(backward_sort(j))>eigval_i(backward_sort(i))))) then
        temp = backward_sort(i)
        backward_sort(i) = backward_sort(j)
        backward_sort(j) = temp
      endif
    enddo
  enddo
! "forward_sort"
! is the map from original eigenvalues to sorted eigenvalues
  call alloc('forward_sort', forward_sort, n)
  do i = 1, n
   forward_sort(backward_sort(i)) = i
  enddo
  if (run_done.or.l_print_every_block) then
  write(6,'(''indexes      ='',10000i6)') (i,i=1,n)
  write(6,'(''backward_sort='',10000i6)') backward_sort(1:n)
  write(6,'(''forward_sort ='',10000i6)') forward_sort(1:n)
  endif

  end subroutine sorting_real_im
         
!===========================================================================
  subroutine analysis_eigval_eigvec(eigvec,sorting,eigval_r,err_r,eigval_i,err_i,n)
!---------------------------------------------------------------------------
! Description : Analysis of eigval/eigvec
! Description : 
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input/output
  real(dp), allocatable           :: eigvec(:,:)
  real(dp), allocatable           :: eigval_r(:)
  real(dp), allocatable           :: err_r(:)
  real(dp), allocatable           :: eigval_i(:)
  real(dp), allocatable           :: err_i(:)
  integer,  allocatable           :: sorting(:)
  integer                         :: n

! local
  integer                         :: i,j,k,nunique,i_eigval,j_eigval,i_pos,jorb
  integer,  allocatable           :: position_in_array(:)
  integer,  allocatable           :: degeneracy(:)
  real(dp), allocatable           :: zero_array(:)
  real(dp)                        :: contrib_max
  real(dp), allocatable           :: contrib(:,:)
  logical, allocatable            :: went_through(:)
  character(len=9), allocatable   :: formt0(:)
  character(len=7)                :: formt1
  character(len=9), allocatable   :: formt2(:)
  character(len=15), allocatable  :: formt3(:)

! begin
  if (header_exe) then
    call object_needed('param_pairs')
    call object_needed('dpsi_av')
    call object_needed('dpsi_dpsi_av')
    
    return
  endif

  if (run_done.or.l_print_every_block) then

! Gather information on unique eigenvalues (real and im part)
! and their original position and degeneracy
  call alloc('went_through',  went_through,  n)
  call alloc('position_in_array', position_in_array, n)
  call alloc('degeneracy', degeneracy, n)
  nunique=0
  went_through=.false.
  do i_eigval=1,n
    if (went_through(i_eigval)) cycle
    if ((n.eq.2*param_nb).and.(eigval_r(i_eigval).lt.0.d0)) cycle
    went_through(i_eigval)=.true.
    nunique=nunique+1
    position_in_array(nunique)=i_eigval
    degeneracy(nunique)=1

    ! search for all repetitions of this newly found unique eigenvalue
    do j_eigval=i_eigval+1,n
      if (went_through(j_eigval)) cycle
      ! real1=real2, im1=im2, X1=X2 Y1=Y2
      if((abs(eigval_r(i_eigval)-eigval_r(j_eigval)).le.0.001).and. &
        &(abs(eigval_i(i_eigval)-eigval_i(j_eigval)).le.0.001).and. &
        &(is_equal(eigvec(:,sorting(i_eigval)),&
                 & eigvec(:,sorting(j_eigval)),1.0d-2).eq.0)) then
        went_through(j_eigval)=.true.
        degeneracy(nunique)=degeneracy(nunique)+1
      !! real1=real2, im1=-im2!=0, this is only an imaginary pair
      !elseif((abs(eigval_r(i_eigval)-eigval_r(j_eigval)).le.0.001).and.     &
      !      &(abs(eigval_i(i_eigval)+eigval_i(j_eigval)).le.0.001).and.     &
      !      &(abs(eigval_i(i_eigval)).ge.0.001)) then
      !  went_through(j_eigval)=.true.
      !  degeneracy(nunique)=degeneracy(nunique)+1
      endif
      if (n.eq.2*param_nb) then
      ! real1=-real2, im1=im2=0, (X1=Y2 Y1=X2 or X1=-Y2 Y1=-X2)
      if    ((abs(eigval_r(i_eigval)+eigval_r(j_eigval)).le.0.001).and.     &
            &(abs(eigval_i(i_eigval)).le.0.001).and.                        &
            &(abs(eigval_i(j_eigval)).le.0.001).and.                        &
            &(((is_equal(eigvec(:n/2  ,sorting(i_eigval)),                  &
                       & eigvec(n/2+1:,sorting(j_eigval)),1.0d-2).eq.0).and.&
            &  (is_equal(eigvec(n/2+1:,sorting(i_eigval)),                  &
                       & eigvec(:n/2  ,sorting(j_eigval)),1.0d-2).eq.0)).or.&
            & ((is_equal(eigvec(:n/2  ,sorting(i_eigval)),                  &
                       &-eigvec(n/2+1:,sorting(j_eigval)),1.0d-2).eq.0).and.&
            &  (is_equal(eigvec(n/2+1:,sorting(i_eigval)),                  &
                       &-eigvec(:n/2  ,sorting(j_eigval)),1.0d-2).eq.0)))) then
        went_through(j_eigval)=.true.
        degeneracy(nunique)=degeneracy(nunique)+1
       !if (eigval_r(i_eigval).lt.0.000) position_in_array(nunique)=j_eigval
      !! real1=-real2, im1=im2!=0, X1=Y2 Y1=X2
      !elseif((abs(eigval_r(i_eigval)+eigval_r(j_eigval)).le.0.001).and.     &
      !      &(abs(eigval_i(i_eigval)-eigval_i(j_eigval)).le.0.001).and.     &
      !      &(abs(eigval_i(i_eigval)).ge.0.001).and.                        &
      !      &(((is_equal(eigvec(:n/2  ,sorting(i_eigval)),                  &
      !                 & eigvec(n/2+1:,sorting(j_eigval)),1.0d-2).eq.0).and.&
      !      &  (is_equal(eigvec(n/2+1:,sorting(i_eigval)),                  &
      !                 & eigvec(:n/2  ,sorting(j_eigval)),1.0d-2).eq.0)).or.&
      !      & ((is_equal(eigvec(:n/2  ,sorting(i_eigval)),                  &
      !                 &-eigvec(n/2+1:,sorting(j_eigval)),1.0d-2).eq.0).and.&
      !      &  (is_equal(eigvec(n/2+1:,sorting(i_eigval)),                  &
      !                 &-eigvec(:n/2  ,sorting(j_eigval)),1.0d-2).eq.0)))) then
      !  went_through(j_eigval)=.true.
      !  degeneracy(nunique)=degeneracy(nunique)+1
      ! !if (eigval_r(i_eigval).lt.0.000) position_in_array(nunique)=j_eigval
      !! real1=-real2, im1=-im2!=0, this is only an imaginary pair
      !elseif((abs(eigval_r(i_eigval)+eigval_r(j_eigval)).le.0.001).and.     &
      !      &(abs(eigval_i(i_eigval)+eigval_i(j_eigval)).le.0.001).and.     &
      !      &(abs(eigval_i(i_eigval)).ge.0.001)) then
      !  went_through(j_eigval)=.true.
      !  degeneracy(nunique)=degeneracy(nunique)+1
      ! !if (eigval_r(i_eigval).lt.0.000) position_in_array(nunique)=j_eigval
      endif
      endif
    enddo
  enddo
  call release('went_through', went_through)
      
! Print the unique eigenvalues...
! ...and (if asked for) the eigenvectors
! - in the case where matrices are of dimension "param_nb":
!   print the eigenvector where:
!     - the type of parameter for each component is written 
!      (in the case of an orbital parameter: additional info on p->q excitation)
!       This is "formt3"
!     - the three major component are highlighted
!       This is "formt2"
! - in the case where matrices are of dimension "2*param_nb":
!       add information on X=Y, X=-Y, X=0, Y=0, or UNKNOWN situation
!       This is "formt1"
  write(6,'(a,i5)') 'Sorted (complex) (unique) eigenvalues:',nunique

  call get_formt0(formt0,n,eigval_r,eigvec,sorting,position_in_array)
  call get_formt3(formt3,n)

  do i = 1, nunique
    write(6,'(a,i8,a,a,4(f12.6,a),i5,a)') 'eigenvalue #',i,': ',formt0(position_in_array(i)),&
      & eigval_r(position_in_array(i)),'(+/-',err_r(position_in_array(i)),') +',&
      & eigval_i(position_in_array(i)),'(+/-',err_i(position_in_array(i)),') i (',&
      & degeneracy(i),')'

    call get_formt1(i,formt1,n,eigvec,sorting,position_in_array)
    call get_formt2(i,formt2,n,eigvec,sorting,position_in_array)

    ! write all out
    if (l_print_eigenvec) then
      if ((n.eq.param_nb).or.(n.eq.param_nb+1)) then
        do j=1,n
          if (abs(eigval_i(position_in_array(i))).eq.0.d0) then
             write(6,'(a,a,i5,a,f11.3,a,a)') &
               & 'eigenvec',formt1,i,':          ',&
               &  eigvec(j, sorting(position_in_array(i))),&
               &  formt2(j),formt3(j)
          elseif (eigval_i(position_in_array(i)).gt.0.d0) then
            write(6,'(a,a,i5,a,f11.3,a,f11.3,a,a)') &
               & 'eigenvec',formt1,i,':          ',&
               &  eigvec(j ,sorting(position_in_array(i)))  ,' + ',eigvec(j,sorting(position_in_array(i))+1),&
               &  formt2(j),formt3(j)
          elseif (eigval_i(position_in_array(i)).lt.0.d0) then
            write(6,'(a,a,i5,a,f11.3,a,f11.3,a,a)') &
               & 'eigenvec',formt1,i,':          ',&
               &  eigvec(j ,sorting(position_in_array(i))-1),' + ',eigvec(j,sorting(position_in_array(i))  ),&
               &  formt2(j),formt3(j)
          endif
        enddo
      elseif (n.eq.2*param_nb) then
        do j=1,n/2
          if (abs(eigval_i(position_in_array(i))).eq.0.d0) then
            write(6,'(a,a,i5,a,2(f11.3,a,a))') &
               & 'eigenvec',formt1,i,':          ',&
               &  eigvec(j     ,sorting(position_in_array(i))),&
               &  formt2(j)    ,formt3(j), &
               &  eigvec(j+n/2 ,sorting(position_in_array(i))),&
               &  formt2(j+n/2),formt3(j+n/2)
          elseif (eigval_i(position_in_array(i)).gt.0.d0) then
            write(6,'(a,a,i5,a,2(f11.3,a,f11.3,a,a))') &
               & 'eigenvec',formt1,i,':          ',&
               &  eigvec(j     ,sorting(position_in_array(i))  ),' + ',eigvec(j    ,sorting(position_in_array(i))+1),&
               &  formt2(j)    ,formt3(j), &
               &  eigvec(j+n/2 ,sorting(position_in_array(i))  ),' + ',eigvec(j+n/2,sorting(position_in_array(i))+1),&
               &  formt2(j+n/2),formt3(j+n/2)
          elseif (eigval_i(position_in_array(i)).lt.0.d0) then
            write(6,'(a,a,i5,a,2(f11.3,a,f11.3,a,a))') &
               & 'eigenvec',formt1,i,':          ',&
               &  eigvec(j     ,sorting(position_in_array(i))-1),' + ',eigvec(j    ,sorting(position_in_array(i))  ),&
               &  formt2(j)    ,formt3(j), &                   
               &  eigvec(j+n/2 ,sorting(position_in_array(i))-1),' + ',eigvec(j+n/2,sorting(position_in_array(i))  ),&
               &  formt2(j+n/2),formt3(j+n/2)
          endif
        enddo
      endif
    endif
  enddo

  if (param_orb_nb.ne.0) then
  write(6,*) ""
  write(6,*) "Here are shown the contributions,"
  write(6,*) "for each orb'th orbital excitation of"
  write(6,*) "  <\Psi_orb|\Psi^n>=  \sum_i \Delta p^n_i <\Psi_orb|\Psi_i>"
  write(6,*) " where \Psi^n is the wavefunction constructed by the n'th eigenvector \Delta p^n_i"
  write(6,*) "Here are only shown the dominant contributions, scaled to the largest component"
  write(6,*) ""

  call alloc('contrib',contrib,n,param_orb_nb)
  contrib_max=0.d0
  jorb=0
  do j=1,param_nb
    if (is_param_type_orb(j)) then
      jorb=jorb+1
      do i=1,n
        contrib(i,jorb)=0
        if (n.eq.param_nb) then
          do k=1,n
            contrib(i,jorb)=contrib(i,jorb)+eigvec(k,i)  *dpsi_dpsi_av(param_pairs(j,k))
          enddo
        elseif (n.eq.param_nb+1) then
          contrib(i,jorb)=eigvec(1,i)
          do k=1,n-1
            contrib(i,jorb)=contrib(i,jorb)+eigvec(k+1,i)  *dpsi_av(k)
          enddo
        elseif (n.eq.2*param_nb) then
          do k=1,n/2
            contrib(i,jorb)=contrib(i,jorb)+eigvec(k,i)    *dpsi_dpsi_av(param_pairs(j,k))
            contrib(i,jorb)=contrib(i,jorb)+eigvec(k+n/2,i)*dpsi_dpsi_av(param_pairs(j,k))
          enddo
        endif
        if (abs(contrib(i,jorb)).ge.contrib_max) contrib_max=abs(contrib(i,jorb))
      enddo
    endif
  enddo
  contrib=contrib/contrib_max

  call alloc('zero_array',zero_array,param_orb_nb)
  do j=1,nunique
  if (abs(eigval_i(position_in_array(j))).eq.0.d0) then
    if (is_equal(contrib(sorting(position_in_array(j)),:),zero_array,0.2d0).ne.0) then
      write(6,'(a,i8,a,4(f12.6,a),i5,a)') &
        & 'eigenvalue #',j,': ',&
        &  eigval_r(position_in_array(j)),'(+/-',err_r(position_in_array(j)),') +',&
        &  eigval_i(position_in_array(j)),'(+/-',err_i(position_in_array(j)),') i (',&
        &  degeneracy(j),')'
      do jorb=1,param_orb_nb
        write(*,'(a,i2,a,i2,a,f11.3)') &
          & '(orbital',ex_orb_1st_lab(ex_orb_ind(jorb)),&
          & '->',ex_orb_2nd_lab(ex_orb_ind(jorb)),')',&
          & contrib(sorting(position_in_array(j)),jorb)
      enddo
    endif
  elseif (eigval_i(position_in_array(j)).gt.0.d0) then
    if((is_equal(contrib(sorting(position_in_array(j))  ,:),zero_array,0.2d0).ne.0).or.&
      &(is_equal(contrib(sorting(position_in_array(j))+1,:),zero_array,0.2d0).ne.0)) then
      write(6,'(a,i8,a,4(f12.6,a),i5,a)') &
        & 'eigenvalue #',j,': ',&
        &  eigval_r(position_in_array(j)),'(+/-',err_r(position_in_array(j)),') +',&
        &  eigval_i(position_in_array(j)),'(+/-',err_i(position_in_array(j)),') i (',&
        &  degeneracy(j),')'
      do jorb=1,param_orb_nb
        write(*,'(a,i2,a,i2,a,2f11.3)') &
          & '(orbital',ex_orb_1st_lab(ex_orb_ind(jorb)),&
          & '->',ex_orb_2nd_lab(ex_orb_ind(jorb)),')',&
          & contrib(sorting(position_in_array(j)),jorb),contrib(sorting(position_in_array(j))+1,jorb)
      enddo
    endif
  elseif (eigval_i(position_in_array(j)).lt.0.d0) then
    if((is_equal(contrib(sorting(position_in_array(j))-1,:),zero_array,0.2d0).ne.0).or.&
      &(is_equal(contrib(sorting(position_in_array(j))  ,:),zero_array,0.2d0).ne.0)) then
      write(6,'(a,i8,a,4(f12.6,a),i5,a)') &
        & 'eigenvalue #',j,': ',&
        & eigval_r(position_in_array(j)),'(+/-',err_r(position_in_array(j)),') +',&
        & eigval_i(position_in_array(j)),'(+/-',err_i(position_in_array(j)),') i (',&
        & degeneracy(j),')'
      do jorb=1,param_orb_nb
        write(*,'(a,i2,a,i2,a,2f11.3)') &
          & '(orbital',ex_orb_1st_lab(ex_orb_ind(jorb)),&
          & '->',ex_orb_2nd_lab(ex_orb_ind(jorb)),')',&
          & contrib(sorting(position_in_array(j))-1,jorb),contrib(sorting(position_in_array(j)),jorb)
      enddo
    endif
  endif
  enddo
  endif

  call release ('position_in_array', position_in_array)
  call release ('degeneracy', degeneracy)

  endif

  end subroutine analysis_eigval_eigvec

!===========================================================================
  subroutine get_formt1(i,formt1,n,eigvec,sorting,position_in_array)
!---------------------------------------------------------------------------
! Description :
! Description : 
!
! Created     : B. Mussard, October 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input/output
  character(len=7), intent(out)              :: formt1
  integer, intent(in)                        :: i,n
  real(dp), allocatable,intent(in)           :: eigvec(:,:)
  integer,  allocatable,intent(in)           :: sorting(:)
  integer,  allocatable,intent(in)           :: position_in_array(:)

! local
  real(dp), allocatable           :: zero_array(:)

  if (l_print_eigenvec) then
    ! information on X=Y, etc...
    if ((n.eq.param_nb).or.(n.eq.param_nb+1)) then
      formt1='       '
    elseif (n.eq.2*param_nb) then
      call alloc('zero_array',zero_array,n/2)
      if     (is_equal(eigvec(:n/2   ,sorting(position_in_array(i))),&
                     & eigvec(n/2+1:n,sorting(position_in_array(i))),1.0d-2).eq.0) then
       formt1='(EQU)  '
      elseif (is_equal(eigvec(:n/2   ,sorting(position_in_array(i))),&
                    & -eigvec(n/2+1:n,sorting(position_in_array(i))),1.0d-2).eq.0) then
       formt1='(OPP)  '
      elseif (is_equal(eigvec(n/2+1:n,sorting(position_in_array(i))),&
                    &  zero_array,1.0d-2).eq.0) then
       formt1='(TDAX) '
      elseif (is_equal(eigvec(:n/2   ,sorting(position_in_array(i))),&
                    &  zero_array,1.0d-2).eq.0) then
       formt1='(TDAY) '
      else
       formt1='(UNK)  '
      endif
    endif
  endif

  end subroutine get_formt1

!===========================================================================
  subroutine get_formt2(i,formt2,n,eigvec,sorting,position_in_array)
!---------------------------------------------------------------------------
! Description :
! Description : 
!
! Created     : B. Mussard, October 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input/output
  character(len=9), allocatable, intent(out) :: formt2(:)
  integer, intent(in)                        :: i,n
  real(dp), allocatable,intent(in)           :: eigvec(:,:)
  integer,  allocatable,intent(in)           :: sorting(:)
  integer,  allocatable,intent(in)           :: position_in_array(:)

! local
  integer                         :: j,maj1pos,maj2pos,maj3pos,maj4pos
  real(dp)                        :: maj1,maj2,maj3,maj4

  if (l_print_eigenvec) then
    call alloc('formt2',formt2,n)
    maj1=0
    do j=1,n
      if    (abs(eigvec(j, sorting(position_in_array(i)))).ge.maj1) then
        maj1=abs(eigvec(j, sorting(position_in_array(i))))
        maj1pos=j
      endif
    enddo
    maj2=0
    do j=1,n
      if   ((abs(eigvec(j, sorting(position_in_array(i)))).ge.maj2).and.&
           &(abs(eigvec(j, sorting(position_in_array(i)))).lt.maj1)) then
        maj2=abs(eigvec(j, sorting(position_in_array(i))))
        maj2pos=j
      endif
    enddo
    maj3=0
    do j=1,n
      if   ((abs(eigvec(j, sorting(position_in_array(i)))).ge.maj3).and.&
           &(abs(eigvec(j, sorting(position_in_array(i)))).lt.maj2)) then
        maj3=abs(eigvec(j, sorting(position_in_array(i))))
        maj3pos=j
      endif
    enddo
    maj4=0
    do j=1,n
      if   ((abs(eigvec(j, sorting(position_in_array(i)))).ge.maj4).and.&
           &(abs(eigvec(j, sorting(position_in_array(i)))).lt.maj3)) then
        maj4=abs(eigvec(j, sorting(position_in_array(i))))
        maj4pos=j
      endif
    enddo
    do j=1,n
      if ((j.eq.maj1pos).or.(j.eq.maj2pos).or.(j.eq.maj3pos).or.(j.eq.maj4pos)) then
        formt2(j)=' [major] '
      else
        formt2(j)='         '
      endif
    enddo
  endif

  end subroutine get_formt2

!===========================================================================
  subroutine get_formt3(formt3,n)
!---------------------------------------------------------------------------
! Description :
! Description : 
!
! Created     : B. Mussard, October 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input/output
  character(len=15), allocatable,intent(out) :: formt3(:)
  integer, intent(in)                        :: n

! local
  integer :: j,jorb

  if (l_print_eigenvec) then
    call alloc('formt3',formt3,n)
    if (n.eq.param_nb) then
      jorb=0
      do j=1,n
        if (is_param_type_orb(j)) then
          jorb=jorb+1
          write(formt3(j),'(a,i2,a,i2,a)') '(orbital',ex_orb_1st_lab(ex_orb_ind(jorb)),&
                                               & '->',ex_orb_2nd_lab(ex_orb_ind(jorb)),')'
        else
          formt3(j)='('//trim(param_type(j))//')'
        endif
      enddo
    elseif (n.eq.param_nb+1) then
      jorb=0
      formt3(1)='(gs)'
      do j=1,n-1
        if (is_param_type_orb(j+1)) then
          jorb=jorb+1
          write(formt3(j+1),'(a,i2,a,i2,a)') '(orbital',ex_orb_1st_lab(ex_orb_ind(jorb)),&
                                                 & '->',ex_orb_2nd_lab(ex_orb_ind(jorb)),')'
        else
          formt3(j+1)='('//trim(param_type(j))//')'
        endif
      enddo
    elseif (n.eq.2*param_nb) then
      jorb=0
      do j=1,n/2
        if (is_param_type_orb(j)) then
          jorb=jorb+1
          write(formt3(j),'(a,i2,a,i2,a)') '(orbital',ex_orb_1st_lab(ex_orb_ind(jorb)),&
                                               & '->',ex_orb_2nd_lab(ex_orb_ind(jorb)),')'
        else
          formt3(j)='('//trim(param_type(j))//')'
        endif
      enddo
      jorb=0
      do j=1,n/2
        if (is_param_type_orb(j)) then
          jorb=jorb+1
          write(formt3(n/2+j),'(a,i2,a,i2,a)') '(orbital',ex_orb_1st_lab(ex_orb_ind(jorb)),&
                                                   & '->',ex_orb_2nd_lab(ex_orb_ind(jorb)),')'
        else
          formt3(n/2+j)='('//trim(param_type(j))//')'
        endif
      enddo
    endif
  endif

  end subroutine get_formt3

!===========================================================================
  subroutine get_formt0(formt0,n,eigval_r,eigvec,sorting,position_in_array)
!---------------------------------------------------------------------------
! Description :
! Description : 
!
! Created     : B. Mussard, October 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! input/output
  character(len=9), allocatable,intent(out) :: formt0(:)
  integer, intent(in)                        :: n
  real(dp), allocatable,intent(in)           :: eigval_r(:)
  real(dp), allocatable,intent(in)           :: eigvec(:,:)
  integer,  allocatable,intent(in)           :: sorting(:)
  integer,  allocatable,intent(in)           :: position_in_array(:)

! local
  integer :: closestpos, largestpos, largestprojpos,i,j
  real(dp):: closest   , largest   , largestproj, contrib

  call alloc('formt0',formt0,n)
  if (n.eq.param_nb+1) then
    closest=1000
    closestpos=0
    largest=0
    largestpos=0
    largestproj=0
    largestprojpos=0
    do j=1,n
      if (abs(eigval_r(j)-etrial).lt.closest) then
        closest=abs(eigval_r(j)-etrial)
        closestpos=j
      endif
      if (abs(eigvec(1, sorting(position_in_array(j)))).gt.largest) then
        largest=abs(eigvec(1, sorting(position_in_array(j))))
        largestpos=j
      endif
      contrib=eigvec(1, sorting(position_in_array(j)))
      do i=1,n-1
        contrib=contrib+eigvec(i+1, sorting(position_in_array(j)))*dpsi_av(i)
      enddo
      if (abs(contrib).gt.largestproj) then
        largestproj=abs(contrib)
        largestprojpos=j
      endif
    enddo
    do j=1,n
      if (j.eq.closestpos)     formt0(j)(1:3)='-->'
      if (j.eq.largestpos)     formt0(j)(4:6)='-->'
      if (j.eq.largestprojpos) formt0(j)(7:9)='-->'
    enddo
  endif

  end subroutine get_formt0

!===========================================================================
  subroutine get_nparmj
!---------------------------------------------------------------------------
! Description : Calculate "nparmj", the number of Jastrow parameters 
! Description : (Copy-Paste of a part of "optimization_menu") 
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! local
  character(len=max_string_len_rout), save :: lhere = 'get_nparmj'
  integer :: nterms4,ia,isp,it,na1,na2,param_i

! begin

  if (use_parser) then
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
      call object_modified ('nparmj')
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
  end subroutine get_nparmj

end module linearresponse_mod
