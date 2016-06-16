module linearresponse_mod

  use all_tools_mod
  use opt_common_mod
  use dets_mod
  use vmc_mod
  use deriv_mod
  use opt_lin_mod

  ! Declaration of global variables and default values

  logical                         :: l_tda=.false.
  logical                         :: do_print=.true.
  real(dp), allocatable           :: amat_av(:,:)
  real(dp), allocatable           :: bmat_av(:,:)
  real(dp), allocatable           :: ovlp_psii_psij_av(:,:)
  real(dp), allocatable           :: linresp_av_eigenval(:)
  real(dp), allocatable           :: linresp_av_eigenval_err(:)

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

    case ('tda')
      call get_next_value (l_tda)

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
  add_diag_mult_exp = 1.d0
  call object_modified ('add_diag_mult_exp')

  !call object_provide ('ncsf')
  !call object_provide ('nparmcsf')
  !call object_provide ('iwcsf')
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
!   routine that wraps the calculation of the number of Jastrow parameters
    write(6,*)
    write(6,'(a,i5)') ' Number of Jastrow parameters:   ', nparmj
  else
    nparmj=0
    call object_modified ('nparmj')
  endif
  write(6,'(a,i5)') ' Number of periodic Jastrow parameters: ', param_pjas_nb

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

  nparm = nparmj+nparmcsf
  call object_modified ('nparm')

! Final messages 
  call object_provide ('nparm')
  call object_provide ('nparmcsf')
  call object_provide ('param_orb_nb')
  call object_provide ('param_exp_nb')
  call object_provide ('param_nb')
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

  call object_average_request('dpsi_av')
  call object_average_request('dpsi_eloc_av')
  call object_average_request('dpsi_dpsi_av')
  call object_average_request('dpsi_dpsi_eloc_av')
  call object_average_request('deloc_av')
  call object_average_request('dpsi_deloc_av')
  if (.not. l_tda) then
    call object_average_request('d2psi_av')
    call object_average_request('d2psi_eloc_av')
  endif

  ! either request for the whole run...
  !call object_error_request('linresp_av_eigenval_err')

  nforce=1
  nwftype=1
  call vmc
  run_done=.true.
  do_print=.true.

  ! ...or ask for it at the end
  !call object_invalidate('amat_av')
  !call object_invalidate('ovlp_psii_psij_av')
  call object_provide('linresp_av_eigenval')

! release
  call object_release ('linresp_av_eigenval',linresp_av_eigenval)
  call object_release ('linresp_av_eigenval_err',linresp_av_eigenval_err)
  call object_release ('amat_av',amat_av)
  call object_release ('bmat_av',bmat_av)
  call object_release ('ovlp_psii_psij_av',ovlp_psii_psij_av)

  end subroutine linearresponse

!===========================================================================
  subroutine  linresp_av_eigenval_bld
!---------------------------------------------------------------------------
! Description : Construct and solve the generalized eigenvalue equation
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
  integer                         :: i,j,info,lwork,temp
  real(dp), allocatable           :: linresp_matrix(:,:)
  real(dp), allocatable           :: ovlp_matrix(:,:)
  real(dp), allocatable           :: work(:)
  real(dp), allocatable           :: eigvec(:,:)
  real(dp), allocatable           :: eigval_r(:)
  real(dp), allocatable           :: eigval_i(:)
  real(dp), allocatable           :: eigval_denom(:)
  integer,  allocatable           :: eigval_srt_ind_to_eigval_ind(:)
  integer,  allocatable           :: eigval_ind_to_eigval_srt_ind(:)

  ! begin
  if (header_exe) then
    call object_create('linresp_av_eigenval')

    call object_needed('param_nb')
    call object_needed('param_pairs')
    call object_needed('diag_stab')

    call object_needed('amat_av')
    call object_needed('ovlp_psii_psij_av')

    return
  endif

  call object_alloc ('linresp_av_eigenval',linresp_av_eigenval,2*param_nb)
  call object_alloc ('linresp_av_eigenval_err',linresp_av_eigenval_err,2*param_nb)

! Construct the ABBA super-matrix from A and B matrices
  call alloc('linresp_matrix',linresp_matrix,2*param_nb,2*param_nb)
  linresp_matrix=0.d0
  if (l_tda) then
    call object_alloc('bmat_av',bmat_av,param_nb,param_nb)
    bmat_av=0.d0
  else
    call object_provide('bmat_av')
  endif
  linresp_matrix(1:param_nb,1:param_nb)=amat_av
  linresp_matrix(param_nb+1:2*param_nb,1:param_nb)=bmat_av
  linresp_matrix(1:param_nb,param_nb+1:2*param_nb)=bmat_av
  linresp_matrix(param_nb+1:2*param_nb,param_nb+1:2*param_nb)=amat_av
  if(l_sym_ham) then
    ham_lin_energy =(ham_lin_energy + transpose(ham_lin_energy))/2.d0
  endif
  do i = 1, param_nb
    if(i > nparmcsf+nparmj .and. i <= nparmcsf+nparmj+param_exp_nb) then
      linresp_matrix(i,i) = linresp_matrix(i,i) + diag_stab * 1.d0
      linresp_matrix(param_nb+i,param_nb+i) = linresp_matrix(param_nb+i,param_nb+i) + diag_stab * 1.d0
    else
      linresp_matrix(i,i) = linresp_matrix(i,i) + diag_stab
      linresp_matrix(i+param_nb,i+param_nb) = linresp_matrix(i+param_nb,i+param_nb) + diag_stab
    endif
  enddo

! Construct the OVERLAP super-matrix from the overlap matrix
  call alloc('ovlp_matrix',ovlp_matrix,2*param_nb,2*param_nb)
  ovlp_matrix=0.d0
  ovlp_matrix(1:param_nb,1:param_nb)=ovlp_psii_psij_av
  ovlp_matrix(param_nb+1:2*param_nb,param_nb+1:2*param_nb)=ovlp_psii_psij_av

! Solve the generalized eigenvalue equation
  if (do_print) then
  write(6,*)
  write(6,'(a,1pd9.1)') 'Solving generalized eigenvalue equation with a_diag =', diag_stab
  endif
! a\ prepare arrays
  call alloc('eigvec',       eigvec,       2*param_nb, 2*param_nb)
  call alloc('eigval_r',     eigval_r,     2*param_nb)
  call alloc('eigval_i',     eigval_i,     2*param_nb)
  call alloc('eigval_denom', eigval_denom, 2*param_nb)

! b\ calculate optimal value of lwork
  lwork = 1
  call alloc('work', work, lwork)
  call dggev('N','V', &
             2*param_nb, linresp_matrix,  &
             2*param_nb, ovlp_matrix,     &
             2*param_nb, eigval_r, eigval_i, eigval_denom, &
             eigvec, 2*param_nb, eigvec, 2*param_nb, &
             work, -1, info)
  if(info /= 0) call die(lhere, 'problem in dggev_alloc: info='+info+' /= 0')
  lwork =  nint(work(1))
  call alloc('work', work, lwork)

! c\ generalized eigenvalue problem
  call dggev('N','V', &
             2*param_nb, linresp_matrix,  &
             2*param_nb, ovlp_matrix,     &
             2*param_nb, eigval_r, eigval_i, eigval_denom, &
             eigvec, 2*param_nb, eigvec, 2*param_nb, &
             work, lwork, info)
  call release ('work', work)
  if(info /= 0) call die(lhere, 'problem in dggev: info='+info+' /= 0 (compare to 2*param_nb='+2*param_nb+')')

! d\ calculate eigenvalue
  do i = 1, 2*param_nb
    eigval_r(i) = eigval_r(i) / eigval_denom(i)
    eigval_i(i) = eigval_i(i) / eigval_denom(i)
  enddo
  call release ('eigval_denom', eigval_denom)

! Sorting out eigenvalues
! ("eigval_srt_ind_to_eigval_ind"
!   is the map from sorted eigenvalues to original eigenvalues)
  call alloc('eigval_srt_ind_to_eigval_ind', eigval_srt_ind_to_eigval_ind, 2*param_nb)
  do i = 1, 2*param_nb
    eigval_srt_ind_to_eigval_ind(i) = i
  enddo
  write(6,'(''eigval_srt_ind_to_eigval_ind='',10000i6)') eigval_srt_ind_to_eigval_ind(1:param_aug_nb)
  do i = 1, 2*param_nb
    do j = i+1, 2*param_nb
      if(eigval_r(eigval_srt_ind_to_eigval_ind(j)) < eigval_r(eigval_srt_ind_to_eigval_ind(i))) then
        temp = eigval_srt_ind_to_eigval_ind(i)
        eigval_srt_ind_to_eigval_ind(i) = eigval_srt_ind_to_eigval_ind(j)
        eigval_srt_ind_to_eigval_ind(j) = temp
      endif
    enddo
  enddo
  write(6,'(''eigval_srt_ind_to_eigval_ind='',10000i6)') eigval_srt_ind_to_eigval_ind(1:param_aug_nb)
! ("eigval_ind_to_eigval_srt_ind"
!   is the map from original eigenvalues to sorted eigenvalues)
  call alloc('eigval_ind_to_eigval_srt_ind', eigval_ind_to_eigval_srt_ind, 2*param_nb)
  do i = 1, 2*param_nb
   eigval_ind_to_eigval_srt_ind(eigval_srt_ind_to_eigval_ind(i)) = i
  enddo
  write(6,'(''eigval_ind_to_eigval_srt_ind='',10000i6)') eigval_ind_to_eigval_srt_ind(1:param_aug_nb)

! Print eigenvalues
  if (do_print) then
  write(6,'(a)') 'Sorted (complex) eigenvalues:'
  do i = 1, 2*param_nb
    write(6,'(a,i5,a,2(f20.6,a))') 'eigenvalue #',i,': ',eigval_r(eigval_srt_ind_to_eigval_ind(i)), ' +', eigval_i(eigval_srt_ind_to_eigval_ind(i)),' i'
  enddo
  endif

! The sorted eigenvalues are "linresp_av_eigenval"
  do i = 1, 2*param_nb
   linresp_av_eigenval (i) = eigval_r(eigval_srt_ind_to_eigval_ind(i))
  enddo

  call release ('eigvec', eigvec)
  call release ('eigval_r', eigval_r)
  call release ('eigval_i', eigval_i)
  call release ('eigval_ind_to_eigval_srt_ind', eigval_ind_to_eigval_srt_ind)
  call release ('eigval_srt_ind_to_eigval_ind', eigval_srt_ind_to_eigval_ind)

  call release('linresp_matrix',linresp_matrix)
  call release('ovlp_matrix',ovlp_matrix)

  end subroutine  linresp_av_eigenval_bld

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

  if (do_print) then
  write(6,*)
  write(6,'(''Adding to ovlp_lin'',es12.4)') add_diag_ovlp
  do i = 1, param_nb
    ovlp_psii_psij_av(i,i) = ovlp_psii_psij_av(i,i)+add_diag_ovlp
  enddo
  call alloc('ovlp_eigvec', ovlp_eigvec, param_nb, param_nb)
  call alloc('ovlp_eigval', ovlp_eigval, param_nb)
  call eigensystem(ovlp_psii_psij_av, ovlp_eigvec, ovlp_eigval, param_nb)
  write(6,*)
  write(6,'(a)') 'Eigenvalues of overlap matrix of current wave function and its first-order derivatives:'
  do i = 1, param_nb
    write(6,'(a,i4,a,es15.8)') 'overlap eigenvalue # ',i,': ',ovlp_eigval(i)
  enddo
  call release('ovlp_eigvec',ovlp_eigvec)
  call release('ovlp_eigval',ovlp_eigval)
  endif

  end subroutine ovlp_psii_psij_av_bld

!===========================================================================
  subroutine  amat_av_bld
!---------------------------------------------------------------------------
! Description : <Psi_i|H|Psi_j>
! Description : (this is inspired by "ham_eigval_av_bld",
! Description :  which does many other things useless in this context)
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
  implicit none

! local
  integer :: i,j,ij

! begin
  if (header_exe) then
    call object_create('amat_av')

    call object_needed('param_nb')
    call object_needed('param_pairs')

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
      amat_av(i,j)=dpsi_dpsi_eloc_av(ij)         &
                 - dpsi_av(j)*dpsi_eloc_av(i)    &
                 - dpsi_av(i)*dpsi_eloc_av(j)    &
                 + dpsi_av(i)*dpsi_av(j)*eloc_av &
                 + dpsi_deloc_covar(i,j)
    enddo
  enddo

  end subroutine amat_av_bld

!===========================================================================
  subroutine  bmat_av_bld
!---------------------------------------------------------------------------
! Description : <Psi_ij|H|Psi_0>
! Description : 
!
! Created     : B. Mussard, June 2016
! Modified    :
!---------------------------------------------------------------------------
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

  end subroutine bmat_av_bld

!===========================================================================
  subroutine  get_nparmj
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
