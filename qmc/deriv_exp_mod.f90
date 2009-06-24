module deriv_exp_mod

  use all_tools_mod
  use basis_mod
  use orbitals_mod
  use determinants_mod
  use csfs_mod
  use electrons_mod
  use psi_mod
  use eloc_mod
  use deriv_orb_mod

! Declaration of global variables and default values
  logical                        :: l_exp_opt_restrict = .true.
  logical                        :: l_deloc_exp_num = .false.
  integer                        :: param_exp_nb = 0
  logical, allocatable           :: is_exp_opt (:)
  logical, allocatable           :: orbital_depends_on_opt_exp (:,:)
  real(dp), allocatable          :: dnorm_basis_dz (:)
  real(dp), allocatable          :: dbasis_ovlp_dz (:,:,:)
  real(dp), allocatable          :: dbasis_ovlp_dz_in_eig_basis (:,:,:)
  real(dp), allocatable          :: dbasis_ovlp_m12_dz (:,:,:)
  real(dp), allocatable          :: dphin_dz (:,:)
  real(dp), allocatable          :: dphin_norm_dz (:,:)
  real(dp), allocatable          :: dphin_ortho_dz (:,:,:)
  real(dp), allocatable          :: grd_dphin_dz (:,:,:)
  real(dp), allocatable          :: grd_dphin_norm_dz (:,:,:)
  real(dp), allocatable          :: grd_dphin_ortho_dz (:,:,:,:)
  real(dp), allocatable          :: lap_dphin_dz (:,:)
  real(dp), allocatable          :: lap_dphin_norm_dz (:,:)
  real(dp), allocatable          :: lap_dphin_ortho_dz (:,:,:)
  real(dp), allocatable          :: dorb_dexp (:,:,:)
  real(dp), allocatable          :: grd_dorb_dexp (:,:,:,:)
  real(dp), allocatable          :: lap_dorb_dexp (:,:,:)

  integer, allocatable                   :: dexp_to_all_bas_nb (:)
  type (type_integer_row), allocatable   :: dexp_to_all_bas (:)
  integer, allocatable                   :: dexp_to_bas_nb (:)
  type (type_integer_row), allocatable   :: dexp_to_bas (:)
  integer, allocatable                   :: bas_to_dexp (:)

  real(dp), allocatable          :: ddet_dexp_unq_up (:,:)
  real(dp), allocatable          :: ddet_dexp_unq_dn (:,:)
  real(dp), allocatable          :: ddet_dexp_col_unq_up (:,:,:)
  real(dp), allocatable          :: ddet_dexp_col_unq_dn (:,:,:)
  real(dp), allocatable          :: dpsid_exp (:)
  real(dp), allocatable          :: dpsi_exp (:)
  real(dp), allocatable          :: dpsi_lnexp (:)
  real(dp), allocatable          :: grd_ddet_dexp_unq_up (:,:,:,:)
  real(dp), allocatable          :: grd_ddet_dexp_unq_dn (:,:,:,:)
  real(dp), allocatable          :: lap_ddet_dexp_unq_up (:,:,:)
  real(dp), allocatable          :: lap_ddet_dexp_unq_dn (:,:,:)
  real(dp), allocatable          :: grd_dpsid_exp_over_dpsid_exp (:,:,:)
  real(dp), allocatable          :: lap_dpsid_exp_over_dpsid_exp (:,:)
  real(dp), allocatable          :: lap_ln_dpsid_exp (:,:)
  real(dp), allocatable          :: sum_lap_ln_dpsid_exp (:)
  real(dp), allocatable          :: grd_dpsi_exp_over_dpsi_exp (:,:,:)
  real(dp), allocatable          :: sum_lap_ln_dpsi_exp (:)
  real(dp), allocatable          :: eloc_kin_exp (:)
  real(dp), allocatable          :: eloc_exp (:)
  real(dp), allocatable          :: deloc_exp (:)
  real(dp), allocatable          :: deloc_exp_num (:)
  real(dp), allocatable          :: deloc_exp1 (:)
  real(dp), allocatable          :: deloc_exp2 (:)
  real(dp), allocatable          :: deloc_lnexp (:)

  real(dp), allocatable          :: slater_mat_exp_trans_inv_up (:,:,:,:,:)
  real(dp), allocatable          :: slater_mat_exp_trans_inv_dn (:,:,:,:,:)

  contains

! ==============================================================================
  subroutine param_exp_nb_bld
! ------------------------------------------------------------------------------
! Description   : number of derivatives wrt to optimized exponents
! Description   : dexp_to_all_bas is the correspondence between exponent parameters
! Description   :                 and all basis exponents to be updated including those 
! Description   :                 not directly involved in the optimization but updated
! Description   :                 because of symmetry constraints
! Description   :                 (for instance, if there is no pz occupied orbitals
! Description   :                  only px and py basis functions will be involved in the optimization,
! Description   :                  the derivative wrt the pz exponent will be zero,
! Description   :                  but pz exponent will also be updated to keep the same exponents for px,py,pz)
! Description   :           
! Description   : dexp_to_bas     is the correspondence between exponent parameters
! Description   :                 and only the basis exponents directly involved in the optimization
!
! Created       : J. Toulouse, 25 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_j, dexp_i, exp_opt_lab_i !, exp_opt_lab_j
  logical, allocatable :: is_basis_func_attributed (:)

! header
  if (header_exe) then

   call object_create ('param_exp_nb', param_exp_nb_index)
   call object_create ('dexp_to_all_bas_nb')
   call object_create ('dexp_to_all_bas')
   call object_create ('dexp_to_bas_nb')
   call object_create ('dexp_to_bas')
   call object_create ('bas_to_dexp')

   call object_needed ('nbasis')
   call object_needed ('n_bas')
   call object_needed ('l_bas')
   call object_needed ('m_bas')
   call object_needed ('ictype_basis')
   call object_needed ('exp_opt_lab_nb')
   call object_needed ('exp_opt_lab')
   call object_needed ('is_exp_opt')
   call object_needed ('basis_fns_name')

   return

  endif

! begin

! allocation
  call alloc ('bas_to_dexp', bas_to_dexp, nbasis)
  call alloc ('is_basis_func_attributed', is_basis_func_attributed, nbasis)

  bas_to_dexp (:) = 0
  is_basis_func_attributed (:)= .false.

! provide zex but do not record dependency
  call object_provide ('zex')
  call object_provide ('coef')

! to improve
  param_exp_nb = 0

  do exp_opt_lab_i = 1, exp_opt_lab_nb

    bas_i = exp_opt_lab (exp_opt_lab_i)

    if (is_basis_func_attributed (bas_i)) cycle

    param_exp_nb = param_exp_nb + 1
    call object_alloc ('dexp_to_all_bas_nb', dexp_to_all_bas_nb, param_exp_nb)
    call object_alloc ('dexp_to_all_bas', dexp_to_all_bas, param_exp_nb)
    call object_alloc ('dexp_to_bas_nb', dexp_to_bas_nb, param_exp_nb)
    call object_alloc ('dexp_to_bas', dexp_to_bas, param_exp_nb)
    dexp_to_all_bas_nb (param_exp_nb) = 0
    dexp_to_bas_nb (param_exp_nb) = 0

!    do exp_opt_lab_j = exp_opt_lab_i, exp_opt_lab_nb
!      bas_j = exp_opt_lab (exp_opt_lab_j)
     do bas_j = 1, nbasis

      if (bas_j == bas_i) then
         dexp_to_all_bas_nb (param_exp_nb) = dexp_to_all_bas_nb (param_exp_nb) + 1
         call append (dexp_to_all_bas (param_exp_nb)%row, bas_j)
         dexp_to_bas_nb (param_exp_nb) = dexp_to_bas_nb (param_exp_nb) + 1
         call append (dexp_to_bas (param_exp_nb)%row, bas_j)
         bas_to_dexp (bas_j) = param_exp_nb
         is_basis_func_attributed (bas_j) = .true.
         cycle
      endif

      if (l_exp_opt_restrict) then

!       restriction on exponent parameters
        if (zex (bas_i, iwf) == zex (bas_j, iwf) .and. ictype_basis(bas_i) == ictype_basis(bas_j) .and. abs(n_bas(bas_i)) == abs(n_bas(bas_j)) .and. abs(l_bas(bas_i)) == abs(l_bas(bas_j))) then

!         exponent parameter involved in the optimization
          if (is_exp_opt (bas_j)) then
             dexp_to_all_bas_nb (param_exp_nb) = dexp_to_all_bas_nb (param_exp_nb) + 1
             call append (dexp_to_all_bas (param_exp_nb)%row, bas_j)
             dexp_to_bas_nb (param_exp_nb) = dexp_to_bas_nb (param_exp_nb) + 1
             call append (dexp_to_bas (param_exp_nb)%row, bas_j)
             bas_to_dexp (bas_j) = param_exp_nb
             is_basis_func_attributed (bas_j) = .true.
             cycle
          
!         exponent parameter not directly involved in the optimization
          else
             dexp_to_all_bas_nb (param_exp_nb) = dexp_to_all_bas_nb (param_exp_nb) + 1
             call append (dexp_to_all_bas (param_exp_nb)%row, bas_j)
             cycle
          endif
      
        endif ! restriction on exponent parameters

      endif ! l_exp_opt_restrict

    enddo ! lab_j

  enddo ! exp_opt_lab_i

  write(6,'(a,i8)') ' Number of exponent parameters =',param_exp_nb
  do dexp_i = 1, param_exp_nb
!    write(6,'(a,i3,a,100i3)') ' Exponent parameter # ',dexp_i,' corresponds to exponents: ', dexp_to_bas(dexp_i)%row (:)
!    write(6,'(a,i3,a,100(a,x))') ' Exponent parameter # ',dexp_i,' corresponds to basis functions: ', (trim(basis_fns_name (dexp_to_bas(dexp_i)%row (bas_i))), bas_i=1, dexp_to_bas_nb(dexp_i))
    write(6,'(a,i3,a,100(a,x))') ' Exponent parameter # ',dexp_i,' corresponds to basis functions: ', (trim(basis_fns_name (dexp_to_all_bas(dexp_i)%row (bas_i))), bas_i=1, dexp_to_all_bas_nb(dexp_i))
  enddo ! dexp_i

  end subroutine param_exp_nb_bld

! ==============================================================================
  subroutine exp_opt_lab_bld
! ------------------------------------------------------------------------------
! Description   : build list of exponents to be optimized if not read
!
! Created       : J. Toulouse, 18 Apr 2008
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, orb_i
  character(len=max_string_len_rout), save :: lhere = 'exp_opt_lab_bld' !fp

! header
  if (header_exe) then

   call object_create ('exp_opt_lab_nb')
   call object_create ('exp_opt_lab')
   call object_create ('is_exp_opt')

   call object_needed ('nbasis')
   call object_needed ('orb_tot_nb')
!   call object_needed ('orb_occ_last_in_wf_lab')
   call object_needed ('orb_occ_in_wf_lab')

   return

  endif

! begin
  call object_alloc ('is_exp_opt', is_exp_opt, nbasis)
  is_exp_opt (:) = .false.

! provide orbital coef
  call object_provide ('coef')

!  do bas_i = 1, nbasis
!   do orb_i = 1, orb_occ_last_in_wf_lab
!     if (orb_occ_in_wf (orb_i) .and. coef (bas_i, orb_i, 1) /= 0.d0) then
!        exp_opt_lab_nb = exp_opt_lab_nb + 1
!        call object_alloc ('exp_opt_lab', exp_opt_lab, exp_opt_lab_nb)
!        exp_opt_lab (exp_opt_lab_nb) = bas_i
!        exit
!      endif
!   enddo ! orb_i
!  enddo ! bas_i

  if (l_opt_orb) then
    call object_provide ('orb_mix_lab')
  endif

  do bas_i = 1, nbasis
   do orb_i = 1, orb_tot_nb

!    occupied orbital with non zero coefficient
     if (orb_occ_in_wf (orb_i) .and. coef (bas_i, orb_i, 1) /= 0.d0) then
        exp_opt_lab_nb = exp_opt_lab_nb + 1
        call object_alloc ('exp_opt_lab', exp_opt_lab, exp_opt_lab_nb)
        exp_opt_lab (exp_opt_lab_nb) = bas_i
        is_exp_opt (bas_i) = .true.
        exit
      endif

!    virtual orbital with non zero coefficient which can be mixed in in orbital optimization
     if (l_opt_orb) then
      if (orb_mix_lab (orb_i) .and. coef (bas_i, orb_i, 1) /= 0.d0) then
        exp_opt_lab_nb = exp_opt_lab_nb + 1
        call object_alloc ('exp_opt_lab', exp_opt_lab, exp_opt_lab_nb)
        exp_opt_lab (exp_opt_lab_nb) = bas_i
        is_exp_opt (bas_i) = .true.
        exit
       endif
      endif

   enddo ! orb_i
  enddo ! bas_i

! check
  call require (lhere, 'exp_opt_lab_nb <= nbasis', exp_opt_lab_nb <= nbasis) !fp

  end subroutine exp_opt_lab_bld

! ==============================================================================
  subroutine orbital_depends_on_opt_exp_bld
! ------------------------------------------------------------------------------
! Description   : orbital_depends_on_opt_exp = true if the orbital
! Description   : depends on the exponents
!
! Created       : J. Toulouse, 19 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer orb_i, dexp_i, dexp_to_bas_i, bas_i

! header
  if (header_exe) then

   call object_create ('orbital_depends_on_opt_exp')

   call object_needed ('orb_occ_last_in_wf_lab')
   call object_needed ('param_exp_nb')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')
   call object_needed ('coef')
   call object_needed ('iwf')

   return

  endif

! begin

! allocation
  call object_alloc ('orbital_depends_on_opt_exp', orbital_depends_on_opt_exp, orb_occ_last_in_wf_lab, param_exp_nb)
  orbital_depends_on_opt_exp (:, :) = .false.

! warning: coef changes during orbital optimization
  do orb_i = 1, orb_occ_last_in_wf_lab
    do dexp_i = 1, param_exp_nb
      do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
        bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
        if (coef (bas_i, orb_i, iwf) /= 0.d0) then
          orbital_depends_on_opt_exp (orb_i, dexp_i) = .true.
          exit
        endif
      enddo ! dexp_to_bas_i
    enddo ! dexp_i
  enddo ! orb_i

  end subroutine orbital_depends_on_opt_exp_bld

! ==============================================================================
  subroutine dnorm_basis_dz_bld
! ------------------------------------------------------------------------------
! Description   : Derivatives of the norm of the basis functions wrt to exponents
! Warning       : implemented only for Slater functions
!
! Created       : J. Toulouse, 17 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, n

! header
  if (header_exe) then

   call object_create ('dnorm_basis_dz')

   call object_needed ('nbasis')
   call object_needed ('norm_basis')
   call object_needed ('zex')
   call object_needed ('n_bas')
   call object_needed ('iwf')

   return

  endif

! begin

! requirements
  if (numr > 0) then
   call die (here, 'numr='+numr+'>0. Implemented only for analytical functions (numr <= 0)')
  endif

! allocation
  call object_alloc ('dnorm_basis_dz', dnorm_basis_dz, nbasis)

  do bas_i = 1, nbasis
    if (zex (bas_i, iwf) == 0.d0) then
      write(6,'(2a,i3,a,i3,a,es15.8)') trim(here), ': zex(',bas_i,',',iwf,')=', zex (bas_i, iwf)
      call die (here, 'exponent must be non zero')
    endif
    n = n_bas (bas_i)
    if (n > 0) then
      dnorm_basis_dz (bas_i) = (n + 0.5d0) * norm_basis (bas_i) / zex (bas_i, iwf)
    else
      call die (here, 'implemented only for Slater functions (n > 0)')
    endif
  enddo

  end subroutine dnorm_basis_dz_bld

! ==============================================================================
  subroutine dbasis_ovlp_dz_bld
! ------------------------------------------------------------------------------
! Description   : derivative of basis overlap matrix wrt exponents
!
! Created       : J. Toulouse, 29 Mai 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_j, dexp_i, dexp_j
  integer n_i, n_j, l_i, l_j, m_i, m_j
  real(dp) exp_i, exp_j

! header
  if (header_exe) then

   call object_create ('dbasis_ovlp_dz')

   call object_needed ('param_exp_nb')
   call object_needed ('nbasis')
   call object_needed ('n_bas')
   call object_needed ('l_bas')
   call object_needed ('m_bas')
   call object_needed ('zex')
   call object_needed ('bas_to_dexp')

   return

  endif

! begin

! requirements
  if (numr /= 0 .and. numr /= -1 .and. numr /= -2 .and. numr /= -3) then
   write(6,*) trim(here),': numr=',numr,' /= 0, - 1, -2, -3'
   write(6,*) trim(here),': implemented only for the case numr = -1 or -2 or -3'
   write(6,*) trim(here),': i.e. analytical slater functions with order: 1s, 2s, 3s, ...., 2p, 3p...'
   call die (here)
  endif
  if (ncent /= 1 ) then
   call die (here, 'ncent='+ncent+' /= 1. Implemented only for one center')
  endif
! end requirements

! allocation
  call object_alloc ('dbasis_ovlp_dz', dbasis_ovlp_dz, nbasis, nbasis, param_exp_nb)

  dbasis_ovlp_dz (:,:,:) = 0.d0

  do bas_i = 1, nbasis

      n_i   = abs(n_bas (bas_i))
      l_i   = l_bas (bas_i)
      m_i   = m_bas (bas_i)
      exp_i = zex (bas_i, 1)
      dexp_i = bas_to_dexp (bas_i)

    do bas_j = 1, nbasis

      n_j   = abs(n_bas (bas_j))
      l_j   = l_bas (bas_j)
      m_j   = m_bas (bas_j)
      exp_j = zex (bas_j, 1)
      dexp_j = bas_to_dexp (bas_j)

      if (dexp_i /= 0) then
       dbasis_ovlp_dz (bas_i, bas_j, dexp_i) = dbasis_ovlp_dz (bas_i, bas_j, dexp_i) - slater_ovlp (n_i+1, l_i, m_i, exp_i, n_j, l_j, m_j, exp_j)
      endif

      if (dexp_j /= 0) then
       dbasis_ovlp_dz (bas_i, bas_j, dexp_j) = dbasis_ovlp_dz (bas_i, bas_j, dexp_j) - slater_ovlp (n_i, l_i, m_i, exp_i, n_j+1, l_j, m_j, exp_j)
      endif

    enddo ! bas_j
  enddo ! bas_i

  end subroutine dbasis_ovlp_dz_bld

! ==============================================================================
  subroutine dbasis_ovlp_dz_in_eig_basis_bld
! ------------------------------------------------------------------------------
! Description   : derivative of overlap matrix wrt exponents
! Description   : in eigenbasis of overlap matrix
!
! Created       : J. Toulouse, 30 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_j, bas_k, bas_l, dexp_i

! header
  if (header_exe) then

   call object_create ('dbasis_ovlp_dz_in_eig_basis')

   call object_needed ('param_exp_nb')
   call object_needed ('nbasis')
   call object_needed ('basis_ovlp_eigvec')
   call object_needed ('dbasis_ovlp_dz')

   return

  endif

! begin

! allocation
  call object_alloc ('dbasis_ovlp_dz_in_eig_basis', dbasis_ovlp_dz_in_eig_basis, nbasis, nbasis, param_exp_nb)

  do dexp_i = 1, param_exp_nb
   do bas_i = 1, nbasis
    do bas_j = 1, nbasis
     dbasis_ovlp_dz_in_eig_basis (bas_i, bas_j, dexp_i) = 0.d0
     do bas_k = 1, nbasis
       do bas_l = 1, nbasis
        dbasis_ovlp_dz_in_eig_basis (bas_i, bas_j, dexp_i) = dbasis_ovlp_dz_in_eig_basis (bas_i, bas_j, dexp_i) + basis_ovlp_eigvec (bas_k, bas_i) * dbasis_ovlp_dz (bas_k, bas_l, dexp_i) * basis_ovlp_eigvec (bas_l, bas_j)
       enddo ! bas_l
     enddo ! bas_k
    enddo ! bas_j
   enddo ! bas_i
  enddo ! dexp_i

  end subroutine dbasis_ovlp_dz_in_eig_basis_bld

! ==============================================================================
  subroutine dbasis_ovlp_m12_dz_bld
! ------------------------------------------------------------------------------
! Description   : derivative of overlap matrix to the power -1/2 wrt exponents
! Description   : according to J. Jorgensen and J. Simons, J. Chem. Phys. 79, 334 (1983)
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_j, bas_k, bas_l, dexp_i
  real(dp) :: dbasis_ovlp_m12_dz_in_eig_basis

! header
  if (header_exe) then

   call object_create ('dbasis_ovlp_m12_dz')

   call object_needed ('param_exp_nb')
   call object_needed ('nbasis')
   call object_needed ('basis_ovlp_eigvec')
   call object_needed ('basis_ovlp_eigval')
   call object_needed ('dbasis_ovlp_dz_in_eig_basis')

   return

  endif

! begin

! allocation
  call object_alloc ('dbasis_ovlp_m12_dz', dbasis_ovlp_m12_dz, nbasis, nbasis, param_exp_nb)

  do dexp_i = 1, param_exp_nb
  do bas_i = 1, nbasis
   do bas_j = 1, nbasis
     dbasis_ovlp_m12_dz (bas_i, bas_j, dexp_i) = 0.d0
     do bas_k = 1, nbasis
       do bas_l = 1, nbasis
       dbasis_ovlp_m12_dz_in_eig_basis = - dbasis_ovlp_dz_in_eig_basis (bas_k, bas_l, dexp_i) /(dsqrt(basis_ovlp_eigval (bas_k)) * dsqrt(basis_ovlp_eigval (bas_l)) * (dsqrt(basis_ovlp_eigval (bas_k)) + dsqrt(basis_ovlp_eigval (bas_l))))
       dbasis_ovlp_m12_dz (bas_i, bas_j, dexp_i) = dbasis_ovlp_m12_dz (bas_i, bas_j, dexp_i) + basis_ovlp_eigvec (bas_i, bas_k) * dbasis_ovlp_m12_dz_in_eig_basis * basis_ovlp_eigvec (bas_j, bas_l)
       enddo ! bas_l
     enddo ! bas_k
   enddo ! bas_j
  enddo ! bas_i
  enddo ! dexp_i

  end subroutine dbasis_ovlp_m12_dz_bld

! ==============================================================================
  subroutine dphin_dz_bld
! ------------------------------------------------------------------------------
! Description   : Derivatives of the unnormalized basis functions wrt to exponents
!
! Created       : J. Toulouse, 25 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, cent_i, dexp_i, dexp_to_bas_i

! header
  if (header_exe) then

   call object_create ('dphin_dz', dphin_dz_index)

   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('basis_fns_cent')
   call object_needed ('r_en')
   call object_needed ('phin')
   call object_needed ('param_exp_nb')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')

   return

  endif

! begin

! requirements
  if (numr > 0) then
   call die (here, 'numr='+numr+'>0. Implemented only for analytical functions (numr <= 0)')
  endif

! allocation
  call object_alloc ('dphin_dz', dphin_dz, nelec, nbasis)

  do dexp_i = 1, param_exp_nb
    do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
      bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
      cent_i = basis_fns_cent (bas_i)
      dphin_dz (1:nelec, bas_i) = -r_en (1:nelec, cent_i) * phin (bas_i, 1:nelec)
    enddo ! dexp_to_bas_i
  enddo ! dexp_i

  end subroutine dphin_dz_bld

! ==============================================================================
  subroutine dphin_norm_dz_bld
! ------------------------------------------------------------------------------
! Description   : Derivatives of normalized basis functions wrt exponents
!
! Created       : J. Toulouse, 17 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, dexp_i, dexp_to_bas_i

! header
  if (header_exe) then

   call object_create ('dphin_norm_dz', dphin_norm_dz_index)

   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('phin')
   call object_needed ('dphin_dz')
   call object_needed ('norm_basis')
   call object_needed ('dnorm_basis_dz')
   call object_needed ('param_exp_nb')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')

   return

  endif

! begin

! allocation
  call object_alloc ('dphin_norm_dz', dphin_norm_dz, nelec, nbasis)

  do dexp_i = 1, param_exp_nb
    do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
      bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
      dphin_norm_dz (1:nelec, bas_i) = norm_basis (bas_i) * dphin_dz (1:nelec, bas_i) + dnorm_basis_dz (bas_i) * phin (bas_i, 1:nelec)
    enddo ! dexp_to_bas_i
  enddo ! dexp_i

  end subroutine dphin_norm_dz_bld

! ==============================================================================
  subroutine dphin_ortho_dz_bld
! ------------------------------------------------------------------------------
! Description   : Derivatives of orthonomalized basis functions wrt exponents
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_k, dexp_i

! header
  if (header_exe) then

   call object_create ('dphin_ortho_dz', dphin_ortho_dz_index)

   call object_needed ('param_exp_nb')
   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('phin')
   call object_needed ('dphin_dz')
   call object_needed ('bas_to_dexp')
   call object_needed ('basis_ovlp_m12')
   call object_needed ('dbasis_ovlp_m12_dz')

   return

  endif

! begin

! allocation
  call object_alloc ('dphin_ortho_dz', dphin_ortho_dz, nelec, nbasis, param_exp_nb)

  do dexp_i = 1, param_exp_nb
     do bas_i = 1, nbasis
      dphin_ortho_dz (1:nelec, bas_i, dexp_i) = 0.d0
        do bas_k = 1, nbasis
         if (bas_to_dexp (bas_k) == dexp_i) then
          dphin_ortho_dz (1:nelec, bas_i, dexp_i) = dphin_ortho_dz (1:nelec, bas_i, dexp_i) + basis_ovlp_m12 (bas_i, bas_k) * dphin_dz (1:nelec, bas_k)
         endif
         dphin_ortho_dz (1:nelec, bas_i, dexp_i) = dphin_ortho_dz (1:nelec, bas_i, dexp_i) + dbasis_ovlp_m12_dz (bas_i, bas_k, dexp_i) * phin (bas_k, 1:nelec)
        enddo ! bas_k
    enddo ! bas_i
  enddo ! dexp_i

  end subroutine dphin_ortho_dz_bld

! ==============================================================================
  subroutine grd_dphin_dz_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of derivatives of the unnormalized basis functions wrt to exponents
!
! Created       : J. Toulouse, 27 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer elec_i, bas_i, dim_i, cent_i, dexp_i, dexp_to_bas_i

! header
  if (header_exe) then

   call object_create ('grd_dphin_dz', grd_dphin_dz_index)

   call object_needed ('nbasis')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('basis_fns_cent')
   call object_needed ('r_en')
   call object_needed ('grd_dist_en')
   call object_needed ('phin')
   call object_needed ('dphin')
   call object_needed ('param_exp_nb')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')

   return

  endif

! begin

! requirements
  if (numr > 0) then
   call die (here, 'numr='+numr+'>0. Implemented only for analytical functions (numr <= 0)')
  endif

! allocation
  call object_alloc ('grd_dphin_dz', grd_dphin_dz, ndim, nelec, nbasis)

  do dexp_i = 1, param_exp_nb
    do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
      bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
      cent_i = basis_fns_cent (bas_i)
      do elec_i = 1, nelec
        do dim_i = 1, ndim
          grd_dphin_dz (dim_i, elec_i, bas_i) = - grd_dist_en (dim_i, elec_i, cent_i) * phin (bas_i, elec_i) &
                                              - r_en (elec_i, cent_i) * dphin (dim_i, bas_i, elec_i)
        enddo ! dim_i
      enddo ! elec_i
    enddo ! dexp_to_bas_i
  enddo ! dexp_i

  end subroutine grd_dphin_dz_bld

! ==============================================================================
  subroutine grd_dphin_norm_dz_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of derivatives of normalized basis functions wrt exponents
!
! Created       : J. Toulouse, 17 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer elec_i, bas_i, dim_i, dexp_i, dexp_to_bas_i

! header
  if (header_exe) then

   call object_create ('grd_dphin_norm_dz', grd_dphin_norm_dz_index)

   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('norm_basis')
   call object_needed ('grd_dphin_dz')
   call object_needed ('dnorm_basis_dz')
   call object_needed ('dphin')
   call object_needed ('param_exp_nb')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')

   return

  endif

! begin

! allocation
  call object_alloc ('grd_dphin_norm_dz', grd_dphin_norm_dz, ndim, nelec, nbasis)

  do dexp_i = 1, param_exp_nb
    do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
      bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
      do elec_i = 1, nelec
        do dim_i = 1, ndim
          grd_dphin_norm_dz (dim_i, elec_i, bas_i) = norm_basis (bas_i) * grd_dphin_dz (dim_i, elec_i, bas_i) + dnorm_basis_dz (bas_i) * dphin (dim_i, bas_i, elec_i)
        enddo ! dim_i
      enddo ! elec_i
    enddo ! dexp_to_bas_i
  enddo ! dexp_i

  end subroutine grd_dphin_norm_dz_bld

! ==============================================================================
  subroutine grd_dphin_ortho_dz_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of derivatives of orthonormalized basis functions wrt exponents
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer elec_i, bas_i, bas_k, dim_i, dexp_i

! header
  if (header_exe) then

   call object_create ('grd_dphin_ortho_dz', grd_dphin_ortho_dz_index)

   call object_needed ('param_exp_nb')
   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('grd_dphin_dz')
   call object_needed ('dphin')
   call object_needed ('bas_to_dexp')
   call object_needed ('basis_ovlp_m12')
   call object_needed ('dbasis_ovlp_m12_dz')

   return

  endif

! begin

! allocation
  call object_alloc ('grd_dphin_ortho_dz', grd_dphin_ortho_dz, ndim, nelec, nbasis, param_exp_nb)

  do dexp_i = 1, param_exp_nb
     do bas_i = 1, nbasis
      do elec_i = 1, nelec
        do dim_i = 1, ndim
          grd_dphin_ortho_dz (dim_i, elec_i, bas_i, dexp_i) = 0.d0
          do bas_k = 1, nbasis
           if (bas_to_dexp (bas_k) == dexp_i) then
           grd_dphin_ortho_dz (dim_i, elec_i, bas_i, dexp_i) =  grd_dphin_ortho_dz (dim_i, elec_i, bas_i, dexp_i) + basis_ovlp_m12 (bas_i, bas_k) * grd_dphin_dz (dim_i, elec_i, bas_k)
           endif
           grd_dphin_ortho_dz (dim_i, elec_i, bas_i, dexp_i) =  grd_dphin_ortho_dz (dim_i, elec_i, bas_i, dexp_i) + dbasis_ovlp_m12_dz (bas_i, bas_k, dexp_i) * dphin (dim_i, bas_k, elec_i)
          enddo ! bas_k
        enddo ! dim_i
      enddo ! elec_i
    enddo ! bas_i
  enddo ! dexp_i

  end subroutine grd_dphin_ortho_dz_bld

! ==============================================================================
  subroutine lap_dphin_dz_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of derivatives of unnormalized basis functions wrt to exponents
!
! Created       : J. Toulouse, 27 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer elec_i, bas_i, cent_i, dim_i, dexp_i, dexp_to_bas_i
  real(dp) dotproduct

! header
  if (header_exe) then

   call object_create ('lap_dphin_dz', lap_dphin_dz_index)

   call object_needed ('nbasis')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('basis_fns_cent')
   call object_needed ('r_en')
   call object_needed ('grd_dist_en')
   call object_needed ('lap_dist_en')
   call object_needed ('grd_dphin_dz')
   call object_needed ('phin')
   call object_needed ('dphin')
   call object_needed ('d2phin')
   call object_needed ('param_exp_nb')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')

   return

  endif

! begin

! requirements
  if (numr > 0) then
   call die (here, 'numr='+numr+'>0. Implemented only for analytical functions (numr <= 0)')
  endif

! allocation
  call object_alloc ('lap_dphin_dz', lap_dphin_dz, nelec, nbasis)

  do dexp_i = 1, param_exp_nb
    do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
      bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
      cent_i = basis_fns_cent (bas_i)

      do elec_i = 1, nelec

      dotproduct = 0.d0
      do dim_i = 1, ndim
        dotproduct = dotproduct + grd_dist_en (dim_i, elec_i, cent_i) * dphin (dim_i, bas_i, elec_i)
      enddo ! dim_i

      lap_dphin_dz (elec_i, bas_i) = - lap_dist_en (elec_i, cent_i) * phin (bas_i, elec_i)          &
                                     - 2.d0 * dotproduct                                            &
                                     - r_en (elec_i, cent_i) * d2phin (bas_i, elec_i)
     enddo ! elec_i

    enddo ! dexp_to_bas_i
  enddo ! dexp_i

  end subroutine lap_dphin_dz_bld

! ==============================================================================
  subroutine lap_dphin_norm_dz_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of derivatives of normalized basis functions wrt exponents
!
! Created       : J. Toulouse, 17 Apr 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, elec_i, dexp_i, dexp_to_bas_i

! header
  if (header_exe) then

   call object_create ('lap_dphin_norm_dz', lap_dphin_norm_dz_index)

   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('norm_basis')
   call object_needed ('lap_dphin_dz')
   call object_needed ('dnorm_basis_dz')
   call object_needed ('d2phin')
   call object_needed ('param_exp_nb')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')

   return

  endif

! begin

! allocation
  call object_alloc ('lap_dphin_norm_dz', lap_dphin_norm_dz, nelec, nbasis)

  do dexp_i = 1, param_exp_nb
    do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
      bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
      do elec_i = 1, nelec
        lap_dphin_norm_dz (elec_i, bas_i) = norm_basis (bas_i) * lap_dphin_dz (elec_i, bas_i) + dnorm_basis_dz (bas_i) * d2phin (bas_i, elec_i)
      enddo ! elec_i
    enddo ! dexp_to_bas_i
  enddo ! dexp_i

  end subroutine lap_dphin_norm_dz_bld

! ==============================================================================
  subroutine lap_dphin_ortho_dz_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of derivatives of orthonormalized basis functions wrt exponents
!
! Created       : J. Toulouse, 29 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer bas_i, bas_k, elec_i, dexp_i

! header
  if (header_exe) then

   call object_create ('lap_dphin_ortho_dz', lap_dphin_ortho_dz_index)

   call object_needed ('nbasis')
   call object_needed ('nelec')
   call object_needed ('lap_dphin_dz')
   call object_needed ('d2phin')
   call object_needed ('param_exp_nb')
   call object_needed ('bas_to_dexp')
   call object_needed ('basis_ovlp_m12')
   call object_needed ('dbasis_ovlp_m12_dz')

   return

  endif

! begin

! allocation
  call object_alloc ('lap_dphin_ortho_dz', lap_dphin_ortho_dz, nelec, nbasis, param_exp_nb)

  do dexp_i = 1, param_exp_nb
    do bas_i = 1, nbasis
      do elec_i = 1, nelec
        lap_dphin_ortho_dz (elec_i, bas_i, dexp_i) = 0.d0
        do bas_k = 1, nbasis
        if ((bas_to_dexp (bas_k) == dexp_i)) then
         lap_dphin_ortho_dz (elec_i, bas_i, dexp_i) = lap_dphin_ortho_dz (elec_i, bas_i, dexp_i) + basis_ovlp_m12 (bas_i, bas_k) * lap_dphin_dz (elec_i, bas_k)
        endif
         lap_dphin_ortho_dz (elec_i, bas_i, dexp_i) = lap_dphin_ortho_dz (elec_i, bas_i, dexp_i) + dbasis_ovlp_m12_dz (bas_i, bas_k, dexp_i) * d2phin (bas_k, elec_i)
        enddo ! bas_k
      enddo ! elec_i
    enddo ! bas_i
  enddo ! dexp_i

  end subroutine lap_dphin_ortho_dz_bld

! ===============================================================================
  subroutine dorb_dexp_bld
! -------------------------------------------------------------------------------
! Description   : Derivatives of occupied orbitals wrt to optimized basis exponents
!
! Created       : J. Toulouse, 25 Jan 2007
! -------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i, bas_i, orb_i, dexp_to_bas_i

! header
  if (header_exe) then

   call object_create ('dorb_dexp')

   call object_needed ('param_exp_nb')
   call object_needed ('nbasis')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')
   call object_needed ('orb_occ_last_in_wf_lab')
   call object_needed ('nelec')
   call object_needed ('orbital_depends_on_opt_exp')

   return

  endif

! begin

! allocation
  call object_alloc ('dorb_dexp', dorb_dexp, nelec, orb_occ_last_in_wf_lab, param_exp_nb)
  dorb_dexp (:,:,:) = 0.d0

  select case (trim(basis_functions_varied))
  case ('unnormalized')
  call object_provide_by_index (dorb_dexp_bld_index, dphin_dz_index)
  call object_provide_by_index (dorb_dexp_bld_index, coef_index)
  do dexp_i = 1, param_exp_nb
   do orb_i = 1, orb_occ_last_in_wf_lab
    if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
     do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
       bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
       dorb_dexp (:, orb_i, dexp_i) = dorb_dexp (:, orb_i, dexp_i) + coef (bas_i, orb_i, iwf) * dphin_dz (:, bas_i)
     enddo ! dexp_to_bas_i
   enddo ! orb_i
  enddo ! dexp_i

  case ('normalized')
  call object_provide_by_index (dorb_dexp_bld_index, dphin_norm_dz_index)
  call object_provide_by_index (dorb_dexp_bld_index, coef_orb_on_norm_basis_index)
  do dexp_i = 1, param_exp_nb
   do orb_i = 1, orb_occ_last_in_wf_lab
    if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
     do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
       bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
       dorb_dexp (:, orb_i, dexp_i) = dorb_dexp (:, orb_i, dexp_i) + coef_orb_on_norm_basis (bas_i, orb_i, iwf) * dphin_norm_dz (:, bas_i)
     enddo ! dexp_to_bas_i
   enddo ! orb_i
  enddo ! dexp_i

  case ('orthonormalized')
  call object_provide_by_index (dorb_dexp_bld_index, dphin_ortho_dz_index)
  call object_provide_by_index (dorb_dexp_bld_index, coef_orb_on_ortho_basis_index)
  do dexp_i = 1, param_exp_nb
   do orb_i = 1, orb_occ_last_in_wf_lab
    if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
     do bas_i = 1, nbasis
       dorb_dexp (:, orb_i, dexp_i) = dorb_dexp (:, orb_i, dexp_i) + coef_orb_on_ortho_basis (bas_i, orb_i, iwf) * dphin_ortho_dz (:, bas_i, dexp_i)
     enddo ! bas_i
   enddo ! orb_i
  enddo ! dexp_i

  case default
   call die (here, 'unknown case >'+trim(basis_functions_varied)+'< for basis_functions_varied.')
  end select

  end subroutine dorb_dexp_bld

! ===============================================================================
  subroutine grd_dorb_dexp_bld
! -------------------------------------------------------------------------------
! Description   : Gradient of derivatives of occupied orbitals wrt to optimized basis exponents
!
! Created       : J. Toulouse, 27 Jan 2007
! -------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i, bas_i, orb_i, dexp_to_bas_i

! header
  if (header_exe) then

   call object_create ('grd_dorb_dexp')

   call object_needed ('param_exp_nb')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('nbasis')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')
   call object_needed ('orb_occ_last_in_wf_lab')
   call object_needed ('orbital_depends_on_opt_exp')

   return

  endif

! begin

! allocation
  call object_alloc ('grd_dorb_dexp', grd_dorb_dexp, ndim, nelec, orb_occ_last_in_wf_lab, param_exp_nb)

  grd_dorb_dexp (:,:,:,:) = 0.d0

  select case (trim(basis_functions_varied))
  case ('unnormalized')
  call object_provide_by_index (grd_dorb_dexp_bld_index, grd_dphin_dz_index)
  call object_provide_by_index (grd_dorb_dexp_bld_index, coef_index)
  do dexp_i = 1, param_exp_nb
   do orb_i = 1, orb_occ_last_in_wf_lab
    if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
      do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
        bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
        grd_dorb_dexp (:, :, orb_i, dexp_i)= grd_dorb_dexp (:, :, orb_i, dexp_i) + coef (bas_i, orb_i, iwf) * grd_dphin_dz (:, :, bas_i)
      enddo ! dexp_to_bas_i
   enddo ! orb_i
  enddo ! dexp_i

  case ('normalized')
  call object_provide_by_index (grd_dorb_dexp_bld_index, grd_dphin_norm_dz_index)
  call object_provide_by_index (grd_dorb_dexp_bld_index, coef_orb_on_norm_basis_index)
  do dexp_i = 1, param_exp_nb
   do orb_i = 1, orb_occ_last_in_wf_lab
    if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
      do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
        bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
        grd_dorb_dexp (:, :, orb_i, dexp_i)= grd_dorb_dexp (:, :, orb_i, dexp_i) + coef_orb_on_norm_basis (bas_i, orb_i, iwf) * grd_dphin_norm_dz (:, :, bas_i)
      enddo ! dexp_to_bas_i
   enddo ! orb_i
  enddo ! dexp_i

  case ('orthonormalized')
  call object_provide_by_index (grd_dorb_dexp_bld_index, grd_dphin_ortho_dz_index)
  call object_provide_by_index (grd_dorb_dexp_bld_index, coef_orb_on_ortho_basis_index)
  do dexp_i = 1, param_exp_nb
   do orb_i = 1, orb_occ_last_in_wf_lab
    if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
      do bas_i = 1, nbasis
        grd_dorb_dexp (:, :, orb_i, dexp_i)= grd_dorb_dexp (:, :, orb_i, dexp_i) + coef_orb_on_ortho_basis (bas_i, orb_i, iwf) * grd_dphin_ortho_dz (:, :, bas_i, dexp_i)
      enddo ! bas_i
   enddo ! orb_i
  enddo ! dexp_i

  case default
   call die (here, 'unknown case >'+trim(basis_functions_varied)+'< for basis_functions_varied.')
  end select

  end subroutine grd_dorb_dexp_bld

! ===============================================================================
  subroutine lap_dorb_dexp_bld
! -------------------------------------------------------------------------------
! Description   : Laplacian of derivatives of occupied orbitals wrt to optimized basis exponents
!
! Created       : J. Toulouse, 27 Jan 2007
! -------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i, bas_i, orb_i, dexp_to_bas_i

! header
  if (header_exe) then

   call object_create ('lap_dorb_dexp')

   call object_needed ('param_exp_nb')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('nbasis')
   call object_needed ('dexp_to_bas_nb')
   call object_needed ('dexp_to_bas')
   call object_needed ('orb_occ_last_in_wf_lab')
   call object_needed ('orbital_depends_on_opt_exp')

   return

  endif

! begin

! allocation
  call object_alloc ('lap_dorb_dexp', lap_dorb_dexp, nelec, orb_occ_last_in_wf_lab, param_exp_nb)
  lap_dorb_dexp (:,:,:) = 0.d0

  select case (trim(basis_functions_varied))
  case ('unnormalized')
  call object_provide_by_index (lap_dorb_dexp_bld_index, lap_dphin_dz_index)
  call object_provide_by_index (lap_dorb_dexp_bld_index, coef_index)
  do dexp_i = 1, param_exp_nb
   do orb_i = 1, orb_occ_last_in_wf_lab
    if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
      do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
        bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
        lap_dorb_dexp (:, orb_i, dexp_i)= lap_dorb_dexp (:, orb_i, dexp_i) + coef (bas_i, orb_i, iwf) * lap_dphin_dz (:, bas_i)
      enddo ! dexp_to_bas_i
   enddo ! orb_i
  enddo ! dexp_i

  case ('normalized')
  call object_provide_by_index (lap_dorb_dexp_bld_index, lap_dphin_norm_dz_index)
  call object_provide_by_index (lap_dorb_dexp_bld_index, coef_orb_on_norm_basis_index)
  do dexp_i = 1, param_exp_nb
   do orb_i = 1, orb_occ_last_in_wf_lab
    if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
      do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
        bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
        lap_dorb_dexp (:, orb_i, dexp_i)= lap_dorb_dexp (:, orb_i, dexp_i) + coef_orb_on_norm_basis (bas_i, orb_i, iwf) * lap_dphin_norm_dz (:, bas_i)
      enddo ! dexp_to_bas_i
   enddo ! orb_i
  enddo ! dexp_i

  case ('orthonormalized')
  call object_provide_by_index (lap_dorb_dexp_bld_index, lap_dphin_ortho_dz_index)
  call object_provide_by_index (lap_dorb_dexp_bld_index, coef_orb_on_ortho_basis_index)
  do dexp_i = 1, param_exp_nb
   do orb_i = 1, orb_occ_last_in_wf_lab
    if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
      do bas_i = 1, nbasis
        lap_dorb_dexp (:, orb_i, dexp_i)= lap_dorb_dexp (:, orb_i, dexp_i) + coef_orb_on_ortho_basis (bas_i, orb_i, iwf) * lap_dphin_ortho_dz (:, bas_i, dexp_i)
      enddo ! bas_i
   enddo ! orb_i
  enddo ! dexp_i

  case default
   call die (here, 'unknown case >'+trim(basis_functions_varied)+'< for basis_functions_varied.')
  end select

  end subroutine lap_dorb_dexp_bld

! ==============================================================================
  subroutine ddet_dexp_unq_bld
! ------------------------------------------------------------------------------
! Description   : derivatives of unique spin-up and spin-down determinants
! Description   : with respect to basis exponents
! Description   : by updating reference determinant via the Sherman-Morrison formula
!
! Created       : J. Toulouse, 25 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer det_unq_up_i, det_unq_dn_i
  integer orb_i, col_i, i
  integer dexp_i
  real(dp) factor_up, factor_dn

! header
  if (header_exe) then

   call object_create ('ddet_dexp_unq_up')
   call object_create ('ddet_dexp_unq_dn')
   call object_create ('ddet_dexp_col_unq_up')
   call object_create ('ddet_dexp_col_unq_dn')

   call object_needed ('param_exp_nb')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('dorb_dexp')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('orbital_depends_on_opt_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('ddet_dexp_unq_up', ddet_dexp_unq_up, ndetup, param_exp_nb)
  call object_alloc ('ddet_dexp_unq_dn', ddet_dexp_unq_dn, ndetdn, param_exp_nb)
  call object_alloc ('ddet_dexp_col_unq_up', ddet_dexp_col_unq_up, nup, ndetup, param_exp_nb)
  call object_alloc ('ddet_dexp_col_unq_dn', ddet_dexp_col_unq_dn, ndn, ndetdn, param_exp_nb)

! loop over optimized exponents
  do dexp_i = 1, param_exp_nb

!  loop over unique spin-up determinants
   do det_unq_up_i = 1, ndetup

     ddet_dexp_unq_up (det_unq_up_i, dexp_i) = 0.d0

!    loop over columns (=orbitals) of determinant
     do col_i = 1, nup

       orb_i = det_unq_orb_lab_srt_up (col_i, det_unq_up_i)
       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle

       factor_up = 0.d0
       do i = 1, nup
        factor_up = factor_up + slater_mat_trans_inv_up (i, col_i, det_unq_up_i) * dorb_dexp (i, orb_i, dexp_i)
       enddo

       ddet_dexp_col_unq_up (col_i, det_unq_up_i, dexp_i) = factor_up * detu (det_unq_up_i)
       ddet_dexp_unq_up (det_unq_up_i, dexp_i) = ddet_dexp_unq_up (det_unq_up_i, dexp_i) + ddet_dexp_col_unq_up (col_i, det_unq_up_i, dexp_i)

     enddo ! col_i

   enddo ! det_unq_up_i

!  loop over unique spin-dn determinants
   do det_unq_dn_i = 1, ndetdn

     ddet_dexp_unq_dn (det_unq_dn_i, dexp_i) = 0.d0

!    loop over columns (=orbitals) of determinant
     do col_i = 1, ndn

       orb_i = det_unq_orb_lab_srt_dn (col_i, det_unq_dn_i)
       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle

       factor_dn = 0.d0
       do i = 1, ndn
        factor_dn = factor_dn + slater_mat_trans_inv_dn (i, col_i, det_unq_dn_i) * dorb_dexp (nup + i, orb_i, dexp_i)
       enddo

       ddet_dexp_col_unq_dn (col_i, det_unq_dn_i, dexp_i) = factor_dn * detd (det_unq_dn_i)
       ddet_dexp_unq_dn (det_unq_dn_i, dexp_i) = ddet_dexp_unq_dn (det_unq_dn_i, dexp_i) + ddet_dexp_col_unq_dn (col_i, det_unq_dn_i, dexp_i)

     enddo ! col_i

   enddo ! det_unq_dn_i

  enddo ! dexp_i

  end subroutine ddet_dexp_unq_bld

! ==============================================================================
  subroutine dpsi_exp_bld
! ------------------------------------------------------------------------------
! Description   : Logarithm derivatives of Psi with respect to basis exponents
!
! Created       : J. Toulouse, 25 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i
  integer csf_i, det_in_csf_i, det_i
  integer det_unq_up_i, det_unq_dn_i

! header
  if (header_exe) then

   call object_create ('dpsi_exp', dpsi_exp_index)
   call object_create ('dpsid_exp')

   call object_needed ('param_exp_nb')
   call object_needed ('ncsf')
   call object_needed ('ndet')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('csf_coef')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('ddet_dexp_unq_up')
   call object_needed ('ddet_dexp_unq_dn')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('psido')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsid_exp', dpsid_exp, param_exp_nb)
  call object_alloc ('dpsi_exp', dpsi_exp, param_exp_nb)

  dpsid_exp (:) = 0.d0

! loop over optimized exponents
  do dexp_i = 1, param_exp_nb

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

        det_i = iwdet_in_csf (det_in_csf_i, csf_i)
        det_unq_up_i = det_to_det_unq_up (det_i)
        det_unq_dn_i = det_to_det_unq_dn (det_i)

        dpsid_exp (dexp_i) = dpsid_exp (dexp_i) + csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i) *  &
        (ddet_dexp_unq_up (det_unq_up_i, dexp_i) * detd (det_unq_dn_i) + detu (det_unq_up_i) * ddet_dexp_unq_dn (det_unq_dn_i, dexp_i))

     enddo ! det_in_csf_i
   enddo ! csf_i

    dpsi_exp (dexp_i) = dpsid_exp (dexp_i) / psido

  enddo ! dexp_i

! tests for He
!  write(6,'(2a,f)') trim(here), ': dpsi_exp =', dpsi_exp (1)
!  write(6,'(2a,f)') trim(here), ': check dpsi_exp=', dorb_dexp(1,1,1)/orb(1,1) + dorb_dexp(2,1,1)/orb(2,1)

  end subroutine dpsi_exp_bld

! ==============================================================================
  subroutine dpsi_lnexp_bld
! ------------------------------------------------------------------------------
! Description   : logarithm derivatives of Psi with respect to logarithm of basis exponents
!
! Created       : J. Toulouse, 23 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i, bas_i

! header
  if (header_exe) then

   call object_create ('dpsi_lnexp', dpsi_lnexp_index)

   call object_needed ('param_exp_nb')
   call object_needed ('dpsi_exp')
   call object_needed ('dexp_to_bas')
   call object_needed ('zex')

   return

  endif

! begin

! allocations
  call object_alloc ('dpsi_lnexp', dpsi_lnexp, param_exp_nb)

  do dexp_i = 1, param_exp_nb
     bas_i = dexp_to_bas (dexp_i)%row (1)
     dpsi_lnexp (dexp_i) = dpsi_exp (dexp_i) * zex (bas_i, iwf)
  enddo ! dexp_i

  end subroutine dpsi_lnexp_bld

! ==============================================================================
  subroutine slater_mat_exp_trans_inv_bld
! ------------------------------------------------------------------------------
! Description   : inverse of transpose of derivatives of Slater matrices wrt to basis exponents
! Description   : using the Sherman-Morison formula
!
! Created       : J. Toulouse, 26 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i
  integer det_unq_up_i, det_unq_dn_i
  integer i, j, k, l
  integer col_i, orb_i
  real(dp) factor_up_inv, factor_dn_inv
  real(dp), allocatable :: ratio_up (:), ratio_dn (:)

! header
  if (header_exe) then

   call object_create ('slater_mat_exp_trans_inv_up')
   call object_create ('slater_mat_exp_trans_inv_dn')

   call object_needed ('param_exp_nb')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('slater_mat_trans_inv_up')
   call object_needed ('slater_mat_trans_inv_dn')
   call object_needed ('dorb_dexp')
   call object_needed ('orbital_depends_on_opt_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('slater_mat_exp_trans_inv_up', slater_mat_exp_trans_inv_up, nup, nup, nup, ndetup, param_exp_nb)
  call object_alloc ('slater_mat_exp_trans_inv_dn', slater_mat_exp_trans_inv_dn, ndn, ndn, ndn, ndetdn, param_exp_nb)
  call alloc ('ratio_up', ratio_up, nup)
  call alloc ('ratio_dn', ratio_dn, ndn)

! loop over optimized exponents
  do dexp_i = 1, param_exp_nb

!  loop over unique spin-up determinants
   do det_unq_up_i = 1, ndetup

!    loop over columns (=orbitals) of determinant
     do col_i = 1, nup

       orb_i = det_unq_orb_lab_srt_up (col_i, det_unq_up_i)
       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle

       do l = 1, nup
        ratio_up (l) = 0.d0
        do i = 1, nup
         ratio_up (l) = ratio_up (l) + slater_mat_trans_inv_up (i, l, det_unq_up_i) * dorb_dexp (i, orb_i, dexp_i)
        enddo  ! i
       enddo ! l

      factor_up_inv = 1.d0/ratio_up (col_i)
      ratio_up (col_i) = ratio_up (col_i) - 1.d0

      do k = 1, nup
        do j = 1 , nup
         slater_mat_exp_trans_inv_up (k, j, col_i, det_unq_up_i, dexp_i) = slater_mat_trans_inv_up (k, j, det_unq_up_i) &
                                 - slater_mat_trans_inv_up (k, col_i, det_unq_up_i) * ratio_up (j) * factor_up_inv
        enddo ! j
      enddo ! k


     enddo ! col_i

   enddo ! det_unq_up_i

!  loop over unique spin-dn determinants
   do det_unq_dn_i = 1, ndetdn

!    loop over columns (=orbitals) of determinant
     do col_i = 1, ndn

       orb_i = det_unq_orb_lab_srt_dn (col_i, det_unq_dn_i)
       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle

       do l = 1, ndn
        ratio_dn (l) = 0.d0
        do i = 1, ndn
         ratio_dn (l) = ratio_dn (l) + slater_mat_trans_inv_dn (i, l, det_unq_dn_i) * dorb_dexp (nup+i, orb_i, dexp_i)
        enddo  ! i
       enddo ! l

       factor_dn_inv = 1.d0/ratio_dn (col_i)
       ratio_dn (col_i) = ratio_dn (col_i) - 1.d0

      do k = 1, ndn
        do j = 1 , ndn
         slater_mat_exp_trans_inv_dn (k, j, col_i, det_unq_dn_i, dexp_i) = slater_mat_trans_inv_dn (k, j, det_unq_dn_i) &
                                 - slater_mat_trans_inv_dn (k, col_i, det_unq_dn_i) * ratio_dn (j) * factor_dn_inv
        enddo ! j
      enddo ! k

     enddo ! col_i

   enddo ! det_unq_dn_i

   enddo ! dexp_i

  end subroutine slater_mat_exp_trans_inv_bld

!! ==============================================================================
!  subroutine slater_mat_exp_trans_inv_2_bld
!! ------------------------------------------------------------------------------
!! Description   : inverse of transpose of derivatives of Slater matrices wrt to basis exponents
!! Description   : using direct inversion
!!
!! Created       : J. Toulouse, 18 May 2007
!! ------------------------------------------------------------------------------
!  include 'modules.h'
!  implicit none
!
!! local
!  integer dexp_i
!  integer det_unq_up_i, det_unq_dn_i
!  integer det_i
!  integer i, j, k, l
!  integer col_i, orb_i
!  real(dp) factor_up_inv, factor_dn_inv
!  real(dp), allocatable :: ratio_up (:), ratio_dn (:)
!
!! header
!  if (header_exe) then
!
!   call object_create ('slater_mat_exp_trans_inv_up_2')
!   call object_create ('slater_mat_exp_trans_inv_dn_2')
!
!   call object_needed ('param_exp_nb')
!   call object_needed ('ndetup')
!   call object_needed ('ndetdn')
!   call object_needed ('nup')
!   call object_needed ('ndn')
!   call object_needed ('slater_mat_trans_inv_up')
!   call object_needed ('slater_mat_trans_inv_dn')
!   call object_needed ('dorb_dexp')
!   call object_needed ('orbital_depends_on_opt_exp')
!
!   return
!
!  endif
!
!! begin
!
!! allocations
!  call object_alloc ('slater_mat_exp_trans_inv_up_2', slater_mat_exp_trans_inv_up_2, nup, nup, nup, ndetup, param_exp_nb)
!  call object_alloc ('slater_mat_exp_trans_inv_dn_2', slater_mat_exp_trans_inv_dn_2, ndn, ndn, ndn, ndetdn, param_exp_nb)
!  call object_alloc ('slater_mat_exp_trans_up', slater_mat_exp_trans_up, nup, nup, nup, ndetup, param_exp_nb)
!  call object_alloc ('slater_mat_exp_trans_dn', slater_mat_exp_trans_dn, ndn, ndn, ndn, ndetdn, param_exp_nb)
!  call object_alloc ('mat_flat_up', mat_flat_up, nup*nup, nup, ndetup, param_exp_nb)
!  call object_alloc ('mat_flat_dn', mat_flat_dn, ndn*ndn, ndn, ndetdn, param_exp_nb)
!
!! loop over optimized exponents
!  do dexp_i = 1, param_exp_nb
!
!!  loop over unique spin-up determinants
!   do det_unq_up_i = 1, ndetup
!
!     det_i = det_unq_up_to_det (det_unq_up_i)
!
!!    loop over columns (=orbitals) of determinant
!     do col_i = 1, nup
!
!       orb_i = det_unq_orb_lab_srt_up (col_i, det_unq_up_i)
!       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
!
!    do orb_i = 1, nup
!      do elec_i = 1, nup
!!        write(6,*) trim(here),': det_i=',det_i
!!        write(6,*) trim(here),': orb_i=',orb_i
!!        write(6,*) trim(here),': elec_i=',elec_i
!!        write(6,*) trim(here),': det_ex_unq_orb_lab_up (orb_i, det_i)=',det_ex_unq_orb_lab_up (orb_i, det_i)
!!        write(6,*) trim(here),': orb=', orb (elec_i, det_ex_unq_orb_lab_up (orb_i, det_i))
!        slater_mat_ex_trans_up (orb_i, elec_i, det_i) = orb (elec_i, det_ex_unq_orb_lab_up (orb_i, det_i))
!      enddo
!    enddo
!
!       slater_mat_exp_trans_up (1:nup, 1:nup, col_i, det_unq_up_i, dexp_i) = orb (1:nup, det_ex_unq_orb_lab_up (orb_i, det_i))
!       slater_mat_exp_trans_up (1:nup, col_i, col_i, det_unq_up_i, dexp_i) = dorb_dexp (1:nup, orb_i, dexp_i)
!
!       call flatten (mat_flat_up (:,:, col_i, det_unq_up_i, dexp_i), slater_mat_exp_trans_up (:,:,col_i, det_unq_up_i, dexp_i), nup, nup)
!       call matinv (mat_flat_up (:,:, col_i, det_unq_up_i, dexp_i), nup, det)
!       call unflatten (mat_flat_up (:,:, col_i, det_unq_up_i, dexp_i), slater_mat_exp_trans_inv_up_2 (:,:,col_i, det_unq_up_i, dexp_i), nup, nup)
!
!     enddo ! col_i
!
!   enddo ! det_unq_up_i
!
!!  loop over unique spin-dn determinants
!   do det_unq_dn_i = 1, ndetdn
!
!     det_i = det_unq_dn_to_det (det_unq_dn_i)
!
!!    loop over columns (=orbitals) of determinant
!     do col_i = 1, ndn
!
!       orb_i = det_unq_orb_lab_srt_dn (col_i, det_unq_dn_i)
!       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle
!
!       slater_mat_exp_trans_dn (1:ndn, 1:ndn, col_i, det_unq_dn_i, dexp_i) = slater_mat_trans_inv_dn (1:ndn, 1:ndn, det_i)
!       slater_mat_exp_trans_dn (1:ndn, col_i, col_i, det_unq_dn_i, dexp_i) = dorb_dexp (1:ndn, orb_i, dexp_i)
!
!       call flatten (mat_flat_up (:,:, col_i, det_unq_dn_i, dexp_i), slater_mat_exp_trans_up (:,:,col_i, det_unq_dn_i, dexp_i), ndn, ndn)
!       call matinv (mat_flat_up (:,:, col_i, det_unq_dn_i, dexp_i), ndn, det)
!       call unflatten (mat_flat_up (:,:, col_i, det_unq_dn_i, dexp_i), slater_mat_exp_trans_inv_dn_2 (:,:,col_i, det_unq_dn_i, dexp_i), ndn, ndn)
!
!     enddo ! col_i
!
!   enddo ! det_unq_dn_i
!
!   enddo ! dexp_i
!
!   write(6,'(a,100f12.6)') 'slater_mat_exp_trans_inv_up=', slater_mat_exp_trans_inv_up
!
!  end subroutine slater_mat_exp_trans_inv_bld

! ==============================================================================
  subroutine grd_ddet_dexp_unq_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of derivatives of unique spin-up and spin-down determinants
! Description   : with respect to basis exponents
! Description   : by updating reference determinant via the Sherman-Morrison formula
!
! Created       : J. Toulouse, 27 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i
  integer det_unq_up_i, det_unq_dn_i
  integer dim_i, elec_up_i, elec_dn_i
  integer col_i, col_j
  integer orb_i, orb_j
  real(dp) factor_up, factor_dn

! header
  if (header_exe) then

   call object_create ('grd_ddet_dexp_unq_up')
   call object_create ('grd_ddet_dexp_unq_dn')

   call object_needed ('param_exp_nb')
   call object_needed ('ndim')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('slater_mat_exp_trans_inv_up')
   call object_needed ('slater_mat_exp_trans_inv_dn')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('dorb')
   call object_needed ('grd_dorb_dexp')
   call object_needed ('ddet_dexp_col_unq_up')
   call object_needed ('ddet_dexp_col_unq_dn')
   call object_needed ('orbital_depends_on_opt_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_ddet_dexp_unq_up', grd_ddet_dexp_unq_up, ndim, nup, ndetup, param_exp_nb)
  call object_alloc ('grd_ddet_dexp_unq_dn', grd_ddet_dexp_unq_dn, ndim, ndn, ndetdn, param_exp_nb)

  grd_ddet_dexp_unq_up (:,:,:,:) = 0.d0
  grd_ddet_dexp_unq_dn (:,:,:,:) = 0.d0

! loop over optimized exponents
  do dexp_i = 1, param_exp_nb

!  loop over unique spin-up determinants
   do det_unq_up_i = 1, ndetup

!   loop over dimensions
    do dim_i = 1, ndim

!    loop over spin-up electrons
     do elec_up_i = 1, nup

!     loop over columns (=orbitals) of determinant
      do col_i = 1, nup

       orb_i = det_unq_orb_lab_srt_up (col_i, det_unq_up_i)
       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle

!      Sherman-Morison formula
       factor_up = 0.d0
       do col_j = 1, nup
        orb_j = det_unq_orb_lab_srt_up (col_j, det_unq_up_i)
        if (col_j == col_i) then
         factor_up = factor_up + slater_mat_exp_trans_inv_up (elec_up_i, col_i, col_i, det_unq_up_i, dexp_i) * grd_dorb_dexp (dim_i, elec_up_i, orb_i, dexp_i)
        else
         factor_up = factor_up + slater_mat_exp_trans_inv_up (elec_up_i, col_j, col_i, det_unq_up_i, dexp_i) * dorb (dim_i, elec_up_i, orb_j)
        endif
       enddo ! col_j

        grd_ddet_dexp_unq_up (dim_i, elec_up_i, det_unq_up_i, dexp_i) = grd_ddet_dexp_unq_up (dim_i, elec_up_i, det_unq_up_i, dexp_i)   &
                                                                     + factor_up * ddet_dexp_col_unq_up (col_i, det_unq_up_i, dexp_i)

      enddo ! col_i

     enddo ! elec_up_i

    enddo ! dim_i

   enddo ! det_unq_up_i

!  loop over unique spin-dn determinants
   do det_unq_dn_i = 1, ndetdn

!   loop over dimensions
    do dim_i = 1, ndim

!    loop over ndn
     do elec_dn_i = 1, ndn

!     loop over columns (=orbitals) of determinant
      do col_i = 1, ndn

       orb_i = det_unq_orb_lab_srt_dn (col_i, det_unq_dn_i)
       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle

!      Sherman-Morison formula
       factor_dn = 0.d0
       do col_j = 1, ndn
        orb_j = det_unq_orb_lab_srt_dn (col_j, det_unq_dn_i)
        if (col_j == col_i) then
         factor_dn = factor_dn + slater_mat_exp_trans_inv_dn (elec_dn_i, col_i, col_i, det_unq_dn_i, dexp_i) * grd_dorb_dexp (dim_i, nup + elec_dn_i, orb_i, dexp_i)
        else
         factor_dn = factor_dn + slater_mat_exp_trans_inv_dn (elec_dn_i, col_j, col_i, det_unq_dn_i, dexp_i) * dorb (dim_i, nup + elec_dn_i, orb_j)
        endif
       enddo ! col_j

        grd_ddet_dexp_unq_dn (dim_i, elec_dn_i, det_unq_dn_i, dexp_i) = grd_ddet_dexp_unq_dn (dim_i, elec_dn_i, det_unq_dn_i, dexp_i)   &
                                                                     + factor_dn * ddet_dexp_col_unq_dn (col_i, det_unq_dn_i, dexp_i)

      enddo ! col_i

     enddo ! elec_dn_i

    enddo ! dim_i

   enddo ! det_unq_dn_i

   enddo ! dexp_i

  end subroutine grd_ddet_dexp_unq_bld

! ==============================================================================
  subroutine lap_ddet_dexp_unq_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of derivatives of unique spin-up and spin-down determinants
! Description   : with respect to basis exponents
! Description   : by updating reference determinant via the Sherman-Morrison formula
!
! Created       : J. Toulouse, 27 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i
  integer det_unq_up_i, det_unq_dn_i
  integer elec_up_i, elec_dn_i
  integer col_i, col_j
  integer orb_i, orb_j
  real(dp) factor_up, factor_dn

! header
  if (header_exe) then

   call object_create ('lap_ddet_dexp_unq_up')
   call object_create ('lap_ddet_dexp_unq_dn')

   call object_needed ('param_exp_nb')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ndetup')
   call object_needed ('ndetdn')
   call object_needed ('slater_mat_exp_trans_inv_up')
   call object_needed ('slater_mat_exp_trans_inv_dn')
   call object_needed ('det_unq_orb_lab_srt_up')
   call object_needed ('det_unq_orb_lab_srt_dn')
   call object_needed ('ddorb')
   call object_needed ('lap_dorb_dexp')
   call object_needed ('ddet_dexp_col_unq_up')
   call object_needed ('ddet_dexp_col_unq_dn')
   call object_needed ('orbital_depends_on_opt_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_ddet_dexp_unq_up', lap_ddet_dexp_unq_up, nup, ndetup, param_exp_nb)
  call object_alloc ('lap_ddet_dexp_unq_dn', lap_ddet_dexp_unq_dn, ndn, ndetdn, param_exp_nb)

  lap_ddet_dexp_unq_up (:,:,:) = 0.d0
  lap_ddet_dexp_unq_dn (:,:,:) = 0.d0

! loop over optimized exponents
  do dexp_i = 1, param_exp_nb

!  loop over unique spin-up determinants
   do det_unq_up_i = 1, ndetup

!    loop over spin-up electrons
     do elec_up_i = 1, nup

!     loop over columns (=orbitals) of determinant
      do col_i = 1, nup

       orb_i = det_unq_orb_lab_srt_up (col_i, det_unq_up_i)
       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle

!      Sherman-Morison formula
       factor_up = 0.d0
       do col_j = 1, nup
        orb_j = det_unq_orb_lab_srt_up (col_j, det_unq_up_i)
        if (col_j == col_i) then
         factor_up = factor_up + slater_mat_exp_trans_inv_up (elec_up_i, col_i, col_i, det_unq_up_i, dexp_i) * lap_dorb_dexp (elec_up_i, orb_i, dexp_i)
        else
         factor_up = factor_up + slater_mat_exp_trans_inv_up (elec_up_i, col_j, col_i, det_unq_up_i, dexp_i) * ddorb (elec_up_i, orb_j)
        endif
       enddo ! col_j

        lap_ddet_dexp_unq_up (elec_up_i, det_unq_up_i, dexp_i) = lap_ddet_dexp_unq_up (elec_up_i, det_unq_up_i, dexp_i)   &
                                                                     + factor_up * ddet_dexp_col_unq_up (col_i, det_unq_up_i, dexp_i)

      enddo ! col_i

     enddo ! elec_up_i

   enddo ! det_unq_up_i


!  loop over unique spin-dn determinants
   do det_unq_dn_i = 1, ndetdn

!    loop over ndn
     do elec_dn_i = 1, ndn

!     loop over columns (=orbitals) of determinant
      do col_i = 1, ndn

       orb_i = det_unq_orb_lab_srt_dn (col_i, det_unq_dn_i)
       if (.not. orbital_depends_on_opt_exp (orb_i, dexp_i)) cycle

!      Sherman-Morison formula
       factor_dn = 0.d0
       do col_j = 1, ndn
        orb_j = det_unq_orb_lab_srt_dn (col_j, det_unq_dn_i)
        if (col_j == col_i) then
         factor_dn = factor_dn + slater_mat_exp_trans_inv_dn (elec_dn_i, col_i, col_i, det_unq_dn_i, dexp_i) * lap_dorb_dexp (nup + elec_dn_i, orb_i, dexp_i)
        else
         factor_dn = factor_dn + slater_mat_exp_trans_inv_dn (elec_dn_i, col_j, col_i, det_unq_dn_i, dexp_i) * ddorb (nup + elec_dn_i, orb_j)
        endif
       enddo ! col_j

        lap_ddet_dexp_unq_dn (elec_dn_i, det_unq_dn_i, dexp_i) = lap_ddet_dexp_unq_dn (elec_dn_i, det_unq_dn_i, dexp_i)   &
                                                                     + factor_dn * ddet_dexp_col_unq_dn (col_i, det_unq_dn_i, dexp_i)

      enddo ! col_i

     enddo ! elec_dn_i

   enddo ! det_unq_dn_i

   enddo ! dexp_i

  end subroutine lap_ddet_dexp_unq_bld

! ==============================================================================
  subroutine grd_dpsid_exp_over_dpsid_exp_bld
! ------------------------------------------------------------------------------
! Description   : Gradient of derivative of determinantal part of Psi with respect to
! Description   : basis exponents over derivative of determinantal part of Psi
! Description   : with respect to basis exponents
!
! Created       : J. Toulouse, 29 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i
  integer csf_i, det_in_csf_i, det_i
  integer det_unq_up_i, det_unq_dn_i
  integer dim_i
  integer elec_i, elec_up_i, elec_dn_i
  real(dp) coefficient


! header
  if (header_exe) then

   call object_create ('grd_dpsid_exp_over_dpsid_exp')

   call object_needed ('param_exp_nb')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ncsf')
   call object_needed ('ndet')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('ddet_dexp_unq_up')
   call object_needed ('ddet_dexp_unq_dn')
   call object_needed ('grd_ddet_dexp_unq_up')
   call object_needed ('grd_ddet_dexp_unq_dn')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('grd_det_unq_up')
   call object_needed ('grd_det_unq_dn')
   call object_needed ('dpsid_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_dpsid_exp_over_dpsid_exp', grd_dpsid_exp_over_dpsid_exp, ndim, nelec, param_exp_nb)

  grd_dpsid_exp_over_dpsid_exp (:,:,:) = 0.d0

! loop over optimized exponents
  do dexp_i = 1, param_exp_nb

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

        det_i = iwdet_in_csf (det_in_csf_i, csf_i)
        det_unq_up_i = det_to_det_unq_up (det_i)
        det_unq_dn_i = det_to_det_unq_dn (det_i)

        coefficient = csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)

!       loop over dimensions
        do dim_i = 1, ndim

         elec_i = 0

!        gradient wrt to spin-up electrons
         do elec_up_i = 1, nup
          elec_i = elec_i + 1

          grd_dpsid_exp_over_dpsid_exp (dim_i, elec_i, dexp_i) = grd_dpsid_exp_over_dpsid_exp (dim_i, elec_i, dexp_i)       &
                      + coefficient * ( grd_ddet_dexp_unq_up (dim_i, elec_up_i, det_unq_up_i, dexp_i) * detd (det_unq_dn_i)     &
                                       + grd_det_unq_up (dim_i, elec_up_i, det_unq_up_i) * ddet_dexp_unq_dn (det_unq_dn_i, dexp_i) )
         enddo ! elec_up_i

!        gradient wrt to spin-dn electrons
         do elec_dn_i = 1, ndn
          elec_i = elec_i + 1

          grd_dpsid_exp_over_dpsid_exp (dim_i, elec_i, dexp_i) = grd_dpsid_exp_over_dpsid_exp (dim_i, elec_i, dexp_i)   &
                       + coefficient * ( ddet_dexp_unq_up (det_unq_up_i, dexp_i) * grd_det_unq_dn (dim_i, elec_dn_i, det_unq_dn_i) &
                                        + detu (det_unq_up_i) * grd_ddet_dexp_unq_dn (dim_i, elec_dn_i, det_unq_dn_i, dexp_i) )
         enddo ! elec_dn_i

        enddo ! dim_i

     enddo ! det_in_csf_i
   enddo ! csf_i

!  write(6,'(2a,f)') trim(here), ': grd_dpsid_exp_over_dpsid_exp(1,1,1)=', grd_dpsid_exp_over_dpsid_exp(1,1,1)
!  write(6,'(2a,f)') trim(here), ': check =', (grd_dorb_dexp(1,1,1,1)*orb(2,1) + dorb(1,1,1)*dorb_dexp(2,1,1))
   grd_dpsid_exp_over_dpsid_exp (:,:,dexp_i) = grd_dpsid_exp_over_dpsid_exp (:,:,dexp_i) / dpsid_exp (dexp_i)

  enddo ! dexp_i

! tests for He
!  write(6,'(2a,f)') trim(here), ': grd_dpsid_exp_over_dpsid_exp(1,1,1)=', grd_dpsid_exp_over_dpsid_exp(1,1,1)
!  write(6,'(2a,f)') trim(here), ': check =', (grd_dorb_dexp(1,1,1,1)*orb(2,1) + dorb(1,1,1)*dorb_dexp(2,1,1))/(dorb_dexp(1,1,1)*orb(2,1) + orb(1,1)*dorb_dexp(2,1,1))

  end subroutine grd_dpsid_exp_over_dpsid_exp_bld

! ==============================================================================
  subroutine lap_dpsid_exp_over_dpsid_exp_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian of derivative of determinantal part of Psi with respect to
! Description   : basis exponents over derivative of determinantal part of Psi
! Description   : with respect to basis exponents
!
! Created       : J. Toulouse, 29 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i
  integer csf_i, det_in_csf_i, det_i
  integer det_unq_up_i, det_unq_dn_i
  integer elec_i, elec_up_i, elec_dn_i
  real(dp) coefficient


! header
  if (header_exe) then

   call object_create ('lap_dpsid_exp_over_dpsid_exp')

   call object_needed ('param_exp_nb')
   call object_needed ('nelec')
   call object_needed ('nup')
   call object_needed ('ndn')
   call object_needed ('ncsf')
   call object_needed ('ndet')
   call object_needed ('ndet_in_csf')
   call object_needed ('iwdet_in_csf')
   call object_needed ('cdet_in_csf')
   call object_needed ('det_to_det_unq_up')
   call object_needed ('det_to_det_unq_dn')
   call object_needed ('ddet_dexp_unq_up')
   call object_needed ('ddet_dexp_unq_dn')
   call object_needed ('lap_ddet_dexp_unq_up')
   call object_needed ('lap_ddet_dexp_unq_dn')
   call object_needed ('detu')
   call object_needed ('detd')
   call object_needed ('lap_det_unq_up')
   call object_needed ('lap_det_unq_dn')
   call object_needed ('dpsid_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_dpsid_exp_over_dpsid_exp', lap_dpsid_exp_over_dpsid_exp, nelec, param_exp_nb)

  lap_dpsid_exp_over_dpsid_exp (:,:) = 0.d0

! loop over optimized exponents
  do dexp_i = 1, param_exp_nb

   do csf_i = 1, ncsf

     do det_in_csf_i = 1, ndet_in_csf (csf_i)

        det_i = iwdet_in_csf (det_in_csf_i, csf_i)
        det_unq_up_i = det_to_det_unq_up (det_i)
        det_unq_dn_i = det_to_det_unq_dn (det_i)

        coefficient = csf_coef (csf_i, 1) * cdet_in_csf (det_in_csf_i, csf_i)

         elec_i = 0

!        laplacian wrt to spin-up electrons
         do elec_up_i = 1, nup
          elec_i = elec_i + 1

          lap_dpsid_exp_over_dpsid_exp (elec_i, dexp_i) = lap_dpsid_exp_over_dpsid_exp (elec_i, dexp_i)       &
                      + coefficient * ( lap_ddet_dexp_unq_up (elec_up_i, det_unq_up_i, dexp_i) * detd (det_unq_dn_i)     &
                                       + lap_det_unq_up (elec_up_i, det_unq_up_i) * ddet_dexp_unq_dn (det_unq_dn_i, dexp_i) )
         enddo ! elec_up_i

!        laplacian wrt to spin-dn electrons
         do elec_dn_i = 1, ndn
          elec_i = elec_i + 1

          lap_dpsid_exp_over_dpsid_exp (elec_i, dexp_i) = lap_dpsid_exp_over_dpsid_exp (elec_i, dexp_i)   &
                       + coefficient * ( ddet_dexp_unq_up (det_unq_up_i, dexp_i) * lap_det_unq_dn (elec_dn_i, det_unq_dn_i) &
                                        + detu (det_unq_up_i) * lap_ddet_dexp_unq_dn (elec_dn_i, det_unq_dn_i, dexp_i) )
         enddo ! elec_dn_i

     enddo ! det_in_csf_i
   enddo ! csf_i

   lap_dpsid_exp_over_dpsid_exp (:,dexp_i) = lap_dpsid_exp_over_dpsid_exp (:,dexp_i) / dpsid_exp (dexp_i)

  enddo ! dexp_i

! tests for He
!  write(6,'(2a,f)') trim(here), ': lap_dpsid_exp_over_dpsid_exp(1,1)=', lap_dpsid_exp_over_dpsid_exp(1,1)
!  write(6,'(2a,f)') trim(here), ': check =', (lap_dorb_dexp(1,1,1)*orb(2,1) + ddorb(1,1)*dorb_dexp(2,1,1))/(dorb_dexp(1,1,1)*orb(2,1) + orb(1,1)*dorb_dexp(2,1,1))

  end subroutine lap_dpsid_exp_over_dpsid_exp_bld

! ==============================================================================
  subroutine lap_ln_dpsid_exp_bld
! ------------------------------------------------------------------------------
! Description   : Laplacian ( ln (psid_exp) )
!
! Created       : J. Toulouse, 29 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dim_i, elec_i, dexp_i
  real(dp) sum

! header
  if (header_exe) then

   call object_create ('lap_ln_dpsid_exp')

   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('param_exp_nb')
   call object_needed ('lap_dpsid_exp_over_dpsid_exp')
   call object_needed ('grd_dpsid_exp_over_dpsid_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('lap_ln_dpsid_exp', lap_ln_dpsid_exp, nelec, param_exp_nb)

  do dexp_i = 1, param_exp_nb
   do elec_i = 1, nelec
    sum = 0.d0
    do dim_i = 1, ndim
     sum = sum + grd_dpsid_exp_over_dpsid_exp (dim_i, elec_i, dexp_i)**2
    enddo ! dim_i
    lap_ln_dpsid_exp (elec_i, dexp_i) = lap_dpsid_exp_over_dpsid_exp (elec_i, dexp_i) - sum
   enddo ! elec_i
  enddo ! dexp_i

  end subroutine lap_ln_dpsid_exp_bld

! ==============================================================================
  subroutine sum_lap_ln_dpsid_exp_bld
! ------------------------------------------------------------------------------
! Description   :  Sum_i lap_ln_dpsid_exp (i)
!
! Created       : J. Toulouse, 29 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer elec_i, dexp_i

! header
  if (header_exe) then

   call object_create ('sum_lap_ln_dpsid_exp')

   call object_needed ('param_exp_nb')
   call object_needed ('nelec')
   call object_needed ('lap_ln_dpsid_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('sum_lap_ln_dpsid_exp', sum_lap_ln_dpsid_exp, param_exp_nb)

  sum_lap_ln_dpsid_exp (:) = 0.d0

  do dexp_i = 1, param_exp_nb
   do elec_i = 1, nelec
     sum_lap_ln_dpsid_exp (dexp_i) = sum_lap_ln_dpsid_exp (dexp_i) + lap_ln_dpsid_exp (elec_i, dexp_i)
   enddo
  enddo

  end subroutine sum_lap_ln_dpsid_exp_bld

! ==============================================================================
  subroutine grd_dpsi_exp_over_dpsi_exp_bld
! ------------------------------------------------------------------------------
! Description   : (Gradient dpsi_exp)/ dpsi_exp
!
! Created       : J. Toulouse, 29 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i

! header
  if (header_exe) then

   call object_create ('grd_dpsi_exp_over_dpsi_exp')

   call object_needed ('param_exp_nb')
   call object_needed ('nelec')
   call object_needed ('ndim')
   call object_needed ('grd_dpsid_exp_over_dpsid_exp')
   call object_needed ('vj')

   return

  endif

! begin

! allocations
  call object_alloc ('grd_dpsi_exp_over_dpsi_exp', grd_dpsi_exp_over_dpsi_exp, ndim, nelec, param_exp_nb)

  do dexp_i = 1, param_exp_nb
     grd_dpsi_exp_over_dpsi_exp (1:ndim, 1:nelec, dexp_i) = grd_dpsid_exp_over_dpsid_exp (1:ndim, 1:nelec, dexp_i) + vj (1:ndim, 1:nelec)
  enddo

  end subroutine grd_dpsi_exp_over_dpsi_exp_bld

! ==============================================================================
  subroutine sum_lap_ln_dpsi_exp_bld
! ------------------------------------------------------------------------------
! Description   : Sum_i lap_ln_dpsi_exp (i) : determinant part + Jastrow
!
! Created       : J. Toulouse, 29 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! header
  if (header_exe) then

   call object_create ('sum_lap_ln_dpsi_exp')

   call object_needed ('param_exp_nb')
   call object_needed ('sum_lap_ln_dpsid_exp')
   call object_needed ('sum_lap_lnj')

   return

  endif

! begin

! allocations
  call object_alloc ('sum_lap_ln_dpsi_exp', sum_lap_ln_dpsi_exp, param_exp_nb)

  sum_lap_ln_dpsi_exp (:) = sum_lap_ln_dpsid_exp (:) + sum_lap_lnj

  end subroutine sum_lap_ln_dpsi_exp_bld

! ==============================================================================
  subroutine eloc_kin_exp_bld
! ------------------------------------------------------------------------------
! Description   : Kinetic local energy of dpsi_exp
!
! Created       : J. Toulouse, 29 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dim_i, elec_i, dexp_i
  real(dp) sum_grd_psi_over_psi_square

! header
  if (header_exe) then

   call object_create ('eloc_kin_exp')

   call object_needed ('param_exp_nb')
   call object_needed ('ndim')
   call object_needed ('nelec')
   call object_needed ('sum_lap_ln_dpsi_exp')
   call object_needed ('grd_dpsi_exp_over_dpsi_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('eloc_kin_exp', eloc_kin_exp, param_exp_nb)

  do dexp_i = 1, param_exp_nb

  sum_grd_psi_over_psi_square = 0.d0

  do dim_i = 1, ndim
    do elec_i = 1, nelec
       sum_grd_psi_over_psi_square = sum_grd_psi_over_psi_square + grd_dpsi_exp_over_dpsi_exp (dim_i, elec_i, dexp_i)**2
    enddo
  enddo

  eloc_kin_exp (dexp_i) =  -0.5d0 * (sum_lap_ln_dpsi_exp (dexp_i) + sum_grd_psi_over_psi_square)

  enddo ! dexp_i

  end subroutine eloc_kin_exp_bld

! ==============================================================================
  subroutine eloc_exp_bld
! ------------------------------------------------------------------------------
! Description   : Total local energy of dpsi_exp
!
! Created       : J. Toulouse, 29 Jan 2006
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! header
  if (header_exe) then

   call object_create ('eloc_exp')

   call object_needed ('param_exp_nb')
   call object_needed ('eloc_kin_exp')
   call object_needed ('eloc_pot')

   return

  endif

! begin

! not for a pseudopotential
  if (nloc > 0) then
    call die (here, 'local energy derivatives wrt exponents not implemented for a pseudopotential')
  endif

! allocations
  call object_alloc ('eloc_exp', eloc_exp, param_exp_nb)

  eloc_exp (:) = eloc_kin_exp (:) + eloc_pot

  end subroutine eloc_exp_bld

! ==============================================================================
  subroutine deloc_exp_bld
! ------------------------------------------------------------------------------
! Description   : derivative of local energy wrt basis exponents
!
! Created       : J. Toulouse, 29 Jan 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! header
  if (header_exe) then

   call object_create ('deloc_exp', deloc_exp_index)

   call object_needed ('param_exp_nb')
   call object_needed ('eloc_exp')
   call object_needed ('eloc')
   call object_needed ('dpsi_exp')

   return

  endif

! begin

! allocations
  call object_alloc ('deloc_exp', deloc_exp, param_exp_nb)

  deloc_exp (:) = (eloc_exp (:)- eloc) * dpsi_exp (:)

  end subroutine deloc_exp_bld

! ==============================================================================
  subroutine move_zex (dexp_i, d)
! ------------------------------------------------------------------------------
! Description   : numerical derivative of local energy wrt basis exponents
! Description   : for checkings
!
! Created       : J. Toulouse, 25 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i, dexp_to_bas_i, bas_i
  real(dp) :: d

! begin
  call object_provide ('param_exp_nb')
  call object_provide ('dexp_to_bas_nb')
  call object_provide ('zex')

  do dexp_to_bas_i = 1, dexp_to_bas_nb (dexp_i)
    bas_i = dexp_to_bas (dexp_i)%row (dexp_to_bas_i)
    zex (bas_i, 1) = zex (bas_i, 1) + d
  enddo ! dexp_to_bas_i

  call object_modified ('zex')
  call distinct_radial_bas

! update orbital coefficients
  select case (trim(basis_functions_varied))
  case ('unnormalized')
   call coef_orb_on_norm_basis_from_coef (1)
  case ('normalized')
   call coef_from_coef_orb_on_norm_basis (1)
  case ('orthonormalized')
   call coef_from_coef_orb_on_ortho_basis (1)
  case default
   call die (here, 'unknown case >'+trim(basis_functions_varied)+'< for basis_functions_varied.')
  end select

  end subroutine move_zex

! ==============================================================================
  subroutine deloc_lnexp_bld
! ------------------------------------------------------------------------------
! Description   : logarithm derivatives of Eloc with respect to logarithm of basis exponents
!
! Created       : J. Toulouse, 23 May 2007
! ------------------------------------------------------------------------------
  include 'modules.h'
  implicit none

! local
  integer dexp_i, bas_i

! header
  if (header_exe) then

   call object_create ('deloc_lnexp', deloc_lnexp_index)

   call object_needed ('param_exp_nb')
   call object_needed ('deloc_exp')
   call object_needed ('dexp_to_bas')
   call object_needed ('zex')

   return

  endif

! begin

! allocations
  call object_alloc ('deloc_lnexp', deloc_lnexp, param_exp_nb)

  do dexp_i = 1, param_exp_nb
     bas_i = dexp_to_bas (dexp_i)%row (1)
     deloc_lnexp (dexp_i) = deloc_exp (dexp_i) * zex (bas_i, iwf)
  enddo ! dexp_i

  end subroutine deloc_lnexp_bld

end module deriv_exp_mod
