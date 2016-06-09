module opt_nwt_mod

  use all_tools_mod
  use deriv_mod
  use opt_lin_mod

! Declaration of global variables and default values
  character(len=max_string_len)  :: hessian_type = 'tu'


  real(dp)                       :: beta_sor = 0.d0
  real(dp), allocatable          :: sh_sor (:,:)
  real(dp), allocatable          :: g_sor (:,:)
  real(dp), allocatable          :: hess_sor (:,:)
  real(dp), allocatable          :: hess_lzr (:,:)
  real(dp), allocatable          :: hess_uf (:,:)
  real(dp), allocatable          :: hess_tu (:,:)
  real(dp), allocatable          :: hess_lin (:,:)
  real(dp), allocatable          :: hess_nwt_energy (:,:)
  real(dp), allocatable          :: hess_nwt_variance (:,:)
  real(dp), allocatable          :: hess_nwt (:,:)
  real(dp), allocatable          :: hess_nwt_err (:,:)
  real(dp), allocatable          :: hess_nwt_sav (:,:)
  real(dp), allocatable          :: hess_nwt_inv (:,:)
  real(dp), allocatable          :: hess_nwt_eigvec (:,:)
  real(dp), allocatable          :: hess_nwt_eigval (:)
  real(dp)                       :: hess_nwt_eigval_min
  real(dp)                       :: hess_nwt_eigval_max
  real(dp), allocatable          :: hess_nwt_stab (:,:)
  real(dp), allocatable          :: hess_nwt_stab_inv (:,:)
  real(dp)                       :: hess_nwt_norm
  real(dp)                       :: hess_nwt_norm_err
  real(dp), allocatable          :: delta_nwt (:)


  contains

!===========================================================================
  subroutine opt_nwt_menu
!---------------------------------------------------------------------------
! Description : menu for newton optimization
!
! Created     : J. Toulouse, 24 Apr 2006
!---------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  character(len=max_string_len_rout), save :: lhere = 'opt_nwt_menu'

! begin

! loop over menu lines
  do
  call get_next_word (word)

  select case (trim(word))
  case ('help')
   write(6,'(a)') 'HELP for Newton optimization menu:'
   write(6,'(a)') 'newton'
   write(6,'(a)') ' hessian = lzr|uf|tu|linear|sorella'
   write(6,'(a)') ' beta = [real] beta parameter for Sorella Hessian'
   write(6,'(a)') 'end'

  case('hessian')
   call get_next_value (hessian_type)

  case ('beta')
   call get_next_value (beta_sor)
   write(6,*) trim(lhere),': beta parameter for Sorella Hessian = ',beta_sor

  case ('end')
   exit

  case default
   call die(lhere, 'unknown word >'+trim(word)+'<.')
  end select

  enddo ! end loop over menu lines

  write(6,'(2a)') trim(lhere),' optimization will be done by Newton method'
  write(6,'(3a)') trim(lhere),' with hessian type = ', trim(hessian_type)

  select case(trim(hessian_type))
   case ('lzr')
   case ('uf')
   case ('tu')
   case ('linear')
   case ('sorella')
   case default
    call die (lhere, 'unknown hessian type >'+trim(hessian_type)+'<.')
  end select

  end subroutine opt_nwt_menu

! ==============================================================================
  subroutine sh_sor_bld
! ------------------------------------------------------------------------------
! Description   : Sh matrix of Sorella, PRB 71 241103 (2005).
!
! Created       : J. Toulouse, 19 Sep 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer i, j

! header
  if (header_exe) then

   call object_create ('sh_sor')

   call object_needed ('param_nb')
   call object_needed ('dpsi_deloc_covar')

   return

  endif

! begin
  call object_alloc ('sh_sor', sh_sor, param_nb, param_nb)

  do i = 1, param_nb
    do j = i, param_nb

      sh_sor (i,j) = dpsi_deloc_covar (i,j) + dpsi_deloc_covar (j,i)

       if (i /= j) then
         sh_sor (j,i) = sh_sor (i,j)
       endif
    enddo
  enddo

  end subroutine sh_sor_bld

! ==============================================================================
  subroutine g_sor_bld
! ------------------------------------------------------------------------------
! Description   : G matrix of Sorella, PRB 71 241103 (2005).
!
! Created       : J. Toulouse, 19 Sep 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer i, j

! header
  if (header_exe) then

   call object_create ('g_sor')

   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('dpsi_dpsi_eloc_av')
   call object_needed ('dpsi_dpsi_av')
   call object_needed ('dpsi_av')
   call object_needed ('eloc_av')
   call object_needed ('gradient')

   return

  endif

! begin
  call object_alloc ('g_sor', g_sor, param_nb, param_nb)

  do i = 1, param_nb
    do j = i, param_nb

      g_sor (i,j) = 2.d0 * ( dpsi_dpsi_eloc_av (param_pairs(i,j)) - dpsi_dpsi_av (param_pairs(i,j)) * eloc_av  &
                            - dpsi_av (i) * gradient (j)/2.d0 - dpsi_av (j) * gradient (i)/2.d0 )

       if (i /= j) then
         g_sor (j,i) = g_sor (i,j)
       endif
    enddo
  enddo

  end subroutine g_sor_bld

! ==============================================================================
  subroutine hess_sor_bld
! ------------------------------------------------------------------------------
! Description   : Hessian matrix of Sorella, PRB 71 241103 (2005).
!
! Created       : J. Toulouse, 19 Sep 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer i, j

! header
  if (header_exe) then

   call object_create ('hess_sor')

   call object_needed ('sh_sor')
   call object_needed ('g_sor')

   return

  endif

! begin
  call object_alloc ('hess_sor', hess_sor, param_nb, param_nb)

  do i = 1, param_nb
    do j = i, param_nb

      hess_sor (i,j) = sh_sor (i,j)  + (1.d0 + beta_sor)* g_sor (i,j)

       if (i /= j) then
         hess_sor (j,i) = hess_sor (i,j)
       endif
    enddo
  enddo

  end subroutine hess_sor_bld

! ==============================================================================
  subroutine hess_lzr_bld
! ------------------------------------------------------------------------------
! Description   : Hessian  of energy from Lin, Zhang, Rappe, JCP 112, 2650 (2000).
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer i, j

! header
  if (header_exe) then

   call object_create ('hess_lzr')

   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('d2psi_eloc_av')
   call object_needed ('d2psi_av')
   call object_needed ('dpsi_dpsi_eloc_av')
   call object_needed ('dpsi_dpsi_av')
   call object_needed ('eloc_av')
   call object_needed ('dpsi_av')
   call object_needed ('gradient')
   call object_needed ('dpsi_deloc_av')

   return

  endif

! begin
  call object_alloc ('hess_lzr', hess_lzr, param_nb, param_nb)


  do i = 1, param_nb
    do j = i, param_nb

       hess_lzr (i,j) = 2.d0 * ( d2psi_eloc_av (param_pairs(i,j)) - d2psi_av (param_pairs(i,j)) * eloc_av   &
             + dpsi_dpsi_eloc_av (param_pairs(i,j)) -  dpsi_dpsi_av (param_pairs(i,j)) * eloc_av    &
             - dpsi_av (i) * gradient (j) - dpsi_av (j) * gradient (i)                                          &
             + dpsi_deloc_av (i,j) )

       if (i /= j) then
         hess_lzr (j,i) = hess_lzr (i,j)
       endif

    enddo
  enddo

  end subroutine hess_lzr_bld

! ==============================================================================
  subroutine hess_uf_bld
! ------------------------------------------------------------------------------
! Description   : Hessian  of energy from Umrigar, Filippi, PRL, 2005.
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer i, j

! header
  if (header_exe) then

   call object_create ('hess_uf')

   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('d2psi_eloc_av')
   call object_needed ('d2psi_av')
   call object_needed ('dpsi_dpsi_eloc_av')
   call object_needed ('dpsi_dpsi_av')
   call object_needed ('eloc_av')
   call object_needed ('dpsi_av')
   call object_needed ('gradient')
   call object_needed ('dpsi_deloc_covar')
!   call object_needed ('dpsi_nwt_eloc_covar')  !instead of gradient
!   call object_needed ('d2eloc_nwt_av') !

   return

  endif

! begin
  call object_alloc ('hess_uf', hess_uf, param_nb, param_nb)

  do i = 1, param_nb
    do j = i, param_nb

       hess_uf (i,j) = 2.d0 * ( d2psi_eloc_av (param_pairs(i,j)) - d2psi_av (param_pairs(i,j)) * eloc_av   &
             + dpsi_dpsi_eloc_av (param_pairs(i,j)) -  dpsi_dpsi_av (param_pairs(i,j)) * eloc_av    &
             - dpsi_av (i) * gradient (j) - dpsi_av (j) * gradient (i))     &
             +  dpsi_deloc_covar (i,j) +  dpsi_deloc_covar (j,i)

! alternative expression Hessian
!       hess_nwt (i,j) = 2.d0 * ( d2psi_eloc_av (param_pairs(i,j)) - d2psi_av (param_pairs(i,j)) * eloc_av   &
!             + dpsi_dpsi_eloc_av (param_pairs(i,j)) -  dpsi_dpsi_av (param_pairs(i,j)) * eloc_av    &
!             - 2.d0 *dpsi_av (i) * dpsi_nwt_eloc_covar (j) - 2.d0 * dpsi_av (j) * dpsi_nwt_eloc_covar (i)      &
!             +  dpsi_deloc_covar (i,j) +  dpsi_deloc_covar (j,i)  )   &
!             + d2eloc_nwt_av (param_pairs(i,j)) !!!!!!!!!!!!!!!

       if (i /= j) then
         hess_uf (j,i) = hess_uf (i,j)
       endif

    enddo
  enddo

  end subroutine hess_uf_bld

! ==============================================================================
  subroutine hess_tu_bld
! ------------------------------------------------------------------------------
! Description   : Hessian  of energy, with rescaling from Toulouse, Umrigar, JCP.
! Description   : to be check, probably not correct
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer i, j
  real(dp) sumD, sumBpD, sumBpD_over_sumD

! header
  if (header_exe) then

   call object_create ('hess_tu')

   call object_needed ('param_nb')
   call object_needed ('param_pairs')
   call object_needed ('d2psi_eloc_av')
   call object_needed ('d2psi_av')
   call object_needed ('dpsi_dpsi_eloc_av')
   call object_needed ('dpsi_dpsi_av')
   call object_needed ('eloc_av')
   call object_needed ('dpsi_av')
   call object_needed ('gradient')
   call object_needed ('dpsi_deloc_covar')

   return

  endif

! begin
  call object_alloc ('hess_tu', hess_tu, param_nb, param_nb)

  sumD = 0.d0
  sumBpD = 0.d0

   do i = 1, param_nb
    do j = i, param_nb
       if (i > nparmcsf .and. i <= nparmj .and. j > nparmcsf .and. j <= nparmj) then
        sumD = sumD + dpsi_deloc_covar (i,j) +  dpsi_deloc_covar (j,i)

        sumBpD = sumBpD + 2.d0 * (2.d0 * (dpsi_dpsi_eloc_av (param_pairs(i,j)) -  dpsi_dpsi_av (param_pairs(i,j)) * eloc_av )  &
                - dpsi_av (i) * gradient (j) - dpsi_av (j) * gradient (i))     &
                +  dpsi_deloc_covar (i,j) +  dpsi_deloc_covar (j,i)
        endif
    enddo
   enddo

   sumBpD_over_sumD = sumBpD/sumD

  do i = 1, param_nb
    do j = i, param_nb

       if (i > nparmcsf .and. i <= nparmj .and. j > nparmcsf .and. j <= nparmj) then
        hess_tu (i,j) = 2.d0 * (d2psi_eloc_av (param_pairs(i,j)) - d2psi_av (param_pairs(i,j)) * eloc_av   &
             - (dpsi_dpsi_eloc_av (param_pairs(i,j)) -  dpsi_dpsi_av (param_pairs(i,j)) * eloc_av) )  &
             + sumBpD_over_sumD * (   &
             +  dpsi_deloc_covar (i,j) +  dpsi_deloc_covar (j,i) )
       else
         hess_tu (i,j) = 2.d0 * (d2psi_eloc_av (param_pairs(i,j)) - d2psi_av (param_pairs(i,j)) * eloc_av   &
             + dpsi_dpsi_eloc_av (param_pairs(i,j)) -  dpsi_dpsi_av (param_pairs(i,j)) * eloc_av    &
             - dpsi_av (i) * gradient (j) - dpsi_av (j) * gradient (i))     &
             +  dpsi_deloc_covar (i,j) +  dpsi_deloc_covar (j,i)
       endif

       if (i /= j) then
         hess_tu (j,i) = hess_tu (i,j)
       endif
    enddo
  enddo

  end subroutine hess_tu_bld

! ==============================================================================
  subroutine hess_lin_bld
! ------------------------------------------------------------------------------
! Description   : Hessian  of energy from the linear method
! Description   : It must essentially equivalent to Sorella's Hessian with beta=0
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

  integer i, j

! header
  if (header_exe) then

   call object_create ('hess_lin')

   call object_needed ('param_nb')
   call object_needed ('eloc_av')
   call object_needed ('ham_lin')
   call object_needed ('ovlp_lin')

   return

  endif

! begin
  call object_alloc ('hess_lin', hess_lin, param_nb, param_nb)

  do i = 1, param_nb
    do j = i, param_nb

       hess_lin (i,j) = 2.d0 * (ham_lin (1+i,1+j) - eloc_av * ovlp_lin (1+i,1+j))

       if (i /= j) then
         hess_lin (j,i) = hess_lin (i,j)
       endif
    enddo
  enddo

  end subroutine hess_lin_bld

! ==============================================================================
  subroutine hess_nwt_energy_bld
! ------------------------------------------------------------------------------
! Description   : Hessian of energy
!
! Created       : J. Toulouse, 15 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('hess_nwt_energy')

   return

  endif

! begin
  call object_alloc ('hess_nwt_energy', hess_nwt_energy, param_nb, param_nb)

  select case(trim(hessian_type))
   case ('lzr')
    call object_provide_by_index (hess_nwt_energy_bld_index, hess_lzr_index)
    hess_nwt_energy (:,:) = hess_lzr (:,:)
   case ('uf')
    call object_provide_by_index (hess_nwt_energy_bld_index, hess_uf_index)
    hess_nwt_energy (:,:) = hess_uf (:,:)
   case ('tu')
    call object_provide_by_index (hess_nwt_energy_bld_index, hess_tu_index)
    hess_nwt_energy (:,:) = hess_tu (:,:)
   case ('linear')
    call object_provide_by_index (hess_nwt_energy_bld_index, hess_lin_index)
    hess_nwt_energy (:,:) = hess_lin (:,:)
   case ('sorella')
    call object_provide_by_index (hess_nwt_energy_bld_index, hess_sor_index)
    hess_nwt_energy (:,:) = hess_sor (:,:)
   case default
    call die (here, 'unknown hessian type >'+trim(hessian_type)+'<.')
  end select

  end subroutine hess_nwt_energy_bld

! ==============================================================================
  subroutine hess_nwt_bld
! ------------------------------------------------------------------------------
! Description   : Hessian of linear combinaison of energy and variance
!
! Created       : J. Toulouse, 25 Apr 2007
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('hess_nwt')

   call object_needed ('hess_nwt_energy')
   call object_needed ('p_var')

   return

  endif

! begin
  call object_alloc ('hess_nwt', hess_nwt, param_nb, param_nb)
  call object_alloc ('hess_nwt_err', hess_nwt_err, param_nb, param_nb)

  if (p_var /= 0.d0) then
    call object_provide_by_index (hess_nwt_bld_index, hessian_variance_index)
    hess_nwt = (1.d0 - p_var) * hess_nwt_energy + p_var * hessian_variance
  else
    hess_nwt = hess_nwt_energy
  endif

  end subroutine hess_nwt_bld

! ==============================================================================
  subroutine hess_nwt_eig_bld
! ------------------------------------------------------------------------------
! Description   : eigenvalues of hess_nwt
!
! Created       : J. Toulouse, 21 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

 ! integer i

! header
  if (header_exe) then

   call object_create ('hess_nwt_eigvec')
   call object_create ('hess_nwt_eigval')
   call object_create ('hess_nwt_eigval_min')
   call object_create ('hess_nwt_eigval_max')

   call object_needed ('param_nb')
   call object_needed ('hess_nwt')

   return

  endif

! begin
  call object_alloc ('hess_nwt_eigvec', hess_nwt_eigvec, param_nb, param_nb)
  call object_alloc ('hess_nwt_eigval', hess_nwt_eigval, param_nb)
  call object_associate ('hess_nwt_eigval_min', hess_nwt_eigval_min)

  call eigensystem (hess_nwt, hess_nwt_eigvec, hess_nwt_eigval, param_nb)

!  write(6,*) trim(here),': eigenvalues of Hessian:'
!  do i = 1, param_nb
!   write(6,'(2a,i3,a,e)') trim(here),': eigenvalue # ',i,'  ', hess_nwt_eigval (i)
!  enddo

! minimal and maximal eigenvalues
  hess_nwt_eigval_min =  minval(hess_nwt_eigval)
  hess_nwt_eigval_max =  maxval(hess_nwt_eigval)

!  write(6,'(2a,1p200d9.1)') trim(here),': hess_nwt_eigval=',hess_nwt_eigval(:)
! write(6,'(2a,1p1d8.1,a,1p1d8.1)') trim(here),': hess_nwt_eigval_min=',hess_nwt_eigval_min,' hess_nwt_eigval_max=',hess_nwt_eigval_max

  end subroutine hess_nwt_eig_bld

! ==============================================================================
  subroutine hess_nwt_stab_bld
! ------------------------------------------------------------------------------
! Description   : stabilized hess_nwt
!
! Created       : J. Toulouse, 21 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

  integer i

! header
  if (header_exe) then

   call object_create ('hess_nwt_stab')

   call object_needed ('param_nb')
   call object_needed ('hess_nwt')
   call object_needed ('hess_nwt_eigval_min')
   call object_needed ('diag_stab')

   return

  endif

! begin
  call object_alloc ('hess_nwt_stab', hess_nwt_stab, param_nb, param_nb)

  hess_nwt_stab (:,:) = hess_nwt (:,:)

! stabilization by modifying the diagonal
  if (l_stab) then
   do i = 1, param_nb
    hess_nwt_stab (i,i) = hess_nwt (i,i) + max(-hess_nwt_eigval_min,0.d0) + diag_stab
   enddo
  endif

  end subroutine hess_nwt_stab_bld

! ==============================================================================
  subroutine hess_nwt_stab_inv_bld
! ------------------------------------------------------------------------------
! Description   :  inverse of stab hess_nwt
!
! Created       : J. Toulouse, 18 Jan 2006
! ------------------------------------------------------------------------------
  implicit none

  real(dp) threshold

! header
  if (header_exe) then

   call object_create ('hess_nwt_stab_inv')

   call object_needed ('param_nb')
   call object_needed ('hess_nwt_stab')

   return

  endif

! begin
  call object_alloc ('hess_nwt_stab_inv', hess_nwt_stab_inv, param_nb, param_nb)

! the negative eigenvalues are already treated, so we can use small threshold
  threshold = 1.d-10 ! to make the inverse stable wrt numerical noise
!  threshold = 1.d-6 ! to make the inverse stable wrt numerical noise
  call inverse_by_svd (hess_nwt_stab, hess_nwt_stab_inv, param_nb, threshold)

  end subroutine hess_nwt_stab_inv_bld

! ==============================================================================
  subroutine hess_nwt_norm_bld
! ------------------------------------------------------------------------------
! Description   :  Frobenius norm of Hessian
! Description   :  norm(H) = sqrt ( trace( transpose(H) * H) )
!
! Created       : J. Toulouse, 20 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('hess_nwt_norm')

   call object_needed ('hess_nwt')

   return

  endif

! begin
  call object_associate ('hess_nwt_norm', hess_nwt_norm)
  call object_associate ('hess_nwt_norm_err', hess_nwt_norm_err)

  hess_nwt_norm = dsqrt(trace(matmul(transpose(hess_nwt),hess_nwt)))

  end subroutine hess_nwt_norm_bld

! ==============================================================================
  subroutine delta_nwt_bld
! ------------------------------------------------------------------------------
! Description   :  variation of parameters for newton method
!
! Created       : J. Toulouse, 10 Jan 2006
! ------------------------------------------------------------------------------
  use all_modules_mod
  implicit none

! header
  if (header_exe) then

   call object_create ('delta_nwt')

   call object_needed ('param_nb')
   call object_needed ('gradient')
   call object_needed ('hess_nwt_stab_inv')

   return

  endif

! begin

! allocations
  call object_alloc ('delta_nwt', delta_nwt, param_nb)

! solve H x = -g
  delta_nwt = matmul(hess_nwt_stab_inv, -gradient)

! check
!  write(6,*) 'delta_nwt=',delta_nwt
!  write(6,*) 'gradient=',gradient
!  write(6,*) 'gradient_calc=',-matmul(hess_nwt_stab, delta_nwt)

  end subroutine delta_nwt_bld

end module opt_nwt_mod
