module linearresponse_mod

  use all_tools_mod
  use opt_common_mod
  use dets_mod
  use vmc_mod
  use deriv_mod
  use opt_lin_mod

  ! Declaration of global variables and default values

  real(dp), allocatable           :: amat_av(:,:)
  real(dp), allocatable           :: amat_av_err(:,:)

  contains

!===========================================================================
  subroutine linearresponse_menu
!---------------------------------------------------------------------------
! Description : menu for linear-response calculations
!
! Created     : J. Toulouse, 18 May 2016
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

  ! local
  character(len=max_string_len_rout), save :: lhere = 'linearresponse_menu'

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
      write(6,'(a)') 'linearresponse'
      write(6,'(a)') 'end'

      !  case ('slater_mat_ex_trans_inv_sm')
      !  call get_next_value (l_slater_mat_ex_trans_inv_sm)

    case ('end')
      exit

    case default
      call die (lhere, 'unknown keyword >'+trim(word)+'<.')
    end select

  enddo ! end loop over menu lines

  l_opt_orb=.true.
  write(6,*)
  write(6,'(3a)') ' Orbital optimization information:'
  call object_provide ('param_orb_nb')
  call object_provide ('det_ex_unq_up_nb')
  call object_provide ('orb_opt_last_lab')
  write(6,'(a,i8)') ' Number of computed orbitals will be ', orb_opt_last_lab

  l_opt_exp=.true.
  write(6,*)
  write(6,'(3a)') ' Exponent optimization information:'
  call object_provide ('param_exp_nb')

  !l_opt_jas=.true.
  !call object_provide ('param_pjas_nb')
  write(6,*)
  !write(6,'(a,i5)') ' Number of Jastrow parameters:   ', nparmj
  !write(6,'(a,i5)') ' Number of periodic Jastrow parameters: ', param_pjas_nb

  l_opt_csf=.true.
  call object_provide ('nparm')
  call object_provide ('nparmcsf')
  call object_provide ('param_nb')
  write(6,'(a,i5)') ' Number of CSF parameters:       ', nparmcsf
  write(6,'(a,i5)') ' Number of orbital parameters:   ', param_orb_nb
  write(6,'(a,i5)') ' Number of exponent parameters:  ', param_exp_nb
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
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none

  ! local
  character(len=max_string_len_rout), save :: lhere = 'linearresponse'

  ! begin
  write(6,*)
  write(6,'(a)') '*********************** LINEAR-RESPONSE CALCULATION **************************'

  ! initializations
  if (.not. l_mode_vmc) then 
    call die (lhere, 'mode must be vmc for linear-response calculations')
  endif 
  call vmc_init

  call object_average_request('dpsi_av')
  call object_average_request('dpsi_eloc_av')
  call object_average_request('dpsi_dpsi_eloc_av')
  call object_average_request('dpsi_dpsi_av')
  call object_average_request('deloc_av')
  call object_average_request('dpsi_deloc_av')

  call object_error_request('amat_av_err')
  call vmc
  run_done=.true.

  end subroutine linearresponse


!===========================================================================
  subroutine  amat_av_bld
!---------------------------------------------------------------------------
! Description : 
!
! Created     : B. Mussard, Mon 06 Jun 2016 11:37:02 AM EDT
! Modified    :
!---------------------------------------------------------------------------
  include 'modules.h'
  implicit none
  
  if (header_exe) then
    call object_create('amat_av')

    call object_needed('param_nb')

    call object_needed('eloc_av')
    call object_needed('dpsi_av')
    call object_needed('dpsi_eloc_av')
    call object_needed('dpsi_dpsi_eloc_av')
    call object_needed('dpsi_dpsi_covar')
    call object_needed('dpsi_deloc_covar')

    return
  endif

  call object_alloc ('amat_av',amat_av,param_nb,param_nb)
  call object_alloc ('amat_av_err',amat_av_err,param_nb,param_nb)
  call object_associate ('amat_av',amat_av,param_nb,param_nb)
  call object_associate ('amat_av_err',amat_av_err,param_nb,param_nb)

  write(6,*) 'THIS IS HOW I WOULD DO A'

  end subroutine  amat_av_bld

end module linearresponse_mod
